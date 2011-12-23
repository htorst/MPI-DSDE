/*
 *  Created on: Dec 22, 2011
 *      Author: Adam Moody
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <signal.h>
#include <string.h>
#include <sys/time.h>
#include "mpi.h"
#include <math.h>

#include "dsde_internal.h"

/* Packs and reduces data through intermediate ranks using k-way dissemination algorithm.
 *
 * Features
 *  - Portable to hetergeneous systems (however only tested on homogeneous systems)
 *  - Destination ranks in send arrays may be provided in any order
 *  - MPI_PROC_NULL is a valid destination
 *  - Any (portable) datatype is supported
 *  - Any reduction op is supported 
 * Limitations
 *  - Inlines reduction data so only efficient with small data
 *  - If too many messages / ranks, then intermediate ranks could be overwhelmed with
 *    too much incoming data
 *  - The same destination can not be listed more than once in send arrays
 *
 * A tunable parameter in this algorithm is the "degree", which must be >= 2.
 * This algorithm runs in ceil(log_degree(N)) rounds, wherein each round, each
 * process sends and receives up to (degree-1) messages.
 *
 * Each message consists of a list of packets, where each packet contains a header
 * and data.  The header is just an integer that lists the destination rank, and the
 * data is the current reduction result for that destination.
 *
 * So a message looks like the following:
 *
 * <dest_rank_1> <reduction_data_1> : <dest_rank_2> <reduction_data_2> : ...
 * 
 * During each step of the algorithm, each process sends and receives lists from other
 * ranks.  It merges each incoming list into a list of data it has yet to send out.
 * During this merge, if more than one list has data for a given rank, these items are
 * reduced to a single item by applying the reduction operation on the data elements.
 * The lists are kept sorted by increasing destination rank so that merges are efficient.
 */

/* function to print error messages, just throw away for now */
#define sparse_abort

typedef struct {
  int rank; /* rank of final destination */
} rsbi_packet_header;

/* qsort integer compare (using first four bytes of structure) */
static int int_cmp_fn(const void* a, const void* b)
{
  return (int) (*(int*)a - *(int*)b);
}

/* pack send data into list, which is allocated and returned in outbuf */
static int sparse_pack(
  const void* sendbuf, int srankcount, const int sranks[], const MPI_Aint sdispls[],
  MPI_Aint count, MPI_Datatype datatype, MPI_Aint true_lb, MPI_Aint true_extent,
  void** outbuf, int* outbuf_size, MPI_Comm comm)
{
  int i;

  /* we'll copy these values over to the output parameters before returning */
  void* buf      = NULL;
  int buf_offset = 0;

  /* prepare our initial data for sending */
  if (srankcount > 0) {
    /* get number of ranks in communicator */
    int ranks;
    MPI_Comm_size(comm, &ranks);

    /* for each message we might send, allocate packet header and space to pack its data */
    int buf_maxsize = (sizeof(rsbi_packet_header) + true_extent) * srankcount;
    buf = malloc(buf_maxsize);
    if (buf == NULL) {
      sparse_abort(1, "Failed to allocate temporary send buffer @ %s:%d",
        __FILE__, __LINE__
      );
    }

    /* get extent of sendtype
     * (extent, not true extent, since we use it with sdispls to find start of each send buffer) */
    MPI_Aint sendlb, sendextent;
    MPI_Type_get_extent(datatype, &sendlb, &sendextent);

    /* allocate space to sort ranks */
    int* sort_list = (int*) malloc(2 * sizeof(int) * srankcount);
    if (sort_list == NULL) {
      sparse_abort(1, "Failed to allocate memory for sort list @ %s:%d",
        __FILE__, __LINE__
      );
    }

    /* prepare list for sorting */
    for (i = 0; i < srankcount; i++) {
      /* record rank and its original index within the sranks array */
      sort_list[i*2+0] = sranks[i];
      sort_list[i*2+1] = i;
    }

    /* sort ranks in ascending order */
    qsort(sort_list, srankcount, 2 * sizeof(int), &int_cmp_fn);

    /* prepare data for sending */
    int last_rank = MPI_PROC_NULL;
    for (i = 0; i < srankcount; i++) {
      /* get the rank this packet is headed to */
      int dest_rank = sort_list[i*2+0];
      if (dest_rank == last_rank && last_rank != MPI_PROC_NULL) {
        sparse_abort(1,"Destination rank %d specified multiple times @ %s:%d",
          dest_rank, __FILE__, __LINE__
        );
      }

      /* check that we have a valid rank */
      if (dest_rank >= 0 && dest_rank < ranks) {
        /* get rank's original index in user send arrays */
        int rank_index = sort_list[i*2+1];

        /* get pointer to start of buffer for this message */
        void* sendptr = (char*)sendbuf + sdispls[rank_index] * sendextent;

        /* allocate space for a packet header */
        rsbi_packet_header* packet = (rsbi_packet_header*) ((char*)buf + buf_offset);
        packet->rank = dest_rank;
        buf_offset += sizeof(rsbi_packet_header);

        /* pack our inlined data */
        DSDE_Memcpy((char*)buf + buf_offset - true_lb, count, datatype, sendptr, count, datatype);
        buf_offset += true_extent;
      } else if (dest_rank != MPI_PROC_NULL) {
        /* error, rank out of range */
        sparse_abort(1, "Invalid destination rank %d @ %s:%d",
          dest_rank, __FILE__, __LINE__
        );
      }
    }

    /* free the sort array */
    if (sort_list != NULL) {
      free(sort_list);
      sort_list = NULL;
    }
  }

  /* update output parameters */
  *outbuf       = buf;
  *outbuf_size  = buf_offset;

  return MPI_SUCCESS;
}

/* given our position within bruck's algorithm, and given our current send buffer
 * and an array of recv buffers, combine received data with data we've yet to send
 * into a single merge buffer in preparation for the next round */
static int sparse_merge_bufs(
  int rank, int ranks, int factor, int degree, /* position within brucks algorithm */
  const void* sendbuf, int sendsize, int nrecv, void* recvbufs[], const int recvsizes[], /* send buf and array of recv bufs */
  MPI_Aint count, MPI_Datatype datatype, MPI_Aint true_lb, MPI_Aint true_extent, MPI_Op op, /* parameters needed for reduction */
  int recvoffsets[], void* payload_bufs[], /* scratch space */
  void* mergebuf, int* out_mergesize) /* output buffer containing merged data */
{
  int i;
  int rc = MPI_SUCCESS;
  rsbi_packet_header* packet = NULL;

  /* initialize our offsets to point to first byte of each input buffer */
  int sendoffset  = 0;
  int mergeoffset = 0;
  for (i = 0; i < nrecv; i++) {
    recvoffsets[i] = 0;
  }

  /* TODO: execute this merge as a k-way merge with min-heap */

  /* merge received data with data we've yet to send */
  int remaining = 1;
  while (remaining) {
    /* scan through and identify lowest rank among our buffer and receive buffers */
    int min_rank = -1;

    /* first check our buffer */
    int send_rank = -1;
    while (sendoffset < sendsize && send_rank == -1) {
      /* get the rank this packet is headed to */
      packet = (rsbi_packet_header*) ((char*)sendbuf + sendoffset);
      int dest_rank = packet->rank;

      /* compute the destination rank relative to our own rank */
      int relative_rank = dest_rank - rank;
      if (relative_rank < 0) {
        relative_rank += ranks;
      }

      int relative_id = (relative_rank / factor) % degree;
      if (relative_id == 0) {
        /* we kept the data for this rank during this round, so consider this rank for merging */
        send_rank = dest_rank;
      } else {
        /* we sent this data for this rank to someone else during this step, so skip it */
        sendoffset += sizeof(rsbi_packet_header) + true_extent;
      }
    }
    if (send_rank != -1) {
      min_rank = send_rank;
    }

    /* now check each of our receive buffers */
    for (i = 0; i < nrecv; i++) {
      if (recvoffsets[i] < recvsizes[i]) {
        packet = (rsbi_packet_header*) ((char*)recvbufs[i] + recvoffsets[i]);
        int dest_rank = packet->rank;
        if (dest_rank < min_rank || min_rank == -1) {
          min_rank = dest_rank;
        }
      }
    }

    /* if we found a rank, merge the data in the new_data buffer */
    if (min_rank != -1) {
      /* initialize our packet header to include merged data */
      rsbi_packet_header* merge_packet = (rsbi_packet_header*) ((char*)mergebuf + mergeoffset);
      merge_packet->rank = min_rank;
      mergeoffset += sizeof(rsbi_packet_header);
      int payload_count = 0;

      /* merge data from the send buffer if the rank matches the min rank */
      if (sendoffset < sendsize) { 
        packet = (rsbi_packet_header*) ((char*)sendbuf + sendoffset);
        int dest_rank = packet->rank;
        if (dest_rank == min_rank) {
          payload_bufs[payload_count] = (char*)sendbuf + sendoffset + sizeof(rsbi_packet_header) - true_lb;
          payload_count++;

          sendoffset += sizeof(rsbi_packet_header) + true_extent;
        }
      }

      /* merge data from each receive buffer that matches the min rank */
      for (i = 0; i < nrecv; i++) {
        if (recvoffsets[i] < recvsizes[i]) {
          packet = (rsbi_packet_header*) ((char*)recvbufs[i] + recvoffsets[i]);
          int dest_rank = packet->rank;
          if (dest_rank == min_rank) {
            payload_bufs[payload_count] = (char*)recvbufs[i] + recvoffsets[i] + sizeof(rsbi_packet_header) - true_lb;
            payload_count++;

            recvoffsets[i] += sizeof(rsbi_packet_header) + true_extent;
          }
        }
      }

      /* merge payloads into a single payload */
      if (payload_count > 0) {
        void* reducebuf = (char*)mergebuf + mergeoffset - true_lb;
        DSDE_Memcpy(reducebuf, count, datatype, payload_bufs[0], count, datatype);
        for (i = 1; i < payload_count; i++) {
          MPI_Reduce_local(payload_bufs[i], reducebuf, count, datatype, op);
        }
        mergeoffset += true_extent;
      }
    } else {
      /* found no rank on this scan, so we must be done */
      remaining = 0;
    }
  }

  /* update output parameters, we just need to tell caller how big the merge buffer is */
  *out_mergesize = mergeoffset;

  return rc;
}

/* execute Bruck's index algorithm to exchange data */
static int sparse_brucks(
  void** inoutbuf, int* inoutbuf_size,
  MPI_Aint count, MPI_Datatype datatype, MPI_Aint true_lb, MPI_Aint true_extent, MPI_Op op,
  MPI_Datatype packed, int degree, MPI_Comm comm)
{
  int i;

  /* record input parameters */
  void* buf    = *inoutbuf;
  int buf_size = *inoutbuf_size;

  /* get our rank and the number of ranks in communicator */
  int rank, ranks;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &ranks);

  /* compute the maximum number of messages we'll send in each round */
  int max_messages = degree - 1;

  /* get true extent of a packed type (should include leading int and count consecutive copies of datatype) */
  MPI_Aint packed_lb, packed_extent;
  MPI_Type_get_true_extent(packed, &packed_lb, &packed_extent);

  /* TODO: for a fixed degree on a given communicator, many of these values won't change,
   * so an optimization would be to cache this on comm after we compute it the first time */

  /* send and destination rank arrays, size arrays, and MPI request and status objects */
  int* dst         = NULL;   /* array of ranks this process will send to during a given round */
  int* src         = NULL;   /* array of ranks this process will receive from in a given round */
  int* send_counts = NULL;   /* array of number of bytes for each rank this process will send to */
  int* recv_counts = NULL;   /* array of number of bytes this process will receive from each src rank */
  int* recv_offs   = NULL;   /* array to hold offset into each of the receive buffers for the merge step */
  void** send_bufs = NULL;   /* array of pointers to send buffer locations */
  void** recv_bufs = NULL;   /* array of pointers to receive buffer locations */
  MPI_Request* req   = NULL; /* array of MPI_Requests for outstanding messages when exchanging number of bytes */
  MPI_Request* req2  = NULL; /* array of MPI_Requests for outstanding messages when exchanging data */
  MPI_Status*  stat  = NULL; /* array of MPI_Status objects for number-of-byte messages */
  MPI_Status*  stat2 = NULL; /* array of MPI_Status objects for data exchange messages */
  void** payload_bufs = NULL; /* array of pointers to payload data to be merged */

  /* to be efficient, we allocate one big section of memory and then setup pointers within this block
   * for each of our internal variables */
  void* scratch = NULL;
  int scratch_size = (5 * sizeof(int) + 2 * sizeof(void*) + 4 * sizeof(MPI_Request) + 4 * sizeof(MPI_Status)) * max_messages + 
                     (1 * sizeof(void*)) * (max_messages + 1);
  if (scratch_size > 0) {
    scratch = (void*) malloc(scratch_size);
    if (scratch == NULL) {
      sparse_abort(1, "Failed to allocate memory for internal data structures @ %s:%d",
        __FILE__, __LINE__
      );
    }

    char* scratch_temp = (char*) scratch;

    dst = (int*) scratch_temp;
    scratch_temp += sizeof(int) * max_messages;

    src = (int*) scratch_temp;
    scratch_temp += sizeof(int) * max_messages;

    send_counts = (int*) scratch_temp;
    scratch_temp += sizeof(int) * max_messages;

    recv_counts = (int*) scratch_temp;
    scratch_temp += sizeof(int) * max_messages;

    recv_offs = (int*) scratch_temp;
    scratch_temp += sizeof(int) * max_messages;

    send_bufs = (void**) scratch_temp;
    scratch_temp += sizeof(void*) * max_messages;

    recv_bufs = (void**) scratch_temp;
    scratch_temp += sizeof(void*) * max_messages;

    payload_bufs = (void**) scratch_temp;
    scratch_temp += sizeof(void*) * (max_messages + 1);

    req = (MPI_Request*) scratch_temp;
    scratch_temp += sizeof(MPI_Request) * max_messages * 2;

    req2 = (MPI_Request*) scratch_temp;
    scratch_temp += sizeof(MPI_Request) * max_messages * 2;

    stat = (MPI_Status*) scratch_temp;
    scratch_temp += sizeof(MPI_Status) * max_messages * 2;

    stat2 = (MPI_Status*) scratch_temp;
    scratch_temp += sizeof(MPI_Status) * max_messages * 2;
  }

  /* execute Bruck's index algorithm to exchange data */
  int factor = 1;
  while (factor < ranks) {
    /* compute the number of messages to send in this round */
    int num_messages = 0;
    for (i = 1; i < degree; i++) {
      if (i * factor < ranks) {
        num_messages++;
      }
    }

    /* allocate a buffer to pack our send data in,
     * assume we need to send the full buffer to each destination in this round */
    void* tmp_send = NULL;
    if (buf_size > 0) {
      tmp_send = malloc(buf_size * num_messages);
      if (tmp_send == NULL) {
        sparse_abort(1, "Failed to allocate temporary pack buffer @ %s:%d",
          __FILE__, __LINE__
        );
      }
    }

    /* determine our source and destination ranks for this round and set our send buffer locations and send size */
    for (i = 0; i < num_messages; i++) {
      dst[i] = rank + (i+1) * factor;
      if (dst[i] >= ranks) {
        dst[i] -= ranks;
      }
      src[i] = rank - (i+1) * factor;
      if (src[i] < 0) {
        src[i] += ranks;
      }

      /* set our send buffer locations, assume we may need to send our entire buffer to each destination */
      send_bufs[i] = (char*)tmp_send + i * buf_size;

      /* initialize our send sizes */
      send_counts[i] = 0;
    }

    /* pack our send messages and count number of bytes in each */
    int buf_offset = 0;
    while (buf_offset < buf_size) {
      /* get the rank this packet is headed to */
      rsbi_packet_header* packet = (rsbi_packet_header*) ((char*)buf + buf_offset);
      int dest_rank = packet->rank;

      /* compute the destination rank relative to our own rank */
      int relative_rank = dest_rank - rank;
      if (relative_rank < 0) {
        relative_rank += ranks;
      }

      int relative_id = (relative_rank / factor) % degree;
      if (relative_id > 0) {
        int index = relative_id - 1;
        memcpy((char*)send_bufs[index] + send_counts[index] * packed_extent, (char*)buf + buf_offset, packed_extent);
        send_counts[index]++;
      }
      buf_offset += packed_extent;
    }

    /* exchange number of bytes for this round */
    for (i = 0; i < num_messages; i++) {
      MPI_Irecv(&recv_counts[i], 1, MPI_INT, src[i], 0, comm, &req[i]);
    }
    for (i = 0; i < num_messages; i++) {
      MPI_Isend(&send_counts[i], 1, MPI_INT, dst[i], 0, comm, &req[num_messages + i]);
    }

    /* eagerly send our non-zero messages (assumed to be relatively small messages) */
    int req2_count = 0;
    for (i = 0; i < num_messages; i++) {
      if (send_counts[i] > 0) {
        MPI_Isend(send_bufs[i], send_counts[i], packed, dst[i], 0, comm, &req2[req2_count]);
        req2_count++;
      }
    }

    /* wait for the number-of-bytes messages to complete */
    if (num_messages > 0) {
      MPI_Waitall(num_messages * 2, req, stat);
    }

    /* count total number of bytes we'll receive in this round */
    int tmp_recv_size = 0;
    for (i = 0; i < num_messages; i++) {
      tmp_recv_size += recv_counts[i] * packed_extent;
    }

    /* allocate buffer to hold all incoming bytes */
    void* tmp_recv = NULL;
    if (tmp_recv_size > 0) {
      tmp_recv = (void*) malloc(tmp_recv_size);
      if (tmp_recv == NULL) {
        sparse_abort(1, "Failed to allocate temporary receive buffer @ %s:%d",
          __FILE__, __LINE__
        );
      }
    }

    /* assign receive buffer addresses based on size of each message */
    int tmp_recv_offset = 0;
    for (i = 0; i < num_messages; i++) {
      recv_bufs[i] = NULL;
      if (tmp_recv_size > 0 && recv_counts[i] > 0) {
        recv_bufs[i] = (char*)tmp_recv + tmp_recv_offset;
        tmp_recv_offset += recv_counts[i] * packed_extent;
      }
    }

    /* finally, recv each non-zero message, and wait on sends and receives of non-zero messages */
    for (i = 0; i < num_messages; i++) {
      if (recv_counts[i] > 0) {
        MPI_Irecv(recv_bufs[i], recv_counts[i], packed, src[i], 0, comm, &req2[req2_count]);
        req2_count++;
      }
    }
    if (req2_count > 0) {
      MPI_Waitall(req2_count, req2, stat2);
    }

    /* free the temporary send buffer */
    if (tmp_send != NULL) {
      free(tmp_send);
      tmp_send = NULL;
    }

    /* allocate space to merge received data with our data */
    int mergebuf_maxsize = buf_size + tmp_recv_size;
    void* mergebuf = malloc(mergebuf_maxsize);
    if (mergebuf == NULL) {
      sparse_abort(1, "Failed to allocate merge buffer @ %s:%d",
        __FILE__, __LINE__
      );
    }

    /* convert recv counts to bytes for merge step */
    for (i = 0; i < num_messages; i++) {
      recv_counts[i] *= packed_extent;
    }

    /* merge data that we didn't send with data we just received */
    int mergebuf_size = 0;
    sparse_merge_bufs(
      rank, ranks, factor, degree,
      buf, buf_size, num_messages, recv_bufs, recv_counts,
      count, datatype, true_lb, true_extent, op,
      recv_offs, payload_bufs,
      mergebuf, &mergebuf_size
    );

    /* we've merged our receive data into the mergebuf, so free the receive buffer */
    if (tmp_recv != NULL) {
      free(tmp_recv);
      tmp_recv = NULL;
    }

    /* free our current buffer and reassign its pointer to point to the merge buffer */
    if (buf != NULL) {
      free(buf);
      buf = NULL;
    }
    buf           = mergebuf;
    buf_size      = mergebuf_size;
    mergebuf      = NULL;
    mergebuf_size = 0;

    /* go on to the next phase of the exchange */
    factor *= degree;
  }

  /* free temporary memory */
  if (scratch != NULL) {
    free(scratch);
    scratch = NULL;
  }

  /* update output parameters */
  *inoutbuf      = buf;
  *inoutbuf_size = buf_size;

  return 0;
}

/* unpack our final list into user receive buffer and free memory holding list */
static int sparse_unpack(
  int* flag, void* recvbuf,
  MPI_Aint count, MPI_Datatype datatype, MPI_Aint true_lb, MPI_Aint true_extent,
  void** inbuf, int* inbuf_size, MPI_Comm comm)
{
  int i;

  /* record input parameters */
  void* buf    = *inbuf;
  int buf_size = *inbuf_size;

  /* check whether any data made it to this rank */
  if (buf_size > 0) {
    /* got something, get our rank in comm */
    int rank;
    MPI_Comm_rank(comm, &rank);

    /* read packet header and check that data arrived at correct rank */
    rsbi_packet_header* packet = (rsbi_packet_header*) buf;
    int dest_rank = packet->rank;
    if (dest_rank != rank) {
      /* error, received data for a different rank (shouldn't happen) */
      sparse_abort(1, "Received data for rank %d @ %s:%d",
        dest_rank, __FILE__, __LINE__
      );
    }

    /* copy data to user's receive buffer and set flag to indicate that it's valid */
    void* ptr = (char*)buf + sizeof(rsbi_packet_header) - true_lb;
    DSDE_Memcpy(recvbuf, count, datatype, ptr, count, datatype);
    *flag = 1;
  } else {
    /* no data reached this process, set flag to 0 to indicate receive buffer is not valid */
    *flag = 0;
  }

  /* free off the result buffer (allocated in sparse_pack) */
  if (buf != NULL) {
    free(buf);
    buf = NULL;
  }

  /* update output parameters */
  *inbuf      = buf;
  *inbuf_size = 0;

  return MPI_SUCCESS;
}

/* Packs and reduces data through intermediate ranks using k-way dissemination algorithm. */
int DSDE_Reduce_scatter_block_brucks(
    const void* sendbuf, int srankcount, const int sranks[], const MPI_Aint sdispls[],
    int* flag, void* recvbuf, MPI_Aint count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
    int degree)
{
  int tmp_rc;
  int rc = MPI_SUCCESS;

  /* TODO: current implementation limits one message per destination from each source */

  /* check that we have a valid value for our degree */
  if (degree < 2) {
    sparse_abort(1, "Degree must be >= 2 @ %s:%d",
      __FILE__, __LINE__
    );
    return !MPI_SUCCESS;
  }

  /* build a new type that is count consecutive entries of user's type */
  MPI_Datatype contig;
  MPI_Type_contiguous(count, datatype, &contig);

  /* get lower bound and extent of count consecutive copies of datatype */
  MPI_Aint contig_lb, contig_extent;
  MPI_Type_get_extent(contig, &contig_lb, &contig_extent);

  /* get true lower bound and true extent of count consecutive copies of datatype */
  MPI_Aint contig_true_lb, contig_true_extent;
  MPI_Type_get_true_extent(contig, &contig_true_lb, &contig_true_extent);

  /* create a new type that represents count copies with true extent from true lower
   * bound, keep user-defined lower bound */
  MPI_Datatype contig_resized;
  MPI_Type_create_resized(contig, contig_lb, contig_true_extent, &contig_resized);

  /* build the datatype representing our packed messsage: <integer> <contig_resized>
   * need to be sure to look out for any leading holes in contig type */
  MPI_Datatype packed;
  int blocks[2] = {1, 1};
  MPI_Aint displs[2];
  MPI_Datatype types[2];
  MPI_Aint int_true_lb, int_true_extent;
  MPI_Type_get_true_extent(MPI_INT, &int_true_lb, &int_true_extent);
  displs[0] = 0;
  displs[1] = int_true_extent - contig_true_lb;
  types[0] = MPI_INT;
  types[1] = contig_resized;
  MPI_Type_create_struct(2, blocks, displs, types, &packed);
  MPI_Type_commit(&packed);

  /* variables to track temporaries between sparse_pack/unpack calls */
  void* buf    = NULL;
  int buf_size = 0;

  /* pack our send data into list,
   * which is allocated and returned by sparse_pack */
  tmp_rc = sparse_pack(
    sendbuf, srankcount, sranks, sdispls,
    count, datatype, contig_true_lb, contig_true_extent,
    &buf, &buf_size, comm
  );
  if (rc == MPI_SUCCESS) {
    rc = tmp_rc;
  }

  /* execute Bruck's index algorithm to exchange data and merge lists */
  tmp_rc = sparse_brucks(
    &buf, &buf_size,
    count, datatype, contig_true_lb, contig_true_extent, op,
    packed, degree, comm
  );
  if (rc == MPI_SUCCESS) {
    rc = tmp_rc;
  }

  /* unpack final list into user receive buffer, and free list */
  tmp_rc = sparse_unpack(
    flag, recvbuf,
    count, datatype, contig_true_lb, contig_true_extent,
    &buf, &buf_size, comm
  );
  if (rc == MPI_SUCCESS) {
    rc = tmp_rc;
  }

  /* TODO: cache these as attributes on the datatype as an optimization */

  /* free off our temporary types */
  MPI_Type_free(&contig_resized);
  MPI_Type_free(&contig);
  MPI_Type_free(&packed);

  return rc;
}
