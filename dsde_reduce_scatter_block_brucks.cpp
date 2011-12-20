/* Code for homogeneous systems */

/* Executes reduction, inlines all data so works best with small messages */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <signal.h>
#include <string.h>
#include <sys/time.h>
#include "mpi.h"
#include <math.h>
#include "dsde.h"

/* function to print error messages, just throw away for now */
#define sparse_abort

/* In the Brucks implementation, we pack and forward data through intermediate ranks.
 * An individual message is packed into an element, and then a list of elements is
 * packed into a packet and forwarded on to the destination.  If the message data is
 * small enough, it is inlined into the element.  Otherwise, only the element header
 * is sent and the data is sent direct as a separate message.
 *
 * The packet format looks like the following:
 *
 * <packet_header> : <element_header> [element_data], <element_header> [element data], ...
 *
 * where the packet header consists of:
 */

typedef struct {
  int rank;    /* rank of final destination */
} rsbi_packet_header;

/* qsort integer compare (using first four bytes of structure) */
static int int_cmp_fn(const void* a, const void* b)
{
  return (int) (*(int*)a - *(int*)b);
}

/* pack send data into buf, which is allocated and returned */
static int sparse_pack(
  const void* sendbuf, int srankcount, const int sranks[], const MPI_Aint sdispls[], MPI_Aint count, MPI_Datatype datatype,
  void** outbuf, int* outbuf_size, MPI_Aint* outextent, MPI_Comm comm)
{
  int i;

  /* we'll copy these values over to the output parameters before returning */
  void* buf      = NULL;
  int buf_offset = 0;

  /* compute memory size of element */
  MPI_Aint true_lb, true_extent;
  MPI_Datatype contig;
  MPI_Type_contiguous(count, datatype, &contig);
  MPI_Type_get_true_extent(contig, &true_lb, &true_extent);
  MPI_Type_free(&contig);

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
        DSDE_Memcpy((char*)buf + buf_offset, count, datatype, sendptr, count, datatype);
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
  *outextent    = true_extent;

  return MPI_SUCCESS;
}

/* given our position within bruck's algorithm, and given our current send buffer
 * and an array of recv buffers, combine received data with data we've yet to send
 * into a single merge buffer in preparation for the next round */
static int sparse_merge_bufs(
  int rank, int ranks, int factor, int degree, /* position within brucks algorithm */
  const void* sendbuf, int sendsize, int nrecv, void* recvbufs[], const int recvsizes[], /* send buf and array of recv bufs */
  MPI_Aint extent, MPI_Aint count, MPI_Datatype datatype, MPI_Op op, /* parameters needed for reduction */
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

  /* merge received data with data we've yet to send */
  int remaining = 1;
  while (remaining) {
    /* TODO: execute this merge as a k-way merge with min-heap */

    /* scan through and identify lowest rank among our buffer and receive buffers */
    int min_rank = -1;

    /* first check our buffer */
    int send_rank = -1;
    while (sendoffset < sendsize && send_rank == -1) {
      packet = (rsbi_packet_header*) ((char*)sendbuf + sendoffset);
      int dest_rank = packet->rank;
      int relative_rank = (dest_rank - rank + ranks) % ranks;
      int relative_id = (relative_rank / factor) % degree;
      if (relative_id == 0) {
        /* we kept the data for this rank during this round, so consider this rank for merging */
        send_rank = dest_rank;
      } else {
        /* we sent this data for this rank to someone else during this step, so skip it */
        sendoffset += sizeof(rsbi_packet_header) + extent;
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
          payload_bufs[payload_count] = ((char*)sendbuf + sendoffset + sizeof(rsbi_packet_header));
          payload_count++;

          sendoffset += sizeof(rsbi_packet_header) + extent;
        }
      }

      /* merge data from each receive buffer that matches the min rank */
      for (i = 0; i < nrecv; i++) {
        if (recvoffsets[i] < recvsizes[i]) {
          packet = (rsbi_packet_header*) ((char*)recvbufs[i] + recvoffsets[i]);
          int dest_rank = packet->rank;
          if (dest_rank == min_rank) {
            payload_bufs[payload_count] = (char*)recvbufs[i] + recvoffsets[i] + sizeof(rsbi_packet_header);
            payload_count++;

            recvoffsets[i] += sizeof(rsbi_packet_header) + extent;
          }
        }
      }

      /* merge payloads into a single payload */
      if (payload_count > 0) {
        void* reducebuf = (char*)mergebuf + mergeoffset;
        DSDE_Memcpy(reducebuf, count, datatype, payload_bufs[0], count, datatype);
        for (i = 1; i < payload_count; i++) {
          MPI_Reduce_local(payload_bufs[i], reducebuf, count, datatype, op);
        }
        mergeoffset += extent;
      }
    } else {
      /* found no rank on this scan, so we must be done */
      remaining = 0;
    }
  }

  /* update outpu parameters */
  *out_mergesize = mergeoffset;

  return rc;
}

/* execute Bruck's index algorithm to exchange data */
static int sparse_brucks(
  void** inoutbuf, int* inoutbuf_size,
  MPI_Aint extent, MPI_Aint count, MPI_Datatype datatype, MPI_Op op,
  int degree, MPI_Comm comm)
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

  /* send and destination rank arrays, size arrays, and MPI request and status objects */
  int* dst         = NULL;   /* array of ranks this process will send to during a given round */
  int* src         = NULL;   /* array of ranks this process will receive from in a given round */
  int* send_bytes  = NULL;   /* array of number of bytes for each rank this process will send to */
  int* recv_bytes  = NULL;   /* array of number of bytes this process will receive from each src rank */
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

    send_bytes = (int*) scratch_temp;
    scratch_temp += sizeof(int) * max_messages;

    recv_bytes = (int*) scratch_temp;
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
      dst[i] = (rank + (i+1) * factor + ranks) % ranks;
      src[i] = (rank - (i+1) * factor + ranks) % ranks;

      /* set our send buffer locations, assume we may need to send our entire buffer to each destination */
      send_bufs[i] = (char*)tmp_send + i * buf_size;

      /* initialize our send sizes */
      send_bytes[i] = 0;
    }

    /* pack our send messages and count number of bytes in each */
    int buf_offset = 0;
    while (buf_offset < buf_size) {
      rsbi_packet_header* packet = (rsbi_packet_header*) ((char*)buf + buf_offset);
      int dest_rank = packet->rank;
      int dest_size = sizeof(rsbi_packet_header) + extent;
      int relative_rank = (dest_rank - rank + ranks) % ranks;
      int relative_id = (relative_rank / factor) % degree;
      if (relative_id > 0) {
        int index = relative_id - 1;
        memcpy((char*)send_bufs[index] + send_bytes[index], (char*)buf + buf_offset, dest_size);
        send_bytes[index] += dest_size;
      }
      buf_offset += dest_size;
    }

    /* exchange number of bytes for this round */
    for (i = 0; i < num_messages; i++) {
      MPI_Irecv(&recv_bytes[i], 1, MPI_INT, src[i], 0, comm, &req[i]);
    }
    for (i = 0; i < num_messages; i++) {
      MPI_Isend(&send_bytes[i], 1, MPI_INT, dst[i], 0, comm, &req[num_messages + i]);
    }

    /* eagerly send our non-zero messages (assumed to be relatively small messages) */
    int req2_count = 0;
    for (i = 0; i < num_messages; i++) {
      if (send_bytes[i] > 0) {
        MPI_Isend(send_bufs[i], send_bytes[i], MPI_BYTE, dst[i], 0, comm, &req2[req2_count]);
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
      tmp_recv_size += recv_bytes[i];
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
      if (tmp_recv_size > 0 && recv_bytes[i] > 0) {
        recv_bufs[i] = (char*)tmp_recv + tmp_recv_offset;
        tmp_recv_offset += recv_bytes[i];
      }
    }

    /* finally, recv each non-zero message, and wait on sends and receives of non-zero messages */
    for (i = 0; i < num_messages; i++) {
      if (recv_bytes[i] > 0) {
        MPI_Irecv(recv_bufs[i], recv_bytes[i], MPI_BYTE, src[i], 0, comm, &req2[req2_count]);
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

    /* merge data that we didn't send with data we just received */
    int mergebuf_size = 0;
    sparse_merge_bufs(
      rank, ranks, factor, degree,
      buf, buf_size, num_messages, recv_bufs, recv_bytes,
      extent, count, datatype, op,
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

static int sparse_unpack(
  int* flag, void* recvbuf, MPI_Aint count, MPI_Datatype datatype,
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
    void* ptr = (char*)buf + sizeof(rsbi_packet_header);
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

/* used to efficiently implement an alltoallv where each process only sends to a handful of other processes,
 * uses the indexing algorithm by Jehoshua Bruck et al, IEEE TPDS, Nov. 97 */
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

  /* variables to track temporaries between sparse_pack/unpack calls */
  void* buf        = NULL;
  int buf_size     = 0;
  MPI_Aint extent  = 0;

  /* pack our send data into buf,
   * which is allocated and returned by sparse_pack */
  tmp_rc = sparse_pack(
    sendbuf, srankcount, sranks, sdispls, count, datatype,
    &buf, &buf_size, &extent, comm
  );
  if (rc == MPI_SUCCESS) {
    rc = tmp_rc;
  }

  /* execute Bruck's index algorithm to exchange data */
  tmp_rc = sparse_brucks(
    &buf, &buf_size,
    extent, count, datatype, op,
    degree, comm
  );
  if (rc == MPI_SUCCESS) {
    rc = tmp_rc;
  }

  /* unpack data into user receive buffers,
   * frees buf and request array which were allocated in sparse_pack,
   * and allocate handle to be returned to caller */
  tmp_rc = sparse_unpack(flag, recvbuf, count, datatype, &buf, &buf_size, comm);
  if (rc == MPI_SUCCESS) {
    rc = tmp_rc;
  }

  return rc;
}
