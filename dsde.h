#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <limits>
#include <iostream>
#include <ctime>
#include <set>
#include <queue>

extern "C" {

/* handle to an object that tracks internal resources allocated during a DSDE call,
 * active handles must be freed with call to DSDE_Free to free internal resources */
typedef void* DSDE_Handle;

/* we define a NULL handle as DSDE_HANDLE_NULL */
extern DSDE_Handle DSDE_HANDLE_NULL;

/* Copies data from one buffer to another on the same process using MPI
 * datatypes.  The number of basic elements specified by dstcount and
 * dsttype must be equal to the number of elements specified by srccount
 * and srctype, and both dsttype and srctype must be committed.
 * Otherwise, the call is erroneous.  */
int DSDE_Memcpy(
  void*        dstbuf,    /* OUT - buffer to copy data to */
  MPI_Aint     dstcount,  /* IN  - number of elements of type dsttype to be stored to dstbuf */
  MPI_Datatype dsttype,   /* IN  - datatype of elements in dstbuf */
  const void*  srcbuf,    /* IN  - source buffer to copy data from */ 
  MPI_Aint     srccount,  /* IN  - number of elements of type srctype to be copied from srcbuf */
  MPI_Datatype srctype    /* IN  - datatype of elements in srcbuf */
);

/* Collective over all procs in specified communicator in which each process specifies
 * which ranks it has data for, and as output, it receives a list of ranks that have
 * sent data to it along with pointers to that data.  Data is stored in internal DSDE
 * resources which must be freed by a call to DSDE_Free to release resources associated
 * with output handle.
 *
 * DSDE_Handle handle;
 * DSDE_Exchange(..., &handle);
 * ... copy received data to proper location ...
 * DSDE_Free(&handle);
 */
int DSDE_Exchange(
  const void*    sendbuf,       /* IN  - starting address of send buffer (choice) */
  int            srankcount,    /* IN  - number of dest. processes (non-negative integer) */
  const int      sranks[],      /* IN  - non-negative integer array (of length srankcount) of dest. processes */
  MPI_Aint       sendcount,     /* IN  - non-negative integer array (of length srankcount) of number of elements to send to sranks[i] */
  const MPI_Aint sdispls[],     /* IN  - integer array (of length srankcount) of displs to send to dests[i] */
  MPI_Datatype   sendtype,      /* IN  - data type of send buffer elements (handle) */
  void**         recvbuf,       /* OUT - starting address of receive buffer, allocated by lib (choice) */
  int*           rrankcount,    /* OUT - number of source processes (non-negative integer) */
  int*           rranks[],      /* OUT - non-negative integer array (of length rrankcount) of source processes */
  MPI_Aint       recvcount,     /* IN  - non-negative integer array (of length rrankcount) of number of elements recevied from rranks[i] */
  MPI_Aint*      rdispls[],     /* OUT - integer array (of length rrankcount) of displs of data received from rranks[i] */
  MPI_Datatype   recvtype,      /* IN  - data type of receive buffer elements (handle) */
  MPI_Comm       comm,          /* IN  - communicator (handle) */
  DSDE_Handle*   handle         /* OUT - DSDE resource (handle) */
);

int DSDE_Exchangev(
  void*        sendbuf,       /* IN  - starting address of send buffer (choice) */
  int          srankcount,    /* IN  - number of dest. processes (non-negative integer) */
  int          sranks[],      /* IN  - non-negative integer array (of length srankcount) of dest. processes */
  MPI_Aint     sendcounts[],  /* IN  - non-negative integer array (of length srankcount) of number of elements to send to sranks[i] */
  MPI_Aint     sdispls[],     /* IN  - integer array (of length srankcount) of displs to send to dests[i] */
  MPI_Datatype sendtype,      /* IN  - data type of send buffer elements (handle) */
  void**       recvbuf,       /* OUT - starting address of receive buffer, allocated by lib (choice) */
  int*         rrankcount,    /* OUT - number of source processes (non-negative integer) */
  int*         rranks[],      /* OUT - non-negative integer array (of length rrankcount) of source processes */
  MPI_Aint*    recvcounts[],  /* OUT - non-negative integer array (of length rrankcount) of number of elements recevied from rranks[i] */
  MPI_Aint*    rdispls[],     /* OUT - integer array (of length rrankcount) of displs of data received from rranks[i] */
  MPI_Datatype recvtype,      /* IN  - data type of receive buffer elements (handle) */
  MPI_Comm     comm,          /* IN  - communicator (handle) */
  DSDE_Handle* handle         /* OUT - DSDE resource (handle) */
);

/* frees resources internally allocated in call to DSDE_Exchange,
 * sets mem to DSDE_HANDLE_NULL */
int DSDE_Free(
  DSDE_Handle* handle /* INOUT - DSDE resource (handle) */
);

/* reduce scatter block in which each process specifies the set of ranks
 * for which it has data, flag is used to indicate whether the receive buffer
 * is valid since there may be 0 procs with data to contribute to the calling proc */
/* TODO: may be useful to change flag to count to record number of processes
 * contributing data */
int DSDE_Reduce_scatter_block(
  const void*    sendbuf,    /* IN  - starting address of send buffer (choice) */
  int            srankcount, /* IN  - number of dest. processes (non-negative integer) */
  const int      sranks[],   /* IN  - integer array (of length srankcount) of dest. processes */
  const MPI_Aint sdispls[],  /* IN  - integer array (of length srankcount) of displs to send to dests[i] */
  int*           flag,       /* OUT - true if data in recvbuf is valid (logical) */
  void*          recvbuf,    /* OUT - starting address of receive buffer (choice) */
  MPI_Aint       count,      /* IN  - element count per block (non-negative integer) */
  MPI_Datatype   datatype,   /* IN  - data type of elements of send and receive buffers (handle) */
  MPI_Op         op,         /* IN  - operation (handle) */
  MPI_Comm       comm        /* IN  - communicator (handle) */
);

/* -----------------------------------------
 * Move this section to an internal header not exposed to user
 * ----------------------------------------- */

/* function pointer to a DSDE_Free implementation that takes a pointer to a handle */
typedef int(*DSDE_Free_fn)(DSDE_Handle*);

/* assumes that handle just points to one big block of memory that must be freed */
int DSDE_Free_single(DSDE_Handle* handle);

/* this may be temp */
int DSDE_Exchange_alltoall(
  void*  sendbuf, int  srankcount, int  sranks[], MPI_Aint  sendcounts[], MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void** recvbuf, int* rrankcount, int* rranks[], MPI_Aint* recvcounts[], MPI_Aint* rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle);
int DSDE_Exchange_reduce_scatter(
  void*  sendbuf, int  srankcount, int  sranks[], MPI_Aint  sendcounts[], MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void** recvbuf, int* rrankcount, int* rranks[], MPI_Aint* recvcounts[], MPI_Aint* rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle);
int DSDE_Exchange_accumulate(
  void*  sendbuf, int  srankcount, int  sranks[], MPI_Aint  sendcounts[], MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void** recvbuf, int* rrankcount, int* rranks[], MPI_Aint* recvcounts[], MPI_Aint* rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle);
int DSDE_Exchange_ibarrier(
  void*  sendbuf, int  srankcount, int  sranks[], MPI_Aint  sendcounts[], MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void** recvbuf, int* rrankcount, int* rranks[], MPI_Aint* recvcounts[], MPI_Aint* rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle);

int DSDE_Exchangev_brucks(
  void*  sendbuf, int  srankcount, int  sranks[], MPI_Aint  sendcounts[], MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void** recvbuf, int* rrankcount, int* rranks[], MPI_Aint* recvcounts[], MPI_Aint* rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle);

int DSDE_Reduce_scatter_block_brucks(
    const void* sendbuf, int srankcount, const int sranks[], const MPI_Aint sdispls[],
    int* flag, void* recvbuf, MPI_Aint count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
    int degree);

} /* extern "C" */
