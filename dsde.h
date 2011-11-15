#include <mpi.h>

extern "C" {

/* handle to an object that tracks internal resources allocated during a DSDE call,
 * active handles must be freed with call to DSDE_Free to free internal resources */
typedef DSDE_Handle void*;

/* we define a NULL handle as DSDE_HANDLE_NULL */
DSDE_Handle DSDE_HANDLE_NULL = NULL;

/* Collective over all procs in specified communicator in which each process specifies
 * which ranks it has data for, and as output, it receives a list of ranks that have
 * sent data to it along with pointers to that data.  Data is stored in internal DSDE
 * resources which must be freed by call to DSDE_Free to release resources assocaited
 * with output handle.
 *
 * DSDE_Handle handle;
 * DSDE_Exchange(..., &handle);
 * ... copy received data to proper location ...
 * DSDE_Free(&handle);
 */
int DSDE_Exchange(
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
  void*        sendbuf,    /* IN  - starting address of send buffer (choice) */
  int          srankcount, /* IN  - number of dest. processes (non-negative integer) */
  int          sranks[],   /* IN  - integer array (of length srankcount) of dest. processes */
  MPI_Aint     sdispls[],  /* IN  - integer array (of length srankcount) of displs to send to dests[i] */
  int*         flag,       /* OUT - true if data in recvbuf is valid (logical) */
  void*        recvbuf,    /* OUT - starting address of receive buffer (choice) */
  MPI_Aint     recvcount,  /* OUT - element count per block (non-negative integer) */
  MPI_Datatype datatype,   /* IN  - data type of elements of send and receive buffers (handle) */
  MPI_Op       op,         /* IN  - operation (handle) */
  MPI_Comm     comm        /* IN  - communicator (handle) */
);

} /* extern "C" */
