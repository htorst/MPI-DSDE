#include "dsde_internal.h"

/* will hold duplicate of MPI_COMM_SELF, needed for DSDE_Memcpy */
/* TODO: this should not be a global, won't work if we have more than one DSDE library instance -- make this an attribute to comm */
static MPI_Comm dsde_comm_self = MPI_COMM_NULL;

/* set a NULL handle to be a NULL pointer */
DSDE_Handle DSDE_HANDLE_NULL = NULL;

/* copy memory from srcbuf to dstbuf using committed MPI datatypes */
int DSDE_Memcpy(
  void* dstbuf,       MPI_Aint dstcount, MPI_Datatype dsttype,
  const void* srcbuf, MPI_Aint srccount, MPI_Datatype srctype)
{
  /* TODO: not thread-safe */
  /* get a dup of MPI_COMM_SELF */
  if (dsde_comm_self == MPI_COMM_NULL) {
    MPI_Comm_dup(MPI_COMM_SELF, &dsde_comm_self);
  }

  /* TODO: need to handle case that Aint value is larger than int */
  /* execute sendrecv to ourself on comm_self */
  MPI_Sendrecv((void*)srcbuf, srccount, srctype, 0, 999, dstbuf, dstcount, dsttype, 0, 999, dsde_comm_self, MPI_STATUS_IGNORE);

  return MPI_SUCCESS;
}

/* free resources associated with handle, if any */
int DSDE_Free(DSDE_Handle* handle)
{
  if (handle != NULL && *handle != DSDE_HANDLE_NULL) {
    DSDE_Free_fn* fn = (DSDE_Free_fn*)(*handle);
    return (*fn)(handle);
  }
  return MPI_SUCCESS;
}

/* assumes that handle just points to one big block of memory that must be freed */
int DSDE_Free_single(DSDE_Handle* handle)
{
  if (handle != NULL && *handle != DSDE_HANDLE_NULL) {
    free(*handle);
    *handle = DSDE_HANDLE_NULL;
  }
  return MPI_SUCCESS;
}

int DSDE_Exchange(
  const void*  sendbuf, int  srankcount, const int  sranks[], MPI_Aint sendcount, const MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void**       recvbuf, int* rrankcount, int*       rranks[], MPI_Aint recvcount, MPI_Aint*       rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle)
{
  int rc = MPI_SUCCESS;

  /* TODO: need to handle case that Aint value is larger than int */

  /* TODO: determine which internal function to call */

  /* TODO: invoke internal function */

  /* for now, just invoke something */
  /* TODO: set degree based on message size, smaller messages may benefit from higher degree values */
  int degree = 2;
  rc = DSDE_Exchange_brucks_inline(
    sendbuf, srankcount, sranks, sendcount, sdispls, sendtype,
    recvbuf, rrankcount, rranks, recvcount, rdispls, recvtype,
    comm, handle,
    degree
  );

  return rc;
}

int DSDE_Exchangev(
  const void*  sendbuf, int  srankcount, const int  sranks[], const MPI_Aint  sendcounts[], const MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void**       recvbuf, int* rrankcount, int*       rranks[], MPI_Aint*       recvcounts[], MPI_Aint*       rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle)
{
  int rc = MPI_SUCCESS;

  /* TODO: need to handle case that Aint value is larger than int */

  /* TODO: determine which internal function to call */

  /* TODO: invoke internal function */

  /* for now, just invoke something */
  /* TODO: set degree based on message size, smaller messages may benefit from higher degree values */
  int degree = 2;
  rc = DSDE_Exchangev_brucks(
    sendbuf, srankcount, sranks, sendcounts, sdispls, sendtype,
    recvbuf, rrankcount, rranks, recvcounts, rdispls, recvtype,
    comm, handle,
    degree
  );

  return rc;
}

int DSDE_Reduce_scatter_block(
  void* sendbuf, int srankcount, int sranks[], MPI_Aint sdispls[],
  int* flag, void* recvbuf, MPI_Aint recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  int rc = MPI_SUCCESS;

  /* assume we don't get any data */
  *flag = 0;

  /* TODO: need to handle case that Aint value is larger than int */

  /* TODO: determine which internal function to call */

  /* TODO: invoke internal function */

  /* right now we only have one implementation so there's no choice */
  /* TODO: set degree based on message size, smaller messages may benefit from higher degree values */
  int degree = 2;
  rc = DSDE_Reduce_scatter_block_brucks(
    sendbuf, srankcount, sranks, sdispls,
    flag, recvbuf, recvcount, datatype, op, comm, degree
  );

  return rc;
}
