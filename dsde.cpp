#include "dsde.h"



/* will hold duplicate of MPI_COMM_SELF, needed for DSDE_Memcpy */
/* TODO: this should not be a global, won't work if we have more than one DSDE library instance -- make this an attribute to comm */
static MPI_Comm dsde_comm_self = MPI_COMM_NULL;

DSDE_Handle DSDE_HANDLE_NULL = NULL;

enum handle_type {
  DSDE_HT_BUF = 1, /* the handle value itself is the pointer to the buffer to be freed */
  DSDE_HT_MAX      /* place all values before this one */
};

/* all handle structs must have a handle_type enum value as their very first field */
/* we don't really need to keep track of all of these fields, but we do so anyway for debugging */
typedef struct handle_dsde_exchange {
  enum handle_type type;
  void*            recvbuf;
  int              rrankcount;
  int*             rranks;
  MPI_Aint*        recvcounts;
  MPI_Aint*        rdispls;
} handle_dsde_exchange_t;


int DSDE_Exchange(
  void*  sendbuf, int  srankcount, int  sranks[], MPI_Aint  sendcounts[], MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void** recvbuf, int* rrankcount, int* rranks[], MPI_Aint* recvcounts[], MPI_Aint* rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle)
{
  /* TODO: need to handle case that Aint value is larger than int */

  /* execute sendrecv to ourself on comm_self */
  int p;
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  std::vector<int> ssizes(p,0), rsizes(p);
  for(int i=0; i < srankcount; i++) ssizes[ sranks[i] ] = sendcounts[i];
  
  MPI_Alltoall(&ssizes[0], 1, MPI_INT, &rsizes[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<MPI_Request> reqs(srankcount);

  MPI_Aint lb, extent;
  MPI_Type_get_extent(sendtype, &lb, &extent);

  for(int i=0; i < srankcount; ++i) {
    MPI_Isend((char*)sendbuf+extent*sdispls[i], sendcounts[i], sendtype, sranks[i], 999, comm, &reqs[i]);
  }

  assert(recvbuf    == NULL); // avoid memory leaks
  assert(rranks     == NULL); // avoid memory leaks
  assert(recvcounts == NULL); // avoid memory leaks
  assert(rdispls    == NULL); // avoid memory leaks

  /* TODO: need to adjust between send and receive types here, since counts and extents may be different? */

  int rbufsize = 0;
  int num_incoming = 0;
  for(int i=0; i<p; i++) if(rsizes[i]) {
    num_incoming++;
    rbufsize+=rsizes[i];
  }

  /* allocate all the memory we need in one shot, then we'll section it off */
  void* memory = (void*) malloc(
    sizeof(handle_dsde_exchange_t) +   /* room to hold the fields for our handle struct */
    rbufsize * extent +                /* space to hold incoming data */
    sizeof(int)      * num_incoming + /* list of source ranks */
    sizeof(MPI_Aint) * num_incoming + /* list of counts from each source */
    sizeof(MPI_Aint) * num_incoming   /* list of displs for each source */
  );

  /* now step through memory and section it off, use tmp as current pointer */
  char* tmp = (char*)memory;

  /* grab space to hold our structure */
  handle_dsde_exchange_t* fields = (handle_dsde_exchange_t*) tmp;
  tmp += sizeof(handle_dsde_exchange_t);

  /* set this handle type to DSDE_HT_BUF */
  fields->type = DSDE_HT_BUF;

  /* grab space to hold incoming data */
  fields->recvbuf = (void*) tmp;
  tmp += rbufsize * extent;

  /* record the number of ranks */
  fields->rrankcount = num_incoming;

  /* grab space to hold list of source ranks */
  fields->rranks = (int*) tmp;
  tmp += sizeof(int) * num_incoming;

  /* grab space to hold list of incoming counts */
  fields->recvcounts = (MPI_Aint*) tmp;
  tmp += sizeof(MPI_Aint) * num_incoming;

  /* grab space to hold list of displacements within recvbuf */
  fields->rdispls = (MPI_Aint*) tmp;
  tmp += sizeof(MPI_Aint) * num_incoming;

  /* set the caller's output parameters */
  *recvbuf    = fields->recvbuf;
  *rrankcount = fields->rrankcount;
  *rranks     = fields->rranks;
  *recvcounts = fields->recvcounts;
  *rdispls    = fields->rdispls;
  *handle     = (DSDE_Handle*) memory;

  int reqoffset = reqs.size() - 1;
  reqs.resize(reqs.size() + num_incoming);
  
  MPI_Aint offset = 0;
  int index = 0;
  for(int i=0; i<p; i++) if(rsizes[i]) {
    MPI_Irecv((char*)(*recvbuf) + extent*offset, rsizes[i], recvtype, i, 999, comm, &reqs[reqoffset++]);
    (*rranks)[index]     = i;
    (*recvcounts)[index] = rsizes[i];
    (*rdispls)[index]    = offset;
    offset += rsizes[i];
    index++;
  }

  MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUSES_IGNORE);

  return MPI_SUCCESS;
}

/* free resources associated with handle, if any */
int DSDE_Free(DSDE_Handle* handle)
{
  if (handle != DSDE_HANDLE_NULL) {
    int type = *(int*) handle;
    if (type == DSDE_HT_BUF) {
      /* in this case, the handle value points directly to the memory that we need to free */
      free(handle);
      *handle = DSDE_HANDLE_NULL;
    }
  }
  return MPI_SUCCESS;
}

/* copy memory from srcbuf to dstbuf using committed MPI datatypes */
int DSDE_Memcpy(
  void* dstbuf,       MPI_Aint dstcount, MPI_Datatype dsttype,
  const void* srcbuf, MPI_Aint srccount, MPI_Datatype srctype)
{
  /* get a dup of MPI_COMM_SELF */

  if (dsde_comm_self == MPI_COMM_NULL) {
    MPI_Comm_dup(MPI_COMM_SELF, &dsde_comm_self);
  }

  /* TODO: need to handle case that Aint value is larger than int */

  /* execute sendrecv to ourself on comm_self */
  MPI_Sendrecv((void*)srcbuf, srccount, srctype, 0, 999, dstbuf, dstcount, dsttype, 0, 999, dsde_comm_self, MPI_STATUS_IGNORE);


  return MPI_SUCCESS;
}

int DSDE_Reduce_scatter_block(
  void* sendbuf, int srankcount, int sranks[], MPI_Aint sdispls[],
  int* flag, void* recvbuf, MPI_Aint recvcount, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
  *flag = 0;
  return MPI_SUCCESS;
}

