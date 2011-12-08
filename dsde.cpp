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

int DSDE_Exchange_alltoall(
  void*  sendbuf, int  srankcount, int  sranks[], MPI_Aint  sendcounts[], MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void** recvbuf, int* rrankcount, int* rranks[], MPI_Aint* recvcounts[], MPI_Aint* rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle) {


  /* TODO: actually, the following should all be centralized: */
  int res = 0, recvsize;
  MPI_Aint sndext, rcvext;
  MPI_Comm libcomm = comm; // we should comm_dup here and play the attribute game
  int p; res = MPI_Comm_size(comm, &p);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Comm_size() (%i)\n", res); return res; }
  int r; res = MPI_Comm_rank(comm, &r);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Comm_size() (%i)\n", res); return res; }
  res = MPI_Type_extent(sendtype, &sndext);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Type_extent() (%i)\n", res); return res; }
  res = MPI_Type_extent(recvtype, &rcvext);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Type_extent() (%i)\n", res); return res; }
  res = MPI_Type_size(recvtype, &recvsize);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Type_size() (%i)\n", res); return res; }
  if(rcvext != recvsize) {
    printf("recvext (%i) is not equal to recvsize (%i), does this make sense?\n", rcvext, recvsize);
  }
  /* end of centralization */

  std::vector<int> ssizes(p,0), rsizes(p);
  for(int i=0; i<srankcount; ++i) { // fill sendsizes
    ssizes[sranks[i]]=sendcounts[i];
  }

  MPI_Alltoall(&ssizes[0], 1, MPI_INT, &rsizes[0], 1, MPI_INT, comm);

  /*std::cout << "["<<r<<"] ";
  for(int i=0; i<rsizes.size(); ++i) { //
    std::cout << rsizes[i] << " ";
  }
  std::cout << "\n";*/

  std::vector<MPI_Request> reqs(srankcount);
  for(int i=0; i<srankcount; ++i) {
    char *sbuf = ((char *) sendbuf) + (sdispls[i]*sndext);
    //std::cout << "["<<r<<"] sending " << sendcounts[i] << " vertices to " << sranks[i] <<" reqs: "<<i+1<<"\n";
    MPI_Isend(sbuf, sendcounts[i], sendtype, sranks[i], 999, comm, &reqs[i]);
  }


  std::vector<int> *lrranks = new std::vector<int>; // local rranks
  std::vector<MPI_Aint> *lrdispls = new std::vector<MPI_Aint>; // local rdispls
  std::vector<MPI_Aint> *lrsizes = new std::vector<MPI_Aint>; // local rsizes
  int cursize = 0; // size of receive buffer in elements
  char *rbuf=NULL;
  *rrankcount = 0;

  for(int i=0; i<p; i++) {
    if(rsizes[i]) {
      lrranks->push_back(i);
      lrdispls->push_back(cursize);
      lrsizes->push_back(rsizes[i]);
      (*rrankcount)++;
      cursize =  lrdispls->back() + rsizes[i];
      rbuf = (char*)realloc(rbuf, cursize);
      reqs.resize(reqs.size()+1);
      MPI_Irecv(rbuf + lrdispls->back()*rcvext, rsizes[i], recvtype, i, 999, comm, &reqs.back());
      //MPI_Recv(&(*recvbuf)[recvbuf->size()-1][0], rsizes[i], MPI_CHAR, i, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //std::cout << "["<<r<<"] receiving " << rsizes[i] << " vertices from " << i <<" reqs: "<<reqs.size()<<"\n";
    }
  }

  MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUSES_IGNORE);

  *recvbuf = rbuf;
  *rranks = &(*lrranks)[0];
  *rdispls = &(*lrdispls)[0];
  *recvcounts = &(*lrsizes)[0];

  /*if(!r && neighbors.size()) {
    std::cout << "["<<r<<"] >" << recvs << "< ";
    for(int i=0; i<neighbors.size(); ++i) {
      std::cout << neighbors[i] << " ";
    }
    std::cout << "\n";
  }*/

  return MPI_SUCCESS;
}

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

