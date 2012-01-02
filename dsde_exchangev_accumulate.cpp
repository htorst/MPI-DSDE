/*
 * dsde_exchange_accumulate.cpp
 *
 *  Created on: Dec 8, 2011
 *      Author: htor
 */
#include "dsde_internal.h"

// TODO: this is horrible and not thread-safe, this should be hung off the comm attribute!
static int mpi2os_recvs;
static char mpi2os_initialized=0;
static MPI_Info info;
static MPI_Win win;

typedef struct {
  DSDE_Free_fn free_fn; ///< the free function for exchangev_alltoall
  char *rbuf; ///< the receive buffer
  std::vector<int> *lrranks; ///< local rranks
  std::vector<MPI_Aint> *lrdispls; ///< local rdispls
  std::vector<MPI_Aint> *lrsizes; ///< local rsizes
} handle_t;

static void prepare_mpi2os(MPI_Comm comm) {
  MPI_Info_create(&info);
  MPI_Info_set(info, (char*)std::string("no_locks").c_str(), (char*)std::string("true").c_str());
  MPI_Win_create(&mpi2os_recvs, sizeof(int), 1, info, comm, &win);
}

/* TODO: this is never called :-( -- should we have a DSDE_Finalize? I'd rater not call this in every DSDE_Free */
static void free_mpi2os() {
  MPI_Win_free(&win);
  MPI_Info_free(&info);
}

static int free(void **handlev) {
  handle_t *handle = (handle_t*)(*handlev);
  free(handle->rbuf);
  delete(handle->lrdispls);
  delete(handle->lrranks);
  delete(handle->lrsizes);
  free(handle);
  *handlev = DSDE_HANDLE_NULL;
}

/* MPI-2 one-sided: uses accumulate instead of reduce_scatter */
int DSDE_Exchangev_accumulate(
  const void*  sendbuf, int  srankcount, const int  sranks[], const MPI_Aint  sendcounts[], const MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void**       recvbuf, int* rrankcount, int*       rranks[], MPI_Aint*       recvcounts[], MPI_Aint*       rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle) {

  handle_t *state_data=(handle_t*)malloc(sizeof(handle_t)); ///< this struct contains the free function and pointers to all data that is returned
  state_data->free_fn = free;
  *handle = (DSDE_Handle*)state_data;

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

  if(!mpi2os_initialized) {
    prepare_mpi2os(libcomm);
  }

  std::vector<MPI_Request> reqs(srankcount);
  for(int i=0; i<srankcount; ++i) {
    char *sbuf = ((char *) sendbuf) + (sdispls[i]*sndext);
    //std::cout << "["<<r<<"] sending " << sendcounts[i] << " vertices to " << sranks[i] <<" reqs: "<<i+1<<"\n";
    //printf("[%i] sending %i elements to %i\n", r, sendcounts[i], sranks[i]);
    MPI_Isend(sbuf, sendcounts[i], sendtype, sranks[i], 999, libcomm, &reqs[i]);
  }

  // TODO: certainly not thread-safe!
  mpi2os_recvs=0;

  MPI_Win_fence(0, win);
  for(int i=0; i<srankcount; i++) { // fill senddests
    int x=1;
    MPI_Accumulate(&x,1,MPI_INT,sranks[i],0,1,MPI_INT,MPI_SUM,win);
  }
  MPI_Win_fence(0, win);

  const int recvs=mpi2os_recvs;

  state_data->lrranks = new std::vector<int>; // local rranks
  state_data->lrdispls = new std::vector<MPI_Aint>; // local rdispls
  state_data->lrsizes = new std::vector<MPI_Aint>; // local rsizes
  int cursize = 0; // size of receive buffer in elements
  state_data->rbuf=NULL; // local receive buffer pointer
  *rrankcount = 0; // input argument - number of receive ranks

  for(int i=0; i<recvs; ++i) {
    MPI_Status stat;
    MPI_Probe(MPI_ANY_SOURCE, 999, libcomm, &stat);
    int count; MPI_Get_count(&stat, recvtype, &count);
    state_data->lrranks->push_back(stat.MPI_SOURCE);
    state_data->lrdispls->push_back(cursize);
    state_data->lrsizes->push_back(count);
    (*rrankcount)++;
    cursize =  state_data->lrdispls->back() + count;
    state_data->rbuf = (char*)realloc(state_data->rbuf, cursize);
    reqs.resize(reqs.size()+1);

    //printf("[%i] receiving %i elements from %i\n", r, count, stat.MPI_SOURCE);
    MPI_Irecv(state_data->rbuf + state_data->lrdispls->back()*rcvext, count, recvtype, stat.MPI_SOURCE, 999, libcomm, &reqs.back());
  }

  MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUSES_IGNORE);

  *recvbuf = state_data->rbuf;
  *rranks = &(*state_data->lrranks)[0];
  *rdispls = &(*state_data->lrdispls)[0];
  *recvcounts = &(*state_data->lrsizes)[0];

  return MPI_SUCCESS;
}
