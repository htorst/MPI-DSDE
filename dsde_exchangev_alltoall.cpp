/*
 * Copyright (c) 2012 The Trustees of University of Illinois.
 *                    All rights reserved.
 *
 * Author(s): Torsten Hoefler <htor@illinois.edu>
 *
 */

/*
 * dsde_exchange_alltoall.c
 *
 *  Created on: Dec 8, 2011
 *      Author: htor
 */
#include "dsde_internal.h"

typedef struct {
  DSDE_Free_fn free_fn; ///< the free function for exchangev_alltoall
  char *rbuf; ///< the receive buffer
  std::vector<int> *lrranks; ///< local rranks
  std::vector<MPI_Aint> *lrdispls; ///< local rdispls
  std::vector<MPI_Aint> *lrsizes; ///< local rsizes
} handle_t;

static int free(void **handlev) {
  handle_t *handle = (handle_t*)(*handlev);
  free(handle->rbuf);
  delete(handle->lrdispls);
  delete(handle->lrranks);
  delete(handle->lrsizes);
  free(handle);
  *handlev = DSDE_HANDLE_NULL;
}

int DSDE_Exchangev_alltoall(
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
  int p; res = MPI_Comm_size(libcomm, &p);
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

  MPI_Alltoall(&ssizes[0], 1, MPI_INT, &rsizes[0], 1, MPI_INT, libcomm);

  /*std::cout << "["<<r<<"] ";
  for(int i=0; i<rsizes.size(); ++i) { //
    std::cout << rsizes[i] << " ";
  }
  std::cout << "\n";*/

  std::vector<MPI_Request> reqs(srankcount);
  for(int i=0; i<srankcount; ++i) {
    char *sbuf = ((char *) sendbuf) + (sdispls[i]*sndext);
    //std::cout << "["<<r<<"] sending " << sendcounts[i] << " vertices to " << sranks[i] <<" reqs: "<<i+1<<"\n";
    MPI_Isend(sbuf, sendcounts[i], sendtype, sranks[i], 999, libcomm, &reqs[i]);
  }


  state_data->lrranks = new std::vector<int>; // local rranks
  state_data->lrdispls = new std::vector<MPI_Aint>; // local rdispls
  state_data->lrsizes = new std::vector<MPI_Aint>; // local rsizes
  int cursize = 0; // size of receive buffer in elements
  state_data->rbuf=NULL;
  *rrankcount = 0;

  for(int i=0; i<p; i++) {
    if(rsizes[i]) {
      state_data->lrranks->push_back(i);
      state_data->lrdispls->push_back(cursize);
      state_data->lrsizes->push_back(rsizes[i]);
      (*rrankcount)++;
      cursize =  state_data->lrdispls->back() + rsizes[i];
      state_data->rbuf = (char*)realloc(state_data->rbuf, cursize);
      reqs.resize(reqs.size()+1);
      MPI_Irecv(state_data->rbuf + state_data->lrdispls->back()*rcvext, rsizes[i], recvtype, i, 999, libcomm, &reqs.back());
      //MPI_Recv(&(*recvbuf)[recvbuf->size()-1][0], rsizes[i], MPI_CHAR, i, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //std::cout << "["<<r<<"] receiving " << rsizes[i] << " vertices from " << i <<" reqs: "<<reqs.size()<<"\n";
    }
  }

  MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUSES_IGNORE);

  *recvbuf = state_data->rbuf;
  *rranks = &(*state_data->lrranks)[0];
  *rdispls = &(*state_data->lrdispls)[0];
  *recvcounts = &(*state_data->lrsizes)[0];

  return MPI_SUCCESS;
}
