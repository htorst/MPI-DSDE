/*
 * dsde_exchange_reduce_scatter.cpp
 *
 *  Created on: Dec 8, 2011
 *      Author: htor
 */
#include "dsde_internal.h"

/* here we get only the processes that we receive from so that we have
 * to probe/malloc/receive the messages */
int DSDE_Exchangev_reduce_scatter(
  const void*  sendbuf, int  srankcount, const int  sranks[], const MPI_Aint  sendcounts[], const MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void**       recvbuf, int* rrankcount, int*       rranks[], MPI_Aint*       recvcounts[], MPI_Aint*       rdispls[], MPI_Datatype recvtype,
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


  std::vector<int> sdests(p,0), rcounts(p,1);
  for(int i=0; i<srankcount; ++i) { // fill senddests
    sdests[ sranks[i] ]=1;
  }

  int recvs;
  MPI_Reduce_scatter(&sdests[0], &recvs, &rcounts[0], MPI_INT, MPI_SUM, libcomm);

  /*std::cout << "["<<r<<"] {" << recvs << "} ";
  for(int i=0; i<sdests.size(); ++i) { //
    std::cout << sdests[i] << " ";
  }
  std::cout << "\n";*/

  std::vector<MPI_Request> reqs(srankcount);
  for(int i=0; i<srankcount; ++i) {
    char *sbuf = ((char *) sendbuf) + (sdispls[i]*sndext);
    //std::cout << "["<<r<<"] sending " << sendcounts[i] << " vertices to " << sranks[i] <<" reqs: "<<i+1<<"\n";
    //printf("[%i] sending %i elements to %i\n", r, sendcounts[i], sranks[i]);
    MPI_Isend(sbuf, sendcounts[i], sendtype, sranks[i], 999, libcomm, &reqs[i]);
  }

  std::vector<int> *lrranks = new std::vector<int>; // local rranks
  std::vector<MPI_Aint> *lrdispls = new std::vector<MPI_Aint>; // local rdispls
  std::vector<MPI_Aint> *lrsizes = new std::vector<MPI_Aint>; // local rsizes
  int cursize = 0; // size of receive buffer in elements
  char *rbuf=NULL; // local receive buffer pointer
  *rrankcount = 0; // input argument - number of receive ranks

  for(int i=0; i<recvs; ++i) {
    MPI_Status stat;
    MPI_Probe(MPI_ANY_SOURCE, 999, libcomm, &stat);
    int count; MPI_Get_count(&stat, recvtype, &count);
    lrranks->push_back(stat.MPI_SOURCE);
    lrdispls->push_back(cursize);
    lrsizes->push_back(count);
    (*rrankcount)++;
    cursize =  lrdispls->back() + count;
    rbuf = (char*)realloc(rbuf, cursize);
    reqs.resize(reqs.size()+1);

    //printf("[%i] receiving %i elements from %i\n", r, count, stat.MPI_SOURCE);
    MPI_Irecv(rbuf + lrdispls->back()*rcvext, count, recvtype, stat.MPI_SOURCE, 999, libcomm, &reqs.back());
  }

  MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUSES_IGNORE);

  *recvbuf = rbuf;
  *rranks = &(*lrranks)[0];
  *rdispls = &(*lrdispls)[0];
  *recvcounts = &(*lrsizes)[0];

  return MPI_SUCCESS;
}
