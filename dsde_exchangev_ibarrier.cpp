/*
 * dsde_exchange_ibarrier.cpp
 *
 *  Created on: Dec 8, 2011
 *      Author: htor
 */
#include "dsde_internal.h"

#ifdef HAVE_DCMF
#include <dcmf_globalcollectives.h>
#else
#include <nbc.h>
#endif

// TODO: this is horrible and not thread-safe, this should be hung off the comm attribute!
static int nbctag=9;
static volatile char done=0;
static char initialized=0;

#ifdef HAVE_DCMF
static void cbfunc(void *clientdata, DCMF_Error_t *error) {
  done=1;
}

DCMF_Protocol_t barr_reg;
DCMF_Request_t barr_req;
#endif


static void prepare_mpi3nbc(MPI_Comm comm) {
#ifdef HAVE_DCMF
  int res;
  MPI_Comm_compare(comm, MPI_COMM_WORLD, &res);
  if(res == MPI_UNEQUAL) {
    printf("DCMF support only works on MPI_COMM_WORLD, bug the developer to extend this lib to work on non-world comms on BG :-)\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  DCMF_GlobalBarrier_Configuration_t barrier_config;
  barrier_config.protocol = DCMF_GI_GLOBALBARRIER_PROTOCOL;
  DCMF_GlobalBarrier_register(&barr_reg, &barrier_config);
#else
  // need to initialize communicator because LibNBC's first call is blocking!
  NBC_Handle hndl;
  NBC_Ibarrier(comm, &hndl);
  NBC_Wait(&hndl);
#endif
}

/* here we get only the processes that we receive from so that we have
 * to probe/malloc/receive the messages */
int DSDE_Exchangev_ibarrier(
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

  if(!initialized) {
    prepare_mpi3nbc(libcomm);
  }
  nbctag++; // the safe tag!

  std::vector<MPI_Request> sreqs(srankcount);
  for(int i=0; i<srankcount; ++i) {
    char *sbuf = ((char *) sendbuf) + (sdispls[i]*sndext);
    //printf("[%i] sending %i elements to %i\n", r, sendcounts[i], sranks[i]);
    MPI_Issend(sbuf, sendcounts[i], sendtype, sranks[i], nbctag, libcomm, &sreqs[i]);
  }

  std::vector<int> *lrranks = new std::vector<int>; // local rranks
  std::vector<MPI_Aint> *lrdispls = new std::vector<MPI_Aint>; // local rdispls
  std::vector<MPI_Aint> *lrsizes = new std::vector<MPI_Aint>; // local rsizes
  int cursize = 0; // size of receive buffer in elements
  char *rbuf=NULL; // local receive buffer pointer
  *rrankcount = 0; // input argument - number of receive ranks
  std::vector<MPI_Request> rreqs;

  int barr_act=0;
  done=0;
#ifndef HAVE_DCMF
  NBC_Handle hndl;
#endif
  while(!done) {
    int flag;
    MPI_Status stat;
    MPI_Iprobe(MPI_ANY_SOURCE, nbctag, libcomm, &flag, &stat);
    if(flag) {
      int count;
      MPI_Get_count(&stat, recvtype, &count);
      lrranks->push_back(stat.MPI_SOURCE);
      lrdispls->push_back(cursize);
      lrsizes->push_back(count);
      (*rrankcount)++;
      cursize =  lrdispls->back() + count;
      rbuf = (char*)realloc(rbuf, cursize);
      rreqs.resize(rreqs.size()+1);

      printf("[%i] receiving %i elements from %i\n", r, count, stat.MPI_SOURCE);
      MPI_Irecv(rbuf + lrdispls->back()*rcvext, count, recvtype, stat.MPI_SOURCE, nbctag, libcomm, &rreqs.back());
    }

    if(barr_act) {
#ifndef HAVE_DCMF
      if(NBC_OK == NBC_Test(&hndl)) {
        NBC_Wait(&hndl); // needed to free request
        //if(!r) std::cout <<"["<<r<<"] barrier finished\n";
        done=1;
      }
#endif
    } else {
      MPI_Testall(sreqs.size(), &sreqs[0], &flag, MPI_STATUSES_IGNORE);
      if(flag || sreqs.empty()) {
#ifndef HAVE_DCMF
        NBC_Ibarrier(libcomm, &hndl);
#else
        DCMF_Callback_t callback={ cbfunc, (void*)NULL };
        DCMF_GlobalBarrier(&barr_reg, &barr_req, callback);
#endif
        //if(!r) std::cout <<"["<<r<<"] starting barrier\n";
        barr_act=1;
      }
    }
  }
  MPI_Waitall(rreqs.size(), &rreqs[0], MPI_STATUSES_IGNORE);

  *recvbuf = rbuf;
  *rranks = &(*lrranks)[0];
  *rdispls = &(*lrdispls)[0];
  *recvcounts = &(*lrsizes)[0];

  return MPI_SUCCESS;
}
