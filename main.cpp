/*
 * main.cpp
 *
 *  Created on: Dec 8, 2011
 *      Author: htor
 */
#include "dsde.h"
#include <vector>

int main() {

  MPI_Init(NULL, NULL);

  DSDE_Handle handle;
  std::vector<int> sbuf, sranks;
  std::vector<MPI_Aint> scounts, sdispls;
  int *rbuf, rrankcount, *rranks;
  MPI_Aint *rdispls, *rcounts;

  int p; MPI_Comm_size(MPI_COMM_WORLD, &p);
  int r; MPI_Comm_rank(MPI_COMM_WORLD, &r);
  if(p < 4) {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // initialize send structure
  sbuf.push_back(r);
  sbuf.push_back(r);
  sbuf.push_back(r);
  scounts.push_back(1);
  scounts.push_back(1);
  scounts.push_back(1);
  sdispls.push_back(0);
  sdispls.push_back(1);
  sdispls.push_back(2);

  sranks.push_back((r+1)%p);
  sranks.push_back((r+2)%p);
  sranks.push_back((r+3)%p);

  DSDE_Exchange_ibarrier(&sbuf[0], sranks.size(), &sranks[0], &scounts[0], &sdispls[0], MPI_INT,
                         (void**)&rbuf, &rrankcount, &rranks, &rcounts, &rdispls, MPI_INT, MPI_COMM_WORLD, &handle);

  if(2==2) {
    printf("[%i] received from %i ranks\n", r, rrankcount);
    for(int i=0; i<rrankcount; ++i) {
      printf("[%i] received %i elements at offset %i (%i) from %i\n", r, rcounts[i], rdispls[i], rbuf[i], rranks[i]);

    }
  }

  MPI_Finalize();
}

