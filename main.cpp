/*
 * Copyright (c) 2012 The Trustees of University of Illinois.
 *                    All rights reserved.
 *
 * Author(s): Torsten Hoefler <htor@illinois.edu>
 *
 */

/*
 * main.cpp
 *
 *  Created on: Dec 8, 2011
 *      Author: htor
 */
#include "dsde.h"
#include <vector>

// real apps should just include dsde.h
#include "dsde_internal.h"

int main() {
  int i;

  MPI_Init(NULL, NULL);

  DSDE_Handle handle;
  std::vector<int> sbuf, sranks;
  std::vector<MPI_Aint> scounts, sdispls;
  int *rbuf, rrankcount, *rranks;
  MPI_Aint *rdispls, *rcounts;

  int p; MPI_Comm_size(MPI_COMM_WORLD, &p);
  int r; MPI_Comm_rank(MPI_COMM_WORLD, &r);
  if(p < 4) {
//    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // initialize send structure
  sbuf.push_back(r+100);
  sbuf.push_back(r+1000);
  sbuf.push_back(r+100);
  sbuf.push_back(r+1000);
  sbuf.push_back(r+100);
  sbuf.push_back(r+1000);
  scounts.push_back(2);
  scounts.push_back(2);
  scounts.push_back(2);
  sdispls.push_back(0);
  sdispls.push_back(2);
  sdispls.push_back(4);

  sranks.push_back((r-1+p)%p);
  sranks.push_back((r+2+p)%p);
//  sranks.push_back((r-3+p)%p);

  int degree = 2;

//  DSDE_Exchangev_brucks(&sbuf[0], sranks.size(), &sranks[0], &scounts[0], &sdispls[0], MPI_INT,
//                         (void**)&rbuf, &rrankcount, &rranks, &rcounts, &rdispls, MPI_INT, MPI_COMM_WORLD, &handle, degree);
  DSDE_Exchangev_accumulate(&sbuf[0], sranks.size(), &sranks[0], &scounts[0], &sdispls[0], MPI_INT,
                         (void**)&rbuf, &rrankcount, &rranks, &rcounts, &rdispls, MPI_INT, MPI_COMM_WORLD, &handle);

  if(2==2) {
    printf("[%i] received from %i ranks\n", r, rrankcount);
    for(int i=0; i<rrankcount; ++i) {
    for (int j=0; j < rcounts[i]; j++) {
      int index = rdispls[i] + j;
      printf("[%i] received %i elements starting at offset %i, offset %d = (%i) from %i\n", r, rcounts[i], rdispls[i], index, rbuf[index], rranks[i]);
    }
    }
  }

  DSDE_Free(&handle);

  /* execute a sparse reduce scatter operation */
  int result_flag;
  std::vector<int> result(scounts[0]);
  DSDE_Reduce_scatter_block_hbrucks(
    &sbuf[0], sranks.size(), &sranks[0], &sdispls[0],
    &result_flag, &result[0],
    scounts[0], MPI_INT, MPI_SUM, MPI_COMM_WORLD, degree
  );

  if (result_flag) {
    for (i=0; i < scounts[0]; i++) {
      printf("[%i] result[%d]=%d\n", r, i, result[i]);
    }
  } else {
    printf("[%i] no result\n", r);
  }

  MPI_Finalize();
}

