/*
 * Copyright (c) 2012 The Trustees of University of Illinois.
 *                    All rights reserved.
 *
 * Author(s): Torsten Hoefler <htor@illinois.edu>
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <limits>
#include <iostream>
#include <ctime>
#include <set>
#include <queue>

#include "mpi.h"
#include "dsde.h"

extern "C" {

/* -----------------------------------------
 * DSDE_Free implementations
 * ----------------------------------------- */

/* function pointer to a DSDE_Free implementation that takes a pointer to a handle */
typedef int(*DSDE_Free_fn)(DSDE_Handle*);

/* assumes that handle just points to one big block of memory that must be freed */
int DSDE_Free_single(DSDE_Handle* handle);

/* -----------------------------------------
 * DSDE_Exchange implementations
 * ----------------------------------------- */

int DSDE_Exchange_brucks_inline(
    const void*  sendbuf, int  srankcount, const int  sranks[], MPI_Aint sendcount, const MPI_Aint  sdispls[], MPI_Datatype sendtype,
    void**       recvbuf, int* rrankcount, int*       rranks[], MPI_Aint recvcount, MPI_Aint*       rdispls[], MPI_Datatype recvtype,
    MPI_Comm comm, DSDE_Handle* handle,
    int degree);

/* -----------------------------------------
 * DSDE_Exchangev implementations
 * ----------------------------------------- */

int DSDE_Exchangev_alltoall(
  const void*  sendbuf, int  srankcount, const int  sranks[], const MPI_Aint  sendcounts[], const MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void**       recvbuf, int* rrankcount, int*       rranks[],       MPI_Aint* recvcounts[],       MPI_Aint* rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle);
int DSDE_Exchangev_reduce_scatter(
  const void*  sendbuf, int  srankcount, const int  sranks[], const MPI_Aint  sendcounts[], const MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void**       recvbuf, int* rrankcount, int*       rranks[], MPI_Aint*       recvcounts[], MPI_Aint*       rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle);
int DSDE_Exchangev_accumulate(
  const void*  sendbuf, int  srankcount, const int  sranks[], const MPI_Aint  sendcounts[], const MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void**       recvbuf, int* rrankcount, int*       rranks[], MPI_Aint*       recvcounts[], MPI_Aint*       rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle);
int DSDE_Exchangev_ibarrier(
  const void*  sendbuf, int  srankcount, const int  sranks[], const MPI_Aint  sendcounts[], const MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void**       recvbuf, int* rrankcount, int*       rranks[], MPI_Aint*       recvcounts[], MPI_Aint*       rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle);
int DSDE_Exchangev_brucks(
  const void*  sendbuf, int  srankcount, const int  sranks[], const MPI_Aint  sendcounts[], const MPI_Aint  sdispls[], MPI_Datatype sendtype,
  void**       recvbuf, int* rrankcount, int*       rranks[], MPI_Aint*       recvcounts[], MPI_Aint*       rdispls[], MPI_Datatype recvtype,
  MPI_Comm comm, DSDE_Handle* handle,
  int degree);

/* -----------------------------------------
 * DSDE_Reduce_scatter_block implementations
 * ----------------------------------------- */

int DSDE_Reduce_scatter_block_brucks(
    const void* sendbuf, int srankcount, const int sranks[], const MPI_Aint sdispls[],
    int* flag, void* recvbuf, MPI_Aint count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
    int degree);

/* avoids creating some datatypes at the cost of being restricted to homogeneous systems */
int DSDE_Reduce_scatter_block_hbrucks(
    const void* sendbuf, int srankcount, const int sranks[], const MPI_Aint sdispls[],
    int* flag, void* recvbuf, MPI_Aint count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
    int degree);

} /* extern "C" */
