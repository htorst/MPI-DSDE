#include <dsde.h>

#include <vector>
#include <stdlib.h>
#include <assert.h>

int DSDE_Exchange(void* sendbuf, int *sendcounts, int *sdispls,
                  MPI_Datatype sendtype, int destcount, int *dests,
                  void** recvbuf, int **recvcounts, MPI_Datatype recvtype, 
                  int *srccount, int **sources, MPI_Comm comm) {

  int p;
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  std::vector<int> ssizes(p,0), rsizes(p);
  for(int i=0; i < destcount; i++) ssizes[ dests[i] ] = sendcounts[i];
  
  MPI_Alltoall(&ssizes[0], 1, MPI_INT, &rsizes[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<MPI_Request> reqs(destcount);

  MPI_Aint lb, extent;
  MPI_Type_get_extent(sendtype, &lb, &extent);

  for(int i=0; i < destcount; ++i) {
    MPI_Isend((char*)sendbuf+extent*sdispls[i], sendcounts[i], sendtype, dests[i], 999, comm, &reqs[i]);
  }

  assert(sources == NULL); // avoid memory leaks
  assert(recvbuf == NULL); // avoid memory leaks
  *srccount = 0;

  int rbufsize = 0;
  for(int i=0; i<p; i++) if(rsizes[i]) {
    (*srccount)++;
    rbufsize+=rsizes[i];
  }

  *recvbuf = malloc(rbufsize*extent);

  int reqoffset = reqs.size()-1;
  reqs.resize(reqs.size()+*srccount);
  *sources = (int*)malloc(sizeof(int)* *srccount);
  *recvcounts = (int*)malloc(sizeof(int)* *srccount);
  
  int offset = 0;
  int index = 0;
  for(int i=0; i<p; i++) if(rsizes[i]) {
    MPI_Irecv((char*)*recvbuf + extent*offset, rsizes[i], recvtype, i, 999, comm, &reqs[reqoffset++]);
    offset += rsizes[i];
    (*sources)[index] = i;
    (*recvcounts)[index] = rsizes[i];
    index++;
  }

  MPI_Waitall(reqs.size(), &reqs[0], MPI_STATUSES_IGNORE);

  return MPI_SUCCESS;

}

