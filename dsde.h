#include <mpi.h>

extern "C" {

int DSDE_Exchange(void* sendbuf, int *sendcounts, int *sdispls,
                  MPI_Datatype sendtype, int destcount, int *dests,
                  void* recvbuf, int *recvcounts, int *rdispls,
                  MPI_Datatype recvtype, int srccount, int *sources,
                  MPI_Comm comm);

}
