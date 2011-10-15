#include <mpi.h>

extern "C" {

int DSDE_Exchange(void* sendbuf, // the sendbuffer -- provided by user
                  int *sendcounts,  // number of elements to send to dests[i]
                  int *sdispls, // displs to send to dests[i]
                  MPI_Datatype sendtype, // datatype
                  int destcount, // number of dest. processes
                  int *dests, // array of dest. processes
                  void** recvbuf, // receive buffer, allocated by lib
                  int **recvcounts, // recvcounts, allocated by lib
                  MPI_Datatype recvtype, // datatype
                  int *srccount, // number of source processes
                  int **sources, // array of source processes
                  MPI_Comm comm); // communicator

}
