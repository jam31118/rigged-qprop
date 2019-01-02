#include <iostream>
#include "mpi.h"

using namespace std;

int main(int argc, char *argv[]) {

  // Declare with optional initialization
  int return_code = -1;
  int num_of_process = -1, rank = -1;

  // MPI Initialization
  return_code = MPI_Init(&argc, &argv);
  if (return_code == MPI_SUCCESS) {
    fprintf(stdout, "[ LOG ] MPI program has been initialized.\n");
  } else {
    fprintf(stderr, "[ERROR] Seomthing got wrong during initialization.\n");
  }

  // Getting MPI communication size and process rank
  return_code = MPI_Comm_size(MPI_COMM_WORLD, &num_of_process);
  if (return_code != MPI_SUCCESS || num_of_process < 0) {
    fprintf(stderr, "[ERROR} Seomthing got wrong during getting number of processes\n");
  }
  return_code = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (return_code != MPI_SUCCESS || rank < 0) {
    fprintf(stderr, "[ERROR} Seomthing got wrong during getting rank\n");
  }

  // Main Body
  fprintf(stdout, "So, I am process with rank %d out of size %d\n", rank, num_of_process);

  // MPI Finalization
  return_code = MPI_Finalize();
  if (return_code == MPI_SUCCESS) {
    fprintf(stdout, "[ LOG ] MPI program has been finalized.\n");
  } else {
    fprintf(stderr, "[ERROR] Seomthing got wrong during finalization.\n");
  }

  // End of the program
  return 0;
}

