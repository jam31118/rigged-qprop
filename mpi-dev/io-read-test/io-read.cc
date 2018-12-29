#include <iostream>
#include <cstdlib>
#include "mpi.h"

using namespace std;

int error_and_exit(int rank, int return_code, const char *func_name) {
  fprintf(stderr, "In process with rank '%d': something got wrong in '%s' with return code: '%d'\n",
      rank, func_name, return_code);
  MPI_Finalize();
  return return_code;
}

int main(int argc, char *argv[]) {

  //// Declare variables with optional initialization
  int return_code = -1;
  int num_of_process = -1, rank = -1;

  //// MPI Initialization
  return_code = MPI_Init(&argc, &argv);
  if (return_code == MPI_SUCCESS) {
    fprintf(stdout, "[ LOG ] MPI program has been initialized.\n");
  } else {
    fprintf(stderr, "[ERROR] Something got wrong during initialization.\n");
    return 1;
  }

  //// Configuration
  // data file path
  const char *file_name = argv[1];
  if (argc >= 2) {
    fprintf(stdout, "file_name: '%s'\n", argv[1]);
  } else {
    fprintf(stderr, "[ERROR] No file name given\n");
    int arg_idx;
    for (arg_idx=0; arg_idx<argc; arg_idx++) {
      fprintf(stdout, "[ LOG ] agrv[%d] = '%s'\n", arg_idx, argv[arg_idx]);
    }
    return 1;
  }

  // set data type
  MPI_Datatype datatype = MPI_DOUBLE;


  //// Getting MPI communication size and process rank
  return_code = MPI_Comm_size(MPI_COMM_WORLD, &num_of_process);
  if (return_code != MPI_SUCCESS || num_of_process < 0) {
    fprintf(stderr, "[ERROR} Something got wrong during getting number of processes\n");
  }
  return_code = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (return_code != MPI_SUCCESS || rank < 0) {
    fprintf(stderr, "[ERROR} Something got wrong during getting rank\n");
  }


  //// Read data from the common file
  fprintf(stdout, "So, I am process with rank %d out of size %d\n", rank, num_of_process);

  // Open file
  MPI_File fh;
  return_code = MPI_File_open(MPI_COMM_WORLD, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_open"); }

  // Get file size
  MPI_Offset file_size;
  return_code = MPI_File_get_size(fh, &file_size);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_get_size"); }
  if (rank == 0) { fprintf(stdout, "[@rank=%d] The size of file '%s' = %lld\n", rank, file_name, file_size); }

  // Construct datatype's size
  int datatype_size;
  return_code = MPI_Type_size(datatype, &datatype_size);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_Type_size"); }

  // Determine number of elements to read for each process
  long num_of_elements_in_file = file_size / datatype_size;
  if (file_size % datatype_size != 0) { return error_and_exit(rank, 1, "The size of given file is not a multiple of elementary datatype's size."); }
  long num_of_elements_per_proc;
  num_of_elements_per_proc = ( num_of_elements_in_file + (num_of_process - 1) ) / num_of_process;
  long num_of_elements_to_read;
  if ( num_of_elements_in_file >= num_of_process ) {
    num_of_elements_to_read = num_of_elements_per_proc;
    if (rank == num_of_process - 1) { 
      num_of_elements_to_read = num_of_elements_in_file - (num_of_process - 1) * num_of_elements_per_proc; 
    }
  } else {
    if ( rank < num_of_elements_in_file ) { num_of_elements_to_read = 1; }
    else { num_of_elements_to_read = 0; } 
  }

  // Determine offset for each process
  MPI_Offset read_offset;
  read_offset = rank * num_of_elements_per_proc * datatype_size;
  
  // Memory allocation
  double *p_buf = (double *) malloc(num_of_elements_to_read * datatype_size);
  if (p_buf == NULL) { return error_and_exit(rank, 1, "malloc"); }

  // Read data collectively from a single data file
  MPI_Status read_status;
  return_code = MPI_File_read_at(fh, read_offset, p_buf, num_of_elements_to_read, datatype, &read_status);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_read_at"); }

  // Print data
  long i;
  for (i=0; i<num_of_elements_to_read; i++) { fprintf(stdout, "[@rank=%d] buf[%ld] = %f\n", rank, i, p_buf[i]); }

  // Close file
  return_code = MPI_File_close(&fh);
  if (return_code != MPI_SUCCESS) { return error_and_exit(rank, return_code, "MPI_File_close"); }


  //// MPI Finalization
  return_code = MPI_Finalize();
  if (return_code == MPI_SUCCESS) {
    fprintf(stdout, "[ LOG ] MPI program has been finalized.\n");
  } else {
    fprintf(stderr, "[ERROR] Something got wrong during finalization.\n");
    return -1;
  }

  // End of the program
  return 0;
}

