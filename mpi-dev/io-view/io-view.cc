#include <iostream>
#include "mpi.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int num_of_process, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &num_of_process);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int N_per_t = 10;
  int N_t = 5;
  
  int num_of_elements_per_proc, num_of_elements_of_me;
  num_of_elements_per_proc = (N_per_t + (num_of_process - 1)) / num_of_process;
  num_of_elements_of_me = num_of_elements_per_proc;
  if (rank == num_of_process - 1) {
    num_of_elements_of_me = N_per_t - (num_of_process-1) * num_of_elements_per_proc;
  }
  fprintf(stdout, "[@rank=%d] num_of_elements_of_me = %d\n", rank, num_of_elements_of_me);
  
  MPI_Datatype element_type, blocK_type;
  int element_type_size;
  element_type = MPI_DOUBLE;
  MPI_Type_size(element_type, &element_type_size);
  MPI_Type_vector(N_t, num_of_elements_of_me, N_per_t, element_type, &blocK_type);
  MPI_Type_commit(&blocK_type);

  MPI_Offset disp;
  disp = rank * num_of_elements_per_proc * element_type_size;
  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD, "arr.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, disp, element_type, blocK_type, "native", MPI_INFO_NULL);

  double buf[num_of_elements_of_me*N_t];
  int t_i, i, offset_per_t, total_index, buf_index;
  for (t_i=0; t_i<N_t; t_i++) {
    offset_per_t = rank * num_of_elements_per_proc;
    for (i=offset_per_t; i<offset_per_t+num_of_elements_of_me; i++) {
      buf_index = t_i * num_of_elements_of_me + (i-offset_per_t);
      total_index = t_i * N_per_t + i;
      buf[buf_index] = 0.1 * total_index;
      fprintf(stdout, "[@rank=%d] buf[%d] = %f\n", rank, buf_index, buf[buf_index]);
    }
  }
  MPI_Status write_status;
  MPI_File_write(fh, buf, num_of_elements_of_me*N_t, MPI_DOUBLE, &write_status);

  MPI_File_close(&fh);

  MPI_Finalize();
}
