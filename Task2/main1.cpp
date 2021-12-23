#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include <sstream>

using namespace std;

int main(int argc, char **argv) {
  int size, rank;
  double t1, t2;
  int buffer_count = (1 << 28);
  int MAX_VALUE = 1e2;
  
  MPI_Init(&argc, &argv); 
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Request request;

  int elements_per_proc = buffer_count / size;
  
  int *A_per_proc = new int[elements_per_proc];

  if (rank == 0) { 
	  int* A = new int[buffer_count];
	  for (int i = 0; i < buffer_count; ++i) {
			 A[i] = (2 * i + 3) % MAX_VALUE;
			//cout << A[i] << " " << endl;
	  }
	  //cout << endl;
	  t1 = MPI_Wtime();
	  for (int i = 0; i < elements_per_proc; ++i)
		  A_per_proc[i] = A[i];

	  for (int i = 1; i < size; ++i)
		MPI_Send(A + i * elements_per_proc, elements_per_proc, MPI_INT, i, 0, MPI_COMM_WORLD);
		//MPI_Isend(A + i * elements_per_proc, elements_per_proc, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
	  
	  //MPI_Wait(&request, MPI_STATUS_IGNORE);
	  delete [] A;
  }
  else {
	  MPI_Recv(A_per_proc, elements_per_proc, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
 
  long long scal_product_local = 0;
  stringstream mystring;
  for(int i = 0; i < elements_per_proc; ++i) {
          scal_product_local += (long long)A_per_proc[i] * A_per_proc[i];
          mystring << scal_product_local << " ";
  }

  long long scal_product_total = 0;
  MPI_Reduce(&scal_product_local, &scal_product_total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
	t2 = MPI_Wtime();
	cout << " = " << scal_product_total << "; clock time = " << t2 - t1 << endl;
  }
  
  delete [] A_per_proc;
  
  MPI_Finalize();
  
  return 0;
}
