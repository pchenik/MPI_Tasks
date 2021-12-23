#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include <random>
#include <fstream>

using namespace std;

int main(int argc, char **argv) {
  int buffer_sizes[] = {64, 128, 1024, 16384, 32768, 131072, 1048576, 16777216, 67108864};
  //int send_count[] = {10, 50, 100, 200, 500, 1000, 2000, 3000, 4000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000};
  int send_count[] = {10, 50, 100, 500, 1000, 2000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000}; 

  int size, rank;
  double t1, t2;
  constexpr int MAX_VALUE = 256;
  char* A = nullptr;
  //char* B = nullptr;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (int t = 0; t < 7; ++t) {
	  int N = buffer_sizes[t];
	  A = new char[N];
	  for (int i = 0; i < N; ++i)
		A[i] = (char)((2 * i + rank) % MAX_VALUE);
	  
	  MPI_Barrier(MPI_COMM_WORLD);
	  
	  for (int i = 0; i < 11; ++i) {
		  t1 = MPI_Wtime();
		  for (int j = 0; j < send_count[i]; ++j) {
			  if (rank == 0) {
				  MPI_Send(A, N, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
				  //for (int i = 0; i < N; ++i)
				//	  A[i] = char(A[i] + 1);
				  MPI_Recv(A, N, MPI_CHAR, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			  }
			  else {
				  MPI_Recv(A, N, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				  //for (int i = 0; i < N; ++i)
				//	  A[i] = char(A[i] % MAX_VALUE);
				  MPI_Send(A, N, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
			  }
		  }
		  if (rank == 0) {
			t2 = MPI_Wtime();
			cout << "Dispatch times = " << send_count[i] << "; N = " << N << "; clock time = " << t2 - t1 << endl;
		  }
	  }
	  
	  if (rank == 0)
	  	cout << t << endl;
	  delete [] A;
	  
	  MPI_Barrier(MPI_COMM_WORLD);
  }
  
  MPI_Finalize();
  return 0;
}

