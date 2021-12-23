#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include <random>
#include <fstream>
#include <sstream>

#define sqr(x) ((x) * (x))

using namespace std;

int main(int argc, char **argv) {
  int buffer_sizes[] = {64, 128, 1024, 16384, 32768, 131072, 1048576, 16777216, 67108864 / 2};
  int size, rank;
  double t1, t2;
  constexpr int MAX_VALUE = 1e9;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (int t = 0; t < 9; ++t) {
	  
	  int total_count = buffer_sizes[t];
	  long long* A = new long long[total_count];
	  long long* A_total = nullptr;
	  t1 = MPI_Wtime(); 

	  for(int i = 0; i < total_count; ++i)
		A[i] = (239 * i + rank) % MAX_VALUE;
	 
	  if (rank == 0)
		A_total = new long long[total_count];
	  MPI_Reduce(A, A_total, total_count, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	  
	  /*stringstream ss;
	  ofstream myfile("output_file.txt");
	  for(int i = 0; i < total_count; ++i) {
                 int k = 0;
                 for(int j = 2; j * j <= A[i]; ++j)
                        if (A[i] % j) ++k;
                 ss << k << " ";
		 myfile << k << " ";
          }
	  myfile.close();*/

	  //cout << "RANK IS " << rank << " SIZE IS " << size << endl;
	  if (rank == 0) {
		t2 = MPI_Wtime();
		//cout << A_total[1] << " " << total_count << "; clock time = " << t2 - t1 << endl;
		cout << t2 - t1 << endl;
		delete [] A_total;
	  }
	  
	  /*MPI_Barrier(MPI_COMM_WORLD);
	  for(int i = 0; i < 10; ++i)
		cout << A[i] << " ";
	  cout << endl << endl;*/

	  delete [] A;
	  
	  MPI_Barrier(MPI_COMM_WORLD);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
	  cout << endl;
  }
  
  MPI_Finalize();
  return 0;
}
