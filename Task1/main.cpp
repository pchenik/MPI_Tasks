#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include <random>
#include <sstream>

using namespace std;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int size, rank;
  double t1, t2;
  int buffer_count = (1 << 28);
  int MAX_VALUE = 1e9;
  int* buffer = nullptr;
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
	  //t1 = MPI_Wtime();
	  //cout << buffer_count << endl;
	  buffer = new int[buffer_count];
	  random_device rd;
          mt19937 gen(rd());
	  uniform_int_distribution<int> distribution(0, MAX_VALUE);
	  for (int i = 0; i < buffer_count; ++i) {
		buffer[i] = (17 * i + 7) % MAX_VALUE;
		//buffer[i] = distribution(gen);
		//cout << buffer[i] << " ";
	  }
	  t1 = MPI_Wtime();
	  //cout << endl;
  }
  //MPI_Barrier(MPI_COMM_WORLD);
  
  int elements_per_proc = buffer_count / size;
  int *array_per_proc = new int[elements_per_proc];
  MPI_Scatter(buffer, elements_per_proc, MPI_INT, array_per_proc, elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
  /*for (int i = 0; i < elements_per_proc; ++i)
	cout << array_per_proc[i] << " ";
  cout << endl;*/

  int min_value = MAX_VALUE;
  stringstream mystring;
  for(int i = 0; i < elements_per_proc; ++i) {
	  min_value = min(min_value, array_per_proc[i]);
	  mystring << min_value << " ";
  }

  int total_min = MAX_VALUE;
  MPI_Reduce(&min_value, &total_min, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

  if (rank == 0) {
	t2 = MPI_Wtime();
	cout << "min = " << total_min << "; clock time = " << t2 - t1 << endl;
  }
  
  delete [] buffer;
  delete [] array_per_proc;
  
  MPI_Finalize();
  
  return 0;
}
