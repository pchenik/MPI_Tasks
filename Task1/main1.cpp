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

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 
 int elements_per_proc = buffer_count / size;
 int *array_per_proc = new int[elements_per_proc];

  if (rank == 0) {
	  int* buffer = new int[buffer_count];
	  random_device rd;
          mt19937 gen(rd());
	  uniform_int_distribution<int> distribution(0, MAX_VALUE);
	  for (int i = 0; i < buffer_count; ++i) {
		buffer[i] = (17 * i + 7) % MAX_VALUE;
		//buffer[i] = distribution(gen);
		//cout << buffer[i] << " ";
	  }

	 t1 = MPI_Wtime();
	 for (int i = 0; i < elements_per_proc; ++i)
                  array_per_proc[i] = buffer[i];

          for (int i = 1; i < size; ++i)
                MPI_Send(buffer + i * elements_per_proc, elements_per_proc, MPI_INT, i, 0, MPI_COMM_WORLD);

	  //cout << endl;
	 delete [] buffer;
  }
  else {
	MPI_Recv(array_per_proc, elements_per_proc, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

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
  
  delete [] array_per_proc;
  
  MPI_Finalize();
  
  return 0;
}
