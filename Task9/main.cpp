#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include <random>
#include <fstream>
#include <sstream>

#define sqr(x) ((x) * (x))

using namespace std;

void compute(int total_count, int my_count, int* my_points[3], stringstream& ss, double& sum_dist) {
  double local_sum[3] = {0.0, 0.0, 0.0};
  for (int j = 0; j < 3; ++j)
    for (int i = 0; i < my_count; ++i)
        local_sum[j] += my_points[j][i];
  //cout << "IN COMPUTE " << total_count << endl; 
  double barycentre[3];
  MPI_Allreduce(local_sum, barycentre, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  for (int j = 0; j < 3; ++j)
    barycentre[j] /= total_count;
  
  for (int i = 0; i < my_count; ++i) {
    double dist = 0.0;
    for (int j = 0; j < 3; ++j)
	dist += (my_points[j][i] - barycentre[j]) * (my_points[j][i] - barycentre[j]);
    dist = sqrt(dist);
    sum_dist += dist;
    //cout << dist << endl;
    ss << dist << endl;
  }
  
}

int main(int argc, char **argv) {
  int buffer_sizes[] = {64, 128, 1024, 16384, 32768, 131072, 1048576, 16777216, 67108864};
  int size, rank;
  double t1, t2;
  constexpr int MAX_VALUE = 1e5;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  int* A[3] = {nullptr, nullptr, nullptr};
  int* A_local[3] = {nullptr, nullptr, nullptr};
  for (int t = 0; t < 9; ++t) {
	  int total_count = buffer_sizes[t];
	  int elements_per_proc = total_count / size;
	  
	  for(int k = 0; k < 3; ++k)
			A_local[k] = new int[elements_per_proc];
	  
	  if (rank == 0) {
		  for(int k = 0; k < 3; ++k)
			A[k] = new int[total_count];
		
		  for(int k = 0; k < 3; ++k)
			for (int i = 0; i < total_count; ++i)
				A[k][i] = (2 * i + k) % MAX_VALUE;
			
		  t1 = MPI_Wtime();
	  }
	  
	  MPI_Barrier(MPI_COMM_WORLD);
	  for (int k = 0; k < 3; ++k)
		MPI_Scatter(A[k], elements_per_proc, MPI_INT, A_local[k], elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
	  
	  stringstream ss;
	  double local_sum = 0, total_sum = 0;

          compute(total_count, elements_per_proc, A_local, ss, local_sum);
	  //MPI_Reduce(&local_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	  if (rank == 0) {
		t2 = MPI_Wtime();
		//cout << total_count << "; clock time = " << t2 - t1 << endl;
		cout << t2 - t1 << endl;
		for(int k = 0; k < 3; ++k)
                        delete [] A[k];
	  }
	  
	  for(int k = 0; k < 3; ++k) {
			delete [] A_local[k];
	  }
	  
	  MPI_Barrier(MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
	cout << endl << endl;
  }

  MPI_Finalize();
  return 0;
}
