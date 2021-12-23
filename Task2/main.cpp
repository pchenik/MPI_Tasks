#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include <random>
#include <fstream>
#include <sstream>

using namespace std;

int main(int argc, char **argv) {
  int size, rank;
  double t1, t2;
  constexpr int buffer_count = (1 << 28);
  constexpr int MAX_VALUE = 1e2;
  int* A = nullptr;
  //int* B = nullptr;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int elements_per_proc = buffer_count / size;
  int *A_per_proc = new int[elements_per_proc];
  //int *B_per_proc = new int[elements_per_proc];

  if (rank == 0) {
          //t1 = MPI_Wtime();
          A = new int[buffer_count];
          //B = new int[buffer_count];
          for (int i = 0; i < buffer_count; ++i) {
                        A[i] = (2 * i + 3) % MAX_VALUE;
                        //B[i] = (3 * i + 1) % MAX_VALUE;
                        //cout << A[i] << " " << B[i] << endl;
          }
          //cout << endl;
          t1 = MPI_Wtime();
  }

  MPI_Scatter(A, elements_per_proc, MPI_INT, A_per_proc, elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
  //MPI_Scatter(B, elements_per_proc, MPI_INT, B_per_proc, elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);

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

  delete [] A;
  //delete [] B;

  delete [] A_per_proc;
  //delete [] B_per_proc;

  MPI_Finalize();
  return 0;
}

