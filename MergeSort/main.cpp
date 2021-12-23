#include <iostream>
#include <random>
#include <algorithm>
#include <mpi.h>
#include <cmath>

using namespace std;

void merge_sort(int *a, int n) {
	if (n <= 1)
		return;
	merge_sort(a, n / 2);
	merge_sort(a + n / 2, n / 2);
	int* res = new int[n];
	merge(a, a + n / 2, a + n / 2, a + n, res);
	for (int i = 0; i < n; ++i)
		a[i] = res[i];
	delete [] res;
}

int* parallel_merge_sort(int rank, int* local_array, int* main_array, int local_N, int depth, MPI_Comm comm) {
	int right_rank;
	int *left, *right, *res;
	merge_sort(local_array, local_N);
	left = local_array;

	//Выполнение слияние отсортированных массивов поднимаясь по двоичному дереву
	for (int cur_depth = 0; cur_depth < depth; ) {
		//Проверка на то, что данный процесс поднялся на уровень вверх по левому ребру
		if ((rank & (~(1ll << cur_depth))) == rank) {
			//Номер процесса, который отвечает за правого потомка
			right_rank = (rank | (1ll << cur_depth));
			right = new int[local_N];
			//Ожидание прихода сообщения с отсортированным массивом от правого потомка
			MPI_Recv(right, local_N, MPI_INT, right_rank, 0,
				MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			res = new int[local_N * 2];
			merge(left, left + local_N, right, right + local_N, res);
			left = res;
			res = nullptr;
			delete [] right;
			local_N <<= 1;
			cur_depth++;
		}
		else {
			MPI_Send(left, local_N, MPI_INT, rank & (~(1ll << cur_depth)), 0, MPI_COMM_WORLD);
              		if (cur_depth) delete [] left;
              		cur_depth = depth;
		}
	}

	if (rank) main_array = left;
	return left;
}

int main(int argc, char** argv) {
	int sizes[] = {64, 128, 1024, 16384, 32768, 131072, 1048576, 16777216, 67108864};
	int N, rank, value, size, depth, local_N;
	int *main_array, *local_array, *array;
	double local_time, start_time, total_time;

        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

	for (int t = 0; t < 9; ++t) {
	N = sizes[t];
	if (rank == 0) {
		main_array = new int[N];
		array = new int[N];
		
		//Генерация массивов
		default_random_engine gen;
        	uniform_int_distribution<int> distribution(0, N);
        	for (int i = 0; i < N; ++i)
                	array[i] = main_array[i] = distribution(gen);

		//Последовательная сортировка
		start_time = MPI_Wtime();
		merge_sort(array, N);
		cout << "For consequtive variant and N = " << N << " time of execution equals to " << 
			MPI_Wtime() - start_time << endl;

		/*for (int i = 0; i < N; ++i)
			cout << main_array[i] << " ";
		cout << endl;*/
	}
	
	//Формирование блоков массивов для каждого процесса и перессылка данных на каждый из них
	local_N = N / size;
	local_array = new int[local_N];
	depth = log2(size);
	MPI_Scatter(main_array, local_N, MPI_INT, local_array, local_N, MPI_INT, 0, MPI_COMM_WORLD);

	start_time = MPI_Wtime();
	if (rank == 0) {
		/*for (int i = 0; i < local_N; ++i)
                        cout << local_array[i] << " ";
                cout << endl;*/
		main_array = parallel_merge_sort(rank, local_array, main_array, local_N, depth, MPI_COMM_WORLD);
	}
	else {
		parallel_merge_sort(rank, local_array, main_array, local_N, depth, MPI_COMM_WORLD);
	}

	//Подсчет времени для каждого процесса и вычисление времени выполнения самого долго из них
	//с помощью редукции
	local_time = MPI_Wtime() - start_time;
	MPI_Reduce(&local_time, &total_time, 1, MPI_DOUBLE,
        	MPI_MAX, 0, MPI_COMM_WORLD);

	if (rank == 0) {

		cout << "For parallel variant and N = " << N << " time of execution equals to " 
			<< total_time  << endl << endl;

		//Проверка на правильность сортировки
		for (int i = 0; i < N - 1; ++i)
			if (main_array[i] > main_array[i + 1] || array[i] > array[i + 1]) {
				cout << "Array is not sorted" << endl;
				break;
			}

	}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}
