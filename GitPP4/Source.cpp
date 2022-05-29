#include <mpi.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>

constexpr int N = 6600;

#define MASTER_RANK 0

#define NEW_BLOCK 1
#define SORTED_BLOCK 2

#ifdef DEBUG_MODE
#define MSEND(y) std::cout << "Master(0): Sending block " << y << std::endl;
#define MRECV(y) std::cout << "Master(0): Received sorted block " << y << std::endl;
#define SSEND(x,y) std::cout << "Slave(" << x << "): Received block " << y << std::endl;
#define SRECV(x,y) std::cout << "Slave(" << x << "): Sending sorted block " << y << std::endl;
#else
#define MSEND(y)
#define MRECV(y)
#define SSEND(x,y)
#define SRECV(x,y)
#endif

void Master(int rank, int comm_size, int block_size);
void Slave(int rank, int comm_size, int block_size);
void GenerateRandomArray(int arr[]);
bool IsSorted(int arr[]);
void InsertionSort(int arr[], int size);
void InsertionSort(int arr_sorted[], int arr_unsorted[], int size1, int size2);

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	
	int rank, comm_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	int block_size = N / (comm_size-1);

	if (rank == MASTER_RANK)
		Master(rank, comm_size, block_size);
	else
		Slave(rank, comm_size, block_size);

	MPI_Finalize();
}

void Master(int rank, int comm_size, int block_size)
{
	int* arr = new int[N];
	GenerateRandomArray(arr);
	for (int block = 0; block < comm_size-1; ++block) {
		MSEND(block)
		MPI_Send(arr + block*block_size, block_size, MPI_INT, rank + 1, NEW_BLOCK, MPI_COMM_WORLD);
	}

	MPI_Status status;
	int last_block_size = N - block_size * (comm_size - 1);
	for (int block = 0; block < comm_size-1; ++block) {
		MPI_Probe(rank+block+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(arr + block*block_size, block_size, MPI_INT, rank+block+1, SORTED_BLOCK, MPI_COMM_WORLD, &status);
		InsertionSort(arr + block*block_size, arr + block_size*(comm_size-1), block_size, last_block_size);
		MRECV(block)
	}
	InsertionSort(arr + block_size*(comm_size-1), last_block_size);

	std::cout << IsSorted(arr) << std::endl;
}

void Slave(int rank, int comm_size, int block_size)
{
	MPI_Status status;
	int* arr_main = new int[block_size];
	MPI_Probe(rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	MPI_Recv(arr_main, block_size, MPI_INT, rank-1, NEW_BLOCK, MPI_COMM_WORLD, &status);
	SRECV(rank, 0)
	InsertionSort(arr_main, block_size);

	int* arr_recv_send = new int[block_size];
	for (int block = 1; block < comm_size - rank; ++block) {
		MPI_Probe(rank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(arr_recv_send, block_size, MPI_INT, rank-1, NEW_BLOCK, MPI_COMM_WORLD, &status);
		SRECV(rank, block)
		InsertionSort(arr_main, arr_recv_send, block_size, block_size);
		SSEND(rank, block)
		MPI_Send(arr_recv_send, block_size, MPI_INT, rank+1, NEW_BLOCK, MPI_COMM_WORLD);
	}
	SSEND(rank, rank-1)
	MPI_Send(arr_main, block_size, MPI_INT, MASTER_RANK, SORTED_BLOCK, MPI_COMM_WORLD);
}

void GenerateRandomArray(int arr[])
{
	srand(time(NULL));
	for(int i = 0; i < N; ++i)
		arr[i] = rand() % N;
}

bool IsSorted(int arr[])
{
	for (int i = 1; i < N; ++i)
		if (arr[i] < arr[i - 1])
			return false;
	return true;
}

void InsertElement(int arr[], int i, int j) {
	int element = arr[i];
	for (i; i > j; --i)
		arr[i] = arr[i-1];

	arr[j] = element;
}

void InsertElement(int arr_sorted[], int arr_unsorted[], int size_sorted, int i, int j) {
	int element = arr_unsorted[i];
	arr_unsorted[i] = arr_sorted[size_sorted-1];

	for (i = size_sorted-1; i > j; --i)
		arr_sorted[i] = arr_sorted[i-1];

	arr_sorted[j] = element;
}

void InsertionSort(int arr[], int size)
{
	for (int i = 1; i < size; ++i)
		for (int j = 0; j < i; ++j)
			if (arr[i] < arr[j])
				InsertElement(arr, i, j);
}

void InsertionSort(int arr_sorted[], int arr_unsorted[], int size_sorted, int size_unsorted)
{
	for (int i = 0; i < size_unsorted; ++i)
		for (int j = 0; j < size_sorted; ++j)
			if (arr_unsorted[i] < arr_sorted[j])
				InsertElement(arr_sorted, arr_unsorted, size_sorted, i, j);
}