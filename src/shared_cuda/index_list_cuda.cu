#include "index_list_cuda.h"

#include <cuda.h>
#include <cub/device/device_select.cuh>
#include <cub/iterator/counting_input_iterator.cuh>

template<typename T>
struct ZeroCmp
{
	const T* conditions;
	const int startid;

	ZeroCmp(const int startid, const T* conditions) :
		startid(startid), conditions(conditions)
	{ }

	__device__ __host__ __forceinline__
	bool operator() (const int &id)
	{
	  return (conditions[ id - startid ] != 0);
	}
};

template <typename T>
static
void c_generate_index_list_cuda_generic_device(
			const T* dev_conditions,
			const int startid, const int endid,
			int* dev_indices,
			int* dev_nvalid, cudaStream_t stream)
{
	static size_t storageSize = 0;
	static char* storage = nullptr;

	const int n = endid - startid + 1;

	// Argument is the offset of the first element
	cub::CountingInputIterator<int> iterator(startid);

	// Determine temporary device storage requirements
	size_t storageRequirement;
	cub::DeviceSelect::Flagged(nullptr, storageRequirement,
			iterator, dev_conditions, dev_indices,
			dev_nvalid, n, stream);

	// Allocate temporary storage (only if not enough)
	if (storageRequirement > storageSize)
	{
		cudaFree(storage);
		cudaMalloc(&storage, storageRequirement);
		storageSize = storageRequirement;
	}

	ZeroCmp<T> select(startid, dev_conditions);
	cub::DeviceSelect::If(storage, storageRequirement,
			iterator, dev_indices,
			dev_nvalid, n,
			select, stream);
}


template <typename T>
static
void c_generate_index_list_cuda_batched_generic(
			const int batch_size,
			const T* dev_conditions, const int cond_stride,
			const int startid, const int endid,
			int* dev_indices, const int idx_stride,
			int* dev_nvalid, cudaStream_t stream)
{
	for (int i = 0; i < batch_size; i++)
		c_generate_index_list_cuda_generic_device(
				dev_conditions + cond_stride*i,
				startid, endid,
				dev_indices + idx_stride*i,
				dev_nvalid + i, stream);
}

template <typename T>
static
void c_generate_index_list_cuda_generic(
			const T* dev_conditions,
			const int startid, const int endid,
			int* dev_indices,
			int& nvalid, cudaStream_t stream)
{
	static int* dev_nvalid = nullptr;
	if (dev_nvalid == nullptr)
			cudaMalloc(&dev_nvalid, sizeof(int));

	c_generate_index_list_cuda_generic_device(
			dev_conditions, startid, endid, dev_indices, dev_nvalid, stream);

	cudaMemcpyAsync(&nvalid, dev_nvalid, sizeof(int), cudaMemcpyDeviceToHost, stream);
	cudaStreamSynchronize(stream);
}


void c_generate_index_list_cuda_i1(
			const char* dev_conditions,
			const int startid, const int endid,
			int* dev_indices,
			int& nvalid, cudaStream_t stream)
{
	c_generate_index_list_cuda_generic(dev_conditions, startid, endid, dev_indices, nvalid, stream);
}

void c_generate_index_list_cuda_i4(
			const int* dev_conditions,
			const int startid, const int endid,
			int* dev_indices,
			int& nvalid, cudaStream_t stream)
{
	c_generate_index_list_cuda_generic(dev_conditions, startid, endid, dev_indices, nvalid, stream);
}

void c_generate_index_list_cuda_batched_i4(
		const int batch_size,
		const int* dev_conditions, const int cond_stride,
		const int startid, const int endid,
		int* dev_indices, const int idx_stride,
		int* dev_nvalid, cudaStream_t stream)
{
	c_generate_index_list_cuda_batched_generic(
			batch_size,
			dev_conditions, cond_stride,
			startid, endid,
			dev_indices, idx_stride,
			dev_nvalid, stream);
}

