#include "index_list_cuda.h"

#include <cuda.h>
#include <cub/device/device_select.cuh>
#include <cub/iterator/counting_input_iterator.cuh>

void c_generate_index_list_cuda_i1(
			const char* dev_conditions,
			const int startid, const int endid,
			int* dev_indices,
			int& nvalid, cudaStream_t stream)
{
	static size_t storageSize = 0;
	static char* storage = nullptr;
	static int*  dev_nvalid = nullptr;

	const int n = endid - startid + 1;

	if (dev_nvalid == nullptr)
		cudaMalloc(&dev_nvalid, sizeof(int));

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

	cub::DeviceSelect::Flagged(storage, storageRequirement,
			iterator, dev_conditions, dev_indices,
			dev_nvalid, n, stream);

	cudaMemcpyAsync(&nvalid, dev_nvalid, sizeof(int), cudaMemcpyDeviceToHost, stream);
	cudaStreamSynchronize(stream);
}



void c_generate_index_list_cuda_batched_i1(
			const char* dev_conditions,
			int nbatches,
			const int startid, const int endid,
			int* dev_indices,
			int& nvalid, cudaStream_t stream)
{
	return;
}

