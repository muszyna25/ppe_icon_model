#pragma once

#ifdef __cplusplus // Are we compiling this with a C++ compiler ?
extern "C"
{
#endif

	void c_generate_index_list_cuda_i1(
			const char* dev_conditions,
			const int startid, const int endid,
			int* dev_indices, int& nvalid, cudaStream_t stream);

#ifdef __cplusplus
}
#endif
