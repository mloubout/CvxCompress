#include <cuda_runtime_api.h>

//
// CUDA kernel that copies a block from one location to another.
//

__global__ 
void cuSimple_Copy(
	char* d_dst,
	char* d_src,
	int nx,
	int ny,
	int nz,
	size_t one_y_size
	)
{
	int thr_Y = threadIdx.y + blockIdx.y * blockDim.y;
	if (thr_Y < ny)
	{
		float* dst = (float*)(d_dst + thr_Y * one_y_size);
		float* src = (float*)(d_src + thr_Y * one_y_size);
		for (int i = threadIdx.x;  i < nz*nx*6;  i+=32)
		{
			dst[i] = src[i];
		}
	}
}

void 
Host_Simple_Copy_Kernel(
	cudaStream_t stream,
	void* d_dst,
	void* d_src,
	int nx,
	int ny,
	int nz,
	size_t one_y_size
	)
{
	dim3 blockShape(32,4,1);
	dim3 gridShape(1,(ny+3)/4,1);
	cuSimple_Copy<<<gridShape,blockShape,0,stream>>>((char*)d_dst,(char*)d_src,nx,ny,nz,one_y_size);
	//gpuErrchk( cudaPeekAtLastError() );
	//gpuErrchk( cudaDeviceSynchronize() );
}

