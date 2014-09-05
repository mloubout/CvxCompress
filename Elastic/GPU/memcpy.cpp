#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <omp.h>
#include <xmmintrin.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <cuda_runtime_api.h>

#define NUM_PAGES 1

void omp_memclear(void* dst, size_t len)
{
	size_t leni = len / 16;
	size_t nn = (len + NUM_PAGES*4096-1) / (NUM_PAGES*4096);
#pragma omp parallel for schedule(static)
	for (int i = 0;  i < nn;  ++i)
	{
		__m128 zero = _mm_set_ps(0.0f, 0.0f, 0.0f, 0.0f);
		size_t i0 = (size_t)i * (NUM_PAGES*256);
		size_t in = leni - i0;
		if (in > (NUM_PAGES*256)) in = (NUM_PAGES*256);
		__m128* d = (__m128*)dst + i0;
		for (int j = 0;  j < in;  ++j)
		{
			_mm_stream_ps((float*)(d+j), zero);
		}
	}
}

void omp_memcpy(void* dst, void* src, size_t len)
{
	size_t leni = len / 16;
	size_t nn = (len + NUM_PAGES*4096-1) / (NUM_PAGES*4096);
#pragma omp parallel for schedule(static)
	for (int i = 0;  i < nn;  ++i)
	{
		size_t i0 = (size_t)i * (NUM_PAGES*256);
		size_t in = leni - i0;
		if (in > (NUM_PAGES*256)) in = (NUM_PAGES*256);
		__m128* d = (__m128*)dst + i0;
		__m128* s = (__m128*)src + i0;
		for (int j = 0;  j < in;  ++j)
		{
			_mm_stream_ps((float*)(d+j),_mm_load_ps((float*)(s+j)));
		}
	}
}

//
// Returns true if test passed.
// 
bool omp_memcpy_module_test()
{
	const size_t blksize_v = 124567;
	const size_t blksize_f = blksize_v * 4;
	const size_t blksize_b = blksize_f * 4;

	const size_t blksize_pad_v = blksize_v + 4;
	const size_t blksize_pad_f = blksize_pad_v * 4;
	const size_t blksize_pad_b = blksize_pad_f * 4;

	void *src, *dst;
	posix_memalign(&src, getpagesize(), blksize_pad_b);
	posix_memalign(&dst, getpagesize(), blksize_pad_b);
	memset(src, 0, blksize_pad_b);
	memset(dst, 0, blksize_pad_b);

	for (int i = 0;  i < blksize_f;  ++i) ((int*)src)[i] = i;
	omp_memcpy(dst, src, blksize_b);

	bool fault = false;
	for (int i = 0;  i < blksize_f && !fault;  ++i)
	{
		int val = ((int*)dst)[i];
		if (val != i)
		{
			printf("Error! Expected %d, found %d\n",i,val);
			fault = true;
		}
	}
	for (int i = blksize_f;  i < blksize_pad_f && !fault;  ++i)
	{
		int val = ((int*)dst)[i];
		if (val != 0)
		{
			printf("Error! Memcpy overwrote trailing region.\n");
			fault = true;
		}
	}

	free(dst);
	free(src);
	
	return !fault;
}

bool Print_Device_Stats(int device_id, double& TFLOPS, double& GB_per_s)
{
	TFLOPS = GB_per_s = 0.0;
	cudaDeviceProp devProps;
	cudaError_t err = cudaGetDeviceProperties(&devProps, device_id);
	if (err == cudaSuccess)
	{
		cudaSetDevice(device_id);
		size_t free,total;
		err = cudaMemGetInfo(&free,&total);
		if (err == cudaSuccess)
		{
			double dFree_MB = (double)free / 1048576.0;
			GB_per_s = (double)devProps.memoryBusWidth * (double)devProps.memoryClockRate / 4e6;
			int Cores_per_SM;
			if (devProps.major == 1)
			{
				Cores_per_SM = 8;
			}
			else if (devProps.major == 2)
			{
				if (devProps.minor == 1)
				{
					Cores_per_SM = 48;
				}
				else
				{
					Cores_per_SM = 32;
				}
			}
			else if (devProps.major == 3)
			{
				Cores_per_SM = 192;
			}
			TFLOPS = (double)devProps.clockRate * (double)devProps.multiProcessorCount * (double)Cores_per_SM / 5e8;
			printf("device_id %d :: %s, CC=%d.%d, Free Mem=%.2f MB, %.3f TFLOPS, %.0f GB/s\n",device_id,devProps.name,devProps.major,devProps.minor,dFree_MB,TFLOPS,GB_per_s);
			return true;
		}
	}
	return false;
}

bool Check_GPUs(int* device_id, int num_devices)
{
	printf("\n");
	int device_count = 0;
	cudaGetDeviceCount(&device_count);
	if (device_count < 1)
	{
		printf("No CUDA capable devices found!\n\n");
		return false;
	}
	double Total_GB_per_s=0.0, Total_TFLOPS=0.0;
	for (int i = 0;  i < num_devices;  ++i)
	{
		double GB_per_s, TFLOPS;
		if (!Print_Device_Stats(device_id[i],TFLOPS,GB_per_s))
		{
			printf("device_id %d not found\n\n",device_id[i]);
			return false;
		}
		Total_GB_per_s += GB_per_s;
		Total_TFLOPS += TFLOPS;
	}
	printf("Aggregate %.3f TFLOPS, %.0f GB/s\n\n",Total_TFLOPS,Total_GB_per_s);
	return true;
}

void cuda_host_memalign(char** p, size_t alignment, size_t len)
{
	posix_memalign((void**)p, alignment, len);
	omp_memclear((void*)(*p), len);
	cudaHostRegister((void*)(*p), len, 0);
}

void cuda_host_free(void* p)
{
	cudaHostUnregister(p);
	free(p);
}

int main(int argc, char* argv[])
{
	//setpriority(PRIO_PROCESS, 0, 19);

	const size_t nx = 1720;
	const size_t ny = 1720;
	const size_t nz = 1024;

	size_t niter = atoi(argv[2]);
	omp_set_num_threads(atoi(argv[1]));

	if (omp_memcpy_module_test())
	{
		printf("omp_mempcy module test PASSED!\n");

		size_t nbX = nx / 4;
		size_t blkSize_pt = (size_t)4 * (size_t)ny * (size_t)nz;
		size_t blkSize_wf_f = blkSize_pt * (size_t)12;
		size_t blkSize_wf_b = blkSize_wf_f * (size_t)4;
		size_t blkSize_wf_half_b = blkSize_wf_b / 2;
		// pad to whole number of vm pages
		size_t ps = getpagesize();
		size_t blkSize_wf_pad_b = ((blkSize_wf_b + ps - 1) / ps) * ps;
		// total size of wavefield buffer
		size_t wf_pad_b = blkSize_wf_pad_b * nbX;
		size_t wf_pad_l = wf_pad_b / 8;
	
		size_t blkSize_em_b = blkSize_pt * (size_t)16;
		size_t blkSize_em_half_b = blkSize_em_b / 2;
		size_t blkSize_em_pad_b = ((blkSize_em_b + ps - 1) / ps) * ps;
		size_t em_pad_b = blkSize_em_pad_b * nbX;
		size_t em_pad_l = em_pad_b / 8;

		// wake up a couple of GPUs
		int* device_id = new int[2];
		device_id[0] = 0;
		device_id[1] = 8;
		Check_GPUs(device_id,2);

		printf("Allocating transfer buffers Host <-> GPU...\n");
		char *h_wf_in_1, *h_wf_in_2, *h_wf_out_1, *h_wf_out_2, *h_em_in_1, *h_em_in_2;
		cuda_host_memalign(&h_wf_in_1, ps, blkSize_wf_pad_b);
		cuda_host_memalign(&h_wf_in_2, ps, blkSize_wf_pad_b);
		cuda_host_memalign(&h_wf_out_1, ps, blkSize_wf_pad_b);
		cuda_host_memalign(&h_wf_out_2, ps, blkSize_wf_pad_b);
		cuda_host_memalign(&h_em_in_1, ps, blkSize_em_pad_b);
		cuda_host_memalign(&h_em_in_2, ps, blkSize_em_pad_b);

		printf("Allocating interleaved wavefields (%.2fGB)...\n",(double)wf_pad_b*1e-9);
		char* wf = 0L;
		posix_memalign((void**)&wf, ps, wf_pad_b);
		// initialize buffers block-by-block to ensure vm pages are pegged to the right socket
		for (size_t i = 0;  i < nbX;  ++i) omp_memclear(wf + i * blkSize_wf_pad_b, blkSize_wf_pad_b);

		//omp_memclear(wf, wf_pad_b);
		printf("Allocating earth model (%.2fGB)...\n",(double)em_pad_b*1e-9);
		char* em = 0L;
		posix_memalign((void**)&em, ps, em_pad_b);
		// initialize buffers block-by-block to ensure vm pages are pegged to the right socket
		for (size_t i = 0;  i < nbX;  ++i) omp_memclear(em + i * blkSize_em_pad_b, blkSize_em_pad_b);
		//omp_memclear(em, em_pad_b);

		printf("Allocating device buffers...\n");
		void *d1_wf_1, *d1_wf_2, *d1_em_1, *d1_em_2;
		cudaSetDevice(device_id[0]);
		cudaMalloc(&d1_wf_1, blkSize_wf_half_b);
		cudaMalloc(&d1_wf_2, blkSize_wf_half_b);
		cudaMalloc(&d1_em_1, blkSize_em_half_b);
		cudaMalloc(&d1_em_2, blkSize_em_half_b);

		cudaStream_t inp1_stream;
		cudaStream_t out1_stream;
		cudaStreamCreate(&inp1_stream);
		cudaStreamCreate(&out1_stream);

		void *d2_wf_1, *d2_wf_2, *d2_em_1, *d2_em_2;
		cudaSetDevice(device_id[1]);
		cudaMalloc(&d2_wf_1, blkSize_wf_half_b);
		cudaMalloc(&d2_wf_2, blkSize_wf_half_b);
		cudaMalloc(&d2_em_1, blkSize_em_half_b);
		cudaMalloc(&d2_em_2, blkSize_em_half_b);

		cudaStream_t inp2_stream;
		cudaStream_t out2_stream;
		cudaStreamCreate(&inp2_stream);
		cudaStreamCreate(&out2_stream);

		// fill buffers for module test
#pragma omp parallel sections
		{
#pragma omp section
			{
				for (size_t i = 0;  i < wf_pad_l;  ++i) ((long*)wf)[i] = i;
			}
#pragma omp section
			{
				for (size_t i = 0;  i < em_pad_l;  ++i) ((long*)em)[i] = i;
			}
		}

		struct timespec before, after;
		clock_gettime(CLOCK_REALTIME, &before);

		size_t i0 = 2, iteration = 1;
		for (size_t i = 0;  i < niter*nbX+3;  ++i)
		{
			//printf("iteration %ld\n",i);

			// rotate buffers
			void* tmp;

			tmp = h_wf_in_1;
			h_wf_in_1 = h_wf_in_2;
			h_wf_in_2 = (char*)tmp;

			tmp = h_wf_out_1;
			h_wf_out_1 = h_wf_out_2;
			h_wf_out_2 = (char*)tmp;

			tmp = h_em_in_1;
			h_em_in_1 = h_em_in_2;
			h_em_in_2 = (char*)tmp;

			tmp = d1_wf_1;
			d1_wf_1 = d1_wf_2;
			d1_wf_2 = tmp;

			tmp = d1_em_1;
			d1_em_1 = d1_em_2;
			d1_em_2 = tmp;

			tmp = d2_wf_1;
			d2_wf_1 = d2_wf_2;
			d2_wf_2 = tmp;

			tmp = d2_em_1;
			d2_em_1 = d2_em_2;
			d2_em_2 = tmp;

			size_t i_inp = i % nbX;

			cudaMemcpyAsync(d1_wf_1, h_wf_in_2, blkSize_wf_half_b, cudaMemcpyHostToDevice, inp1_stream);
			cudaMemcpyAsync(d1_em_1, h_em_in_2, blkSize_em_half_b, cudaMemcpyHostToDevice, inp1_stream);
			cudaMemcpyAsync(h_wf_out_1, d1_wf_2, blkSize_wf_half_b, cudaMemcpyDeviceToHost, out1_stream);
			
			cudaMemcpyAsync(d2_wf_1, h_wf_in_2 + blkSize_wf_half_b, blkSize_wf_half_b, cudaMemcpyHostToDevice, inp2_stream);
			cudaMemcpyAsync(d2_em_1, h_em_in_2 + blkSize_em_half_b, blkSize_em_half_b, cudaMemcpyHostToDevice, inp2_stream);
			cudaMemcpyAsync(h_wf_out_1 + blkSize_wf_half_b, d2_wf_2, blkSize_wf_half_b, cudaMemcpyDeviceToHost, out2_stream);

			omp_memcpy(h_wf_in_1, wf + i_inp * blkSize_wf_pad_b, blkSize_wf_pad_b);
			omp_memcpy(h_em_in_1, em + i_inp * blkSize_em_pad_b, blkSize_em_pad_b);
	
			if (i > 2)
			{
				size_t i_out = (i - 3) % nbX;
				omp_memcpy(wf + i_out * blkSize_wf_pad_b, h_wf_out_2, blkSize_wf_pad_b);
			}

			cudaStreamSynchronize(inp1_stream);
			cudaStreamSynchronize(out1_stream);

			cudaStreamSynchronize(inp2_stream);
			cudaStreamSynchronize(out2_stream);

			if (i >= i0 && (i-i0) >= nbX)
			{
				clock_gettime(CLOCK_REALTIME, &after);
				double elapsed_time = (double)after.tv_sec + (double)after.tv_nsec * 1e-9 - (double)before.tv_sec - (double)before.tv_nsec * 1e-9;
				double BW_MCells = (double)nx * (double)ny * (double)nz / (elapsed_time * 1e6);
				printf("Iteration %04ld :: Elapsed time = %.2fs, BW = %.2f MCells/s\n",iteration,elapsed_time,BW_MCells);
				before = after;
				i0 += nbX;
				++iteration;
			}
		}

		// check memory for correctness
#pragma omp parallel sections
		{
#pragma omp section
			{
				bool fault = false;
				for (size_t i = 0;  i < wf_pad_l && !fault;  ++i)
				{
					if (((long*)wf)[i] != i)
					{
						printf("WF fault found at index %ld (found %ld, expected %ld)\n",i,((long*)wf)[i],i);
						for (int j = 0;  j < 4;  ++j) printf("%ld\n",((long*)wf)[i+j]);
						fault = true;
					}
				}
				if (!fault) printf("WF test PASSED!\n");
			}
#pragma omp section
			{
				bool fault = false;
				for (size_t i = 0;  i < em_pad_l && !fault;  ++i)
				{
					if (((long*)em)[i] != i)
					{
						printf("EM fault found at index %ld\n",i);
						fault = true;
					}
				}
				if (!fault) printf("EM test PASSED!\n");
			}
		}

		// die
		cudaFree(d1_em_2);
		cudaFree(d1_em_1);
		cudaFree(d1_wf_2);
		cudaFree(d1_wf_1);
		
		cudaFree(d2_em_2);
		cudaFree(d2_em_1);
		cudaFree(d2_wf_2);
		cudaFree(d2_wf_1);

		cuda_host_free(h_em_in_2);
		cuda_host_free(h_em_in_1);
		cuda_host_free(h_wf_out_2);
		cuda_host_free(h_wf_out_1);
		cuda_host_free(h_wf_in_2);
		cuda_host_free(h_wf_in_1);
		
		free(em);
		free(wf);

		/*
		const size_t blksize = 160 * 1024 * 1024;
		int niter = atoi(argv[2]);

		omp_set_num_threads(atoi(argv[1]));

		void *blk1, *blk2;
		posix_memalign(&blk1, getpagesize(), blksize);
		posix_memalign(&blk2, getpagesize(), blksize);
		omp_memclear(blk1, blksize);  // these two calls bind the pages
		omp_memclear(blk2, blksize);
		for (int i = 0;  i < blksize;  ++i)
		{
			((char*)blk2)[i] = (char)(i%127);
		}

		struct timespec before, after;
		clock_gettime(CLOCK_REALTIME, &before);

		for (int i = 0;  i < niter;  ++i)
		{
			omp_memcpy(blk1, blk2, blksize);
		}

		clock_gettime(CLOCK_REALTIME, &after);
		double elapsed_time = (double)after.tv_sec + (double)after.tv_nsec * 1e-9 - (double)before.tv_sec - (double)before.tv_nsec * 1e-9;
		double BW_GB = (double)blksize * (double)niter * 2 / (elapsed_time * 1e9);
		printf("Elapsed time = %.2fs, BW = %.2f GB/s\n",elapsed_time,BW_GB);

		size_t acc = 0;
		for (int i = 0;  i < blksize;  ++i) acc += ((char*)blk1)[i];
		printf("acc = %ld\n",acc);
		*/
	}
	
	return 0;
}

