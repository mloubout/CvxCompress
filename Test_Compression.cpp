#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <chrono>
#include "CvxCompress.hxx"
#include "Read_Raw_Volume.hxx"

// #define PAPI
#ifndef __INTEL_COMPILER
#undef PAPI  // PAPI only works with intel compiler (for now).
#endif

#ifdef PAPI
#include "papi.h"
#endif

using std::chrono::system_clock;
typedef std::chrono::system_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<double> fsec;

void XZ_Slice(const char* filename, int nx, int ny, int nz, float* vol, int iy)
{
	FILE* fp8 = fopen(filename, "w");
	if (fp8 != 0L)
	{
		for (long iz = 0;  iz < nz;  ++iz)
		{
			for (long ix = 0;  ix < nx;  ++ix)
			{
				fprintf(fp8,"%ld %ld %e\n",ix,iz,vol[iz*ny*nx+iy*nx+ix]);
			}
			fprintf(fp8,"\n");
		}
		fclose(fp8);
	}
}

float Compute_RMS(float* vol, int nx, int ny, int nz)
{
	double acc = 0.0;
#pragma omp parallel for reduction(+:acc)
	for (int iz = 0;  iz < nz;  ++iz)
	{
		float* p = vol + (long)iz * (long)nx * (long)ny;
		for (int i = 0;  i < nx*ny;  ++i) acc += p[i] * p[i];
	}
	return sqrt(acc/((double)nx*(double)ny*(double)nz));
}

int main(int argc, char* argv[])
{
	if (argc != 8)
	{
		printf("Usage : %s <raw-volume-file> <scale> <bx> <by> <bz> <compressed-volume-file> <num-loops>\n",argv[0]);
		return -1;
	}

	const bool use_local_RMS = false;
	printf("Using %s RMS.\n",use_local_RMS?"local":"global");

	int num_threads = 0;
#pragma omp parallel 
	{
		num_threads = omp_get_num_threads();
	}

#ifdef PAPI
        int retval = PAPI_library_init( PAPI_VER_CURRENT );
        assert(retval == PAPI_VER_CURRENT);
	int sp_ops_events[1]{PAPI_SP_OPS};
	bool sp_ops_counter_is_available = (PAPI_query_event(sp_ops_events[0]) == PAPI_OK);
	if (sp_ops_counter_is_available)
	{
		assert(PAPI_thread_init((unsigned long (*)(void))(omp_get_num_threads)) == PAPI_OK);
	}
	else
	{
		printf("PAPI does not support PAPI_VEC_SP counter on this machine.\n");
		PAPI_shutdown();
	}
#endif

	printf("Reading raw volume from file %s...\n",argv[1]);
	float scale = atof(argv[2]);
	printf("Scale is %.3e\n",scale);
	int bx = atoi(argv[3]);
	int by = atoi(argv[4]);
	int bz = atoi(argv[5]);
	int num_loops = atoi(argv[7]);

	int nx,ny,nz;
	float* vol = 0L;
	Read_Raw_Volume(argv[1],nx,ny,nz,vol);

#ifdef PAPI
	long long sp_ops1 = 0, sp_ops2 = 0;
	double total_elapsed_time1 = 0.0, total_elapsed_time2 = 0.0;
#endif
	for (int iLoop = 0;  iLoop < num_loops;  ++iLoop)
	{
		unsigned int* compressed = 0L;
		posix_memalign((void**)&compressed, 64, ((long)sizeof(float)*(long)nx*(long)ny*(long)nz*5l)>>2);

		CvxCompress* compressor = new CvxCompress();

		printf("Compressing.\n");
#ifdef PAPI
		if (sp_ops_counter_is_available)
		{
#pragma omp parallel for schedule(static,1)
			for (int iThr = 0;  iThr < num_threads;  ++iThr)
			{
				PAPI_start_counters(sp_ops_events,1);
			}
		}
#endif
		long compressed_length = 0;
		auto start = Time::now();
		float ratio = compressor->Compress(scale,vol,nx,ny,nz,bx,by,bz,use_local_RMS,compressed,compressed_length);
		auto stop = Time::now();
		fsec elapsed = stop-start;
		double mcells_per_sec = (double)nx * (double)ny * (double)nz / (elapsed.count() * 1e6);
#ifdef PAPI
		if (sp_ops_counter_is_available)
		{
#pragma omp parallel for schedule(static,1) reduction(+:sp_ops1)
			for (int iThr = 0;  iThr < num_threads;  ++iThr)
			{
				long long curr_sp_ops = 0L;
				PAPI_stop_counters(&curr_sp_ops,1);
				sp_ops1 += curr_sp_ops;
			}
			total_elapsed_time1 += elapsed.count();
		}
#endif
		printf("Compression throughput was %.0f MC/s - Compression ratio is %.2f:1\n",mcells_per_sec,ratio);

		printf("Decompressing.\n");
		int nx2,ny2,nz2;
#ifdef PAPI
		if (sp_ops_counter_is_available)
		{
#pragma omp parallel for schedule(static,1)
			for (int iThr = 0;  iThr < num_threads;  ++iThr)
			{
				PAPI_start_counters(sp_ops_events,1);
			}
		}
#endif
		start = Time::now();
		float* vol2 = compressor->Decompress(nx2,ny2,nz2,compressed,compressed_length);
		stop = Time::now();
		elapsed = stop-start; 
		mcells_per_sec = (double)nx2 * (double)ny2 * (double)nz2 / (elapsed.count() * 1e6);
#ifdef PAPI
		if (sp_ops_counter_is_available)
		{
#pragma omp parallel for schedule(static,1) reduction(+:sp_ops2)
			for (int iThr = 0;  iThr < num_threads;  ++iThr)
			{
				long long curr_sp_ops = 0L;
				PAPI_stop_counters(&curr_sp_ops,1);
				sp_ops2 += curr_sp_ops;
			}
			total_elapsed_time2 += elapsed.count();
		}
#endif
		printf("Decompression throughput was %.0f MC/s\n",mcells_per_sec);

		printf("Generating error volume.\n");
		float *vol3;
		posix_memalign((void**)&vol3,64,(long)sizeof(float)*(long)nx*(long)ny*(long)nz);
		double acc = 0.0, acc2 = 0.0;
#pragma omp parallel for reduction(+:acc,acc2)
		for (int iz = 0;  iz < nz;  ++iz)
		{
			float* p1 = vol + (long)iz * (long)nx * (long)ny;
			float* p2 = vol2 + (long)iz * (long)nx * (long)ny;
			float* p3 = vol3 + (long)iz * (long)nx * (long)ny;
			for (int i = 0;  i < nx*ny;  ++i)
			{
				p3[i] = p2[i] - p1[i];
				acc += p1[i] * p1[i];
				acc2 += p3[i] * p3[i];
			}
		}
		float rms_orig = sqrt(acc/((double)nx*(double)ny*(double)nz));
		float rms_diff = sqrt(acc2/((double)nx*(double)ny*(double)nz));
		printf("Error is %.3e\n",rms_diff/rms_orig);

		if (iLoop == (num_loops-1))
		{
			printf("Write XZ slice original.\n");
			XZ_Slice("original.txt",nx,ny,nz,vol,ny-1);
			XZ_Slice("compressed.txt",nx,ny,nz,vol2,ny-1);
			XZ_Slice("difference.txt",nx,ny,nz,vol3,ny-1);

			FILE* fp = fopen(argv[6],"w");
			if (fp != 0L)
			{
				printf("Writing compressed volume to %s...\n",argv[6]);
				fwrite(compressed,1,compressed_length,fp);
				fclose(fp);
			}
		}

		delete compressor;
		free(vol3);
		free(vol2);
		free(compressed);
	}
#ifdef PAPI
	if (sp_ops_counter_is_available)
	{
		// PAPI only counted FLOPS for master thread.
		// For simplicity, we assume all threads did the same amount of work,
		// thus we just multiply by number of threads to estimate GFLOPS for the whole machine.
		double GFLOPS_compress = (double)sp_ops1 / (total_elapsed_time1 * 1e9);
		double GFLOPS_decompress = (double)sp_ops2 / (total_elapsed_time2 * 1e9);
		printf("PAPI says we averaged %.2f GFLOPS during compression, which took a total of %.2f seconds.\n",GFLOPS_compress,total_elapsed_time1);
		printf("PAPI says we averaged %.2f GFLOPS during decompression, which took a total of %.2f seconds.\n",GFLOPS_decompress,total_elapsed_time2);
	}
#endif

	free(vol);
	return 0;
}
