#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "CvxCompress.hxx"

#define PAPI
#ifndef __INTEL_COMPILER
#undef PAPI
#endif

#define VERBOSE

#ifdef PAPI
#include "papi.h"
//#define PAPI_GFLOPS
#endif

using namespace std;

int main(int argc, char* argv[])
{
	double target_snr = atof(argv[1]);
	double target_error = pow(10.0,-(target_snr/20.0));
	printf("Target SNR is %.0fdB, target error is %e\n",target_snr,target_error);

	const int nx = 800;
	const int ny = 828;
	const int nz = 573;
	const int nn = 1;
	const string filepath = "/cpfs/lfs01/ESDRD/guos/compressiontest/cube140.bin";
	const string outpath = "compressed_field_guojian.bin";
	const string decomp_outpath = "decompressed_field_guojian.bin";

	const int bx = 32;
	const int by = 32;
	const int bz = 32;
	float scale = 1e-2f;
	const bool use_local_RMS = false;
	printf("Using %s RMS.\n",use_local_RMS?"local":"global");

	int num_threads = 0;
#pragma omp parallel
	{
		num_threads = omp_get_num_threads();
	}

	long volsize = (long)nz * (long)ny * (long)nx;
	long totsize = (long)nn * volsize;
	long totsize_b = totsize * 4l;
	float *vol,*vol2,*vol3;
	posix_memalign((void**)&vol, 64, totsize_b);
	posix_memalign((void**)&vol2, 64, totsize_b);
	posix_memalign((void**)&vol3, 64, volsize*4l);
	assert(vol != 0L);
	assert(vol2 != 0L);
	assert(vol3 != 0L);
#pragma omp parallel for
	for (long ix = 0;  ix < (long)nx;  ++ix)
	{
		long nyz = (long)ny * (long)nz;
		long idx = ix * nyz;
		memset((void*)(vol+idx),0,nyz*4);
		memset((void*)(vol2+idx),0,nyz*4);
	}
	memset((void*)vol3,0,volsize*4l);

	FILE* fp = fopen(filepath.c_str(),"rb");
	assert(fp != 0L);
	fread(vol,sizeof(float),totsize,fp);
	fclose(fp);

	long length;
	unsigned int* compressed;
	posix_memalign((void**)&compressed,64,5*totsize_b/4);
	assert(compressed != 0L);
#pragma omp parallel for
	for (int ix = 0;  ix < ((5l*(long)nx)>>2);  ++ix)
	{
		long nyz = (long)ny * (long)nz;
		long idx = ix * nyz;
		memset((void*)(compressed+idx),0,nyz*4);
	}

#ifdef PAPI
	long long fpops;
	int retval = PAPI_library_init( PAPI_VER_CURRENT );
	assert(retval == PAPI_VER_CURRENT);

	int fip = 0;
	assert ( PAPI_query_event( PAPI_TOT_CYC ) == PAPI_OK );

	PAPI_shutdown(  );

#ifdef PAPI_GFLOPS
	int vec_sp_events[1]{PAPI_SP_OPS};
#else
	int vec_sp_events[1]{PAPI_TOT_CYC};
#endif
	long long vec_sp_ops = 0;
#endif
	double tot_elapsed_time = 0.0;

	FILE* fp2 = fopen("results.txt","w");
	assert(fp2 != 0L);
	printf("Starting compression test\n");
	double achieved_snr = 0.0;
	int number_of_successes = 0;
	CvxCompress* compressor = new CvxCompress();
	while (number_of_successes < 5)
	{
		int cnt_snr = 0;
		double acc_snr = 0.0;
		for (int i = 0;  i < nn;  ++i)
		{
			float* p1 = vol + (long)i * volsize;
			float* p2 = vol2 + (long)i * volsize;

			// check if wavefield has NaN's in it
			/*
			   bool has_NaN = false;
			   for (long j = 0;  j < volsize && !has_NaN;  ++j)
			   {
			   if (isnan(p1[j])) has_NaN = true;
			   }
			   assert(!has_NaN);
			   */

#ifdef PAPI
			assert( PAPI_start_counters(vec_sp_events,1) == PAPI_OK );
#endif
			struct timespec before, after;
			clock_gettime(CLOCK_REALTIME,&before);

			float ratio = compressor->Compress(scale,p1,nx,ny,nz,bx,by,bz,use_local_RMS,compressed,length);

			clock_gettime(CLOCK_REALTIME,&after);
			double elapsed1 = (double)after.tv_sec + (double)after.tv_nsec * 1e-9 - (double)before.tv_sec - (double)before.tv_nsec * 1e-9;
			double mcells_per_sec1 = (double)nx * (double)ny * (double)nz / (elapsed1 * 1e6);
			tot_elapsed_time += elapsed1;

			clock_gettime(CLOCK_REALTIME,&before);
			compressor->Decompress(p2,nx,ny,nz,compressed,length);

			clock_gettime(CLOCK_REALTIME,&after);
			double elapsed2 = (double)after.tv_sec + (double)after.tv_nsec * 1e-9 - (double)before.tv_sec - (double)before.tv_nsec * 1e-9;
			double mcells_per_sec2 = (double)nx * (double)ny * (double)nz / (elapsed2 * 1e6);
			tot_elapsed_time += elapsed2;
#ifdef PAPI
			long long curr_vec_sp_ops = 0L;
			assert( PAPI_stop_counters(&curr_vec_sp_ops,1) == PAPI_OK );
			vec_sp_ops += curr_vec_sp_ops;
#endif

#ifdef VERBOSE
			for (long j = 0;  j < volsize;  ++j)
			{
				vol3[j] += (p1[j] - p2[j]);
			}

			char str[1024];
			sprintf(str,"XZ_%04d",i+1);
			FILE* fp4 = fopen(str,"w");
			sprintf(str,"XZ_compressed_%04d",i+1);
			FILE* fp5 = fopen(str,"w");
			sprintf(str,"XZ_difference_%04d",i+1);
			FILE* fp6 = fopen(str,"w");
			sprintf(str,"XZ_rundiff_%04d",i+1);
			FILE* fp7 = fopen(str,"w");
			int iy = ny / 2;
			for (int iz = 0;  iz < nz;  ++iz)
			{
				for (int ix = 0;  ix < nx;  ++ix)
				{
					int idx = (iz * ny + iy ) * nx + ix;
					fprintf(fp4,"%d %d %e\n",ix,iz,p1[idx]);
					fprintf(fp5,"%d %d %e\n",ix,iz,p2[idx]);
					fprintf(fp6,"%d %d %e\n",ix,iz,p1[idx]-p2[idx]);
					fprintf(fp7,"%d %d %e\n",ix,iz,vol3[idx]/(double)(i+1));
				}
				fprintf(fp4,"\n");
				fprintf(fp5,"\n");
				fprintf(fp6,"\n");
				fprintf(fp7,"\n");
			}
			fclose(fp4);
			fclose(fp5);
			fclose(fp6);
			fclose(fp7);
#endif

			double acc1=0.0, acc2=0.0;
#pragma omp parallel for reduction(+:acc1,acc2)
			for (long idx = 0;  idx < volsize;  ++idx)
			{
				double val1 = p1[idx];
				double val2 = p1[idx] - p2[idx];
				acc1 += val1 * val1;
				acc2 += val2 * val2;
			}
			acc1 = sqrt(acc1/(double)volsize);
			acc2 = sqrt(acc2/(double)volsize);
			double error = acc2 / acc1;
			double snr = -20.0 * log10(error);
			acc_snr += snr;
			++cnt_snr;

			printf("vol %d, scale = %e, compression ratio = %.2f:1, compression throughput = %.0f MC/s, decompression throughput = %.0f MC/s, error = %.6e, SNR = %.0fdB\n",i+1,scale,ratio,mcells_per_sec1,mcells_per_sec2,error,snr);
			fprintf(fp2,"%d, %.2f, %.6e, %.5f, %.5f, %.0f, %.0f\n",i+1,ratio,error,elapsed1,elapsed2,mcells_per_sec1,mcells_per_sec2);
		}
		achieved_snr = acc_snr / (double)cnt_snr;
		double achieved_error = pow(10.0,-(achieved_snr/20.0));
		double scale_mulfac = target_error / achieved_error;
		scale *= scale_mulfac;
		if (abs(achieved_snr-target_snr) < 0.1)
			++number_of_successes;
		else
			number_of_successes = 0;
	}

#ifdef VERBOSE
	FILE* fp3 = fopen(outpath.c_str(),"wb");
	assert(fp3 != 0L);
	fwrite(compressed,1,length,fp3);
	fclose(fp3);

	FILE* fp4 = fopen(decomp_outpath.c_str(),"wb");
	assert(fp4 != 0L);
	fwrite(vol2,sizeof(float),volsize,fp4);
	fclose(fp4);
#endif

#ifdef PAPI
	double papi_gflops = (double)num_threads * (double)vec_sp_ops / ((double)tot_elapsed_time * 1e9);
#ifdef PAPI_GFLOPS
	printf("PAPI says we averaged %.2f GLOPS.\nTotal compression and decompression times were %.2f seconds\n",papi_gflops,tot_elapsed_time);
#else
	printf("PAPI says we averaged %.2f GCYCLES.\nTotal compression and decompression times were %.2f seconds\n",papi_gflops,tot_elapsed_time);
#endif
#else
	printf("Total compression and decompression times were %.2f seconds\n",tot_elapsed_time);
#endif

	fclose(fp2);

	return 0;
}
