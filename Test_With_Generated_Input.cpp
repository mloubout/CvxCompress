#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <time.h>
#include <chrono>
#include "CvxCompress.hxx"
#include <iostream>
using namespace std;

using std::chrono::high_resolution_clock;
typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<double> fsec;

int main(int argc, char* argv[])
{
  int bx = 32;
  int by = 32;
  int bz = 32;
  const float scale = 1e-2f;
  const bool use_local_RMS = false;
  const bool include_block_lengths = false;
  printf("Using %s RMS.\n",use_local_RMS?"local":"global");
  for(int itries=1; itries<=3; itries++){
    int nx = 320*itries; //slow
    int ny = 416*itries; //medium
    int nz = 352*itries; //fast
    printf("Using size of input (nslow, middle, fast) = ( %d, %d, %d ) \n", nx, ny, nz);

    long overall_compressed = 0;

    long volsize = (long)nz * (long)ny * (long)nx;
    long totsize = volsize;               // size in # of floats
    long totsize_b = totsize * 4l;        // size in bytes
    long compressed_array_size = totsize; // same size as input should generally be safe

    float *vol,*vol2;
    posix_memalign((void**)&vol,  64, totsize_b);   assert(vol != 0L);
    posix_memalign((void**)&vol2, 64, totsize_b);   assert(vol2 != 0L);
    memset((void*)vol2,0,totsize_b);
    for (long ix = 0;  ix < (long)nx;  ++ix) {
      long nyz = (long)ny * (long)nz;
      long idx = ix * nyz;
      // Set input here
      float xval = sin(ix*M_PI/(nx/itries)*10); // 10 perios of sin, constant x-slices
      for(int ii = 0; ii<nyz; ii++) vol[idx+ii] = xval; 
    }
    unsigned int* compressed;
    posix_memalign((void**)&compressed, 64, compressed_array_size);  assert(compressed != 0L);
    memset((void*)(compressed), 0, compressed_array_size);

    double tot_elapsed_time = 0.0;

    printf("Starting compression test\n");
    CvxCompress* compressor = new CvxCompress();
    long compressed_length;

    // check if wavefield has NaN's in it
    bool has_NaN = false;
    for (long j = 0;  j < volsize && !has_NaN;  ++j) if (isnan(vol[j])) has_NaN = true;
    assert(!has_NaN);
  
    // struct timespec before, after;
    // clock_gettime(CLOCK_REALTIME,&before);
    auto start = Time::now();
    // **********************  COMPRESSING **********************  
    // float ratio = compressor->Compress_Safe(scale,vol,nz,ny,nx,bz,by,bx,
    // 					      use_local_RMS, include_block_lengths,
    // 					      (char *)compressed, compressed_array_size,  compressed_length);

    float ratio = compressor->Compress(scale, vol,
				       nz,ny,nx,
				       bz,by,bx,
				       use_local_RMS,
				       compressed,
				       compressed_length);

    // clock_gettime(CLOCK_REALTIME,&after);
    auto stop = Time::now();
    // double elapsed1 = (double)after.tv_sec + (double)after.tv_nsec * 1e-9 - (double)before.tv_sec - (double)before.tv_nsec * 1e-9;
    // double mcells_per_sec1 = (double) volsize / (elapsed1 * 1e6);
    fsec elapsed1 = stop-start;
    double mcells_per_sec1 = (double) volsize / (elapsed1.count() * 1e6);
    tot_elapsed_time += elapsed1.count();
    overall_compressed = (long)compressed_length;
	
    // clock_gettime(CLOCK_REALTIME,&before);
    start = Time::now();
    // **********************  DECOMPRESSING **********************  
    //  compressor->Decompress_Safe(vol2,nz,ny,nx, (char*) compressed, compressed_length);
    compressor->Decompress(vol2, nz,ny,nx, compressed, compressed_length);
  
    // clock_gettime(CLOCK_REALTIME,&after);
    stop = Time::now();
    // double elapsed2 = (double)after.tv_sec + (double)after.tv_nsec * 1e-9 - (double)before.tv_sec - (double)before.tv_nsec * 1e-9;
    // double mcells_per_sec2 = (double) volsize / (elapsed2 * 1e6);
    fsec elapsed2 = stop-start;
    double mcells_per_sec2 = (double) volsize / (elapsed2.count() * 1e6);
    tot_elapsed_time += elapsed2.count();

    double acc1=0.0, acc2=0.0, acc3=0.0;
    for (long idx = 0;  idx < volsize;  ++idx)
      {
	double val1 = vol[idx];
	double val2 = vol[idx] - vol2[idx];
	double val3 = vol2[idx];
	acc1 += val1 * val1;
	acc2 += val2 * val2;
	acc3 += val3 * val3;
      }
    acc1 = sqrt(acc1/(double)volsize);
    acc2 = sqrt(acc2/(double)volsize);
    acc3 = sqrt(acc3/(double)volsize);

    printf("RMS:\n input :     %f\n output:     %f\n Difference: %f\n",acc1,acc3,acc2);
    double error = acc2 / acc1;
    double snr = -20.0 * log10(error);

    printf("compression ratio (return value) = %.2f:1, compression throughput = %.0f MC/s, decompression throughput = %.0f MC/s, error = %.6e, SNR = %.1f dB\n",ratio,mcells_per_sec1,mcells_per_sec2,error,snr);
    printf("Total compression and decompression times were %.2f seconds\n",tot_elapsed_time);
    double overall_ratio = ((double)volsize * 4.0) / (double)overall_compressed;
    printf("Total compression ratio (based on compress_length) was %.2f:1, compressed length in bytes = %d \n",overall_ratio, overall_compressed);
  }// itries
  return 0;
}
