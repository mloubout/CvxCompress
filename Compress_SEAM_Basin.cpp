#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "CvxCompress.hxx"
#include <iostream>
using namespace std;

//#define VERBOSE

using namespace std;

int main(int argc, char* argv[])
{
  int nx = 320; //slow
  int ny = 416; //medium
  int nz = 352; //fast
  const string filepath = "/cpfs/lfs02/ESDRD/akep/CompressionTest/SEAM1_le.bin";
  const string outpath = "compressed_seam_basin.bin";

  int bx = 32;
  int by = 32;
  int bz = 32;
  const float scale = 1e-2f;
  const bool use_local_RMS = false;
  const bool include_block_lengths = false;
  printf("Using %s RMS.\n",use_local_RMS?"local":"global");

  long overall_compressed = 0;

  int num_threads = 0;
#pragma omp parallel
  {
    num_threads = omp_get_num_threads();
  }

  long volsize = (long)nz * (long)ny * (long)nx;
  long totsize = volsize;
  long totsize_b = totsize * 4l;
  long compressed_array_size = 5l*totsize_b/4l;

  float *vol,*vol2;
  posix_memalign((void**)&vol,  64, totsize_b);
  posix_memalign((void**)&vol2, 64, totsize_b);
  assert(vol != 0L);
  assert(vol2 != 0L);
#pragma omp parallel for
  for (long ix = 0;  ix < (long)nx;  ++ix)
    {
      long nyz = (long)ny * (long)nz;
      long idx = ix * nyz;
      memset((void*)(vol+idx),0,nyz*4);
      memset((void*)(vol2+idx),0,nyz*4);
    }

  FILE* fp = fopen(filepath.c_str(),"rb");
  assert(fp != 0L);
  size_t items_read = fread(vol,sizeof(float),totsize,fp);
  fclose(fp);

  // DEBUG
  double inMS=0;
  for(long i=0; i<totsize; i++) {
    inMS += vol[i]*vol[i];
  }
  double inRMS = sqrt(inMS/totsize);
  cout << "  DEBUG: Input RMS is " << inRMS << endl;
  double norm = 1.0/inRMS;
  for(long i=0; i<totsize; i++) {
    vol[i] *= norm;
  }
  inMS = 0;
  for(long i=0; i<totsize; i++) {
    inMS += vol[i]*vol[i];
  }
  cout << "  DEBUG: Scaled input RMS is " << sqrt(inMS/totsize) << endl;
  cout << "  DEBUG: Reading from file: "  <<  filepath.c_str() <<endl;
  cout << "  DEBUG: Number of floats read from input                                          : " << items_read << endl; 
  cout << "  DEBUG: Expected number of floats from dimensions (n1, n2, n3) = ( " << nz << ", " << ny << ", " << nz << " ): " << volsize << endl;   

  unsigned int* compressed;
  posix_memalign((void**)&compressed,64, compressed_array_size);
  assert(compressed != 0L);
#pragma omp parallel for
  for (int ix = 0;  ix < ((5l*(long)nx)>>2);  ++ix)
    {
      long nyz = (long)ny * (long)nz;
      long idx = ix * nyz;
      memset((void*)(compressed+idx),0,nyz*4);
    }

  double tot_elapsed_time = 0.0;

  FILE* fp2 = fopen("results.txt","w");
  assert(fp2 != 0L);
  printf("Starting compression test\n");
  CvxCompress* compressor = new CvxCompress();
  long compressed_length;

  // check if wavefield has NaN's in it
  bool has_NaN = false;
  for (long j = 0;  j < volsize && !has_NaN;  ++j)
    {
      if (isnan(vol[j])) has_NaN = true;
    }
  assert(!has_NaN);
  
  struct timespec before, after;
  clock_gettime(CLOCK_REALTIME,&before);
  // **********************  COMPRESSING **********************  
  // float ratio = compressor->Compress_Safe(scale,vol,nz,ny,nx,bz,by,bx,
  // 					      use_local_RMS,include_block_lengths,
  // 					      compressed,

  float ratio = compressor->Compress(scale, vol,
				     nz,ny,nx,
				     bz,by,bx,
				     use_local_RMS,
				     compressed,
				     compressed_length);

  clock_gettime(CLOCK_REALTIME,&after);
  double elapsed1 = (double)after.tv_sec + (double)after.tv_nsec * 1e-9 - (double)before.tv_sec - (double)before.tv_nsec * 1e-9;
  double mcells_per_sec1 = (double)nx * (double)ny * (double)nz / (elapsed1 * 1e6);
  tot_elapsed_time += elapsed1;
  overall_compressed += (long)compressed_length;
	
  clock_gettime(CLOCK_REALTIME,&before);
  // **********************  DECOMPRESSING **********************  
  //      compressor->Decompress_Safe(vol2,nz,ny,nx,compressed,compressed_length);
  compressor->Decompress(vol2, nz,ny,nx, compressed, compressed_length);
  
  // DEBUG
  unsigned long sum_compressed;
  for(int ii=0; ii<compressed_array_size/sizeof(unsigned int) ; ii++) sum_compressed = compressed[ii]; 
  printf("  DEBUG: Sum of all elements in unsigned int array containing compressed data: %d\n", sum_compressed);

  clock_gettime(CLOCK_REALTIME,&after);
  double elapsed2 = (double)after.tv_sec + (double)after.tv_nsec * 1e-9 - (double)before.tv_sec - (double)before.tv_nsec * 1e-9;
  double mcells_per_sec2 = (double)nx * (double)ny * (double)nz / (elapsed2 * 1e6);
  tot_elapsed_time += elapsed2;

  double acc1=0.0, acc2=0.0, acc3=0.0;
#pragma omp parallel for reduction(+:acc1,acc2,acc3)
  for (long idx = 0;  idx < volsize;  ++idx)
    {
      double val1 = vol[idx];
      double val2 = vol[idx] - vol2[idx];
      double val3 = vol2[idx];
      acc1 += val1 * val1;
      acc2 += val2 * val2;
      acc2 += val3 * val3;
    }
  acc1 = sqrt(acc1/(double)volsize);
  acc2 = sqrt(acc2/(double)volsize);
  acc3 = sqrt(acc3/(double)volsize);
  printf("RMS:\n input:     %f\n output:     %f\n Difference: %f\n",acc1,acc3,acc2);
  double error = acc2 / acc1;
  double snr = -20.0 * log10(error);

  printf("compression ratio = %.2f:1, compression throughput = %.0f MC/s, decompression throughput = %.0f MC/s, error = %.6e, SNR = %.1f dB\n",ratio,mcells_per_sec1,mcells_per_sec2,error,snr);
  fprintf(fp2,"%.2f, %.6e, %.5f, %.5f, %.0f, %.0f\n",ratio,error,elapsed1,elapsed2,mcells_per_sec1,mcells_per_sec2);

  printf("Total compression and decompression times were %.2f seconds\n",tot_elapsed_time);
  double overall_ratio = ((double)nx * (double)ny * (double)nz * 4.0) / (double)overall_compressed;
  printf("Total compression ratio (all snapshots) was %.2f:1\n",overall_ratio);

  fclose(fp2);

  return 0;
}
