#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "CvxCompress.hxx"
#include "Read_Raw_Volume.hxx"

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

	for (int iLoop = 0;  iLoop < num_loops;  ++iLoop)
	{
		unsigned int* compressed = 0L;
		posix_memalign((void**)&compressed, 64, ((long)sizeof(float)*(long)nx*(long)ny*(long)nz*5l)>>2);

		CvxCompress* compressor = new CvxCompress();

		printf("Compressing.\n");	
		long compressed_length = 0;
		struct timespec before, after;
		clock_gettime(CLOCK_REALTIME,&before);
		float ratio = compressor->Compress(scale,vol,nx,ny,nz,bx,by,bz,compressed,compressed_length);
		clock_gettime(CLOCK_REALTIME,&after);
		double elapsed = (double)after.tv_sec + (double)after.tv_nsec * 1e-9 - (double)before.tv_sec - (double)before.tv_nsec * 1e-9;
		double mcells_per_sec = (double)nx * (double)ny * (double)nz / (elapsed * 1e6);
		printf("Compression throughput was %.0f MC/s - Compression ratio is %.2f:1\n",mcells_per_sec,ratio);

		printf("Decompressing.\n");
		int nx2,ny2,nz2;
		clock_gettime(CLOCK_REALTIME,&before);
		float* vol2 = compressor->Decompress(nx2,ny2,nz2,compressed,compressed_length);
		clock_gettime(CLOCK_REALTIME,&after);
		elapsed = (double)after.tv_sec + (double)after.tv_nsec * 1e-9 - (double)before.tv_sec - (double)before.tv_nsec * 1e-9;
		mcells_per_sec = (double)nx2 * (double)ny2 * (double)nz2 / (elapsed * 1e6);
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
	free(vol);
	return 0;
}
