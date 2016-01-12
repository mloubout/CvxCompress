#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void PrintUsage(const char* cmd)
{
	printf("Usage: %s <nx> <ny> <nz> <output-file>\n",cmd);
}

int main(int argc, char* argv[])
{
	if (argc != 5)
	{
		PrintUsage(argv[0]);
		return -1;
	}

	int nx = atoi(argv[1]);
	int ny = atoi(argv[2]);
	int nz = atoi(argv[3]);
	if (nx <= 0 || ny <= 0 || nz <= 0)
	{
		printf("Bad volume dimensions %d by %d by %d.\n",nx,ny,nz);
		PrintUsage(argv[0]);
		return -2;
	}

	FILE* fp = fopen(argv[4], "wb");
	if (fp == 0L)
	{
		printf("Unable to open %s for writing.\n",argv[4]);
		PrintUsage(argv[0]);
		return -3;
	}

	fwrite(&nx,sizeof(int),1,fp);
	fwrite(&ny,sizeof(int),1,fp);
	fwrite(&nz,sizeof(int),1,fp);
	float* arr = new float[nx];
	//srand48(time(0L));
	for (int i = 0;  i < nx;  ++i) arr[i] = 0.0f; //drand48();
	for (int i = 0;  i < ny*nz;  ++i) fwrite(arr,sizeof(float),nx,fp);

	fclose(fp);
	return 0;
}
