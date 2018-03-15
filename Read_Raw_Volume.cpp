#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

/*!
 * Read uncompressed volume from a binary file.
 *
 */
void Read_Raw_Volume(const char* filename, int& nx, int& ny, int& nz, float*& vol)
{
        FILE* fp = fopen(filename, "rb");
        if (fp != 0L)
        {
                fread(&nx,sizeof(int),1,fp);
                fread(&ny,sizeof(int),1,fp);
                fread(&nz,sizeof(int),1,fp);
                printf("nx=%d, ny=%d, nz=%d\n",nx,ny,nz);
                size_t nn = (size_t)nx * (size_t)ny * (size_t)nz;
                posix_memalign((void**)&vol, 64, nn * (size_t)sizeof(float));
#pragma omp parallel for schedule(static,1)
		for (long iz = 0;  iz < nz;  ++iz)
		{
			memset((void*)(vol+iz*(long)nx*(long)ny),0,(long)sizeof(float)*(long)nx*(long)ny);
		}
                fread(vol,sizeof(float),nn,fp);
                fclose(fp);
        }
	else
	{
		printf("Error! Unable to open file %s for reading.\nAborting\n",filename);
		exit(-1);
	}
}

