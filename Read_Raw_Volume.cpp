#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

/*!
* Read uncompressed volume from a binary file.
* 2024.10.27: hijacked this method to replace reading a binary with generating a volume from sinusoids + random
*/
void Read_Raw_Volume(const char* filename, int& nx, int& ny, int& nz, float*& vol) {
        nx = 151;
        ny = 101;
        nz =  51;
        size_t nn = (size_t)nx * (size_t)ny * (size_t)nz;
        posix_memalign((void**)&vol, 64, nn * (size_t)sizeof(float));
        int x0 = (nx-1) / 2;
        int y0 = (ny-1) / 2;
        int z0 = (nz-1) / 2;

#pragma omp parallel for schedule(static,1)
        for (long iz = 0;  iz < nz;  ++iz) 
        {
                memset((void*)(vol+iz*(long)nx*(long)ny),0,(long)sizeof(float)*(long)nx*(long)ny);
        }

#pragma omp parallel for schedule(static,1)
        for (long iz = 0;  iz < nz;  ++iz) 
        {
                for (long iy = 0;  iy < ny;  ++iy) 
                {
                        for (long ix = 0;  ix < nx;  ++ix) 
                        {
                                double x = ix - x0;
                                double y = iy - y0;
                                double z = iz - z0;
                                double r = sqrt(x*x + y*y + z*z);
                                vol[iz*nx*ny + ix*ny + iy] = sin(r / 10) + drand48() / 100;
                        }
                }
        }
}
