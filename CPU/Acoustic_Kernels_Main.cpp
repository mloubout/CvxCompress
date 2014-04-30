/*
** ISO, VTI and TTI with variable density and inverse Q.
** Codes provided for SSE and AVX equipped CPUs.
** 12/06/2012.
*/

/****************************************************************************
  TTI3MOD(): CALCULATE TTI pseudo-PRESSURE WAVEFIELD IN 3D VIA FINITE DIFFERENCE
             4th order time, 10th order space, variable density, C44>0
             OUTPUT directly to SEGY
   STAGGERED DIFFERENCES THROUGHOUT
   questions: joestefani@chevron.com   Mar 2010/2011

   Differential Equations System, uses 2 dynamic variables P and Q:
   (subscripts refer to differentiation with either t(ime), x, y, z;
    C11,C33,C13,C44 are the earth model parameters) density included

   Ptt =      C11* [{G1}(B{g1}P)+{G2}(B{g2}P)] + (C13+C44)*{G3}(B{g3}Q) + C44* {G3}(B{g3}P)
   Qtt = (C13+C44)*[{G1}(B{g1}P)+{G2}(B{g2}P)] +      C33 *{G3}(B{g3}Q) + C44*[{G1}(B{g1}Q)+{G2}(B{g2}Q)]

   where  {Gi}(B{gi}(*)) ==> d/dXm[Rim B Rin d(*)/dXn]  (No sum on i), with B = 1/Dn = Buoyancy
   and R = 3D rotation matrix:
       ( cosdip*cosazm  cosdip*sinazm  -sindip )
   R = (       -sinazm         cosazm     0    )
       ( sindip*cosazm  sindip*sinazm   cosdip )

   In this FD implementation, there are 13 3D volumes:

   2 compressed table-lookup earth parameter 3D volumes: VelAnis & DenAng
   instead of lookup can calc some parms directly
   VelAnis: CDDDDDDDDEEEEEEEEVVVVVVVVVVVVVVV
   Variable Nbits shift mask TableEntries
   Vp       15      0  32767   32768
   Eps       8     15    255     256
   Del       8     23    255     256
   C44/C33   1     31      1       2

   DenAng:  ........DDDDDDDDAAAAAAAADDDDDDDD
   Variable Nbits shift mask TableEntries
   Dip        8     0    255     256
   Azm        8     8    255     256
   Den        8    16    255     256

   4 dynamic pressure 3D volumes p,q and r,s;
      p & q correspond to the P and Q above at the current time step t, and
      r & s are the p,q counterparts at t-1 (and also t+1)

   7 intermediate variable 3D volumes V1, V2, V3, V4, V5; and  Ap, Aq for OT4;
  
   Input Parameters:
   seisname:    root name of seismic output files
   id:          identifying integer tag for this shot
   vpname epsetaname etaflag deltaname: filenames of vpvertical, epseta, epsetaflag and delta 3D volumes
   (need to be in ZslowYmediXfast order, and with correct byte endian for the host machine)
   VsoVp0: nominal Vs/Vp0 ratio (goes to 0 in isotropic)
   dnname, dipdxname, azmdyname, dipxdipyflag, degreesflag: filenames of density, dipdipx and azmdipy angles & dipflag
   swapin: 1 to swap big/little endian input
   Velspecial 0,1,2, Dnspecial 0, 1, 2   SedDenmin: avoids the gap in compression
   smoothanisflag: (=1) smooths anisotropy & dip (to avoid sharp anis boundaries)
   isorad: (about 15) puts an isotropic zone around source
   absorbz0:    0=normal reflecting free surface, 1=absorbing top surface
   srcghost:    1=turn on effective src ghosting when absorbz0=1
   recghost:    1=turn on effective rec ghosting when absorbz0=1
   dh dz:       grid spacing in horizontal and vertical direction, respectively
   stretchfacz: multiplicative factor to stretch z axis per cell (~1.003)
   nx ny nz:    number of cells in x, y, z, respectively
   srcx srcy srcz: floating pt values of source location in x,y,z
   sourcetype:  0=put source on nearest grid point, 1=interpolate source location
   (sourcetype=1 will put an effective homogeneous halo around the source location
   one wavelength in radius)
   xrecstart xrecend xrecstride: integer grid start, end & stride in X
   yrecstart yrecend yrecstride: integer grid start, end & stride in Y
   zrecstart zrecend zrecstride: integer grid start, end & stride in Z
   maxtime:     maximum recording time (seconds)
   timestartrec:time at which to begin recording
   dtout:       output time sample interval (seconds)
   newfmax:     = 0.: use default maximum frequency
                > 1.: resets maximum frequency to  newfmax
                < 0.: scales maximum frequency to |newfmax|*default
   gamfac: multiplier of Courant number gamma (typically <= 1)
   OTflag: 2=OT2, 4=OT4
   movieframestride: if>0: internal time stride between movie frames

*/

//#define TMJ_PROFILE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <unistd.h>

//#include "tti3mod_propagator.hh"

#define min(a,b) (a < b ? a : b)

//
// This is a Visual Studio 2010 consession.
// Passing __m128 args by value does not work correctly.
//
__m128 A5h_A4h;
__m128 A4h_A3h;
__m128 A3h_A2h;
__m128 A2h_A1h;
__m128 A5h_mA5h;
__m128 mA1h_A1h;

__m128 A5h_A4h_A3h_A2h;
__m128 A1h_A1h_A1h_A1h;

__m128 A5_A4_A3_A2;
__m128 A1_A1_A1_A1;

float A1h, A2h, A3h, A4h, A5h;
float B1h, B2h, B3h, B4h;
float C1h, C2h, C3h;
float D1h, D2h;
float E1h;

float A1, A2, A3, A4, A5;
float B1, B2, B3, B4;
float C1, C2, C3;
float D1, D2;
float E1;

__m128 A2h_A1h_zzz_zzz_1st_11pt;
__m128 zzz_A5h_A4h_A3h_1st_11pt;
__m128 A2h_A1h_A0h_zzz_2nd_11pt;
__m128 zzz_A5h_A4h_A3h_2nd_11pt;

float A1h_1st_11pt, A2h_1st_11pt, A3h_1st_11pt, A4h_1st_11pt, A5h_1st_11pt;
float A0h_2nd_11pt, A1h_2nd_11pt, A2h_2nd_11pt, A3h_2nd_11pt, A4h_2nd_11pt, A5h_2nd_11pt;

__m128 A2_A1_zz_zz_1st_11pt;
__m128 zz_A5_A4_A3_1st_11pt;
__m128 A2_A1_A0_zz_2nd_11pt;
__m128 zz_A5_A4_A3_2nd_11pt;

float A1_1st_11pt, A2_1st_11pt, A3_1st_11pt, A4_1st_11pt, A5_1st_11pt;
float A0_2nd_11pt, A1_2nd_11pt, A2_2nd_11pt, A3_2nd_11pt, A4_2nd_11pt, A5_2nd_11pt;

float *lutVp2, *lutEps, *lutDel, *lutDen, *lutBuoy, *lutsDip, *lutcDip, *lutsAzm, *lutcAzm, *lutQ, *lutc44c33;

float _mm_Velmin;
float _mm_Velbinsize;
float _mm_dt2;
float _mm_Epsmin;
float _mm_Epsbinsize;
float _mm_Delmin;
float _mm_Delbinsize;
float _mm_Denmin;
float _mm_Denbinsize;
float _mm_Qmin;
float _mm_Qbinsize;
float _mm_C44C33min;
float _mm_C44C33binsize;

void Init_Earth_Model_Compression(
	float Velmin,
	float Velbinsize,
	float dt2,
	float Epsmin,
	float Epsbinsize,
	float Delmin,
	float Delbinsize,
	float Denmin,
	float Denbinsize,
	float Dipmin,
	float Dipbinsize,
	float Azmmin,
	float Azmbinsize,
	float Qmin,
	float Qbinsize,
	float C44C33min,
	float C44C33binsize
	)
{
	_mm_Velmin = Velmin;
	_mm_Velbinsize = Velbinsize;
	_mm_dt2 = dt2;
	_mm_Epsmin = Epsmin;
	_mm_Epsbinsize = Epsbinsize;
	_mm_Delmin = Delmin;
	_mm_Delbinsize = Delbinsize;
	_mm_Denmin = Denmin;
	_mm_Denbinsize = Denbinsize;
	_mm_Qmin = Qmin;
	_mm_Qbinsize = Qbinsize;
	_mm_C44C33min = C44C33min;
	_mm_C44C33binsize = C44C33binsize;
}

extern void ISODenQ_TimeStep(
        int logLevel,
        int stage,
        __m128* pq,
        __m128* rs,
        __m128* Apq,
        int* VelAnis,
        float* spgx,
        float* spgy,
        float* spgz,
        int dimx,
        int dimy,
        int dimz,
        int xh,
        int yh,
        int zh,
        int bsX,
        int bsY,
        int x0,
        int x1,
        int y0,
        int y1,
        int z0,
        int z1
        );

extern void VTIDenQ_TimeStep(
        int logLevel,
        int stage,
        __m128* pq,
        __m128* rs,
        __m128* Apq,
        int* VelAnis,
        int* DenAng,
        float* spgx,
        float* spgy,
        float* spgz,
        int dimx,
        int dimy,
        int dimz,
        int xh,
        int yh,
        int zh,
        int bsX,
        int bsY,
        int x0,
        int x1,
        int y0,
        int y1,
        int z0,
        int z1
        );

extern void TTIDenQ_TimeStep(
        int logLevel,
        int stage,
        __m128* pq,
        __m128* rs,
        __m128* Apq,
        int* VelAnis,
        int* DenAng,
        float* spgx,
        float* spgy,
        float* spgz,
        int dimx,
        int dimy,
        int dimz,
        int xh,
        int yh,
        int zh,
        int bsX,
        int bsY,
        int x0,
        int x1,
        int y0,
        int y1,
        int z0,
        int z1
        );

extern void AVX_VTIDenQ_TimeStep(
        int logLevel,
        int stage,
        __m128* pq,
        __m128* rs,
        __m128* Apq,
        int* VelAnis,
        int* DenAng,
        float* spgx,
        float* spgy,
        float* spgz,
        int dimx,
        int dimy,
        int dimz,
        int xh,
        int yh,
        int zh,
        int bsX,
        int bsY,
        int x0,
        int x1,
        int y0,
        int y1,
        int z0,
        int z1
        );

extern void AVX_TTIDenQ_TimeStep(
        int logLevel,
        int stage,
        __m128* pq,
        __m128* rs,
        __m128* Apq,
        int* VelAnis,
        int* DenAng,
        float* spgx,
        float* spgy,
        float* spgz,
        int dimx,
        int dimy,
        int dimz,
        int xh,
        int yh,
        int zh,
        int bsX,
        int bsY,
        int x0,
        int x1,
        int y0,
        int y1,
        int z0,
        int z1
        );

// everything within anonymous namespace has file scope only.
namespace {

/* global variables */
int nx,ny,nz, nx1,ny1,nz1, halfOperLen, halfOperLen1;
int xstart,xend, ystart,yend, zend;
int xsrc, ysrc, zsrc, absorbz0, srcghost, recghost, zsrcghost;
int xrec, xrecstart, xrecend, xrecstride, nxrec;
int yrec, yrecstart, yrecend, yrecstride;
int zrec, zrecstart, zrecend, zrecstride, zrecstartgrid;
long shotsamp=0;

float dh, dz, stretchfacz;
int ***DenAng, ***VelAnis;
float ***V1, ***V2, ***V3;
#ifndef STANDALONE
float *recvec;
#endif
float *recsponge;

#include "Acoustic_Kernels_Common.cpp"
#include "Acoustic_Kernels_Util.cpp"

#define swapout 1
#define zero 1e-20f
/* spongecoeff can vary between 0.00005 to 0.00010 for width of 30 */
//#define spongecoeff 0.0003f  // optimal for spongewidth 30
// #define spongecoeff 0.000075f
#define r2d 57.29578f
#define d2r 0.0174532925f
#define MICRO 1.e-6f
#define OperLen 10

//
// new propagator
//

int xh = 12;
int yh = 9;
int zh = 0;

int dimxh, dimyh, dimzh;

float *spgx=0L, *spgy=0L, *spgz=0L;

int *PadDenAng = 0L, *PadVelAnis = 0L;

__m128* pq = 0L;
__m128* rs = 0L;
__m128* Apq = 0L;

int initDone = 0;

// Determine number of physical CPU cores.
// Hyper-threaded logical cores are not counted.
// Cache_Size is per core and in whole KB.
int Get_Physical_Core_Count(int& Cache_Size_Per_Core_KB)
{
	FILE* fp = fopen("/proc/cpuinfo","r");
	if (fp != 0L)
	{
		char buf[256];
		int max_phys_id = -1, max_cores = -1, max_cache_size = -1;
		while (fgets(buf, 256, fp) != 0L)
		{
			int phys_id = -1;
			if (sscanf(buf, "physical id : %d", &phys_id) == 1)
			{
				if (phys_id > max_phys_id) max_phys_id = phys_id;
			}
			
			int cpu_cores = -1;
			if (sscanf(buf, "cpu cores : %d", &cpu_cores) == 1)
			{
				if (cpu_cores > max_cores) max_cores = cpu_cores;
			}

			int cache_size = -1;
			if (sscanf(buf, "cache size : %d", &cache_size) == 1)
			{
				if (cache_size > max_cache_size) max_cache_size = cache_size;
			}
		}
		fclose(fp);
		if (max_phys_id >= 0 && max_cores > 0)
		{
			if (max_cache_size >= 0)
                	{
				Cache_Size_Per_Core_KB = max_cache_size / max_cores;
			}
			else
			{
				Cache_Size_Per_Core_KB = -1;
			}
			return (max_phys_id+1) * max_cores;
		}
		else
		{	
			return -1;
		}
	}
	return -1;
}

//
// TILE SHAPES USED FOR OPTIMAL TILE ADJUSTMENT
//

int num_tile_shapes = 0;
int* pos = 0L;
int* tile_shape = 0L;
float* avgPerf = 0L;

void Generate_Tile_Shapes(
	int logLevel,
	int is_dry_run,
	int KernelType,
	int Min_X,
	int Max_X,
	int Inc_X,
	int Min_Y,
	int Max_Y,
	int Inc_Y,
	int num_threads,
	unsigned long Avail_Cache_KB
	)
{
	if (is_dry_run == 3)
	{
		int max_num_tile_shapes = 3;
		pos = (int*)malloc(max_num_tile_shapes*sizeof(int));
		tile_shape = (int*)malloc(2*max_num_tile_shapes*sizeof(int));
		avgPerf = (float*)malloc(max_num_tile_shapes*sizeof(float));

		for (int tsCnt = 0;  tsCnt < max_num_tile_shapes;  ++tsCnt)
		{
			pos[tsCnt] = tsCnt;
			avgPerf[tsCnt] = 0.0f;
			tile_shape[tsCnt*2] = 32 + tsCnt*8;
			tile_shape[tsCnt*2+1] = 32 + tsCnt*8;
		}

		num_tile_shapes = max_num_tile_shapes;
	}
	else
	{
		if (logLevel >= 4) printf("Min_X = %d, Max_X = %d, Min_Y = %d, Max_Y = %d\n",Min_X,Max_X,Min_Y,Max_Y);

		int max_num_tile_shapes = (((Max_X-Min_X+Inc_X-1)/Inc_X)+1) * (((Max_Y-Min_Y+Inc_Y-1)/Inc_Y)+1);
		pos = (int*)malloc(max_num_tile_shapes*sizeof(int));
		tile_shape = (int*)malloc(2*max_num_tile_shapes*sizeof(int));
		avgPerf = (float*)malloc(max_num_tile_shapes*sizeof(float));

		// first eliminate tile shapes that won't fit into cache and tile shapes where bsX < bsY
		int tsCnt = 0;
		unsigned long Min_Avail_Cache_KB = 25 * Avail_Cache_KB / 100;  // don't use tile shapes with < 25% cache usage
		unsigned long Max_Avail_Cache_KB = 88 * Avail_Cache_KB / 100;  // don't use tile shapes with > 88% cache usage
		for (int X = Min_X;  X <= Max_X;  X += Inc_X)
		{
			for (int Y = Min_Y;  Y <= Max_Y;  Y += Inc_Y)
			{
				if (X >= Y)
				{
					unsigned long Cache_Usage_Per_Thread_KB = 0;
					if (KernelType == 0)
					{
						ISODenQ_Comp_Cache_Usage(X,Y,Cache_Usage_Per_Thread_KB);
					}
					else if (KernelType == 1)
					{
						VTIDenQ_Comp_Cache_Usage(X,Y,Cache_Usage_Per_Thread_KB);
					}
					else if (KernelType == 2)
					{
						TTIDenQ_Comp_Cache_Usage(X,Y,Cache_Usage_Per_Thread_KB);
					}
					unsigned long Total_Cache_Usage_KB = Cache_Usage_Per_Thread_KB * (unsigned long)num_threads;
					if (Total_Cache_Usage_KB >= Min_Avail_Cache_KB && Total_Cache_Usage_KB <= Max_Avail_Cache_KB)
					{
						double cache_usage_percent = 100.0 * (double)Total_Cache_Usage_KB / (double)Avail_Cache_KB;
						if (logLevel >= 4)
						{
							printf("%3d x %3d :: Total cache usage %5ld / %5ld KB (%.0f%%).\n",X,Y,Total_Cache_Usage_KB,Avail_Cache_KB,cache_usage_percent);
						}
						pos[tsCnt] = tsCnt;
						int ii = tsCnt << 1;
						tile_shape[ii] = X;
						tile_shape[ii+1] = Y;
						avgPerf[tsCnt] = 0.0f;
						++tsCnt;
					}
				}
			}
		}
		num_tile_shapes = tsCnt;
	}
	if (logLevel >= 3) printf("%d tile shapes generated.\n",num_tile_shapes);
}

void Sort_TileSize_Set(int* pos, float* avgPerf, int CurrSize)
{
	// bubble sort according to performance.
	int nSwap;
	do
	{
		nSwap = 0;
		for (int i = 0;  i < CurrSize-1;  ++i)
		{
			if (avgPerf[pos[i]] < avgPerf[pos[i+1]])
			{
				int tmp = pos[i];
				pos[i] = pos[i+1];
				pos[i+1] = tmp;
				++nSwap;
			}
		}
	} while (nSwap > 0);
}

void Get_TileSize(int bsIdx, int& bsX, int& bsY)
{
	if (bsIdx < num_tile_shapes)
	{
		int ii = bsIdx << 1;
		bsX = tile_shape[ii];	
		bsY = tile_shape[ii+1];
	}
	else
	{
		printf("ERROR! :: Get_TileSize, bsIdx argument of %d is invalid!\n",bsIdx);
		exit(-1);
	}
}

void omp_memset(__m128* p, unsigned long wflen)
{
	unsigned long wflen_32 = wflen >> 5;
#pragma omp parallel for
	for (int i = 0;  i < 32;  ++i)
	{
		unsigned long pos = wflen_32 * (unsigned long)i;
		memset((void*)(((char*)p)+pos), 0, wflen_32);
	}
}

void omp_setconst(
	int KernelType,
	__m128* pq, 
	unsigned long dimxh,
	unsigned long dimyh,
	unsigned long dimzh,
	float val
	)
{
	unsigned long stride_z_m128 = dimyh * (dimxh >> 1);
	if (KernelType == 0) stride_z_m128 = stride_z_m128 >> 1;
	__m128 vecval = _mm_set1_ps(val);
#pragma omp parallel for schedule(dynamic)
	for (unsigned long iZ = 0;  iZ < dimzh;  ++iZ)
	{
		__m128* curr_pq = pq + stride_z_m128 * iZ;
		for (unsigned long ii = 0;  ii < stride_z_m128;  ++ii)
		{
			curr_pq[ii] = vecval;
		}
	}
}

unsigned long Allocate_Propagation_Fields(
	int logLevel,
	int is_dry_run,
	int OTflag,
	int KernelType,
	__m128*& pq,
	__m128*& rs,
	__m128*& Apq,
	int dimx,
	int dimy,
	int dimz,
	int xh,
	int yh,
	int zh,
	float***& V1,
	float***& V2,
	float***& V3
	)
{
	unsigned long max_allocated_memory = 0L;

	dimxh = dimx + 2*xh;
	dimyh = dimy + 2*yh;
	dimzh = dimz + 2*zh;

	unsigned long cell_size = KernelType == 0 ? 4 : 8;

	unsigned long wflen = (unsigned long)dimxh * (unsigned long)dimyh * (unsigned long)dimzh * cell_size;
	if (logLevel >= 4) printf("dimx=%d, dimy=%d, dimz=%d, xh=%d, yh=%d, zh=%d, wflen = %ld bytes\n",dimx,dimy,dimz,xh,yh,zh,wflen);

	max_allocated_memory += wflen;
	if (is_dry_run == 0 || is_dry_run == 3)
	{
		posix_memalign((void**)&pq, 64, wflen);
		if (pq != 0L)
		{
			if (logLevel >= 4) printf("Allocated %ld bytes for pq\n",wflen);
			omp_memset(pq, wflen);
		}
		else
		{
			printf("FAILED to allocate pq!\n");
			exit(-1);
		}
	}

	max_allocated_memory += wflen;
	if (is_dry_run == 0 || is_dry_run == 3)
	{
		posix_memalign((void**)&rs, 64, wflen);
		if (rs != 0L)
		{
			if (logLevel >= 4) printf("Allocated %ld bytes for rs\n",wflen);
			omp_memset(rs, wflen);
		}
		else
		{
			printf("FAILED to allocate rs!\n");
			exit(-1);
		}
	}

	if (OTflag == 4)
	{
		max_allocated_memory += wflen;
		if (is_dry_run == 0 || is_dry_run == 3)
		{
			posix_memalign((void**)&Apq, 64, wflen);
			if (Apq != 0L)
			{
				if (logLevel >= 4) printf("Allocated %ld bytes for Apq\n",wflen);
				omp_memset(Apq, wflen);
			}
			else
			{
				printf("FAILED to allocate Apq!\n");
				exit(-1);
			}
		}
	}
	else
	{
		Apq = 0L;
	}

	// allocate old structures on top of new padded blocks
	unsigned long offset = (unsigned long)dimxh * (unsigned long)dimyh * (unsigned long)dimzh;
	max_allocated_memory += 3 * (unsigned long)dimz * (unsigned long)sizeof(float**);
	if (is_dry_run == 0 || is_dry_run == 3)
	{
		V1 = (float***)malloc(dimz*sizeof(float**));
		V2 = (float***)malloc(dimz*sizeof(float**));
		if (KernelType > 0)
		{
			V3 = (float***)malloc(dimz*sizeof(float**));
		}
		else
		{
			V3 = 0L;
		}
	}
	for (int iZ = 0;  iZ < dimz;  ++iZ)
	{
		max_allocated_memory += 3 * (unsigned long)dimy * (unsigned long)sizeof(float*);
		if (is_dry_run == 0 || is_dry_run == 3)
		{
			V1[iZ] = (float**)malloc(dimy*sizeof(float*));
			V2[iZ] = (float**)malloc(dimy*sizeof(float*));
			if (KernelType > 0)
			{
				V3[iZ] = (float**)malloc(dimy*sizeof(float*));
			}
			for (int iY = 0;  iY < dimy;  ++iY)
			{
				unsigned long idx = (unsigned long)dimxh * ( (unsigned long)dimyh * (unsigned long)(iZ + zh) + (unsigned long)(iY + yh) ) + (unsigned long)xh;
				V1[iZ][iY] = ((float*)pq) + idx;
				if (KernelType > 0)
				{
					V2[iZ][iY] = ((float*)pq) + offset + idx;
					V3[iZ][iY] = ((float*)rs) + idx;
				}
				else
				{
					V2[iZ][iY] = ((float*)rs) + idx;
				}
			}
		}
	}

	return max_allocated_memory;
}

void Wipe_Propagation_Wavefields(
	int KernelType,
	__m128*& pq,
	__m128*& rs,
	__m128*& Apq,
	int dimx,
	int dimy,
	int dimz,
	int xh,
	int yh,
	int zh
	)
{
	dimxh = dimx + 2*xh;
	dimyh = dimy + 2*yh;
	dimzh = dimz + 2*zh;
	
	unsigned long wflen = (unsigned long)dimxh * (unsigned long)dimyh * (unsigned long)dimzh * (unsigned long)(KernelType == 0 ? 4 : 8);
	if (pq != 0L) omp_memset(pq, wflen);
	if (rs != 0L) omp_memset(rs, wflen);
	if (Apq != 0L) omp_memset(Apq, wflen);
}

unsigned long Allocate_Padded_Earth_Model(
	int logLevel,
	int is_dry_run,
	int*& PadDenAng,
	int nx,
	int ny,
	int nz,
	int xh,
	int yh,
	int zh,
	int***& DenAng
	)
{
	unsigned long max_allocated_memory = 0;

	dimxh = nx + 2*xh;
	dimyh = ny + 2*yh;
	dimzh = nz + 2*zh;
	unsigned long wflen = (unsigned long)dimxh * (unsigned long)dimyh * (unsigned long)dimzh * (unsigned long)4;
	max_allocated_memory += wflen;
	if (is_dry_run == 0 || is_dry_run == 3)
	{
		posix_memalign((void**)&PadDenAng, 64, wflen);
		if (PadDenAng != 0L)
		{
			if (logLevel >= 4) printf("Allocated %ld bytes for PadDenAng!\n",wflen);
			//memset((void*)PadDenAng, 0, wflen);
			omp_memset((__m128*)PadDenAng, wflen);
		}
		else
		{
			printf("FAILED to allocate PadDenAng!\n");
			exit(-1);
		}
	}

	// map old structures onto new padded blocks
	max_allocated_memory += (unsigned long)nz * (unsigned long)sizeof(int**);
	if (is_dry_run == 0 || is_dry_run == 3) DenAng = (int***)malloc(nz*sizeof(int**));
	for (int iZ = 0;  iZ < nz;  ++iZ)
	{
		max_allocated_memory += (unsigned long)ny * (unsigned long)sizeof(int*);
		if (is_dry_run == 0 || is_dry_run == 3)
		{
			 DenAng[iZ] = (int**)malloc(ny*sizeof(int*));
			 for (int iY = 0;  iY < ny;  ++iY)
			 {
				 unsigned long idx = (unsigned long)dimxh * ( (unsigned long)dimyh * (unsigned long)(iZ + zh) + (unsigned long)(iY + yh) ) + (unsigned long)xh;
				 DenAng[iZ][iY] = PadDenAng + idx; 
			 }
		}
	}

	return max_allocated_memory;
}

void Free_Padded_Earth_Model(
	int*& PadDenAng,
	int nz,
	int***& DenAng
	)
{
	if (DenAng != 0L)
	{
		for (int iZ = 0;  iZ < nz;  ++iZ)
		{
			free((void*)DenAng[iZ]);
		}
		free((void*)DenAng);
	}
	if (PadDenAng != 0L) free((void*)PadDenAng);
}

inline void AddToPQ(
	int KernelType,
	__m128* pq,
	int iX,
	int iY,
	int iZ,
	float val_P,
	float val_Q
	)
{
	if (KernelType > 0)
	{
		unsigned long idx = ( ( (unsigned long)(iZ+zh) * (unsigned long)dimyh + (unsigned long)(iY+yh) ) * (unsigned long)dimxh + (unsigned long)(iX+xh) ) << 1;
		float* fpq = (float*)pq;
		fpq[idx] += val_P;
		fpq[idx+1] += val_Q;
	}
	else
	{
		unsigned long idx = ( ( (unsigned long)(iZ+zh) * (unsigned long)dimyh + (unsigned long)(iY+yh) ) * (unsigned long)dimxh + (unsigned long)(iX+xh) );
		float* fpq = (float*)pq;
		fpq[idx] += val_P;
	}
}

inline float GetP(
	int KernelType,
	__m128* pq,
	int iX,
	int iY,
	int iZ
	)
{
	unsigned long idx = ( ( (unsigned long)(iZ+zh) * (unsigned long)dimyh + (unsigned long)(iY+yh) ) * (unsigned long)dimxh + (unsigned long)(iX+xh) );
	if (KernelType > 0) idx = idx << 1;
	float* fpq = (float*)pq;
	return fpq[idx];
}

inline float GetQ(
	int KernelType,
	__m128* pq,
	int iX,
	int iY,
	int iZ
	)
{
	unsigned long idx = ( ( (unsigned long)(iZ+zh) * (unsigned long)dimyh + (unsigned long)(iY+yh) ) * (unsigned long)dimxh + (unsigned long)(iX+xh) );
	if (KernelType > 0)
	{
		// return Q
		idx = idx << 1;
		float* fpq = (float*)pq;
		return fpq[idx+1];
	}
	else
	{
		// ISO kernel, Q == P
		// return P
		float* fpq = (float*)pq;
		return fpq[idx];
	}
}

/***** ALTER_DTOUT: calculate a "round" dtout just > or < input delt *****/
/* e.g. will change delt from .00314159 to .004 or .003 dep. on direction */
float alter_dtout(int direction, float delt)
{
	int i, tencount=0;
	while(delt < 1.0f)  
	{
		delt *= 10.0f;  
		++tencount; 
	}

	if(direction == 1)        
		delt = (float)((int)(delt) + 1);
	else if(direction == -1)  
		delt = (float)((int)(delt));

	for(i=0; i<tencount; i++)  
		delt /= 10.0f;

	return(delt);
}


/***** SRC: set up source time function stf[] and related parameters. *****/

void src(int logLevel, float dt, float fmax, int type, char* stfname, int* tsrc, float* stf)
{ 
	if (logLevel >= 4)
	{
		printf("src(logLevel=%d, dt=%e, fmax=%e, type=%d, stfname=%s, tsrc=%s, stf=%s)\n",logLevel,dt,fmax,type,stfname!=0L?stfname:"nil",tsrc!=0L?"ok":"nil",stf!=0L?"ok":"nil");
		fflush(stdout);
	}
	if (type == 1)
	{
		/* type==1: first derivative of a Gaussian, with linear extension tapers */
		int i, imax, imaxhalf, ntap;
		float w0, wmax, ts,t,rt, wt, wt2, diff;

		wmax = 2.0f*M_PI*fmax;
		/* Note:   tsrc = ts/dt = 2/(gam*khmax)*(rt*rt)*(Vmax/Vmin) */

		/* if(type==1) */  /* only one type for now */
		{
			rt = 3.571625f; /* guarantees SourceAmplitude(tmax)   = 0.01*MaxAmplitude
					   and guarantees SpectralAmplitude(fmax) = 0.01*MaxSpectrum
					   and: wmax/w0 = rt = w0*ts/2  */
			w0 = wmax/rt;  /* w0i = 1./w0; */
			ts = 2.0f*rt/w0;  /* total source time */
			imax = (int)(ts/dt) + 1;

			for(i=0;i<imax;i++)
			{ t=i*dt-0.5f*ts;
				stf[i]  = -sqrtf(M_E)*w0*t*exp(-0.5f*w0*w0*t*t);
			}

			/* taper (linearly extend) front and back ends */
			/* front end */
			diff = stf[1]-stf[0];
			ntap = (int)(fabs(stf[0]/diff));
			for(i=imax-1; i>=0; i--) stf[i+ntap] = stf[i];  /* shift */
			for(i=ntap-1; i>=0; i--) stf[i] = stf[i+1] - diff; /* taper */
			imax += ntap;

			/* back end: */
			diff = stf[imax-1]-stf[imax-2];
			ntap = (int)(fabs(stf[imax-1]/diff));
			for(i=0; i<ntap; i++)  stf[imax+i] = stf[imax+i-1] + diff; /* taper */
			imax += ntap;
		}

		*tsrc = imax;

		if (logLevel >= 4) printf("SOURCE TYPE 1 : imax=%d\n",imax);

		// for(i=0; i<imax; i++) stf[i]=0.0f; stf[0] = 1.0f;

		// for(i=0; i<imax; i++) printf("%d  %f\n", i,stf[i]);
	}
	else if (type == 2)
	{
		/* type==2: user provided source function from file */
		int nfine;
		float dtfine;
		FILE* stffile = fopen(stfname,"r");
		if (stffile == 0L)
		{
			printf("ERROR! src(...) - Source wavelet file '%s' cannot be read.\n",stfname);
			fflush(stdout);
			exit(-1);
		}
		fscanf(stffile,"%d %f", &nfine, &dtfine);
		float* stffine = (float*)malloc((nfine+1)*sizeof(float));
		for(int i=0; i<nfine; i++) fscanf(stffile,"%f", &stffine[i]);
		stffine[nfine] = 0.;

		int imax = (int)((nfine-1)*dtfine/dt) + 1;
		float absmax = -1e37f;
		for(int i=0; i<imax; i++)
		{
			float t = i*dt;
			int tfine = (int)(t/dtfine);
			float frac = t/dtfine - tfine;
			float val = (1.-frac)*stffine[tfine] + frac*stffine[tfine+1];
			stf[i] = val;
			float absval = val < 0.0f ? -val : val;
			if (absval > absmax) absmax = absval;
		}
		*tsrc = imax;

		if (logLevel >= 4) printf("SOURCE TYPE 2 : nfine=%d, dtfine=%e, imax=%d, absmax=%e\n",nfine,dtfine,imax,absmax);
		
		for(int i=0; i<imax; i++) 
		{
			stf[i] /= absmax;
			if (logLevel >= 4) printf("stf[%d] = %e\n",i,stf[i]);
		}
	}
}


/***** SAVETIMESLICE: save seismic data at time interval dtout
   at several x,y locations at a common zrec or
   at several z locations at a common xrec,yrec
   Pressure = (2p+q)/3, (r+2s)/3, omit the "3";  in an isotropic medium,
   p =q and (2p+q)/3 is just = p, but for receivers in anisotropic medium,
   p!=q and need full expression *****/
#ifdef STANDALONE
void saveTimeslice(int KernelType, float tfrac, int skylayer, FILE* fp_recvec, float* write_buf)
{
	long nwr = 0;
	float* curr_recsponge = recsponge;
	for(int z=zrecstartgrid; z<=zrecend; z+=zrecstride)
	{
		for(int y=yrecstart; ((yrecstride > 0 && y<=yrecend) || (yrecstride < 0 && y>=yrecend)); y+=yrecstride)
		{ 
			int write_buf_idx = 0;
			int zrecghost = 2*skylayer - z;
			if (absorbz0 && recghost && zrecend==zrecstartgrid && zrecghost > 0)
			{
				// add receiver ghost
				for(int x=xrecstart; ((xrecstride > 0 && x<=xrecend) || (xrecstride < 0 && x>=xrecend)); x+=xrecstride, ++write_buf_idx, ++shotsamp)
				{
					float sponge_factor = *(curr_recsponge++);
					write_buf[write_buf_idx] = 
						(      tfrac * sponge_factor * (2.0f*GetP(KernelType,pq,x,y,z)         + GetQ(KernelType,pq,x,y,z))         + 
						       (1.0f-tfrac) *          (2.0f*GetP(KernelType,rs,x,y,z)         + GetQ(KernelType,rs,x,y,z)))        - 
						(      tfrac * sponge_factor * (2.0f*GetP(KernelType,pq,x,y,zrecghost) + GetQ(KernelType,pq,x,y,zrecghost)) +
                                                       (1.0f-tfrac) *          (2.0f*GetP(KernelType,rs,x,y,zrecghost) + GetQ(KernelType,rs,x,y,zrecghost)));
				}
			}
			else
			{	
				// no receiver ghost
				for(int x=xrecstart; ((xrecstride > 0 && x<=xrecend) || (xrecstride < 0 && x>=xrecend)); x+=xrecstride, ++write_buf_idx, ++shotsamp)
				{
					float sponge_factor = *(curr_recsponge++);
					write_buf[write_buf_idx] = 
						(      tfrac * sponge_factor * (2.0f*GetP(KernelType,pq,x,y,z) + GetQ(KernelType,pq,x,y,z)) + 
						       (1.0f-tfrac) *          (2.0f*GetP(KernelType,rs,x,y,z) + GetQ(KernelType,rs,x,y,z)));
				}
			}
			fwrite(write_buf, sizeof(float), write_buf_idx, fp_recvec);
			nwr += (long)sizeof(float) * (long)write_buf_idx;
		}
	}
//	printf("Wrote %ld bytes\n",nwr);
}
#else
void saveTimeslice(int KernelType, float tfrac, int skylayer)
{
	float* curr_recsponge = recsponge;
	for(int z=zrecstartgrid; z<=zrecend; z+=zrecstride)
	{
		for(int y=yrecstart; ((yrecstride > 0 && y<=yrecend) || (yrecstride < 0 && y>=yrecend)); y+=yrecstride)
		{ 
			for(int x=xrecstart; ((xrecstride > 0 && x<=xrecend) || (xrecstride < 0 && x>=xrecend)); x+=xrecstride, shotsamp++)
			{
				float sponge_factor = *(curr_recsponge++);
				recvec[shotsamp] = 
					(      tfrac * sponge_factor * (2.0f*GetP(KernelType,pq,x,y,z) + GetQ(KernelType,pq,x,y,z)) + 
					(1.0f-tfrac) *                 (2.0f*GetP(KernelType,rs,x,y,z) + GetQ(KernelType,rs,x,y,z)));
			}
			//    recvec[shotsamp] =       tfrac*(q[z][y][x]) +
			//                       (1.0f-tfrac)*(s[z][y][x]);

			if(absorbz0 && recghost && zrecend==zrecstartgrid) // excludes receiver ghost in vsp 
			{ 
				int zrecghost = 2*skylayer - z;
				if(zrecghost > 0)    // back up nx and superpose
				{ 
					int x;
					for(x=xrecstart, shotsamp -= nxrec; ((xrecstride > 0 && x<=xrecend) || (xrecstride < 0 && x>=xrecend)); x+=xrecstride, shotsamp++)
					{
						float sponge_factor = *(curr_recsponge++);
						recvec[shotsamp] -= 
							(      tfrac * sponge_factor * (2.0f*GetP(KernelType,pq,x,y,zrecghost) + GetQ(KernelType,pq,x,y,zrecghost)) +
							(1.0f-tfrac) *                 (2.0f*GetP(KernelType,rs,x,y,zrecghost) + GetQ(KernelType,rs,x,y,zrecghost)));
					}
				}
			}
		}
	}
}
#endif

void saveSponge(int KernelType, __m128* sponge_3D, int skylayer)
{
	float* curr_recsponge = recsponge;
	for(int z=zrecstartgrid; z<=zrecend; z+=zrecstride)
	{
		for(int y=yrecstart; ((yrecstride > 0 && y<=yrecend) || (yrecstride < 0 && y>=yrecend)); y+=yrecstride)
		{ 
			for(int x=xrecstart; ((xrecstride > 0 && x<=xrecend) || (xrecstride < 0 && x>=xrecend)); x+=xrecstride)
			{
				*(curr_recsponge++) = GetP(KernelType,sponge_3D,x,y,z);
			}

			if(absorbz0 && recghost && zrecend==zrecstartgrid) /* excludes receiver ghost in vsp */
			{ 
				int zrecghost = 2*skylayer - z;
				if(zrecghost > 0)    /* back up nx and superpose */
				{ 
					int x;
					for(x=xrecstart; ((xrecstride > 0 && x<=xrecend) || (xrecstride < 0 && x>=xrecend)); x+=xrecstride)
					{
						*(curr_recsponge++) = GetP(KernelType,sponge_3D,x,y,zrecghost);
					}
				}
			}
		}
	}
}

void genSEGYFilename(char* segybasename, int ffid, char* segyname)
{
	if     (ffid<10)    sprintf(segyname,"%s0000%d.segy", segybasename, ffid);
	else if(ffid<100)   sprintf(segyname,"%s000%d.segy", segybasename, ffid);
	else if(ffid<1000)  sprintf(segyname,"%s00%d.segy", segybasename, ffid);
	else if(ffid<10000) sprintf(segyname,"%s0%d.segy", segybasename, ffid);
	else                sprintf(segyname,"%s%d.segy", segybasename, ffid);
}

/***** CONVERT ENTIRE SHOT PROFILE TO SEGY & WRITE TO DISK
       Transpose 3D shot profile output from p[t][y][x] to p[y][x][t]
       (i.e., from planes of constant time to trace columns);
       Optionally swap data & header bytes from linux little endian
       (for ProMAX readability if data and/or headers need byte swapping); *****/
void writeSEGY(float* seisdata, char *segyfilename, char *hdrstring, int swapflag, int leaky_integrate_flag,
		int ffid, float srcx, float srcy, float srcz,
		float dt, float dtout, float timestartrec, int nsamp, int nrec,
		float recxstart, int nx_inline, float dx_inline,
		float recystart, int ny_crossline, float dy_crossline,
		float reczstart, int nzrec, float dzrec)
{
	/* SEGY DATA TYPES ***/
	char reel_id_hdr1[3200];
	memset((void*)reel_id_hdr1, 0, 3200);

	struct
	{ 
		int jobid;  
		int lineid;  
		int reelid;
		short ntrc_per_record; 
		short nauxtrc;
		short dtreel; 
		short dtfield;
		short nsampreel;
		short nsampfield;
		short datafmt; 
		short cmpfold;  
		short sortcode;
		char skip[370];
	} reel_id_hdr2;
	memset((void*)&reel_id_hdr2, 0, sizeof(reel_id_hdr2));

	struct
	{ 
		int trcseqno;  
		int skip0;  
		int isrc;  
		int ichan;
		int skip1;  
		int cmpbin; 
		int trcensemb; 
		short code;
		char skip3[6];
		int offset;  
		int recelev;
		int elevatsrc;
		int srcdepth; 
		char skip4[16];
		short scalar1;  
		short scalar2;
		int srcx;  
		int srcy;  
		int recx;    
		int recy;  
		short lenunit;
		char skip5[18]; 
		short tstartrec; 
		char skip6[4]; 
		short nsamp;  
		short dtmicro; 
		char skip7[82];
		float cmp_x;  
		float cmp_y; 
		int iline_no; 
		int xline_no;
		float xoff;  
		float yoff; 
		float azim;
		char skip8[12];
	} trc_id_hdr;
	memset((void*)&trc_id_hdr, 0, sizeof(trc_id_hdr));
	/* cmp_x starts at byte position 201 */

	char suffix[16];
	short one2=1, five2=5, nrec2, dtmicro2, nsamp2, fold2, trc_sortcode2, tstartrec2;
	short neg100 = -100;
	int one=1, one4=1, elevatsrc, recelev, srcdepth, xsrc, ysrc, xrec, yrec, izrec;
	int ichn, ichan, trcseq, trcens;
	int xline, iline, xinline, ycrossline, offset;
	int t, trcseqno, trcensemb, dtmicro;
	long zxyline;
	long ynx, currtrace;
	long ptr, nsegyvol;
	float recx, recy, recz, cmpx, cmpy, xoff, yoff, azim;
	float *segyvol, *trace;
	void swap2bytes(short*, int), swap4bytes(int*, int);
	FILE *segyfile;
	reel_id_hdr2.jobid = 0;
	trc_id_hdr.trcseqno = 0;

	/*** OPEN OUTPUT FILE ***/
	if( !(segyfile = fopen(segyfilename,"w")) )
	{ 
		fprintf(stderr,"Cannot open %s\n",segyfilename);
		exit(0);
	}

	/*** ALLOCATE MEM ***/
	trace = (float *)malloc(nsamp*sizeof(float));
	nsegyvol = nrec*(nsamp*sizeof(float) + 240) + 3600;  /* 240=trc header, 3600=main+reel */
	segyvol = (float *)malloc(nsegyvol);  /* Rec array */

	/*** FILL REEL ID HEADER 1 ***/
	/* Write Time of Day */
	/*status = gettimeofday(&time, &tz);
	  secs = time.tv_sec;
	  timenow = ctime(&secs);
	  strcat(reel_id_hdr1, timenow);
	 */
	strcpy(reel_id_hdr1, hdrstring);
	memcpy(segyvol, &reel_id_hdr1, 3200);
	ptr = 800;  /* 800 = number of 4-byte words in 3200 byte master header */

	/*** FILL REEL ID HEADER 2 ***/
	dtmicro = (int)(1000000.*dtout + 0.5);
	trc_sortcode2 = 1;  /* as recorded, no sorting */
	fold2 = 1;
	one2 = one;       one4 = one;
	nrec2 = nrec;     dtmicro2 = dtmicro;        nsamp2 = nsamp;

	if(swapflag)
	{
		swap2bytes(&one2, 1);        swap2bytes(&five2, 1);  swap4bytes(&one4, 1);
		swap2bytes(&nrec2, 1);       swap2bytes(&dtmicro2, 1);
		swap2bytes(&nsamp2, 1);      swap2bytes(&trc_sortcode2, 1);
		swap2bytes(&fold2, 1);
	}
	reel_id_hdr2.jobid = reel_id_hdr2.lineid = reel_id_hdr2.reelid = one4;
	reel_id_hdr2.ntrc_per_record = nrec2;
	reel_id_hdr2.dtreel = reel_id_hdr2.dtfield = dtmicro2;
	reel_id_hdr2.nsampreel = reel_id_hdr2.nsampfield = nsamp2;
	reel_id_hdr2.datafmt = five2;
	reel_id_hdr2.cmpfold = fold2;
	reel_id_hdr2.sortcode = trc_sortcode2;
	memcpy(&segyvol[ptr], &reel_id_hdr2, 400);
	ptr += 100;  /* 100 = number of 4-byte words in 400 byte master header */


	/*** FILL SOURCE-RELATED PART OF TRACE HEADER ***/
	elevatsrc = 0;  
	srcdepth = (int)(100.*srcz);
	xsrc = (int)(100.*srcx);  
	ysrc = (int)(100.*srcy);
	tstartrec2 = (int)(timestartrec*1000. + 0.5);
	if(swapflag)
	{
		swap4bytes(&ffid, 1);
		swap4bytes(&elevatsrc, 1);     swap4bytes(&srcdepth, 1);
		swap4bytes(&xsrc, 1);          swap4bytes(&ysrc, 1);
		swap2bytes(&tstartrec2, 1);
		swap2bytes(&neg100, 1);
	}
	trc_id_hdr.isrc = ffid;
	trc_id_hdr.elevatsrc = elevatsrc; trc_id_hdr.srcdepth = srcdepth;
	trc_id_hdr.srcx = xsrc;           trc_id_hdr.srcy = ysrc;
	trc_id_hdr.nsamp = nsamp2;
	trc_id_hdr.tstartrec = tstartrec2;
	trc_id_hdr.dtmicro = dtmicro2;
	trc_id_hdr.scalar1 = neg100;
	trc_id_hdr.scalar2 = neg100;


	/*** READ IN SEISMIC DATA (slow axis is Time, med axis is Y, fast axis is X)
	  AND WRITE OUT TRACE HEADER + DATA TRACE ***/
	for(izrec = 0, trcseqno=ichn=0; izrec<nzrec; izrec++)
	{
		zxyline = (long)izrec *  (long)nx_inline * (long)ny_crossline;
		recz = reczstart + izrec*dzrec;
		recelev = -(int)(100.*recz);  if(swapflag) swap4bytes(&recelev, 1);
		trc_id_hdr.recelev = recelev;

		for(ycrossline=0; ycrossline<ny_crossline; ycrossline++)
		{
			ynx = zxyline + (long)ycrossline * (long)nx_inline;
			recy = recystart + ycrossline*dy_crossline;
			yrec = (int)(100.*recy);  xline = ycrossline+1;
			if(swapflag) { swap4bytes(&yrec, 1); swap4bytes(&xline, 1); }
			trc_id_hdr.recy = yrec;
			trc_id_hdr.iline_no = xline; /* yes, this is correct */

			for(xinline=0, trcensemb=1;    xinline<nx_inline;   xinline++, trcensemb++)
			{ 
				currtrace = (long)ynx + (long)xinline;

				trcseq = trcseqno++;       ichan = ichn++;      trcens = trcensemb;
				recx = recxstart + xinline*dx_inline;
				xrec = (int)(100.*recx);
				xoff = recx - srcx;        yoff = recy - srcy;
				cmpx = 0.5*(srcx + recx);  cmpy = 0.5*(srcy + recy);
				iline = xinline+1;
				offset = (int)(sqrtf(yoff*yoff + xoff*xoff) + 0.5);
				azim = r2d*atan2f(yoff, xoff);

				if(swapflag)
				{ 
					swap4bytes(&trcseq, 1);  swap4bytes(&ichan, 1);  swap4bytes(&trcens, 1);
					swap4bytes(&xrec, 1);
					swap4bytes((int*)(&cmpx), 1); swap4bytes((int*)(&cmpy), 1);
					swap4bytes(&iline, 1);
					swap4bytes((int*)(&xoff), 1); swap4bytes((int*)(&yoff), 1);
					swap4bytes(&offset, 1);       swap4bytes((int*)(&azim), 1);
				}

				/* Assign & Write Trace Header */
				trc_id_hdr.trcseqno = trcseq;
				trc_id_hdr.ichan = ichan;
				trc_id_hdr.trcensemb = trcens;
				trc_id_hdr.offset = offset;
				trc_id_hdr.recx = xrec;
				trc_id_hdr.cmp_x = cmpx;
				trc_id_hdr.cmp_y = cmpy;
				trc_id_hdr.xline_no = iline; /* yes, this is correct */
				trc_id_hdr.xoff = xoff;
				trc_id_hdr.yoff = yoff;
				trc_id_hdr.azim = azim;
				memcpy(&segyvol[ptr], &trc_id_hdr, 240);
				ptr += 60;  /* 60 = number of 4-byte words in 240 byte trace header */

				/* Read one trace into trace[] array and swapbytes */
				for(t=0; t<nsamp; t++)  trace[t] = seisdata[currtrace + (long)t * (long)nrec];
				if (leaky_integrate_flag)
				{
					float leaky_fac = 0.98f;  // constant provided by Joe Stefani
					float Vz = 0.0f;
					for (t=1; t<nsamp; ++t)
					{
						Vz = leaky_fac * Vz + 0.5f * (trace[t-1] + trace[t]) * dtout;
						trace[t-1] = Vz;
					}
					// shift trace one sample
					for (t=nsamp-1; t>0; --t) trace[t] = trace[t-1];
				}
				if(swapflag)  swap4bytes((int*)trace,nsamp);

				/* Write Trace Vector to memory */
				memcpy(&segyvol[ptr], trace, nsamp*sizeof(float));
				ptr += nsamp;
			}
		}
	}

	/* WRITE EVERYTHING TO DISK IN ONE GO */
	fwrite(segyvol, nsegyvol, 1, segyfile);

	fclose(segyfile);

	free((void*)trace);
	free((void*)segyvol);
}


/***** swap2bytes ba --> ab *****/
void swap2bytes(short *i2, int n)
{
	int i;
	short a,b;
	for (i=0; i<n; i++)
	{ 
		a = i2[i] << 8;    
		b = (i2[i] >> 8) & 255; 
		i2[i] = a | b; 
	}
}


/***** swap4bytes:  dcba --> abcd *****/
void swap4bytes(int *i4, int n)
{
  int k, i, a, b, c, d, bmask = 16711680, cmask = 65280, dmask = 255;
  for(k=0; k<n; k++)
  { i = i4[k];
    a =  i << 24;          b = (i << 8)  & bmask;
    c = (i >> 8) & cmask;  d = (i >> 24) & dmask;
    i4[k] = a | b | c | d ;
  }
}

/* Allocate a 2D float matrix */
float **matrix2(int nslow, int nfast)
{ int z;   float **m;
  m=(float **) malloc((unsigned) nslow*sizeof(float*));
  if (!m)  { fprintf(stderr,"Memory allocation failure 1 in matrix2.\n");
             exit(0); }

  for(z=0; z<nslow; z++)
   {  m[z]=(float *) malloc((unsigned) nfast*sizeof(float));
      if (!m[z]) { fprintf(stderr,"Memory allocation failure 2 in matrix2.\n");
                   exit(0); }
   }
  return m;
}

/* Allocate a 3D float matrix */
float ***matrix3(int nslow, int nmed, int nfast)
{ int z,y;   float ***m;
  m=(float ***) malloc((unsigned) nslow*sizeof(float**));
  if (!m)  { fprintf(stderr,"Memory allocation failure 1 in matrix3.\n");
             exit(0); }

  for(z=0; z<nslow; z++)
   {  m[z]=(float **) malloc((unsigned) nmed*sizeof(float*));
      if (!m[z]) { fprintf(stderr,"Memory allocation failure 2 in matrix3.\n");
                   exit(0); }

      for(y=0; y<nmed; y++)
       {  m[z][y]=(float *) malloc((unsigned) nfast*sizeof(float));
          if (!m[z][y])
          { fprintf(stderr,"Memory allocation failure 3 in matrix3.\n");
            exit(0); }
       }
   }
  return m;
}

void free_matrix3(float*** m, int nslow, int nmed, int nfast)
{
	for (int iZ = 0;  iZ < nslow;  ++iZ)
	{
		for (int iY = 0;  iY < nmed;  ++iY)
		{
			free((void*)m[iZ][iY]);
		}
		free((void*)m[iZ]);
	}
	free((void*)m);
}

void free_imatrix3(int*** m, int nslow, int nmed, int nfast)
{
	free_matrix3((float***)m, nslow, nmed, nfast);
}

/* Allocate a 3D int matrix */
int ***imatrix3(int nslow, int nmed, int nfast)
{ int z,y;   int ***m;
  m=(int ***) malloc((unsigned) nslow*sizeof(int**));
  if (!m)  { fprintf(stderr,"Memory allocation failure 1 in imatrix3.\n");
             exit(0); }

  for(z=0; z<nslow; z++)
   {  m[z]=(int **) malloc((unsigned) nmed*sizeof(int*));
      if (!m[z]) { fprintf(stderr,"Memory allocation failure 2 in imatrix3.\n");
                   exit(0); }

      for(y=0; y<nmed; y++)
       {  m[z][y]=(int *) malloc((unsigned) nfast*sizeof(int));
          if (!m[z][y])
          { fprintf(stderr,"Memory allocation failure 3 in imatrix3.\n");
            exit(0); }
       }
   }
  return m;
}

const char* getAxisSymbol(int axis_index)
{
	switch (axis_index)
	{
		case 0:
			return "X";
		case 1:
			return "Y";
		case 2:
			return "Z";
		default:
			return "Invalid!";
	}
}

void compUVWRangesFromXYZ(
	int arg_fastAxis, int arg_medAxis, int arg_slowAxis,
	int arg_xbeg, int arg_xend, int arg_ybeg, int arg_yend, int arg_zbeg, int arg_zend,
	int arg_nx, int arg_ny, int arg_nz,
	int& dimu, int& actual_dimu, int& dimv, int& dimw,
	int& ubeg, int& uend, int& vbeg, int& vend, int& wbeg, int& wend
	)
{
	switch (arg_fastAxis)
	{
	case 0:
		dimu = arg_nx;
		ubeg = arg_xbeg;
		uend = arg_xend;
		break;
	case 1:
		dimu = arg_ny;
		ubeg = arg_ybeg;
		uend = arg_yend;
		break;
	case 2:
		dimu = arg_nz;
		ubeg = arg_zbeg;
		uend = arg_zend;
		break;
	default:
		fprintf(stderr, "ERROR! Invalid value for arg_fastAxis\n");
		exit(-1);
	}
	actual_dimu = uend - ubeg + 1;
	switch (arg_medAxis)
	{
	case 0:
		dimv = arg_nx;
		vbeg = arg_xbeg;
		vend = arg_xend;
		break;
	case 1:
		dimv = arg_ny;
		vbeg = arg_ybeg;
		vend = arg_yend;
		break;
	case 2:
		dimv = arg_nz;
		vbeg = arg_zbeg;
		vend = arg_zend;
		break;
	default:
		fprintf(stderr, "ERROR! Invalid value for arg_medAxis\n");
		exit(-1);
	}
	switch (arg_slowAxis)
	{
	case 0:
		dimw = arg_nx;
		wbeg = arg_xbeg;
		wend = arg_xend;
		break;
	case 1:
		dimw = arg_ny;
		wbeg = arg_ybeg;
		wend = arg_yend;
		break;
	case 2:
		dimw = arg_nz;
		wbeg = arg_zbeg;
		wend = arg_zend;
		break;
	default:
		fprintf(stderr, "ERROR! Invalid value for arg_slowAxis\n");
		exit(-1);
	}
}

void compXYZFromUVW(
	int arg_fastAxis, int arg_medAxis, int arg_slowAxis,
	int u, int v, int w,
	int& x, int& y, int& z
	)
{
	switch (arg_fastAxis)
	{
	case 0:
		x = u;
		break;
	case 1:
		y = u;
		break;
	case 2:
		z = u;
		break;
	default:
		fprintf(stderr, "ERROR! Bad value for U!\n");
		exit(-1);
	}
	switch (arg_medAxis)
	{
	case 0:
		x = v;
		break;
	case 1:
		y = v;
		break;
	case 2:
		z = v;
		break;
	default:
		fprintf(stderr, "ERROR! Bad value for V!\n");
		exit(-1);
	}
	switch (arg_slowAxis)
	{
	case 0:
		x = w;
		break;
	case 1:
		y = w;
		break;
	case 2:
		z = w;
		break;
	default:
		fprintf(stderr, "ERROR! Bad value for W!\n");
		exit(-1);
	}
}

float ISO_Determine_Vpmax(
	int logLevel,
	char* filename_vp,
	float scalarval_vp,
	int swapin,
	int ubeg,
	int uend,
	int vbeg,
	int vend,
	int wbeg,
	int wend,
	int dimu,
	int dimv,
	int dimw,
	int arg_fastAxis,
	int arg_medAxis,
	int arg_slowAxis
	)
{
	if (ubeg < 0 || uend >= dimu || vbeg < 0 || vend >= dimv || wbeg < 0 || wend >= dimw)
	{
		printf("ERROR! Read_Single_Earth_Model_Attribute - Sub volume extends beyond dimensions of original volume!\n");
		printf("u = [%d, %d], dimu = %d\n",ubeg,uend,dimu);
		printf("v = [%d, %d], dimv = %d\n",vbeg,vend,dimv);
		printf("w = [%d, %d], dimw = %d\n",wbeg,wend,dimw);
	}	

	// Read from file or use constant value.
	bool constflag_vp = false; 
	float constval_vp = 0.0f;
	FILE* file_vp = 0L;
	if (filename_vp == 0L)
	{
		constval_vp = scalarval_vp;
		constflag_vp = true;
		if (logLevel >= 3) printf("nil filename: Using const val = %f\n", constval_vp);
	}
	else
	{
		if( !(file_vp = fopen(filename_vp,"r")))
		{ 
			constval_vp = atof(filename_vp);  
			constflag_vp=true;
			printf("%s not found: Using const val = %f\n", filename_vp,constval_vp); 
		}
	}

	int actual_dimu = uend - ubeg + 1;

	// loop over u,v,w.
	// read, convert to right endian-ness and transpose to x-y-z ordering.
	int nthreads = 0;
#pragma omp parallel
	{
		nthreads = omp_get_num_threads();
	}
	if (nthreads < 1) nthreads = 1;  // just in case...
	float* rawinpbuf = new float[actual_dimu*nthreads*3];
	int nwpu = (wend-wbeg)/100;
	if (nwpu < 1) nwpu = 1;
	float vpmax = -1e37f;
#pragma omp parallel for schedule(dynamic)
	for (int w = wbeg;  w <= wend;  ++w)
	{
		int curr_thread = omp_get_thread_num();
		float* inpbuf_vp = rawinpbuf + actual_dimu * curr_thread * 3;
		float* inpbuf_eta = inpbuf_vp + actual_dimu;
		float* inpbuf_dta = inpbuf_eta + actual_dimu;
		if (logLevel >= 3 && ((w-wbeg) % nwpu) == 0)
		{
#pragma omp critical
			{
				printf("\r%s (%.0f%%)","Scanning Vp, Eta and Dta",100.0f*(float)(w-wbeg)/(float)(wend-wbeg));
				fflush(stdout);
			}
		}
		for (int v = vbeg;  v <= vend;  ++v)
		{
			long offset = (long)sizeof(float) * ((long)w * (long)dimu * (long)dimv + (long)v * (long)dimu + (long)ubeg);
			//fprintf(stderr, "Input offset = %ld u,v,w = [%d,%d,%d]\n", offset, ubeg, v, w);

			// read from files. we do this sequentially, since hitting the filesystem with many parallel requests most likely will trash performance.
			if (constflag_vp)
			{
				for (int u = ubeg;  u <= uend;  ++u) inpbuf_vp[u-ubeg] = constval_vp;
			}
#pragma omp critical
			{
				if (!constflag_vp)
				{
					fseek(file_vp, offset, SEEK_SET);
					fread(inpbuf_vp, sizeof(float), actual_dimu, file_vp);
				}
			}
			if(swapin)
			{
				if (!constflag_vp ) swap4bytes((int*)inpbuf_vp , actual_dimu);
			}
			for (int u = ubeg;  u <= uend;  ++u)
			{
				int x, y, z;
				compXYZFromUVW(arg_fastAxis,arg_medAxis,arg_slowAxis,u-ubeg,v-vbeg,w-wbeg,x,y,z);
				float vpZ = inpbuf_vp[u-ubeg];
				if (vpZ > vpmax)
				{
#pragma omp critical
					{
						vpmax = vpZ;
					}
				}
			}
		}
	}
	delete [] rawinpbuf;

	if (logLevel >= 3)
	{
		printf("\rvpmax = %.2f                                                                                   \n",vpmax);
		fflush(stdout);
	}

	// clean up
	if (file_vp  != 0L) fclose(file_vp );

	return vpmax;
}

float Determine_Vpmax(
	int logLevel,
	char* filename_vp,
	float scalarval_vp,
	char* filename_eta,
	float scalarval_eta,
	int etaflag,
	char* filename_dta,
	float scalarval_dta,
	int swapin,
	int ubeg,
	int uend,
	int vbeg,
	int vend,
	int wbeg,
	int wend,
	int dimu,
	int dimv,
	int dimw,
	int arg_fastAxis,
	int arg_medAxis,
	int arg_slowAxis
	)
{
	if (ubeg < 0 || uend >= dimu || vbeg < 0 || vend >= dimv || wbeg < 0 || wend >= dimw)
	{
		printf("ERROR! Read_Single_Earth_Model_Attribute - Sub volume extends beyond dimensions of original volume!\n");
		printf("u = [%d, %d], dimu = %d\n",ubeg,uend,dimu);
		printf("v = [%d, %d], dimv = %d\n",vbeg,vend,dimv);
		printf("w = [%d, %d], dimw = %d\n",wbeg,wend,dimw);
	}	

	// Read from file or use constant value.
	bool constflag_vp = false; 
	float constval_vp = 0.0f;
	FILE* file_vp = 0L;
	if (filename_vp == 0L)
	{
		constval_vp = scalarval_vp;
		constflag_vp = true;
		if (logLevel >= 3) printf("nil filename: Using const val = %f\n", constval_vp);
	}
	else
	{
		if( !(file_vp = fopen(filename_vp,"r")))
		{ 
			constval_vp = atof(filename_vp);  
			constflag_vp=true;
			printf("%s not found: Using const val = %f\n", filename_vp,constval_vp); 
		}
	}

	bool constflag_eta = false; 
	float constval_eta = 0.0f;
	FILE* file_eta = 0L;
	if (filename_eta == 0L)
	{
		constval_eta = scalarval_eta;
		constflag_eta = true;
		if (logLevel >= 3) printf("nil filename: Using const val = %f\n", constval_eta);
	}
	else
	{
		if( !(file_eta = fopen(filename_eta,"r")))
		{ 
			constval_eta = atof(filename_eta);  
			constflag_eta=true;
			printf("%s not found: Using const val = %f\n", filename_eta,constval_eta); 
		}
	}

	bool constflag_dta = false; 
	float constval_dta = 0.0f;
	FILE* file_dta = 0L;
	if (etaflag)
	{
		if (filename_dta == 0L)
		{
			constval_dta = scalarval_dta;
			constflag_dta = true;
			if (logLevel >= 3) printf("nil filename: Using const val = %f\n", constval_dta);
		}
		else
		{
			if( !(file_dta = fopen(filename_dta,"r")))
			{ 
				constval_dta = atof(filename_dta);  
				constflag_dta=true;
				printf("%s not found: Using const val = %f\n", filename_dta,constval_dta); 
			}
		}
	}

	int actual_dimu = uend - ubeg + 1;

	// loop over u,v,w.
	// read, convert to right endian-ness and transpose to x-y-z ordering.
	int nthreads = 0;
#pragma omp parallel
	{
		nthreads = omp_get_num_threads();
	}
	if (nthreads < 1) nthreads = 1;  // just in case...
	float* rawinpbuf = new float[actual_dimu*nthreads*3];
	int nwpu = (wend-wbeg)/100;
	if (nwpu < 1) nwpu = 1;
	float vpmax = -1e37f;
#pragma omp parallel for schedule(dynamic)
	for (int w = wbeg;  w <= wend;  ++w)
	{
		int curr_thread = omp_get_thread_num();
		float* inpbuf_vp = rawinpbuf + actual_dimu * curr_thread * 3;
		float* inpbuf_eta = inpbuf_vp + actual_dimu;
		float* inpbuf_dta = inpbuf_eta + actual_dimu;
		if (logLevel >= 3 && ((w-wbeg) % nwpu) == 0)
		{
#pragma omp critical
			{
				printf("\r%s (%.0f%%)","Scanning Vp, Eta and Dta",100.0f*(float)(w-wbeg)/(float)(wend-wbeg));
				fflush(stdout);
			}
		}
		for (int v = vbeg;  v <= vend;  ++v)
		{
			long offset = (long)sizeof(float) * ((long)w * (long)dimu * (long)dimv + (long)v * (long)dimu + (long)ubeg);
			//fprintf(stderr, "Input offset = %ld u,v,w = [%d,%d,%d]\n", offset, ubeg, v, w);

			// read from files. we do this sequentially, since hitting the filesystem with many parallel requests most likely will trash performance.
			if (constflag_vp)
			{
				for (int u = ubeg;  u <= uend;  ++u) inpbuf_vp[u-ubeg] = constval_vp;
			}
			if (constflag_eta)
			{
				for (int u = ubeg;  u <= uend;  ++u) inpbuf_eta[u-ubeg] = constval_eta;
			}
			if (etaflag && constflag_dta)
			{
				for (int u = ubeg;  u <= uend;  ++u) inpbuf_dta[u-ubeg] = constval_dta;
			}
#pragma omp critical
			{
				if (!constflag_vp)
				{
					fseek(file_vp, offset, SEEK_SET);
					fread(inpbuf_vp, sizeof(float), actual_dimu, file_vp);
				}
				if (!constflag_eta)
				{
					fseek(file_eta, offset, SEEK_SET);
					fread(inpbuf_eta, sizeof(float), actual_dimu, file_eta);
				}
				if (etaflag && !constflag_dta)
				{
					fseek(file_dta, offset, SEEK_SET);
					fread(inpbuf_dta, sizeof(float), actual_dimu, file_dta);
				}
			}
			if(swapin)
			{
				if (!constflag_vp ) swap4bytes((int*)inpbuf_vp , actual_dimu);
				if (!constflag_eta) swap4bytes((int*)inpbuf_eta, actual_dimu);
				if (etaflag && !constflag_dta) swap4bytes((int*)inpbuf_dta, actual_dimu);
			}
			for (int u = ubeg;  u <= uend;  ++u)
			{
				int x, y, z;
				compXYZFromUVW(arg_fastAxis,arg_medAxis,arg_slowAxis,u-ubeg,v-vbeg,w-wbeg,x,y,z);
				float vpZ = inpbuf_vp[u-ubeg];
				float eps = inpbuf_eta[u-ubeg];
				if (etaflag)
				{
					// input is really eta(in EpsIn), so convert to eps
					float dta = inpbuf_dta[u-ubeg];
					eps = 0.5f*((1.0f+2.0f*eps)*(1.0f+2.0f*dta) - 1.0f);
				}
				float vpX = sqrtf(1.0f+2.0f*eps)*vpZ;
				float vp = vpX > vpZ ? vpX : vpZ;
				if (vp > vpmax)
				{
#pragma omp critical
					{
						vpmax = vp;
					}
				}
			}
		}
	}
	delete [] rawinpbuf;

	if (logLevel >= 3)
	{
		printf("\rvpmax = %.2f                                                                                   \n",vpmax);
		fflush(stdout);
	}

	// clean up
	if (file_vp  != 0L) fclose(file_vp );
	if (file_eta != 0L) fclose(file_eta);
	if (file_dta != 0L) fclose(file_dta);

	return vpmax;
}

void Read_Single_Earth_Model_Attribute(
	int logLevel,
	char* filename,
	float scalarval,
	int swapin,
	int ubeg,
	int uend,
	int vbeg,
	int vend,
	int wbeg,
	int wend,
	int dimu,
	int dimv,
	int dimw,
	int arg_fastAxis,
	int arg_medAxis,
	int arg_slowAxis,
	int zpad,
	float*** attr,
	float& attrmin,
	float& attrmax
	)
{
	if (ubeg < 0 || uend >= dimu || vbeg < 0 || vend >= dimv || wbeg < 0 || wend >= dimw)
	{
		printf("ERROR! Read_Single_Earth_Model_Attribute - Sub volume extends beyond dimensions of original volume!\n");
		printf("u = [%d, %d], dimu = %d\n",ubeg,uend,dimu);
		printf("v = [%d, %d], dimv = %d\n",vbeg,vend,dimv);
		printf("w = [%d, %d], dimw = %d\n",wbeg,wend,dimw);
	}	

	// Read from file or use constant value.
	bool constflag = false; 
	float constval = 0.0f;
	FILE* file = 0L;
	if (filename == 0L)
	{
		constval = scalarval;
		constflag = true;
		if (logLevel >= 3) printf("nil filename: Using const val = %f\n", constval);
	}
	else
	{
		if( !(file = fopen(filename,"r")))
		{ 
			constval = atof(filename);  
			constflag=true;
			printf("%s not found: Using const val = %f\n", filename,constval); 
		}
	}

	int actual_dimu = uend - ubeg + 1;

	// loop over u,v,w.
	// read, convert to right endian-ness and transpose to x-y-z ordering.
	if (constflag)
	{
#pragma omp parallel for schedule(dynamic)
		for (int w = wbeg;  w <= wend;  ++w)
		{
			for (int v = vbeg;  v <= vend;  ++v)
			{
				for (int u = ubeg;  u <= uend;  ++u)
				{
					int x, y, z;
					compXYZFromUVW(arg_fastAxis,arg_medAxis,arg_slowAxis,u-ubeg,v-vbeg,w-wbeg,x,y,z);
					attr[z+zpad][y][x] = constval;
				}
			}
		}
		attrmin = constval;
		attrmax = constval;
	}
	else
	{
		attrmin = 1e37f;
		attrmax = -1e37f;
		int nthreads = 0;
#pragma omp parallel
		{
			nthreads = omp_get_num_threads();
		}
		if (nthreads < 1) nthreads = 1;  // just in case...
		float* rawinpbuf = new float[actual_dimu*nthreads];
		int nwpu = (wend-wbeg)/100;
		if (nwpu < 1) nwpu = 1;
#pragma omp parallel for schedule(dynamic)
		for (int w = wbeg;  w <= wend;  ++w)
		{
			int curr_thread = omp_get_thread_num();
			float* inpbuf = rawinpbuf + actual_dimu * curr_thread;
			if (logLevel >= 3 && ((w-wbeg) % nwpu) == 0)
			{
#pragma omp critical
				{
					printf("\r%s (%.0f%%)",filename,100.0f*(float)(w-wbeg)/(float)(wend-wbeg));
					fflush(stdout);
				}
			}
			for (int v = vbeg;  v <= vend;  ++v)
			{
				long offset = (long)sizeof(float) * ((long)w * (long)dimu * (long)dimv + (long)v * (long)dimu + (long)ubeg);
				//fprintf(stderr, "Input offset = %ld u,v,w = [%d,%d,%d]\n", offset, ubeg, v, w);

				// read from files. we do this sequentially, since hitting the filesystem with many parallel requests most likely will trash performance.
#pragma omp critical
				{
					fseek(file, offset, SEEK_SET);
					fread(inpbuf, sizeof(float), actual_dimu, file);
				}
				if(swapin) swap4bytes((int*)inpbuf, actual_dimu);
				for (int u = ubeg;  u <= uend;  ++u)
				{
					int x, y, z;
					compXYZFromUVW(arg_fastAxis,arg_medAxis,arg_slowAxis,u-ubeg,v-vbeg,w-wbeg,x,y,z);
					float val = inpbuf[u-ubeg];
					attr[z+zpad][y][x] = val;
					if (val < attrmin)
					{
#pragma omp critical
						{
							attrmin = val;
							//printf("new value for attrmin = %f at offset %ld\n",attrmin,offset);
						}
					}
					if (val > attrmax)
					{
#pragma omp critical
						{
							attrmax = val;
							//printf("new value for attrmax = %f at offset %ld\n",attrmax,offset);
						}
					}
				}
			}
		}
		delete [] rawinpbuf;
	}

	if (logLevel >= 3 && !constflag)
	{
		printf("\r%s (min = %f, max = %f)\n",filename,attrmin,attrmax);
		fflush(stdout);
	}

	// clean up
	if (file != 0L) fclose(file);
}

void Read_Single_Earth_Model_Attribute(
	int logLevel,
	char* filename,
	float scalarval,
	int swapin,
	int ubeg,
	int uend,
	int vbeg,
	int vend,
	int wbeg,
	int wend,
	int dimu,
	int dimv,
	int dimw,
	int arg_fastAxis,
	int arg_medAxis,
	int arg_slowAxis,
	int zpad,
	float*** attr
	)
{
	float attrmin, attrmax;
	Read_Single_Earth_Model_Attribute(logLevel,filename,scalarval,swapin,ubeg,uend,vbeg,vend,wbeg,wend,dimu,dimv,dimw,arg_fastAxis,arg_medAxis,arg_slowAxis,zpad,attr,attrmin,attrmax);
}

// normalize angle to range [0,2*pi>
float Normalize_Angle_Radians(float angle, float center)
{
	const float pi = 3.141592653589793116f;
	const float two_pi = 6.283185307179586232f;
	return angle - two_pi * floor((angle + pi - center) / two_pi);
}
} // end of anonymous namespace

//
// Returns recout data.
// Contains receiver traces stored in this order X-Y-Z-T.
// T is time.
//
#ifndef STANDALONE
float* TTIDenQ_GetRecout()
{
	//return recvec;
	// TO-DO: Add back this feature for SeisSpace module. Need temporary file somewhere.
	return recvec;
}
#endif

char* getSourceVertInterpStr(int source_vert_interp)
{
	switch (source_vert_interp)
	{
	case 0:
		return "None";
	case 1:
		return "Reciprocal difference";
	case 2:
		return "Reciprocal average";
	default:
		return "???";
	}
}

float Get_Density(
		int KernelType,
		int*** VelAnis,
		int*** DenAng,
		float Denbinsize,
		float Denmin,
		int xsrc,
		int ysrc,
		int zsrc
		)
{
	if (KernelType == 0)
	{
		int iVelAnis = VelAnis[zsrc][ysrc][xsrc];
		return (float)((iVelAnis >> SHIFTDen) & DENMASK) * Denbinsize + Denmin;
	}
	else
	{
		int iDenAng = DenAng[zsrc][ysrc][xsrc];
		return (float)((iDenAng >> SHIFTDen) & DENMASK) * Denbinsize + Denmin;
	}
}

float Get_Vp(
		int KernelType,
		int*** VelAnis,
		int*** DenAng,
		float Velbinsize,
		float Velmin,
		int xsrc,
		int ysrc,
		int zsrc
	    )
{
	int iVelAnis = VelAnis[zsrc][ysrc][xsrc];
	return (float)(iVelAnis & VELMASK) * Velbinsize + Velmin;
}

float Get_Bulk_Modulus(
		int KernelType,
		int*** VelAnis,
		int*** DenAng,
		float Denbinsize,
                float Denmin,
		float Velbinsize,
                float Velmin,
		int xsrc,
		int ysrc,
		int zsrc
		)
{
	float Vp = Get_Vp(KernelType,VelAnis,DenAng,Velbinsize,Velmin,xsrc,ysrc,zsrc);
	float Dn = Get_Density(KernelType,VelAnis,DenAng,Denbinsize,Denmin,xsrc,ysrc,zsrc);
	return Dn * Vp * Vp;
}

float Compute_Source_Vertical_Interpolation_Term(
		int i,
		int source_vert_interp,
		float bsrcavg,
		float bulk_modulus,
		float KOBNrec
		)
{
	if (source_vert_interp == 1)
	{
		float stencil_coeff = 0.0f;
		switch (i)
		{
			case -4:
				stencil_coeff = -A5;
				break;
			case -3:
				stencil_coeff = -A4;
				break;
			case -2:
				stencil_coeff = -A3;
				break;
			case -1:
				stencil_coeff = -A2;
				break;
			case  0:
				stencil_coeff = -A1;
				break;
			case  1:
				stencil_coeff =  A1;
				break;
			case  2:
				stencil_coeff =  A2;
				break;
			case  3:
				stencil_coeff =  A3;
				break;
			case  4:
				stencil_coeff =  A4;
				break;
			case  5:
				stencil_coeff =  A5;
				break;
			default:
				stencil_coeff = 0.0f;
				break;
		}
		return bsrcavg * bulk_modulus/KOBNrec * stencil_coeff;
	}
	else if (source_vert_interp == 2)
	{
		float lagr1 = 0.605621338f;
		float lagr2 = -0.13458252f;
		float lagr3 = 0.034606934f;
		float lagr4 = -0.00617981f;
		float lagr5 = 0.000534058f;

		float stencil_coeff = 0.0f;
		switch (i)
		{
			case -4:
				stencil_coeff =  lagr5;
				break;
			case -3:
				stencil_coeff =  lagr4;
				break;
			case -2:
				stencil_coeff =  lagr3;
				break;
			case -1:
				stencil_coeff =  lagr2;
				break;
			case  0:
				stencil_coeff =  lagr1;
				break;
			case  1:
				stencil_coeff =  lagr1;
				break;
			case  2:
				stencil_coeff =  lagr2;
				break;
			case  3:
				stencil_coeff =  lagr3;
				break;
			case  4:
				stencil_coeff =  lagr4;
				break;
			case  5:
				stencil_coeff =  lagr5;
				break;
			default:
				stencil_coeff = 0.0f;
					break;
		}
		return bulk_modulus/KOBNrec * stencil_coeff;
	}
	return 0.0f;
}

void VarDenQ_ComputeShot(
	int num_obn_nodes,
	float* obn_node_locations,
	int max_num_threads,
	int dump_xz_itvl,
	int is_dry_run,
	unsigned long& max_allocated_memory,
	double& total_runtime,
	int KernelType,		// 0->ISO, 1->VTI, 2->TTI
	char* seisname,
	int id,
	char* vpname,
	char* epsetaname,
	int etaflag,
	char* deltaname,
	float VsoVp0,
	char* dnname,
	char* dipdxname,
	char* azmdyname,
	float rotationAngleWorldToLocalDegree,
	int dipxdipyflag,
	int degreesflag,
	int swapin,
	char* Qname,
	float scalarVp,
	float scalarEta,
	float scalarDelta,
	float scalarDen,
	float scalarTilt,
	float scalarAzm,
	float scalarQ,
	int smoothanisflag,
	int isorad,
	int arg_absorbz0,
	int arg_srcghost,
	int arg_recghost,
	int spongewidth_x,
	int spongewidth_y,
	int spongewidth_z_lo,
	int spongewidth_z_hi,
	float spongecoeff_x,
	float spongecoeff_y,
	float spongecoeff_z_lo,
	float spongecoeff_z_hi,
	float arg_dh,
	float arg_dz,
	float arg_stretchfacz,
	int arg_nx,
	int arg_ny,
	int arg_nz,
	int arg_fastAxis,	// applies to input files only. For all axis arguments, the following values are valid:
	int arg_medAxis,	// 0->x, 1->y, 2->z
	int arg_slowAxis,	// ex.: For the "regular" X-Y-Z ordering (X being the fast axis, Z being the slow, these 3 args will be 0,1,2 respectively.
	int arg_xbeg,
	int arg_xend,
	int arg_ybeg,
	int arg_yend,
	int arg_zbeg,
	int arg_zend,
	float srcx,
	float srcy,
	float srcz,
	int sourcetype,
	int source_vert_interp,
	char* stfname,
	int arg_xrecstart,
	int arg_xrecend,
	int arg_xrecstride,
	int arg_yrecstart,
	int arg_yrecend,
	int arg_yrecstride,
	int arg_zrecstart,
	int arg_zrecend,
	int arg_zrecstride,
	float maxtime,
	float timestartrec,
	float dtout,
	float newfmax,
	float gamfac,
	int OTflag,
	int logLevel
	)
{
	char hdrstring[2000],infostring[200];

	int i,j,k, x,y,z, nynzmax, zstart, nx2,ny2,nz2;
	int t, tsrc,tmax,tstartrec,tstartsim, srcboxH, srcboxZ, tvolsrc;
	int nrec, nyrec, nzrec, ntout, itout;
	int isohalfrad, isozstart, isoystart, isoxstart, isotapersrc=1, smthanisflag;
	int dnconstflag=0, vpconstflag=0, epsetaconstflag=0, deltaconstflag=0, dipdxconstflag=0, azmdyconstflag=0, Qconstflag=0;
	int velanis, ***VelAnisOrigsrc, velanisIsosrc;
	// int tid, nthread;  placeholder for threading, not currently used

	float dt, dtX, dtZ, dt2, gam, khmax, minfac;
	float eps, delta, C44, newC44ref, newVsoVp0, newVsoVp0max=0.f, C13thresh;
	float dnconst, vpconst, epsetaconst, deltaconst, dipdxconst, azmdyconst, Qconst;
	float tout, dtsrcfine, tfrac;
	float fmax, fmaxX, fmaxZ;
	float velsrc, velsrcX, velsrcZ, initsrc=1, stf[4000], smthnorm, cum, isofac;
	float val, vp,vpmin,vpmax,vpminX,vpmaxX,vpminZ,vpmaxZ;
	float Epsmax, Delmax, Denmax, Dipmax, Azmmax, Qmax;
	float Epsmin, Delmin, Denmin, Dipmin, Azmmin, Qmin;
	float invVelrange, invEpsrange, invDelrange, invDenrange, invDiprange, invAzmrange, invQrange;
	float Velbinsize,  Epsbinsize,  Delbinsize,  Denbinsize,  Dipbinsize,  Azmbinsize, Qbinsize;
	float tandx, tandy;
	float ***VpIn, ***EpsIn, ***DelIn;
	float ***swap;

	FILE *seisfile;
	FILE *vpfile, *epsetafile, *deltafile, *dnfile, *dipdxfile, *azmdyfile, *Qfile;

	initDone = 0;

	if (num_obn_nodes > 0)
	{
		// only need to load density for this run mode, so set everything else to (isotropic) constants.
		vpname = "4500.0";
		epsetaname = "0.0";
		etaflag = 0;
		deltaname = "0.0";
		VsoVp0 = 0.0f;
		dipdxname = "0.0";
		azmdyname = "0.0";
		rotationAngleWorldToLocalDegree = 0.0f;
		dipxdipyflag = 0;
		degreesflag = 0;
		Qname = "1e9";		
	}

	if (logLevel >= 4)
	{
		printf("num_obn_nodes = %d\n", num_obn_nodes);
		printf("max_num_threads = %d\n", max_num_threads);
		printf("dump_xz_itvl = %d\n", dump_xz_itvl);
		printf("seisname = %s\n", seisname != 0L ? seisname : "nil");
		printf("id = %d\n", id);
		printf("vpname = %s\n", vpname != 0L ? vpname : "nil");
		printf("epsetaname = %s\n", epsetaname != 0L ? epsetaname : "nil");
		printf("etaflag = %d\n", etaflag);
		printf("deltaname = %s\n", deltaname != 0L ? deltaname : "nil");
		printf("VsoVp0 = %f\n", VsoVp0);
		printf("dnname = %s\n", dnname != 0L ? dnname : "nil");
		printf("dipdxname = %s\n", dipdxname != 0L ? dipdxname : "nil");
		printf("azmdyname = %s\n", azmdyname != 0L ? azmdyname : "nil");
		printf("rotationAngleWorldToLocalDegree = %.2f\n",rotationAngleWorldToLocalDegree);
		printf("dipxdipyflag = %d\n", dipxdipyflag);
		printf("degreesflag = %d\n", degreesflag);
		printf("swapin = %d\n", swapin);
		printf("Qname = %s\n", Qname != 0L ? Qname : "nil");
		printf("scalarVp = %f\n", scalarVp);
		printf("scalarEta = %f\n", scalarEta);
		printf("scalarDelta = %f\n", scalarDelta);
		printf("scalarDen = %f\n", scalarDen);
		printf("scalarTilt = %f\n", scalarTilt);
		printf("scalarAzm = %f\n", scalarAzm);
		printf("scalarQ = %f\n", scalarQ);
		printf("smoothanisflag = %d\n", smoothanisflag);
		printf("isorad = %d\n", isorad);
		printf("arg_absorbz0 = %s\n", arg_absorbz0 ? "Absorbing" : "Free Surface");
		printf("arg_srcghost = %s\n", arg_srcghost ? "Yes" : "No");
		printf("arg_recghost = %s\n", arg_recghost ? "Yes" : "No");
		printf("spongewidth_x = %d\n", spongewidth_x);
		printf("spongewidth_y = %d\n", spongewidth_y);
		printf("spongewidth_z_lo = %d\n", spongewidth_z_lo);
		printf("spongewidth_z_hi = %d\n", spongewidth_z_hi);
		printf("spongecoeff_x = %e\n", spongecoeff_x);
		printf("spongecoeff_y = %e\n", spongecoeff_y);
		printf("spongecoeff_z_lo = %e\n", spongecoeff_z_lo);
		printf("spongecoeff_z_hi = %e\n", spongecoeff_z_hi);
		printf("arg_dh = %f\n", arg_dh);
		printf("arg_dz = %f\n", arg_dz);
		printf("arg_stretchfacz = %f\n", arg_stretchfacz);
		printf("arg_nx = %d\n", arg_nx);
		printf("arg_ny = %d\n", arg_ny);
		printf("arg_nz = %d\n", arg_nz);
		printf("arg_fastAxis = %s\n", getAxisSymbol(arg_fastAxis)); 
		printf("arg_medAxis = %s\n", getAxisSymbol(arg_medAxis)); 
		printf("arg_slowAxis = %s\n", getAxisSymbol(arg_slowAxis)); 
		printf("arg_xbeg = %d\n", arg_xbeg);
		printf("arg_xend = %d\n", arg_xend);
		printf("arg_ybeg = %d\n", arg_ybeg);
		printf("arg_yend = %d\n", arg_yend);
		printf("arg_zbeg = %d\n", arg_zbeg);
		printf("arg_zend = %d\n", arg_zend);
		printf("srcx = %f\n", srcx);
		printf("srcy = %f\n", srcy);
		printf("srcz = %f\n", srcz);
		printf("sourcetype = %d\n",  sourcetype);
		printf("source_vert_interp = %d (%s)\n", source_vert_interp, getSourceVertInterpStr(source_vert_interp));
		printf("stfname = %s\n", stfname);
		printf("arg_xrecstart = %d\n", arg_xrecstart);
		printf("arg_xrecend = %d\n", arg_xrecend);
		printf("arg_xrecstride = %d\n", arg_xrecstride);
		printf("arg_yrecstart = %d\n", arg_yrecstart);
		printf("arg_yrecend = %d\n", arg_yrecend);
		printf("arg_yrecstride = %d\n", arg_yrecstride);
		printf("arg_zrecstart = %d\n", arg_zrecstart);
		printf("arg_zrecend = %d\n", arg_zrecend);
		printf("arg_zrecstride = %d\n", arg_zrecstride);
		printf("maxtime = %f\n", maxtime);
		printf("timestartrec = %f\n", timestartrec);
		printf("dtout = %f\n", dtout);
		printf("newfmax = %f\n", newfmax);
		printf("gamfac = %f\n", gamfac);
		printf("OTflag = %d\n", OTflag);
		printf("logLevel = %d\n", logLevel);
	}

	// KernelType == 0 -> Isotropic propagator. No need for isorad.
	if (KernelType == 0) isorad = 0;

#ifdef STANDALONE
	char segyfilename[1024];
	genSEGYFilename(seisname, id, segyfilename);
	char tmpfilename[1024];
	sprintf(tmpfilename,"%s.tmp",segyfilename);

	// verify that SEGY file can be written before we do X hours of propagation work!
	FILE* fp_tmp = 0L;
	if (is_dry_run == 0)
	{
		fp_tmp = fopen(tmpfilename, "w");
		if (fp_tmp == 0L)
		{
			fprintf(stderr, "ERROR! Cannot open %s for writing!",tmpfilename);
			exit(-1);
		}
	}
#endif

	// Enable FTZ and DAZ modes.
	unsigned int old_csr = _mm_getcsr();
	_mm_setcsr(old_csr|0x8940);

	absorbz0 = arg_absorbz0;
	srcghost = arg_srcghost;
	recghost = arg_recghost;

	dh = arg_dh;
	dz = arg_dz;
	stretchfacz = arg_stretchfacz;

	/*** S E T  C O N S T A N T  C O E F F I C I E N T S ***/
	if(OTflag==4)  
	{ 
		khmax = 2.41f;  gam = 0.6601f*gamfac; 
	}
	else /* OT2 */ 
	{
		khmax = 1.76f;  gam = 0.3801f*gamfac;
	}

	int actual_nx = arg_xend - arg_xbeg + 1;
	nx = ((actual_nx + 7) >> 3) << 3; 		// arg_nx;
	ny = arg_yend - arg_ybeg + 1; 			// arg_ny;
	nz = arg_zend - arg_zbeg + 1; 			// arg_nz;
	int num_x_zero_pad = nx - actual_nx;
	if (logLevel >= 4)
	{
		printf("nx = %d, actual_nx = %d, num_x_zero_pad = %d\n",nx,actual_nx,num_x_zero_pad);
	}

	xrecstart = arg_xrecstart - arg_xbeg;
	xrecend = arg_xrecend - arg_xbeg;
	xrecstride = arg_xrecstride;

	yrecstart = arg_yrecstart - arg_ybeg;
	yrecend = arg_yrecend - arg_ybeg;
	yrecstride = arg_yrecstride;

	zrecstart = arg_zrecstart - arg_zbeg;
	zrecend = arg_zrecend - arg_zbeg;
	zrecstride = arg_zrecstride;

	int cache_size_per_core_KB;
	int numCPU = Get_Physical_Core_Count(cache_size_per_core_KB);
	// initial tile size. used until full volume is reached.
	// ISO cache < VTI cache < TTI cache.
	int bsX, bsY;
	switch (KernelType)
	{
		case 0:  // ISO
			bsX = 96;
			bsY = 32;
			break;
		case 1:  // VTI
			bsX = 64;
			bsY = 32;
			break;
		case 2:  // TTI
			bsX = 48;
			bsY = 32;
			break;
		default:
			bsX = 32;
			bsY = 32;
			break;
	}

	int expanding_box = 1;
	if (max_num_threads < 0)
	{
		expanding_box = 0;
		max_num_threads = -max_num_threads;
	}
	if (!expanding_box)
	{
		printf("EXPANDING BOX DISABLED! This will increase the run time.\n");
	}
	if (max_num_threads > numCPU)
	{
		max_num_threads = numCPU;
	}
	omp_set_num_threads(max_num_threads);

	if (is_dry_run) dump_xz_itvl = 0;  // disable this for dry runs.
	int Supports_AVX = SupportsAVX();
	if (KernelType == 0) Supports_AVX = 0;  // No ISO AVX kernel exists.
	if (dump_xz_itvl < 0)
	{
		Supports_AVX = 0;
	}
	if (dump_xz_itvl > 1)
	{
		printf("OUTPUT OF XZ CROSS SECTION EVERY %d TIMESTEPS ENABLED!\n",dump_xz_itvl);
	}

	if (is_dry_run == 3)
	{
		// disable these fields in the earth model, since they have no impact on run times.
		dnname = 0L;
		dipdxname = 0L;
		azmdyname = 0L;
		Qname = 0L;
		scalarDen = 1.0f;
		scalarTilt = 0.0f;
		scalarAzm = 0.0f;
		scalarQ = 1e9f;
	}

	int num_threads = 0;
#pragma omp parallel
	{
		num_threads = omp_get_num_threads();
	}

	Generate_Tile_Shapes(logLevel,is_dry_run,KernelType, min(32,nx), min(128,nx), 8, min(16,ny), min(128,ny), 8, num_threads, (unsigned long)(cache_size_per_core_KB*num_threads));

	if (logLevel >= 4)
	{
		printf("OPENMP will use %d threads. Node has %d physical CPU cores.\n",num_threads,numCPU);
		if (cache_size_per_core_KB > 0) printf("Total available cache is %d KB.\n",cache_size_per_core_KB*numCPU);
	}

	max_allocated_memory = 0L;

	struct timespec glob_start;
	clock_gettime(CLOCK_REALTIME, &glob_start);

	/*
	if(xrecend < xrecstart) { xrecend=xrecstart; xrecstride=1; }
	if(yrecend < yrecstart) { yrecend=yrecstart; yrecstride=1; }
	if(zrecend < zrecstart) { zrecend=zrecstart; zrecstride=1; }
	*/

	/** S O U R C E  G E O M E T R Y **/
	xsrc = (int)(srcx/dh+0.5f) - arg_xbeg;
	ysrc = (int)(srcy/dh+0.5f) - arg_ybeg;
	zsrc = (int)(srcz/dz+0.5f) - arg_zbeg;
	if(zsrc<1) zsrc=1;

	// Adjust skylayer to include deepest source and receiver.
	int skylayer = 0;
	if (absorbz0)
	{
		int deep_z = 0;
		if (arg_srcghost && zsrc > deep_z) deep_z = zsrc;
		if (arg_recghost && zrecend > deep_z) deep_z = zrecend;
		skylayer = spongewidth_z_lo + 9 + deep_z;
		if (logLevel >= 2)
		{
			printf("Absorbing z surface. SKYLAYER = %d\n",skylayer);
		}
	}
	else
	{
		if (logLevel >= 2)
		{
			printf("Reflecting free surface.\n");
		}
	}

	if(absorbz0)  nz += skylayer;
	nx1 = nx-1;  ny1 = ny-1;  nz1 = nz-1;
	nx2 = nx-2;  ny2 = ny-2;  nz2 = nz-2;
	if(OTflag != 2  &&  OTflag != 4)
	{ 
		fprintf(stderr,"OTflag must = 2 or 4\n"); exit(0); 
	}

	// don't need DenAng for ISO kernel
	if (KernelType > 0)
	{
		max_allocated_memory += Allocate_Padded_Earth_Model(logLevel,is_dry_run,PadDenAng,nx,ny,nz,xh,yh,zh,DenAng);
	}
	else
	{
		DenAng = 0L;
		PadDenAng = 0L;
	}
	max_allocated_memory += Allocate_Padded_Earth_Model(logLevel,is_dry_run,PadVelAnis,nx,ny,nz,xh,yh,zh,VelAnis);

	max_allocated_memory += Allocate_Propagation_Fields(logLevel,is_dry_run,OTflag,KernelType,pq,rs,Apq,nx,ny,nz,xh,yh,zh,V1,V2,V3);

	switch (KernelType)
	{
		case 0:
			if (logLevel >= 2) printf("\nI S O   T %d\n\n",OTflag);
			break;
		case 1:
			if (logLevel >= 2) printf("\nV T I   T %d\n\n",OTflag);
			break;
		case 2:
			if (logLevel >= 2) printf("\nT T I   T %d\n\n",OTflag);
			break;
		default:
			printf("\nUNKNOWN PROPAGATOR TYPE\n\n");
			exit(-1);
	}

	if (is_dry_run)
	{
		if (KernelType == 0)
		{
			printf("\nIsotropic kernel not supported yet!\nNo information is available.\n\n");
			return;
		}

		if (logLevel >= 2) printf("\nD R Y   R U N   L E V E L   %d\n\n",is_dry_run);

		if (is_dry_run == 1)
		{
			float appr_vpmax = 15000.0f;
			float appr_dt = gam * ( (dz/appr_vpmax < dh/appr_vpmax) ? dz/appr_vpmax : dh/appr_vpmax );	
			int appr_tmax = (int)(maxtime/appr_dt) + 1;

			double per_core_throughput = 0.0;
			if (KernelType == 1)
			{
				// VTI
				per_core_throughput = 30.0;
			}
			else if (KernelType == 2)
			{
				// TTI
				per_core_throughput = 20.0;
			}
			if (OTflag == 4) per_core_throughput /= 2.0;
			double total_throughput = per_core_throughput * (double)num_threads;

			double mcells = (double)nx * (double)ny * (double)nz * 1e-6;
			double seconds_per_timestep = mcells / total_throughput;
			total_runtime = (double)appr_tmax * seconds_per_timestep * 0.9;

			if (logLevel >= 1) 
			{
				printf("In order to run this job, a total of %.2f GB of free memory is required.\n",(double)max_allocated_memory / 1073741824.0);
				printf("Approximately %d timestep are required.\n",appr_tmax);
				printf("Approximate throughput is %.0f MCells/s.\n",total_throughput);
				printf("Approximate run time is %.2f hours\n\n",total_runtime/3600.0);
			}

			return;
		}
	}

	/*** OPTIMIZED CONVOLUTION COEFFS from global optimization over k = w/v ***/
	/* FIRST DERIVATIVE COEFFS: */
	/*
	   if(OTflag==4)
	   { 
	   a1 =  1.250394163714f;  a2 = -0.119656543874f;  a3 = 0.031206223579f;
	   a4 = -0.009128136972f;  a5 =  0.001882183398f;
	   }
	   else
	   { 
	   a1 =  1.248489029341f;  a2 = -0.120133754290f;  a3 = 0.031688119039f;
	   a4 = -0.008048796917f;  a5 =  0.001090357653f;
	   }
	   b1 =  1.231650129521f;  b2 = -0.103861125624f;  b3 = 0.020166542235f;
	   b4 = -0.002985637689f;
	   c1 =  1.199634495725f;  c2 = -0.080370339530f;  c3 = 0.008295304573f;
	   d1 =  1.134389630713f;  d2 = -0.044796543571f;
	   e1 =  0.5f;

	   A1h = a1/dh;  A2h = a2/dh;  A3h = a3/dh;  A4h = a4/dh;  A5h = a5/dh;
	   B1h = b1/dh;  B2h = b2/dh;  B3h = b3/dh;  B4h = b4/dh;
	   C1h = c1/dh;  C2h = c2/dh;  C3h = c3/dh;
	   D1h = d1/dh;  D2h = d2/dh;
	   E1h = e1/dh;
	 */
	halfOperLen=OperLen/2; halfOperLen1=halfOperLen-1;

	/*** M E M O R Y  A L L O C A T I O N  ***/
	if (is_dry_run == 0 || is_dry_run == 3)
	{
		if (logLevel >= 4) fprintf(stdout,"Allocating memory for lookup tables...\n");

		lutVp2 = (float*)malloc((VELMASK+1)*sizeof(float));
		lutDen = (float*)malloc((DENMASK+1)*sizeof(float));
		lutBuoy = (float*)malloc((DENMASK+1)*sizeof(float));
		lutQ = (float*)malloc((QMASK+1)*sizeof(float));
		if (
				lutVp2 == 0L ||
				lutDen == 0L ||
				lutBuoy == 0L ||
				lutQ == 0L 
		   )
		{
			printf("FAILED TO ALLOCATE LOOKUP TABLES!\n");
			exit(-1);
		}

		if (KernelType > 0)
		{
			lutEps = (float*)malloc((EPSMASK+1)*sizeof(float));
			lutDel = (float*)malloc((DELMASK+1)*sizeof(float));
			lutc44c33 = (float*)malloc((C44C33MASK+1)*sizeof(float));
			if ( 
					lutEps == 0L ||
					lutDel == 0L ||
					lutc44c33 == 0L

			   )
			{
				printf("FAILED TO ALLOCATE LOOKUP TABLES!\n");
				exit(-1);
			}
		}
		else
		{
			lutEps = 0L;
			lutDel = 0L;
			lutc44c33 = 0L;
		}

		if (KernelType > 1)
		{
			lutsDip = (float*)malloc((DIPMASK+1)*sizeof(float));
			lutcDip = (float*)malloc((DIPMASK+1)*sizeof(float));
			lutsAzm = (float*)malloc((AZMMASK+1)*sizeof(float));
			lutcAzm = (float*)malloc((AZMMASK+1)*sizeof(float));
			if (
					lutsDip == 0L ||
					lutcDip == 0L ||
					lutsAzm == 0L ||
					lutcAzm == 0L
			   )
			{
				printf("FAILED TO ALLOCATE LOOKUP TABLES!\n");
				exit(-1);
			}
		}
		else
		{
			lutsDip = 0L;
			lutcDip = 0L;
			lutsAzm = 0L;
			lutcAzm = 0L;
		}
	}

	/*** R E A D  Vpvertical, Eps, Delta, Den, Dip, Azm  F I L E S (X=fast axis) ***/
	if (logLevel >= 2) fprintf(stdout,"Reading Earth Model Files...\n");
	// alias pointer names for readability:
	VpIn = V1; EpsIn = V2; DelIn = V3;

	vpmaxX = vpmaxZ = Epsmax = Delmax = Denmax = Dipmax = Azmmax = Qmax = -10.0f;
	vpminX = vpminZ = Epsmin = Delmin = Denmin = Dipmin = Azmmin = Qmin =  1.e6f;

	// Read subvolume from input file(s).
	// Note that any axis ordering is allowed in input files, not just X-Y-Z.
	// Read in U-V-W order, then tranpose to X-Y-Z on-the-fly.
	int dimu, actual_dimu, dimv, dimw;
	int ubeg, uend, vbeg, vend, wbeg, wend;
	compUVWRangesFromXYZ(
			arg_fastAxis,arg_medAxis,arg_slowAxis,arg_xbeg,arg_xend,arg_ybeg,arg_yend,arg_zbeg,arg_zend,arg_nx,arg_ny,arg_nz,
			dimu,actual_dimu,dimv,dimw,ubeg,uend,vbeg,vend,wbeg,wend);
	if (logLevel >= 4)
	{
		printf("dimu = %d, actual_dimu = %d\n",dimu,actual_dimu);
		printf("dimv = %d\n",dimv);
		printf("dimw = %d\n",dimw);
		printf("ubeg = %d\n",ubeg);
		printf("uend = %d\n",uend);
		printf("vbeg = %d\n",vbeg);
		printf("vend = %d\n",vend);
		printf("wbeg = %d\n",wbeg);
		printf("wend = %d\n",wend);
	}

	if(absorbz0) zstart = skylayer;  else zstart = 0;

	// read Vp and Eps first.

	if (is_dry_run == 2)
	{
		float vpmax = 4500.0f;
		if (KernelType > 0)
		{
			vpmax = Determine_Vpmax(logLevel,vpname,scalarVp,epsetaname,scalarEta,etaflag,deltaname,scalarDelta,swapin,ubeg,uend,vbeg,vend,wbeg,wend,dimu,dimv,dimw,arg_fastAxis,arg_medAxis,arg_slowAxis);
		}
		else
		{
			vpmax = ISO_Determine_Vpmax(logLevel,vpname,scalarVp,swapin,ubeg,uend,vbeg,vend,wbeg,wend,dimu,dimv,dimw,arg_fastAxis,arg_medAxis,arg_slowAxis);
		}
		
		float dt = gam * ( (dz/vpmax < dh/vpmax) ? dz/vpmax : dh/vpmax );	
		int tmax = (int)(maxtime/dt) + 1;

		double per_core_throughput = 0.0;
		if (KernelType == 1)
		{
			// VTI
			per_core_throughput = 30.0;
		}
		else if (KernelType == 2)
		{
			// TTI
			per_core_throughput = 20.0;
		}
		if (OTflag == 4) per_core_throughput /= 2.0;
		double total_throughput = per_core_throughput * (double)num_threads;

		double mcells = (double)nx * (double)ny * (double)nz * 1e-6;
		double seconds_per_timestep = mcells / total_throughput;
		total_runtime = (double)tmax * seconds_per_timestep * 0.9;

		if (logLevel >= 1)
		{
			printf("\nIn order to run this job, a total of %.2f GB of free memory is required.\n",(double)max_allocated_memory / 1073741824.0);
			printf("%d timestep are required.\n",tmax);
			printf("Approximate throughput is %.0f MCells/s.\n",total_throughput);
			printf("Approximate run time is %.2f hours\n\n",total_runtime/3600.0);
		}

		return;
	}
	else
	{
		if (KernelType > 0)
		{
			Read_Single_Earth_Model_Attribute(logLevel,vpname,scalarVp,swapin,ubeg,uend,vbeg,vend,wbeg,wend,dimu,dimv,dimw,arg_fastAxis,arg_medAxis,arg_slowAxis,zstart,VpIn,vpminZ,vpmaxZ);
			Read_Single_Earth_Model_Attribute(logLevel,deltaname,scalarDelta,swapin,ubeg,uend,vbeg,vend,wbeg,wend,dimu,dimv,dimw,arg_fastAxis,arg_medAxis,arg_slowAxis,zstart,DelIn,Delmin,Delmax);
			Read_Single_Earth_Model_Attribute(logLevel,epsetaname,scalarEta,swapin,ubeg,uend,vbeg,vend,wbeg,wend,dimu,dimv,dimw,arg_fastAxis,arg_medAxis,arg_slowAxis,zstart,EpsIn);

			vpminX = Epsmin =  1e37;
			vpmaxX = Epsmax = -1e37;
#pragma omp parallel for schedule(dynamic)
			for(int z=zstart; z<nz; z++)
			{
				float _vpmaxX = vpmaxX;
				float _vpminX = vpminX;

				float _Epsmax = Epsmax;
				float _Epsmin = Epsmin;

				for(int y=0; y<ny; y++)
				{
					// calc min & max values
					for(int x=0; x<actual_nx; x++)
					{ 
						if(etaflag) // input is really eta(in EpsIn), so convert to eps
							EpsIn[z][y][x] = 0.5f*((1.0f+2.0f*EpsIn[z][y][x])*(1.0f+2.0f*DelIn[z][y][x]) - 1.0f);
						float val = EpsIn[z][y][x];
						if (val > _Epsmax) _Epsmax = val;
						if (val < _Epsmin) _Epsmin = val;

						val = sqrtf(1.0f+2.0f*EpsIn[z][y][x])*VpIn[z][y][x];
						if (val > _vpmaxX) _vpmaxX = val;  
						if (val < _vpminX) _vpminX = val;
					}
				}
#pragma omp critical
				{
					if (_Epsmin < Epsmin) Epsmin = _Epsmin;
					if (_Epsmax > Epsmax) Epsmax = _Epsmax;

					if (_vpminX < vpminX) vpminX = _vpminX;
					if (_vpmaxX > vpmaxX) vpmaxX = _vpmaxX;
				}
			}

			// compute vminX, vmaxX, vminZ and vmaxZ.
			vpmin = (vpminZ<vpminX)?vpminZ:vpminX;
			vpmax = (vpmaxX>vpmaxZ)?vpmaxX:vpmaxZ;
			// joe: need to adjust minfac according to actual rotated vels and expanding dz
			minfac = (vpminX/dh < vpminZ/dz)?vpminX/dh:vpminZ/dz;
			fmax = khmax*minfac/(2.0f*M_PI);
			// dt = gam * ( (dz/vpmaxZ < dh/vpmaxX) ? dz/vpmaxZ : dh/vpmaxX );
			dt = gam * ( (dz/vpmax < dh/vpmax) ? dz/vpmax : dh/vpmax );

			//vpmaxZ += 0.1f;
			//vpmaxX += 0.1f;
			//Epsmax += 0.00001f;
			invVelrange = (vpminZ == vpmaxZ) ? 0.0f : 1.0f / (vpmaxZ-vpminZ);
			Velbinsize = (vpminZ == vpmaxZ) ? 0.0f : (vpmaxZ-vpminZ)/VELMASK;

			invEpsrange = (Epsmin == Epsmax) ? 0.0f : 1.0f / (Epsmax-Epsmin);
			Epsbinsize = (Epsmin == Epsmax) ? 0.0f : (Epsmax-Epsmin)/EPSMASK;

			invDelrange = (Delmin == Delmax) ? 0.0f : 1.0f / (Delmax-Delmin);
			Delbinsize = (Delmin == Delmax) ? 0.0f : (Delmax-Delmin)/DELMASK;

			if (logLevel >= 4)
			{
				printf("vpminZ = %f, vpmaxZ = %f, invVelrange = %e, Velbinsize = %e\n",vpminZ,vpmaxZ,invVelrange,Velbinsize);
				printf("Epsmin = %f, Epsmax = %f, invEpsrange = %e, Epsbinsize = %e\n",Epsmin,Epsmax,invEpsrange,Epsbinsize);
				printf("Delmin = %f, Delmax = %f, invDelrange = %e, Delbinsize = %e\n",Delmin,Delmax,invDelrange,Delbinsize);
			}

			// add Vp, Eps to compressed earth model
#pragma omp parallel for schedule(dynamic)
			for(int z=zstart; z<nz; z++)
			{
				for(int y=0; y<ny; y++)
				{
					for(int x=0; x<actual_nx; x++)
					{ 
						// Pack VelAnis index array
						unsigned int indexVpdt2 = (unsigned int)(((VpIn[z][y][x]-vpminZ)*invVelrange*VELMASK) + 0.5f);
						unsigned int indexEps = (unsigned int)(((EpsIn[z][y][x]-Epsmin)*invEpsrange*EPSMASK) + 0.5f);
						unsigned int indexDel = (unsigned int)(((DelIn[z][y][x]-Delmin)*invDelrange*DELMASK) + 0.5f);

						unsigned int indexC44C33;
						if(indexEps == 0 && indexDel == 0) indexC44C33 = 0; // isotropic
						else                               indexC44C33 = (unsigned int)1;

						VelAnis[z][y][x] = (indexVpdt2 & VELMASK) | ((indexEps & EPSMASK) << SHIFTEps) | ((indexDel & DELMASK) << SHIFTDel) | ((indexC44C33 & C44C33MASK) << SHIFTC44C33);
					}
				}
			}

			// compute fmax. need dt for that.
			if     (newfmax >  1.1f) fmax  =  newfmax;
			else if(newfmax < -0.1f) fmax *= -newfmax;
			tmax = (int)(maxtime/dt) + 1;
			float fq = fmax / 3.0f;

			Read_Single_Earth_Model_Attribute(logLevel,dnname,scalarDen,swapin,ubeg,uend,vbeg,vend,wbeg,wend,dimu,dimv,dimw,arg_fastAxis,arg_medAxis,arg_slowAxis,zstart,VpIn,Denmin,Denmax);

			invDenrange = (Denmin == Denmax) ? 0.0f : 1.0f / (Denmax-Denmin);
			Denbinsize = (Denmin == Denmax) ? 0.0f : (Denmax-Denmin)/DENMASK;

			if (logLevel >= 4)
			{
				printf("Denmin = %f, Denmax = %f, invDenrange = %e, Denbinsize = %e\n",Denmin,Denmax,invDenrange,Denbinsize);
			}

#pragma omp parallel for schedule(dynamic)
			for(int z=zstart; z<nz; z++)
			{
				for(int y=0; y<ny; y++)
				{
					for(int x=0; x<actual_nx; x++)
					{ 
						// Pack DenAng index array
						unsigned int indexDen = (unsigned int)(((VpIn[z][y][x]-Denmin)*invDenrange*DENMASK) + 0.5f);

						DenAng[z][y][x] = ((indexDen & DENMASK) << SHIFTDen);
					}
				}
			}

			if (KernelType > 1)
			{
				float rotationAngleWorldToLocal = rotationAngleWorldToLocalDegree * d2r;
				Read_Single_Earth_Model_Attribute(logLevel,dipdxname,scalarTilt,swapin,ubeg,uend,vbeg,vend,wbeg,wend,dimu,dimv,dimw,arg_fastAxis,arg_medAxis,arg_slowAxis,zstart,VpIn);
				Read_Single_Earth_Model_Attribute(logLevel,azmdyname,scalarAzm,swapin,ubeg,uend,vbeg,vend,wbeg,wend,dimu,dimv,dimw,arg_fastAxis,arg_medAxis,arg_slowAxis,zstart,EpsIn);

				Dipmin = Azmmin =  1e37;
				Dipmax = Azmmax = -1e37;
#pragma omp parallel for schedule(dynamic)
				for(int z=zstart; z<nz; z++)
				{
					float _Dipmax = Dipmax;
					float _Dipmin = Dipmin;

					float _Azmmax = Azmmax;
					float _Azmmin = Azmmin;

					for(int y=0; y<ny; y++)
					{
						// calc min & max values
						for(int x=0; x<actual_nx; x++)
						{
							if(degreesflag)
							{
								VpIn[z][y][x] *= d2r; EpsIn[z][y][x] *= d2r;
							}
							if(dipxdipyflag) // input is dipx(in VpIn) & dipy(in EpsIn), so convert to dip,azm 
							{
								tandx = tanf(VpIn[z][y][x]);  tandy = tanf(EpsIn[z][y][x]);
								VpIn[z][y][x] = sqrtf(tandx*tandx + tandy*tandy);
								EpsIn[z][y][x] = atan2f(tandy,tandx);
							}

							// TMJ 07/25/12 - Rotate azimuth angles to local coordinate system.
							EpsIn[z][y][x] = Normalize_Angle_Radians(EpsIn[z][y][x] + rotationAngleWorldToLocal, 0.0f);
							VpIn[z][y][x] = Normalize_Angle_Radians(VpIn[z][y][x], 0.0f);

							val = VpIn[z][y][x];
							if (val > _Dipmax) _Dipmax=val;
							if (val < _Dipmin) _Dipmin=val;

							val = EpsIn[z][y][x];
							if (val > _Azmmax) _Azmmax=val;
							if (val < _Azmmin) _Azmmin=val;
						}
					}
#pragma omp critical
					{
						if (_Dipmin < Dipmin) Dipmin = _Dipmin;
						if (_Dipmax > Dipmax) Dipmax = _Dipmax;

						if (_Azmmin < Azmmin) Azmmin = _Azmmin;
						if (_Azmmax > Azmmax) Azmmax = _Azmmax;
					}
				}

				// demote kernel if possible.
				/*
				if (is_dry_run == 0 && KernelType == 2 && Dipmin == 0.0f && Dipmax == 0.0f)
				{
					if (Epsmin == 0.0f && Epsmax == 0.0f && Delmin == 0.0f && Delmax == 0.0f)
					{
						printf("\nTilt angles and anisotropic parameters are zero everywhere!\nDemoting kernel type from TTI to ISO.\n\n");
						Azmmin = 0.0f;
						Azmmax = 0.0f;
						KernelType = 0;
					}
					else
					{
						printf("\nTilt angles are zero everywhere!\nDemoting kernel type from TTI to VTI.\n\n");
						Azmmin = 0.0f;
						Azmmax = 0.0f;
						KernelType = 1;
					}
				}
				*/

				//Dipmax += 0.01f;
				invDiprange = (Dipmin == Dipmax) ? 0.0f : 1.0f / (Dipmax-Dipmin);
				Dipbinsize = (Dipmin == Dipmax) ? 0.0f : (Dipmax-Dipmin)/DIPMASK;

				//Azmmax += 0.01f;
				invAzmrange = (Azmmin == Azmmax) ? 0.0f : 1.0f / (Azmmax-Azmmin);
				Azmbinsize = (Azmmin == Azmmax) ? 0.0f : (Azmmax-Azmmin)/AZMMASK;

				if (logLevel >= 4)
				{
					printf("Dipmin = %f, Dipmax = %f, invDiprange = %e, Dipbinsize = %e\n",Dipmin,Dipmax,invDiprange,Dipbinsize);
					printf("Azmmin = %f, Azmmax = %f, invAzmrange = %e, Azmbinsize = %e\n",Azmmin,Azmmax,invAzmrange,Azmbinsize);
				}

				if (KernelType > 1)
				{
#pragma omp parallel for schedule(dynamic)
					for(int z=zstart; z<nz; z++)
					{
						for(int y=0; y<ny; y++)
						{
							for(int x=0; x<actual_nx; x++)
							{ 
								// Pack DenAng index array
								unsigned int indexDip = (unsigned int)(((VpIn[z][y][x]-Dipmin)*invDiprange*DIPMASK) + 0.5f);
								unsigned int indexAzm = (unsigned int)(((EpsIn[z][y][x]-Azmmin)*invAzmrange*AZMMASK) + 0.5f);

								DenAng[z][y][x] |= (indexDip & DIPMASK) | ((indexAzm & AZMMASK) << SHIFTAzm);
							}
						}
					}
				}
			}

			Read_Single_Earth_Model_Attribute(logLevel,Qname,scalarQ,swapin,ubeg,uend,vbeg,vend,wbeg,wend,dimu,dimv,dimw,arg_fastAxis,arg_medAxis,arg_slowAxis,zstart,VpIn);

			Qmin = 1e37f;
			Qmax = -1e37f;
#pragma omp parallel for schedule(dynamic)
			for(int z=zstart; z<nz; z++)
			{
				float _min, _max;
#pragma omp critical
				{
					_min = Qmin;
					_max = Qmax;
				}
				for(int y=0; y<ny; y++)
				{ 
					// calc min & max values
					for(int x=0; x<actual_nx; x++)
					{
						float val = VpIn[z][y][x];
						// transform Qatten into a scalar multiplier at each node location: this makes Q proportional to freq,
						// because only one dominant freq is used to represent all frequencies
						// Qatten = 1-(PI/Q)*dt/Tdom = 1-PI/Q*dt*fmax/3 = 1-0.333*PI*fmax*dt/Q (applied to amps, not energy)
						//val = (val > 1e8f) ? 1.0f : 1.0f - 0.333f*M_PI*fmax*dt/val;
						val = (val > 1e8f) ? 1.0f : exp(-M_PI*fq*dt/val);
						VpIn[z][y][x] = val;
						if (val > _max) _max=val;  
						if (val < _min) _min=val;
					}
				}
#pragma omp critical
				{
					if (_min < Qmin) Qmin = _min;
					if (_max > Qmax) Qmax = _max;
				}
			}

			if (Qmax > Qmin)
			{
				invQrange = (float)((double)QMASK / ((double)Qmax-(double)Qmin));
				Qbinsize = (float)(((double)Qmax-(double)Qmin)/(double)QMASK);
			}
			else
			{
				invQrange = Qbinsize = 0.0f;
			}

			if (logLevel >= 4)
			{
				printf("Qmin = %e, Qmax = %e, invQrange = %e, Qbinsize = %e\n",Qmin,Qmax,invQrange,Qbinsize);
			}

#pragma omp parallel for schedule(dynamic) 
			for(int z=zstart; z<nz; z++)
			{
				for(int y=0; y<ny; y++)
				{
					for(int x=0; x<actual_nx; x++)
					{ 
						// Pack DenAng index array
						unsigned int indexQ = (unsigned int)((VpIn[z][y][x]-Qmin)*invQrange + 0.5f);

						DenAng[z][y][x] |= ((indexQ & QMASK) << SHIFTQ);
					}
				}
			}

			if (KernelType == 0)
			{
				// kernel was demoted from TTI to ISO.
				// must repack earth model for ISO.
#pragma omp parallel for schedule(dynamic)
				for(int z=zstart; z<nz; z++)
				{
					for(int y=0; y<ny; y++)
					{
						for(int x=0; x<actual_nx; x++)
						{
							unsigned int indexVpdt2 = VelAnis[z][y][x] & VELMASK;
							unsigned int indexDen = (DenAng[z][y][x] >> SHIFTDen) & DENMASK;
							unsigned int indexQ = (DenAng[z][y][x] >> SHIFTQ) & QMASK;
							VelAnis[z][y][x] = (indexVpdt2 & VELMASK) | ((indexDen & DENMASK) << SHIFTDen) | ((indexQ & QMASK) << SHIFTQ);
						}
					}
				}

				// release denang, since it is not used by isotropic code.
				Free_Padded_Earth_Model(PadDenAng,nz,DenAng);
			}
		}
		else
		{
			// ISO
			Read_Single_Earth_Model_Attribute(logLevel,vpname,scalarVp,swapin,ubeg,uend,vbeg,vend,wbeg,wend,dimu,dimv,dimw,arg_fastAxis,arg_medAxis,arg_slowAxis,zstart,VpIn,vpminZ,vpmaxZ);

			vpminX = vpminZ;
			vpmaxX = vpmaxZ;

			// compute vminX, vmaxX, vminZ and vmaxZ.
			vpmin = (vpminZ<vpminX)?vpminZ:vpminX;
			vpmax = (vpmaxX>vpmaxZ)?vpmaxX:vpmaxZ;
			// joe: need to adjust minfac according to actual rotated vels and expanding dz
			minfac = (vpminX/dh < vpminZ/dz)?vpminX/dh:vpminZ/dz;
			fmax = khmax*minfac/(2.0f*M_PI);
			// dt = gam * ( (dz/vpmaxZ < dh/vpmaxX) ? dz/vpmaxZ : dh/vpmaxX );
			dt = gam * ( (dz/vpmax < dh/vpmax) ? dz/vpmax : dh/vpmax );

			//vpmaxZ += 0.1f;
			//vpmaxX += 0.1f;
			//Epsmax += 0.00001f;
			invVelrange = (vpminZ == vpmaxZ) ? 0.0f : 1.0f / (vpmaxZ-vpminZ);
			Velbinsize = (vpminZ == vpmaxZ) ? 0.0f : (vpmaxZ-vpminZ)/VELMASK;

			if (logLevel >= 4)
			{
				printf("vpminZ = %f, vpmaxZ = %f, invVelrange = %e, Velbinsize = %e\n",vpminZ,vpmaxZ,invVelrange,Velbinsize);
			}

			// add Vp, Eps to compressed earth model
#pragma omp parallel for schedule(dynamic)
			for(int z=zstart; z<nz; z++)
			{
				for(int y=0; y<ny; y++)
				{
					for(int x=0; x<actual_nx; x++)
					{ 
						// Pack VelAnis index array
						unsigned int indexVpdt2 = (unsigned int)(((VpIn[z][y][x]-vpminZ)*invVelrange*VELMASK) + 0.5f);

						VelAnis[z][y][x] = (indexVpdt2 & VELMASK);
					}
				}
			}

			// compute fmax. need dt for that.
			if     (newfmax >  1.1f) fmax  =  newfmax;
			else if(newfmax < -0.1f) fmax *= -newfmax;
			tmax = (int)(maxtime/dt) + 1;
			float fq = fmax / 3.0f;

			Read_Single_Earth_Model_Attribute(logLevel,dnname,scalarDen,swapin,ubeg,uend,vbeg,vend,wbeg,wend,dimu,dimv,dimw,arg_fastAxis,arg_medAxis,arg_slowAxis,zstart,VpIn,Denmin,Denmax);

			invDenrange = (Denmin == Denmax) ? 0.0f : 1.0f / (Denmax-Denmin);
			Denbinsize = (Denmin == Denmax) ? 0.0f : (Denmax-Denmin)/DENMASK;

			if (logLevel >= 4)
			{
				printf("Denmin = %f, Denmax = %f, invDenrange = %e, Denbinsize = %e\n",Denmin,Denmax,invDenrange,Denbinsize);
			}

#pragma omp parallel for schedule(dynamic)
			for(int z=zstart; z<nz; z++)
			{
				for(int y=0; y<ny; y++)
				{
					for(int x=0; x<actual_nx; x++)
					{ 
						// Pack DenAng index array
						unsigned int indexDen = (unsigned int)(((VpIn[z][y][x]-Denmin)*invDenrange*DENMASK) + 0.5f);

						VelAnis[z][y][x] |= ((indexDen & DENMASK) << SHIFTDen);
					}
				}
			}

			Read_Single_Earth_Model_Attribute(logLevel,Qname,scalarQ,swapin,ubeg,uend,vbeg,vend,wbeg,wend,dimu,dimv,dimw,arg_fastAxis,arg_medAxis,arg_slowAxis,zstart,VpIn);

			Qmin = 1e37f;
			Qmax = -1e37f;
#pragma omp parallel for schedule(dynamic)
			for(int z=zstart; z<nz; z++)
			{
				float _min, _max;
#pragma omp critical
				{
					_min = Qmin;
					_max = Qmax;
				}
				for(int y=0; y<ny; y++)
				{ 
					// calc min & max values
					for(int x=0; x<actual_nx; x++)
					{
						float val = VpIn[z][y][x];
						// transform Qatten into a scalar multiplier at each node location: this makes Q proportional to freq,
						// because only one dominant freq is used to represent all frequencies
						// Qatten = 1-(PI/Q)*dt/Tdom = 1-PI/Q*dt*fmax/3 = 1-0.333*PI*fmax*dt/Q (applied to amps, not energy)
						//val = (val > 1e8f) ? 1.0f : 1.0f - 0.333f*M_PI*fmax*dt/val;
						val = (val > 1e8f) ? 1.0f : exp(-M_PI*fq*dt/val);
						VpIn[z][y][x] = val;
						if (val > _max) _max=val;  
						if (val < _min) _min=val;
					}
				}
#pragma omp critical
				{
					if (_min < Qmin) Qmin = _min;
					if (_max > Qmax) Qmax = _max;
				}
			}

			if (Qmax > Qmin)
			{
				invQrange = (float)((double)QMASK / ((double)Qmax-(double)Qmin));
				Qbinsize = (float)(((double)Qmax-(double)Qmin)/(double)QMASK);
			}
			else
			{
				invQrange = Qbinsize = 0.0f;
			}

			if (logLevel >= 4)
			{
				printf("Qmin = %e, Qmax = %e, invQrange = %e, Qbinsize = %e\n",Qmin,Qmax,invQrange,Qbinsize);
			}

#pragma omp parallel for schedule(dynamic) 
			for(int z=zstart; z<nz; z++)
			{
				for(int y=0; y<ny; y++)
				{
					for(int x=0; x<actual_nx; x++)
					{ 
						// Pack DenAng index array
						unsigned int indexQ = (unsigned int)((VpIn[z][y][x]-Qmin)*invQrange + 0.5f);

						VelAnis[z][y][x] |= ((indexQ & QMASK) << SHIFTQ);
					}
				}
			}
		}
	}

	// Extend earth model in X if actual_nx differs from nx
	if (actual_nx < nx)
	{
		for (int z = zstart;  z < nz;  ++z)
		{
			for (int y = 0;  y < ny;  ++y)
			{
				for (int x = actual_nx;  x < nx;  ++x)
				{
					VelAnis[z][y][x] = VelAnis[z][y][actual_nx-1];
					if (KernelType > 0)
					{
						DenAng[z][y][x] = DenAng[z][y][actual_nx-1];
					}
				}
			}
		}
	}

	// Extend earth model if we have a skylayer
	if (zstart > 0)
	{
		for (int z = 0;  z < zstart;  ++z)
		{
			for (int y = 0;  y < ny;  ++y)
			{
				for (int x = 0;  x < nx;  ++x)
				{
					VelAnis[z][y][x] = VelAnis[zstart][y][x];
					if (KernelType > 0)
					{
						DenAng[z][y][x] = DenAng[zstart][y][x];
					}
				}
			}
		}
	}

	strcpy(hdrstring,"TTI3: TTI Pressure Finite Difference 3D in SEGY Format\n");
	sprintf(infostring,"dh=%.2f  dz=%.2f  maxtime=%.3f  vpminX=%.2f vpmaxX=%.2f vpminZ=%.2f vpmaxZ=%.2f\n",
			dh,dz,maxtime,vpminX,vpmaxX,vpminZ,vpmaxZ);
	strcat(hdrstring, infostring);
	fprintf(stdout,"%.2f < VpX < %.2f    %.2f < VpZ < %.2f\n", vpminX,vpmaxX-1.,vpminZ,vpmaxZ-1.);

	/* dt = dh*gam/Vmax = khmax*Vmin/(2pi*fmax) * gam/Vmax ~~ 0.111(Vmin/Vmax)/fmax
	   < 0.111/fmax < dtout < 0.5/Fnyquist <= 0.5/fmax  */
	if(dtout < dt)        dtout = alter_dtout( 1, dt);
	if(dtout > 0.5f/fmax) dtout = alter_dtout(-1, 0.5f/fmax);
	sprintf(infostring,"fmax=%.2f (default fmax=%.2f) dtinternal=%f  tmax=%d  dtexternal=%f\n",
			fmax, khmax*minfac/(2.0f*M_PI), dt,tmax,dtout);
	strcat(hdrstring, infostring);

	/* POPULATE THE VALUE TABLES */
	dt2 = dt*dt;  // Veltable: note that Vp has dt2 included
	for (int i = 0;  i <= VELMASK;  ++i)
	{
		float Vp = vpminZ + (float)i * Velbinsize;
		lutVp2[i] = Vp * Vp * dt2;
	}

	if (KernelType > 0)
	{
		for (int i = 0;  i <= EPSMASK;  ++i) lutEps[i] = Epsmin + (float)i * Epsbinsize;
		for (int i = 0;  i <= DELMASK;  ++i) lutDel[i] = Delmin + (float)i * Delbinsize;

		lutc44c33[0] = 0.0f;  lutc44c33[1] = VsoVp0*VsoVp0;
		// PosDef check: need C11*C33 > C13*C13 --> C44/C33 >~ -2*(eps-del)/(del*del)
	}

	for (int i = 0;  i <= DENMASK;  ++i)
	{
		float Den = Denmin + (float)i * Denbinsize;
		lutDen[i] = Den;
		lutBuoy[i] = 1.0f / Den;
	}

	if (KernelType > 1)
	{
		for (int i = 0;  i <= DIPMASK;  ++i)
		{
			float angle = Dipmin + (float)i * Dipbinsize;
			lutsDip[i] = sinf(angle);
			lutcDip[i] = cosf(angle);
		}

		for (int i = 0;  i <= AZMMASK;  ++i)
		{
			float angle = Azmmin + (float)i * Azmbinsize;
			lutsAzm[i] = sinf(angle);
			lutcAzm[i] = cosf(angle);
		}
	}

	for (int i=0;  i<=QMASK;  ++i)
	{
		lutQ[i] = Qmin + Qbinsize * (float)i;
	}

	Init_Earth_Model_Compression(
		vpminZ,Velbinsize,dt2,
		Epsmin,Epsbinsize,
		Delmin,Delbinsize,
		Denmin,Denbinsize,
		Dipmin,Dipbinsize,
		Azmmin,Azmbinsize,
		Qmin,Qbinsize,
		0.0f,VsoVp0*VsoVp0
		);

	/** R E C E I V E R  G E O M E T R Y **/
	if (logLevel >= 2) fprintf(stdout,"Assigning Geometry...\n");
	nyrec = (yrecend-yrecstart)/yrecstride + 1;
	yrecend = yrecstart + (nyrec-1)*yrecstride;
	if(yrecend >= ny) 
	{ 
		fprintf(stderr,"yrec too close to edge\n"); exit(0); 
	}

	nxrec = (xrecend-xrecstart)/xrecstride + 1;
	xrecend = xrecstart + (nxrec-1)*xrecstride;
	if(xrecend >= actual_nx) 
	{ 
		fprintf(stderr,"xrec too close to edge\n"); exit(0); 
	}

	nzrec = (zrecend-zrecstart)/zrecstride + 1;
	if(absorbz0) zrecstartgrid = zrecstart + skylayer;  else zrecstartgrid = zrecstart;
	zrecend = zrecstartgrid + (nzrec-1)*zrecstride;
	if(zrecend >= nz) 
	{ 
		fprintf(stderr,"zrec too close to edge\n"); exit(0); 
	}

	nrec = nxrec*nyrec*nzrec;

	/* Rec 3D float array */
	int nsampout = (int)ceil(maxtime/dtout) + 1;
#ifdef STANDALONE
	float* tmp_write_buf = (float*)malloc(nxrec*sizeof(float));
#else
	long recveclen = (long)nrec * (long)nsampout * (long)sizeof(float);
	recvec = (float *)malloc(recveclen);
	memset((void*)recvec, 0, recveclen);
	if (logLevel >= 2) printf("Allocated %ld MB for recvec\n",recveclen/(long)1048576);
#endif

	long recspongelen = (long)nrec * (long)2 * (long)sizeof(float);
	recsponge = (float*)malloc(recspongelen);
	memset((void*)recsponge, 0, recspongelen);
	if (logLevel >= 4) printf("Allocated %ld KB for recsponge\n",recspongelen/(long)1024);

	if(absorbz0)
	{ 
		zsrc += skylayer;
		if(srcghost)
		{ 
			zsrcghost = 2*skylayer - zsrc;
			if(zsrcghost<1) 
			{ 
				fprintf(stderr,"ghost src out of zbounds\n"); exit(0); 
			}
		}

		sprintf(infostring,"Absorbing Top Surface. \n");
		strcat(hdrstring, infostring);
		float fZSRC = (zsrc+arg_zbeg-skylayer)*dz;
		if (source_vert_interp > 0) fZSRC += (dz/2.0f);
		if(sourcetype==0 || sourcetype==1) sprintf(infostring,"srcx=%.2f  srcy=%.2f  srcz=%.2f\n",
				(xsrc+arg_xbeg)*dh, (ysrc+arg_ybeg)*dh, fZSRC);
		else              sprintf(infostring,"srcx=%.2f  srcy=%.2f  srcz=%.2f\n",
				srcx, srcy, srcz);
		strcat(hdrstring, infostring);
	}
	else
	{ 
		sprintf(infostring,"Free Reflecting Surface.\n");
		strcat(hdrstring, infostring);
		float fZSRC = (zsrc+arg_zbeg-skylayer)*dz;
		if (source_vert_interp > 0) fZSRC += (dz/2.0f);
		if(sourcetype==0 || sourcetype==1) sprintf(infostring,"srcx=%.2f  srcy=%.2f  srcz=%.2f\n",
				(xsrc+arg_xbeg)*dh, (ysrc+arg_ybeg)*dh, fZSRC);
		else              sprintf(infostring,"srcx=%.2f  srcy=%.2f  srcz=%.2f\n",
				srcx, srcy, srcz);
		strcat(hdrstring, infostring);
	}

	/*** C O M P U T E   S O U R C E   T I M E   F U N C T I O N ***/
	velanis = VelAnis[zsrc][ysrc][xsrc];
	velsrcZ = sqrtf(lutVp2[velanis & VELMASK]/dt2);
	if (KernelType > 0)
	{
		velsrcX = velsrcZ * sqrtf(1.0f + 2.0f*lutEps[(velanis >> SHIFTEps) & EPSMASK]);
	}
	else
	{
		velsrcX = velsrcZ;
	}
	if (logLevel >= 2)
	{
		int denang = KernelType > 0 ? DenAng[zsrc][ysrc][xsrc] : VelAnis[zsrc][ysrc][xsrc];
		fprintf(stdout,"after parm compression: VelsrcZ = %f   Densrc = %f\n", velsrcZ, 
				lutDen[(denang >> SHIFTDen)& DENMASK]);
	}

	src(logLevel, dt, fmax, sourcetype+1, stfname, &tsrc, stf);
	tstartsim = 0;
	tstartrec = (int)(timestartrec/dt);
	timestartrec = tstartrec*dt;
	srcboxH = srcboxZ = 0;

	if (logLevel >= 5)
	{
		printf("Source function computed\n");
		fflush(stdout);
	}

	/*** PUT TEMPORARY ISOTROPIC & Vs=0 HALO AROUND SOURCE ***/
	// preserve original source velanis values for later re-insertion
	velsrc  = (velsrcZ < velsrcX) ? velsrcZ : velsrcX; // min velsrc
	tvolsrc = (int)((isorad*dh/velsrc + 4.0f/fmax)/dt);
	VelAnisOrigsrc = imatrix3(2*isorad+1, 2*isorad+1, 2*isorad+1);

	velanisIsosrc = (int)(((velsrc-vpminZ)*invVelrange*VELMASK) + 0.5f);

	isozstart = (zsrc-isorad >= 0) ? zsrc-isorad : 0;
	isoystart = ysrc-isorad;
	isoxstart = xsrc-isorad;

	if (logLevel >= 5)
	{
		printf("PUT TEMPORARY ISOTROPIC & Vs=0 HALO AROUND SOURCE\n");
		fflush(stdout);
	}

	/*
	// check for negative eta in source region; if found, do not taper iso region
	for(z=isozstart; z<=zsrc+isorad; z++)
	for(y=isoystart; y<=ysrc+isorad; y++)
	for(x=isoxstart; x<=xsrc+isorad; x++)
	if(C11[z][y][x]*C33[z][y][x] < C44ref*(C11[z][y][x] + C33[z][y][x] + 2.0f*C13[z][y][x])
	+ C13[z][y][x]*C13[z][y][x]) isotapersrc = 0;
	if(isotapersrc)
	{ // (this block tapers the isotropic source region but requires eta nonnegative)
	isohalfrad = isorad/2;

	for(z=isozstart; z<=zsrc+isorad; z++)
	{ for(y=isoystart; y<=ysrc+isorad; y++)
	for(x=isoxstart; x<=xsrc+isorad; x++)
	{ C11origsrc[z-isozstart][y-isoystart][x-isoxstart] = C11[z][y][x];
	C13origsrc[z-isozstart][y-isoystart][x-isoxstart] = C13[z][y][x];

	isofac = (sqrtf((float)((x-xsrc)*(x-xsrc)+(y-ysrc)*(y-ysrc)+(z-zsrc)*(z-zsrc)))
	- isohalfrad) / (isorad - isohalfrad);
	if(isofac >=1.0f) continue;  // full aniso, modify nothing
	else if(isofac < 0.0f)  C11[z][y][x] = C13[z][y][x] = C33[z][y][x]; // pure iso
	else  // iso/aniso penumbra
	{ C11[z][y][x] = isofac*C11[z][y][x] + (1.0f-isofac)*C33[z][y][x];
	C13[z][y][x] = isofac*C13[z][y][x] + (1.0f-isofac)*C33[z][y][x];
	if(C13[z][y][x]*C13[z][y][x] > C11[z][y][x]*C33[z][y][x])
	C13[z][y][x] = sqrtf(C11[z][y][x]*C33[z][y][x]);
	}
	}
	}
	}
	else
	 */  // this next part OK
	if (isorad > 0)
	{ 
		for(z=isozstart; z<=zsrc+isorad; z++)
		{ 
			for(y=isoystart; y<=ysrc+isorad; y++)
			{
				for(x=isoxstart; x<=xsrc+isorad; x++)
				{ 
					VelAnisOrigsrc[z-isozstart][y-isoystart][x-isoxstart] = VelAnis[z][y][x];
					if(sqrtf((x-xsrc)*(x-xsrc)+(y-ysrc)*(y-ysrc)+(z-zsrc)*(z-zsrc)) < isorad)
						VelAnis[z][y][x] = velanisIsosrc; // pure iso
				}
			}
		}
	}

	/*** INITIALIZATION OF DYNAMIC VARIABLES TO "0" (faster with 1.e-20 than hard 0.) ***/
	/*
	for(z=0; z<nz; z++)
		for(y=0; y<ny; y++)
			for(x=0; x<nx; x++)
				p[z][y][x]  = q[z][y][x]  = r[z][y][x]  = s[z][y][x]  = Ap[z][y][x] = Aq[z][y][x] =
					V1[z][y][x] = V2[z][y][x] = V3[z][y][x] = V4[z][y][x] = V5[z][y][x] = zero;
	*/
	
	Compute_Stencils(OTflag,dh,dz);

	if (logLevel >= 5)
	{
		printf("Compute_Stencils\n");
		fflush(stdout);
	}

	posix_memalign((void**)&spgx, 64, (nx+ny+nz)*sizeof(float));
	spgy = spgx + nx;
	spgz = spgy + ny;
	Compute_Sponges(spongecoeff_x,spongecoeff_y,spongecoeff_z_lo,spongecoeff_z_hi,spongewidth_x,spongewidth_y,spongewidth_z_lo,spongewidth_z_hi,absorbz0,nx,num_x_zero_pad,ny,nz,spgx,spgy,spgz);

	if (logLevel >= 5)
	{
		printf("Compute_Sponges\n");
		fflush(stdout);
	}

	if (KernelType > 0) Precompute_Anisotropy_Parameters();
	Wipe_Propagation_Wavefields(KernelType,pq,rs,Apq,nx,ny,nz,xh,yh,zh);

	if (dump_xz_itvl > 0)
	{
		Write_XZ_Slice_EM(logLevel,"~/slices/em",PadDenAng,PadVelAnis,KernelType,nx,ny,nz,xh,yh,zh,ysrc);
	}

	if (logLevel >= 5)
	{
		printf("Wipe_Propagation_Wavefields\n");
		fflush(stdout);
	}

	// Get sponge factors for all receiver locations by calling ABCsponge.
	omp_setconst(KernelType,pq,dimxh,dimyh,dimzh,1.0f);
	ABCsponge(KernelType,pq,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,0,nx-1,0,ny-1,0,nz-1);
	if (dump_xz_itvl > 0)
	{
		char str2[256];
		sprintf(str2,"~/slices/sponge_y=%d.dat",ysrc);
		Write_XZ_Slice_P(logLevel,KernelType,str2,pq,nx,ny,nz,xh,yh,zh,ysrc);
	}
	saveSponge(KernelType,pq,skylayer);
	omp_setconst(KernelType,pq,dimxh,dimyh,dimzh,0.0f);

	if (logLevel >= 5)
	{
		printf("saveSponge\n");
		fflush(stdout);
	}

	int CurrTileSetSize = num_tile_shapes;
	int CurrTileIdx = 0;
	double mcells_full = (double)(nx1-1) * (double)(ny1-1) * (double)(nz-halfOperLen+1) * 1e-6;
	
	int ETA_Estimate_Cnt = 0;
	double acc_time = 0.0;

	// free up V1, V2 and V3 since they are no longer needed.
	for (int iZ = 0;  iZ < nz;  ++iZ)
	{
		// V* points to pq and rs, so no need to free Y's.
		if (V1 != 0L) free((void*)V1[iZ]);
		if (V2 != 0L) free((void*)V2[iZ]);
		if (V3 != 0L) free((void*)V3[iZ]);
	}
	if (V1 != 0L) free((void*)V1);
	if (V2 != 0L) free((void*)V2);
	if (V3 != 0L) free((void*)V3);

	if (logLevel >= 5)
	{
		printf("free up V1, V2 and V3 since they are no longer needed.\n");
		fflush(stdout);
	}

	float bsrcavg = 0.0f, bsrcavg_ghost = 0.0f;
	float *Kfac = 0L, *Kfac_ghost = 0L;
	if (source_vert_interp > 0)
	{
		float KOBNrec = Get_Bulk_Modulus(KernelType,VelAnis,DenAng,Denbinsize,Denmin,Velbinsize,vpminZ,nx/2,ny/2,zrecstartgrid);

		// store bulk modulus for vertical stencil
		// TO-DO: Absorbing z0 - how is ghost source handled?
		Kfac = new float[10];
		float bsrcavg = -((0.5f / Get_Density(KernelType,VelAnis,DenAng,Denbinsize,Denmin,xsrc,ysrc,zsrc)) + (0.5f / Get_Density(KernelType,VelAnis,DenAng,Denbinsize,Denmin,xsrc,ysrc,zsrc+1)));
		for (int i = -4;  i <= 5;  ++i)
		{
			int abs_zsrc = zsrc + i;
                        if (abs_zsrc < 0) abs_zsrc = -abs_zsrc;
			float bulk_modulus = Get_Bulk_Modulus(KernelType,VelAnis,DenAng,Denbinsize,Denmin,Velbinsize,vpminZ,xsrc,ysrc,abs_zsrc);
			Kfac[i+4] = Compute_Source_Vertical_Interpolation_Term(i,source_vert_interp,bsrcavg,bulk_modulus,KOBNrec);
		}
		
		if (absorbz0 && srcghost)
		{
			Kfac_ghost = new float[10];
			float bsrcavg_ghost = -((0.5f / Get_Density(KernelType,VelAnis,DenAng,Denbinsize,Denmin,xsrc,ysrc,zsrcghost)) + (0.5f / Get_Density(KernelType,VelAnis,DenAng,Denbinsize,Denmin,xsrc,ysrc,zsrcghost+1)));
			for (int i = -4;  i <= 5;  ++i)
			{
				int abs_zsrc = zsrcghost + i;
				if (abs_zsrc < 0) abs_zsrc = -abs_zsrc;
				float bulk_modulus = Get_Bulk_Modulus(KernelType,VelAnis,DenAng,Denbinsize,Denmin,Velbinsize,vpminZ,xsrc,ysrc,abs_zsrc);
				Kfac_ghost[i+4] = Compute_Source_Vertical_Interpolation_Term(i,source_vert_interp,bsrcavg_ghost,bulk_modulus,KOBNrec);
			}
		}
	}

	initDone = 1;

        if (num_obn_nodes > 0)
        {
                // verify obn node locations against density field and exit
                FILE* delta_z_file = fopen("delta_z", "w");
                FILE* den_z_file = fopen("den_z_node_locations.txt", "w");
                float min_delta_z = 1e37f;
                float max_delta_z = -1e37f;
                for (int i = 0;  i < num_obn_nodes;  ++i)
                {
                        float* p = obn_node_locations + i * 3;
                        float x = p[0];
                        float y = p[1];
                        float z = -p[2];

                        int ixsrc = (int)(x/dh+0.5f) - arg_xbeg;
                        int iysrc = (int)(y/dh+0.5f) - arg_ybeg;
                        int izsrc = 0;
                        int done = 0;
                        do
                        {
                                int denang = KernelType > 0 ? DenAng[izsrc+zstart][iysrc][ixsrc] : VelAnis[izsrc+zstart][iysrc][ixsrc];
                                int iden = (denang >> SHIFTDen) & DENMASK;
                                float den = (float)iden * Denbinsize + Denmin;
                                if (den < 1.1f)
                                {
                                        ++izsrc;
                                }
                                else
                                {
                                        done = 1;
                                }
                        } while (!done);
                        --izsrc;  // move back up since the cell we end up on is first cell of sediment.
                        float snap_z = floor(z/dz) * dz;
                        float den_z = (float)izsrc * dz;
                        printf("node %f %f : z = %f, snap_z = %f, den_z = %f\n",x,y,z,snap_z,den_z);

                        float delta_z = snap_z - den_z;
                        if (delta_z < min_delta_z) min_delta_z = delta_z;
                        if (delta_z > max_delta_z) max_delta_z = delta_z;

			float grid_x = (float)(ixsrc + arg_xbeg) * dh;		// grid_x = src_x after snapping to grid. Actual source location used in propagation.
			float grid_y = (float)(iysrc + arg_ybeg) * dh;		// grid_y = src_y after snapping to grid.

                        fprintf(delta_z_file, "%f %f %f\n", grid_x, grid_y, delta_z);
                        fprintf(den_z_file, "%f %f %f\n", grid_x, grid_y, -den_z);
                }
                fclose(delta_z_file);
                fclose(den_z_file);
                printf("delta_z = [%f, %f]. delta_z = z_from_node_file - z_from_density_field\n",min_delta_z,max_delta_z);
                exit(0);
        }

	/*** L O O P   O V E R   S I M U L A T I O N   T I M E ***/
	int FullVolume = 0;
	if (logLevel >= 1)
	{
		fprintf(stdout,"fmax=%.2f (default fmax=%.2f) dtinternal=%f  tmax=%d  tsrc=%d  tvolsrc=%d\n",
				fmax, khmax*minfac/(2.0f*M_PI), dt, tmax, tsrc, tvolsrc);
		fprintf(stdout,"Starting Time Loop...\n");
	}

	printf("Running %s Kernel!\n",(Supports_AVX ? "AVX" : "SSE"));

	int tdry = 0;
	float max_net_comp_throughput = 0.0f;
	double subvolsizeacc = 0.0;
	for(t=tstartsim, ntout=0, itout=0; ntout < nsampout; ++t)
	{
		if (!FullVolume)
		{
			/* EXPANDING BOX fractional (gam) cell per time step */
			if (expanding_box)
			{
				xstart = xsrc - srcboxH - 20 - gam*t; if(xstart < 1)   xstart=1;
				xend   = xsrc + srcboxH + 20 + gam*t; if(xend >= nx1) xend=nx1-1;
				ystart = ysrc - srcboxH - 20 - gam*t; if(ystart < 1)   ystart=1;
				yend   = ysrc + srcboxH + 20 + gam*t; if(yend >= ny1) yend=ny1-1;
				zend   = zsrc + srcboxZ + 20 + (int)(dh/dz*gam*t);
				if(zend > nz - halfOperLen) zend = nz - halfOperLen;
			}
			else
			{
				xstart = 1;
				xend = nx1-1;
				ystart = 1;
				yend = ny1-1;
				zend = nz - halfOperLen;
			}
			if ((xstart == 1) && (xend == nx1-1) && (ystart == 1) && (yend == ny1-1) && (zend == nz - halfOperLen))
			{
				FullVolume = 1;
				if (logLevel >= 1)
				{
					printf("EXPANDING BOX - Reached full volume!\n");
				}
			}
		}

		int bsIdx = -1;
		if (FullVolume)
		{
			if (CurrTileSetSize > 1)
			{
				bsIdx = pos[CurrTileIdx];
				Get_TileSize(bsIdx,bsX,bsY);
				if (logLevel >= 1)
				{
					if (CurrTileIdx == 0) printf("Profiling -- %d tile shapes remaining\n",CurrTileSetSize);
					if (logLevel >= 3) printf("%3d x %3d :: ",bsX,bsY);
				}
			}
		}

		struct timespec before, after;
		clock_gettime(CLOCK_REALTIME, &before);
		//fprintf(stderr,"%d ", t);

		/* REPLACE HALO WITH ORIGINAL PROPERTIES */
		if(isorad > 0 && t == tvolsrc)
		{
			for(z=isozstart; z<=zsrc+isorad; z++)
			{
				for(y=isoystart; y<=ysrc+isorad; y++)
				{
					for(x=isoxstart; x<=xsrc+isorad; x++)
					{
						VelAnis[z][y][x] = VelAnisOrigsrc[z-isozstart][y-isoystart][x-isoxstart];
					}
				}
			}
			//printf("Replaced earth model in range x=[%d,%d] y=[%d,%d] z=[%d,%d]\n",isoxstart,xsrc+isorad,isoystart,ysrc+isorad,isozstart,zsrc+isorad);
			//Write_XZ_Slice_EM("slices/em_orig",PadDenAng,PadVelAnis,nx,ny,nz,xh,yh,zh,ny/2);
		}

		/* ADD SOURCE TERM */
		if(t<tsrc && (sourcetype==0 || sourcetype ==1))
		{
			if (source_vert_interp == 0)
			{
				/* addSource(stf[t]); */
				//p[zsrc][ysrc][xsrc] += stf[t];
				//q[zsrc][ysrc][xsrc] += stf[t];
				AddToPQ(KernelType,pq,xsrc,ysrc,zsrc,stf[t],stf[t]);
				if (logLevel >= 5)
				{
					printf("Added %e to pq[%d,%d,%d]\n",stf[t],xsrc,ysrc,zsrc);
					fflush(stdout);
				}

				if(absorbz0 && srcghost)
				{ 
					//p[zsrcghost][ysrc][xsrc] -= stf[t];
					//q[zsrcghost][ysrc][xsrc] -= stf[t];
					AddToPQ(KernelType,pq,xsrc,ysrc,zsrcghost,-stf[t],-stf[t]);
					if (logLevel >= 5)
					{
						printf("Added %e to pq[%d,%d,%d]\n",-stf[t],xsrc,ysrc,zsrcghost);
						fflush(stdout);
					}
				}
			}
			else
			{
				for (int i = -4;  i <= 5;  ++i) 
				{
					int abs_zsrc = zsrc + i;
					float source_term = stf[t] * Kfac[i+4];
					AddToPQ(KernelType,pq,xsrc,ysrc,abs_zsrc,source_term,source_term);
				}
				if(absorbz0 && srcghost)
				{
					for (int i = -4;  i <= 5;  ++i)
					{
						int abs_zsrc = zsrcghost + i;
						float source_term = stf[t] * Kfac_ghost[i+4];
						AddToPQ(KernelType,pq,xsrc,ysrc,abs_zsrc,-source_term,-source_term);
					}
				}
			}
		}

		int skip_propagation = (is_dry_run == 0 || (is_dry_run == 3 && FullVolume)) ? 0 : 1;
		if (!skip_propagation)
		{
			if (KernelType == 0)
			{
				// ISO
				// No AVX iso kernel yet
				if (OTflag == 2)
				{
					ISODenQ_TimeStep(logLevel,0,pq,rs,pq,PadVelAnis,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
				}
				else
				{
					ISODenQ_TimeStep(logLevel,1,pq,rs,Apq,PadVelAnis,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
					ISODenQ_TimeStep(logLevel,2,pq,rs,Apq,PadVelAnis,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
				}
			}
			else if (KernelType == 1)
			{
				// VTI
				if (OTflag == 2)
				{
					if (Supports_AVX)
					{
						AVX_VTIDenQ_TimeStep(logLevel,0,pq,rs,pq,PadVelAnis,PadDenAng,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
					}
					else
					{
						VTIDenQ_TimeStep(logLevel,0,pq,rs,pq,PadVelAnis,PadDenAng,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
					}
				}
				else
				{
					if (Supports_AVX)
					{
						AVX_VTIDenQ_TimeStep(logLevel,1,pq,rs,Apq,PadVelAnis,PadDenAng,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
						AVX_VTIDenQ_TimeStep(logLevel,2,pq,rs,Apq,PadVelAnis,PadDenAng,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
					}
					else
					{
						VTIDenQ_TimeStep(logLevel,1,pq,rs,Apq,PadVelAnis,PadDenAng,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
						VTIDenQ_TimeStep(logLevel,2,pq,rs,Apq,PadVelAnis,PadDenAng,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
					}
				}
			}
			else if (KernelType == 2)
			{
				// TTI
				if (OTflag == 2)
				{
					if (Supports_AVX)
					{
						AVX_TTIDenQ_TimeStep(logLevel,0,pq,rs,pq,PadVelAnis,PadDenAng,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
					}
					else
					{
						TTIDenQ_TimeStep(logLevel,0,pq,rs,pq,PadVelAnis,PadDenAng,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
					}
				}
				else
				{
					if (Supports_AVX)
					{
						AVX_TTIDenQ_TimeStep(logLevel,1,pq,rs,Apq,PadVelAnis,PadDenAng,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
						AVX_TTIDenQ_TimeStep(logLevel,2,pq,rs,Apq,PadVelAnis,PadDenAng,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
					}
					else
					{
						TTIDenQ_TimeStep(logLevel,1,pq,rs,Apq,PadVelAnis,PadDenAng,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
						TTIDenQ_TimeStep(logLevel,2,pq,rs,Apq,PadVelAnis,PadDenAng,spgx,spgy,spgz,nx,ny,nz,xh,yh,zh,bsX,bsY,xstart,xend,ystart,yend,0,zend);
					}
				}
			}
		}

		// REMOVE SOURCE TERM FROM WHAT IS NOW PREVIOUS TIMESTEP
		if(t<tsrc && (sourcetype==0 || sourcetype ==1))
		{ 
			if (source_vert_interp == 0)
			{
				// addSource(stf[t]);
				//p[zsrc][ysrc][xsrc] += stf[t];
				//q[zsrc][ysrc][xsrc] += stf[t];
				AddToPQ(KernelType,pq,xsrc,ysrc,zsrc,-stf[t],-stf[t]);
				if (logLevel >= 5)
				{
					printf("Subtracted %e from pq[%d,%d,%d]\n",-stf[t],xsrc,ysrc,zsrc);
					fflush(stdout);
				}

				if(absorbz0 && srcghost)
				{ 
					//p[zsrcghost][ysrc][xsrc] -= stf[t];
					//q[zsrcghost][ysrc][xsrc] -= stf[t];
					AddToPQ(KernelType,pq,xsrc,ysrc,zsrcghost,stf[t],stf[t]);
					if (logLevel >= 5)
					{
						printf("Subtracted %e from pq[%d,%d,%d]\n",stf[t],xsrc,ysrc,zsrcghost);
						fflush(stdout);
					}
				}
			}
			else
			{
				for (int i = -4;  i <= 5;  ++i) 
				{
					int abs_zsrc = zsrc + i;
					float source_term = stf[t] * Kfac[i+4];
					AddToPQ(KernelType,pq,xsrc,ysrc,abs_zsrc,-source_term,-source_term);
				}
				if(absorbz0 && srcghost)
				{
					for (int i = -4;  i <= 5;  ++i)
					{
						int abs_zsrc = zsrcghost + i;
						float source_term = stf[t] * Kfac_ghost[i+4];
						AddToPQ(KernelType,pq,xsrc,ysrc,abs_zsrc,source_term,source_term);
					}
				}
			}
		}

		/* ABSORBING BOUNDARY CONDITIONS */
		//ABCsponge(pq,spongecoeff,spongewidth,absorbz0,nx,ny,nz,xh,yh,zh,xstart,xend,ystart,yend,zstart,zend);
		//ABCsponge();

		clock_gettime(CLOCK_REALTIME, &after);
		double elapsed_time = (double)after.tv_sec + (double)after.tv_nsec * 1e-9 - ((double)before.tv_sec + (double)before.tv_nsec * 1e-9);

		double mcells = (double)(xend-xstart+1) * (double)(yend-ystart+1) * (double)(zend+1) * 1e-6;
		double subvolsize = mcells/mcells_full;
		if (!FullVolume) subvolsizeacc += subvolsize;
		if (!skip_propagation)
		{	
			double overcomp, flops_per_cell, gflops, eff_freq, net_comp_throughput;
			TTIDenQ_Compute_Performance(Supports_AVX,mcells,elapsed_time,num_threads,bsX,bsY,KernelType,(OTflag==2)?0:3,overcomp,flops_per_cell,gflops,eff_freq,net_comp_throughput);
			if (net_comp_throughput > max_net_comp_throughput) max_net_comp_throughput = net_comp_throughput;
			if (logLevel >= 3)
			{
				if (FullVolume)
				{
					printf("Iter %03d/%03d;  %.2f sec - %.0f MCells/s - %.0f GFLOPS - %.2f GHz\n",t+1,tmax,elapsed_time,net_comp_throughput,gflops,eff_freq); 
				}
				else
				{
					printf("Iter %03d/%03d;  %.1f%% of full volume - %.2f sec - %.0f MCells/s - %.0f GFLOPS - %.2f GHz\n",t+1,tmax,100.0*subvolsize,elapsed_time,net_comp_throughput,gflops,eff_freq);
				}
			}
		}

		if (FullVolume)
		{
			if (is_dry_run == 3)
			{
				++tdry;
				if (tdry >= 3)
				{
					double expbox_fraction = (subvolsizeacc + (double)(tmax - t)) / (double)(tmax - tstartsim);
					double seconds_per_timestep = mcells / max_net_comp_throughput;
					total_runtime = seconds_per_timestep*(double)(tmax-tstartsim)*expbox_fraction;

					if (logLevel >= 1)
					{
						printf("\n");
						printf("In order to run this job, a total of %.2f GB of free memory is required.\n",(double)max_allocated_memory / 1073741824.0);
						printf("%d timestep are required, of whom %d are expanding box timesteps.\n",tmax,t-tstartsim);
						printf("Expanding box reduces computations by %.2f%%\n",100.0*(1.0-expbox_fraction));
						printf("Throughput is %.0f MCells/s.\n",max_net_comp_throughput);
						printf("Approximate run time is %.2f hours\n\n",total_runtime/3600.0);
					}

					return;
				}
			}


			if (CurrTileSetSize > 1)
			{
				avgPerf[bsIdx] += (mcells/elapsed_time);
				++CurrTileIdx;
				if (CurrTileIdx == CurrTileSetSize)
				{
					Sort_TileSize_Set(pos,avgPerf,CurrTileSetSize);
					int NewCurrTileSetSize = CurrTileSetSize >> 1;
					if (CurrTileSetSize > 4 && NewCurrTileSetSize <= 4)
					{
						// slow down tile shape set reduction
						CurrTileSetSize = 4;
					}
					else if (CurrTileSetSize > 1 && CurrTileSetSize <= 4)
					{
						--CurrTileSetSize;
					}
					else
					{
						CurrTileSetSize = NewCurrTileSetSize;
					}
					CurrTileIdx = 0;
				}
			}
			else
			{
				if (ETA_Estimate_Cnt > 20)
				{
					double avg_time = acc_time / (double)ETA_Estimate_Cnt;
					acc_time = elapsed_time;
					ETA_Estimate_Cnt = 1;

					double elapsed_so_far = (double)after.tv_sec + (double)after.tv_nsec * 1e-9 - ((double)glob_start.tv_sec + (double)glob_start.tv_nsec * 1e-9);
					double ETA = avg_time * (double)(tmax - t - 1);
					double TOT = elapsed_so_far + ETA;
					total_runtime = TOT;  // output to caller.
					
					int iETA = (int)ETA;
					int iETA_Hr = iETA / 3600;
					int iETA_Mn = (iETA-(3600*iETA_Hr))/60;
					int iETA_Sc = iETA-(3600*iETA_Hr)-(60*iETA_Mn);
					
					int iTOT = (int)TOT;
					int iTOT_Hr = iTOT / 3600;
					int iTOT_Mn = (iTOT-(3600*iTOT_Hr))/60;
					int iTOT_Sc = iTOT-(3600*iTOT_Hr)-(60*iTOT_Mn);

					if (logLevel >= 2)
					{
						if (iETA_Hr > 0)
						{
							printf("Estimated completion in %d hours %d minutes (total run time %d hours %d minutes)\n",iETA_Hr,iETA_Mn,iTOT_Hr,iTOT_Mn);
						}
						else
						{
							if (iTOT_Hr > 0)
							{
								printf("Estimated completion in %d minutes (total run time %d hours %d minutes)\n",iETA_Mn,iTOT_Hr,iTOT_Mn);
							}
							else
							{
								printf("Estimated completion in %d minutes (total run time %d minutes)\n",iETA_Mn,iTOT_Mn);
							}
						}
					}
				}
				else
				{
					acc_time += elapsed_time;
					++ETA_Estimate_Cnt;
				}
			}
			if (CurrTileSetSize == 1)
			{
				CurrTileSetSize = 0;
				Get_TileSize(pos[0],bsX,bsY);
				if (logLevel >= 1) printf("Optimal tile size is %3d x %3d.\n",bsX,bsY);
			}
		}

		if (is_dry_run == 0)
		{
			if (dump_xz_itvl > 0 && ((t%dump_xz_itvl) == 0))
			{
				//printf("Finished time step %d\n",t);
				char str1[256];
				sprintf(str1,"~/slices/slice_z=%d_t=%d.dat",zsrc,t);
				char str2[256];
				sprintf(str2,"~/slices/slice_y=%d_t=%d.dat",ysrc,t);
				char str3[256];
				sprintf(str3,"~/slices/slice_r_z=%d_t=%d.dat",zsrc,t);
				char str4[256];
				sprintf(str4,"~/slices/slice_r_y=%d_t=%d.dat",ysrc,t);

				//Range_P(pq,nx,ny,nz,xh,yh,zh);
				Write_XY_Slice_P(logLevel,KernelType,str1,rs,nx,ny,nz,xh,yh,zh,zsrc);
				Write_XZ_Slice_P(logLevel,KernelType,str2,rs,nx,ny,nz,xh,yh,zh,ysrc);
				Write_XY_Slice_P(logLevel,KernelType,str3,pq,nx,ny,nz,xh,yh,zh,zsrc);
				Write_XZ_Slice_P(logLevel,KernelType,str4,pq,nx,ny,nz,xh,yh,zh,ysrc);
			}

			/* SAVE SEISMIC DATA AT TIME INTERVAL DTOUT */
			//tfrac = ((t+1)*dt-tout)/dt; if(t==tstartsim) tfrac = 1.0f;

			float tout = tstartsim*dt + itout * dtout;
			float tfrac = t + 1.0f - tout/dt; 
			if (t == tstartsim) tfrac = 1.0f;
			if (tfrac > 0.0f && tfrac <= 1.0f)
			{ 
				if (t >= tstartrec)  
				{ 
#ifdef STANDALONE
					saveTimeslice(KernelType,tfrac,skylayer,fp_tmp,tmp_write_buf); 
#else
					saveTimeslice(KernelType,tfrac,skylayer); 
#endif
					//printf("NEW :: t = %d, tfrac = %f, tout = %f\n",t,tfrac,tout);
					++ntout; 
				}
				++itout;
			}
		}

		/* ADD SOURCE TERM */
		if(t<tsrc && (sourcetype==0 || sourcetype ==1))
		{
			if (source_vert_interp == 0)
			{
				/* addSource(stf[t]); */
				//p[zsrc][ysrc][xsrc] += stf[t];
				//q[zsrc][ysrc][xsrc] += stf[t];
				AddToPQ(KernelType,pq,xsrc,ysrc,zsrc,stf[t],stf[t]);
				if (logLevel >= 5)
				{
					printf("Added %e to pq[%d,%d,%d]\n",stf[t],xsrc,ysrc,zsrc);
					fflush(stdout);
				}

				if(absorbz0 && srcghost)
				{ 
					//p[zsrcghost][ysrc][xsrc] -= stf[t];
					//q[zsrcghost][ysrc][xsrc] -= stf[t];
					AddToPQ(KernelType,pq,xsrc,ysrc,zsrcghost,-stf[t],-stf[t]);
					if (logLevel >= 5)
					{
						printf("Added %e to pq[%d,%d,%d]\n",-stf[t],xsrc,ysrc,zsrcghost);
						fflush(stdout);
					}
				}
			}
			else
			{
				for (int i = -4;  i <= 5;  ++i) 
				{
					int abs_zsrc = zsrc + i;
					float source_term = stf[t] * Kfac[i+4];
					AddToPQ(KernelType,pq,xsrc,ysrc,abs_zsrc,source_term,source_term);
				}
				if(absorbz0 && srcghost)
				{
					for (int i = -4;  i <= 5;  ++i)
					{
						int abs_zsrc = zsrcghost + i;
						float source_term = stf[t] * Kfac_ghost[i+4];
						AddToPQ(KernelType,pq,xsrc,ysrc,abs_zsrc,-source_term,-source_term);
					}
				}
			}
		}

		/* SWAP POINTERS ( p<->r and q<->s ) */
		//swap = p;  p = r;  r = swap;
		//swap = q;  q = s;  s = swap;
		__m128* tmp = pq;
		pq = rs;
		rs = tmp;

	}  // END OF LOOP OVER TIME

#ifdef STANDALONE
	if (tmp_write_buf != 0L) free((void*)tmp_write_buf);	
#endif
	if (recsponge != 0L) free((void*)recsponge);

	if (Kfac != 0L) delete [] Kfac;
	if (Kfac_ghost != 0L) delete [] Kfac_ghost;

	if (lutVp2    != 0L) free((void*)lutVp2);
	if (lutEps    != 0L) free((void*)lutEps);
	if (lutDel    != 0L) free((void*)lutDel);
	if (lutDen    != 0L) free((void*)lutDen);
	if (lutBuoy   != 0L) free((void*)lutBuoy);
	if (lutsDip   != 0L) free((void*)lutsDip);
	if (lutcDip   != 0L) free((void*)lutcDip);
	if (lutsAzm   != 0L) free((void*)lutsAzm);
	if (lutcAzm   != 0L) free((void*)lutcAzm);
	if (lutQ      != 0L) free((void*)lutQ);
	if (lutc44c33 != 0L) free((void*)lutc44c33);

	if (pos != 0L) free((void*)pos);
	if (tile_shape != 0L) free((void*)tile_shape);
	if (avgPerf != 0L) free((void*)avgPerf);

	if (spgx != 0L) free((void*)spgx);
	// spgy and spgz are never allocated as separate arrays, they are piggy backed onto spgx.

	free((void*)pq);	pq = 0L;
	free((void*)rs);	rs = 0L;
	if (Apq != 0L)
	{
		free((void*)Apq);
		Apq = 0L;
	}

	Free_Padded_Earth_Model(PadVelAnis,nz,VelAnis);
	if (KernelType > 0) Free_Padded_Earth_Model(PadDenAng,nz,DenAng);

	_mm_setcsr(old_csr);

	struct timespec glob_end;
	clock_gettime(CLOCK_REALTIME, &glob_end);
	total_runtime = (double)glob_end.tv_sec + (double)glob_end.tv_nsec * 1e-9 - ((double)glob_start.tv_sec + (double)glob_start.tv_nsec * 1e-9);

#ifdef STANDALONE
	sprintf(infostring,"timestartrec=%.3f  nsamp=%d\n", timestartrec, ntout);
	strcat(hdrstring, infostring);

	if(shotsamp != nrec*ntout)
		fprintf(stderr,"Warning: shotsamp = %ld  #floats written to shotprof = %d\n",
				shotsamp,nx*nyrec*ntout);


	/*** CONVERT ENTIRE SHOT PROFILE TO SEGY (first free most memory)
	  AND WRITE TO DISK ***/

	if (logLevel >= 1) fprintf(stdout,"Flushing remainder of samples to tmp file...\n");
	fclose(fp_tmp);

	if (logLevel >= 1) fprintf(stdout,"Reading MUX'ed SEGY data from tmp file...\n");
	fp_tmp = fopen(tmpfilename, "r");
	long seisdataitems = (long)nrec * (long)nsampout;
	long seisdatalen = seisdataitems * (long)sizeof(float);
	float* seisdata = (float *)malloc(seisdatalen);
	fread(seisdata, sizeof(float), seisdataitems, fp_tmp);
	fclose(fp_tmp);
	remove(tmpfilename);  // clean up 

	if (logLevel >= 1) fprintf(stdout,"Writing SEGY Data...\n");
	float zsrc_segy = (absorbz0 ? (zsrc+arg_zbeg-skylayer)*dz : (zsrc+arg_zbeg)*dz);
	if (source_vert_interp > 0) zsrc_segy += 0.5f * dz;
	int leaky_integrate_flag = source_vert_interp == 1 ? 1 : 0;
	writeSEGY(
			seisdata, segyfilename, hdrstring, swapout, leaky_integrate_flag, id, 
			(xsrc+arg_xbeg)*dh, (ysrc+arg_ybeg)*dh, zsrc_segy, 
			dt, dtout, timestartrec, ntout, nrec,
			(xrecstart+arg_xbeg)*dh, nxrec, xrecstride*dh,
			(yrecstart+arg_ybeg)*dh, nyrec, yrecstride*dh,
			(zrecstart+arg_zbeg)*dz, nzrec, zrecstride*dz);

	/*** ALT: WRITE FLAT BINARY TO DISK ***/
	// swap4bytes((int*)recvec, shotsamp);
	// fwrite(recvec, sizeof(float), shotsamp, seisfile);
#endif
}  

#ifdef STANDALONE
extern void Call_Cmp_11pt_DDX(
        __m128* pq,
        __m128& ddx_01,
        __m128& ddx_23
        );
extern void Call_Cmp_11pt_DDX_DDX2(
        __m128* pq,
        __m128& ddx_01,
        __m128& ddx_23,
        __m128& ddx2_01,
        __m128& ddx2_23
        );

int main(int argc, char** argv)
{
	printf("\nTTI 3D FD Modeling tool. v1.9.1\n\n");

	/*
	__m128 zzA5A4A3 = _mm_setr_ps(0.0f, 5.0f, 4.0f, 3.0f);
	__m128 A2A1zzzz = _mm_setr_ps(2.0f, 1.0f, 0.0f, 0.0f);

	__m128 zz_A5 = _mm_shuffle_ps(zzA5A4A3, zzA5A4A3, 0x50);
	__m128 A5_A4 = _mm_shuffle_ps(zzA5A4A3, zzA5A4A3, 0xA5);
	__m128 A4_A3 = _mm_shuffle_ps(zzA5A4A3, zzA5A4A3, 0xFA);

	__m128 A3_A2 = _mm_shuffle_ps(zzA5A4A3, A2A1zzzz, 0x0F);
	__m128 A2_A1 = _mm_shuffle_ps(A2A1zzzz, A2A1zzzz, 0x50);
	__m128 A1_zz = _mm_shuffle_ps(A2A1zzzz, A2A1zzzz, 0xA5);

	__m128 A5_zz = _mm_shuffle_ps(zzA5A4A3, zzA5A4A3, 0x05);
	__m128 A4_A5 = _mm_shuffle_ps(zzA5A4A3, zzA5A4A3, 0x5A);
	__m128 A3_A4 = _mm_shuffle_ps(zzA5A4A3, zzA5A4A3, 0xAF);

	__m128 A2_A3 = _mm_shuffle_ps(A2A1zzzz ,zzA5A4A3, 0xF0);
	__m128 A1_A2 = _mm_shuffle_ps(A2A1zzzz, A2A1zzzz, 0x05);
	__m128 zz_A1 = _mm_shuffle_ps(A2A1zzzz, A2A1zzzz, 0x5A);

	float* f_zz_A5 = (float*)&zz_A5;
	float* f_A5_A4 = (float*)&A5_A4;
	float* f_A4_A3 = (float*)&A4_A3;

	float* f_A3_A2 = (float*)&A3_A2;
	float* f_A2_A1 = (float*)&A2_A1;
	float* f_A1_zz = (float*)&A1_zz;

	printf("zz_A5 = [%f, %f, %f, %f]\n",f_zz_A5[0], f_zz_A5[1], f_zz_A5[2], f_zz_A5[3]);
	printf("A5_A4 = [%f, %f, %f, %f]\n",f_A5_A4[0], f_A5_A4[1], f_A5_A4[2], f_A5_A4[3]);
	printf("A4_A3 = [%f, %f, %f, %f]\n",f_A4_A3[0], f_A4_A3[1], f_A4_A3[2], f_A4_A3[3]);

	printf("A3_A2 = [%f, %f, %f, %f]\n",f_A3_A2[0], f_A3_A2[1], f_A3_A2[2], f_A3_A2[3]);
	printf("A2_A1 = [%f, %f, %f, %f]\n",f_A2_A1[0], f_A2_A1[1], f_A2_A1[2], f_A2_A1[3]);
	printf("A1_zz = [%f, %f, %f, %f]\n",f_A1_zz[0], f_A1_zz[1], f_A1_zz[2], f_A1_zz[3]);

	float* f_A5_zz = (float*)&A5_zz;
	float* f_A4_A5 = (float*)&A4_A5;
	float* f_A3_A4 = (float*)&A3_A4;

	float* f_A2_A3 = (float*)&A2_A3;
	float* f_A1_A2 = (float*)&A1_A2;
	float* f_zz_A1 = (float*)&zz_A1;

	printf("A5_zz = [%f, %f, %f, %f]\n",f_A5_zz[0], f_A5_zz[1], f_A5_zz[2], f_A5_zz[3]);
	printf("A4_A5 = [%f, %f, %f, %f]\n",f_A4_A5[0], f_A4_A5[1], f_A4_A5[2], f_A4_A5[3]);
	printf("A3_A4 = [%f, %f, %f, %f]\n",f_A3_A4[0], f_A3_A4[1], f_A3_A4[2], f_A3_A4[3]);

	printf("A2_A3 = [%f, %f, %f, %f]\n",f_A2_A3[0], f_A2_A3[1], f_A2_A3[2], f_A2_A3[3]);
	printf("A1_A2 = [%f, %f, %f, %f]\n",f_A1_A2[0], f_A1_A2[1], f_A1_A2[2], f_A1_A2[3]);
	printf("zz_A1 = [%f, %f, %f, %f]\n",f_zz_A1[0], f_zz_A1[1], f_zz_A1[2], f_zz_A1[3]);

	Compute_Stencils(2,1.0f,1.0f);
	printf("Compute_Stencils done!\n");
	fflush(stdout);

	float *buffer1 = 0L, *buffer2 = 0L, *buffer3 = 0L;
	posix_memalign((void**)&buffer1, 64, 4096);
	posix_memalign((void**)&buffer2, 64, 4096);
	posix_memalign((void**)&buffer3, 64, 4096);
	memset((void*)buffer1, 0, 4096);
	memset((void*)buffer2, 0, 4096);
	memset((void*)buffer3, 0, 4096);
	buffer1[520] = 1.0f;
	printf("Buffers allocated.\n");
	fflush(stdout);

	for (int i = 2;  i < (1024-16)/8;  ++i)
	{
		__m128 ddx_01, ddx_23, ddx2_01, ddx2_23;
		Call_Cmp_11pt_DDX_DDX2((__m128*)&(buffer1[i*8]), ddx_01, ddx_23, ddx2_01, ddx2_23);
		_mm_store_ps(&(buffer2[i*8]), ddx_01);
		_mm_store_ps(&(buffer2[i*8+4]), ddx_23);
		_mm_store_ps(&(buffer3[i*8]), ddx2_01);
		_mm_store_ps(&(buffer3[i*8+4]), ddx2_23);
	}
	printf("Call_Cmp_11pt_DDX done!\n");
	fflush(stdout);

	for (int i = 0;  i < 1024;  ++i)
	{
		float val = buffer2[i];
		if (val != 0.0f)
		{
			printf("buffer2[%04d] = %f\n",i,val);
		}
	}
	for (int i = 0;  i < 1024;  ++i)
	{
		float val = buffer3[i];
		if (val != 0.0f)
		{
			printf("buffer3[%04d] = %f\n",i,val);
		}
	}
	printf("Printout done!\n");
	fflush(stdout);

	if (buffer3 != 0L) free(buffer3);
	if (buffer2 != 0L) free(buffer2);
	if (buffer1 != 0L) free(buffer1);
	exit(0);
	*/

	if(argc < 2 || argc > 6) 
	{ 
		printf("Usage: <<Modeling         >> tti3mod parmfile [num_threads] [dump_xz_itvl] [-info0|1|2|3] [-log0|1|2|3|4]\n"); 
		printf(" OR    <<QC Node Locations>> tti3mod parmfile node_locations_file.\n\n");
		printf("       node_locations: program will read node locations from file,\n");
		printf("       determine Z from density profile and output delta_Z and Z from density to two separate files.\n");
		printf("       These files are named 'delta_z' and 'den_z_node_locations.txt'.\n");
		printf("       num_threads is an optional argument. If provided, openmp will only use num_threads threads.\n");
		printf("       num_threads < 0 disables expanding box. This is meant for testing ONLY. \n\n");
		printf("       dump_xz_itvl is an optional argument. If provided, X-Z cross section will be dumped\n");
		printf("       to ASCII file every <dump_xz_itvl> shots. The format is suitable for viewing with gnuplot.\n");
		printf("       dump_xz_itvl must be greater than 1.\n\n");
		printf("       -info is an optional switch. If provided, the executable will load the parmfile and calculate\n");
		printf("       all run time parameters, but will not actually run propagation. Run time parameters are printed\n");
		printf("       out instead. These include how much memory is needed to run the job.\n");
		printf("       -log is an optional switch. Sets log level to 0, 1, 2, 3 or 4. Higher level means more verbosity.\n\n");
		exit(0); 
	}

	int argument_mode = 0;
        int num_obn_nodes = 0;
        float* obn_node_locations = 0L;
	if (argc == 3)
	{
		// check if first argument is parmfile and second argument is node locations file.
		FILE* obn_node_file = fopen(argv[2], "r");
		if (obn_node_file != 0L)
		{
			argument_mode = 1;

			// count number of lines to determine size of node locations array
			char buf[1024];
			num_obn_nodes = 0;
			int done = 0;
			do
			{
				char* str = fgets(buf, 1024, obn_node_file);
				if (str != 0L)
				{
					++num_obn_nodes;
				}
				else
				{
					done = 1;
				}
			}  while (!done);
			// allocate space
			printf("Node locations file has %d locations\n",num_obn_nodes);
			if (num_obn_nodes > 0)
			{
				obn_node_locations = new float[num_obn_nodes*3];
			}
			// rewind and read in nodes
			rewind(obn_node_file);
			for (int i = 0;  i < num_obn_nodes;  ++i)
			{
				char* str = fgets(buf, 1024, obn_node_file);
				float x,y,z,id;
				sscanf(str, "%f %f %f %f", &x, &y, &z, &id);
				printf("Node %d : %f %f %f\n",i+1,x,y,z);
				float* p = obn_node_locations + 3 * i;
				p[0] = x;
				p[1] = y;
				p[2] = z;
			}

			fclose(obn_node_file);
		}
	}

	int max_num_threads = 9999;
	int dump_xz_itvl = 0;
	int is_dry_run = 0;
	int logLevel = 3;

	if (argument_mode == 0)
	{
		if (argc > 2) sscanf(argv[2], "%d", &max_num_threads);

		if (argc > 3) sscanf(argv[3], "%d", &dump_xz_itvl);

		if (argc > 4)
		{
			if (strcmp(argv[4], "-info0") == 0) is_dry_run = 0;
			else if (strcmp(argv[4], "-info1") == 0) is_dry_run = 1;
			else if (strcmp(argv[4], "-info2") == 0) is_dry_run = 2;
			else if (strcmp(argv[4], "-info3") == 0) is_dry_run = 3;
			else is_dry_run = 0;
		}

		if (argc > 5)
		{
			if (strcmp(argv[5], "-log0") == 0) logLevel = 0;
			else if (strcmp(argv[5], "-log1") == 0) logLevel = 1;
			else if (strcmp(argv[5], "-log2") == 0) logLevel = 2;  // Normal
			else if (strcmp(argv[5], "-log3") == 0) logLevel = 3;  // Fine
			else if (strcmp(argv[5], "-log4") == 0) logLevel = 4;  // ALL
			else if (strcmp(argv[5], "-log5") == 0) logLevel = 5;  // Most debug messages
			else if (strcmp(argv[5], "-log6") == 0) logLevel = 6;  // All debug messages
		}
	}

	/*** O P E N  A N D  R E A D  P A R A M E T E R  F I L E ***/
	FILE* parmfile = fopen(argv[1], "r");
	if (parmfile == 0L)
	{
		fprintf(stderr,"Cannot find parmfile\n"); 
		exit(0); 
	}
	else
	{
		char seisname[256], vpname[256], epsetaname[256], deltaname[256];
		char dnname[256], dipdxname[256], azmdyname[256], Qname[256];
		char stfname[256];
		int id, etaflag, dipxdipyflag, degreesflag, swapin, smoothanisflag, isorad, arg_absorbz0, arg_srcghost, arg_recghost, arg_nx, arg_ny, arg_nz;
		int arg_sub_xoff, arg_sub_yoff, arg_sub_zoff, arg_sub_nx, arg_sub_ny, arg_sub_nz;
		float VsoVp0, arg_dh, arg_dz, arg_stretchfacz;
		float srcx, srcy, srcz, maxtime, timestartrec, dtout, newfmax, gamfac;
		int sourcetype, source_vert_interp, OTflag;
		int arg_xrecstart, arg_xrecend, arg_xrecstride;
		int arg_yrecstart, arg_yrecend, arg_yrecstride;
		int arg_zrecstart, arg_zrecend, arg_zrecstride;
		int arg_fastAxis, arg_medAxis, arg_slowAxis;
		int KernelType;	// 0->ISO, 1->VTI, 2->TTI
		int spongewidth_x, spongewidth_y, spongewidth_z_lo, spongewidth_z_hi;
		float spongeendval;

		fscanf(parmfile,"%s %d %d %d %d %d %s %s %d %s %f %s %s %s %d %d %d %s", seisname, &id, &KernelType, &arg_fastAxis, &arg_medAxis, &arg_slowAxis,
				vpname, epsetaname, &etaflag, deltaname, &VsoVp0, dnname, dipdxname, azmdyname, &dipxdipyflag, &degreesflag, &swapin, Qname);
		fscanf(parmfile,"%d %d %d %d %d", &smoothanisflag, &isorad, &arg_absorbz0, &arg_srcghost, &arg_recghost);
		fscanf(parmfile,"%d %f", &spongewidth_x, &spongeendval);
		fscanf(parmfile,"%f %f %f %d %d %d", &arg_dh, &arg_dz, &arg_stretchfacz, &arg_nx, &arg_ny, &arg_nz);
		fscanf(parmfile,"%d %d %d %d %d %d", &arg_sub_xoff, &arg_sub_yoff, &arg_sub_zoff, &arg_sub_nx, &arg_sub_ny, &arg_sub_nz);
		fscanf(parmfile,"%f %f %f %d %d %s", &srcx, &srcy, &srcz, &sourcetype, &source_vert_interp, stfname);
		fscanf(parmfile,"%d %d %d", &arg_xrecstart, &arg_xrecend, &arg_xrecstride);
		fscanf(parmfile,"%d %d %d", &arg_yrecstart, &arg_yrecend, &arg_yrecstride);
		fscanf(parmfile,"%d %d %d", &arg_zrecstart, &arg_zrecend, &arg_zrecstride);
		fscanf(parmfile,"%f %f %f %f %f %d", &maxtime, &timestartrec, &dtout,
				&newfmax, &gamfac, &OTflag /* &nthread */ );
		fclose(parmfile);

		spongewidth_y = spongewidth_x;
		spongewidth_z_hi = (int)(((float)spongewidth_x * arg_dh / arg_dz) + 0.5f);
		spongewidth_z_lo = arg_absorbz0 ? spongewidth_z_hi : spongewidth_z_hi;  // lo sponge may be longer than hi, although in this case isn't

		// compute sponge coefficient
		
		float spongecoeff_x = (1.0f - spongeendval) / (float)((spongewidth_x-1) * (spongewidth_x-1));
		float spongecoeff_y = (1.0f - spongeendval) / (float)((spongewidth_y-1) * (spongewidth_y-1));
		float spongecoeff_z_lo = (1.0f - spongeendval) / (float)((spongewidth_z_lo-1) * (spongewidth_z_lo-1));
		float spongecoeff_z_hi = (1.0f - spongeendval) / (float)((spongewidth_z_hi-1) * (spongewidth_z_hi-1));

		// Specifies how data is ordered in input files.
		// The default setup is X-Y-Z, which means that X is the fast axis and Z is the slow.
		// The parmfiles do not allow this ordering to be specified, so all jobs run via parmfiles assume X-Y-Z.
		// TMJ 07/09/12 - Added support for axis spec in parmfile.
		//arg_fastAxis = 0;
		//arg_medAxis = 1;
		//arg_slowAxis = 2;
		if (
			arg_fastAxis < 0 || arg_fastAxis > 2 || 
			arg_medAxis  < 0 || arg_medAxis  > 2 || 
			arg_slowAxis < 0 || arg_slowAxis > 2 ||
			(arg_fastAxis + arg_medAxis + arg_slowAxis) != 3
			)
		{
			printf("Invalid axis ordering : %d %d %d\nAborting!\n",arg_fastAxis,arg_medAxis,arg_slowAxis);
			return -1;
		}

		// stand alone always operates in local coordinates, so this parameter is always zero.
		float rotationAngleWorldToLocalDegree = 0.0f;

		unsigned long max_allocated_memory = 0L;
		double elapsed_time = 0.0;
		VarDenQ_ComputeShot(
                                num_obn_nodes, obn_node_locations, 
				max_num_threads,dump_xz_itvl,is_dry_run,max_allocated_memory,elapsed_time,KernelType,
				seisname,id,vpname,epsetaname,etaflag,deltaname,VsoVp0,dnname,dipdxname,azmdyname,rotationAngleWorldToLocalDegree,dipxdipyflag,degreesflag,swapin,Qname,
				0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,
				smoothanisflag,isorad,arg_absorbz0,arg_srcghost,arg_recghost,
				spongewidth_x,spongewidth_y,spongewidth_z_lo,spongewidth_z_hi,
				spongecoeff_x,spongecoeff_y,spongecoeff_z_lo,spongecoeff_z_hi,
				arg_dh,arg_dz,arg_stretchfacz,arg_nx,arg_ny,arg_nz,arg_fastAxis,arg_medAxis,arg_slowAxis,
				arg_sub_xoff, arg_sub_xoff + arg_sub_nx - 1, 
				arg_sub_yoff, arg_sub_yoff + arg_sub_ny - 1, 
				arg_sub_zoff, arg_sub_zoff + arg_sub_nz - 1,
				srcx,srcy,srcz,sourcetype,source_vert_interp,stfname,
				arg_xrecstart,arg_xrecend,arg_xrecstride,
				arg_yrecstart,arg_yrecend,arg_yrecstride,
				arg_zrecstart,arg_zrecend,arg_zrecstride,
				maxtime,timestartrec,dtout,newfmax,gamfac,OTflag,logLevel
				);
		return 0;
	}
	return -1;
} /*** E N D  O F  M A I N ***/
#endif

