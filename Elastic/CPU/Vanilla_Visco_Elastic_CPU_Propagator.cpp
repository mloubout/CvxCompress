#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <emmintrin.h>
#include <tmmintrin.h>
#include <immintrin.h>

#include "Vanilla_Visco_Elastic_CPU_Propagator.hxx"

//
// D E S I G N   C H O I C E S
//
// This is a CPU implementation of an elastic orthorhombic seismic forward modeling propagator.
// 8th order in space, 2nd order in time.
// 
// The design is a typical "pencil" design, the X-Y face is tiled. Tile size can vary, but X tile size must
// be a multiple of 8. Optimal tile size is hard to determine, generally speaking it is the largest tile
// where all the inputs will fit into cache. Tile size can be controlled with method Set_TileSize(...).
// Default tile size is 96 x 16. 
//
// The main routine is Timestep(...), which moves the wavefields one timestep forward in time.
// Timestep(...) has two main passes, the first one propagates stresses, the second propagates particle velocities.
// Each pass has the following nested loops:
//
// Loop over tiles
//   Loop over Z
//     Loop over Y dimension of tile
//       Loop over X dimension of tile
//
// The outermost loop over "tiles" is parallelized with openmp, so each tile is propagated by a separate thread.
// The tile is pushed from one end of the volume to the other by the loop over Z.
// The "X dimension of tile" loop is done in increments of 8, i.e. 8 consecutive cells are computed simultaneously.
// All computations are done with AVX and SSE intrinsics, i.e. explicit vectorization.
// Explicit prefetch instructions are inserted into the code to minimize global memory latency.
// The prefetching attemps to load all data into cache for the next iteration of the "Y dimension of tile" loop.
// 

float Vanilla_Visco_Elastic_CPU_Propagator::WLUP(float* wf, int x, int y, int z)
{
	if (x >= 0 && x < _nx && y >= 0 && y < _ny && z >= 0 && z < _nz)
		return wf[CmpIdx(x,y,z)];
	else
		return 0.0f;
}

void Vanilla_Visco_Elastic_CPU_Propagator::Unit_Test_Comp_DX_DY_DZ(
		int nx,
		int ny,
		int nz,
		float dx,
		float dy,
		float dz,
		int xoff,
		int yoff,
		int zoff
		)
{
	float inv_dx = 1.0f / dx;
	float inv_dy = 1.0f / dy;
	float inv_dz = 1.0f / dz;
	long idx = CmpIdx(nx/2,ny/2,nz/2);
	Get_TXY_ptr()[idx] = 1.0f;
	for (int z = 0;  z < nz;  ++z)
	{
		for (int y = 0;  y < ny;  ++y)
		{
			for (int x = 0;  x < nx;  x+=8)
			{
				__m256 d1dx, d1dy, d1dz;
				Compute_DX_DY_DZ(Get_TXY_ptr(),x,y,z,xoff,yoff,zoff,inv_dx,inv_dy,inv_dz,d1dx,d1dy,d1dz);
				*((__m256*)(Get_TXX_ptr() + CmpIdx(x,y,z))) = d1dx;
				*((__m256*)(Get_TYY_ptr() + CmpIdx(x,y,z))) = d1dy;
				*((__m256*)(Get_TZZ_ptr() + CmpIdx(x,y,z))) = d1dz;
			}
		}
	}
	for (int z = 0;  z < nz;  ++z)
        {
                for (int y = 0;  y < ny;  ++y)
                {
                        for (int x = 0;  x < nx;  ++x)
                        {
				// compute d1dx the slow way
				float d1dx = (
					_C0 * ( WLUP(Get_TXY_ptr(), x+xoff+1, y, z) - WLUP(Get_TXY_ptr(), x+xoff  , y, z) ) + 
					_C1 * ( WLUP(Get_TXY_ptr(), x+xoff+2, y, z) - WLUP(Get_TXY_ptr(), x+xoff-1, y, z) ) + 
					_C2 * ( WLUP(Get_TXY_ptr(), x+xoff+3, y, z) - WLUP(Get_TXY_ptr(), x+xoff-2, y, z) ) + 
					_C3 * ( WLUP(Get_TXY_ptr(), x+xoff+4, y, z) - WLUP(Get_TXY_ptr(), x+xoff-3, y, z) ) ) * inv_dx;
				float d1dy = (
					_C0 * ( WLUP(Get_TXY_ptr(), x, y+yoff+1, z) - WLUP(Get_TXY_ptr(), x, y+yoff  , z) ) + 
					_C1 * ( WLUP(Get_TXY_ptr(), x, y+yoff+2, z) - WLUP(Get_TXY_ptr(), x, y+yoff-1, z) ) + 
					_C2 * ( WLUP(Get_TXY_ptr(), x, y+yoff+3, z) - WLUP(Get_TXY_ptr(), x, y+yoff-2, z) ) + 
					_C3 * ( WLUP(Get_TXY_ptr(), x, y+yoff+4, z) - WLUP(Get_TXY_ptr(), x, y+yoff-3, z) ) ) * inv_dy;
				float d1dz = (
					_C0 * ( WLUP(Get_TXY_ptr(), x, y, z+zoff+1) - WLUP(Get_TXY_ptr(), x, y, z+zoff  ) ) + 
					_C1 * ( WLUP(Get_TXY_ptr(), x, y, z+zoff+2) - WLUP(Get_TXY_ptr(), x, y, z+zoff-1) ) + 
					_C2 * ( WLUP(Get_TXY_ptr(), x, y, z+zoff+3) - WLUP(Get_TXY_ptr(), x, y, z+zoff-2) ) + 
					_C3 * ( WLUP(Get_TXY_ptr(), x, y, z+zoff+4) - WLUP(Get_TXY_ptr(), x, y, z+zoff-3) ) ) * inv_dz;
				float d1dx_fast = WLUP(Get_TXX_ptr(), x, y, z);
				float d1dy_fast = WLUP(Get_TYY_ptr(), x, y, z);
				float d1dz_fast = WLUP(Get_TZZ_ptr(), x, y, z);
				if (d1dx != d1dx_fast)
				{
					printf("d1dx_SLOW[%d,%d,%d]=%e != d1dx_FAST=%e\n",x,y,z,d1dx,d1dx_fast);
				}
				if (d1dy != d1dy_fast)
				{
					printf("d1dy_SLOW[%d,%d,%d]=%e != d1dy_FAST=%e\n",x,y,z,d1dy,d1dy_fast);
				}
				if (d1dz != d1dz_fast)
				{
					printf("d1dz_SLOW[%d,%d,%d]=%e != d1dz_FAST=%e\n",x,y,z,d1dz,d1dz_fast);
				}
			}
		}
	}
}

void Vanilla_Visco_Elastic_CPU_Propagator::Run_Unit_Tests()
{
	const int nx = 400;
	const int ny = 400;
	const int nz = 400;
	const float dx = 50.0f;
	const float dy = 50.0f;
	const float dz = 50.0f;
	const float fq = 10.0f;

	Vanilla_Visco_Elastic_CPU_Propagator* prop = new Vanilla_Visco_Elastic_CPU_Propagator(nx,ny,nz,dx,dy,dz,fq);
	prop->Unit_Test_Comp_DX_DY_DZ(nx,ny,nz,dx,dy,dz,0,0,0);
	delete prop;

	prop = new Vanilla_Visco_Elastic_CPU_Propagator(nx,ny,nz,dx,dy,dz,fq);
	prop->Unit_Test_Comp_DX_DY_DZ(nx,ny,nz,dx,dy,dz,-1,-1,0);
	delete prop;

	prop = new Vanilla_Visco_Elastic_CPU_Propagator(nx,ny,nz,dx,dy,dz,fq);
	prop->Unit_Test_Comp_DX_DY_DZ(nx,ny,nz,dx,dy,dz,-1,0,-1);
	delete prop;
}

#define NUM_PAGES 1
void Vanilla_Visco_Elastic_CPU_Propagator::omp_memclear(void* dst, size_t len)
{
        size_t leni = len / 16;
        size_t nn = (len + NUM_PAGES*getpagesize()-1) / (NUM_PAGES*getpagesize());
        size_t One_Thread_Full = NUM_PAGES * getpagesize() / sizeof(__m128);
        //printf("omp_memclear len=%ld, leni=%ld, nn=%ld, One_Thread_Full=%ld\n",len,leni,nn,One_Thread_Full);
#pragma omp parallel for schedule(static)
        for (int i = 0;  i < nn;  ++i)
        {
                __m128 zero = _mm_set_ps(0.0f, 0.0f, 0.0f, 0.0f);
                size_t i0 = (size_t)i * One_Thread_Full;
                size_t in = leni - i0;
                if (in > One_Thread_Full) in = One_Thread_Full;
                __m128* d = (__m128*)dst + i0;
                for (int j = 0;  j < in;  ++j)
                {
                        _mm_stream_ps((float*)(d+j), zero);
                }
        }
}

void Vanilla_Visco_Elastic_CPU_Propagator::omp_memcpy(void* dst, void* src, size_t len)
{
        size_t leni = len / 16;
        size_t nn = (len + NUM_PAGES*getpagesize()-1) / (NUM_PAGES*getpagesize());
        size_t One_Thread_Full = NUM_PAGES * getpagesize() / sizeof(__m128);
#pragma omp parallel for schedule(static)
        for (int i = 0;  i < nn;  ++i)
        {
                size_t i0 = (size_t)i * One_Thread_Full;
                size_t in = leni - i0;
                if (in > One_Thread_Full) in = One_Thread_Full;
                __m128* d = (__m128*)dst + i0;
                __m128* s = (__m128*)src + i0;
                for (int j = 0;  j < in;  ++j)
                {
                        if ((j&3) == 0) _mm_prefetch((char*)(s+j+16),_MM_HINT_T0);
                        _mm_stream_ps((float*)(d+j),_mm_load_ps((float*)(s+j)));
                }
        }
}
#undef NUM_PAGES

Vanilla_Visco_Elastic_CPU_Propagator::Vanilla_Visco_Elastic_CPU_Propagator(
		int nx,
		int ny,
		int nz,
		float dx,
		float dy,
		float dz,
		float fq
		)
{
	// make nx a multiple of 8
	_nx = ((nx + 7) >> 3) << 3;
	_ny = ny;
	_nz = nz;

	_nabc_top = 0;
	_nabc_bot = 10;
	_nabc_sdx = 10;
	_nabc_sdy = 10;
	
	_stride_y = (long)_nx;
	_stride_z = (long)_ny * _stride_y;

	_mm256_stride_y = _stride_y / 8;
	_mm256_stride_z = _stride_z / 8;

	_dx = dx;
	_dy = dy;
	_dz = dz;

	_fq = fq;

	_c[0] = -_C3;
	_c[1] = -_C2;
	_c[2] = -_C1;
	_c[3] = -_C0;
	_c[4] = _C0;
	_c[5] = _C1;
	_c[6] = _C2;
	_c[7] = _C3;

	// default tile size is 128 by 32
	// user can overwrite this parameter.
	_bx = 128;
	_by = 32;

	_vpvert_avtop = 0.0f;
	_vpvert_avbot = 0.0f;

	// calculate size of one flattened array in bytes
	long size = (long)_nx * (long)_ny * (long)_nz * (long)sizeof(float);

	// align on 64 byte boundaries
	posix_memalign((void**)&_txx, 64, size);
	posix_memalign((void**)&_tyy, 64, size);
	posix_memalign((void**)&_tzz, 64, size);
	posix_memalign((void**)&_txy, 64, size);
	posix_memalign((void**)&_txz, 64, size);
	posix_memalign((void**)&_tyz, 64, size);

	posix_memalign((void**)&_vx, 64, size);
	posix_memalign((void**)&_vy, 64, size);
	posix_memalign((void**)&_vz, 64, size);
	posix_memalign((void**)&_sx, 64, size);
	posix_memalign((void**)&_sy, 64, size);
	posix_memalign((void**)&_sz, 64, size);

	posix_memalign((void**)&_Vp, 64, size);
	posix_memalign((void**)&_Vs, 64, size);
	posix_memalign((void**)&_Rho, 64, size);
	posix_memalign((void**)&_Qp, 64, size);

	// set all arrays to zero.
	// memory pages are assigned to the socket that touches them first,
	// using a multi-threaded method to zero out the pages initially ensures
	// the pages are evenly distributed among both sockets, which increases
	// memory performance.
	omp_memclear(_txx, size);
	omp_memclear(_tyy, size);
	omp_memclear(_tzz, size);
	omp_memclear(_txy, size);
	omp_memclear(_txz, size);
	omp_memclear(_tyz, size);

	omp_memclear(_vx, size);
	omp_memclear(_vy, size);
	omp_memclear(_vz, size);
	omp_memclear(_sx, size);
	omp_memclear(_sy, size);
	omp_memclear(_sz, size);

	omp_memclear(_Vp, size);
	omp_memclear(_Vs, size);
	omp_memclear(_Rho, size);
	omp_memclear(_Qp, size);
}

Vanilla_Visco_Elastic_CPU_Propagator::~Vanilla_Visco_Elastic_CPU_Propagator()
{
	if (_txx != 0L) free((void*)_txx);
	if (_tyy != 0L) free((void*)_tyy);
	if (_tzz != 0L) free((void*)_tzz);
	if (_txy != 0L) free((void*)_txy);
	if (_txz != 0L) free((void*)_txz);
	if (_tyz != 0L) free((void*)_tyz);

	if (_vx != 0L) free((void*)_vx);
	if (_vy != 0L) free((void*)_vy);
	if (_vz != 0L) free((void*)_vz);
	if (_sx != 0L) free((void*)_sx);
	if (_sy != 0L) free((void*)_sy);
	if (_sz != 0L) free((void*)_sz);

	if (_Vp != 0L) free((void*)_Vp);
	if (_Vs != 0L) free((void*)_Vs);
	if (_Rho != 0L) free((void*)_Rho);
	if (_Qp != 0L) free((void*)_Qp);
}

//
// Compute first order derivative along x axis for "wf" wavefield.
// This code uses SSE intrinsics for speed.
// 
inline void Vanilla_Visco_Elastic_CPU_Propagator::Compute_DX(
		float* wf,
		int x,
		int y,
		int z,
		int xoff,
		float inv_dx,
		__m256& dx
		)
{
	// we use sse intrinsics for this because of the palignr instruction
	// ..load a window of 16 consecutive floats containing the 8 cells we are interested in plus 4 halo cells on either side.
	__m128i* p = (__m128i*)(wf + CmpIdx(x-4,y,z));
	__m128i w0 = (x > 0) ? p[0] : _mm_castps_si128(_mm_setzero_ps());
	__m128i w1 = p[1];
	__m128i w2 = p[2];
	__m128i w3 = ((x+8) < _nx) ? p[3] : _mm_castps_si128(_mm_setzero_ps());
	_mm_prefetch((char*)(p+1+_stride_y), _MM_HINT_T0);
	// reverse order of floats in vector because of endian issue
	w0 = _mm_shuffle_epi32(w0, 0x1B);
	w1 = _mm_shuffle_epi32(w1, 0x1B);
	w2 = _mm_shuffle_epi32(w2, 0x1B);
	w3 = _mm_shuffle_epi32(w3, 0x1B);
	if (xoff == 0)
	{
		// shift 16 float window by 1 float to the left (4 bytes).
		// palignr shifts to the right, shifting 4 bytes to the left is the same as shifting 12 bytes to the right.
		w0 = _mm_alignr_epi8(w0, w1, 12);
		w1 = _mm_alignr_epi8(w1, w2, 12);
		w2 = _mm_alignr_epi8(w2, w3, 12);
		w3 = _mm_alignr_epi8(w3, w3, 12);
	}
	__m128 acc0 = _mm_mul_ps(_mm_set1_ps(_c[0]), _mm_castsi128_ps(w0));
	__m128 acc1 = _mm_mul_ps(_mm_set1_ps(_c[0]), _mm_castsi128_ps(w1));
	for (int i = 1;  i < 8;  ++i)
	{
		// shift 16 float window by 1 float (4 bytes).
		w0 = _mm_alignr_epi8(w0, w1, 12);
		w1 = _mm_alignr_epi8(w1, w2, 12);
		w2 = _mm_alignr_epi8(w2, w3, 12);
		w3 = _mm_alignr_epi8(w3, w3, 12);
		acc0 = _mm_add_ps(acc0, _mm_mul_ps(_mm_set1_ps(_c[i]), _mm_castsi128_ps(w0)));
		acc1 = _mm_add_ps(acc1, _mm_mul_ps(_mm_set1_ps(_c[i]), _mm_castsi128_ps(w1)));
	}
	//acc0 = _mm_shuffle_ps(acc0, acc0, 0x1B);
	//acc1 = _mm_shuffle_ps(acc1, acc1, 0x1B);
	// concatenate the two 4 float vectors into one 8 float vector and output
	dx = _mm256_insertf128_ps(_mm256_castps128_ps256(acc0), acc1, 1);
	dx = _mm256_shuffle_ps(dx, dx, 0x1B);
	dx = _mm256_mul_ps(dx, _mm256_set1_ps(inv_dx));
}

//
// Compute first order derivative along Y axis for "wf" wavefield.
// This code uses AVX intrinsics for speed.
// 
inline void Vanilla_Visco_Elastic_CPU_Propagator::Compute_DY(
		float* wf,
		int x,
		int y,
		int z,
		int yoff,
		float inv_dy,
		__m256& dy
		)
{
	int y4dy = y + yoff;
	__m256* q = (__m256*)(wf + CmpIdx(x,y4dy,z));
	__m256 yval_m3 = (y4dy-3) >= 0 ? q[-3*_mm256_stride_y] : _mm256_setzero_ps();
	__m256 yval_p4 = (y4dy+4) < _ny ? q[4*_mm256_stride_y] : _mm256_setzero_ps();
	dy = _mm256_mul_ps(_mm256_set1_ps(_C3), _mm256_sub_ps(yval_p4, yval_m3));
	__m256 yval_m2 = (y4dy-2) >= 0 ? q[-2*_mm256_stride_y] : _mm256_setzero_ps();
	__m256 yval_p3 = (y4dy+3) < _ny ? q[3*_mm256_stride_y] : _mm256_setzero_ps();
	dy = _mm256_add_ps(dy, _mm256_mul_ps(_mm256_set1_ps(_C2), _mm256_sub_ps(yval_p3, yval_m2)));
	__m256 yval_m1 = (y4dy-1) >= 0 ? q[-1*_mm256_stride_y] : _mm256_setzero_ps();
	__m256 yval_p2 = (y4dy+2) < _ny ? q[2*_mm256_stride_y] : _mm256_setzero_ps();
	dy = _mm256_add_ps(dy, _mm256_mul_ps(_mm256_set1_ps(_C1), _mm256_sub_ps(yval_p2, yval_m1)));
	__m256 yval_p0 = y4dy >= 0 ? q[0] : _mm256_setzero_ps();
	__m256 yval_p1 = (y4dy+1) < _ny ? q[1*_mm256_stride_y] : _mm256_setzero_ps();
	dy = _mm256_add_ps(dy, _mm256_mul_ps(_mm256_set1_ps(_C0), _mm256_sub_ps(yval_p1, yval_p0)));
	dy = _mm256_mul_ps(dy, _mm256_set1_ps(inv_dy));
	_mm_prefetch((char*)(q+5*_mm256_stride_y), _MM_HINT_T0);
}

//
// Compute first order derivative along Z axis for "wf" wavefield.
// This code uses AVX intrinsics for speed.
// 
inline void Vanilla_Visco_Elastic_CPU_Propagator::Compute_DZ(
		float* wf,
		int x,
		int y,
		int z,
		int zoff,
		float inv_dz,
		__m256& dz
		)
{
	int z4dz = z + zoff;
	__m256* q = (__m256*)(wf + CmpIdx(x,y,z4dz));
	__m256 zval_m3 = (z4dz-3) >= 0 ? q[-3*_mm256_stride_z] : _mm256_setzero_ps();
	__m256 zval_p4 = (z4dz+4) < _nz ? q[4*_mm256_stride_z] : _mm256_setzero_ps();
	dz = _mm256_mul_ps(_mm256_set1_ps(_C3), _mm256_sub_ps(zval_p4, zval_m3));
	__m256 zval_m2 = (z4dz-2) >= 0 ? q[-2*_mm256_stride_z] : _mm256_setzero_ps();
	__m256 zval_p3 = (z4dz+3) < _nz ? q[3*_mm256_stride_z] : _mm256_setzero_ps();
	dz = _mm256_add_ps(dz, _mm256_mul_ps(_mm256_set1_ps(_C2), _mm256_sub_ps(zval_p3, zval_m2)));
	__m256 zval_m1 = (z4dz-1) >= 0 ? q[-1*_mm256_stride_z] : _mm256_setzero_ps();
	__m256 zval_p2 = (z4dz+2) < _nz ? q[2*_mm256_stride_z] : _mm256_setzero_ps();
	dz = _mm256_add_ps(dz, _mm256_mul_ps(_mm256_set1_ps(_C1), _mm256_sub_ps(zval_p2, zval_m1)));
	__m256 zval_p0 = z4dz >= 0 ? q[0] : _mm256_setzero_ps();
	__m256 zval_p1 = (z4dz+1) < _nz ? q[1*_mm256_stride_z] : _mm256_setzero_ps();
	dz = _mm256_add_ps(dz, _mm256_mul_ps(_mm256_set1_ps(_C0), _mm256_sub_ps(zval_p1, zval_p0)));
	dz = _mm256_mul_ps(dz, _mm256_set1_ps(-inv_dz));  // minus sign is for positive Z axis up.
	_mm_prefetch((char*)(q+4*_mm256_stride_z+_mm256_stride_y), _MM_HINT_T0);
}

//
// Computes first order derivative along X, Y and Z axis for "wf" wavefield.
// 
inline void Vanilla_Visco_Elastic_CPU_Propagator::Compute_DX_DY_DZ(
		float* wf,
		int x,
		int y,
		int z,
		int xoff,
		int yoff,
		int zoff,
		float inv_dx,
		float inv_dy,
		float inv_dz,
		__m256& dx,
		__m256& dy,
		__m256& dz
		)
{
	Compute_DX(wf,x,y,z,xoff,inv_dx,dx);
	Compute_DY(wf,x,y,z,yoff,inv_dy,dy);
	Compute_DZ(wf,x,y,z,zoff,inv_dz,dz);
}

void Vanilla_Visco_Elastic_CPU_Propagator::Compute_On_Kite_CIJs(
		int x,
		int y,
		int z,
		__m256& c11,
		__m256& c22,
		__m256& c33,
		__m256& c44,
		__m256& c55,
		__m256& c66,
		__m256& c12,
		__m256& c13,
		__m256& c23
		)
{
	long CurrIdx = CmpIdx(x,y,z);

	__m256 Vp = *((__m256*)(Get_Vp_ptr() + CurrIdx));
	__m256 Vs = *((__m256*)(Get_Vs_ptr() + CurrIdx));
	__m256 Density = *((__m256*)(Get_Rho_ptr() + CurrIdx));
	__m256 Q = *((__m256*)(Get_Q_ptr() + CurrIdx));
	__m256 Epsilon1 = _mm256_setzero_ps();
	__m256 Epsilon2 = _mm256_setzero_ps();
	__m256 Gamma1 = _mm256_setzero_ps();
	__m256 Gamma2 = _mm256_setzero_ps();
	__m256 Delta1 = _mm256_setzero_ps();
	__m256 Delta2 = _mm256_setzero_ps();
	__m256 Delta3 = _mm256_setzero_ps();
	_mm_prefetch((char*)(Get_Vp_ptr() + CurrIdx + _stride_y), _MM_HINT_NTA);
	_mm_prefetch((char*)(Get_Vs_ptr() + CurrIdx + _stride_y), _MM_HINT_NTA);
	_mm_prefetch((char*)(Get_Rho_ptr() + CurrIdx + _stride_y), _MM_HINT_NTA);
	_mm_prefetch((char*)(Get_Q_ptr() + CurrIdx + _stride_y), _MM_HINT_NTA);

	c33 = _mm256_mul_ps(_mm256_mul_ps(Vp, Vp), Density);
	c11 = _mm256_mul_ps(c33, _mm256_add_ps(_mm256_set1_ps(1.0f), _mm256_add_ps(Epsilon2, Epsilon2)));
	c22 = _mm256_mul_ps(c33, _mm256_add_ps(_mm256_set1_ps(1.0f), _mm256_add_ps(Epsilon1, Epsilon1)));
	c55 = _mm256_mul_ps(_mm256_mul_ps(Vs, Vs), Density);
	c66 = _mm256_mul_ps(c55, _mm256_add_ps(_mm256_set1_ps(1.0f), _mm256_add_ps(Gamma1, Gamma1)));
	c44 = _mm256_div_ps(c66, _mm256_add_ps(_mm256_set1_ps(1.0f), _mm256_add_ps(Gamma2, Gamma2)));
	c12 = _mm256_sub_ps(_mm256_sqrt_ps(_mm256_mul_ps(_mm256_sub_ps(c11, c66), _mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(2.0f), _mm256_mul_ps(Delta3, c11)), _mm256_sub_ps(c11, c66)))), c66);
	c13 = _mm256_sub_ps(_mm256_sqrt_ps(_mm256_mul_ps(_mm256_sub_ps(c33, c55), _mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(2.0f), _mm256_mul_ps(Delta2, c33)), _mm256_sub_ps(c33, c55)))), c55);
	c23 = _mm256_sub_ps(_mm256_sqrt_ps(_mm256_mul_ps(_mm256_sub_ps(c33, c44), _mm256_add_ps(_mm256_mul_ps(_mm256_set1_ps(2.0f), _mm256_mul_ps(Delta1, c33)), _mm256_sub_ps(c33, c44)))), c44);
}

float Vanilla_Visco_Elastic_CPU_Propagator::_Compute_ABC(
        int x,
        int y,
        int z,
        int vol_nx,
        int vol_ny,
        int vol_nz,
        int nabc_top,
        int nabc_bot,
        int nabc_sdx,
        int nabc_sdy,
        float vpvert_avtop,
        float vpvert_avbot,
        float inv_DX,
        float inv_DY,
        float inv_DZ
        )
{
        bool Is_West = x < nabc_sdx;
        bool Is_East = x >= vol_nx - nabc_sdx;

        bool Is_South = y < nabc_sdy;
        bool Is_North = y >= vol_ny - nabc_sdy;

        bool Is_Top = z < nabc_top;
        bool Is_Bot = z >= vol_nz - nabc_bot;

        bool Is_X = Is_West || Is_East;
        bool Is_Y = Is_South || Is_North;
        bool Is_Z = Is_Top || Is_Bot;

        if ( !Is_X && !Is_Y && !Is_Z )
        {
		// no dampening
		return 0.0f;
        }
        else
        {
                if (Is_Z)
                {
			float zr = Is_Top ? (float)(nabc_top - z - 1) / (float)nabc_top : (float)(nabc_bot - vol_nz + z) / (float)nabc_bot;
			float deta_max = Is_Top ? vpvert_avtop * 18.0f * inv_DZ / (float)nabc_top : vpvert_avbot * 18.0f * inv_DZ / (float)nabc_bot;
			if (Is_X && Is_Y)
			{
				// cube
				float xr = Is_West ? (float)(nabc_sdx - x - 1) / (float)nabc_sdx : (float)(nabc_sdx - vol_nx + x) / (float)nabc_sdx;
				float yr = Is_South ? (float)(nabc_sdy - y - 1) / (float)nabc_sdy : (float)(nabc_sdy - vol_ny + y) / (float)nabc_sdy;
				xr = xr * xr * xr;
				xr = xr * xr;
				yr = yr * yr * yr;
				yr = yr * yr;
				zr = zr * zr * zr;
				zr = zr * zr;
				return powf(deta_max*deta_max*(xr+yr+zr),1.0f/3.0f);
			}
			else if (Is_Y)
			{
				// south or north beam
				float yr = Is_South ? (float)(nabc_sdy - y - 1) / (float)nabc_sdy : (float)(nabc_sdy - vol_ny + y) / (float)nabc_sdy;
				yr = yr * yr;
				yr = yr * yr;
				zr = zr * zr;
				zr = zr * zr;
				return sqrtf(deta_max*deta_max*(yr+zr));
			}
			else if (Is_X)
			{
				// west or east beam
				float xr = Is_West ? (float)(nabc_sdx - x - 1) / (float)nabc_sdx : (float)(nabc_sdx - vol_nx + x) / (float)nabc_sdx;
				xr = xr * xr;
				xr = xr * xr;
				zr = zr * zr;
				zr = zr * zr;
				return sqrtf(deta_max*deta_max*(xr+zr));
			}
			else
			{
				// top or bottom plate
				return deta_max*zr*zr;
			}
		}
		else if (Is_Y) // south or north plate
		{
			float yr = Is_South ? (float)(nabc_sdy - y - 1) / (float)nabc_sdy : (float)(nabc_sdy - vol_ny + y) / (float)nabc_sdy;
			float deta_maxsdy_top = vpvert_avtop * 18.0f * inv_DY / (float)nabc_sdy;
			float deta_maxsdy_bot = vpvert_avbot * 18.0f * inv_DY / (float)nabc_sdy;
			float deta_grady = (deta_maxsdy_bot - deta_maxsdy_top) / (float)(vol_nz-1+nabc_bot);
			float deta_maxsdy = deta_maxsdy_top + deta_grady * (float)(z-nabc_top);
			return deta_maxsdy*yr*yr;
		}
		else // west or east plate
		{
			float xr = Is_West ? (float)(nabc_sdx - x - 1) / (float)nabc_sdx : (float)(nabc_sdx - vol_nx + x) / (float)nabc_sdx;
			float deta_maxsdx_top = vpvert_avtop * 18.0f * inv_DX / (float)nabc_sdx;
			float deta_maxsdx_bot = vpvert_avbot * 18.0f * inv_DX / (float)nabc_sdx;
			float deta_gradx = (deta_maxsdx_bot - deta_maxsdx_top) / (float)(vol_nz-1+nabc_bot);
			float deta_maxsdx = deta_maxsdx_top + deta_gradx * (float)(z-nabc_top);
			return deta_maxsdx*xr*xr;
		}
	}
}

__m256 Vanilla_Visco_Elastic_CPU_Propagator::Compute_ABC(
        int x,
        int y,
        int z,
        int vol_nx,
        int vol_ny,
        int vol_nz,
        int nabc_top,
        int nabc_bot,
        int nabc_sdx,
        int nabc_sdy,
        float vpvert_avtop,
        float vpvert_avbot,
        float inv_DX,
        float inv_DY,
        float inv_DZ
        )
{
	return _mm256_set_ps(
		_Compute_ABC(x+7,y,z,vol_nx,vol_ny,vol_nz,nabc_top,nabc_bot,nabc_sdx,nabc_sdy,vpvert_avtop,vpvert_avbot,inv_DX,inv_DY,inv_DZ),
		_Compute_ABC(x+6,y,z,vol_nx,vol_ny,vol_nz,nabc_top,nabc_bot,nabc_sdx,nabc_sdy,vpvert_avtop,vpvert_avbot,inv_DX,inv_DY,inv_DZ),
		_Compute_ABC(x+5,y,z,vol_nx,vol_ny,vol_nz,nabc_top,nabc_bot,nabc_sdx,nabc_sdy,vpvert_avtop,vpvert_avbot,inv_DX,inv_DY,inv_DZ),
		_Compute_ABC(x+4,y,z,vol_nx,vol_ny,vol_nz,nabc_top,nabc_bot,nabc_sdx,nabc_sdy,vpvert_avtop,vpvert_avbot,inv_DX,inv_DY,inv_DZ),
		_Compute_ABC(x+3,y,z,vol_nx,vol_ny,vol_nz,nabc_top,nabc_bot,nabc_sdx,nabc_sdy,vpvert_avtop,vpvert_avbot,inv_DX,inv_DY,inv_DZ),
		_Compute_ABC(x+2,y,z,vol_nx,vol_ny,vol_nz,nabc_top,nabc_bot,nabc_sdx,nabc_sdy,vpvert_avtop,vpvert_avbot,inv_DX,inv_DY,inv_DZ),
		_Compute_ABC(x+1,y,z,vol_nx,vol_ny,vol_nz,nabc_top,nabc_bot,nabc_sdx,nabc_sdy,vpvert_avtop,vpvert_avbot,inv_DX,inv_DY,inv_DZ),
		_Compute_ABC(x,y,z,vol_nx,vol_ny,vol_nz,nabc_top,nabc_bot,nabc_sdx,nabc_sdy,vpvert_avtop,vpvert_avbot,inv_DX,inv_DY,inv_DZ)
		);
}

void Vanilla_Visco_Elastic_CPU_Propagator::Propagate_Stresses_Orthorhombic_Tile(
	float dti,
	int x0,
	int y0,
	int iz,
	int bx,
	int by,
	float inv_dx,
	float inv_dy,
	float inv_dz
	)
{
	for (int iy = 0;  iy < by;  ++iy)
	{
		for (int ix = 0;  ix < bx;  ix+=8)
		{
			Propagate_Stresses_Orthorhombic_8x(dti,ix+x0,iy+y0,iz,inv_dx,inv_dy,inv_dz);
		}
	}
}

void Vanilla_Visco_Elastic_CPU_Propagator::Propagate_Stresses_Orthorhombic_8x(
		float dti,
		int x,
		int y,
		int z,
		float inv_dx,
		float inv_dy,
		float inv_dz
		)
{
	//
	// compute derivatives
	//
	__m256 dxVx, dyVx, dzVx;
	Compute_DX_DY_DZ(Get_Vx_ptr(), x, y, z,  0,  0,  0, inv_dx, inv_dy, inv_dz, dxVx, dyVx, dzVx);

	__m256 dxVy, dyVy, dzVy;
	Compute_DX_DY_DZ(Get_Vy_ptr(), x, y, z, -1, -1,  0, inv_dx, inv_dy, inv_dz, dxVy, dyVy, dzVy);

	__m256 dxVz, dyVz, dzVz;
	Compute_DX_DY_DZ(Get_Vz_ptr(), x, y, z, -1,  0, -1, inv_dx, inv_dy, inv_dz, dxVz, dyVz, dzVz);

	__m256 dtexx = dxVx;
	__m256 dteyy = dyVy;
	__m256 dtezz = dzVz;
	__m256 dteyz2 = _mm256_add_ps(dzVy, dyVz);
	__m256 dtexz2 = _mm256_add_ps(dzVx, dxVz);
	__m256 dtexy2 = _mm256_add_ps(dyVx, dxVy);

	//
	// On-kite CIJs
	//
	__m256 c11, c22, c33, c44, c55, c66, c12, c13, c23;	
	Compute_On_Kite_CIJs(x,y,z,c11,c22,c33,c44,c55,c66,c12,c13,c23);

	//
	// PDE
	// 

	long CurrIdx = CmpIdx(x,y,z);

	__m256 mm_dti = _mm256_set1_ps(dti);

	__m256 deta = Compute_ABC(x,y,z,_nx,_ny,_nz,_nabc_top,_nabc_bot,_nabc_sdx,_nabc_sdy,_vpvert_avtop,_vpvert_avbot,inv_dx,inv_dy,inv_dz);
	__m256 dabc = _mm256_sub_ps(_mm256_set1_ps(1.0f), _mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(0.5f), deta), mm_dti));
	dabc = _mm256_div_ps(dabc, _mm256_add_ps(_mm256_set1_ps(1.0f), _mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(0.5f), deta), mm_dti)));

	__m256 old_txx = *((__m256*)(Get_TXX_ptr() + CurrIdx));
	_mm_prefetch((char*)(Get_TXX_ptr() + CurrIdx + _stride_y), _MM_HINT_T0);
	__m256 txx = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(c11, dtexx), _mm256_mul_ps(c12, dteyy)), _mm256_mul_ps(c13, dtezz));
	txx = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_add_ps(dabc, _mm256_set1_ps(1.0f)), mm_dti), txx), _mm256_mul_ps(_mm256_mul_ps(dabc, dabc), old_txx));

	__m256 old_tyy = *((__m256*)(Get_TYY_ptr() + CurrIdx));
	_mm_prefetch((char*)(Get_TYY_ptr() + CurrIdx + _stride_y), _MM_HINT_T0);
	__m256 tyy = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(c12, dtexx), _mm256_mul_ps(c22, dteyy)), _mm256_mul_ps(c23, dtezz));
	tyy = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_add_ps(dabc, _mm256_set1_ps(1.0f)), mm_dti), tyy), _mm256_mul_ps(_mm256_mul_ps(dabc, dabc), old_tyy));

	__m256 old_tzz = *((__m256*)(Get_TZZ_ptr() + CurrIdx));
	_mm_prefetch((char*)(Get_TZZ_ptr() + CurrIdx + _stride_y), _MM_HINT_T0);
	__m256 tzz = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(c13, dtexx), _mm256_mul_ps(c23, dteyy)), _mm256_mul_ps(c33, dtezz));
	tzz = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_add_ps(dabc, _mm256_set1_ps(1.0f)), mm_dti), tzz), _mm256_mul_ps(_mm256_mul_ps(dabc, dabc), old_tzz));

	__m256 old_txy = *((__m256*)(Get_TXY_ptr() + CurrIdx));
	_mm_prefetch((char*)(Get_TXY_ptr() + CurrIdx + _stride_y), _MM_HINT_T0);
	__m256 txy = _mm256_mul_ps(c66, dtexy2);
	txy = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_add_ps(dabc, _mm256_set1_ps(1.0f)), mm_dti), txy), _mm256_mul_ps(_mm256_mul_ps(dabc, dabc), old_txy));

	if (z == 0)
	{
		__m256 c13_ = _mm256_sub_ps(c33, _mm256_add_ps(c55, c55));
		__m256 dum1 = _mm256_div_ps(_mm256_mul_ps(c13_, c13_), c33);
		__m256 dum2 = _mm256_sub_ps(c13_, dum1);
		dum1 = _mm256_sub_ps(c11, dum1);
		txx = _mm256_add_ps(old_txx, _mm256_mul_ps(mm_dti, _mm256_add_ps(_mm256_mul_ps(dum1, dxVx), _mm256_mul_ps(dum2, dyVy))));
		tyy = _mm256_add_ps(old_tyy, _mm256_mul_ps(mm_dti, _mm256_add_ps(_mm256_mul_ps(dum2, dxVx), _mm256_mul_ps(dum1, dyVy))));
		txy = _mm256_setzero_ps();
		tzz = _mm256_setzero_ps();
	}
	*((__m256*)(Get_TXX_ptr() + CurrIdx)) = txx;
	*((__m256*)(Get_TYY_ptr() + CurrIdx)) = tyy;
	*((__m256*)(Get_TZZ_ptr() + CurrIdx)) = tzz;
	*((__m256*)(Get_TXY_ptr() + CurrIdx)) = txy;

	__m256 old_txz = *((__m256*)(Get_TXZ_ptr() + CurrIdx));
	_mm_prefetch((char*)(Get_TXZ_ptr() + CurrIdx + _stride_y), _MM_HINT_T0);
	__m256 txz = _mm256_mul_ps(c55, dtexz2);
	txz = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_add_ps(dabc, _mm256_set1_ps(1.0f)), mm_dti), txz), _mm256_mul_ps(_mm256_mul_ps(dabc, dabc), old_txz));
	*((__m256*)(Get_TXZ_ptr() + CurrIdx)) = txz;

	__m256 old_tyz = *((__m256*)(Get_TYZ_ptr() + CurrIdx));
	_mm_prefetch((char*)(Get_TYZ_ptr() + CurrIdx + _stride_y), _MM_HINT_T0);
	__m256 tyz = _mm256_mul_ps(c44, dteyz2);
	tyz = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(_mm256_add_ps(dabc, _mm256_set1_ps(1.0f)), mm_dti), tyz), _mm256_mul_ps(_mm256_mul_ps(dabc, dabc), old_tyz));
	*((__m256*)(Get_TYZ_ptr() + CurrIdx)) = tyz;
}

void Vanilla_Visco_Elastic_CPU_Propagator::Propagate_Particle_Velocities_Tile(
	float dti,
	float fq,
	int x0,
	int y0,
	int iz,
	int bx,
	int by,
	float inv_dx,
	float inv_dy,
	float inv_dz
	)
{
	for (int iy = 0;  iy < by;  ++iy)
	{
		for (int ix = 0;  ix < bx;  ix+=8)
		{
			Propagate_Particle_Velocities_8x(dti,fq,ix+x0,iy+y0,iz,inv_dx,inv_dy,inv_dz);
		}
	}
}

void Vanilla_Visco_Elastic_CPU_Propagator::Propagate_Particle_Velocities_8x(
		float dti,
		float fq,
		int x,
		int y,
		int z,
		float inv_dx,
		float inv_dy,
		float inv_dz
		)
{
	// 
	// Compute derivatives
	//
	__m256 dxTxx;
	Compute_DX(Get_TXX_ptr(), x, y, z, -1, inv_dx, dxTxx);

	__m256 dyTyy;
	Compute_DY(Get_TYY_ptr(), x, y, z,  0, inv_dy, dyTyy);

	__m256 dzTzz;
	Compute_DZ(Get_TZZ_ptr(), x, y, z,  0, inv_dz, dzTzz);

	__m256 dxTxy;
	Compute_DX(Get_TXY_ptr(), x, y, z,  0, inv_dx, dxTxy);

	__m256 dyTxy;
	Compute_DY(Get_TXY_ptr(), x, y, z, -1, inv_dy, dyTxy);

	__m256 dxTxz;
	Compute_DX(Get_TXZ_ptr(), x, y, z,  0, inv_dx, dxTxz);

	__m256 dzTxz;
	Compute_DZ(Get_TXZ_ptr(), x, y, z, -1, inv_dz, dzTxz);

	__m256 dyTyz;
	Compute_DY(Get_TYZ_ptr(), x, y, z, -1, inv_dy, dyTyz);

	__m256 dzTyz;
	Compute_DZ(Get_TYZ_ptr(), x, y, z, -1, inv_dz, dzTyz);

	//
	// PDE
	//

	long CurrIdx = CmpIdx(x,y,z);

	__m256 mm_dti = _mm256_set1_ps(dti);
	
	// ..compute itausig and difitau
	__m256 Q = *((__m256*)(Get_Q_ptr() + CurrIdx));
	//
	// prefetch earth model properties for next block-Y loop iteration into  non-temporal cache area.
	// this helps reduce cache usage since the NTA is cleared out after the first usage.
	//
	_mm_prefetch((char*)(Get_Q_ptr() + CurrIdx + _stride_y), _MM_HINT_NTA);
	__m256 wq = _mm256_set1_ps(fq*6.283185307179586476925286766559f);
	__m256 te = _mm256_add_ps(_mm256_set1_ps(1.0f), _mm256_sqrt_ps(_mm256_add_ps(_mm256_set1_ps(1.0f), _mm256_mul_ps(Q,Q))));
	te = _mm256_div_ps(te, _mm256_mul_ps(Q,wq));
	__m256 itausig = _mm256_mul_ps(te, _mm256_mul_ps(wq,wq));
	__m256 difitau = _mm256_sub_ps(_mm256_div_ps(_mm256_set1_ps(1.0f), te), itausig);

	// Update viscoelastic(SLS) vector field:
	__m256 itausig_term = _mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(0.5f), mm_dti), itausig);
	__m256 const1 = _mm256_div_ps(_mm256_set1_ps(1.0f), _mm256_add_ps(_mm256_set1_ps(1.0f), itausig_term));
	__m256 const2 = _mm256_sub_ps(_mm256_set1_ps(1.0f), itausig_term);
	__m256 const3 = _mm256_mul_ps(mm_dti, difitau);

	__m256 old_sx = *((__m256*)(Get_Sx_ptr() + CurrIdx));
	__m256 old_sy = *((__m256*)(Get_Sy_ptr() + CurrIdx));
	__m256 old_sz = *((__m256*)(Get_Sz_ptr() + CurrIdx));
	//
	// prefetch for next iteration of block-Y loop into L2 cache
	//
	_mm_prefetch((char*)(Get_Sx_ptr() + CurrIdx + _stride_y), _MM_HINT_T0);
	_mm_prefetch((char*)(Get_Sy_ptr() + CurrIdx + _stride_y), _MM_HINT_T0);
	_mm_prefetch((char*)(Get_Sz_ptr() + CurrIdx + _stride_y), _MM_HINT_T0);

	__m256 sx = _mm256_mul_ps(const3, _mm256_add_ps(dxTxx, _mm256_add_ps(dyTxy, dzTxz)));
	sx = _mm256_add_ps(sx, _mm256_mul_ps(const2, old_sx));
	sx = _mm256_mul_ps(const1, sx);

	__m256 sy = _mm256_mul_ps(const3, _mm256_add_ps(dxTxy, _mm256_add_ps(dyTyy, dzTyz)));
	sy = _mm256_add_ps(sy, _mm256_mul_ps(const2, old_sy));
	sy = _mm256_mul_ps(const1, sy);

	__m256 sz = _mm256_mul_ps(const3, _mm256_add_ps(dxTxz, _mm256_add_ps(dyTyz, dzTzz)));
	sz = _mm256_add_ps(sz, _mm256_mul_ps(const2, old_sz));
	sz = _mm256_mul_ps(const1, sz);

	*((__m256*)(Get_Sx_ptr() + CurrIdx)) = sx;
	*((__m256*)(Get_Sy_ptr() + CurrIdx)) = sy;
	*((__m256*)(Get_Sz_ptr() + CurrIdx)) = sz;

	// Absorbing boundary decay funct (for Maxwell viscoelastic model):
	__m256 deta = Compute_ABC(x,y,z,_nx,_ny,_nz,_nabc_top,_nabc_bot,_nabc_sdx,_nabc_sdy,_vpvert_avtop,_vpvert_avbot,inv_dx,inv_dy,inv_dz);
	__m256 dabc = _mm256_sub_ps(_mm256_set1_ps(1.0f), _mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(0.5f), deta), mm_dti));
	dabc = _mm256_div_ps(dabc, _mm256_add_ps(_mm256_set1_ps(1.0f), _mm256_mul_ps(_mm256_mul_ps(_mm256_set1_ps(0.5f), deta), mm_dti)));

	// Update viscoelastic particle velocities:
	__m256 old_vx = *((__m256*)(Get_Vx_ptr() + CurrIdx));
	__m256 old_vy = *((__m256*)(Get_Vy_ptr() + CurrIdx));
	__m256 old_vz = *((__m256*)(Get_Vz_ptr() + CurrIdx));
	//
	// prefetch for next iteration of block-Y loop into L2 cache
	//
	_mm_prefetch((char*)(Get_Vx_ptr() + CurrIdx + _stride_y), _MM_HINT_T0);
	_mm_prefetch((char*)(Get_Vy_ptr() + CurrIdx + _stride_y), _MM_HINT_T0);
	_mm_prefetch((char*)(Get_Vz_ptr() + CurrIdx + _stride_y), _MM_HINT_T0);

	__m256 Density = *((__m256*)(Get_Rho_ptr() + CurrIdx));
	//
	// prefetch earth model properties for next block-Y loop iteration into  non-temporal cache area.
	// this helps reduce cache usage since the NTA is cleared out after the first usage.
	//
	_mm_prefetch((char*)(Get_Rho_ptr() + CurrIdx + _stride_y), _MM_HINT_NTA);
	__m256 factor = _mm256_div_ps(mm_dti, Density);

	__m256 vx = _mm256_add_ps(_mm256_mul_ps(dabc,old_vx), _mm256_mul_ps(factor, _mm256_add_ps(dxTxx, _mm256_add_ps(dyTxy, _mm256_add_ps(dzTxz, sx)))));
	__m256 vy = _mm256_add_ps(_mm256_mul_ps(dabc,old_vy), _mm256_mul_ps(factor, _mm256_add_ps(dxTxy, _mm256_add_ps(dyTyy, _mm256_add_ps(dzTyz, sy)))));
	__m256 vz = _mm256_add_ps(_mm256_mul_ps(dabc,old_vz), _mm256_mul_ps(factor, _mm256_add_ps(dxTxz, _mm256_add_ps(dyTyz, _mm256_add_ps(dzTzz, sz)))));

	*((__m256*)(Get_Vx_ptr() + CurrIdx)) = vx;
	*((__m256*)(Get_Vy_ptr() + CurrIdx)) = vy;
	*((__m256*)(Get_Vz_ptr() + CurrIdx)) = vz;
}

void Vanilla_Visco_Elastic_CPU_Propagator::Prepare_For_Propagation()
{
	_vpvert_avtop = 0.0f;
	_vpvert_avbot = 0.0f;
	for (int y = 0;  y < Get_NY();  ++y)
	{
		for (int x = 0;  x < Get_NX();  ++x)
		{
			_vpvert_avtop += Vp(x,y,0);
			_vpvert_avbot += Vp(x,y,Get_NZ()-1);
		}
	}
	_vpvert_avtop /= (float)(Get_NY() * Get_NX());
	_vpvert_avbot /= (float)(Get_NY() * Get_NX());
}

double Vanilla_Visco_Elastic_CPU_Propagator::Timestep(float dti)
{
	struct timespec ts1;
	clock_gettime(CLOCK_REALTIME, &ts1);

	//
	// tiling is done on the X-Y face of the volume.
	//

	// figure out how many tiles we have
	int Nbx = (_nx + _bx - 1) / _bx;
	int Nby = (_ny + _by - 1) / _by;
	int NNN = Nbx * Nby;

	float inv_dx = 1.0f / _dx;
	float inv_dy = 1.0f / _dy;
	float inv_dz = 1.0f / _dz;

	// propagate stresses one timestep
#pragma omp parallel for
	for (int iTile = 0;  iTile < NNN;  ++iTile)
	{
		// calculate tile indexes.
		// these range from 0 to Nbx-1 for ibx and 0 to Nby-1 for iby
		int iby = iTile / Nbx;
		int ibx = (iTile - iby * Nbx);

		// calculate X and Y coordinates for this tile
		int x0 = ibx * _bx;
		int y0 = iby * _by;
		// int z0 = 0;  // by definition
		
		// clip tile
		int bx = x0 + _bx > _nx ? _nx - x0 : _bx;
		int by = y0 + _by > _ny ? _ny - y0 : _by;

		// push the tile along Z from one end of the volume to the other
		for (int iz = 0;  iz < _nz;  ++iz)
		{
			Propagate_Stresses_Orthorhombic_Tile(dti,x0,y0,iz,bx,by,inv_dx,inv_dy,inv_dz);
		}
	}

	// propagate particle velocities one timestep
#pragma omp parallel for
	for (int iTile = 0;  iTile < NNN;  ++iTile)
	{
		// calculate tile indexes.
		// these range from 0 to Nbx-1 for ibx and 0 to Nby-1 for iby
		int iby = iTile / Nbx;
		int ibx = (iTile - iby * Nbx);

		// calculate X and Y coordinates for this tile
		int x0 = ibx * _bx;
		int y0 = iby * _by;
		// int z0 = 0;  // by definition
		
		// clip tile
		int bx = x0 + _bx > _nx ? _nx - x0 : _bx;
		int by = y0 + _by > _ny ? _ny - y0 : _by;

		// push the tile along Z from one end of the volume to the other
		for (int iz = 0;  iz < _nz;  ++iz)
		{
			Propagate_Particle_Velocities_Tile(dti,_fq,x0,y0,iz,bx,by,inv_dx,inv_dy,inv_dz);
		}
	}

	struct timespec ts2;
	clock_gettime(CLOCK_REALTIME, &ts2);

	double elapsed_time = (double)ts2.tv_sec + (double)ts2.tv_nsec * 1e-9 - (double)ts1.tv_sec - (double)ts1.tv_nsec * 1e-9;
	double mpt_per_second = (double)Get_NX() * (double)Get_NY() * (double)Get_NZ() * 1e-6 / elapsed_time;

	return mpt_per_second;
}

