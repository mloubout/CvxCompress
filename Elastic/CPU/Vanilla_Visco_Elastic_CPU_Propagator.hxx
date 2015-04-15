#ifndef VANILLA_VISCO_ELASTIC_CPU_PROPAGATOR_HXX
#define VANILLA_VISCO_ELASTIC_CPU_PROPAGATOR_HXX

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <emmintrin.h>
#include <tmmintrin.h>
#include <immintrin.h>

//
// Vanilla implementation of visco-elastic isotropic finite difference propagator.
// 8th order in space, 2nd order in time.
//
// Original algorithm provided by Kurt Nihei.
// Questions can be directed to Thor Johnsen, thor.johnsen@chevron.com
//

class Vanilla_Visco_Elastic_CPU_Propagator
{
public:
	Vanilla_Visco_Elastic_CPU_Propagator(
		int nx,
		int ny,
		int nz,
		float dx,
		float dy,
		float dz,
		float fq
		);
	~Vanilla_Visco_Elastic_CPU_Propagator();

	static void Run_Unit_Tests();

	inline int Get_NX() {return _nx;}
	inline int Get_NY() {return _ny;}
	inline int Get_NZ() {return _nz;}

	static float Compute_Max_Internal_Timestep(float max_V, float Courant_Multiplier);

	//
	// Compute internal constants needed for propagation.
	//
	void Prepare_For_Propagation();

	//
	// Advance all wavefields "Internal_Timestep" seconds.
	// Returns throughput in million points per second.
	//
	double Timestep(float Internal_Timestep);

	//
	// Run properties
	//
	
	inline void Set_Tile_Size(int bx, int by) {_bx = bx;  _by = by;}

	//
	// Modeling properties
	//
	
	inline float Get_FQ() {return _fq;}
	inline void Set_FQ(float fq) {_fq = fq;}

	//
	// "Get" methods for the 12 wavefields.
	// Should be reasonably fast, but rely heavily on the compiler being able to eliminate common sub expressions.
	// x must be within [0,Get_NX()-1]
	// y must be within [0,Get_NY()-1]
	// z must be within [0,Get_NZ()-1]
	// 

	inline float Vx(int x, int y, int z) {return Get_Vx_ptr()[CmpIdx(x,y,z)];}
	inline float Vy(int x, int y, int z) {return Get_Vy_ptr()[CmpIdx(x,y,z)];}
	inline float Vz(int x, int y, int z) {return Get_Vz_ptr()[CmpIdx(x,y,z)];}

	inline float Sx(int x, int y, int z) {return Get_Sx_ptr()[CmpIdx(x,y,z)];}
	inline float Sy(int x, int y, int z) {return Get_Sy_ptr()[CmpIdx(x,y,z)];}
	inline float Sz(int x, int y, int z) {return Get_Sz_ptr()[CmpIdx(x,y,z)];}

	inline float Txx(int x, int y, int z) {return Get_TXX_ptr()[CmpIdx(x,y,z)];}
	inline float Tyy(int x, int y, int z) {return Get_TYY_ptr()[CmpIdx(x,y,z)];}
	inline float Tzz(int x, int y, int z) {return Get_TZZ_ptr()[CmpIdx(x,y,z)];}
	inline float Txy(int x, int y, int z) {return Get_TXY_ptr()[CmpIdx(x,y,z)];}
	inline float Txz(int x, int y, int z) {return Get_TXZ_ptr()[CmpIdx(x,y,z)];}
	inline float Tyz(int x, int y, int z) {return Get_TYZ_ptr()[CmpIdx(x,y,z)];}

	//
	// "Set" methods for the 12 wavefields.
	// Should be reasonably fast, but rely heavily on the compiler being able to eliminate common sub expressions.
	// x must be within [0,Get_NX()-1]
	// y must be within [0,Get_NY()-1]
	// z must be within [0,Get_NZ()-1]
	//

	inline void Set_Vx(int x, int y, int z, float new_Vx) {Get_Vx_ptr()[CmpIdx(x,y,z)] = new_Vx;}
	inline void Set_Vy(int x, int y, int z, float new_Vy) {Get_Vy_ptr()[CmpIdx(x,y,z)] = new_Vy;}
	inline void Set_Vz(int x, int y, int z, float new_Vz) {Get_Vz_ptr()[CmpIdx(x,y,z)] = new_Vz;}

	inline void Set_Sx(int x, int y, int z, float new_Sx) {Get_Sx_ptr()[CmpIdx(x,y,z)] = new_Sx;}
	inline void Set_Sy(int x, int y, int z, float new_Sy) {Get_Sy_ptr()[CmpIdx(x,y,z)] = new_Sy;}
	inline void Set_Sz(int x, int y, int z, float new_Sz) {Get_Sz_ptr()[CmpIdx(x,y,z)] = new_Sz;}

	inline void Set_Txx(int x, int y, int z, float new_Txx) {Get_TXX_ptr()[CmpIdx(x,y,z)] = new_Txx;}
	inline void Set_Tyy(int x, int y, int z, float new_Tyy) {Get_TYY_ptr()[CmpIdx(x,y,z)] = new_Tyy;}
	inline void Set_Tzz(int x, int y, int z, float new_Tzz) {Get_TZZ_ptr()[CmpIdx(x,y,z)] = new_Tzz;}
	inline void Set_Txy(int x, int y, int z, float new_Txy) {Get_TXY_ptr()[CmpIdx(x,y,z)] = new_Txy;}
	inline void Set_Txz(int x, int y, int z, float new_Txz) {Get_TXZ_ptr()[CmpIdx(x,y,z)] = new_Txz;}
	inline void Set_Tyz(int x, int y, int z, float new_Tyz) {Get_TYZ_ptr()[CmpIdx(x,y,z)] = new_Tyz;}

	//
	// "Get" methods for the earth model properties.
	// Should be reasonably fast, but rely heavily on the compiler being able to eliminate common sub expressions.
	// x must be within [0,Get_NX()-1]
	// y must be within [0,Get_NY()-1]
	// z must be within [0,Get_NZ()-1]
	//
	
	inline float Vp(int x, int y, int z) {return Get_Vp_ptr()[CmpIdx(x,y,z)];}
	inline float Vs(int x, int y, int z) {return Get_Vs_ptr()[CmpIdx(x,y,z)];}
	inline float Rho(int x, int y, int z) {return Get_Rho_ptr()[CmpIdx(x,y,z)];}
	inline float Q(int x, int y, int z) {return Get_Q_ptr()[CmpIdx(x,y,z)];}

	//
	// "Set" methods for earth model properties.
	// Should be reasonably fast, but rely heavily on the compiler being able to eliminate common sub expressions.
	// x must be within [0,Get_NX()-1]
	// y must be within [0,Get_NY()-1]
	// z must be within [0,Get_NZ()-1]
	//
	
	inline void Set_Vp(int x, int y, int z, float new_Vp) {Get_Vp_ptr()[CmpIdx(x,y,z)] = new_Vp;}
	inline void Set_Vs(int x, int y, int z, float new_Vs) {Get_Vs_ptr()[CmpIdx(x,y,z)] = new_Vs;}
	inline void Set_Rho(int x, int y, int z, float new_Rho) {Get_Rho_ptr()[CmpIdx(x,y,z)] = new_Rho;}
	inline void Set_Q(int x, int y, int z, float new_Q) {Get_Q_ptr()[CmpIdx(x,y,z)] = new_Q;}

private:
	// Get earth model properties as flat arrays.
	inline float* Get_Vp_ptr() {return _Vp;}
	inline float* Get_Vs_ptr() {return _Vs;}
	inline float* Get_Rho_ptr() {return _Rho;}
	inline float* Get_Q_ptr() {return _Qp;}

	// Get wavefields as flat arrays.
	inline float* Get_TXX_ptr() {return _txx;}
	inline float* Get_TYY_ptr() {return _tyy;}
	inline float* Get_TZZ_ptr() {return _tzz;}
	inline float* Get_TXY_ptr() {return _txy;}
	inline float* Get_TXZ_ptr() {return _txz;}
	inline float* Get_TYZ_ptr() {return _tyz;}

	inline float* Get_Vx_ptr() {return _vx;}
	inline float* Get_Vy_ptr() {return _vy;}
	inline float* Get_Vz_ptr() {return _vz;}
	
	inline float* Get_Sx_ptr() {return _sx;}
	inline float* Get_Sy_ptr() {return _sy;}
	inline float* Get_Sz_ptr() {return _sz;}

	// Use this method to compute flat array index for x,y,z coordinate
	inline long CmpIdx(int x, int y, int z) {return ((long)z * (long)_ny + (long)y) * (long)_nx + (long)x;}

private:
	// compute 1st order dx, dy and dz derivative of wavefield for 8 consecutive points (along X, fast dimension).
	inline void Compute_DX(
			float* wf,
			int x,
			int y,
			int z,
			int xoff,
			float inv_dx,
			__m256& dx
			);
	inline void Compute_DY(
			float* wf,
			int x,
			int y,
			int z,
			int yoff,
			float inv_dy,
			__m256& dy
			);
	inline void Compute_DZ(
			float* wf,
			int x,
			int y,
			int z,
			int zoff,
			float inv_dz,
			__m256& dz
			);
	inline void Compute_DX_DY_DZ(
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
			);
	inline void Compute_On_Kite_CIJs(
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
			);
	inline float _Compute_ABC(
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
			);
	__m256 Compute_ABC(
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
                        );
	inline void Propagate_Stresses_Orthorhombic_8x(
			float dti,
			int x,
			int y,
			int z,
			float inv_dx,
			float inv_dy,
			float inv_dz
			);
	inline void Propagate_Stresses_Orthorhombic_Tile(
			float dti,
			int x0,
			int y0,
			int iz,
			int bx,
			int by,
			float inv_dx,
			float inv_dy,
			float inv_dz
			);
	inline void Propagate_Particle_Velocities_8x(
			float dti,
			float fq,
			int x,
			int y,
			int z,
			float inv_dx,
			float inv_dy,
			float inv_dz
			);
	inline void Propagate_Particle_Velocities_Tile(
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
			);
	float WLUP(float* wf, int x, int y, int z);
	void Unit_Test_Comp_DX_DY_DZ(
			int nx,
			int ny,
			int nz,
			float dx,
			float dy,
			float dz,
			int xoff,
			int yoff,
			int zoff
			);
	void omp_memclear(void* dst, size_t len);
	void omp_memcpy(void* dst, void* src, size_t len);

private:
	int _nx;
	int _ny;
	int _nz;

	int _nabc_top;
	int _nabc_bot;
	int _nabc_sdx;
	int _nabc_sdy;

	float _dx;
	float _dy;
	float _dz;

	float _fq;

	float* _txx;
	float* _tyy;
	float* _tzz;
	float* _txy;
	float* _txz;
	float* _tyz;

	float* _vx;
	float* _vy;
	float* _vz;
	float* _sx;
	float* _sy;
	float* _sz;

	float* _Vp;
	float* _Vs;
	float* _Rho;
	float* _Qp;

	float _vpvert_avtop;
	float _vpvert_avbot;

        // Optimized (for lower dispersion) coefficients
        static const float _C0 = 1.1850912100109303f;
        static const float _C1 = -0.0835270299926924f;
        static const float _C2 = 0.016837760894350576f;
        static const float _C3 = -0.0027386181103430177f;

	float _c[8];

	// tile (or block) size
	int _bx;
	int _by;

	// strides for flattened arrays.
	// _stride_x is equal to 1 by definition.
	long _stride_y;
	long _stride_z;

	// strides for __m256 vector type.
	// __m256 avx vector type is simply 8 consecutive floats,
	// so moving one __m256 value ahead is the same as moving 8 floats.
	long _mm256_stride_y;
	long _mm256_stride_z;
};

#endif

