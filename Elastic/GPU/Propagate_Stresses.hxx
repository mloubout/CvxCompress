#ifndef CVX_ESDRD_MI_TMJ_PROPAGATE_STRESSES_HXX
#define CVX_ESDRD_MI_TMJ_PROPAGATE_STRESSES_HXX

#include "Elastic_Interpolation.hxx"

void Host_Propagate_Stresses_Orthorhombic_Kernel(
	int timestep,
	cudaStream_t stream,
	int x0,			// x coordinate of westernmost coordinate in block
	int y0,			// y coordinate of southernmost coordinate in block
	int y1,
	int m1_y0,
	int m1_y1,
	int vol_nx,		// dimensions of global volume
	int vol_ny,
	int vol_nz,
	float dti,
	unsigned int* em,	// earth model, 4 interleaved integers. y(0)
	float* cmp,		// txx, tyy, tzz, txy, txz and tyz, middle, t(1), y(0)
	float* m1L,		// Vx, Vy, Vz, Sx, Sy, Sz in that order. left halo, t(0), y(0)
        float* m1C,		// ..middle, t(0), y(0)
        float* m1R,		// ..right halo, t(0), y(0)
        float* m2C,		// txx, tyy, tzz, txy, txz and tyz. middle, t(0), y(0)
	float C0,
	float C1,
	float C2,
	float C3,
	float inv_DX,		// 1 / DX
	float inv_DY,		// 1 / DY
	float inv_DZ,		// 1 / DZ
	float vpvert_avtop,
	float vpvert_avbot,
	int nabc_sdx,
	int nabc_sdy,
	int nabc_top,
	int nabc_bot,
	float Vp_min,
	float Vp_range,
	float Vs_min,
	float Vs_range,
	float Density_min,
	float Density_range,
	float Dip_min, 
	float Dip_range, 
	float Azimuth_min, 
	float Azimuth_range, 
	float Rake_min, 
	float Rake_range, 
	float Delta1_min,
	float Delta1_range,
	float Delta2_min,
	float Delta2_range,
	float Delta3_min,
	float Delta3_range,
	float Epsilon1_min,
	float Epsilon1_range,
	float Epsilon2_min,
	float Epsilon2_range,
	float Gamma1_min,
	float Gamma1_range,
	float Gamma2_min,
	float Gamma2_range,
	int one_y_size,
	bool inject_source,
	Elastic_Interpolation_t source_interpolation_method,
	bool is_pressure,
	float ampl1,
	float svaw_sample,
	float xsou,
	float ysou,
	float zsou
	);

#endif

