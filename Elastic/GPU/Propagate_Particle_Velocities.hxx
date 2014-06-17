#ifndef CVX_ESDRD_MI_TMJ_PROPAGATE_PARTICLE_VELOCITIES
#define CVX_ESDRD_MI_TMJ_PROPAGATE_PARTICLE_VELOCITIES

#include "Elastic_Interpolation.hxx"

void 
Host_Propagate_Particle_Velocities_Kernel(
	int timestep,
	cudaStream_t stream,
	int num_z,
	int x0,
	int y0,
	int y1,
	int m1_y0,
	int m1_y1,
	int vol_nx,
	int vol_ny,
	int vol_nz,
	float dti,
	void* em,		// earth model
	void* cmp,		// newly computed values should be stored here
	void* m1L,		// strain rates, left halo
	void* m1C,		// strain rates, middle
	void* m1R,		// strain rates, right halo
	void* m2C,		// particle velocities from previous timestep, middle
	float C0,
        float C1,
        float C2,
        float C3,
        float inv_DX,           // 1 / DX
        float inv_DY,           // 1 / DY
        float inv_DZ,           // 1 / DZ
	float vpvert_avtop,
	float vpvert_avbot,
	int nabc_sdx,
	int nabc_sdy,
	int nabc_top,
	int nabc_bot,
	float Q_min,
	float Q_range,
	float fq,
	float Density_min,
	float Density_range,
	int one_y_size,
	bool inject_source,
	Elastic_Interpolation_t source_interpolation_method,
	bool is_force,
	bool is_velocity,
	float ampl1,
	float ampl2,
	float ampl3,
	float svaw_sample,
	float xsou,
	float ysou,
	float zsou
	);

#endif

