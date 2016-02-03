#ifndef CVX_ESDRD_MI_TMJ_PROPAGATE_PARTICLE_VELOCITIES
#define CVX_ESDRD_MI_TMJ_PROPAGATE_PARTICLE_VELOCITIES

#include "Elastic_Interpolation.hxx"

void 
Host_Propagate_Particle_Velocities_Kernel(
	int timestep,
	cudaStream_t stream,
	int spatial_order,	// either 8 or 16
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
	float C8_0,
	float C8_1,
	float C8_2,
	float C8_3,
	float C16_0,
	float C16_1,
	float C16_2,
	float C16_3,
	float C16_4,
	float C16_5,
	float C16_6,
	float C16_7,
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
	float Vp_min,
	float Vp_range,
	float Vs_min,
	float Vs_range,
	float Density_min,
	float Density_range,
	int one_y_size,
	bool inject_source,
	bool source_ghost_enabled,
	int ghost_sea_surface_z,
	Elastic_Interpolation_t source_interpolation_method,
	bool is_force,
	bool is_velocity,
	float ampl1,
	float ampl2,
	float ampl3,
	float svaw_sample,
	float xsou,
	float ysou,
	float zsou,
	bool is_p_reciprocity,
	float bmod_ref,
	float rho_ref
	);

#endif

