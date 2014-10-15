#include <cuda_runtime_api.h>
#include "Elastic_Interpolation.hxx"

//
// CUDA kernel that propagates particle velocity wavefield and memory variable.
//
// 
//

__device__ 
int cuCompTXXIdx(int offset)
{
	int abs_offset = offset + (threadIdx.x & 3) + 4;
	int quotient = abs_offset / 4;
	int remainder = abs_offset & 3;
	return threadIdx.y*96 + quotient*32 + (threadIdx.x&28) + remainder;
}

__device__ 
int cuCompTXXIdx_2(int offset)
{
	int abs_offset = offset + (threadIdx.x & 3) + 4;
	int quotient = abs_offset / 4;
	int remainder = abs_offset & 3;
	return (threadIdx.y+4)*96 + quotient*32 + (threadIdx.x&28) + remainder;
}

__device__ 
int cuCompTYYIdx(int offset)
{
	return (offset + 4 + threadIdx.y) * 32 + threadIdx.x;
}

__device__ 
int cuCompTYYIdx_2(int offset)
{
	return (offset + 8 + threadIdx.y) * 32 + threadIdx.x;
}

__device__ 
float cuTransposeXZY2XYZ(float* buf, float v)
{
	__syncthreads();  // wait for previous step to finish using buf
	buf[threadIdx.x+((threadIdx.x&28)*8)+threadIdx.y*4] = v;
	__syncthreads();  // wait for all threads to finish writing to buf
	float retval = buf[threadIdx.x+threadIdx.y*36];
	__syncthreads();  // wait for all threads to finish reading from buf
	return retval;
}

__device__ 
float cuBessi0(float X)
{
	// Modified Bessel Function of zero order.
	// From Numerical Recipes, Press et al. (1986), pp. 177

	const float P1 = 1.0f;
	const float P2 = 3.5156229f;
	const float P3 = 3.0899424f;
	const float P4 = 1.2067492f;
	const float P5 = 0.2659732f;
	const float P6 = 0.360768e-1f;
	const float P7 = 0.45813e-2f;

	const float Q1 = 0.39894228f;
	const float Q2 = 0.1328592e-1f;
	const float Q3 = 0.225319e-2f;
	const float Q4 = -0.157565e-2f;
	const float Q5 = 0.916281e-2f;
	const float Q6 = -0.2057706e-1f;
	const float Q7 = 0.2635537e-1f;
	const float Q8 = -0.1647633e-1f;
	const float Q9 = 0.392377e-2f;

	if (fabsf(X) < 3.75f)
	{
		float Y = (X*X) / (3.75f*3.75f);
		return P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))));
	}
	else
	{
		float AX = fabsf(X);
		float Y = 3.75f / AX;
		return (expf(AX)/sqrtf(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))));
	}
}

__device__ 
float cuGen_Single_Sinc_Weight(
	int tx,
	float dx_frac
	)
{
	const float b = 4.14f; // optimal Kaiser window param for kmax = 2pi/3
        const float r = 4.0f;  // half-width of sinc interpolator
        const float pi = 3.1415926535897932384626433832795f;

	int ix = tx + 1;

	// cells at which to sample sinc func [normalized]
	float x_cell = (float)ix - r - dx_frac;

	// compute Kaiser window:
	float b_x = (fabsf(x_cell) <= r) ? b*sqrtf(1.0f - ((x_cell*x_cell)/(r*r))) : 0.0f;
	float win_x = cuBessi0(b_x) / cuBessi0(b);

	// compute sinc interpolation function:
	float fsinc_x = (x_cell == 0.0f) ? 1.0f : win_x * sin(x_cell*pi)/(x_cell*pi);

	return fsinc_x;
}

__device__ 
float cuGen_Sinc_Weight(
	int tx,
	int ty,
	int tz,
	float dx_frac,
	float dy_frac,
	float dz_frac
	)
{
	return cuGen_Single_Sinc_Weight(tx,dx_frac) * cuGen_Single_Sinc_Weight(ty,dy_frac) * cuGen_Single_Sinc_Weight(tz,dz_frac);

	/*
	const float b = 4.14f; // optimal Kaiser window param for kmax = 2pi/3
	const float r = 4.0f;  // half-width of sinc interpolator
	const float pi = 3.1415926535897932384626433832795f;

	int ix = tx + 1;
	int iy = ty + 1;
	int iz = tz + 1;

	// cells at which to sample sinc func [normalized]
	float z_cell = (float)iz - r - dz_frac;

	// compute Kaiser window:
	float b_z = (fabsf(z_cell) <= r) ? b * sqrtf(1.0f - ((z_cell*z_cell)/(r*r))) : 0.0f;
	float win_z = cuBessi0(b_z) / cuBessi0(b);

	// compute sinc interpolation function:
	float fsinc_z = (z_cell == 0.0f) ? 1.0f : win_z * sin(z_cell*pi)/(z_cell*pi);

	// cells at which to sample sinc func [normalized]
	float y_cell = (float)iy - r - dy_frac;

	// compute Kaiser window:
	float b_y = (fabsf(y_cell) <= r) ? b * sqrtf(1.0f - ((y_cell*y_cell)/(r*r))) : 0.0f;
	float win_y = cuBessi0(b_y) / cuBessi0(b);

	// compute sinc interpolation function:
	float fsinc_y = (y_cell == 0.0f) ? 1.0f : win_y * sin(y_cell*pi)/(y_cell*pi);

	// cells at which to sample sinc func [normalized]
	float x_cell = (float)ix - r - dx_frac;

	// compute Kaiser window:
	float b_x = (fabsf(x_cell) <= r) ? b*sqrtf(1.0f - ((x_cell*x_cell)/(r*r))) : 0.0f;
	float win_x = cuBessi0(b_x) / cuBessi0(b);

	// compute sinc interpolation function:
	float fsinc_x = (x_cell == 0.0f) ? 1.0f : win_x * sin(x_cell*pi)/(x_cell*pi);

	return fsinc_x * fsinc_y * fsinc_z;
	*/
}

__device__ 
void _cuApply_Source_Term_To_VxVyVz(
	int thr_z,
	unsigned int* em,
	float Q_min,
	float Q_range,
	float Density_min,
	float Density_range,
	float* cmp,
	int x0,
        int y0,
        int z0,
        int nx,
        int ny,
        int nz,
	float dti,
        bool is_force,
	float ampl1,
	float ampl2,
	float ampl3,
        float xs,
        float ys,
        float zs,
        float val,
	int icell,
	int jcell,
	int kcell
	)
{
	int my_x = icell + threadIdx.x - 3 - x0;
	int my_y = jcell + threadIdx.y - 3 - y0;
	int my_z = kcell + thr_z       - 3;
	
	if (
			(my_x >= 0 && my_x < nx) && 
			(my_y >= 0 && my_y < ny) && 
			(my_z > -4 && my_z < nz) 
	   )
	{
		// ..fractional distance from grid pt to sou:
		float dx_frac = xs - (float)(icell - 1);
		float dy_frac = ys - (float)(jcell - 1);
                float dz_frac = zs - (float)(kcell - 1);

		// (fx/vx sou needs to be shifted +0.5icell to colloc w/ pr)
		float vx_dx_frac = xs + 0.5 - (float)(icell - 1);
		// (fy/vy sou needs to be shifted -0.5dy to colloc w/ pr)
		float vy_dy_frac = ys - 0.5 - (float)(jcell - 1);
		// (fz/vz sou need to be shifted -0.5kcell to colloc w/ pr)
		float vz_dz_frac = zs - 0.5 - (float)(kcell - 1);

		float vx_fsinc = cuGen_Sinc_Weight(threadIdx.x,threadIdx.y,thr_z,vx_dx_frac,dy_frac,dz_frac);
		float vy_fsinc = cuGen_Sinc_Weight(threadIdx.x,threadIdx.y,thr_z,dx_frac,vy_dy_frac,dz_frac);
		float vz_fsinc = cuGen_Sinc_Weight(threadIdx.x,threadIdx.y,thr_z,dx_frac,dy_frac,vz_dz_frac);

		if (vx_fsinc != 0.0 || vy_fsinc != 0.0 || vz_fsinc != 0.0)
		{
			// mirror source if necessary
			my_z = my_z < 0 ? -my_z : my_z;

			int one_wf_size_f = nx * nz;
			int one_y_size_f = one_wf_size_f * 6;
			int idx = my_x + my_y * one_y_size_f + my_z * 4;

			if (is_force)
			{
				int em_one_word_size_f = one_wf_size_f;
				int em_one_y_size_f = em_one_word_size_f * 4;
				int em_word3 = em[my_x+my_y*em_one_y_size_f+my_z*4+3*em_one_word_size_f];
				float Q, Density;
				cuUnpack_Q_Density(em_word3,Q_min,Q_range,&Q,Density_min,Density_range,&Density);

				cmp[idx                ] = cmp[idx                ] + vx_fsinc * dti * (val * ampl1) / Density;
				cmp[idx+  one_wf_size_f] = cmp[idx+  one_wf_size_f] + vy_fsinc * dti * (val * ampl2) / Density;
				cmp[idx+2*one_wf_size_f] = cmp[idx+2*one_wf_size_f] + vz_fsinc * dti * (val * ampl3) / Density;

				//printf("Adding source term (%f * [%f-%f-%f] * %f) to Vx,Vy,Vz at %d,%d,%d\n",dti,vx_fsinc,vy_fsinc,vz_fsinc,val/Density,my_x+x0,my_y+y0,my_z);
			}
			else
			{
				cmp[idx                ] = cmp[idx                ] + vx_fsinc * dti * val * ampl1;
				cmp[idx+  one_wf_size_f] = cmp[idx+  one_wf_size_f] + vy_fsinc * dti * val * ampl2;
				cmp[idx+2*one_wf_size_f] = cmp[idx+2*one_wf_size_f] + vz_fsinc * dti * val * ampl3;

				//printf("Adding source term (%f * [%f-%f-%f] * %f) to Vx,Vy,Vz at %d,%d,%d\n",dti,vx_fsinc,vy_fsinc,vz_fsinc,val,my_x+x0,my_y+y0,my_z);
			}
		}
	}
}

__global__ 
void cuApply_Source_Term_To_VxVyVz(
	void* em,
	float Q_min,
	float Q_range,
	float Density_min,
	float Density_range,
	void* cmp,
	int x0,
        int y0,
        int z0,
        int nx,
        int ny,
        int nz,
	float dti,
        bool is_force,
	float ampl1,
	float ampl2,
	float ampl3,
        float xs,
        float ys,
        float zs,
        float val
	)
{
	// fx/vx contribution:

	// nearest grid point:
	int icell = (int)lrintf(xs) + 1; // to left of extrap pt
	int jcell = (int)lrintf(ys) + 1;
	int kcell = (int)lrintf(zs) + 1; // above interp pt:

	for (int thr_z = 0;  thr_z < 8;  ++thr_z)
	{
		_cuApply_Source_Term_To_VxVyVz(thr_z,(unsigned int*)em,Q_min,Q_range,Density_min,Density_range,(float*)cmp,x0,y0,z0,nx,ny,nz,dti,is_force,ampl1,ampl2,ampl3,xs,ys,zs,val,icell,jcell,kcell);
	}
}

__device__ 
void cuAdd_Source_Term_To_Single_VxVyVz(
	int* em,
        float Q_min,
        float Q_range,
        float Density_min,
        float Density_range,
        float* cmp,
        int nx,
        int ny,
        int nz,
        bool is_force,
        float delta_Vx,
        float delta_Vy,
        float delta_Vz,
	int ix,
	int iy,
	int iz,
	int Vx_ix,
	int Vy_iy,
	int Vz_iz,
	int one_wf_size_f,
	int one_y_size_f,
	int em_one_word_size_f,
	int em_one_y_size_f
        )
{
	if (delta_Vz != 0.0f && Vx_ix >= 0 && Vx_ix < nx && iy >= 0 && iy < ny && iz >= 0 && iz < nz)
	{
		int idx_Vx = Vx_ix +    iy * one_y_size_f +    iz * 4;

		if (is_force)
		{
			int em_word3 = em[Vx_ix+iy*em_one_y_size_f+iz*4+3*em_one_word_size_f];
			float Q, Density;
			cuUnpack_Q_Density(em_word3,Q_min,Q_range,&Q,Density_min,Density_range,&Density);
			cmp[idx_Vx] = cmp[idx_Vx] + delta_Vx / Density;
		}
		else
		{
			cmp[idx_Vx] = cmp[idx_Vx] + delta_Vx;
		}
	}
	if (delta_Vy != 0.0f && ix >= 0 && ix < nx && Vy_iy >= 0 && Vy_iy < ny && iz >= 0 && iz < nz)
	{
		int idx_Vy =    ix + Vy_iy * one_y_size_f +    iz * 4 +     one_wf_size_f;

		if (is_force)
		{
			int em_word3 = em[ix+Vy_iy*em_one_y_size_f+iz*4+3*em_one_word_size_f];
			float Q, Density;
			cuUnpack_Q_Density(em_word3,Q_min,Q_range,&Q,Density_min,Density_range,&Density);
			cmp[idx_Vy] = cmp[idx_Vy] + delta_Vy / Density;
		}
		else
		{
			cmp[idx_Vy] = cmp[idx_Vy] + delta_Vy;
		}
	}
	if (delta_Vz != 0.0f && ix >= 0 && ix < nx && iy >= 0 && iy < ny && Vz_iz >= 0 && Vz_iz < nz)
	{
		int idx_Vz =    ix +    iy * one_y_size_f + Vz_iz * 4 + 2 * one_wf_size_f;

		if (is_force)
		{
			int em_word3 = em[ix+iy*em_one_y_size_f+Vz_iz*4+3*em_one_word_size_f];
			float Q, Density;
			cuUnpack_Q_Density(em_word3,Q_min,Q_range,&Q,Density_min,Density_range,&Density);
			cmp[idx_Vz] = cmp[idx_Vz] + delta_Vz / Density;
		}
		else
		{
			cmp[idx_Vz] = cmp[idx_Vz] + delta_Vz;
			//printf("delta_Vz = %e, ix=%d, iy=%d, Vz_iz=%d\n",delta_Vz,ix,iy,Vz_iz);
		}
	}
}

__global__ 
void cuApply_Point_Source_To_VxVyVz(
	int* em,
	float Q_min,
        float Q_range,
        float Density_min,
        float Density_range,
        float* cmp,
        int x0,
        int y0,
        int nx,
        int ny,
        int nz,
        float dti,
        bool is_force,
        float ampl1,
        float ampl2,
        float ampl3,
        float xs,
        float ys,
        float zs,
        float val
        )
{
	int ix = (int)lrintf(xs) - x0;
	int iy = (int)lrintf(ys) - y0;
	int iz = (int)lrintf(zs);

	int one_wf_size_f = nx * nz;
	int one_y_size_f = one_wf_size_f * 6;

	int em_one_word_size_f = one_wf_size_f;
	int em_one_y_size_f = em_one_word_size_f * 4;

	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Q_min,Q_range,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,dti*val*ampl1,dti*val*ampl2,dti*val*ampl3,ix,iy,iz,ix,iy,iz,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f);
}

__global__ 
void cuApply_Trilinear_Source_To_VxVyVz(
	int* em,
	float Q_min,
        float Q_range,
        float Density_min,
        float Density_range,
        float* cmp,
        int x0,
        int y0,
        int nx,
        int ny,
        int nz,
        float dti,
        bool is_force,
        float ampl1,
        float ampl2,
        float ampl3,
        float xs,
        float ys,
        float zs,
        float val
        )
{
	int ix = (int)truncf(xs) - x0;
	int iy = (int)truncf(ys) - y0;
	int iz = (int)truncf(zs);

	float xd = 1.0f - (xs - (float)(ix+x0));
	float yd = 1.0f - (ys - (float)(iy+y0));
	float zd = 1.0f - (zs - (float)iz);

	int Vx_ix = (int)truncf(xs+0.5f) - x0;
	int Vy_iy = (int)truncf(ys-0.5f) - y0;
	int Vz_iz = (int)truncf(zs-0.5f);

	float Vx_xd = 1.0f - (xs - (float)(Vx_ix+x0) + 0.5f);
	float Vy_yd = 1.0f - (ys - (float)(Vy_iy+y0) - 0.5f);
	float Vz_zd = 1.0f - (zs - (float)Vz_iz - 0.5f);

	int one_wf_size_f = nx * nz;
	int one_y_size_f = one_wf_size_f * 6;

	int em_one_word_size_f = one_wf_size_f;
	int em_one_y_size_f = em_one_word_size_f * 4;

	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Q_min,Q_range,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*Vx_xd*   yd*   zd,
		dti*val*ampl2*   xd*Vy_yd*   zd,
		dti*val*ampl3*   xd*   yd*Vz_zd,
		   ix,   iy,   iz,
		Vx_ix,Vy_iy,Vz_iz,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f);
	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Q_min,Q_range,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*(1.0f-Vx_xd)*   yd*   zd,
		dti*val*ampl2*(1.0f-   xd)*Vy_yd*   zd,
		dti*val*ampl3*(1.0f-   xd)*   yd*Vz_zd,
		   ix+1,   iy,   iz,
		Vx_ix+1,Vy_iy,Vz_iz,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f);
	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Q_min,Q_range,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*(1.0f-Vx_xd)*(1.0f-   yd)*   zd,
		dti*val*ampl2*(1.0f-   xd)*(1.0f-Vy_yd)*   zd,
		dti*val*ampl3*(1.0f-   xd)*(1.0f-   yd)*Vz_zd,
		   ix+1,   iy+1,   iz,
		Vx_ix+1,Vy_iy+1,Vz_iz,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f);
	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Q_min,Q_range,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*Vx_xd*(1.0f-   yd)*   zd,
		dti*val*ampl2*   xd*(1.0f-Vy_yd)*   zd,
		dti*val*ampl3*   xd*(1.0f-   yd)*Vz_zd,
		   ix,   iy+1,   iz,
		Vx_ix,Vy_iy+1,Vz_iz,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f);
	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Q_min,Q_range,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*Vx_xd*   yd*(1.0f-   zd),
		dti*val*ampl2*   xd*Vy_yd*(1.0f-   zd),
		dti*val*ampl3*   xd*   yd*(1.0f-Vz_zd),
		   ix,   iy,   iz+1,
		Vx_ix,Vy_iy,Vz_iz+1,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f);
	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Q_min,Q_range,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*(1.0f-Vx_xd)*   yd*(1.0f-   zd),
		dti*val*ampl2*(1.0f-   xd)*Vy_yd*(1.0f-   zd),
		dti*val*ampl3*(1.0f-   xd)*   yd*(1.0f-Vz_zd),
		   ix+1,   iy,   iz+1,
		Vx_ix+1,Vy_iy,Vz_iz+1,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f);
	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Q_min,Q_range,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*(1.0f-Vx_xd)*(1.0f-   yd)*(1.0f-   zd),
		dti*val*ampl2*(1.0f-   xd)*(1.0f-Vy_yd)*(1.0f-   zd),
		dti*val*ampl3*(1.0f-   xd)*(1.0f-   yd)*(1.0f-Vz_zd),
		   ix+1,   iy+1,   iz+1,
		Vx_ix+1,Vy_iy+1,Vz_iz+1,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f);
	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Q_min,Q_range,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*Vx_xd*(1.0f-   yd)*(1.0f-   zd),
		dti*val*ampl2*   xd*(1.0f-Vy_yd)*(1.0f-   zd),
		dti*val*ampl3*   xd*(1.0f-   yd)*(1.0f-Vz_zd),
		   ix,   iy+1,   iz+1,
		Vx_ix,Vy_iy+1,Vz_iz+1,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f);
}

__global__ 
#if __CUDA_ARCH__ >= 300
__launch_bounds__(1280)
#elif __CUDA_ARCH__ >= 200
__launch_bounds__(768)
#endif
void 
cuPropagate_Particle_Velocities_Kernel(
	int timestep,
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
	float* cmp,		// Vx, Vy, Vz, Sx, Sy and Sz, middle, t(1), y(0)
	float* m1L,		// txx, tyy, tzz, txy, txz and tyz in that order. left halo, t(0), y(0)
        float* m1C,		// ..middle, t(0), y(0)
        float* m1R,		// ..right halo, t(0), y(0)
        float* m2C,		// Vx, Vy, Vz, Sx, Sy and Sz in that order. middle, t(-1), y(0)
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
	float Q_min,
	float Q_range,
	float fq,
	float Density_min,
	float Density_range,
	int one_wf_size_f,
	int one_y_size_f,
	int em_one_word_size_f,
	int em_one_y_size_f,
	int tyy_off,
	int tzz_off,
	int txy_off,
	int txz_off,
	int tyz_off
	)
{
	__shared__ float buf[768];	// NON-persistent buffer

	__shared__ float tzzbuf[384];	// persistent buffers
	__shared__ float txzbuf[384];   // some values are transferred from one iZ to the next
	__shared__ float tyzbuf[384];

	int z_per_block = (((vol_nz/8) + gridDim.z - 1) / gridDim.z) * 8;
	int z0 = z_per_block * blockIdx.z;
	int z1 = z0 + z_per_block - 1;
	if (z1 >= vol_nz) z1 = vol_nz - 1;
	int nz = z1 - z0 + 1;
	if (nz <= 0) return;

	int offset = (threadIdx.y + blockIdx.y * 8) * one_y_size_f + threadIdx.x + z0 * 4;

	// populate persistent buffers
	int y = y0 + (threadIdx.y + blockIdx.y * 8);
	if (z0 == 0)
	{
		tzzbuf[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset+tzz_off] : 0.0f);
		txzbuf[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset+txz_off] : 0.0f);
		tyzbuf[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset+tyz_off] : 0.0f);
		if (threadIdx.y < 4)
		{
			tzzbuf[threadIdx.x+(3-threadIdx.y)*32] = -tzzbuf[threadIdx.x+(5+threadIdx.y)*32];
			txzbuf[threadIdx.x+(3-threadIdx.y)*32] = -txzbuf[threadIdx.x+(4+threadIdx.y)*32];
			tyzbuf[threadIdx.x+(3-threadIdx.y)*32] = -tyzbuf[threadIdx.x+(4+threadIdx.y)*32];
		}
		if (threadIdx.y == 4)
		{
			tzzbuf[threadIdx.x+4*32] = 0.0f;
		}
	}
	else
	{
		tzzbuf[threadIdx.x+threadIdx.y*32] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset-16+tzz_off] : 0.0f);
		txzbuf[threadIdx.x+threadIdx.y*32] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset-16+txz_off] : 0.0f);
		tyzbuf[threadIdx.x+threadIdx.y*32] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset-16+tyz_off] : 0.0f);
		tzzbuf[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset+tzz_off] : 0.0f);
		txzbuf[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset+txz_off] : 0.0f);
		tyzbuf[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset+tyz_off] : 0.0f);
	}

	for (int iZ = 0;  iZ < nz/8;  ++iZ)
	{
		int x = x0 + (threadIdx.x & 3);
		int z = z0 + iZ * 8 + (threadIdx.x / 4);

		float tmp3, tmp7, tmp8;
		if (z < vol_nz-8 && y <= m1_y1)
		{
			tmp3 = m1C[offset+tzz_off+32];
			tmp7 = m1C[offset+txz_off+32];
			tmp8 = m1C[offset+tyz_off+32];
		}
		else
		{
			tmp3 = tmp7 = tmp8 = 0.0f;
		}

		float tmp4, tmp5, txx_m4;
		if (m1L != 0L && y <= m1_y1)
		{
			tmp4 = m1L[offset+txy_off];
			tmp5 = m1L[offset+txz_off];
			txx_m4 = m1L[offset];
		}
		else
		{
			tmp4 = tmp5 = txx_m4 = 0.0f;
		}

		float tmp6, txy_p4;
		if (m1R != 0L && y <= m1_y1)
		{
			tmp6 = m1R[offset+txz_off];
			txy_p4 = m1R[offset+txy_off];
		}
		else
		{
			tmp6 = txy_p4 = 0.0f;
		}

		unsigned int em_word3 = (y <= y1) ? em[(threadIdx.y + blockIdx.y*8) * em_one_y_size_f + (iZ*32) + (z0*4) + threadIdx.x + 3*em_one_word_size_f] : 0;

		float txx_p0 = y <= m1_y1 ? m1C[offset] : 0.0f;
                float txy_p0 = y <= m1_y1 ? m1C[offset+txy_off] : 0.0f;

		float tmp2 = y <= m1_y1 ? m1C[offset+tyy_off] : 0.0f;

		float tmp1, tmp9, tmp10;
		if (threadIdx.y < 4)
		{
			if (y-4 >= m1_y0)
			{
				tmp1 = m1C[offset+tyy_off-4*one_y_size_f];
				tmp9 = m1C[offset+txy_off-4*one_y_size_f];
				tmp10 = m1C[offset+tyz_off-4*one_y_size_f];
			}
			else
			{
				tmp1 = tmp9 = tmp10 = 0.0f;
			}
		}
		else
		{
			if (y+4 <= m1_y1)
			{
				tmp1 = m1C[offset+tyy_off+4*one_y_size_f];
				tmp9 = m1C[offset+txy_off+4*one_y_size_f];
				tmp10 = m1C[offset+tyz_off+4*one_y_size_f];
			}
			else
			{
				tmp1 = tmp9 = tmp10 = 0.0f;
			}
		}

		// compute dxtxx
		buf[threadIdx.x+threadIdx.y*96] = txx_m4;
		buf[threadIdx.x+threadIdx.y*96+32] = txx_p0;
		buf[threadIdx.x+threadIdx.y*96+64] = m1R != 0L && y <= m1_y1 ? m1R[offset] : 0.0f;
		__syncthreads();
		float dxtxx = ( C0 * (txx_p0               - buf[cuCompTXXIdx(-1)]) + 
				C1 * (buf[cuCompTXXIdx(1)] - buf[cuCompTXXIdx(-2)]) + 
				C2 * (buf[cuCompTXXIdx(2)] - buf[cuCompTXXIdx(-3)]) +
				C3 * (buf[cuCompTXXIdx(3)] - txx_m4               ) ) * inv_DX;
		__syncthreads();  // wait for computes before reusing buf
		// TMJ 02/11/14 - Works
		//if (timestep == 1 && dxtxx != 0.0f) printf("TIMESTEP 1 :: DXTXX ( %d,%d,%d ) = %e\n",x,y,z,dxtxx);

		// compute dytyy
		buf[threadIdx.x+threadIdx.y*32+128] = tmp2;  // deposit middle section
		buf[threadIdx.x+threadIdx.y*32+64*(threadIdx.y&4)] = tmp1;
		__syncthreads();
		float dytyy = ( C0 * (buf[cuCompTYYIdx(1)] - buf[cuCompTYYIdx( 0)]) + 
				C1 * (buf[cuCompTYYIdx(2)] - buf[cuCompTYYIdx(-1)]) +
				C2 * (buf[cuCompTYYIdx(3)] - buf[cuCompTYYIdx(-2)]) +
				C3 * (buf[cuCompTYYIdx(4)] - buf[cuCompTYYIdx(-3)]) ) * inv_DY;
		// TMJ 02/11/14 - Works
		//if (timestep == 1 && dytyy != 0.0f) printf("TIMESTEP 1 :: DYTYY ( %d,%d,%d ) = %e\n",x,y,z,dytyy);

		// compute dztzz
		// ..load 8 next z and transpose to XYZ
		//float tmp3 = (iZ < ((nz/8)-1)) ? m1C[offset+tzz_off+32] : 0.0f;
		float v2 = cuTransposeXZY2XYZ(buf, tmp3);
		buf[threadIdx.x+threadIdx.y*32] = tzzbuf[threadIdx.x+threadIdx.y*32];  // copy 8 deepest z from tzz buf
		if (threadIdx.y < 4)
		{
			float v3 = tzzbuf[threadIdx.x+threadIdx.y*32+256];
			buf[threadIdx.x+threadIdx.y*32+256] = v3;  // copy 4 shallowest z from tzz buf
			buf[threadIdx.x+threadIdx.y*32+384] = v2;  // copy 4 deepest z from next block of tzz
			tzzbuf[threadIdx.x+threadIdx.y*32] = v3;  // shift tzzbuf by 8 z
		}
		// ..store next 8 z in tzzbuf
		__syncthreads(); // needed to prevent race condition
		tzzbuf[threadIdx.x+threadIdx.y*32+128] = v2;
		// note that we can use cuCompTYYIdx in place of cuCompTZZIdx after the transpose
		float dztzz = -( C0 * (buf[cuCompTYYIdx(1)] - buf[cuCompTYYIdx( 0)]) + 
				C1 * (buf[cuCompTYYIdx(2)] - buf[cuCompTYYIdx(-1)]) +
				C2 * (buf[cuCompTYYIdx(3)] - buf[cuCompTYYIdx(-2)]) +
				C3 * (buf[cuCompTYYIdx(4)] - buf[cuCompTYYIdx(-3)]) ) * inv_DZ;
		//float tmp4 = m1L != 0L ? m1L[offset+txy_off] : 0.0f;
		dztzz = cuTransposeXZY2XYZ(buf,dztzz);  // this actually transposes back from XYZ to XZY.
		// TMJ 02/11/14 - Works
		//if (timestep == 1 && dztzz != 0.0f) printf("TIMESTEP 1 :: DZTZZ ( %d,%d,%d ) = %e\n",x,y,z,dztzz);

		// compute dxtxy
                buf[threadIdx.x+threadIdx.y*96] = tmp4;
                buf[threadIdx.x+threadIdx.y*96+32] = txy_p0;
		buf[threadIdx.x+threadIdx.y*96+64] = txy_p4;
                __syncthreads();
                float dxtxy = ( C0 * (buf[cuCompTXXIdx(1)] - txy_p0               ) +
                                C1 * (buf[cuCompTXXIdx(2)] - buf[cuCompTXXIdx(-1)]) +
                                C2 * (buf[cuCompTXXIdx(3)] - buf[cuCompTXXIdx(-2)]) +
                                C3 * (txy_p4               - buf[cuCompTXXIdx(-3)]) ) * inv_DX;
		// TMJ 02/11/14 - Works
		//if (timestep == 1 && dxtxy != 0.0f) printf("TIMESTEP 1 :: DXTXY ( %d,%d,%d ) = %e\n",x,y,z,dxtxy);

		// ..compute dytxy
		float v4 = buf[threadIdx.x+threadIdx.y*96+32];  // read middle section for dytxy from shared memory
		__syncthreads();
		buf[threadIdx.x+threadIdx.y*32+128] = v4;  // deposit middle section
		buf[threadIdx.x+threadIdx.y*32+64*(threadIdx.y&4)] = tmp9;
		__syncthreads();
		float dytxy = ( C0 * (buf[cuCompTYYIdx(0)] - buf[cuCompTYYIdx(-1)]) +
                                C1 * (buf[cuCompTYYIdx(1)] - buf[cuCompTYYIdx(-2)]) +
                                C2 * (buf[cuCompTYYIdx(2)] - buf[cuCompTYYIdx(-3)]) +
                                C3 * (buf[cuCompTYYIdx(3)] - buf[cuCompTYYIdx(-4)]) ) * inv_DY;
		// TMJ 02/11/14 - Works
		//if (timestep == 1 && dytxy != 0.0f) printf("TIMESTEP 1 :: DYTXY ( %d,%d,%d ) = %e\n",x,y,z,dytxy);

		// compute dxtxz
		//float tmp5 = m1L != 0L ? m1L[offset+txz_off] : 0.0f;
		//float tmp6 = m1R != 0L ? m1R[offset+txz_off] : 0.0f;
                float txz_p0 = cuTransposeXZY2XYZ(buf, txzbuf[threadIdx.x+threadIdx.y*32+128]);  // read middle section from persistent txz buffer
		buf[threadIdx.x+threadIdx.y*96] = tmp5;
                buf[threadIdx.x+threadIdx.y*96+32] = txz_p0;
                float txz_p4 = tmp6;
		buf[threadIdx.x+threadIdx.y*96+64] = txz_p4;
                __syncthreads();
                float dxtxz = ( C0 * (buf[cuCompTXXIdx(1)] - txz_p0               ) +
                                C1 * (buf[cuCompTXXIdx(2)] - buf[cuCompTXXIdx(-1)]) +
                                C2 * (buf[cuCompTXXIdx(3)] - buf[cuCompTXXIdx(-2)]) +
                                C3 * (txz_p4               - buf[cuCompTXXIdx(-3)]) ) * inv_DX;
		// TMJ 02/11/14 - Works
		//if (timestep == 1 && dxtxz != 0.0f) printf("TIMESTEP 1 :: DXTXZ ( %d,%d,%d ) = %e\n",x,y,z,dxtxz);

		// ..compute dztxz
		//float tmp7 = (iZ < ((nz/8)-1)) ? m1C[offset+txz_off+32] : 0.0f;
		float v5 = cuTransposeXZY2XYZ(buf, tmp7);  // read next 8 z from gmem
		buf[threadIdx.x+threadIdx.y*32] = txzbuf[threadIdx.x+threadIdx.y*32];  // copy 8 deepest z from txz buf
                if (threadIdx.y < 4)
                {
                        float v6 = txzbuf[threadIdx.x+threadIdx.y*32+256];
                        buf[threadIdx.x+threadIdx.y*32+256] = v6;  // copy 4 shallowest z from txz buf
                        buf[threadIdx.x+threadIdx.y*32+384] = v5;  // copy 4 deepest z from next block of txz
                        txzbuf[threadIdx.x+threadIdx.y*32] = v6;  // shift txzbuf by 8 z
                }
                // ..store next 8 z in txzbuf
		__syncthreads();
                txzbuf[threadIdx.x+threadIdx.y*32+128] = v5;
                // note that we can use cuCompTYYIdx in place of cuCompTZZIdx after the transpose
                float dztxz = -( C0 * (buf[cuCompTYYIdx(0)] - buf[cuCompTYYIdx(-1)]) +
                                C1 * (buf[cuCompTYYIdx(1)] - buf[cuCompTYYIdx(-2)]) +
                                C2 * (buf[cuCompTYYIdx(2)] - buf[cuCompTYYIdx(-3)]) +
                                C3 * (buf[cuCompTYYIdx(3)] - buf[cuCompTYYIdx(-4)]) ) * inv_DZ;
                dztxz = cuTransposeXZY2XYZ(buf,dztxz);  // this actually transposes back from XYZ to XZY.
		// TMJ 02/11/14 - Works
		//if (timestep == 1 && dztxz != 0.0f) printf("TIMESTEP 1 :: DZTXZ ( %d,%d,%d ) = %e\n",x,y,z,dztxz);

		// compute dytyz
		float tyz_p0 = cuTransposeXZY2XYZ(buf, tyzbuf[threadIdx.x+threadIdx.y*32+128] );  // read middle section from persistent tyz buffer
		buf[threadIdx.x+threadIdx.y*32+128] = tyz_p0;
		buf[threadIdx.x+threadIdx.y*32+64*(threadIdx.y&4)] = tmp10;
		__syncthreads();
		float dytyz = ( C0 * (tyz_p0               - buf[cuCompTYYIdx(-1)]) + 
				C1 * (buf[cuCompTYYIdx(1)] - buf[cuCompTYYIdx(-2)]) +
				C2 * (buf[cuCompTYYIdx(2)] - buf[cuCompTYYIdx(-3)]) +
				C3 * (buf[cuCompTYYIdx(3)] - buf[cuCompTYYIdx(-4)]) ) * inv_DY;
		// TMJ 02/11/14 - Works
		//if (timestep == 1 && dytyz != 0.0f) printf("TIMESTEP 1 :: DYTYZ ( %d,%d,%d ) = %e\n",x,y,z,dytyz);
		/*
		{
			int x = x0 + (threadIdx.x & 3);
			int y = y0 + (threadIdx.y + blockIdx.y * 8);
			int z = iZ * 8 + (threadIdx.x / 4);
			if (x == 200 && y == 200 && z == 200)
			{
				for (int i = -4;  i <= 3;  ++i)
				{
					printf("tyz[200,%d,200] = %e\n",200+i,buf[cuCompTYYIdx(i)]);
				}
				printf("tyz_p0 = %e\n",tyz_p0);
			}
		}
		*/

		// ..compute dztyz
		//float tmp8 = (iZ < ((nz/8)-1)) ? m1C[offset+tyz_off+32] : 0.0f;
		float v8 = cuTransposeXZY2XYZ(buf, tmp8);  // read next 8 z from gmem
		buf[threadIdx.x+threadIdx.y*32] = tyzbuf[threadIdx.x+threadIdx.y*32];  // copy 8 deepest z from tyz buf
                if (threadIdx.y < 4)
                {
                        float v9 = tyzbuf[threadIdx.x+threadIdx.y*32+256];
                        buf[threadIdx.x+threadIdx.y*32+256] = v9;  // copy 4 shallowest z from tyz buf
                        buf[threadIdx.x+threadIdx.y*32+384] = v8;  // copy 4 deepest z from next block of tyz
                        tyzbuf[threadIdx.x+threadIdx.y*32] = v9;  // shift tyzbuf by 8 z
                }
                // ..store next 8 z in tyzbuf
                __syncthreads();
                tyzbuf[threadIdx.x+threadIdx.y*32+128] = v8;
                // note that we can use cuCompTYYIdx in place of cuCompTZZIdx after the transpose
                float dztyz = -( C0 * (buf[cuCompTYYIdx(0)] - buf[cuCompTYYIdx(-1)]) +
                                C1 * (buf[cuCompTYYIdx(1)] - buf[cuCompTYYIdx(-2)]) +
                                C2 * (buf[cuCompTYYIdx(2)] - buf[cuCompTYYIdx(-3)]) +
                                C3 * (buf[cuCompTYYIdx(3)] - buf[cuCompTYYIdx(-4)]) ) * inv_DZ;
                dztyz = cuTransposeXZY2XYZ(buf,dztyz);  // this actually transposes back from XYZ to XZY.
		// TMJ 02/11/14 - Works
		//if (timestep == 1 && dztyz != 0.0f) printf("TIMESTEP 1 :: DZTYZ ( %d,%d,%d ) = %e\n",x,y,z,dztyz);

		if (y <= y1)
		{
			// get word3 from earth model
			float Q, Density;
			cuUnpack_Q_Density(em_word3,Q_min,Q_range,&Q,Density_min,Density_range,&Density);
			Q = 1.0f / Q;  // compressed model actually stores inverse of Q.

			// ..compute itausig and difitau
			float wq = 6.2831853072f * fq;
			float te = (1.0f + sqrtf(1.0f + Q*Q)) / (Q*wq);
			float tt = 1.0f / (te * wq * wq);
			float itausig = 1.0f / tt;
			float difitau = ((1.0f / te) - itausig);

			// Update viscoelastic(SLS) vector field:
			float const1 = 1.0f / (1.0f + 0.5f*dti*itausig);
			float const2 = (1.0f - 0.5f*dti*itausig);
			float const3 = dti*difitau;

			float old_sx = m2C[offset+3*one_wf_size_f];
			float old_sy = m2C[offset+4*one_wf_size_f];
			float old_sz = m2C[offset+5*one_wf_size_f];

			float sx = const3*(dxtxx + dytxy + dztxz);
			sx = sx + const2*old_sx;
			sx = const1*sx;

			float sy = const3*(dxtxy + dytyy + dztyz);
			sy = sy + const2*old_sy;
			sy = const1*sy;

			float sz = const3*(dxtxz + dytyz + dztzz);
			sz = sz + const2*old_sz;
			sz = const1*sz;

			//if (timestep == 1 && sx != 0.0f) printf("TIMESTEP 1 :: SX ( %d,%d,%d ) = %e\n",x,y,z,sx);
			//if (timestep == 1 && sy != 0.0f) printf("TIMESTEP 1 :: SY ( %d,%d,%d ) = %e\n",x,y,z,sy);
			//if (timestep == 1 && sz != 0.0f) printf("TIMESTEP 1 :: SZ ( %d,%d,%d ) = %e\n",x,y,z,sz);

			cmp[offset+3*one_wf_size_f] = sx;
			cmp[offset+4*one_wf_size_f] = sy;
			cmp[offset+5*one_wf_size_f] = sz;

			// Absorbing boundary decay funct (for Maxwell viscoelastic model):
			//int x = x0 + (threadIdx.x & 3);
			//int y = y0 + (threadIdx.y + blockIdx.y * 8);
			//int z = iZ * 8 + (threadIdx.x / 4);
			float deta = Compute_ABC(x,y,z,vol_nx,vol_ny,vol_nz,nabc_top,nabc_bot,nabc_sdx,nabc_sdy,vpvert_avtop,vpvert_avbot,inv_DX,inv_DY,inv_DZ);
			float dabc = (1.0f - 0.5f*deta*dti) / (1.0f + 0.5f*deta*dti);

			//if (x == 400 && y == 400 && z == 400)
			//{
			//	printf("Q=%e, Density=%e, wq=%e, te=%e, tt=%e, itausig=%e, difitau=%e, const1=%e, const2=%e, const3=%e, deta=%e, dabc=%e\n",
			//		Q,Density,wq,te,tt,itausig,difitau,const1,const2,const3,deta,dabc);
			//}

			// Update viscoelastic particle velocities:
			float old_vx = m2C[offset];
			float old_vy = m2C[offset+one_wf_size_f];
			float old_vz = m2C[offset+2*one_wf_size_f];

			float factor = dti / Density;

			float vx = dabc*old_vx + factor*dxtxx + factor*dytxy + factor*dztxz + factor*sx;
			float vy = dabc*old_vy + factor*dxtxy + factor*dytyy + factor*dztyz + factor*sy;
			float vz = dabc*old_vz + factor*dxtxz + factor*dytyz + factor*dztzz + factor*sz;

			//if (timestep == 1 && vx != 0.0f) printf("TIMESTEP 1 :: VX ( %d,%d,%d ) = %e\n",x,y,z,vx);
			//if (timestep == 1 && vy != 0.0f) printf("TIMESTEP 1 :: VY ( %d,%d,%d ) = %e\n",x,y,z,vy);
			//if (timestep == 1 && vz != 0.0f) printf("TIMESTEP 1 :: VZ ( %d,%d,%d ) = %e\n",x,y,z,vz);

			cmp[offset] = vx;
			cmp[offset+one_wf_size_f] = vy;
			cmp[offset+2*one_wf_size_f] = vz;

			/*
			   if (x == 501 && y == 401 && z == 401)
			   {
			   printf("\nPropagate_Particle_Velocities\n");
			   printf("-----------------------------\n");
			   printf("timestep = %d\n",timestep);
			   printf("dti=%e\n",dti);
			   printf("dxtxx=%e\n",dxtxx);
			   printf("dytxy=%e\n",dytxy);
			   printf("dztxz=%e\n",dztxz);
			   printf("dxtxy=%e\n",dxtxy);
			   printf("dytyy=%e\n",dytyy);
			   printf("dztyz=%e\n",dztyz);
			   printf("dxtxz=%e\n",dxtxz);
			   printf("dytyz=%e\n",dytyz);
			   printf("dztzz=%e\n",dztzz);
			   printf("wq=%e\n",wq);
			   printf("te=%e\n",te);
			   printf("tt=%e\n",tt);
			   printf("itausig=%e\n",itausig);
			   printf("difitau=%e\n",difitau);
			   printf("const1=%e\n",const1);
			   printf("const2=%e\n",const2);
			   printf("const3=%e\n",const3);

			   printf("rho=%e\n",Density);
			   printf("dabc=%e\n",dabc);
			   printf("deta=%e\n",deta);

			   printf("sx=%e\n",sx);
			   printf("sy=%e\n",sy);
			   printf("sz=%e\n",sz);
			   printf("vx=%e\n",vx);
			   printf("vy=%e\n",vy);
			   printf("vz=%e\n",vz);

			   printf("\n");
			   }
			 */
		}

		// increase offsets
		offset += 32;
	}
}

#ifdef GPU_DEBUG
__global__ 
void 
cuNon_Zeros_Kernel(
	int x0,			// x coordinate of westernmost coordinate in block
	int y0,			// y coordinate of southernmost coordinate in block
	float* m1L,		// txx, tyy, tzz, txy, txz and tyz in that order. left halo, t(0), y(0)
        float* m1C,		// ..middle, t(0), y(0)
        float* m1R,		// ..right halo, t(0), y(0)
	int nx,
	int ny,
	int nz
	)
{
	int one_wf_size_f = 4 * nz;
	int one_y_size_f = one_wf_size_f * 6;
	float* p = m1C + one_y_size_f * (threadIdx.y + blockIdx.y * 8) + threadIdx.x;
	for (int iZ = 0;  iZ < nz/8;  ++iZ)
	{
		float txx = p[iZ*32];
		if (txx != 0.0f)
		{
			int x = x0 + (threadIdx.x & 3);
			int y = threadIdx.y + blockIdx.y * 8;
			int z = iZ * 8 + (threadIdx.x >> 2);
			printf("txx[%d,%d,%d] = %f\n",x,y,z,txx);
		}
		float tyy = p[iZ*32+one_wf_size_f];
		if (tyy != 0.0f)
		{
			int x = x0 + (threadIdx.x & 3);
			int y = threadIdx.y + blockIdx.y * 8;
			int z = iZ * 8 + (threadIdx.x >> 2);
			printf("tyy[%d,%d,%d] = %f\n",x,y,z,tyy);
		}
		float tzz = p[iZ*32+2*one_wf_size_f];
		if (tzz != 0.0f)
		{
			int x = x0 + (threadIdx.x & 3);
			int y = threadIdx.y + blockIdx.y * 8;
			int z = iZ * 8 + (threadIdx.x >> 2);
			printf("tzz[%d,%d,%d] = %f\n",x,y,z,tzz);
		}
		float txy = p[iZ*32+3*one_wf_size_f];
		if (txy != 0.0f)
		{
			int x = x0 + (threadIdx.x & 3);
			int y = threadIdx.y + blockIdx.y * 8;
			int z = iZ * 8 + (threadIdx.x >> 2);
			printf("txy[%d,%d,%d] = %f\n",x,y,z,txy);
		}
		float txz = p[iZ*32+4*one_wf_size_f];
		if (txz != 0.0f)
		{
			int x = x0 + (threadIdx.x & 3);
			int y = threadIdx.y + blockIdx.y * 8;
			int z = iZ * 8 + (threadIdx.x >> 2);
			printf("txz[%d,%d,%d] = %f\n",x,y,z,txz);
		}
		float tyz = p[iZ*32+5*one_wf_size_f];
		if (tyz != 0.0f)
		{
			int x = x0 + (threadIdx.x & 3);
			int y = threadIdx.y + blockIdx.y * 8;
			int z = iZ * 8 + (threadIdx.x >> 2);
			printf("tyz[%d,%d,%d] = %f\n",x,y,z,tyz);
		}
	}
}
#endif

//
// Relative Y is the Y current Y coordinate relative to the first Y position in block.
// Relative Y can be negative if the block lacks Y halo on low side.
// 

//
// Wavefields are interleaved in the following order: X-Z-WF-Y
// The particle velocities are stored in this order: Vx, Vy, Vz, Sx, Sy, Sz
// The strain rates are stored in this order: txx, tyy, tzz, txy, txz, tyz
//

void 
Host_Propagate_Particle_Velocities_Kernel(
	int timestep,
	cudaStream_t stream,
	int num_z,		// number of blocks along z axis
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
	)
{
	int one_wf_size = one_y_size / 6;
	int em_one_word_size = one_wf_size;
	int em_one_y_size = em_one_word_size * 4;

	const int tyy_off = one_wf_size / 4;
	const int tzz_off = 2 * tyy_off;
	const int txy_off = 3 * tyy_off;
	const int txz_off = 4 * tyy_off;
	const int tyz_off = 5 * tyy_off;

	int nx = 4;
	int ny = y1 - y0 + 1;
	int nz = vol_nz;

	dim3 blockShape(32,8,1);
	dim3 gridShape(1,(ny+7)/8,num_z);

	//cuNon_Zeros_Kernel<<<gridShape,blockShape,0,stream>>>(x0,y0,(float*)m1L,(float*)m1C,(float*)m1R,nx,ny,nz);

	cuPropagate_Particle_Velocities_Kernel<<<gridShape,blockShape,0,stream>>>(
		timestep,
		x0,y0,y1,m1_y0,m1_y1,vol_nx,vol_ny,vol_nz,dti,
		(unsigned int*)em,(float*)cmp,(float*)m1L,(float*)m1C,(float*)m1R,(float*)m2C,
		C0,C1,C2,C3,inv_DX,inv_DY,inv_DZ,
		vpvert_avtop,vpvert_avbot,nabc_sdx,nabc_sdy,nabc_top,nabc_bot,Q_min,Q_range/255.0f,fq,Density_min,Density_range/255.0f,
		one_wf_size/4,one_y_size/4,em_one_word_size/4,em_one_y_size/4,
		tyy_off,tzz_off,txy_off,txz_off,tyz_off);
#ifdef GPU_DEBUG
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
#endif

	//
	// add source term(s)
	//
	if (inject_source && (is_force || is_velocity))
	{
		if (source_interpolation_method == Point)
		{
			dim3 blockShape3(1,1,1);
			dim3 gridShape3(1,1,1);
			cuApply_Point_Source_To_VxVyVz<<<gridShape3,blockShape3,0,stream>>>((int*)em,Q_min,Q_range,Density_min,Density_range,(float*)cmp,x0,y0,nx,ny,nz,dti,is_force,ampl1,ampl2,ampl3,xsou,ysou,zsou,svaw_sample);
		}
		else if(source_interpolation_method == Trilinear)
		{
			dim3 blockShape3(1,1,1);
			dim3 gridShape3(1,1,1);
			cuApply_Trilinear_Source_To_VxVyVz<<<gridShape3,blockShape3,0,stream>>>((int*)em,Q_min,Q_range,Density_min,Density_range,(float*)cmp,x0,y0,nx,ny,nz,dti,is_force,ampl1,ampl2,ampl3,xsou,ysou,zsou,svaw_sample);
		}
		else if(source_interpolation_method == Sinc)
		{
			// use only one thread along z to prevent possible race condition
			dim3 blockShape2(8,8,1);
			dim3 gridShape2(1,1,1);
			cuApply_Source_Term_To_VxVyVz<<<gridShape2,blockShape2,0,stream>>>(em,Q_min,Q_range,Density_min,Density_range,cmp,x0,y0,0,nx,ny,nz,dti,is_force,ampl1,ampl2,ampl3,xsou,ysou,zsou,svaw_sample);
		}
	}
#ifdef GPU_DEBUG
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
#endif
}

