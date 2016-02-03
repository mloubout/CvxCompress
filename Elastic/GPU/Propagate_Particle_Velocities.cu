#include <cuda_runtime_api.h>
#include "Elastic_Interpolation.hxx"

//
// CUDA kernel that propagates particle velocity wavefield and memory variable.
//
// 
//

__device__ __forceinline__
int cuCompTXXIdx(int offset)
{
        int abs_offset = offset + (threadIdx.x & 3) + 4;
        int quotient = abs_offset / 4;
        int remainder = abs_offset & 3;
        return threadIdx.y*96 + quotient*32 + (threadIdx.x&28) + remainder;
}

__device__ __forceinline__ 
int cuCompTXXIdx_16(int offset)
{
	int abs_offset = offset + (threadIdx.x & 7) + 8;
	int quotient = abs_offset / 8;
	int remainder = abs_offset & 7;
	return threadIdx.y*192 + quotient*64 + (threadIdx.x&56) + remainder;
}

__device__ __forceinline__
int cuCompTYYIdx(int offset)
{
        return (offset + 4 + threadIdx.y) * 32 + threadIdx.x;
}

__device__ __forceinline__
int cuCompTYYIdx_16(int offset)
{
	return (offset + 8 + threadIdx.y) * 64 + threadIdx.x;
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
float cuTransposeXZY2XYZ_16(float* buf, float v)
{
	__syncthreads();  // wait for previous step to finish using buf
	buf[threadIdx.x+((threadIdx.x&56)*8)+threadIdx.y*8] = v;
	__syncthreads();  // wait for all threads to finish writing to buf
	float retval = buf[threadIdx.x+threadIdx.y*72];
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

	float AX = fabsf(X);
	if (AX < 3.75f)
	{
		float Y = (X*X) / (3.75f*3.75f);
		return P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))));
	}
	else
	{
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
	float fsinc_x = (x_cell == 0.0f) ? 1.0f : win_x * sinf(x_cell*pi)/(x_cell*pi);

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
}

__device__ 
void __cuApply_Source_Term_To_VxVyVz(
	int thr_z,
	unsigned int* em,
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
	int kcell,
	bool is_p_reciprocity,
	float bmod_ref,
	float rho_ref,
	bool source_ghost_enabled,
	int ghost_sea_surface_z
	)
{
	int my_x = icell + threadIdx.x - 3 - x0;
	int my_y = jcell + threadIdx.y - 3 - y0;
	int my_z = kcell + thr_z       - 3;
	
	if (
			(my_x >= 0 && my_x < nx) && 
			(my_y >= 0 && my_y < ny) && 
			(my_z >= ghost_sea_surface_z && my_z < nz) // TMJ 04/08/15 Turned off mirroring, it's not in Kurt Nihei's original code
	   )
	{
		// ..fractional distance from grid pt to sou:
		float dx_frac = xs - (float)icell;
		float dy_frac = ys - (float)jcell;
                float dz_frac = zs - (float)kcell;

		// (fx/vx sou needs to be shifted +0.5icell to colloc w/ pr)
		float vx_dx_frac = xs + 0.5f - (float)icell;
		// (fy/vy sou needs to be shifted -0.5dy to colloc w/ pr)
		float vy_dy_frac = ys - 0.5f - (float)jcell;
		// (fz/vz sou need to be shifted -0.5kcell to colloc w/ pr)
		float vz_dz_frac = zs - 0.5f - (float)kcell;

		float vx_fsinc = cuGen_Sinc_Weight(threadIdx.x,threadIdx.y,thr_z,vx_dx_frac,dy_frac,dz_frac);
		float vy_fsinc = cuGen_Sinc_Weight(threadIdx.x,threadIdx.y,thr_z,dx_frac,vy_dy_frac,dz_frac);
		float vz_fsinc = cuGen_Sinc_Weight(threadIdx.x,threadIdx.y,thr_z,dx_frac,dy_frac,vz_dz_frac);

		if (vx_fsinc != 0.0f || vy_fsinc != 0.0f || vz_fsinc != 0.0f)
		{
			// mirror source if necessary
			//my_z = my_z < 0 ? -my_z : my_z;

			int one_wf_size_f = nx * nz;
			int one_y_size_f = one_wf_size_f * 6;
			int idx = my_x + my_y * one_y_size_f + my_z * nx;

			int my_z_ghost = 2 * ghost_sea_surface_z - my_z;
			int idx_ghost = my_x + my_y * one_y_size_f + my_z_ghost * nx;
			int my_z_vz_ghost = my_z_ghost - 1;
			int idx_vz_ghost = my_x + my_y * one_y_size_f + my_z_vz_ghost * nx;

			int em_one_word_size_f = one_wf_size_f;
			int em_one_y_size_f = em_one_word_size_f * 4;

			if (is_force)
			{
				int em_word3 = em[my_x+my_y*em_one_y_size_f+my_z*nx+3*em_one_word_size_f];
				float Density;
				cuUnpack_Density(em_word3,Density_min,Density_range,&Density);

				cmp[idx                ] = cmp[idx                ] + vx_fsinc * dti * (val * ampl1) / Density;
				cmp[idx+  one_wf_size_f] = cmp[idx+  one_wf_size_f] + vy_fsinc * dti * (val * ampl2) / Density;
				cmp[idx+2*one_wf_size_f] = cmp[idx+2*one_wf_size_f] + vz_fsinc * dti * (val * ampl3) / Density;

				//printf("Adding source term (%f * [%f-%f-%f] * %f) to Vx,Vy,Vz at %d,%d,%d\n",dti,vx_fsinc,vy_fsinc,vz_fsinc,val/Density,my_x+x0,my_y+y0,my_z);

				if (source_ghost_enabled)
				{
					int em_word3_ghost = em[my_x+my_y*em_one_y_size_f+my_z_ghost*nx+3*em_one_word_size_f];
					float Density_ghost;
					cuUnpack_Density(em_word3_ghost,Density_min,Density_range,&Density_ghost);
				
					cmp[idx_ghost                   ] = cmp[idx_ghost                   ] - vx_fsinc * dti * (val * ampl1) / Density_ghost;
					cmp[idx_ghost   +  one_wf_size_f] = cmp[idx_ghost   +  one_wf_size_f] - vy_fsinc * dti * (val * ampl2) / Density_ghost;
					cmp[idx_vz_ghost+2*one_wf_size_f] = cmp[idx_vz_ghost+2*one_wf_size_f] + vz_fsinc * dti * (val * ampl3) / Density_ghost;
					// TMJ, Vz requires a POSITIVE ghost, hence + sign. Not a bug.
				}
			}
			else
			{
				int em_word3 = em[my_x+my_y*em_one_y_size_f+my_z*nx+3*em_one_word_size_f];
				float rho;
				cuUnpack_Density(em_word3,Density_min,Density_range,&rho);
				float scale_sou = 1.0f;
				if (is_p_reciprocity)
				{
					scale_sou = -1.0f / (rho * bmod_ref);
				}
				else
				{
					scale_sou = rho_ref / rho;
				}

				//if (threadIdx.x == 3 && threadIdx.y == 3)
				//{
				//	printf("\n*****\nmy_z= %d, thr_z=%d :: rho = %e, bmod_ref = %e, scale_sou = %e, vz_fsinc = %e\n*****\n\n",my_z,thr_z,rho,bmod_ref,scale_sou,vz_fsinc);
				//}

				cmp[idx                ] = cmp[idx                ] + vx_fsinc * dti * scale_sou * val * ampl1;
				cmp[idx+  one_wf_size_f] = cmp[idx+  one_wf_size_f] + vy_fsinc * dti * scale_sou * val * ampl2;
				cmp[idx+2*one_wf_size_f] = cmp[idx+2*one_wf_size_f] + vz_fsinc * dti * scale_sou * val * ampl3;

				//printf("Adding source term (%f * [%f-%f-%f] * %f) to Vx,Vy,Vz at %d,%d,%d\n",dti,vx_fsinc,vy_fsinc,vz_fsinc,val,my_x+x0,my_y+y0,my_z);

				if (source_ghost_enabled)
				{
					int em_word3_ghost = em[my_x+my_y*em_one_y_size_f+my_z_ghost*nx+3*em_one_word_size_f];
					float rho_ghost;
					cuUnpack_Density(em_word3_ghost,Density_min,Density_range,&rho_ghost);
					float scale_sou_ghost = 1.0f;
					if (is_p_reciprocity)
					{
						scale_sou_ghost = -1.0f / (rho_ghost * bmod_ref);
					}
					else
					{
						scale_sou_ghost = rho_ref / rho_ghost;
					}

					cmp[idx_ghost                   ] = cmp[idx_ghost                   ] - vx_fsinc * dti * scale_sou_ghost * val * ampl1;
					cmp[idx_ghost   +  one_wf_size_f] = cmp[idx_ghost   +  one_wf_size_f] - vy_fsinc * dti * scale_sou_ghost * val * ampl2;
					cmp[idx_vz_ghost+2*one_wf_size_f] = cmp[idx_vz_ghost+2*one_wf_size_f] + vz_fsinc * dti * scale_sou_ghost * val * ampl3;
					// TMJ, Vz requires a POSITIVE ghost, hence + sign. Not a bug.
				}
			}
		}
	}
}

__device__ 
void _cuApply_Source_Term_To_VxVyVz(
	void* em,
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
        float val,
	bool is_p_reciprocity,
	float bmod_ref,
	float rho_ref,
	bool source_ghost_enabled,
	int ghost_sea_surface_z
	)
{
	// fx/vx contribution:

	// nearest grid point:
	int icell = (int)lrintf(xs);
	int jcell = (int)lrintf(ys);
	int kcell = (int)lrintf(zs);

	for (int thr_z = 0;  thr_z < 8;  ++thr_z)
	{
		__cuApply_Source_Term_To_VxVyVz(
				thr_z,(unsigned int*)em,Density_min,Density_range,
				(float*)cmp,x0,y0,z0,nx,ny,nz,dti,is_force,ampl1,ampl2,ampl3,
				xs,ys,zs,val,icell,jcell,kcell,is_p_reciprocity,bmod_ref,rho_ref,source_ghost_enabled,ghost_sea_surface_z);
	}
}

__global__ 
void cuApply_Source_Term_To_VxVyVz(
	void* em,
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
        float val,
	bool source_ghost_enabled,
	int ghost_sea_surface_z,
	bool is_p_reciprocity,
	float bmod_ref,
	float rho_ref
	)
{
	_cuApply_Source_Term_To_VxVyVz(em,Density_min,Density_range,cmp,x0,y0,z0,nx,ny,nz,dti,is_force,ampl1,ampl2,ampl3,xs,ys,zs,val,is_p_reciprocity,bmod_ref,rho_ref,source_ghost_enabled,ghost_sea_surface_z);
}

__device__ 
void cuAdd_Source_Term_To_Single_VxVyVz(
	int* em,
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
	int em_one_y_size_f,
	bool is_p_reciprocity,
	float bmod_ref,
	float rho_ref
        )
{
	if (delta_Vx != 0.0f && Vx_ix >= 0 && Vx_ix < nx && iy >= 0 && iy < ny && iz >= 0 && iz < nz)
	{
		int idx_Vx = Vx_ix +    iy * one_y_size_f +    iz * nx;

		int em_word3 = em[Vx_ix+iy*em_one_y_size_f+iz*nx+3*em_one_word_size_f];
		if (is_force)
		{
			float Density;
			cuUnpack_Density(em_word3,Density_min,Density_range,&Density);
			cmp[idx_Vx] = cmp[idx_Vx] + delta_Vx / Density;
		}
		else
		{
			float rho;
			cuUnpack_Density(em_word3,Density_min,Density_range,&rho);
			float scale_sou = 1.0f;
			if (is_p_reciprocity)
			{
				scale_sou = -1.0f / (rho * bmod_ref);
			}
			else
			{
				scale_sou = rho_ref / rho;
			}
			cmp[idx_Vx] = cmp[idx_Vx] + delta_Vx * scale_sou;
		}
	}
	if (delta_Vy != 0.0f && ix >= 0 && ix < nx && Vy_iy >= 0 && Vy_iy < ny && iz >= 0 && iz < nz)
	{
		int idx_Vy =    ix + Vy_iy * one_y_size_f +    iz * nx +     one_wf_size_f;

		int em_word3 = em[ix+Vy_iy*em_one_y_size_f+iz*nx+3*em_one_word_size_f];
		if (is_force)
		{
			float Density;
			cuUnpack_Density(em_word3,Density_min,Density_range,&Density);
			cmp[idx_Vy] = cmp[idx_Vy] + delta_Vy / Density;
		}
		else
		{
			float rho;
			cuUnpack_Density(em_word3,Density_min,Density_range,&rho);
			float scale_sou = 1.0f;
			if (is_p_reciprocity)
			{
				scale_sou = -1.0f / (rho * bmod_ref);
			}
			else
			{
				scale_sou = rho_ref / rho;
			}
			cmp[idx_Vy] = cmp[idx_Vy] + delta_Vy * scale_sou;
		}
	}
	if (delta_Vz != 0.0f && ix >= 0 && ix < nx && iy >= 0 && iy < ny && Vz_iz >= 0 && Vz_iz < nz)
	{
		int idx_Vz =    ix +    iy * one_y_size_f + Vz_iz * nx + 2 * one_wf_size_f;

		int em_word3 = em[ix+iy*em_one_y_size_f+Vz_iz*nx+3*em_one_word_size_f];
		if (is_force)
		{
			float Density;
			cuUnpack_Density(em_word3,Density_min,Density_range,&Density);
			cmp[idx_Vz] = cmp[idx_Vz] + delta_Vz / Density;
		}
		else
		{
			float rho;
			cuUnpack_Density(em_word3,Density_min,Density_range,&rho);
			float scale_sou = 1.0f;
			if (is_p_reciprocity)
			{
				scale_sou = -1.0f / (rho * bmod_ref);
			}
			else
			{
				scale_sou = rho_ref / rho;
			}
			cmp[idx_Vz] = cmp[idx_Vz] + delta_Vz * scale_sou;
			//printf("delta_Vz = %e, ix=%d, iy=%d, Vz_iz=%d\n",delta_Vz,ix,iy,Vz_iz);
		}
	}
}

__device__ 
void _cuApply_Point_Source_To_VxVyVz(
	int* em,
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
        float val,
	bool is_p_reciprocity,
	float bmod_ref,
	float rho_ref
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
		em,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,dti*val*ampl1,dti*val*ampl2,dti*val*ampl3,ix,iy,iz,ix,iy,iz,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f,
		is_p_reciprocity,bmod_ref,rho_ref);
}

__global__ 
void cuApply_Point_Source_To_VxVyVz(
	int* em,
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
        float val,
	bool source_ghost_enabled,
	float ghost_sea_surface,
	bool is_p_reciprocity,
	float bmod_ref,
	float rho_ref
        )
{
	_cuApply_Point_Source_To_VxVyVz(em,Density_min,Density_range,cmp,x0,y0,nx,ny,nz,dti,is_force,ampl1,ampl2,ampl3,xs,ys,zs,val,is_p_reciprocity,bmod_ref,rho_ref);
	if (source_ghost_enabled)
	{
		// TMJ, Vz requires a POSITIVE ghost, which I make happen by flipping the sign of ampl3.
		_cuApply_Point_Source_To_VxVyVz(em,Density_min,Density_range,cmp,x0,y0,nx,ny,nz,dti,is_force,ampl1,ampl2,-ampl3,xs,ys,2.0f*ghost_sea_surface-zs,-val,is_p_reciprocity,bmod_ref,rho_ref);
	}
}

__device__ 
void _cuApply_Trilinear_Source_To_VxVyVz(
	int* em,
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
        float val,
	bool is_p_reciprocity,
	float bmod_ref,
	float rho_ref
        )
{
	int ix = (int)truncf(xs) - x0;
	int iy = (int)truncf(ys) - y0;
	int iz = (int)truncf(zs);

	float xd = 1.0f - (xs - (float)(ix+x0));
	float yd = 1.0f - (ys - (float)(iy+y0));
	float zd = 1.0f - (zs - (float)iz);

	float Vx_xs = xs + 0.5f;
	float Vy_ys = ys - 0.5f;
	float Vz_zs = zs - 0.5f;

	int Vx_ix = (int)truncf(Vx_xs) - x0;
	int Vy_iy = (int)truncf(Vy_ys) - y0;
	int Vz_iz = (int)truncf(Vz_zs);

	float Vx_xd = 1.0f - (Vx_xs - (float)(Vx_ix+x0));
	float Vy_yd = 1.0f - (Vy_ys - (float)(Vy_iy+y0));
	float Vz_zd = 1.0f - (Vz_zs - (float)Vz_iz);

	int one_wf_size_f = nx * nz;
	int one_y_size_f = one_wf_size_f * 6;

	int em_one_word_size_f = one_wf_size_f;
	int em_one_y_size_f = em_one_word_size_f * 4;

	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*Vx_xd*   yd*   zd,
		dti*val*ampl2*   xd*Vy_yd*   zd,
		dti*val*ampl3*   xd*   yd*Vz_zd,
		   ix,   iy,   iz,
		Vx_ix,Vy_iy,Vz_iz,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f,
		is_p_reciprocity,bmod_ref,rho_ref);
	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*(1.0f-Vx_xd)*   yd*   zd,
		dti*val*ampl2*(1.0f-   xd)*Vy_yd*   zd,
		dti*val*ampl3*(1.0f-   xd)*   yd*Vz_zd,
		   ix+1,   iy,   iz,
		Vx_ix+1,Vy_iy,Vz_iz,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f,
		is_p_reciprocity,bmod_ref,rho_ref);
	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*(1.0f-Vx_xd)*(1.0f-   yd)*   zd,
		dti*val*ampl2*(1.0f-   xd)*(1.0f-Vy_yd)*   zd,
		dti*val*ampl3*(1.0f-   xd)*(1.0f-   yd)*Vz_zd,
		   ix+1,   iy+1,   iz,
		Vx_ix+1,Vy_iy+1,Vz_iz,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f,
		is_p_reciprocity,bmod_ref,rho_ref);
	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*Vx_xd*(1.0f-   yd)*   zd,
		dti*val*ampl2*   xd*(1.0f-Vy_yd)*   zd,
		dti*val*ampl3*   xd*(1.0f-   yd)*Vz_zd,
		   ix,   iy+1,   iz,
		Vx_ix,Vy_iy+1,Vz_iz,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f,
		is_p_reciprocity,bmod_ref,rho_ref);
	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*Vx_xd*   yd*(1.0f-   zd),
		dti*val*ampl2*   xd*Vy_yd*(1.0f-   zd),
		dti*val*ampl3*   xd*   yd*(1.0f-Vz_zd),
		   ix,   iy,   iz+1,
		Vx_ix,Vy_iy,Vz_iz+1,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f,
		is_p_reciprocity,bmod_ref,rho_ref);
	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*(1.0f-Vx_xd)*   yd*(1.0f-   zd),
		dti*val*ampl2*(1.0f-   xd)*Vy_yd*(1.0f-   zd),
		dti*val*ampl3*(1.0f-   xd)*   yd*(1.0f-Vz_zd),
		   ix+1,   iy,   iz+1,
		Vx_ix+1,Vy_iy,Vz_iz+1,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f,
		is_p_reciprocity,bmod_ref,rho_ref);
	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*(1.0f-Vx_xd)*(1.0f-   yd)*(1.0f-   zd),
		dti*val*ampl2*(1.0f-   xd)*(1.0f-Vy_yd)*(1.0f-   zd),
		dti*val*ampl3*(1.0f-   xd)*(1.0f-   yd)*(1.0f-Vz_zd),
		   ix+1,   iy+1,   iz+1,
		Vx_ix+1,Vy_iy+1,Vz_iz+1,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f,
		is_p_reciprocity,bmod_ref,rho_ref);
	cuAdd_Source_Term_To_Single_VxVyVz(
		em,Density_min,Density_range,
		cmp,nx,ny,nz,
		is_force,
		dti*val*ampl1*Vx_xd*(1.0f-   yd)*(1.0f-   zd),
		dti*val*ampl2*   xd*(1.0f-Vy_yd)*(1.0f-   zd),
		dti*val*ampl3*   xd*(1.0f-   yd)*(1.0f-Vz_zd),
		   ix,   iy+1,   iz+1,
		Vx_ix,Vy_iy+1,Vz_iz+1,
		one_wf_size_f,one_y_size_f,em_one_word_size_f,em_one_y_size_f,
		is_p_reciprocity,bmod_ref,rho_ref);
}

__global__ 
void cuApply_Trilinear_Source_To_VxVyVz(
	int* em,
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
        float val,
	bool source_ghost_enabled,
	float ghost_sea_surface,
	bool is_p_reciprocity,
	float bmod_ref,
	float rho_ref
        )
{
	_cuApply_Trilinear_Source_To_VxVyVz(em,Density_min,Density_range,cmp,x0,y0,nx,ny,nz,dti,is_force,ampl1,ampl2,ampl3,xs,ys,zs,val,is_p_reciprocity,bmod_ref,rho_ref);
	if (source_ghost_enabled)
	{
		// TMJ, Vz require a POSITIVE ghost, which I make happen by flipping the sign of ampl3 for this call.
		_cuApply_Trilinear_Source_To_VxVyVz(em,Density_min,Density_range,cmp,x0,y0,nx,ny,nz,dti,is_force,ampl1,ampl2,-ampl3,xs,ys,2.0f*ghost_sea_surface-zs,-val,is_p_reciprocity,bmod_ref,rho_ref);
	}
}

__global__
#if __CUDA_ARCH__ >= 370
__launch_bounds__(256,8)
#elif __CUDA_ARCH__ >= 300
__launch_bounds__(256,5)
#elif __CUDA_ARCH__ >= 200
__launch_bounds__(256,3)
#endif
void
cuPropagate_Particle_Velocities_Kernel(
        int timestep,
        int x0,                 // x coordinate of westernmost coordinate in block
        int y0,                 // y coordinate of southernmost coordinate in block
        int y1,
        int m1_y0,
        int m1_y1,
        int vol_nx,             // dimensions of global volume
        int vol_ny,
        int vol_nz,
        float dti,
        unsigned int* em,       // earth model, 4 interleaved integers. y(0)
        float* cmp,             // Vx, Vy, Vz, Sx, Sy and Sz, middle, t(1), y(0)
        float* m1L,             // txx, tyy, tzz, txy, txz and tyz in that order. left halo, t(0), y(0)
        float* m1C,             // ..middle, t(0), y(0)
        float* m1R,             // ..right halo, t(0), y(0)
        float* m2C,             // Vx, Vy, Vz, Sx, Sy and Sz in that order. middle, t(-1), y(0)
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
        __shared__ float buf[768];      // NON-persistent buffer

        __shared__ float tzzbuf[384];   // persistent buffers
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

                // compute dytyy
                buf[threadIdx.x+threadIdx.y*32+128] = tmp2;  // deposit middle section
                buf[threadIdx.x+threadIdx.y*32+64*(threadIdx.y&4)] = tmp1;
                __syncthreads();
                float dytyy = ( C0 * (buf[cuCompTYYIdx(1)] - buf[cuCompTYYIdx( 0)]) +
                                C1 * (buf[cuCompTYYIdx(2)] - buf[cuCompTYYIdx(-1)]) +
                                C2 * (buf[cuCompTYYIdx(3)] - buf[cuCompTYYIdx(-2)]) +
                                C3 * (buf[cuCompTYYIdx(4)] - buf[cuCompTYYIdx(-3)]) ) * inv_DY;

                // compute dztzz
                // ..load 8 next z and transpose to XYZ
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
                dztzz = cuTransposeXZY2XYZ(buf,dztzz);  // this actually transposes back from XYZ to XZY.

                // compute dxtxy
                buf[threadIdx.x+threadIdx.y*96] = tmp4;
                buf[threadIdx.x+threadIdx.y*96+32] = txy_p0;
                buf[threadIdx.x+threadIdx.y*96+64] = txy_p4;
                __syncthreads();
                float dxtxy = ( C0 * (buf[cuCompTXXIdx(1)] - txy_p0               ) +
                                C1 * (buf[cuCompTXXIdx(2)] - buf[cuCompTXXIdx(-1)]) +
                                C2 * (buf[cuCompTXXIdx(3)] - buf[cuCompTXXIdx(-2)]) +
                                C3 * (txy_p4               - buf[cuCompTXXIdx(-3)]) ) * inv_DX;

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

                // compute dxtxz
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

                // compute dytyz
                float tyz_p0 = cuTransposeXZY2XYZ(buf, tyzbuf[threadIdx.x+threadIdx.y*32+128] );  // read middle section from persistent tyz buffer
                buf[threadIdx.x+threadIdx.y*32+128] = tyz_p0;
                buf[threadIdx.x+threadIdx.y*32+64*(threadIdx.y&4)] = tmp10;
                __syncthreads();
                float dytyz = ( C0 * (tyz_p0               - buf[cuCompTYYIdx(-1)]) +
                                C1 * (buf[cuCompTYYIdx(1)] - buf[cuCompTYYIdx(-2)]) +
                                C2 * (buf[cuCompTYYIdx(2)] - buf[cuCompTYYIdx(-3)]) +
                                C3 * (buf[cuCompTYYIdx(3)] - buf[cuCompTYYIdx(-4)]) ) * inv_DY;

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

                if (y <= y1)
                {
                        // get word3 from earth model
                        float Q, Density;
                        cuUnpack_Q_Density(em_word3,Q_min,Q_range,&Q,Density_min,Density_range,&Density);
                        Q = 1.0f / Q;  // compressed model actually stores inverse of Q.

                        float deta = Compute_ABC(x,y,z,vol_nx,vol_ny,vol_nz,nabc_top,nabc_bot,nabc_sdx,nabc_sdy,vpvert_avtop,vpvert_avbot,inv_DX,inv_DY,inv_DZ);
                        float dabc = (1.0f - 0.5f*deta*dti) / (1.0f + 0.5f*deta*dti);

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

                        cmp[offset+3*one_wf_size_f] = sx;
                        cmp[offset+4*one_wf_size_f] = sy;
                        cmp[offset+5*one_wf_size_f] = sz;

                        // Update viscoelastic particle velocities:
                        float old_vx = m2C[offset];
                        float old_vy = m2C[offset+one_wf_size_f];
                        float old_vz = m2C[offset+2*one_wf_size_f];

                        float factor = dti / Density;

                        float vx = factor * ( sx + dxtxx + dytxy + dztxz );
                        float vy = factor * ( sy + dxtxy + dytyy + dztyz );
                        float vz = factor * ( sz + dxtxz + dytyz + dztzz );
                        vx = vx + dabc * old_vx;
                        vy = vy + dabc * old_vy;
                        vz = vz + dabc * old_vz;

                        cmp[offset] = vx;
                        cmp[offset+one_wf_size_f] = vy;
                        cmp[offset+2*one_wf_size_f] = vz;
                }

                // increase offsets
                offset += 32;
        }
}

__global__ 
#if __CUDA_ARCH__ >= 370
__launch_bounds__(512,4)
#elif __CUDA_ARCH__ >= 300
__launch_bounds__(512,2)
#elif __CUDA_ARCH__ >= 200
__launch_bounds__(512,1)
#endif
void 
cuPropagate_Particle_Velocities_Kernel_16(
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
	float C4,
	float C5,
	float C6,
	float C7,
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
	__shared__ float buf[1536];	// NON-persistent buffer

	__shared__ float tzzbuf[1024];	// persistent buffers
	__shared__ float txzbuf[1024];   // some values are transferred from one iZ to the next
	__shared__ float tyzbuf[1024];

	const int nx = 8;

	int z_per_block = (((vol_nz/8) + gridDim.z - 1) / gridDim.z) * 8;
	int z0 = z_per_block * blockIdx.z;
	int z1 = z0 + z_per_block - 1;
	if (z1 >= vol_nz) z1 = vol_nz - 1;
	int nz = z1 - z0 + 1;
	if (nz <= 0) return;

	int offset = (threadIdx.y + blockIdx.y * 8) * one_y_size_f + threadIdx.x + z0 * 8;

	// populate persistent buffers
	int y = y0 + (threadIdx.y + blockIdx.y * 8);
	if (z0 == 0)
	{
		tzzbuf[threadIdx.x+threadIdx.y*64+512] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset+tzz_off] : 0.0f);
		txzbuf[threadIdx.x+threadIdx.y*64+512] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset+txz_off] : 0.0f);
		tyzbuf[threadIdx.x+threadIdx.y*64+512] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset+tyz_off] : 0.0f);
		tzzbuf[threadIdx.x+(7-threadIdx.y)*64] = threadIdx.y == 7 ? 0.0f : -tzzbuf[threadIdx.x+(9+threadIdx.y)*64];
		txzbuf[threadIdx.x+(7-threadIdx.y)*64] = threadIdx.y == 7 ? 0.0f : -txzbuf[threadIdx.x+(9+threadIdx.y)*64];
		tyzbuf[threadIdx.x+(7-threadIdx.y)*64] = -tyzbuf[threadIdx.x+(8+threadIdx.y)*64];
		if (threadIdx.y == 0) tzzbuf[threadIdx.x+512] = 0.0f;
	}
	else
	{
		tzzbuf[threadIdx.x+threadIdx.y*64] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset-64+tzz_off] : 0.0f);
		txzbuf[threadIdx.x+threadIdx.y*64] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset-64+txz_off] : 0.0f);
		tyzbuf[threadIdx.x+threadIdx.y*64] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset-64+tyz_off] : 0.0f);
		tzzbuf[threadIdx.x+threadIdx.y*64+512] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset+tzz_off] : 0.0f);
		txzbuf[threadIdx.x+threadIdx.y*64+512] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset+txz_off] : 0.0f);
		tyzbuf[threadIdx.x+threadIdx.y*64+512] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset+tyz_off] : 0.0f);
	}

	for (int iZ = 0;  iZ < nz/8;  ++iZ)
	{
		int x = x0 + (threadIdx.x & 7);
		int z = z0 + iZ * 8 + (threadIdx.x / 8);

		float tmp3, tmp7, tmp8;
		if (z < vol_nz-8 && y <= m1_y1)
		{
			tmp3 = m1C[offset+tzz_off+64];
			tmp7 = m1C[offset+txz_off+64];
			tmp8 = m1C[offset+tyz_off+64];
		}
		else
		{
			tmp3 = tmp7 = tmp8 = 0.0f;
		}

		float tmp4, tmp5, txx_m8;
		if (m1L != 0L && y <= m1_y1)
		{
			tmp4 = m1L[offset+txy_off];
			tmp5 = m1L[offset+txz_off];
			txx_m8 = m1L[offset];
		}
		else
		{
			tmp4 = tmp5 = txx_m8 = 0.0f;
		}

		float tmp6, txy_p8;
		if (m1R != 0L && y <= m1_y1)
		{
			tmp6 = m1R[offset+txz_off];
			txy_p8 = m1R[offset+txy_off];
		}
		else
		{
			tmp6 = txy_p8 = 0.0f;
		}

		unsigned int em_word3 = (y <= y1) ? em[(threadIdx.y + blockIdx.y*8) * em_one_y_size_f + (iZ*64) + (z0*nx) + threadIdx.x + 3*em_one_word_size_f] : 0;

		float txx_p0 = y <= m1_y1 ? m1C[offset] : 0.0f;
                float txy_p0 = y <= m1_y1 ? m1C[offset+txy_off] : 0.0f;

		float tmp2 = y <= m1_y1 ? m1C[offset+tyy_off] : 0.0f;

		float tmp1a, tmp1b, tmp9a, tmp9b, tmp10a, tmp10b;
		if (y-8 >= m1_y0)
		{
			tmp1a = m1C[offset+tyy_off-8*one_y_size_f];
			tmp9a = m1C[offset+txy_off-8*one_y_size_f];
			tmp10a = m1C[offset+tyz_off-8*one_y_size_f];
		}
		else
		{
			tmp1a = tmp9a = tmp10a = 0.0f;
		}
		if (y+8 <= m1_y1)
		{
			tmp1b = m1C[offset+tyy_off+8*one_y_size_f];
			tmp9b = m1C[offset+txy_off+8*one_y_size_f];
			tmp10b = m1C[offset+tyz_off+8*one_y_size_f];
		}
		else
		{
			tmp1b = tmp9b = tmp10b = 0.0f;
		}

		// compute dxtxx
		buf[threadIdx.x+threadIdx.y*192] = txx_m8;
		buf[threadIdx.x+threadIdx.y*192+64] = txx_p0;
		buf[threadIdx.x+threadIdx.y*192+128] = m1R != 0L && y <= m1_y1 ? m1R[offset] : 0.0f;
		__syncthreads();
		float dxtxx = ( C0 * (txx_p0                  - buf[cuCompTXXIdx_16(-1)]) + 
				C1 * (buf[cuCompTXXIdx_16(1)] - buf[cuCompTXXIdx_16(-2)]) + 
				C2 * (buf[cuCompTXXIdx_16(2)] - buf[cuCompTXXIdx_16(-3)]) +
				C3 * (buf[cuCompTXXIdx_16(3)] - buf[cuCompTXXIdx_16(-4)]) + 
			 	C4 * (buf[cuCompTXXIdx_16(4)] - buf[cuCompTXXIdx_16(-5)]) +
			 	C5 * (buf[cuCompTXXIdx_16(5)] - buf[cuCompTXXIdx_16(-6)]) +
			 	C6 * (buf[cuCompTXXIdx_16(6)] - buf[cuCompTXXIdx_16(-7)]) +
			 	C7 * (buf[cuCompTXXIdx_16(7)] - txx_m8                  ) ) * inv_DX;
		__syncthreads();  // wait for computes before reusing buf

		// compute dytyy
		buf[threadIdx.x+threadIdx.y*64] = tmp1a;
		buf[threadIdx.x+threadIdx.y*64+512] = tmp2;
		buf[threadIdx.x+threadIdx.y*64+1024] = tmp1b;
		__syncthreads();
		float dytyy = ( C0 * (buf[cuCompTYYIdx_16(1)] - buf[cuCompTYYIdx_16( 0)]) + 
				C1 * (buf[cuCompTYYIdx_16(2)] - buf[cuCompTYYIdx_16(-1)]) +
				C2 * (buf[cuCompTYYIdx_16(3)] - buf[cuCompTYYIdx_16(-2)]) +
				C3 * (buf[cuCompTYYIdx_16(4)] - buf[cuCompTYYIdx_16(-3)]) +
				C4 * (buf[cuCompTYYIdx_16(5)] - buf[cuCompTYYIdx_16(-4)]) +
				C5 * (buf[cuCompTYYIdx_16(6)] - buf[cuCompTYYIdx_16(-5)]) +
				C6 * (buf[cuCompTYYIdx_16(7)] - buf[cuCompTYYIdx_16(-6)]) +
				C7 * (buf[cuCompTYYIdx_16(8)] - buf[cuCompTYYIdx_16(-7)]) ) * inv_DY;

		// compute dztzz
		// ..load 8 next z and transpose to XYZ
		//float tmp3 = (iZ < ((nz/8)-1)) ? m1C[offset+tzz_off+32] : 0.0f;
		float v2 = cuTransposeXZY2XYZ_16(buf, tmp3);
		buf[threadIdx.x+threadIdx.y*64] = tzzbuf[threadIdx.x+threadIdx.y*64];
		buf[threadIdx.x+threadIdx.y*64+512] = tzzbuf[threadIdx.x+threadIdx.y*64+512];
		buf[threadIdx.x+threadIdx.y*64+1024] = v2;
		// ..store next 16*z in tzzbuf
		tzzbuf[threadIdx.x+threadIdx.y*64] = tzzbuf[threadIdx.x+threadIdx.y*64+512];
		tzzbuf[threadIdx.x+threadIdx.y*64+512] = v2;
		__syncthreads(); // needed to prevent race condition
		// note that we can use cuCompTYYIdx in place of cuCompTZZIdx after the transpose
		float dztzz = -( C0 * (buf[cuCompTYYIdx_16(1)] - buf[cuCompTYYIdx_16( 0)]) + 
				C1 * (buf[cuCompTYYIdx_16(2)] - buf[cuCompTYYIdx_16(-1)]) +
				C2 * (buf[cuCompTYYIdx_16(3)] - buf[cuCompTYYIdx_16(-2)]) +
				C3 * (buf[cuCompTYYIdx_16(4)] - buf[cuCompTYYIdx_16(-3)]) +
				C4 * (buf[cuCompTYYIdx_16(5)] - buf[cuCompTYYIdx_16(-4)]) +
				C5 * (buf[cuCompTYYIdx_16(6)] - buf[cuCompTYYIdx_16(-5)]) +
				C6 * (buf[cuCompTYYIdx_16(7)] - buf[cuCompTYYIdx_16(-6)]) +
				C7 * (buf[cuCompTYYIdx_16(8)] - buf[cuCompTYYIdx_16(-7)]) ) * inv_DZ;
		//float tmp4 = m1L != 0L ? m1L[offset+txy_off] : 0.0f;
		dztzz = cuTransposeXZY2XYZ_16(buf,dztzz);  // this actually transposes back from XYZ to XZY.

		// compute dxtxy
                buf[threadIdx.x+threadIdx.y*192] = tmp4;
                buf[threadIdx.x+threadIdx.y*192+64] = txy_p0;
		buf[threadIdx.x+threadIdx.y*192+128] = txy_p8;
                __syncthreads();
                float dxtxy = ( C0 * (buf[cuCompTXXIdx_16(1)] - txy_p0               ) +
                                C1 * (buf[cuCompTXXIdx_16(2)] - buf[cuCompTXXIdx_16(-1)]) +
                                C2 * (buf[cuCompTXXIdx_16(3)] - buf[cuCompTXXIdx_16(-2)]) +
                                C3 * (buf[cuCompTXXIdx_16(4)] - buf[cuCompTXXIdx_16(-3)]) +
				C4 * (buf[cuCompTXXIdx_16(5)] - buf[cuCompTXXIdx_16(-4)]) +
				C5 * (buf[cuCompTXXIdx_16(6)] - buf[cuCompTXXIdx_16(-5)]) +
				C6 * (buf[cuCompTXXIdx_16(7)] - buf[cuCompTXXIdx_16(-6)]) +
				C7 * (txy_p8                  - buf[cuCompTXXIdx_16(-7)]) ) * inv_DX;

		// ..compute dytxy
		float v4 = buf[threadIdx.x+threadIdx.y*192+64];  // read middle section for dytxy from shared memory
		__syncthreads();
		buf[threadIdx.x+threadIdx.y*64] = tmp9a;
		buf[threadIdx.x+threadIdx.y*64+512] = v4;
		buf[threadIdx.x+threadIdx.y*64+1024] = tmp9b;
		__syncthreads();
		float dytxy = ( C0 * (buf[cuCompTYYIdx_16(0)] - buf[cuCompTYYIdx_16(-1)]) +
                                C1 * (buf[cuCompTYYIdx_16(1)] - buf[cuCompTYYIdx_16(-2)]) +
                                C2 * (buf[cuCompTYYIdx_16(2)] - buf[cuCompTYYIdx_16(-3)]) +
                                C3 * (buf[cuCompTYYIdx_16(3)] - buf[cuCompTYYIdx_16(-4)]) +
                                C4 * (buf[cuCompTYYIdx_16(4)] - buf[cuCompTYYIdx_16(-5)]) +
                                C5 * (buf[cuCompTYYIdx_16(5)] - buf[cuCompTYYIdx_16(-6)]) +
                                C6 * (buf[cuCompTYYIdx_16(6)] - buf[cuCompTYYIdx_16(-7)]) +
                                C7 * (buf[cuCompTYYIdx_16(7)] - buf[cuCompTYYIdx_16(-8)]) ) * inv_DY;

		// compute dxtxz
		//float tmp5 = m1L != 0L ? m1L[offset+txz_off] : 0.0f;
		//float tmp6 = m1R != 0L ? m1R[offset+txz_off] : 0.0f;
                float txz_p0 = cuTransposeXZY2XYZ_16(buf, txzbuf[threadIdx.x+threadIdx.y*64+512]);  // read middle section from persistent txz buffer
		buf[threadIdx.x+threadIdx.y*192] = tmp5;
                buf[threadIdx.x+threadIdx.y*192+64] = txz_p0;
                float txz_p8 = tmp6;
		buf[threadIdx.x+threadIdx.y*192+128] = txz_p8;
                __syncthreads();
                float dxtxz = ( C0 * (buf[cuCompTXXIdx_16(1)] - txz_p0               ) +
                                C1 * (buf[cuCompTXXIdx_16(2)] - buf[cuCompTXXIdx_16(-1)]) +
                                C2 * (buf[cuCompTXXIdx_16(3)] - buf[cuCompTXXIdx_16(-2)]) +
                                C3 * (buf[cuCompTXXIdx_16(4)] - buf[cuCompTXXIdx_16(-3)]) +
				C4 * (buf[cuCompTXXIdx_16(5)] - buf[cuCompTXXIdx_16(-4)]) +
				C5 * (buf[cuCompTXXIdx_16(6)] - buf[cuCompTXXIdx_16(-5)]) +
				C6 * (buf[cuCompTXXIdx_16(7)] - buf[cuCompTXXIdx_16(-6)]) +
				C7 * (txz_p8                  - buf[cuCompTXXIdx_16(-7)]) ) * inv_DX;

		// ..compute dztxz
		//float tmp7 = (iZ < ((nz/8)-1)) ? m1C[offset+txz_off+32] : 0.0f;
		float v5 = cuTransposeXZY2XYZ_16(buf, tmp7);  // read next 8 z from gmem
		buf[threadIdx.x+threadIdx.y*64] = txzbuf[threadIdx.x+threadIdx.y*64];
		buf[threadIdx.x+threadIdx.y*64+512] = txzbuf[threadIdx.x+threadIdx.y*64+512];
		buf[threadIdx.x+threadIdx.y*64+1024] = v5;
                // ..store next 16 z in txzbuf
		txzbuf[threadIdx.x+threadIdx.y*64] = txzbuf[threadIdx.x+threadIdx.y*64+512];
		txzbuf[threadIdx.x+threadIdx.y*64+512] = v5;
		__syncthreads();
                // note that we can use cuCompTYYIdx in place of cuCompTZZIdx after the transpose
                float dztxz = -( C0 * (buf[cuCompTYYIdx_16(0)] - buf[cuCompTYYIdx_16(-1)]) +
                                C1 * (buf[cuCompTYYIdx_16(1)] - buf[cuCompTYYIdx_16(-2)]) +
                                C2 * (buf[cuCompTYYIdx_16(2)] - buf[cuCompTYYIdx_16(-3)]) +
                                C3 * (buf[cuCompTYYIdx_16(3)] - buf[cuCompTYYIdx_16(-4)]) + 
                                C4 * (buf[cuCompTYYIdx_16(4)] - buf[cuCompTYYIdx_16(-5)]) + 
                                C5 * (buf[cuCompTYYIdx_16(5)] - buf[cuCompTYYIdx_16(-6)]) + 
                                C6 * (buf[cuCompTYYIdx_16(6)] - buf[cuCompTYYIdx_16(-7)]) + 
                                C7 * (buf[cuCompTYYIdx_16(7)] - buf[cuCompTYYIdx_16(-8)]) ) * inv_DZ;
                dztxz = cuTransposeXZY2XYZ_16(buf,dztxz);  // this actually transposes back from XYZ to XZY.

		// compute dytyz
		float tyz_p0 = cuTransposeXZY2XYZ_16(buf, tyzbuf[threadIdx.x+threadIdx.y*64+512] );  // read middle section from persistent tyz buffer
		buf[threadIdx.x+threadIdx.y*64] = tmp10a;
		buf[threadIdx.x+threadIdx.y*64+512] = tyz_p0;
		buf[threadIdx.x+threadIdx.y*64+1024] = tmp10b;
		__syncthreads();
		float dytyz = ( C0 * (tyz_p0                  - buf[cuCompTYYIdx_16(-1)]) + 
				C1 * (buf[cuCompTYYIdx_16(1)] - buf[cuCompTYYIdx_16(-2)]) +
				C2 * (buf[cuCompTYYIdx_16(2)] - buf[cuCompTYYIdx_16(-3)]) +
				C3 * (buf[cuCompTYYIdx_16(3)] - buf[cuCompTYYIdx_16(-4)]) +
				C4 * (buf[cuCompTYYIdx_16(4)] - buf[cuCompTYYIdx_16(-5)]) +
				C5 * (buf[cuCompTYYIdx_16(5)] - buf[cuCompTYYIdx_16(-6)]) +
				C6 * (buf[cuCompTYYIdx_16(6)] - buf[cuCompTYYIdx_16(-7)]) +
				C7 * (buf[cuCompTYYIdx_16(7)] - buf[cuCompTYYIdx_16(-8)]) ) * inv_DY;

		// ..compute dztyz
		//float tmp8 = (iZ < ((nz/8)-1)) ? m1C[offset+tyz_off+32] : 0.0f;
		float v8 = cuTransposeXZY2XYZ_16(buf, tmp8);  // read next 8 z from gmem
		buf[threadIdx.x+threadIdx.y*64] = tyzbuf[threadIdx.x+threadIdx.y*64];
		buf[threadIdx.x+threadIdx.y*64+512] = tyzbuf[threadIdx.x+threadIdx.y*64+512];
		buf[threadIdx.x+threadIdx.y*64+1024] = v8;
                // ..store next 16 z in txzbuf
		tyzbuf[threadIdx.x+threadIdx.y*64] = tyzbuf[threadIdx.x+threadIdx.y*64+512];
		tyzbuf[threadIdx.x+threadIdx.y*64+512] = v8;
                __syncthreads();
                // note that we can use cuCompTYYIdx in place of cuCompTZZIdx after the transpose
                float dztyz = -( C0 * (buf[cuCompTYYIdx_16(0)] - buf[cuCompTYYIdx_16(-1)]) +
                                C1 * (buf[cuCompTYYIdx_16(1)] - buf[cuCompTYYIdx_16(-2)]) +
                                C2 * (buf[cuCompTYYIdx_16(2)] - buf[cuCompTYYIdx_16(-3)]) +
                                C3 * (buf[cuCompTYYIdx_16(3)] - buf[cuCompTYYIdx_16(-4)]) +
                                C4 * (buf[cuCompTYYIdx_16(4)] - buf[cuCompTYYIdx_16(-5)]) +
                                C5 * (buf[cuCompTYYIdx_16(5)] - buf[cuCompTYYIdx_16(-6)]) +
                                C6 * (buf[cuCompTYYIdx_16(6)] - buf[cuCompTYYIdx_16(-7)]) +
                                C7 * (buf[cuCompTYYIdx_16(7)] - buf[cuCompTYYIdx_16(-8)]) ) * inv_DZ;
                dztyz = cuTransposeXZY2XYZ_16(buf,dztyz);  // this actually transposes back from XYZ to XZY.

		if (y <= y1)
		{
			// get word3 from earth model
			float Q, Density;
			cuUnpack_Q_Density(em_word3,Q_min,Q_range,&Q,Density_min,Density_range,&Density);
			Q = 1.0f / Q;  // compressed model actually stores inverse of Q.

			float deta = Compute_ABC(x,y,z,vol_nx,vol_ny,vol_nz,nabc_top,nabc_bot,nabc_sdx,nabc_sdy,vpvert_avtop,vpvert_avbot,inv_DX,inv_DY,inv_DZ);
			float dabc = (1.0f - 0.5f*deta*dti) / (1.0f + 0.5f*deta*dti);

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

			cmp[offset+3*one_wf_size_f] = sx;
			cmp[offset+4*one_wf_size_f] = sy;
			cmp[offset+5*one_wf_size_f] = sz;

			// Update viscoelastic particle velocities:
			float old_vx = m2C[offset];
			float old_vy = m2C[offset+one_wf_size_f];
			float old_vz = m2C[offset+2*one_wf_size_f];

			float factor = dti / Density;

			float vx = factor * ( sx + dxtxx + dytxy + dztxz );
			float vy = factor * ( sy + dxtxy + dytyy + dztyz );
			float vz = factor * ( sz + dxtxz + dytyz + dztzz );
			vx = vx + dabc * old_vx;
			vy = vy + dabc * old_vy;
			vz = vz + dabc * old_vz;

			cmp[offset] = vx;
			cmp[offset+one_wf_size_f] = vy;
			cmp[offset+2*one_wf_size_f] = vz;
		}

		// increase offsets
		offset += 64;
	}
}

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
	int spatial_order,
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
	bool is_p_reciprocity,	// only valid for common receiver gather: pressure receivers, a.k.a. reciprocal sources
	float bmod_ref,
	float rho_ref
	)
{
	assert(spatial_order == 8 || spatial_order == 16);
	if (!source_ghost_enabled) ghost_sea_surface_z = 0;

	int one_wf_size = one_y_size / 6;
	int em_one_word_size = one_wf_size;
	int em_one_y_size = em_one_word_size * 4;

	const int tyy_off = one_wf_size / 4;
	const int tzz_off = 2 * tyy_off;
	const int txy_off = 3 * tyy_off;
	const int txz_off = 4 * tyy_off;
	const int tyz_off = 5 * tyy_off;

	int nx = spatial_order / 2;
	int ny = y1 - y0 + 1;
	int nz = vol_nz;

	if (spatial_order == 8)
	{
		dim3 blockShape(32,8,1);
		dim3 gridShape(1,(ny+7)/8,num_z);
		cuPropagate_Particle_Velocities_Kernel<<<gridShape,blockShape,0,stream>>>(
				timestep,
				x0,y0,y1,m1_y0,m1_y1,vol_nx,vol_ny,vol_nz,dti,
				(unsigned int*)em,(float*)cmp,(float*)m1L,(float*)m1C,(float*)m1R,(float*)m2C,
				C8_0,C8_1,C8_2,C8_3,
				inv_DX,inv_DY,inv_DZ,
				vpvert_avtop,vpvert_avbot,nabc_sdx,nabc_sdy,nabc_top,nabc_bot,Q_min,Q_range/255.0f,fq,Density_min,Density_range/255.0f,
				one_wf_size/4,one_y_size/4,em_one_word_size/4,em_one_y_size/4,
				tyy_off,tzz_off,txy_off,txz_off,tyz_off);
	}
	else if (spatial_order == 16)
	{
		dim3 blockShape(64,8,1);
		dim3 gridShape(1,(ny+7)/8,num_z);
		cuPropagate_Particle_Velocities_Kernel_16<<<gridShape,blockShape,0,stream>>>(
				timestep,
				x0,y0,y1,m1_y0,m1_y1,vol_nx,vol_ny,vol_nz,dti,
				(unsigned int*)em,(float*)cmp,(float*)m1L,(float*)m1C,(float*)m1R,(float*)m2C,
				C16_0,C16_1,C16_2,C16_3,C16_4,C16_5,C16_6,C16_7,
				inv_DX,inv_DY,inv_DZ,
				vpvert_avtop,vpvert_avbot,nabc_sdx,nabc_sdy,nabc_top,nabc_bot,Q_min,Q_range/255.0f,fq,Density_min,Density_range/255.0f,
				one_wf_size/4,one_y_size/4,em_one_word_size/4,em_one_y_size/4,
				tyy_off,tzz_off,txy_off,txz_off,tyz_off);
	}
#ifdef GPU_DEBUG
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
#endif

	//
	// add source term(s)
	//
	if (inject_source && (is_force || is_velocity))
	{
		float ghost_sea_surface = (float)ghost_sea_surface_z;
		if (source_interpolation_method == Point)
		{
			dim3 blockShape3(1,1,1);
			dim3 gridShape3(1,1,1);
			cuApply_Point_Source_To_VxVyVz<<<gridShape3,blockShape3,0,stream>>>(
					(int*)em,Density_min,Density_range/255.0f,(float*)cmp,x0,y0,nx,ny,nz,dti,is_force,
					ampl1,ampl2,ampl3,xsou,ysou,zsou,svaw_sample,
					source_ghost_enabled,ghost_sea_surface,
					is_p_reciprocity,bmod_ref,rho_ref);
		}
		else if(source_interpolation_method == Trilinear)
		{
			dim3 blockShape3(1,1,1);
			dim3 gridShape3(1,1,1);
			cuApply_Trilinear_Source_To_VxVyVz<<<gridShape3,blockShape3,0,stream>>>(
					(int*)em,Density_min,Density_range/255.0f,(float*)cmp,x0,y0,nx,ny,nz,dti,is_force,
					ampl1,ampl2,ampl3,xsou,ysou,zsou,svaw_sample,
					source_ghost_enabled,ghost_sea_surface,
					is_p_reciprocity,bmod_ref,rho_ref);
		}
		else if(source_interpolation_method == Sinc)
		{
			// use only one thread along z to prevent possible race condition
			dim3 blockShape2(8,8,1);
			dim3 gridShape2(1,1,1);
			cuApply_Source_Term_To_VxVyVz<<<gridShape2,blockShape2,0,stream>>>(
					em,Density_min,Density_range/255.0f,cmp,x0,y0,0,nx,ny,nz,dti,is_force,
					ampl1,ampl2,ampl3,xsou,ysou,zsou,svaw_sample,
					source_ghost_enabled,ghost_sea_surface_z,
					is_p_reciprocity,bmod_ref,rho_ref);
		}
	}
#ifdef GPU_DEBUG
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
#endif
}

