//
// This kernel propagates strain rates one timestep.
// Compile with L1 caching turned off (-Xptxas -dlcm=cg)
//

__device__ 
void __cuApply_Source_Term_To_TxxTyyTzz(
	int thr_z,
	int timestep,
	float* cmp,
	int x0,
        int y0,
        int z0,
        int nx,
        int ny,
        int nz,
	float dti,
	float ampl1,
        float xs,
        float ys,
        float zs,
        float val,
	int icell,
	int jcell,
	int kcell,
	unsigned int* em,
	float Vp_min,
	float Vp_range,
	float Vs_min,
	float Vs_range,
	float Density_min,
	float Density_range,
	float bmod_ref,
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
			(my_z >= ghost_sea_surface_z && my_z < nz)  // TMJ 04/08/15 Turned off mirroring, it's not in Kurt Nihei's original code
	   )
	{
		// ..fractional distance from grid pt to sou:
		float dx_frac = xs - (float)icell;
		float dy_frac = ys - (float)jcell;
                float dz_frac = zs - (float)kcell;

                float fsinc = cuGen_Sinc_Weight(threadIdx.x,threadIdx.y,thr_z,dx_frac,dy_frac,dz_frac);
                if (fsinc != 0.0f)
                {
			// mirror source if necessary
			//my_z = my_z < 0 ? -my_z : my_z;
		
			// TMJ 05/06/14
			// Mirroring introduces a potential race condition, two threads will update same cell
			// The way to solve this is to let same thread handle all the z indexes

                        int one_wf_size_f = nx * nz;
                        int one_y_size_f = one_wf_size_f * 6;

			int em_one_word_size_f = one_wf_size_f;
			int em_one_y_size_f = em_one_word_size_f * 4;

			int emIdx = my_x + my_y * em_one_y_size_f + my_z * nx;
			unsigned int em_word0 = em[emIdx];
			unsigned int em_word3 = em[emIdx + 3*em_one_word_size_f];
			float rho, bmod;
			cuUnpack_And_Compute_Bulk_Modulus(em_word0, em_word3, Vp_min, Vp_range, Vs_min, Vs_range, Density_min, Density_range, &rho, &bmod);
			float scale_sou = bmod / bmod_ref;

                        int idx = my_x + my_y * one_y_size_f + my_z * nx;
			cmp[idx                ] = cmp[idx                ] - fsinc * scale_sou * dti * val * ampl1;
			cmp[idx+  one_wf_size_f] = cmp[idx+  one_wf_size_f] - fsinc * scale_sou * dti * val * ampl1;
			cmp[idx+2*one_wf_size_f] = cmp[idx+2*one_wf_size_f] - fsinc * scale_sou * dti * val * ampl1;

			if (source_ghost_enabled)
			{
				int my_z_ghost = 2 * ghost_sea_surface_z - my_z;
				int emIdx_ghost = my_x + my_y * em_one_y_size_f + my_z_ghost * nx;
				int idx_ghost = my_x + my_y * one_y_size_f + my_z_ghost * nx;

				unsigned int em_word0_ghost = em[emIdx_ghost];
				unsigned int em_word3_ghost = em[emIdx_ghost + 3*em_one_word_size_f];
				float rho_ghost, bmod_ghost;
				cuUnpack_And_Compute_Bulk_Modulus(em_word0_ghost, em_word3_ghost, Vp_min, Vp_range, Vs_min, Vs_range, Density_min, Density_range, &rho_ghost, &bmod_ghost);
				float scale_sou_ghost = bmod_ghost / bmod_ref;

				cmp[idx_ghost                ] = cmp[idx_ghost                ] + fsinc * scale_sou_ghost * dti * val * ampl1;
				cmp[idx_ghost+  one_wf_size_f] = cmp[idx_ghost+  one_wf_size_f] + fsinc * scale_sou_ghost * dti * val * ampl1;
				cmp[idx_ghost+2*one_wf_size_f] = cmp[idx_ghost+2*one_wf_size_f] + fsinc * scale_sou_ghost * dti * val * ampl1;
			}
                }
        }
}

__device__
void _cuApply_Source_Term_To_TxxTyyTzz(
	int timestep,
        float* cmp,
        int x0,
        int y0,
        int z0,
        int nx,
        int ny,
        int nz,
        float dti,
        float ampl1,
        float xs,
        float ys,
        float zs,
        float val,
	unsigned int* em,
	float Vp_min,
	float Vp_range,
	float Vs_min,
	float Vs_range,
	float Density_min,
	float Density_range,
	float bmod_ref,
	bool source_ghost_enabled,
	int ghost_sea_surface_z
        )
{
        // nearest grid point:
        int icell = (int)lrintf(xs);
        int jcell = (int)lrintf(ys);
        int kcell = (int)lrintf(zs);

	for (int thr_z = 0;  thr_z < 8;  ++thr_z)
	{
        	__cuApply_Source_Term_To_TxxTyyTzz(thr_z,timestep,cmp,x0,y0,z0,nx,ny,nz,dti,ampl1,xs,ys,zs,val,icell,jcell,kcell,em,Vp_min,Vp_range,Vs_min,Vs_range,Density_min,Density_range,bmod_ref,source_ghost_enabled,ghost_sea_surface_z);
	}
}

__global__
void cuApply_Source_Term_To_TxxTyyTzz(
	int timestep,
        float* cmp,
        int x0,
        int y0,
        int z0,
        int nx,
        int ny,
        int nz,
        float dti,
        float ampl1,
        float xs,
        float ys,
        float zs,
        float val,
	bool source_ghost_enabled,
	int ghost_sea_surface_z,
	unsigned int* em,
	float Vp_min,
	float Vp_range,
	float Vs_min,
	float Vs_range,
	float Density_min,
	float Density_range,
	float bmod_ref
        )
{
	_cuApply_Source_Term_To_TxxTyyTzz(timestep,cmp,x0,y0,z0,nx,ny,nz,dti,ampl1,xs,ys,zs,val,em,Vp_min,Vp_range,Vs_min,Vs_range,Density_min,Density_range,bmod_ref,source_ghost_enabled,ghost_sea_surface_z);
}

__device__ 
void cuApply_Point_Source_To_Single_TxxTyyTzz(
	float* cmp,
        int nx,
        int ny,
        int nz,
        int ix,
	int iy,
	int iz,
	float delta,
	int one_wf_size_f,
	int one_y_size_f,
	unsigned int* em,
	float Vp_min,
	float Vp_range,
	float Vs_min,
	float Vs_range,
	float Density_min,
	float Density_range,
	float bmod_ref
        )
{
	if (delta != 0.0f && ix >= 0 && ix < nx && iy >= 0 && iy < ny && iz >= 0 && iz < nz)
	{
		int idx = ix + iy * one_y_size_f + iz * nx;

		int em_one_word_size_f = one_wf_size_f;
		int em_one_y_size = em_one_word_size_f * 4;
		int emIdx = ix + iy * em_one_y_size + iz * nx;
		unsigned int em_word0 = em[emIdx];
		unsigned int em_word3 = em[emIdx + 3*em_one_word_size_f];
		float rho, bmod;
		cuUnpack_And_Compute_Bulk_Modulus(em_word0, em_word3, Vp_min, Vp_range, Vs_min, Vs_range, Density_min, Density_range, &rho, &bmod);
		float scale_sou = bmod / bmod_ref;

		//printf("\n*****\nmy_x=%d, my_y=%d, my_z=%d :: rho = %e, bmod = %e, scale_sou = %e\n*****\n\n",my_x,my_y,my_z,rho,bmod,scale_sou);
		
		cmp[idx                ] = cmp[idx                ] + delta * scale_sou;
		cmp[idx+  one_wf_size_f] = cmp[idx+  one_wf_size_f] + delta * scale_sou;
		cmp[idx+2*one_wf_size_f] = cmp[idx+2*one_wf_size_f] + delta * scale_sou;
	
		//printf("Source Trilinear ix=%d, iy=%d, iz=%d, delta=%e\n",ix,iy,iz,delta);
	}
}

__global__ 
void cuApply_Point_Source_To_TxxTyyTzz(
	float* cmp,
	int x0,
        int y0,
        int nx,
        int ny,
        int nz,
        float dti,
        float ampl1,
        float xs,
        float ys,
        float zs,
        float val,
	bool source_ghost_enabled,
	int ghost_sea_surface_z,
	unsigned int* em,
	float Vp_min,
	float Vp_range,
	float Vs_min,
	float Vs_range,
	float Density_min,
	float Density_range,
	float bmod_ref
        )
{
	int ix = (int)lrintf(xs) - x0;
	int iy = (int)lrintf(ys) - y0;
	int iz = (int)lrintf(zs);

	int one_wf_size_f = nx * nz;
	int one_y_size_f = one_wf_size_f * 6;

	cuApply_Point_Source_To_Single_TxxTyyTzz(
			cmp,nx,ny,nz,
			ix,iy,iz,
			-dti*val*ampl1,
			one_wf_size_f,one_y_size_f,
			em,Vp_min,Vp_range,Vs_min,Vs_range,Density_min,Density_range,bmod_ref
			);
	if (source_ghost_enabled)
	{
		int izd = 2 * ghost_sea_surface_z - iz;
		cuApply_Point_Source_To_Single_TxxTyyTzz(
                        cmp,nx,ny,nz,
                        ix,iy,izd,
                        dti*val*ampl1,
                        one_wf_size_f,one_y_size_f,
			em,Vp_min,Vp_range,Vs_min,Vs_range,Density_min,Density_range,bmod_ref
                        );
	}
}

__device__ 
void _cuApply_Trilinear_Source_To_TxxTyyTzz(
	float* cmp,
	int x0,
        int y0,
        int nx,
        int ny,
        int nz,
        float dti,
        float ampl1,
        float xs,
        float ys,
        float zs,
        float val,
	unsigned int* em,
	float Vp_min,
	float Vp_range,
	float Vs_min,
	float Vs_range,
	float Density_min,
	float Density_range,
	float bmod_ref
        )
{
	int ix = (int)truncf(xs) - x0;
	int iy = (int)truncf(ys) - y0;
	int iz = (int)truncf(zs);

	float xd = 1.0f - (xs - (float)(ix+x0));
	float yd = 1.0f - (ys - (float)(iy+y0));
	float zd = 1.0f - (zs - (float)iz);

	int one_wf_size_f = nx * nz;
	int one_y_size_f = one_wf_size_f * 6;

	cuApply_Point_Source_To_Single_TxxTyyTzz(
			cmp,nx,ny,nz,
			ix,iy,iz,
			-dti*val*ampl1*xd*yd*zd,
			one_wf_size_f,one_y_size_f,
			em,Vp_min,Vp_range,Vs_min,Vs_range,Density_min,Density_range,bmod_ref
			);
	cuApply_Point_Source_To_Single_TxxTyyTzz(
			cmp,nx,ny,nz,
			ix+1,iy,iz,
			-dti*val*ampl1*(1.0f-xd)*yd*zd,
			one_wf_size_f,one_y_size_f,
			em,Vp_min,Vp_range,Vs_min,Vs_range,Density_min,Density_range,bmod_ref
			);
	cuApply_Point_Source_To_Single_TxxTyyTzz(
			cmp,nx,ny,nz,
			ix+1,iy+1,iz,
			-dti*val*ampl1*(1.0f-xd)*(1.0f-yd)*zd,
			one_wf_size_f,one_y_size_f,
			em,Vp_min,Vp_range,Vs_min,Vs_range,Density_min,Density_range,bmod_ref
			);
	cuApply_Point_Source_To_Single_TxxTyyTzz(
			cmp,nx,ny,nz,
			ix,iy+1,iz,
			-dti*val*ampl1*xd*(1.0f-yd)*zd,
			one_wf_size_f,one_y_size_f,
			em,Vp_min,Vp_range,Vs_min,Vs_range,Density_min,Density_range,bmod_ref
			);
	cuApply_Point_Source_To_Single_TxxTyyTzz(
			cmp,nx,ny,nz,
			ix,iy,iz+1,
			-dti*val*ampl1*xd*yd*(1.0f-zd),
			one_wf_size_f,one_y_size_f,
			em,Vp_min,Vp_range,Vs_min,Vs_range,Density_min,Density_range,bmod_ref
			);
	cuApply_Point_Source_To_Single_TxxTyyTzz(
			cmp,nx,ny,nz,
			ix+1,iy,iz+1,
			-dti*val*ampl1*(1.0f-xd)*yd*(1.0f-zd),
			one_wf_size_f,one_y_size_f,
			em,Vp_min,Vp_range,Vs_min,Vs_range,Density_min,Density_range,bmod_ref
			);
	cuApply_Point_Source_To_Single_TxxTyyTzz(
			cmp,nx,ny,nz,
			ix+1,iy+1,iz+1,
			-dti*val*ampl1*(1.0f-xd)*(1.0f-yd)*(1.0f-zd),
			one_wf_size_f,one_y_size_f,
			em,Vp_min,Vp_range,Vs_min,Vs_range,Density_min,Density_range,bmod_ref
			);
	cuApply_Point_Source_To_Single_TxxTyyTzz(
			cmp,nx,ny,nz,
			ix,iy+1,iz+1,
			-dti*val*ampl1*xd*(1.0f-yd)*(1.0f-zd),
			one_wf_size_f,one_y_size_f,
			em,Vp_min,Vp_range,Vs_min,Vs_range,Density_min,Density_range,bmod_ref
			);
}

__global__ 
void cuApply_Trilinear_Source_To_TxxTyyTzz(
	float* cmp,
	int x0,
        int y0,
        int nx,
        int ny,
        int nz,
        float dti,
        float ampl1,
        float xs,
        float ys,
        float zs,
        float val,
	bool source_ghost_enabled,
	float ghost_sea_surface,
	unsigned int* em,
	float Vp_min,
	float Vp_range,
	float Vs_min,
	float Vs_range,
	float Density_min,
	float Density_range,
	float bmod_ref
        )
{
	_cuApply_Trilinear_Source_To_TxxTyyTzz(cmp,x0,y0,nx,ny,nz,dti,ampl1,xs,ys,zs,val,em,Vp_min,Vp_range,Vs_min,Vs_range,Density_min,Density_range,bmod_ref);
	if (source_ghost_enabled)
	{
		_cuApply_Trilinear_Source_To_TxxTyyTzz(cmp,x0,y0,nx,ny,nz,dti,ampl1,xs,ys,2.0f*ghost_sea_surface-zs,-val,em,Vp_min,Vp_range,Vs_min,Vs_range,Density_min,Density_range,bmod_ref);
	}
}

__device__
void cuCompute_DXDYDZ_LoBuf(
        int xoff,
        int yoff,
        int zoff,
        float* dx,
        float* dy,
        float* dz,
        float tmp1,
        float tmp2,
        float tmp3,
        float tmp4,
        float C0,
        float C1,
        float C2,
        float C3,
        float inv_DX,           // 1 / DX
        float inv_DY,           // 1 / DY
        float inv_DZ,           // 1 / DZ
        int one_y_size_f,
        float* buf,
        float* vx_prev,
        int offset
        )
{
        // compute dx
        float Vx_p0 = cuTransposeXZY2XYZ(buf, vx_prev[threadIdx.x+threadIdx.y*32+128]);  // read middle section from persistent Vx buffer
        buf[threadIdx.x+threadIdx.y*96] = tmp1;
        buf[threadIdx.x+threadIdx.y*96+32] = Vx_p0;
        float Vx_p4 = tmp2;
        buf[threadIdx.x+threadIdx.y*96+64] = Vx_p4;
        __syncthreads();
        *dx = (         C0 * (buf[cuCompTXXIdx(1+xoff)] - buf[cuCompTXXIdx(   xoff)]) +
                        C1 * (buf[cuCompTXXIdx(2+xoff)] - buf[cuCompTXXIdx(-1+xoff)]) +
                        C2 * (buf[cuCompTXXIdx(3+xoff)] - buf[cuCompTXXIdx(-2+xoff)]) +
                        C3 * (buf[cuCompTXXIdx(4+xoff)] - buf[cuCompTXXIdx(-3+xoff)]) ) * inv_DX;

        // ..compute dyVx
        float v4 = buf[threadIdx.x+threadIdx.y*96+32];  // read middle section for dyVx from shared memory
        __syncthreads();
        buf[threadIdx.x+threadIdx.y*32+128] = v4;  // deposit middle section
        buf[threadIdx.x+threadIdx.y*32+64*(threadIdx.y&4)] = tmp4;
        __syncthreads();
        *dy = (         C0 * (buf[cuCompTYYIdx(1+yoff)] - buf[cuCompTYYIdx(   yoff)]) +
                        C1 * (buf[cuCompTYYIdx(2+yoff)] - buf[cuCompTYYIdx(-1+yoff)]) +
                        C2 * (buf[cuCompTYYIdx(3+yoff)] - buf[cuCompTYYIdx(-2+yoff)]) +
                        C3 * (buf[cuCompTYYIdx(4+yoff)] - buf[cuCompTYYIdx(-3+yoff)]) ) * inv_DY;

        // ..compute dzVx
        float v5 = cuTransposeXZY2XYZ(buf, tmp3);  // read next 8 z from gmem
        buf[threadIdx.x+threadIdx.y*32] = vx_prev[threadIdx.x+threadIdx.y*32];  // copy 8 deepest z from txz buf
        if (threadIdx.y < 4)
        {
                float v6 = vx_prev[threadIdx.x+threadIdx.y*32+256];
                buf[threadIdx.x+threadIdx.y*32+256] = v6;  // copy 4 shallowest z from txz buf
                buf[threadIdx.x+threadIdx.y*32+384] = v5;  // copy 4 deepest z from next block of txz
                vx_prev[threadIdx.x+threadIdx.y*32] = v6;  // shift txzbuf by 8 z
        }
        // ..store next 8 z in txzbuf
        __syncthreads();
        vx_prev[threadIdx.x+threadIdx.y*32+128] = v5;
        // note that we can use cuCompTYYIdx in place of cuCompTZZIdx after the transpose
        *dz = -(        C0 * (buf[cuCompTYYIdx(1+zoff)] - buf[cuCompTYYIdx(   zoff)]) +
                        C1 * (buf[cuCompTYYIdx(2+zoff)] - buf[cuCompTYYIdx(-1+zoff)]) +
                        C2 * (buf[cuCompTYYIdx(3+zoff)] - buf[cuCompTYYIdx(-2+zoff)]) +
                        C3 * (buf[cuCompTYYIdx(4+zoff)] - buf[cuCompTYYIdx(-3+zoff)]) ) * inv_DZ;
        *dz = cuTransposeXZY2XYZ(buf,*dz);  // this actually transposes back from XYZ to XZY.
}

__device__ 
void cuCompute_DXDYDZ_LoBuf_16(
	int xoff,
	int yoff,
	int zoff,
	float* dx,
        float* dy,
        float* dz,
	float tmp1,
	float tmp2,
	float tmp3,
	float tmp4a,
	float tmp4b,
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
	int one_y_size_f,
	float* buf,
	float* vx_prev,
	int offset
	)
{
	// compute dx
	float Vx_p0 = cuTransposeXZY2XYZ_16(buf, vx_prev[threadIdx.x+threadIdx.y*64+512]);  // read middle section from persistent Vx buffer
	buf[threadIdx.x+threadIdx.y*192] = tmp1;
	buf[threadIdx.x+threadIdx.y*192+64] = Vx_p0;
	float Vx_p8 = tmp2;
	buf[threadIdx.x+threadIdx.y*192+128] = Vx_p8;
	__syncthreads();
	*dx = (         C0 * (buf[cuCompTXXIdx_16(1+xoff)] - buf[cuCompTXXIdx_16(   xoff)]) +
                        C1 * (buf[cuCompTXXIdx_16(2+xoff)] - buf[cuCompTXXIdx_16(-1+xoff)]) +
                        C2 * (buf[cuCompTXXIdx_16(3+xoff)] - buf[cuCompTXXIdx_16(-2+xoff)]) +
                        C3 * (buf[cuCompTXXIdx_16(4+xoff)] - buf[cuCompTXXIdx_16(-3+xoff)]) +
			C4 * (buf[cuCompTXXIdx_16(5+xoff)] - buf[cuCompTXXIdx_16(-4+xoff)]) +
			C5 * (buf[cuCompTXXIdx_16(6+xoff)] - buf[cuCompTXXIdx_16(-5+xoff)]) +
			C6 * (buf[cuCompTXXIdx_16(7+xoff)] - buf[cuCompTXXIdx_16(-6+xoff)]) +
			C7 * (Vx_p8                        - buf[cuCompTXXIdx_16(-7+xoff)]) ) * inv_DX;

	// ..compute dyVx
	float v4 = buf[threadIdx.x+threadIdx.y*192+64];  // read middle section for dyVx from shared memory
	__syncthreads();
	buf[threadIdx.x+threadIdx.y*64] = tmp4a;
	buf[threadIdx.x+threadIdx.y*64+512] = v4;  // deposit middle section
	buf[threadIdx.x+threadIdx.y*64+1024] = tmp4b;
	__syncthreads();
	*dy = (         C0 * (buf[cuCompTYYIdx_16(1+yoff)] - buf[cuCompTYYIdx_16(   yoff)]) +
                        C1 * (buf[cuCompTYYIdx_16(2+yoff)] - buf[cuCompTYYIdx_16(-1+yoff)]) +
                        C2 * (buf[cuCompTYYIdx_16(3+yoff)] - buf[cuCompTYYIdx_16(-2+yoff)]) +
                        C3 * (buf[cuCompTYYIdx_16(4+yoff)] - buf[cuCompTYYIdx_16(-3+yoff)]) +
                        C4 * (buf[cuCompTYYIdx_16(5+yoff)] - buf[cuCompTYYIdx_16(-4+yoff)]) +
                        C5 * (buf[cuCompTYYIdx_16(6+yoff)] - buf[cuCompTYYIdx_16(-5+yoff)]) +
                        C6 * (buf[cuCompTYYIdx_16(7+yoff)] - buf[cuCompTYYIdx_16(-6+yoff)]) +
                        C7 * (buf[cuCompTYYIdx_16(8+yoff)] - buf[cuCompTYYIdx_16(-7+yoff)]) ) * inv_DY;

	// ..compute dzVx
	float v5 = cuTransposeXZY2XYZ_16(buf, tmp3);  // read next 8 z from gmem
	buf[threadIdx.x+threadIdx.y*64] = vx_prev[threadIdx.x+threadIdx.y*64];  // copy 16 deepest z from txz buf
	buf[threadIdx.x+threadIdx.y*64+512] = vx_prev[threadIdx.x+threadIdx.y*64+512];
	buf[threadIdx.x+threadIdx.y*64+1024] = v5;
	// ..shift vx_prev by 8*z
	vx_prev[threadIdx.x+threadIdx.y*64] = vx_prev[threadIdx.x+threadIdx.y*64+512];
	vx_prev[threadIdx.x+threadIdx.y*64+512] = v5;
	__syncthreads();
	// note that we can use cuCompTYYIdx_16 in place of cuCompTZZIdx after the transpose
	*dz = -(        C0 * (buf[cuCompTYYIdx_16(1+zoff)] - buf[cuCompTYYIdx_16(   zoff)]) +
                        C1 * (buf[cuCompTYYIdx_16(2+zoff)] - buf[cuCompTYYIdx_16(-1+zoff)]) +
                        C2 * (buf[cuCompTYYIdx_16(3+zoff)] - buf[cuCompTYYIdx_16(-2+zoff)]) +
                        C3 * (buf[cuCompTYYIdx_16(4+zoff)] - buf[cuCompTYYIdx_16(-3+zoff)]) +
                        C4 * (buf[cuCompTYYIdx_16(5+zoff)] - buf[cuCompTYYIdx_16(-4+zoff)]) +
                        C5 * (buf[cuCompTYYIdx_16(6+zoff)] - buf[cuCompTYYIdx_16(-5+zoff)]) +
                        C6 * (buf[cuCompTYYIdx_16(7+zoff)] - buf[cuCompTYYIdx_16(-6+zoff)]) +
                        C7 * (buf[cuCompTYYIdx_16(8+zoff)] - buf[cuCompTYYIdx_16(-7+zoff)]) ) * inv_DZ;
	*dz = cuTransposeXZY2XYZ_16(buf,*dz);  // this actually transposes back from XYZ to XZY.
}

__global__
#if __CUDA_ARCH__ >= 370
__launch_bounds__(256,8)
#elif __CUDA_ARCH__ >= 300
__launch_bounds__(256,5)
#elif __CUDA_ARCH__ >= 200
__launch_bounds__(256,3)
#endif
void cuPropagate_Stresses_Orthorhombic_Kernel(
        int timestep,
        int x0,                 // x coordinate of westernmost coordinate in block
        int y0,                 // y coordinate of southernmost coordinate in block
        int y1,
        int m1_y0,
        int m1_y1,
        int vol_nx,             // dimensions of global volume
        int vol_ny,
        int vol_nz,
        float dti_half,
        unsigned int* em,       // earth model, 4 interleaved integers. y(0)
        float* cmp,             // txx, tyy, tzz, txy, txz and tyz, middle, t(1), y(0)
        float* m1L,             // Vx, Vy, Vz, Sx, Sy, Sz in that order. left halo, t(0), y(0)
        float* m1C,             // ..middle, t(0), y(0)
        float* m1R,             // ..right halo, t(0), y(0)
        float* m2C,             // txx, tyy, tzz, txy, txz and tyz. middle, t(0), y(0)
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
        int one_wf_size_f,
        int one_y_size_f,
        int em_one_word_size_f,
        int em_one_y_size_f
        )
{
        // work buffer
        __shared__ float buf[768];

        // two persistent buffers used to hold Vx, Vy or Vz values from previous loop iteration
        __shared__ float vx_prev[384];
        __shared__ float vy_prev[384];
        __shared__ float vz_prev[384];

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
                vx_prev[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset] : 0.0f);
                vy_prev[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset+one_wf_size_f] : 0.0f);
                vz_prev[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset+2*one_wf_size_f] : 0.0f);
                if (threadIdx.y < 4)
                {
                        // mirror
                        vx_prev[threadIdx.x+(3-threadIdx.y)*32] = vx_prev[threadIdx.x+(5+threadIdx.y)*32];
                        vy_prev[threadIdx.x+(3-threadIdx.y)*32] = vy_prev[threadIdx.x+(5+threadIdx.y)*32];
                        vz_prev[threadIdx.x+(3-threadIdx.y)*32] = vz_prev[threadIdx.x+(4+threadIdx.y)*32];
                }
        }
        else
        {
                 vx_prev[threadIdx.x+threadIdx.y*32] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset-16] : 0.0f);
                vy_prev[threadIdx.x+threadIdx.y*32] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset-16+one_wf_size_f] : 0.0f);
                vz_prev[threadIdx.x+threadIdx.y*32] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset-16+2*one_wf_size_f] : 0.0f);
                vx_prev[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset] : 0.0f);
                vy_prev[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset+one_wf_size_f] : 0.0f);
                vz_prev[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,y<=m1_y1 ? m1C[offset+2*one_wf_size_f] : 0.0f);
        }
        __syncthreads();

        for (int iZ = 0;  iZ < nz/8;  ++iZ)
        {
                int x = x0 + (threadIdx.x & 3);
                int z = z0 + iZ * 8 + (threadIdx.x / 4);

                float tmp1, tmp5, tmp9;
                if (m1L != 0L && y <= m1_y1)
                {
                        tmp1 = m1L[offset];
                        tmp5 = m1L[offset+one_wf_size_f];
                        tmp9 = m1L[offset+2*one_wf_size_f];
                }
                else
                {
                        tmp1 = tmp5 = tmp9 = 0.0f;
                }

                float tmp2, tmp6, tmpA;
                if (m1R != 0L && y <= m1_y1)
                {
                        tmp2 = m1R[offset];
                        tmp6 = m1R[offset+one_wf_size_f];
                        tmpA = m1R[offset+2*one_wf_size_f];
                }
                else
                {
                        tmp2 = tmp6 = tmpA = 0.0f;
                }

                float tmp3, tmp7, tmpB;
                if (z < vol_nz-8 && y <= m1_y1)
                {
                        tmp3 = m1C[offset+32];
                        tmp7 = m1C[offset+one_wf_size_f+32];
                        tmpB = m1C[offset+2*one_wf_size_f+32];
                }
                else
                {
                        tmp3 = tmp7 = tmpB = 0.0f;
                }

                float tmp4, tmp8, tmpC;
                if (threadIdx.y < 4)
                {
                        if (y-4 >= m1_y0)
                        {
                                tmp4 = m1C[offset-4*one_y_size_f];
                                tmp8 = m1C[offset+one_wf_size_f-4*one_y_size_f];
                                tmpC = m1C[offset+2*one_wf_size_f-4*one_y_size_f];
                        }
                        else
                        {
                                tmp4 = tmp8 = tmpC = 0.0f;
                        }
                }
                else
                {
                        if (y+4 <= m1_y1)
                        {
                                tmp4 = m1C[offset+4*one_y_size_f];
                                tmp8 = m1C[offset+one_wf_size_f+4*one_y_size_f];
                                tmpC = m1C[offset+2*one_wf_size_f+4*one_y_size_f];
                        }
                        else
                        {
                                tmp4 = tmp8 = tmpC = 0.0f;
                        }
                }

                float dxVx, dyVx, dzVx;
                cuCompute_DXDYDZ_LoBuf(0,0,0,&dxVx,&dyVx,&dzVx,tmp1,tmp2,tmp3,tmp4,C0,C1,C2,C3,inv_DX,inv_DY,inv_DZ,one_y_size_f,buf,vx_prev,offset);

                float dxVy, dyVy, dzVy;
                cuCompute_DXDYDZ_LoBuf(-1,-1,0,&dxVy,&dyVy,&dzVy,tmp5,tmp6,tmp7,tmp8,C0,C1,C2,C3,inv_DX,inv_DY,inv_DZ,one_y_size_f,buf,vy_prev,offset);

                float dxVz, dyVz, dzVz;
                cuCompute_DXDYDZ_LoBuf(-1,0,-1,&dxVz,&dyVz,&dzVz,tmp9,tmpA,tmpB,tmpC,C0,C1,C2,C3,inv_DX,inv_DY,inv_DZ,one_y_size_f,buf,vz_prev,offset);

                float dtexx = dxVx;
                float dteyy = dyVy;
                float dtezz = dzVz;
                float dteyz2 = dzVy + dyVz;
                float dtexz2 = dzVx + dxVz;
                float dtexy2 = dyVx + dxVy;

                if (y <= y1)
                {
                        int emIdx = (threadIdx.y + blockIdx.y * 8) * em_one_y_size_f + (iZ*32) + (z0*4) + threadIdx.x;
                        float c11, c22, c33, c44, c55, c66, c12, c13, c23;
                        float dip, azimuth, rake;
                        cuUnpack_And_Compute_On_Kite_CIJs(
                                        em[emIdx], em[emIdx+em_one_word_size_f], em[emIdx+2*em_one_word_size_f], em[emIdx+3*em_one_word_size_f],
                                        Vp_min, Vp_range,
                                        Vs_min, Vs_range,
                                        Density_min, Density_range,
                                        Dip_min, Dip_range,
                                        Azimuth_min, Azimuth_range,
                                        Rake_min, Rake_range,
                                        Delta1_min, Delta1_range,
                                        Delta2_min, Delta2_range,
                                        Delta3_min, Delta3_range,
                                        Epsilon1_min, Epsilon1_range,
                                        Epsilon2_min, Epsilon2_range,
                                        Gamma1_min, Gamma1_range,
                                        Gamma2_min, Gamma2_range,
                                        &dip, &azimuth, &rake,
                                        &c11, &c22, &c33, &c44, &c55, &c66, &c12, &c13, &c23
                                        );

                        // Absorbing boundary decay funct (for Maxwell viscoelastic model):
                        float deta = Compute_ABC(x,y,z,vol_nx,vol_ny,vol_nz,nabc_top,nabc_bot,nabc_sdx,nabc_sdy,vpvert_avtop,vpvert_avbot,inv_DX,inv_DY,inv_DZ);
                        float dabc = (1.0f - 0.5f*deta*dti_half) / (1.0f + 0.5f*deta*dti_half);

                        float old_txx = m2C[offset];
                        float txx = c11 * dtexx + c12 * dteyy + c13 * dtezz;
			txx = (1.0f + dabc) * dti_half * txx + dabc * dabc * old_txx;

                        float old_tyy = m2C[offset+one_wf_size_f];
                        float tyy = c12 * dtexx + c22 * dteyy + c23 * dtezz;
			tyy = (1.0f + dabc) * dti_half * tyy + dabc * dabc * old_tyy;

                        float old_tzz = m2C[offset+2*one_wf_size_f];
                        float tzz = c13 * dtexx + c23 * dteyy + c33 * dtezz;
			tzz = (1.0f + dabc) * dti_half * tzz + dabc * dabc * old_tzz;

                        float old_txy = m2C[offset+3*one_wf_size_f];
                        float txy = c66 * dtexy2;
			txy = (1.0f + dabc) * dti_half * txy + dabc * dabc * old_txy;

                        if (z == 0)
                        {
                                float c13_ = c33 - 2.0f * c55;
                                float dum1 = (c13_ * c13_) / c33;
                                float dum2 = c13_ - dum1;
                                dum1 = c11 - dum1;
                                txx = old_txx + 2.0f * dti_half * (dum1 * dxVx + dum2 * dyVy);
                                tyy = old_tyy + 2.0f * dti_half * (dum2 * dxVx + dum1 * dyVy);
                                // txy = 0.0f;  04/08/15 This was a bug, according to Kurt Nihei
                                tzz = 0.0f;
                        }

                        cmp[offset] = txx;
                        cmp[offset+one_wf_size_f] = tyy;
                        cmp[offset+2*one_wf_size_f] = tzz;
                        cmp[offset+3*one_wf_size_f] = txy;

                        float old_txz = m2C[offset+4*one_wf_size_f];
                        float txz = c55 * dtexz2;
			txz = (1.0f + dabc) * dti_half * txz + dabc * dabc * old_txz;
                        cmp[offset+4*one_wf_size_f] = txz;

                        float old_tyz = m2C[offset+5*one_wf_size_f];
                        float tyz = c44 * dteyz2;
			tyz = (1.0f + dabc) * dti_half * tyz + dabc * dabc * old_tyz;
                        cmp[offset+5*one_wf_size_f] = tyz;
                }

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
void cuPropagate_Stresses_Orthorhombic_Kernel_16(
	int timestep,
	int x0,			// x coordinate of westernmost coordinate in block
	int y0,			// y coordinate of southernmost coordinate in block
	int y1,
	int m1_y0,
	int m1_y1,
	int vol_nx,		// dimensions of global volume
	int vol_ny,
	int vol_nz,
	float dti_half,
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
	int one_wf_size_f,
	int one_y_size_f,
	int em_one_word_size_f,
	int em_one_y_size_f
	)
{
	// work buffer
	__shared__ float buf[1536];
	
	// two persistent buffers used to hold Vx, Vy or Vz values from previous loop iteration
	__shared__ float vx_prev[1024];
	__shared__ float vy_prev[1024];
	__shared__ float vz_prev[1024];

	const int bsX = 8;
	const int bsX_Mask = 7;
	//const int bsX_Shift = 3;
	const int bsZ = 8;

	int z_per_block = (((vol_nz/bsZ) + gridDim.z - 1) / gridDim.z) * bsZ;
	int z0 = z_per_block * blockIdx.z;
	int z1 = z0 + z_per_block - 1;
	if (z1 >= vol_nz) z1 = vol_nz - 1;
	int nz = z1 - z0 + 1;
	if (nz <= 0) return;

        int offset = (threadIdx.y + blockIdx.y * 8) * one_y_size_f + threadIdx.x + z0 * bsX;

        // populate persistent buffers
	int y = y0 + (threadIdx.y + blockIdx.y * 8);
	if (z0 == 0)
	{
		vx_prev[threadIdx.x+threadIdx.y*64+512] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset] : 0.0f);
		vy_prev[threadIdx.x+threadIdx.y*64+512] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset+one_wf_size_f] : 0.0f);
		vz_prev[threadIdx.x+threadIdx.y*64+512] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset+2*one_wf_size_f] : 0.0f);
		vx_prev[threadIdx.x+(7-threadIdx.y)*64] = threadIdx.y == 7 ? 0.0f : vx_prev[threadIdx.x+(9+threadIdx.y)*64];  // vx_prev[z=-8] is never used.
		vy_prev[threadIdx.x+(7-threadIdx.y)*64] = threadIdx.y == 7 ? 0.0f : vy_prev[threadIdx.x+(9+threadIdx.y)*64];  // vy_prev[z=-8] is never used.
		vz_prev[threadIdx.x+(7-threadIdx.y)*64] = vz_prev[threadIdx.x+(8+threadIdx.y)*64];
	}
	else
	{
		vx_prev[threadIdx.x+threadIdx.y*64] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset-64] : 0.0f);
		vy_prev[threadIdx.x+threadIdx.y*64] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset-64+one_wf_size_f] : 0.0f);
		vz_prev[threadIdx.x+threadIdx.y*64] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset-64+2*one_wf_size_f] : 0.0f);
		vx_prev[threadIdx.x+threadIdx.y*64+512] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset] : 0.0f);
		vy_prev[threadIdx.x+threadIdx.y*64+512] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset+one_wf_size_f] : 0.0f);
		vz_prev[threadIdx.x+threadIdx.y*64+512] = cuTransposeXZY2XYZ_16(buf,y<=m1_y1 ? m1C[offset+2*one_wf_size_f] : 0.0f);
	}
	__syncthreads();

	for (int iZ = 0;  iZ < nz/bsZ;  ++iZ)
	{
		int x = x0 + (threadIdx.x & bsX_Mask);
		int z = z0 + iZ * bsZ + (threadIdx.x / bsX);

		float tmp1, tmp5, tmp9;
		if (m1L != 0L && y <= m1_y1)
		{
			tmp1 = m1L[offset];
			tmp5 = m1L[offset+one_wf_size_f];
			tmp9 = m1L[offset+2*one_wf_size_f];
		}
		else
		{
			tmp1 = tmp5 = tmp9 = 0.0f;
		}

		float tmp2, tmp6, tmpA;
		if (m1R != 0L && y <= m1_y1)
		{
			tmp2 = m1R[offset];
			tmp6 = m1R[offset+one_wf_size_f];
			tmpA = m1R[offset+2*one_wf_size_f];
		}
		else
		{
			tmp2 = tmp6 = tmpA = 0.0f;
		}

		float tmp3, tmp7, tmpB;
		if (z < vol_nz-bsZ && y <= m1_y1)
		{
			tmp3 = m1C[offset+64];
			tmp7 = m1C[offset+one_wf_size_f+64];
			tmpB = m1C[offset+2*one_wf_size_f+64];
		}
		else
		{
			tmp3 = tmp7 = tmpB = 0.0f;
		}

		float tmp4a, tmp4b, tmp8a, tmp8b, tmpCa, tmpCb;
		if (y-8 >= m1_y0)
		{
			tmp4a = m1C[offset-8*one_y_size_f];
			tmp8a = m1C[offset+one_wf_size_f-8*one_y_size_f];
			tmpCa = m1C[offset+2*one_wf_size_f-8*one_y_size_f];
		}
		else
		{
			tmp4a = tmp8a = tmpCa = 0.0f;
		}
		if (y+8 <= m1_y1)
		{
			tmp4b = m1C[offset+8*one_y_size_f];
			tmp8b = m1C[offset+one_wf_size_f+8*one_y_size_f];
			tmpCb = m1C[offset+2*one_wf_size_f+8*one_y_size_f];
		}
		else
		{
			tmp4b = tmp8b = tmpCb = 0.0f;
		}

		float dxVx, dyVx, dzVx;
		cuCompute_DXDYDZ_LoBuf_16(0,0,0,&dxVx,&dyVx,&dzVx,tmp1,tmp2,tmp3,tmp4a,tmp4b,C0,C1,C2,C3,C4,C5,C6,C7,inv_DX,inv_DY,inv_DZ,one_y_size_f,buf,vx_prev,offset);

		float dxVy, dyVy, dzVy;
		cuCompute_DXDYDZ_LoBuf_16(-1,-1,0,&dxVy,&dyVy,&dzVy,tmp5,tmp6,tmp7,tmp8a,tmp8b,C0,C1,C2,C3,C4,C5,C6,C7,inv_DX,inv_DY,inv_DZ,one_y_size_f,buf,vy_prev,offset);

		float dxVz, dyVz, dzVz;
		cuCompute_DXDYDZ_LoBuf_16(-1,0,-1,&dxVz,&dyVz,&dzVz,tmp9,tmpA,tmpB,tmpCa,tmpCb,C0,C1,C2,C3,C4,C5,C6,C7,inv_DX,inv_DY,inv_DZ,one_y_size_f,buf,vz_prev,offset);

		float dtexx = dxVx;
		float dteyy = dyVy;
		float dtezz = dzVz;
		float dteyz2 = dzVy + dyVz;
		float dtexz2 = dzVx + dxVz;
		float dtexy2 = dyVx + dxVy;

		if (y <= y1)
		{
			int emIdx = (threadIdx.y + blockIdx.y * 8) * em_one_y_size_f + (iZ*bsZ*bsX) + (z0*bsX) + threadIdx.x;
			float c11, c22, c33, c44, c55, c66, c12, c13, c23;
			float dip, azimuth, rake;
			cuUnpack_And_Compute_On_Kite_CIJs(
					em[emIdx], em[emIdx+em_one_word_size_f], em[emIdx+2*em_one_word_size_f], em[emIdx+3*em_one_word_size_f],
					Vp_min, Vp_range,
					Vs_min, Vs_range,
					Density_min, Density_range,
					Dip_min, Dip_range,
					Azimuth_min, Azimuth_range,
					Rake_min, Rake_range,
					Delta1_min, Delta1_range,
					Delta2_min, Delta2_range,
					Delta3_min, Delta3_range,
					Epsilon1_min, Epsilon1_range,
					Epsilon2_min, Epsilon2_range,
					Gamma1_min, Gamma1_range,
					Gamma2_min, Gamma2_range,
					&dip, &azimuth, &rake,
					&c11, &c22, &c33, &c44, &c55, &c66, &c12, &c13, &c23
					);

			// Absorbing boundary decay funct (for Maxwell viscoelastic model):
			float deta = Compute_ABC(x,y,z,vol_nx,vol_ny,vol_nz,nabc_top,nabc_bot,nabc_sdx,nabc_sdy,vpvert_avtop,vpvert_avbot,inv_DX,inv_DY,inv_DZ);
			float dabc = (1.0f - 0.5f*deta*dti_half) / (1.0f + 0.5f*deta*dti_half);

			float old_txx = m2C[offset];
			float txx = c11 * dtexx + c12 * dteyy + c13 * dtezz;
			txx = (1.0f + dabc) * dti_half * txx + dabc * dabc * old_txx;

			float old_tyy = m2C[offset+one_wf_size_f];
			float tyy = c12 * dtexx + c22 * dteyy + c23 * dtezz;
			tyy = (1.0f + dabc) * dti_half * tyy + dabc * dabc * old_tyy;

			float old_tzz = m2C[offset+2*one_wf_size_f];
			float tzz = c13 * dtexx + c23 * dteyy + c33 * dtezz;
			tzz = (1.0f + dabc) * dti_half * tzz + dabc * dabc * old_tzz;

			float old_txy = m2C[offset+3*one_wf_size_f];
			float txy = c66 * dtexy2;
			txy = (1.0f + dabc) * dti_half * txy + dabc * dabc * old_txy;

			if (z == 0)
			{
				float c13_ = c33 - 2.0f * c55;
				float dum1 = (c13_ * c13_) / c33;
				float dum2 = c13_ - dum1;
				dum1 = c11 - dum1;
				txx = old_txx + 2.0f * dti_half * (dum1 * dxVx + dum2 * dyVy);
				tyy = old_tyy + 2.0f * dti_half * (dum2 * dxVx + dum1 * dyVy);
				// txy = 0.0f;  04/08/15 This was a bug, according to Kurt Nihei
				tzz = 0.0f;
			}

			cmp[offset] = txx;
			cmp[offset+one_wf_size_f] = tyy;
			cmp[offset+2*one_wf_size_f] = tzz;
			cmp[offset+3*one_wf_size_f] = txy;

			float old_txz = m2C[offset+4*one_wf_size_f];
			float txz = c55 * dtexz2;
			txz = (1.0f + dabc) * dti_half * txz + dabc * dabc * old_txz;
			cmp[offset+4*one_wf_size_f] = txz;

			float old_tyz = m2C[offset+5*one_wf_size_f];
			float tyz = c44 * dteyz2;
			tyz = (1.0f + dabc) * dti_half * tyz + dabc * dabc * old_tyz;
			cmp[offset+5*one_wf_size_f] = tyz;
		}

		offset += 64;
	}
}

void Host_Propagate_Stresses_Orthorhombic_Kernel(
	int timestep,
	cudaStream_t stream,
	int spatial_order,	// either 8 or 16
	int num_z,		// number of thread blocks along z axis. must be 1 or more.
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
	bool source_ghost_enabled,
	int ghost_sea_surface_z,
	Elastic_Interpolation_t source_interpolation_method,
	bool is_pressure,
	float ampl1,
	float svaw_sample,
	float xsou,
	float ysou,
	float zsou,
	float bmod_ref
	)
{
	assert(spatial_order == 8 || spatial_order == 16);
	if (!source_ghost_enabled) ghost_sea_surface_z = 0;

	//printf("inject_source=%s, is_pressure=%s, ampl1=%e, svaw_sample=%e, xsou=%f, ysou=%f, zsou=%f\n",inject_source?"Y":"N",is_pressure?"Y":"N",ampl1,svaw_sample,xsou,ysou,zsou);

	int one_wf_size = one_y_size / 6;
	int em_one_word_size = one_wf_size;
	int em_one_y_size = em_one_word_size * 4;

	//printf("nx=%d, ny=%d, nz=%d\n",nx,ny,nz);
	//printf("vpvert_avtop = %f, vpvert_avbot = %f\n",vpvert_avtop,vpvert_avbot);
	//printf("has_low_YHalo = %s, has_high_YHalo = %s\n",has_low_YHalo?"Y":"N",has_high_YHalo?"Y":"N");

	int nx = spatial_order / 2;
	int ny = y1 - y0 + 1;
	int nz = vol_nz;

	if (spatial_order == 8)
	{
		dim3 blockShape(32,8,1);
		dim3 gridShape(1,(ny+7)/8,num_z);
		cuPropagate_Stresses_Orthorhombic_Kernel<<<gridShape,blockShape,0,stream>>>(
				timestep,
				x0,y0,y1,m1_y0,m1_y1,vol_nx,vol_ny,vol_nz,dti/2.0f,
				em,cmp,m1L,m1C,m1R,m2C,
				C8_0,C8_1,C8_2,C8_3,
				inv_DX,inv_DY,inv_DZ,
				vpvert_avtop,vpvert_avbot,
				nabc_sdx,nabc_sdy,nabc_top,nabc_bot,
				Vp_min,Vp_range/65535.0f,
				Vs_min,Vs_range/65535.0f,
				Density_min,Density_range/255.0f,
				Dip_min,Dip_range/255.0f,
				Azimuth_min,Azimuth_range/255.0f,
				Rake_min,Rake_range/255.0f,
				Delta1_min,Delta1_range/255.0f,
				Delta2_min,Delta2_range/255.0f,
				Delta3_min,Delta3_range/255.0f,
				Epsilon1_min,Epsilon1_range/255.0f,
				Epsilon2_min,Epsilon2_range/255.0f,
				Gamma1_min,Gamma1_range/255.0f,
				Gamma2_min,Gamma2_range/255.0f,
				one_wf_size/4,one_y_size/4,
				em_one_word_size/4,em_one_y_size/4
					);
	}
	else if (spatial_order == 16)
	{
		dim3 blockShape(64,8,1);
		dim3 gridShape(1,(ny+7)/8,num_z);
		cuPropagate_Stresses_Orthorhombic_Kernel_16<<<gridShape,blockShape,0,stream>>>(
				timestep,
				x0,y0,y1,m1_y0,m1_y1,vol_nx,vol_ny,vol_nz,dti/2.0f, 
				em,cmp,m1L,m1C,m1R,m2C,
				C16_0,C16_1,C16_2,C16_3,C16_4,C16_5,C16_6,C16_7,
				inv_DX,inv_DY,inv_DZ,
				vpvert_avtop,vpvert_avbot,
				nabc_sdx,nabc_sdy,nabc_top,nabc_bot,
				Vp_min,Vp_range/65535.0f,
				Vs_min,Vs_range/65535.0f,
				Density_min,Density_range/255.0f,
				Dip_min,Dip_range/255.0f,
				Azimuth_min,Azimuth_range/255.0f,
				Rake_min,Rake_range/255.0f,
				Delta1_min,Delta1_range/255.0f,
				Delta2_min,Delta2_range/255.0f,
				Delta3_min,Delta3_range/255.0f,
				Epsilon1_min,Epsilon1_range/255.0f,
				Epsilon2_min,Epsilon2_range/255.0f,
				Gamma1_min,Gamma1_range/255.0f,
				Gamma2_min,Gamma2_range/255.0f,
				one_wf_size/4,one_y_size/4,
				em_one_word_size/4,em_one_y_size/4
					);
	}
#ifdef GPU_DEBUG
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
#endif

	//
	// add source term(s)
	//
	if (inject_source && is_pressure && ampl1 != 0.0f)
	{
		float ghost_sea_surface = (float)ghost_sea_surface_z;
		if (source_interpolation_method == Point)
		{
			dim3 blockShape3(1,1,1);
			dim3 gridShape3(1,1,1);
			cuApply_Point_Source_To_TxxTyyTzz<<<gridShape3,blockShape3,0,stream>>>(
					cmp,x0,y0,nx,ny,nz,dti,ampl1,xsou,ysou,zsou,svaw_sample,source_ghost_enabled,ghost_sea_surface_z,
					em,Vp_min,Vp_range/65535.0f,Vs_min,Vs_range/65535.0f,Density_min,Density_range/255.0f,bmod_ref);
		}
		else if(source_interpolation_method == Trilinear)
		{
			dim3 blockShape3(1,1,1);
			dim3 gridShape3(1,1,1);
			cuApply_Trilinear_Source_To_TxxTyyTzz<<<gridShape3,blockShape3,0,stream>>>(
					cmp,x0,y0,nx,ny,nz,dti,ampl1,xsou,ysou,zsou,svaw_sample,source_ghost_enabled,ghost_sea_surface,
					em,Vp_min,Vp_range/65535.0f,Vs_min,Vs_range/65535.0f,Density_min,Density_range/255.0f,bmod_ref);
		}
		else if(source_interpolation_method == Sinc)
		{
			// use only one thread along z to prevent possible race condition
			dim3 blockShape2(8,8,1);
			dim3 gridShape2(1,1,1);
			cuApply_Source_Term_To_TxxTyyTzz<<<gridShape2,blockShape2,0,stream>>>(
				timestep,cmp,x0,y0,0,nx,ny,nz,dti,ampl1,xsou,ysou,zsou,svaw_sample,source_ghost_enabled,ghost_sea_surface_z,
				em,Vp_min,Vp_range/65535.0f,Vs_min,Vs_range/65535.0f,Density_min,Density_range/255.0f,bmod_ref);
		}
	}
#ifdef GPU_DEBUG
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
#endif
}

