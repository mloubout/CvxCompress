//
// This kernel propagates strain rates one timestep.
// Compile with L1 caching turned off (-Xptxas -dlcm=cg)
//

__device__ 
void _cuApply_Source_Term_To_TxxTyyTzz(
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
		float dx_frac = (float)xs - (float)(icell - 1);
		float dy_frac = (float)ys - (float)(jcell - 1);
                float dz_frac = (float)zs - (float)(kcell - 1);

                float fsinc = cuGen_Sinc_Weight(threadIdx.x,threadIdx.y,thr_z,dx_frac,dy_frac,dz_frac);
		/*
		if (timestep == 3) printf("TIMESTEP %d :: dx_frac ( %d, %d, %d ) = %e\n",timestep,my_x+x0,my_y+y0,my_z,dx_frac);
		if (timestep == 3) printf("TIMESTEP %d :: dy_frac ( %d, %d, %d ) = %e\n",timestep,my_x+x0,my_y+y0,my_z,dy_frac);
		if (timestep == 3) printf("TIMESTEP %d :: dz_frac ( %d, %d, %d ) = %e\n",timestep,my_x+x0,my_y+y0,my_z,dz_frac);
		if (timestep == 3) printf("TIMESTEP %d :: xs ( %d, %d, %d ) = %e\n",timestep,my_x+x0,my_y+y0,my_z,xs);
		if (timestep == 3) printf("TIMESTEP %d :: ys ( %d, %d, %d ) = %e\n",timestep,my_x+x0,my_y+y0,my_z,ys);
		if (timestep == 3) printf("TIMESTEP %d :: zs ( %d, %d, %d ) = %e\n",timestep,my_x+x0,my_y+y0,my_z,zs);
		*/
                if (fsinc != 0.0)
                {
			// mirror source if necessary
			my_z = my_z < 0 ? -my_z : my_z;
		
			// TMJ 05/06/14
			// Mirroring introduces a potential race condition, two threads will update same cell
			// The way to solve this is to let same thread handle all the z indexes

                        int one_wf_size_f = nx * nz;
                        int one_y_size_f = one_wf_size_f * 6;
                        int idx = my_x + my_y * one_y_size_f + my_z * 4;

			cmp[idx                ] = cmp[idx                ] - fsinc * dti * val * ampl1;
			cmp[idx+  one_wf_size_f] = cmp[idx+  one_wf_size_f] - fsinc * dti * val * ampl1;
			cmp[idx+2*one_wf_size_f] = cmp[idx+2*one_wf_size_f] - fsinc * dti * val * ampl1;

			//if (timestep == 1) printf("TIMESTEP 1 :: FSINC ( %d, %d, %d ) = %e\n",my_x+x0,my_y+y0,my_z,fsinc);
                }
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
        float val
        )
{
        // nearest grid point:
        int icell = (int)lrintf(xs) + 1; // to left of extrap pt
        int jcell = (int)lrintf(ys) + 1;
        int kcell = (int)lrintf(zs) + 1; // above interp pt:

	for (int thr_z = 0;  thr_z < 8;  ++thr_z)
	{
        	_cuApply_Source_Term_To_TxxTyyTzz(thr_z,timestep,cmp,x0,y0,z0,nx,ny,nz,dti,ampl1,xs,ys,zs,val,icell,jcell,kcell);
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
	float inv_DX,		// 1 / DX
	float inv_DY,		// 1 / DY
	float inv_DZ,		// 1 / DZ
	int ny,
	int nz,
	int iZ,
        bool do_Lo_YHalo,
	bool has_high_YHalo,
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

__global__
#if __CUDA_ARCH__ >= 300
__launch_bounds__(1280)
#elif __CUDA_ARCH__ >= 200
__launch_bounds__(768)
#endif
void cuPropagate_Stresses_Orthorhombic_Kernel(
	int timestep,
	int x0,			// x coordinate of westernmost coordinate in block
	int y0,			// y coordinate of southernmost coordinate in block
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
	bool has_low_YHalo,	// true if m1 has low yhalo
	bool has_high_YHalo,	// true if m1 has high yhalo
	int nx,
	int ny,
	int nz,
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

        const bool do_Lo_YHalo = (blockIdx.y > 0 || has_low_YHalo) ? true : false;

        int offset = (threadIdx.y + blockIdx.y * 8) * one_y_size_f + threadIdx.x;

        // populate persistent buffers
        vx_prev[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,m1C[offset]);
        vy_prev[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,m1C[offset+one_wf_size_f]);
        vz_prev[threadIdx.x+threadIdx.y*32+128] = cuTransposeXZY2XYZ(buf,m1C[offset+2*one_wf_size_f]);
        if (threadIdx.y < 4)
        {
                vx_prev[threadIdx.x+(3-threadIdx.y)*32] = vx_prev[threadIdx.x+(5+threadIdx.y)*32];
                vy_prev[threadIdx.x+(3-threadIdx.y)*32] = vy_prev[threadIdx.x+(5+threadIdx.y)*32];
                vz_prev[threadIdx.x+(3-threadIdx.y)*32] = vz_prev[threadIdx.x+(4+threadIdx.y)*32];
        }
	__syncthreads();

	for (int iZ = 0;  iZ < nz/8;  ++iZ)
	{
		//int x = x0 + (threadIdx.x & 3);
		//int y = y0 + (threadIdx.y + blockIdx.y * 8);
		//int z = iZ * 8 + (threadIdx.x / 4);

		float tmp1, tmp5, tmp9;
		if (m1L != 0L)
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
		if (m1R != 0L)
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
		if (iZ < ((nz/8)-1))
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
			if (do_Lo_YHalo)
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
			if ((blockIdx.y*8+threadIdx.y+4) < (has_high_YHalo ? ny+4 : ny))
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
		cuCompute_DXDYDZ_LoBuf(0,0,0,&dxVx,&dyVx,&dzVx,tmp1,tmp2,tmp3,tmp4,C0,C1,C2,C3,inv_DX,inv_DY,inv_DZ,ny,nz,iZ,do_Lo_YHalo,has_high_YHalo,one_y_size_f,buf,vx_prev,offset);
                //if (dxVx != 0.0f) printf("x=%d,y=%d,z=%d - dxVx = %f\n",x0+(threadIdx.x&3),y0+threadIdx.y+blockIdx.y*8,iZ*8+(threadIdx.x>>2),dxVx);
                //if (dyVx != 0.0f) printf("x=%d,y=%d,z=%d - dyVx = %f\n",x0+(threadIdx.x&3),y0+threadIdx.y+blockIdx.y*8,iZ*8+(threadIdx.x>>2),dyVx);
                //if (dzVx != 0.0f) printf("x=%d,y=%d,z=%d - dzVx = %f\n",x0+(threadIdx.x&3),y0+threadIdx.y+blockIdx.y*8,iZ*8+(threadIdx.x>>2),dzVx);

		float dxVy, dyVy, dzVy;
		cuCompute_DXDYDZ_LoBuf(-1,-1,0,&dxVy,&dyVy,&dzVy,tmp5,tmp6,tmp7,tmp8,C0,C1,C2,C3,inv_DX,inv_DY,inv_DZ,ny,nz,iZ,do_Lo_YHalo,has_high_YHalo,one_y_size_f,buf,vy_prev,offset);
                //if (dxVy != 0.0f) printf("x=%d,y=%d,z=%d - dxVy = %f\n",x0+(threadIdx.x&3),y0+threadIdx.y+blockIdx.y*8,iZ*8+(threadIdx.x>>2),dxVy);
                //if (dyVy != 0.0f) printf("x=%d,y=%d,z=%d - dyVy = %f\n",x0+(threadIdx.x&3),y0+threadIdx.y+blockIdx.y*8,iZ*8+(threadIdx.x>>2),dyVy);
                //if (dzVy != 0.0f) printf("x=%d,y=%d,z=%d - dzVy = %f\n",x0+(threadIdx.x&3),y0+threadIdx.y+blockIdx.y*8,iZ*8+(threadIdx.x>>2),dzVy);

		float dxVz, dyVz, dzVz;
		cuCompute_DXDYDZ_LoBuf(-1,0,-1,&dxVz,&dyVz,&dzVz,tmp9,tmpA,tmpB,tmpC,C0,C1,C2,C3,inv_DX,inv_DY,inv_DZ,ny,nz,iZ,do_Lo_YHalo,has_high_YHalo,one_y_size_f,buf,vz_prev,offset);
                //if (dxVz != 0.0f) printf("x=%d,y=%d,z=%d - dxVz = %f\n",x0+(threadIdx.x&3),y0+threadIdx.y+blockIdx.y*8,iZ*8+(threadIdx.x>>2),dxVz);
                //if (dyVz != 0.0f) printf("x=%d,y=%d,z=%d - dyVz = %f\n",x0+(threadIdx.x&3),y0+threadIdx.y+blockIdx.y*8,iZ*8+(threadIdx.x>>2),dyVz);
                //if (dzVz != 0.0f) printf("x=%d,y=%d,z=%d - dzVz = %f\n",x0+(threadIdx.x&3),y0+threadIdx.y+blockIdx.y*8,iZ*8+(threadIdx.x>>2),dzVz);
		
		float dtexx = dxVx;
		float dteyy = dyVy;
		float dtezz = dzVz;
		float dteyz2 = dzVy + dyVz;
		float dtexz2 = dzVx + dxVz;
		float dtexy2 = dyVx + dxVy;

		//if (dtexx != 0.0f && timestep == 1) printf("TIMESTEP 1 :: DTEXX (%d, %d, %d) = %e\n",x,y,z,dtexx);
		//if (dteyy != 0.0f && timestep == 1) printf("TIMESTEP 1 :: DTEYY (%d, %d, %d) = %e\n",x,y,z,dteyy);
		//if (dtezz != 0.0f && timestep == 1) printf("TIMESTEP 1 :: DTEZZ (%d, %d, %d) = %e\n",x,y,z,dtezz);
		//if (dteyz2 != 0.0f && timestep == 1) printf("TIMESTEP 1 :: DTEYZ2 (%d, %d, %d) = %e\n",x,y,z,dteyz2);
		//if (dtexz2 != 0.0f && timestep == 1) printf("TIMESTEP 1 :: DTEXZ2 (%d, %d, %d) = %e\n",x,y,z,dtexz2);
		//if (dtexy2 != 0.0f && timestep == 1) printf("TIMESTEP 1 :: DTEXY2 (%d, %d, %d) = %e\n",x,y,z,dtexy2);

		int emIdx = (threadIdx.y + blockIdx.y * 8) * em_one_y_size_f + (iZ*32) + threadIdx.x;
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
		int x = x0 + (threadIdx.x & 3);
		int y = y0 + (threadIdx.y + blockIdx.y * 8);
		int z = iZ * 8 + (threadIdx.x / 4);

		float deta = Compute_ABC(x,y,z,vol_nx,vol_ny,vol_nz,nabc_top,nabc_bot,nabc_sdx,nabc_sdy,vpvert_avtop,vpvert_avbot,inv_DX,inv_DY,inv_DZ);
		float dabc = (1.0f - 0.5f*deta*dti) / (1.0f + 0.5f*deta*dti);

		float old_txx = m2C[offset];
		float txx = c11 * dtexx + c12 * dteyy + c13 * dtezz;
		txx = (dabc + 1.0f) * dti * txx + dabc * dabc * old_txx;

		float old_tyy = m2C[offset+one_wf_size_f];
		float tyy = c12 * dtexx + c22 * dteyy + c23 * dtezz;
		tyy = (dabc + 1.0f) * dti * tyy + dabc * dabc * old_tyy;

		float old_tzz = m2C[offset+2*one_wf_size_f];
		float tzz = c13 * dtexx + c23 * dteyy + c33 * dtezz;
		tzz = (dabc + 1.0f) * dti * tzz + dabc * dabc * old_tzz;

		float old_txy = m2C[offset+3*one_wf_size_f];
		float txy = c66 * dtexy2;
		txy = (dabc + 1.0f) * dti * txy + dabc * dabc * old_txy;

		if (z == 0)
		{
			float c13_ = c33 - 2.0f * c55;
			float dum1 = (c13_ * c13_) / c33;
			float dum2 = c13_ - dum1;
			dum1 = c11 - dum1;
			txx = old_txx + dti * (dum1 * dxVx + dum2 * dyVy);
			tyy = old_tyy + dti * (dum2 * dxVx + dum1 * dyVy);
			txy = 0.0f;
			tzz = 0.0f;
		}

		//float FAKE = c12;
		//FAKE = -FAKE / 3.0f;

		cmp[offset] = txx;
		cmp[offset+one_wf_size_f] = tyy;
		cmp[offset+2*one_wf_size_f] = tzz;
		cmp[offset+3*one_wf_size_f] = txy;

		float old_txz = m2C[offset+4*one_wf_size_f];
		float txz = c55 * dtexz2;
		txz = (dabc + 1.0f) * dti * txz + dabc * dabc * old_txz;
		cmp[offset+4*one_wf_size_f] = txz;

		float old_tyz = m2C[offset+5*one_wf_size_f];
		float tyz = c44 * dteyz2;
		tyz = (dabc + 1.0f) * dti * tyz + dabc * dabc * old_tyz;
		cmp[offset+5*one_wf_size_f] = tyz;

		//if (txx != 0.0f && timestep == 1) printf("TIMESTEP 1 :: TXX (%d, %d, %d) = %e\n",x,y,z,txx);
		//if (tyy != 0.0f && timestep == 1) printf("TIMESTEP 1 :: TYY (%d, %d, %d) = %e\n",x,y,z,tyy);
		//if (tzz != 0.0f && timestep == 1) printf("TIMESTEP 1 :: TZZ (%d, %d, %d) = %e\n",x,y,z,tzz);
		//if (txy != 0.0f && timestep == 1) printf("TIMESTEP 1 :: TXY (%d, %d, %d) = %e\n",x,y,z,txy);
		//if (txz != 0.0f && timestep == 1) printf("TIMESTEP 1 :: TXZ (%d, %d, %d) = %e\n",x,y,z,txz);
		//if (tyz != 0.0f && timestep == 1) printf("TIMESTEP 1 :: TYZ (%d, %d, %d) = %e\n",x,y,z,tyz);

		/*
		if (x == 501 && y == 401 && z == 401)
		{
			printf("\nPropagate_Stress_Isotropic\n");
			printf("--------------------------\n");
			printf("timestep = %d\n",timestep);
			printf("dti = %e\n",dti);
			printf("dtexx = %e\n",dtexx);
			printf("dteyy = %e\n",dteyy);
			printf("dtezz = %e\n",dtezz);
			printf("dteyz2 = %e\n",dteyz2);
			printf("dtexz2 = %e\n",dtexz2);
			printf("dtexy2 = %e\n",dtexy2);
			printf("c11 = %e\n",c11);
			printf("c22 = %e\n",c22);
			printf("c33 = %e\n",c33);
			printf("c44 = %e\n",c44);
			printf("c55 = %e\n",c55);
			printf("c66 = %e\n",c66);
			printf("c12 = %e\n",c12);
			printf("c13 = %e\n",c13);
			printf("c23 = %e\n",c23);
			printf("dabc = %e\n",dabc);
			printf("deta = %e\n",deta);
			printf("txx = %e\n",txx);
			printf("tyy = %e\n",tyy);
			printf("tzz = %e\n",tzz);
			printf("txy = %e\n",txy);
			printf("txz = %e\n",txz);
			printf("tyz = %e\n",tyz);
			printf("\n");
		}
		*/

		offset += 32;
	}
}

void Host_Propagate_Stresses_Orthorhombic_Kernel(
	int timestep,
	cudaStream_t stream,
	int x0,			// x coordinate of westernmost coordinate in block
	int y0,			// y coordinate of southernmost coordinate in block
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
	bool has_low_YHalo,	// true if m1 has low yhalo
	bool has_high_YHalo,	// true if m1 has high yhalo
	int nx,
	int ny,
	int nz,
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
	bool is_pressure,
	float ampl1,
	float svaw_sample,
	float xsou,
	float ysou,
	float zsou
	)
{
	//printf("inject_source=%s, is_pressure=%s, ampl1=%e, svaw_sample=%e, xsou=%f, ysou=%f, zsou=%f\n",inject_source?"Y":"N",is_pressure?"Y":"N",ampl1,svaw_sample,xsou,ysou,zsou);

	int one_wf_size = one_y_size / 6;
	int em_one_word_size = one_wf_size;
	int em_one_y_size = em_one_word_size * 4;

	//printf("nx=%d, ny=%d, nz=%d\n",nx,ny,nz);
	//printf("vpvert_avtop = %f, vpvert_avbot = %f\n",vpvert_avtop,vpvert_avbot);
	//printf("has_low_YHalo = %s, has_high_YHalo = %s\n",has_low_YHalo?"Y":"N",has_high_YHalo?"Y":"N");

	dim3 blockShape(32,8,1);
	dim3 gridShape(1,(ny+7)/8,1);

	cuPropagate_Stresses_Orthorhombic_Kernel<<<gridShape,blockShape,0,stream>>>(
		timestep,
		x0,y0,vol_nx,vol_ny,vol_nz,dti/2.0f,  // TMJ - divide dti by two for the abc
		em,cmp,m1L,m1C,m1R,m2C,
		C0,C1,C2,C3,
		inv_DX,inv_DY,inv_DZ,
		has_low_YHalo,has_high_YHalo,
		nx,ny,nz,
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
#ifdef GPU_DEBUG
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
#endif

	//
	// add source term(s)
	//
	if (inject_source && is_pressure && ampl1 != 0.0f)
	{
		// use only one thread along z to prevent possible race condition
		dim3 blockShape2(8,8,1);
		dim3 gridShape2(1,1,1);
		cuApply_Source_Term_To_TxxTyyTzz<<<gridShape2,blockShape2,0,stream>>>(timestep,cmp,x0,y0,0,nx,ny,nz,dti,ampl1,xsou,ysou,zsou,svaw_sample);
	}
#ifdef GPU_DEBUG
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
#endif
}

