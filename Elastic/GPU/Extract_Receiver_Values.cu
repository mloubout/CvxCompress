#include <cuda_runtime_api.h>
#include "Elastic_Interpolation.hxx"

__device__ 
float cuTrilinear_Interpolation(
	int ix,
	int x0,
	int x1,
	int iy,
	int y0,
	int y1,
	int iz,
	int z1,
	int one_y_size_f,
	float xd,
	float yd,
	float zd,
	float* cmp
	)
{
	bool xlo = ix >= x0;
	bool xhi = ix < x1;
	bool ylo = iy >= y0;
	bool yhi = iy < y1;
	bool zlo = iz >= 0;
	bool zhi = iz < z1;

	float V000 = xlo && ylo && zlo ? cmp[(ix  -x0) + (iy  -y0) * one_y_size_f +  iz    * 4] : 0.0f;
	float V100 = xhi && ylo && zlo ? cmp[(ix+1-x0) + (iy  -y0) * one_y_size_f +  iz    * 4] : 0.0f;
	float V010 = xlo && yhi && zlo ? cmp[(ix  -x0) + (iy+1-y0) * one_y_size_f +  iz    * 4] : 0.0f;
	float V110 = xhi && yhi && zlo ? cmp[(ix+1-x0) + (iy+1-y0) * one_y_size_f +  iz    * 4] : 0.0f;
	float V001 = xlo && ylo && zhi ? cmp[(ix  -x0) + (iy  -y0) * one_y_size_f + (iz+1) * 4] : 0.0f;
	float V101 = xhi && ylo && zhi ? cmp[(ix+1-x0) + (iy  -y0) * one_y_size_f + (iz+1) * 4] : 0.0f;
	float V011 = xlo && yhi && zhi ? cmp[(ix  -x0) + (iy+1-y0) * one_y_size_f + (iz+1) * 4] : 0.0f;
	float V111 = xhi && yhi && zhi ? cmp[(ix+1-x0) + (iy+1-y0) * one_y_size_f + (iz+1) * 4] : 0.0f;

	float c00 = V000 * ( 1.0f - xd ) + V100 * xd;
	float c10 = V010 * ( 1.0f - xd ) + V110 * xd;
	float c01 = V001 * ( 1.0f - xd ) + V101 * xd;
	float c11 = V011 * ( 1.0f - xd ) + V111 * xd;

	float c0 = c00 * ( 1.0f - yd ) + c10 * yd;
	float c1 = c01 * ( 1.0f - yd ) + c11 * yd;

	float c = c0 * ( 1.0f - zd ) + c1 * zd;

	return c;
}

__device__ 
float cuComp_Sum_Of_TXX_TYY_TZZ(
	int ix,
	int iy,
	int iz,
	int x0,
	int y0,
	int one_y_size_f,
	int one_wf_size_f,
	float* cmp
	)
{
	int idx = (ix - x0) + (iy -y0) * one_y_size_f +  iz * 4;
	return cmp[idx] + cmp[idx+one_wf_size_f] + cmp[idx+2*one_wf_size_f];
}

__device__ 
float cuTrilinear_Interpolation_TXX_TYY_TZZ(
	int ix,
	int x0,
	int x1,
	int iy,
	int y0,
	int y1,
	int iz,
	int z1,
	int one_y_size_f,
	int one_wf_size_f,
	float xd,
	float yd,
	float zd,
	float* cmp
	)
{
	bool xlo = ix >= x0;
	bool xhi = ix < x1;
	bool ylo = iy >= y0;
	bool yhi = iy < y1;
	bool zlo = iz >= 0;
	bool zhi = iz < z1;

	float V000 = xlo && ylo && zlo ? cuComp_Sum_Of_TXX_TYY_TZZ(ix  ,iy  ,iz  ,x0,y0,one_y_size_f,one_wf_size_f,cmp) : 0.0f;
	float V100 = xhi && ylo && zlo ? cuComp_Sum_Of_TXX_TYY_TZZ(ix+1,iy  ,iz  ,x0,y0,one_y_size_f,one_wf_size_f,cmp) : 0.0f;
	float V010 = xlo && yhi && zlo ? cuComp_Sum_Of_TXX_TYY_TZZ(ix  ,iy+1,iz  ,x0,y0,one_y_size_f,one_wf_size_f,cmp) : 0.0f;
	float V110 = xhi && yhi && zlo ? cuComp_Sum_Of_TXX_TYY_TZZ(ix+1,iy+1,iz  ,x0,y0,one_y_size_f,one_wf_size_f,cmp) : 0.0f;
	float V001 = xlo && ylo && zhi ? cuComp_Sum_Of_TXX_TYY_TZZ(ix  ,iy  ,iz+1,x0,y0,one_y_size_f,one_wf_size_f,cmp) : 0.0f;
	float V101 = xhi && ylo && zhi ? cuComp_Sum_Of_TXX_TYY_TZZ(ix+1,iy  ,iz+1,x0,y0,one_y_size_f,one_wf_size_f,cmp) : 0.0f;
	float V011 = xlo && yhi && zhi ? cuComp_Sum_Of_TXX_TYY_TZZ(ix  ,iy+1,iz+1,x0,y0,one_y_size_f,one_wf_size_f,cmp) : 0.0f;
	float V111 = xhi && yhi && zhi ? cuComp_Sum_Of_TXX_TYY_TZZ(ix+1,iy+1,iz+1,x0,y0,one_y_size_f,one_wf_size_f,cmp) : 0.0f;

	float c00 = V000 * ( 1.0f - xd ) + V100 * xd;
	float c10 = V010 * ( 1.0f - xd ) + V110 * xd;
	float c01 = V001 * ( 1.0f - xd ) + V101 * xd;
	float c11 = V011 * ( 1.0f - xd ) + V111 * xd;

	float c0 = c00 * ( 1.0f - yd ) + c10 * yd;
	float c1 = c01 * ( 1.0f - yd ) + c11 * yd;

	float c = c0 * ( 1.0f - zd ) + c1 * zd;

	return c;
}

__device__
void cuExtract_Receiver_Particle_Velocity_Values(
        float* cmp,
        int x0,
        int y0,
        int nx,
        int ny,
        int nz,
	int ghost_sea_surface_z,
        float* rcv_binned,
        int num_rx,
        float* res
        )
{
        // linearized thread index
	int tid = blockIdx.x * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x; 
        if (tid < num_rx)
        {
                // determine receiver type flags
                int* iloc = (int*)rcv_binned;
                int num_files = iloc[0];

		// determine my receiver location and selection flags
		bool found_it = false;
		bool rcv_ghost = false;
		float recx, recy, recz;
		int recflags=0, interp=0, tid0=0, outidx=0, offloc=1+2*num_files, offres=0;
		for (int iFile = 0;  iFile < num_files;  ++iFile)
		{	
			int nn = iloc[1+2*iFile];
			int flags = iloc[2+2*iFile] & 14;
			if (flags)
			{
				// this file wants particle velocities
				int num_wf = 0;
				if (flags & 2) ++num_wf;
				if (flags & 4) ++num_wf;
				if (flags & 8) ++num_wf;
				int tid1 = tid0 + nn - 1;
				if (tid >= tid0 && tid <= tid1)
				{
					recflags = flags;
					interp = (iloc[2+2*iFile] >> 16) & 7;
					rcv_ghost = (iloc[2+2*iFile] & (1 << 31)) != 0 ? true : false;
					recx = rcv_binned[offloc+3*(tid-tid0)];
					recy = rcv_binned[offloc+3*(tid-tid0)+1];
					recz = rcv_binned[offloc+3*(tid-tid0)+2];
					outidx = offres + num_wf * (tid-tid0);
					found_it = true;
					break;
				}
				tid0 = tid1 + 1;
				offres += num_wf * nn;
			}
			offloc += 3 * nn;
		}

		if (found_it)
		{
			Elastic_Interpolation_t interpolation_method = (Elastic_Interpolation_t)interp;

			//printf("tid=%d, outidx=%d, offloc=%d, offres=%d, flags=%d, rcv=[%.2f,%.2f,%.2f]\n",tid,outidx,offloc,offres,recflags,recx,recy,recz);
			int one_wf_size_f = nx * nz;
			int one_y_size_f = one_wf_size_f * 6;

			if (interpolation_method == Point)
			{
				// snap to nearest cell
				int ix = (int)lrintf(recx);
				int iy = (int)lrintf(recy);
				int iz = (int)lrintf(recz);
				int idx = (ix-x0) + (iy-y0) * one_y_size_f + iz * 4;
				if (rcv_ghost)
				{
					int izg = 2 * ghost_sea_surface_z - iz;
					int idxg = idx - (iz - izg) * 4;
					if (recflags & 2) res[outidx++] = cmp[idx] - cmp[idxg];
                                        if (recflags & 4) res[outidx++] = cmp[idx+one_wf_size_f] - cmp[idxg+one_wf_size_f];
                                        if (recflags & 8) res[outidx++] = cmp[idx+2*one_wf_size_f] - cmp[idxg+2*one_wf_size_f];
				}
				else
				{
					if (recflags & 2) res[outidx++] = cmp[idx];
					if (recflags & 4) res[outidx++] = cmp[idx+one_wf_size_f];
					if (recflags & 8) res[outidx++] = cmp[idx+2*one_wf_size_f];
				}
			}
			else if (interpolation_method == Trilinear)
			{
				int ix = (int)truncf(recx);
				int iy = (int)truncf(recy);
				int iz = (int)truncf(recz);
				float xd = recx - (float)ix;
				float yd = recy - (float)iy;
				float zd = recz - (float)iz;
				if (rcv_ghost)
				{
					int izg = 2 * ghost_sea_surface_z - iz - 1;
					float zdg = 1.0f - zd;
					if (recflags & 2) res[outidx++] = 
						cuTrilinear_Interpolation(ix,x0,x0+nx-1,iy,y0,y0+ny-1,iz,nz-1,one_y_size_f,xd,yd,zd,cmp)
							+ cuTrilinear_Interpolation(ix,x0,x0+nx-1,iy,y0,y0+ny-1,izg,nz-1,one_y_size_f,xd,yd,zdg,cmp);
					if (recflags & 4) res[outidx++] = 
						cuTrilinear_Interpolation(ix,x0,x0+nx-1,iy,y0,y0+ny-1,iz,nz-1,one_y_size_f,xd,yd,zd,cmp+one_wf_size_f)
							+ cuTrilinear_Interpolation(ix,x0,x0+nx-1,iy,y0,y0+ny-1,izg,nz-1,one_y_size_f,xd,yd,zdg,cmp+one_wf_size_f);
					if (recflags & 8) res[outidx++] = 
						cuTrilinear_Interpolation(ix,x0,x0+nx-1,iy,y0,y0+ny-1,iz,nz-1,one_y_size_f,xd,yd,zd,cmp+2*one_wf_size_f)
							+ cuTrilinear_Interpolation(ix,x0,x0+nx-1,iy,y0,y0+ny-1,izg,nz-1,one_y_size_f,xd,yd,zdg,cmp+2*one_wf_size_f);
				}
				else
				{
					if (recflags & 2) res[outidx++] = cuTrilinear_Interpolation(ix,x0,x0+nx-1,iy,y0,y0+ny-1,iz,nz-1,one_y_size_f,xd,yd,zd,cmp);
					if (recflags & 4) res[outidx++] = cuTrilinear_Interpolation(ix,x0,x0+nx-1,iy,y0,y0+ny-1,iz,nz-1,one_y_size_f,xd,yd,zd,cmp+one_wf_size_f);
					if (recflags & 8) res[outidx++] = cuTrilinear_Interpolation(ix,x0,x0+nx-1,iy,y0,y0+ny-1,iz,nz-1,one_y_size_f,xd,yd,zd,cmp+2*one_wf_size_f);
				}
			}
			else if (interpolation_method == Sinc)
			{
				// nearest grid point:
				int icell = (int)lrintf(recx) + 1; // to left of extrap pt
				int jcell = (int)lrintf(recy) + 1;
				int kcell = (int)lrintf(recz) + 1; // above interp pt:

				// ..fractional distance from grid pt to sou:
				float dx_frac = recx - (float)(icell - 1);
				float dy_frac = recy - (float)(jcell - 1);
				float dz_frac = recz - (float)(kcell - 1);

				int ix_min = icell - 3;
				int ix_max = icell + 4;
				if (ix_min < x0) ix_min = x0;
				if (ix_max >= nx+x0) ix_max = nx + x0 - 1;

				int iy_min = jcell - 3;
				int iy_max = jcell + 4;
				if (iy_min < y0) iy_min = y0;
				if (iy_max >= ny+y0) iy_max = ny + y0 - 1;

				int iz_min = kcell - 3;
				int iz_max = kcell + 4;
				if (iz_min < 0) iz_min = 0;
				if (iz_max >= nz) iz_max = nz - 1;

				float sinc_Dx[8];
				float sinc_Dy[8];
				float sinc_Dz[8];
				for (int ix = ix_min;  ix <= ix_max;  ++ix) sinc_Dx[ix-icell+3] = cuGen_Single_Sinc_Weight(ix-icell+3,dx_frac);
				for (int iy = iy_min;  iy <= iy_max;  ++iy) sinc_Dy[iy-jcell+3] = cuGen_Single_Sinc_Weight(iy-jcell+3,dy_frac);
				for (int iz = iz_min;  iz <= iz_max;  ++iz) sinc_Dz[iz-kcell+3] = cuGen_Single_Sinc_Weight(iz-kcell+3,dz_frac);

				float acc_Vx = 0.0f;
				if (recflags & 2)
				{
					float sinc_Vx[8];
					float vx_dx_frac = recx + 0.5f - (float)(icell - 1);
					for (int ix = ix_min;  ix <= ix_max;  ++ix) sinc_Vx[ix-icell+3] = cuGen_Single_Sinc_Weight(ix-icell+3,vx_dx_frac);
					for (int iz = iz_min;  iz <= iz_max;  ++iz)
					{
						for (int iy = iy_min;  iy <= iy_max;  ++iy)
						{
							for (int ix = ix_min;  ix <= ix_max;  ++ix)
							{
								int idx = (ix-x0) + (iy-y0) * one_y_size_f + iz * 4;
								float fsinc = sinc_Vx[ix-icell+3] * sinc_Dy[iy-jcell+3] * sinc_Dz[iz-kcell+3];
								acc_Vx += fsinc * cmp[idx];
							}
						}
					}
				}

				float acc_Vy = 0.0f;
				if (recflags & 4)
				{
					float sinc_Vy[8];
					float vy_dy_frac = recy - 0.5f - (float)(jcell - 1);
					for (int iy = iy_min;  iy <= iy_max;  ++iy) sinc_Vy[iy-jcell+3] = cuGen_Single_Sinc_Weight(iy-jcell+3,vy_dy_frac);
					for (int iz = iz_min;  iz <= iz_max;  ++iz)
					{
						for (int iy = iy_min;  iy <= iy_max;  ++iy)
						{
							for (int ix = ix_min;  ix <= ix_max;  ++ix)
							{
								int idx = (ix-x0) + (iy-y0) * one_y_size_f + iz * 4;
								float fsinc = sinc_Dx[ix-icell+3] * sinc_Vy[iy-jcell+3] * sinc_Dz[iz-kcell+3];
								acc_Vy += fsinc * cmp[idx+one_wf_size_f];
							}
						}
					}
				}

				float acc_Vz = 0.0f;
				if (recflags & 8)
				{
					float sinc_Vz[8];
					float vz_dz_frac = recz - 0.5f - (float)(kcell - 1);
					for (int iz = iz_min;  iz <= iz_max;  ++iz) sinc_Vz[iz-kcell+3] = cuGen_Single_Sinc_Weight(iz-kcell+3,vz_dz_frac);
					for (int iz = iz_min;  iz <= iz_max;  ++iz)
					{
						for (int iy = iy_min;  iy <= iy_max;  ++iy)
						{
							for (int ix = ix_min;  ix <= ix_max;  ++ix)
							{
								int idx = (ix-x0) + (iy-y0) * one_y_size_f + iz * 4;
								float fsinc = sinc_Dx[ix-icell+3] * sinc_Dy[iy-jcell+3] * sinc_Vz[iz-kcell+3];
								acc_Vz += fsinc * cmp[idx+2*one_wf_size_f];
							}
						}
					}
				}

				if (rcv_ghost)
				{
					int kcell_ghost = 2 * ghost_sea_surface_z - kcell - 1;
					float dz_frac_ghost = 1.0f - dz_frac;

					int iz_min_ghost = kcell_ghost - 3;
					int iz_max_ghost = kcell_ghost + 4;
					if (iz_min_ghost < 0) iz_min_ghost = 0;
					if (iz_max_ghost >= nz) iz_max_ghost = nz - 1;

					float sinc_Dz_ghost[8];
					for (int iz = iz_min_ghost;  iz <= iz_max_ghost;  ++iz) sinc_Dz_ghost[iz-kcell_ghost+3] = cuGen_Single_Sinc_Weight(iz-kcell_ghost+3,dz_frac_ghost);
					
					if (recflags & 2)
					{
						float sinc_Vx[8];
						float vx_dx_frac = recx + 0.5f - (float)(icell - 1);
						for (int ix = ix_min;  ix <= ix_max;  ++ix) sinc_Vx[ix-icell+3] = cuGen_Single_Sinc_Weight(ix-icell+3,vx_dx_frac);
						for (int iz = iz_min_ghost;  iz <= iz_max_ghost;  ++iz)
						{
							for (int iy = iy_min;  iy <= iy_max;  ++iy)
							{
								for (int ix = ix_min;  ix <= ix_max;  ++ix)
								{
									int idx = (ix-x0) + (iy-y0) * one_y_size_f + iz * 4;
									float fsinc = sinc_Vx[ix-icell+3] * sinc_Dy[iy-jcell+3] * sinc_Dz_ghost[iz-kcell_ghost+3];
									acc_Vx -= fsinc * cmp[idx];
								}
							}
						}
					}

					if (recflags & 4)
					{
						float sinc_Vy[8];
						float vy_dy_frac = recy - 0.5f - (float)(jcell - 1);
						for (int iy = iy_min;  iy <= iy_max;  ++iy) sinc_Vy[iy-jcell+3] = cuGen_Single_Sinc_Weight(iy-jcell+3,vy_dy_frac);
						for (int iz = iz_min_ghost;  iz <= iz_max_ghost;  ++iz)
						{
							for (int iy = iy_min;  iy <= iy_max;  ++iy)
							{
								for (int ix = ix_min;  ix <= ix_max;  ++ix)
								{
									int idx = (ix-x0) + (iy-y0) * one_y_size_f + iz * 4;
									float fsinc = sinc_Dx[ix-icell+3] * sinc_Vy[iy-jcell+3] * sinc_Dz[iz-kcell_ghost+3];
									acc_Vy -= fsinc * cmp[idx+one_wf_size_f];
								}
							}
						}
					}

					if (recflags & 8)
					{
						float sinc_Vz_ghost[8];
						float vz_dz_frac = recz - 0.5f - (float)(kcell - 1);
						float vz_dz_frac_ghost = 1.0f - vz_dz_frac;
						for (int iz = iz_min_ghost;  iz <= iz_max_ghost;  ++iz) sinc_Vz_ghost[iz-kcell_ghost+3] = cuGen_Single_Sinc_Weight(iz-kcell_ghost+3,vz_dz_frac_ghost);
						for (int iz = iz_min_ghost;  iz <= iz_max_ghost;  ++iz)
						{
							for (int iy = iy_min;  iy <= iy_max;  ++iy)
							{
								for (int ix = ix_min;  ix <= ix_max;  ++ix)
								{
									int idx = (ix-x0) + (iy-y0) * one_y_size_f + iz * 4;
									float fsinc = sinc_Dx[ix-icell+3] * sinc_Dy[iy-jcell+3] * sinc_Vz_ghost[iz-kcell_ghost+3];
									acc_Vz -= fsinc * cmp[idx+2*one_wf_size_f];
								}
							}
						}
					}
				}

				if (recflags & 2) res[outidx++] = acc_Vx;
				if (recflags & 4) res[outidx++] = acc_Vy;
				if (recflags & 8) res[outidx++] = acc_Vz;
			}
		}
	}
}

__device__
void cuExtract_Receiver_Pressure_Values(
        float* cmp,
        int x0,
        int y0,
        int nx,
        int ny,
        int nz,
	int ghost_sea_surface_z,
        float* rcv_binned,
        int num_rx,
        float* res
        )
{
        // linearized thread index
	int tid = blockIdx.x * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;
        if (tid < num_rx)
        {
                // determine receiver type flags
                int* iloc = (int*)rcv_binned;
                int num_files = iloc[0];

		// determine my receiver location and selection flags
		bool found_it = false;
		bool rcv_ghost = false;
		float recx, recy, recz;
		int tid0=0, interp=0, outidx=0, offloc=1+2*num_files, offres=0;
		for (int iFile = 0;  iFile < num_files;  ++iFile)
		{	
			int nn = iloc[1+2*iFile];
			int flags = iloc[2+2*iFile] & 1;
			if (flags)
			{
				// this file wants receiver pressure
				int tid1 = tid0 + nn - 1;
				if (tid >= tid0 && tid <= tid1)
				{
					interp = (iloc[2+2*iFile] >> 16) & 7;
					rcv_ghost = (iloc[2+2*iFile] & (1 << 31)) != 0 ? true : false;
					recx = rcv_binned[offloc+3*(tid-tid0)];
					recy = rcv_binned[offloc+3*(tid-tid0)+1];
					recz = rcv_binned[offloc+3*(tid-tid0)+2];
					outidx = offres + (tid-tid0);
					found_it = true;
					break;
				}
				tid0 = tid1 + 1;
				offres += nn;
			}
			offloc += 3 * nn;
		}

		if (found_it)
		{
			Elastic_Interpolation_t interpolation_method = (Elastic_Interpolation_t)interp;
			//printf("tid=%d, outidx=%d, offloc=%d, offres=%d, rcv=[%.2f,%.2f,%.2f]\n",tid,outidx,offloc,offres,recx,recy,recz);

			int one_wf_size_f = nx * nz;
			int one_y_size_f = one_wf_size_f * 6;

			if (interpolation_method == Point)
			{
				// snap to nearest cell
				int ix = (int)lrintf(recx);
				int iy = (int)lrintf(recy);
				int iz = (int)lrintf(recz);
				int idx = (ix-x0) + (iy-y0) * one_y_size_f + iz * 4;
				if (rcv_ghost)
				{
					int izg = 2* ghost_sea_surface_z - iz;
					int idxg = idx - (iz - izg) * 4;
					res[outidx++] = -(cmp[idx] - cmp[idxg] + cmp[idx+one_wf_size_f] - cmp[idxg+one_wf_size_f] + cmp[idx+2*one_wf_size_f] - cmp[idxg+2*one_wf_size_f]) / 3.0f;
				}
				else
				{
					res[outidx++] = -(cmp[idx] + cmp[idx+one_wf_size_f] + cmp[idx+2*one_wf_size_f]) / 3.0f;
				}
			}
			else if (interpolation_method == Trilinear)
			{
				int ix = (int)truncf(recx);
				int iy = (int)truncf(recy);
				int iz = (int)truncf(recz);
				float xd = recx - (float)ix;
				float yd = recy - (float)iy;
				float zd = recz - (float)iz;
				if (rcv_ghost)
				{
					int izg = 2 * ghost_sea_surface_z - iz - 1;
					float zdg = 1.0f - zd;
					res[outidx++] = 
						-(cuTrilinear_Interpolation_TXX_TYY_TZZ(ix,x0,x0+nx-1,iy,y0,y0+ny-1,iz,nz-1,one_y_size_f,one_wf_size_f,xd,yd,zd,cmp)
								+cuTrilinear_Interpolation_TXX_TYY_TZZ(ix,x0,x0+nx-1,iy,y0,y0+ny-1,izg,nz-1,one_y_size_f,one_wf_size_f,xd,yd,zdg,cmp)) / 3.0f;
				}
				else
				{
					res[outidx++] = -(cuTrilinear_Interpolation_TXX_TYY_TZZ(ix,x0,x0+nx-1,iy,y0,y0+ny-1,iz,nz-1,one_y_size_f,one_wf_size_f,xd,yd,zd,cmp)) / 3.0f;
				}
			}
			else if (interpolation_method == Sinc)
			{
				// nearest grid point:
				int icell = (int)lrintf(recx) + 1; // to left of extrap pt
				int jcell = (int)lrintf(recy) + 1;
				int kcell = (int)lrintf(recz) + 1; // above interp pt:

				// ..fractional distance from grid pt to sou:
				float dx_frac = recx - (float)(icell - 1);
				float dy_frac = recy - (float)(jcell - 1);
				float dz_frac = recz - (float)(kcell - 1);

				int ix_min = icell - 3;
				int ix_max = icell + 4;
				if (ix_min < x0) ix_min = x0;
				if (ix_max >= nx+x0) ix_max = nx + x0 - 1;

				int iy_min = jcell - 3;
				int iy_max = jcell + 4;
				if (iy_min < y0) iy_min = y0;
				if (iy_max >= ny+y0) iy_max = ny + y0 - 1;

				int iz_min = kcell - 3;
				int iz_max = kcell + 4;
				if (iz_min < 0) iz_min = 0;
				if (iz_max >= nz) iz_max = nz - 1;

				float sinc_Dx[8];
				float sinc_Dy[8];
				float sinc_Dz[8];
				for (int ii = 0;  ii < 8;  ++ii) sinc_Dx[ii] = sinc_Dy[ii] = sinc_Dz[ii] = 0.0f;
				for (int ix = ix_min;  ix <= ix_max;  ++ix) sinc_Dx[ix-icell+3] = cuGen_Single_Sinc_Weight(ix-icell+3,dx_frac);
				for (int iy = iy_min;  iy <= iy_max;  ++iy) sinc_Dy[iy-jcell+3] = cuGen_Single_Sinc_Weight(iy-jcell+3,dy_frac);
				for (int iz = iz_min;  iz <= iz_max;  ++iz) sinc_Dz[iz-kcell+3] = cuGen_Single_Sinc_Weight(iz-kcell+3,dz_frac);

				float acc = 0.0f;
				for (int iz = iz_min;  iz <= iz_max;  ++iz)
				{
					for (int iy = iy_min;  iy <= iy_max;  ++iy)
					{
						for (int ix = ix_min;  ix <= ix_max;  ++ix)
						{
							int idx = (ix-x0) + (iy-y0) * one_y_size_f + iz * 4;
							float fsinc = sinc_Dx[ix-icell+3] * sinc_Dy[iy-jcell+3] * sinc_Dz[iz-kcell+3];
							acc += fsinc * (cmp[idx] + cmp[idx+one_wf_size_f] + cmp[idx+2*one_wf_size_f]);
							//printf("rec=[%f,%f,%f] = %f, x0=%d, y0=%d, ix=%d, iy=%d, iz=%d, fsinc_x=%f, fsinc_y=%f, fsinc_z=%f, txx=%f, tyy=%f, tzz=%f\n",recx,recy,recz,acc,x0,y0,ix,iy,iz,sinc_Dx[ix-icell+3],sinc_Dy[iy-jcell+3],sinc_Dz[iz-kcell+3],cmp[idx],cmp[idx+one_wf_size_f],cmp[idx+2*one_wf_size_f]);
						}
					}
				}
				if (rcv_ghost)
				{
					int kcell_ghost = 2 * ghost_sea_surface_z - kcell - 1;
					float dz_frac_ghost = 1.0f - dz_frac;

					int iz_min_ghost = kcell_ghost - 3;
					int iz_max_ghost = kcell_ghost + 4;
					if (iz_min_ghost < 0) iz_min_ghost = 0;
					if (iz_max_ghost >= nz) iz_max_ghost = nz - 1;

					float sinc_Dz_ghost[8];
					for (int iz = iz_min_ghost;  iz <= iz_max_ghost;  ++iz) sinc_Dz_ghost[iz-kcell_ghost+3] = cuGen_Single_Sinc_Weight(iz-kcell_ghost+3,dz_frac_ghost);

					for (int iz = iz_min_ghost;  iz <= iz_max_ghost;  ++iz)
					{
						for (int iy = iy_min;  iy <= iy_max;  ++iy)
						{
							for (int ix = ix_min;  ix <= ix_max;  ++ix)
							{
								int idx = (ix-x0) + (iy-y0) * one_y_size_f + iz * 4;
								float fsinc = sinc_Dx[ix-icell+3] * sinc_Dy[iy-jcell+3] * sinc_Dz_ghost[iz-kcell_ghost+3];
								acc -= fsinc * (cmp[idx] + cmp[idx+one_wf_size_f] + cmp[idx+2*one_wf_size_f]);
							}
						}
					}
				}
				acc = (-acc) / 3.0f;
				//printf("sinc_Dx = %f,%f,%f,%f,%f,%f,%f,%f\n",sinc_Dx[0],sinc_Dx[1],sinc_Dx[2],sinc_Dx[3],sinc_Dx[4],sinc_Dx[5],sinc_Dx[6],sinc_Dx[7]);
				//printf("sinc_Dy = %f,%f,%f,%f,%f,%f,%f,%f\n",sinc_Dy[0],sinc_Dy[1],sinc_Dy[2],sinc_Dy[3],sinc_Dy[4],sinc_Dy[5],sinc_Dy[6],sinc_Dy[7]);
				//printf("sinc_Dz = %f,%f,%f,%f,%f,%f,%f,%f\n",sinc_Dz[0],sinc_Dz[1],sinc_Dz[2],sinc_Dz[3],sinc_Dz[4],sinc_Dz[5],sinc_Dz[6],sinc_Dz[7]);
				//printf("rec=[%f,%f,%f] = %f, x0=%d, y0=%d, ix=[%d,%d], iy[%d,%d], iz=[%d,%d]\n",recx,recy,recz,acc,x0,y0,ix_min,ix_max,iy_min,iy_max,iz_min,iz_max);
				res[outidx++] = acc;
			}
		}
	}
}

__global__
void cuExtract_Receiver_Values(
	int is_vp_0, int is_vp_1, int is_vp_2, int is_vp_3, int is_vp_4, int is_vp_5, int is_vp_6, int is_vp_7, 
	int is_vp_8, int is_vp_9, int is_vp_10, int is_vp_11, int is_vp_12, int is_vp_13, int is_vp_14, int is_vp_15, 
	float* cmp_0, float* cmp_1, float* cmp_2, float* cmp_3, float* cmp_4, float* cmp_5, float* cmp_6, float* cmp_7, 
	float* cmp_8, float* cmp_9, float* cmp_10, float* cmp_11, float* cmp_12, float* cmp_13, float* cmp_14, float* cmp_15, 
	int x0_0, int x0_1, int x0_2, int x0_3, int x0_4, int x0_5, int x0_6, int x0_7, 
	int x0_8, int x0_9, int x0_10, int x0_11, int x0_12, int x0_13, int x0_14, int x0_15, 
	int y0_0, int y0_1, int y0_2, int y0_3, int y0_4, int y0_5, int y0_6, int y0_7, 
	int y0_8, int y0_9, int y0_10, int y0_11, int y0_12, int y0_13, int y0_14, int y0_15, 
	int nx_0, int nx_1, int nx_2, int nx_3, int nx_4, int nx_5, int nx_6, int nx_7, 
	int nx_8, int nx_9, int nx_10, int nx_11, int nx_12, int nx_13, int nx_14, int nx_15, 
	int ny_0, int ny_1, int ny_2, int ny_3, int ny_4, int ny_5, int ny_6, int ny_7, 
	int ny_8, int ny_9, int ny_10, int ny_11, int ny_12, int ny_13, int ny_14, int ny_15, 
	int nz_0, int nz_1, int nz_2, int nz_3, int nz_4, int nz_5, int nz_6, int nz_7, 
	int nz_8, int nz_9, int nz_10, int nz_11, int nz_12, int nz_13, int nz_14, int nz_15, 
	int ghost_sea_surface_z,
	float* rcv_binned_0, float* rcv_binned_1, float* rcv_binned_2, float* rcv_binned_3, float* rcv_binned_4, float* rcv_binned_5, float* rcv_binned_6, float* rcv_binned_7, 
	float* rcv_binned_8, float* rcv_binned_9, float* rcv_binned_10, float* rcv_binned_11, float* rcv_binned_12, float* rcv_binned_13, float* rcv_binned_14, float* rcv_binned_15, 
	int num_rx_0, int num_rx_1, int num_rx_2, int num_rx_3, int num_rx_4, int num_rx_5, int num_rx_6, int num_rx_7, 
	int num_rx_8, int num_rx_9, int num_rx_10, int num_rx_11, int num_rx_12, int num_rx_13, int num_rx_14, int num_rx_15, 
	float* res_0, float* res_1, float* res_2, float* res_3, float* res_4, float* res_5, float* res_6, float* res_7, 
	float* res_8, float* res_9, float* res_10, float* res_11, float* res_12, float* res_13, float* res_14, float* res_15
	)
{
	if (blockIdx.y == 0)
	{
		if (is_vp_0)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_0,x0_0,y0_0,nx_0,ny_0,nz_0,ghost_sea_surface_z,rcv_binned_0,num_rx_0,res_0);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_0,x0_0,y0_0,nx_0,ny_0,nz_0,ghost_sea_surface_z,rcv_binned_0,num_rx_0,res_0);
		}
	}
	else if (blockIdx.y == 1)
	{
		if (is_vp_1)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_1,x0_1,y0_1,nx_1,ny_1,nz_1,ghost_sea_surface_z,rcv_binned_1,num_rx_1,res_1);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_1,x0_1,y0_1,nx_1,ny_1,nz_1,ghost_sea_surface_z,rcv_binned_1,num_rx_1,res_1);
		}
	}
	else if (blockIdx.y == 2)
	{
		if (is_vp_2)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_2,x0_2,y0_2,nx_2,ny_2,nz_2,ghost_sea_surface_z,rcv_binned_2,num_rx_2,res_2);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_2,x0_2,y0_2,nx_2,ny_2,nz_2,ghost_sea_surface_z,rcv_binned_2,num_rx_2,res_2);
		}
	}
	else if (blockIdx.y == 3)
	{
		if (is_vp_3)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_3,x0_3,y0_3,nx_3,ny_3,nz_3,ghost_sea_surface_z,rcv_binned_3,num_rx_3,res_3);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_3,x0_3,y0_3,nx_3,ny_3,nz_3,ghost_sea_surface_z,rcv_binned_3,num_rx_3,res_3);
		}
	}
	else if (blockIdx.y == 4)
	{
		if (is_vp_4)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_4,x0_4,y0_4,nx_4,ny_4,nz_4,ghost_sea_surface_z,rcv_binned_4,num_rx_4,res_4);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_4,x0_4,y0_4,nx_4,ny_4,nz_4,ghost_sea_surface_z,rcv_binned_4,num_rx_4,res_4);
		}
	}
	else if (blockIdx.y == 5)
	{
		if (is_vp_5)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_5,x0_5,y0_5,nx_5,ny_5,nz_5,ghost_sea_surface_z,rcv_binned_5,num_rx_5,res_5);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_5,x0_5,y0_5,nx_5,ny_5,nz_5,ghost_sea_surface_z,rcv_binned_5,num_rx_5,res_5);
		}
	}
	else if (blockIdx.y == 6)
	{
		if (is_vp_6)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_6,x0_6,y0_6,nx_6,ny_6,nz_6,ghost_sea_surface_z,rcv_binned_6,num_rx_6,res_6);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_6,x0_6,y0_6,nx_6,ny_6,nz_6,ghost_sea_surface_z,rcv_binned_6,num_rx_6,res_6);
		}
	}
	else if (blockIdx.y == 7)
	{
		if (is_vp_7)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_7,x0_7,y0_7,nx_7,ny_7,nz_7,ghost_sea_surface_z,rcv_binned_7,num_rx_7,res_7);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_7,x0_7,y0_7,nx_7,ny_7,nz_7,ghost_sea_surface_z,rcv_binned_7,num_rx_7,res_7);
		}
	}
	else if (blockIdx.y == 8)
	{
		if (is_vp_8)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_8,x0_8,y0_8,nx_8,ny_8,nz_8,ghost_sea_surface_z,rcv_binned_8,num_rx_8,res_8);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_8,x0_8,y0_8,nx_8,ny_8,nz_8,ghost_sea_surface_z,rcv_binned_8,num_rx_8,res_8);
		}
	}
	else if (blockIdx.y == 9)
	{
		if (is_vp_9)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_9,x0_9,y0_9,nx_9,ny_9,nz_9,ghost_sea_surface_z,rcv_binned_9,num_rx_9,res_9);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_9,x0_9,y0_9,nx_9,ny_9,nz_9,ghost_sea_surface_z,rcv_binned_9,num_rx_9,res_9);
		}
	}
	else if (blockIdx.y == 10)
	{
		if (is_vp_10)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_10,x0_10,y0_10,nx_10,ny_10,nz_10,ghost_sea_surface_z,rcv_binned_10,num_rx_10,res_10);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_10,x0_10,y0_10,nx_10,ny_10,nz_10,ghost_sea_surface_z,rcv_binned_10,num_rx_10,res_10);
		}
	}
	else if (blockIdx.y == 11)
	{
		if (is_vp_11)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_11,x0_11,y0_11,nx_11,ny_11,nz_11,ghost_sea_surface_z,rcv_binned_11,num_rx_11,res_11);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_11,x0_11,y0_11,nx_11,ny_11,nz_11,ghost_sea_surface_z,rcv_binned_11,num_rx_11,res_11);
		}
	}
	else if (blockIdx.y == 12)
	{
		if (is_vp_12)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_12,x0_12,y0_12,nx_12,ny_12,nz_12,ghost_sea_surface_z,rcv_binned_12,num_rx_12,res_12);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_12,x0_12,y0_12,nx_12,ny_12,nz_12,ghost_sea_surface_z,rcv_binned_12,num_rx_12,res_12);
		}
	}
	else if (blockIdx.y == 13)
	{
		if (is_vp_13)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_13,x0_13,y0_13,nx_13,ny_13,nz_13,ghost_sea_surface_z,rcv_binned_13,num_rx_13,res_13);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_13,x0_13,y0_13,nx_13,ny_13,nz_13,ghost_sea_surface_z,rcv_binned_13,num_rx_13,res_13);
		}
	}
	else if (blockIdx.y == 14)
	{
		if (is_vp_14)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_14,x0_14,y0_14,nx_14,ny_14,nz_14,ghost_sea_surface_z,rcv_binned_14,num_rx_14,res_14);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_14,x0_14,y0_14,nx_14,ny_14,nz_14,ghost_sea_surface_z,rcv_binned_14,num_rx_14,res_14);
		}
	}
	else if (blockIdx.y == 15)
	{
		if (is_vp_15)
		{
			cuExtract_Receiver_Particle_Velocity_Values(cmp_15,x0_15,y0_15,nx_15,ny_15,nz_15,ghost_sea_surface_z,rcv_binned_15,num_rx_15,res_15);
		}
		else
		{
			cuExtract_Receiver_Pressure_Values(cmp_15,x0_15,y0_15,nx_15,ny_15,nz_15,ghost_sea_surface_z,rcv_binned_15,num_rx_15,res_15);
		}
	}
}

void 
Host_Extract_Receiver_Values(
	cudaStream_t stream,
	int* is_vp,
        float** cmp,
        int* x0,
        int* y0,
        int* nx,
        int* ny,
        int* nz,
	int ghost_sea_surface_z,
        float** rcv_binned,
        int* num_rx,
        float** res,
	int num_kernels
        )
{
	// find max(num_rx)
	int max_num_rx = num_rx[0];
	for (int i = 1;  i < num_kernels;  ++i) if (num_rx[i] > max_num_rx) max_num_rx = num_rx[i];	

	dim3 blockShape(32,4,1);
        dim3 gridShape((max_num_rx+127)/128,num_kernels,1);
	cuExtract_Receiver_Values<<<gridShape,blockShape,0,stream>>>(
		is_vp[0], is_vp[1], is_vp[2], is_vp[3], is_vp[4], is_vp[5], is_vp[6], is_vp[7], is_vp[8], is_vp[9], is_vp[10], is_vp[11], is_vp[12], is_vp[13], is_vp[14], is_vp[15],
		cmp[0], cmp[1], cmp[2], cmp[3], cmp[4], cmp[5], cmp[6], cmp[7], cmp[8], cmp[9], cmp[10], cmp[11], cmp[12], cmp[13], cmp[14], cmp[15],
		x0[0], x0[1], x0[2], x0[3], x0[4], x0[5], x0[6], x0[7], x0[8], x0[9], x0[10], x0[11], x0[12], x0[13], x0[14], x0[15],
		y0[0], y0[1], y0[2], y0[3], y0[4], y0[5], y0[6], y0[7], y0[8], y0[9], y0[10], y0[11], y0[12], y0[13], y0[14], y0[15],
		nx[0], nx[1], nx[2], nx[3], nx[4], nx[5], nx[6], nx[7], nx[8], nx[9], nx[10], nx[11], nx[12], nx[13], nx[14], nx[15],
		ny[0], ny[1], ny[2], ny[3], ny[4], ny[5], ny[6], ny[7], ny[8], ny[9], ny[10], ny[11], ny[12], ny[13], ny[14], ny[15],
		nz[0], nz[1], nz[2], nz[3], nz[4], nz[5], nz[6], nz[7], nz[8], nz[9], nz[10], nz[11], nz[12], nz[13], nz[14], nz[15],
		ghost_sea_surface_z,
		rcv_binned[0], rcv_binned[1], rcv_binned[2], rcv_binned[3], rcv_binned[4], rcv_binned[5], rcv_binned[6], rcv_binned[7], 
		rcv_binned[8], rcv_binned[9], rcv_binned[10], rcv_binned[11], rcv_binned[12], rcv_binned[13], rcv_binned[14], rcv_binned[15],
		num_rx[0], num_rx[1], num_rx[2], num_rx[3], num_rx[4], num_rx[5], num_rx[6], num_rx[7], num_rx[8], num_rx[9], num_rx[10], num_rx[11], num_rx[12], num_rx[13], num_rx[14], num_rx[15],
		res[0], res[1], res[2], res[3], res[4], res[5], res[6], res[7], res[8], res[9], res[10], res[11], res[12], res[13], res[14], res[15]
		);
	//gpuErrchk( cudaStreamSynchronize(stream) );  // DEBUG!!!
}

