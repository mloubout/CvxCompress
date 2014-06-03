#include <stdio.h>
#include <cuda_runtime.h>

#define VELMASK 32767
#define EPSMASK 255
#define DELMASK 255
#define C44C33MASK 1
#define QMASK 255
#define DENMASK 255
#define DIPMASK 255
#define AZMMASK 255

#define SHIFTEps 15
#define SHIFTDel 23
#define SHIFTC44C33 31

#define SHIFTQ 24
#define SHIFTDen 16
#define SHIFTAzm  8

__device__ 
float Decode_Index(
	int Index,
	int shift,
	int mask,
	float min,
	float scaler
	)
{
	return (float)((Index >> shift) & mask) * scaler + min;
}

__device__
float Decode_Index_No_Shift(
	int Index,
	int mask,
	float min,
	float scaler
	)
{
	return (float)(Index & mask) * scaler + min;
}

__device__ 
void Decode_DenAng(
	int DenAng,
	float& inv_Q,
	float& Den,
	float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler
	)
{
	inv_Q = Decode_Index(DenAng,SHIFTQ,QMASK,inv_Q_min,inv_Q_scaler);
	Den = Decode_Index(DenAng,SHIFTDen,DENMASK,Den_min,Den_scaler);
}

__device__ 
void Decode_DenAng(
	int DenAng,
	float& inv_Q,
	float& Den,
	float& Dip,
	float& Azm,
	float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
        float Dip_min,
        float Dip_scaler,
        float Azm_min,
        float Azm_scaler
	)
{
	inv_Q = Decode_Index(DenAng,SHIFTQ,QMASK,inv_Q_min,inv_Q_scaler);
	Den = Decode_Index(DenAng,SHIFTDen,DENMASK,Den_min,Den_scaler);
	Dip = Decode_Index_No_Shift(DenAng,DIPMASK,Dip_min,Dip_scaler);
	Azm = Decode_Index(DenAng,SHIFTAzm,AZMMASK,Azm_min,Azm_scaler);
}

__device__ 
void Decode_DenAng_And_VelAnis(
	int DenAng,
	int VelAnis,
	float& inv_Q,
	float& Den,
	float& C44C33,
	float& Vel,
	float& Del,
	float& Eps,
	float inv_Q_min,
	float inv_Q_scaler,
	float Den_min,
	float Den_scaler,
        float C44C33_min,
        float C44C33_scaler,
        float Vel_min,
        float Vel_scaler,
        float Del_min,
        float Del_scaler,
        float Eps_min,
        float Eps_scaler
	)
{
	inv_Q = Decode_Index(DenAng,SHIFTQ,QMASK,inv_Q_min,inv_Q_scaler);
        Den = Decode_Index(DenAng,SHIFTDen,DENMASK,Den_min,Den_scaler);
	C44C33 = Decode_Index(VelAnis,SHIFTC44C33,C44C33MASK,C44C33_min,C44C33_scaler);
	Vel = Decode_Index_No_Shift(VelAnis,VELMASK,Vel_min,Vel_scaler);
	Del = Decode_Index(VelAnis,SHIFTDel,DELMASK,Del_min,Del_scaler);
	Eps = Decode_Index(VelAnis,SHIFTEps,EPSMASK,Eps_min,Eps_scaler);
}

__global__
#if __CUDA_ARCH__ >= 300
__launch_bounds__(1536)
#elif __CUDA_ARCH__ >= 130
__launch_bounds__(960)
#endif
void Compute_DX_DY_Main(
	const float* d_pq0,
	const float* d_pq1,
	const float* d_pq2,
	const float* d_pq3,
	float* d_dx,
	float* d_dy,
	int dimy,
	int dimz,
	const float a1h,
	const float a2h,
	const float a3h,
	const float a4h,
	const float a5h
	)
{
	int stride_z = 16;
	int stride_y = stride_z * dimz;

	const int thr_z = blockIdx.y * 4; // blockIdx.y is intentional, CC < 2.0 support only 2D blocks

	const float2 *pq0;
	if (threadIdx.y == 0)
	{
		pq0 = ((const float2*)d_pq0) + threadIdx.x + thr_z * stride_z;
	}
	else if (threadIdx.y == 1)
	{
		pq0 = ((const float2*)d_pq1) + threadIdx.x + thr_z * stride_z;
	}
	else if (threadIdx.y == 2)
	{
		pq0 = ((const float2*)d_pq1) + threadIdx.x + (thr_z + 2) * stride_z;
	}
	else if (threadIdx.y == 3)
	{
		pq0 = ((const float2*)d_pq2) + threadIdx.x + thr_z * stride_z;
	}
	else if (threadIdx.y == 4)
	{
		pq0 = ((const float2*)d_pq3) + threadIdx.x + thr_z * stride_z;
	}
	else if (threadIdx.y == 5)
	{
		pq0 = ((const float2*)d_pq3) + threadIdx.x + (thr_z + 2) * stride_z;
	}

	float part_dpdy[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
	float part_dqdy[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

	__shared__ float pq[640];

	const int lidx = (threadIdx.x & 15) + (threadIdx.x / 16) * 80;

	const int tid = threadIdx.x + threadIdx.y * blockDim.x;
	const int px = tid % 48;
	const int pz = tid / 48;
	
	float* cp = pq + pz * 80 + px + 8;

	int stride_z_dy = 64;
	int stride_y_dy = stride_z_dy * dimz;

	float2* dx = ((float2*)d_dx) + tid + thr_z * stride_z_dy;
	float2* dy = ((float2*)d_dy) + tid + thr_z * stride_z_dy + 4 * stride_y_dy;

#ifdef NAN_DESU_KA
	int fu1 = 0;
	int fu2 = 0;
#endif

	for (int iY = 0;  iY < dimy;  ++iY)
	{
		if (threadIdx.y == 0)
		{
			float2 v = pq0[0];
			pq[lidx] = v.x;
			pq[320+lidx] = v.y;

			v = pq0[stride_z*2];
			pq[lidx+160] = v.x;
			pq[320+lidx+160] = v.y;
		}
		else if (threadIdx.y == 1)
		{
			float2 v = pq0[0];
			pq[lidx+16] = v.x;
			pq[320+lidx+16] = v.y;
		}
		else if (threadIdx.y == 2)
		{
			float2 v = pq0[0];
			pq[lidx+160+16] = v.x;
			pq[320+lidx+160+16] = v.y;
		}
		else if (threadIdx.y == 3)
		{
			float2 v = pq0[0];
			pq[lidx+32] = v.x;
			pq[320+lidx+32] = v.y;

			v = pq0[stride_z*2];
			pq[lidx+160+32] = v.x;
			pq[320+lidx+160+32] = v.y;
		}
		else if (threadIdx.y == 4)
		{
			float2 v = pq0[0];
			pq[lidx+48] = v.x;
			pq[320+lidx+48] = v.y;
		}
		else if (threadIdx.y == 5)
		{
			float2 v = pq0[0];
			pq[lidx+160+48] = v.x;
			pq[320+lidx+160+48] = v.y;
		}
		
		pq0 += stride_y;
		__syncthreads();

		float dpdx = a5h * cp[5];
		float dqdx = a5h * cp[320+5];
		dpdx += a4h * cp[4];
		dqdx += a4h * cp[320+4];
		dpdx += a3h * cp[3];
		dqdx += a3h * cp[320+3];
		dpdx += a2h * cp[2];
		dqdx += a2h * cp[320+2];
		dpdx += a1h * cp[1];
		dqdx += a1h * cp[320+1];
		dpdx -= a1h * cp[0];
		dqdx -= a1h * cp[320+0];
		dpdx -= a2h * cp[-1];
		dqdx -= a2h * cp[320-1];
		dpdx -= a3h * cp[-2];
		dqdx -= a3h * cp[320-2];
		dpdx -= a4h * cp[-3];
		dqdx -= a4h * cp[320-3];
		dpdx -= a5h * cp[-4];
		dqdx -= a5h * cp[320-4];
		
		float pval = *cp;
		float dpdy = part_dpdy[0] + a5h * pval; 
		part_dpdy[0] = part_dpdy[1] + a4h * pval;
		part_dpdy[1] = part_dpdy[2] + a3h * pval;
		part_dpdy[2] = part_dpdy[3] + a2h * pval;
		part_dpdy[3] = part_dpdy[4] + a1h * pval;
		part_dpdy[4] = part_dpdy[5] - a1h * pval;
		part_dpdy[5] = part_dpdy[6] - a2h * pval;
		part_dpdy[6] = part_dpdy[7] - a3h * pval;
		part_dpdy[7] = part_dpdy[8] - a4h * pval;
		part_dpdy[8] = -a5h * pval;

		float qval = cp[320];
		float dqdy = part_dqdy[0] + a5h * qval;
		part_dqdy[0] = part_dqdy[1] + a4h * qval;
		part_dqdy[1] = part_dqdy[2] + a3h * qval;
		part_dqdy[2] = part_dqdy[3] + a2h * qval;
		part_dqdy[3] = part_dqdy[4] + a1h * qval;
		part_dqdy[4] = part_dqdy[5] - a1h * qval;
		part_dqdy[5] = part_dqdy[6] - a2h * qval;
		part_dqdy[6] = part_dqdy[7] - a3h * qval;
		part_dqdy[7] = part_dqdy[8] - a4h * qval;
		part_dqdy[8] = -a5h * qval;

		float2 v2;
		v2.x = dpdx;
		v2.y = dqdx;
		*dx = v2;
		dx += stride_y_dy;
#ifdef NAN_DESU_KA
		if (!fu1 && (dpdx != 0.0f || dqdx != 0.0f))
		{
			fu1 = 1;
			printf("Compute_DX_DY_Main ; threadIdx.x=[%d,%d,%d] blockIdx=[%d,%d] :: dpdx=%e, dqdx=%e\n",threadIdx.x,threadIdx.y,threadIdx.z,blockIdx.x,blockIdx.y,dpdx,dqdx);
		}
#endif

		if (iY >= 9)
		{
#ifdef NAN_DESU_KA
			if (!fu2 && (dpdy != 0.0f || dqdy != 0.0f))
			{
				fu2 = 1;
				printf("Compute_DX_DY_Main ; threadIdx.x=[%d,%d,%d] blockIdx=[%d,%d] :: dpdy=%e, dqdy=%e\n",threadIdx.x,threadIdx.y,threadIdx.z,blockIdx.x,blockIdx.y,dpdy,dqdy);
			}
#endif
			float2 v3;
			v3.x = dpdy;
			v3.y = dqdy;
			*dy = v3;
			dy += stride_y_dy;
		}

		__syncthreads();
	}
}

__global__
#if __CUDA_ARCH__ >= 300
__launch_bounds__(2048)
#elif __CUDA_ARCH__ >= 130
__launch_bounds__(1344)
#endif
void Compute_DZ_Main(
	const float* d_pq0,
        const float* d_pq1,
        const float* d_pq2,
        const float* d_pq3,
	float* d_dz,
	int dimy,
	int dimz,
	const float a1,
	const float a2,
	const float a3,
	const float a4,
	const float a5,
	const float e1
	)
{
	const int stride_z = 32;
	const int stride_y = stride_z * dimz;

	const int thr_y = blockIdx.y + 4; 

	const float *pq0;
	int lidx;
	if (threadIdx.y == 0)
	{
		pq0 = d_pq0 + threadIdx.x + thr_y * stride_y;
		lidx = (threadIdx.x & 1) * 320 + (threadIdx.x / 2);
	}
	else if (threadIdx.y == 1)
	{
		pq0 = d_pq1 + threadIdx.x + thr_y * stride_y;
		lidx = (threadIdx.x & 1) * 320 + (threadIdx.x / 2) + 16;
	}
	else if (threadIdx.y == 2)
	{
		pq0 = d_pq1 + threadIdx.x + thr_y * stride_y + 2 * stride_z;
		lidx = (threadIdx.x & 1) * 320 + (threadIdx.x / 2) + 16 + 80 * 2;
	}
	else if (threadIdx.y == 3)
	{
		pq0 = d_pq2 + threadIdx.x + thr_y * stride_y;
		lidx = (threadIdx.x & 1) * 320 + (threadIdx.x / 2) + 32;
	}
	else if (threadIdx.y == 4)
	{
		pq0 = d_pq3 + threadIdx.x + thr_y * stride_y;
		lidx = (threadIdx.x & 1) * 320 + (threadIdx.x / 2) + 48;
	}
	else if (threadIdx.y == 5)
	{
		pq0 = d_pq3 + threadIdx.x + thr_y * stride_y + 2 * stride_z;
		lidx = (threadIdx.x & 1) * 320 + (threadIdx.x / 2) + 48 + 80 * 2;
	}

	float part_dpdz[3] = {0.0f, 0.0f, 0.0f};
	float part_dqdz[3] = {0.0f, 0.0f, 0.0f};

	__shared__ float pq[768];

	const int tid = threadIdx.x + threadIdx.y * blockDim.x;
	const int px = tid % 48;
	const int pz = tid / 48;

	const int stride_z_dz = 48;
	const int stride_y_dz = stride_z_dz * dimz;

	float2* dz = (float2*)d_dz + thr_y * stride_y_dz;
	
	float* cp = pq + px + 8;

	// DZ lead-in - Z=0
	if (threadIdx.y != 2 && threadIdx.y != 5) pq[lidx] = *pq0;
	pq0 += stride_z;
	__syncthreads();

	float pval = cp[0];
	float qval=  cp[320];
	if (pz == 0)
	{
		part_dpdz[0] = -a1 * pval;
		part_dqdz[0] = -a1 * qval;
		
		part_dpdz[1] = -a5 * pval;
		part_dqdz[1] = -a5 * qval;
	}
	else if (pz == 1)
	{
		part_dpdz[0] = -a2 * pval;
		part_dqdz[0] = -a2 * qval;		
	}
	else if (pz == 2)
	{
		part_dpdz[0] = -a3 * pval;
		part_dqdz[0] = -a3 * qval;
	}
	else if (pz == 3)
	{
		part_dpdz[0] = -a4 * pval;
		part_dqdz[0] = -a4 * qval;
	}
	__syncthreads();

	// DZ lead-in - Z=[1,4]
	if (threadIdx.y == 0 || threadIdx.y == 3)
	{
		pq[lidx] = pq0[0];
		pq[lidx+80] = pq0[stride_z];
		pq[lidx+160] = pq0[2*stride_z];
		pq[lidx+240] = pq0[3*stride_z];
	}
	else
	{
		pq[lidx] = pq0[0];
		pq[lidx+80] = pq0[stride_z];
	}
	pq0 += 4*stride_z;
	__syncthreads();

	if (pz == 0)
	{
		part_dpdz[0] += (a2 + a1) * cp[  0] + (a3 + a2) * cp[    80] + (a4 + a3) * cp[    160] + (a5 + a4) * cp[    240];
		part_dqdz[0] += (a2 + a1) * cp[320] + (a3 + a2) * cp[320+80] + (a4 + a3) * cp[320+160] + (a5 + a4) * cp[320+240];

		part_dpdz[1] +=    -a4    * cp[  0]      -a3    * cp[    80]      -a2    * cp[    160]      -a1    * cp[    240];
		part_dqdz[1] +=    -a4    * cp[320]      -a3    * cp[320+80]      -a2    * cp[320+160]      -a1    * cp[320+240];

		part_dpdz[2]  =                                                                             -a5    * cp[    240];
		part_dqdz[2]  =                                                                             -a5    * cp[320+240];
	}
	else if (pz == 1)
	{
		part_dpdz[0] += (a3 - a1) * cp[  0] + (a4 + a1) * cp[    80] + (a5 + a2) * cp[    160]      +a3    * cp[    240];
		part_dqdz[0] += (a3 - a1) * cp[320] + (a4 + a1) * cp[320+80] + (a5 + a2) * cp[320+160]      +a3    * cp[320+240];
		
		part_dpdz[1] +=    -a5    * cp[  0]      -a4    * cp[    80]      -a3    * cp[    160]      -a2    * cp[    240];
		part_dqdz[1] +=    -a5    * cp[320]      -a4    * cp[320+80]      -a3    * cp[320+160]      -a2    * cp[320+240];
	}
	else if (pz == 2)
	{
		part_dpdz[0] += (a4 - a2) * cp[  0] + (a5 - a1) * cp[    80]      +a1    * cp[    160]      +a2    * cp[    240];
		part_dqdz[0] += (a4 - a2) * cp[320] + (a5 - a1) * cp[320+80]      +a1    * cp[320+160]      +a2    * cp[320+240];

		part_dpdz[1] +=                          -a5    * cp[    80]      -a4    * cp[    160]      -a3    * cp[    240];
		part_dqdz[1] +=                          -a5    * cp[320+80]      -a4    * cp[320+160]      -a3    * cp[320+240];
	}
	else if (pz == 3)
	{
		part_dpdz[0] += (a5 - a3) * cp[  0]      -a2    * cp[    80]      -a1    * cp[    160]      +a1    * cp[    240];
		part_dqdz[0] += (a5 - a3) * cp[320]      -a2    * cp[320+80]      -a1    * cp[320+160]      +a1    * cp[320+240];

		part_dpdz[1] +=                                                   -a5    * cp[    160]      -a4    * cp[    240];
		part_dqdz[1] +=                                                   -a5    * cp[320+160]      -a4    * cp[320+240];
	}
	__syncthreads();

	// regular loop
	for (int iZ = 0;  iZ < dimz;  iZ+=4)
	{
		if (iZ < dimz-8)
		{
			if (threadIdx.y == 0 || threadIdx.y == 3)
			{
				pq[lidx    ] = pq0[         0];
				pq[lidx+ 80] = pq0[  stride_z];
				pq[lidx+160] = pq0[2*stride_z];
				pq[lidx+240] = pq0[3*stride_z];
			}
			else
			{
				pq[lidx   ] = pq0[       0];
				pq[lidx+80] = pq0[stride_z];
			}
		}
		else if (iZ < dimz-4)
		{
			if (threadIdx.y == 0 || threadIdx.y == 3)
			{
				pq[lidx    ] = pq0[         0];
				pq[lidx+ 80] = pq0[  stride_z];
				pq[lidx+160] = pq0[2*stride_z];
			}
			else if (threadIdx.y == 1 || threadIdx.y == 4)
			{
				pq[lidx   ] = pq0[       0];
				pq[lidx+80] = pq0[stride_z];
			}
			else
			{
				pq[lidx   ] = pq0[       0];
			}
		}
		pq0 += 4*stride_z;
		__syncthreads();

		float dpdz, dqdz;
		if (iZ < dimz-12)					// four full stencils
		{
			if (pz == 0)
			{
				dpdz         = part_dpdz[0] + a5 * cp[  0];
				part_dpdz[0] = part_dpdz[1] + a1 * cp[  0] + a2 * cp[    80] + a3 * cp[    160] + a4 * cp[    240];
				part_dpdz[1] = part_dpdz[2] - a4 * cp[  0] - a3 * cp[    80] - a2 * cp[    160] - a1 * cp[    240];
				part_dpdz[2] =                                                                  - a5 * cp[    240];

				dqdz         = part_dqdz[0] + a5 * cp[320];
				part_dqdz[0] = part_dqdz[1] + a1 * cp[320] + a2 * cp[320+80] + a3 * cp[320+160] + a4 * cp[320+240];
				part_dqdz[1] = part_dqdz[2] - a4 * cp[320] - a3 * cp[320+80] - a2 * cp[320+160] - a1 * cp[320+240];
				part_dqdz[2] =                                                                  - a5 * cp[320+240];
			}
			else if (pz == 1)
			{
				dpdz         = part_dpdz[0] + a4 * cp[  0] + a5 * cp[    80];
				part_dpdz[0] = part_dpdz[1] - a1 * cp[  0] + a1 * cp[    80] + a2 * cp[    160] + a3 * cp[    240];
				part_dpdz[1] =              - a5 * cp[  0] - a4 * cp[    80] - a3 * cp[    160] - a2 * cp[    240];

				dqdz         = part_dqdz[0] + a4 * cp[320] + a5 * cp[320+80];
				part_dqdz[0] = part_dqdz[1] - a1 * cp[320] + a1 * cp[320+80] + a2 * cp[320+160] + a3 * cp[320+240];
				part_dqdz[1] =              - a5 * cp[320] - a4 * cp[320+80] - a3 * cp[320+160] - a2 * cp[320+240];
			}
			else if (pz == 2)
			{
				dpdz         = part_dpdz[0] + a3 * cp[  0] + a4 * cp[    80] + a5 * cp[    160];
				part_dpdz[0] = part_dpdz[1] - a2 * cp[  0] - a1 * cp[    80] + a1 * cp[    160] + a2 * cp[    240];
				part_dpdz[1] =                             - a5 * cp[    80] - a4 * cp[    160] - a3 * cp[    240];

				dqdz         = part_dqdz[0] + a3 * cp[320] + a4 * cp[320+80] + a5 * cp[320+160];
				part_dqdz[0] = part_dqdz[1] - a2 * cp[320] - a1 * cp[320+80] + a1 * cp[320+160] + a2 * cp[320+240];
				part_dqdz[1] =                             - a5 * cp[320+80] - a4 * cp[320+160] - a3 * cp[320+240];
			}
			else if (pz == 3)
			{
				dpdz         = part_dpdz[0] + a2 * cp[  0] + a3 * cp[    80] + a4 * cp[    160] + a5 * cp[    240];
				part_dpdz[0] = part_dpdz[1] - a3 * cp[  0] - a2 * cp[    80] - a1 * cp[    160] + a1 * cp[    240];
				part_dpdz[1] =                                               - a5 * cp[    160] - a4 * cp[    240];

				dqdz         = part_dqdz[0] + a2 * cp[320] + a3 * cp[320+80] + a4 * cp[320+160] + a5 * cp[320+240];
				part_dqdz[0] = part_dqdz[1] - a3 * cp[320] - a2 * cp[320+80] - a1 * cp[320+160] + a1 * cp[320+240];
				part_dqdz[1] =                                               - a5 * cp[320+160] - a4 * cp[320+240];
			}
		}
		else if (iZ < dimz-8)					// four full stencils, but switch pz==3 to end pt stencil
		{
			if (pz == 0)
			{
				dpdz         = part_dpdz[0] + a5 * cp[  0];
				part_dpdz[0] = part_dpdz[1] + a1 * cp[  0] + a2 * cp[    80] + a3 * cp[    160] + a4 * cp[    240];
				part_dpdz[1] =                                                                  - e1 * cp[    240];

				dqdz         = part_dqdz[0] + a5 * cp[320];
				part_dqdz[0] = part_dqdz[1] + a1 * cp[320] + a2 * cp[320+80] + a3 * cp[320+160] + a4 * cp[320+240];
				part_dqdz[1] =                                                                  - e1 * cp[320+240];
			}
			else if (pz == 1)
			{
				dpdz         = part_dpdz[0] + a4 * cp[  0] + a5 * cp[    80];
				part_dpdz[0] = part_dpdz[1] - a1 * cp[  0] + a1 * cp[    80] + a2 * cp[    160] + a3 * cp[    240];

				dqdz         = part_dqdz[0] + a4 * cp[320] + a5 * cp[320+80];
				part_dqdz[0] = part_dqdz[1] - a1 * cp[320] + a1 * cp[320+80] + a2 * cp[320+160] + a3 * cp[320+240];
			}
			else if (pz == 2)
			{
				dpdz         = part_dpdz[0] + a3 * cp[  0] + a4 * cp[    80] + a5 * cp[    160];
				part_dpdz[0] = part_dpdz[1] - a2 * cp[  0] - a1 * cp[    80] + a1 * cp[    160] + a2 * cp[    240];

				dqdz         = part_dqdz[0] + a3 * cp[320] + a4 * cp[320+80] + a5 * cp[320+160];
				part_dqdz[0] = part_dqdz[1] - a2 * cp[320] - a1 * cp[320+80] + a1 * cp[320+160] + a2 * cp[320+240];
			}
			else if (pz == 3)
			{
				dpdz         = part_dpdz[0] + a2 * cp[  0] + a3 * cp[    80] + a4 * cp[    160] + a5 * cp[    240];
				part_dpdz[0] =                                               - e1 * cp[    160] + e1 * cp[    240];

				dqdz         = part_dqdz[0] + a2 * cp[320] + a3 * cp[320+80] + a4 * cp[320+160] + a5 * cp[320+240];
				part_dqdz[0] =                                               - e1 * cp[320+160] + e1 * cp[320+240];
			}
		}
		else if (iZ < dimz-4)					// three full and one end pt stencils
		{
			if (pz == 0)
			{
				dpdz         = part_dpdz[0] + a5 * cp[  0];
				dqdz         = part_dqdz[0] + a5 * cp[320];
			}
			else if (pz == 1)
			{
				dpdz         = part_dpdz[0] + a4 * cp[  0] + a5 * cp[    80];
				dqdz         = part_dqdz[0] + a4 * cp[320] + a5 * cp[320+80];
			}
			else if (pz == 2)
			{
				dpdz         = part_dpdz[0] + a3 * cp[  0] + a4 * cp[    80] + a5 * cp[    160];
				dqdz         = part_dqdz[0] + a3 * cp[320] + a4 * cp[320+80] + a5 * cp[320+160];
			}
			else if (pz == 3)
			{
				dpdz         = part_dpdz[0];
				dqdz         = part_dqdz[0];
			}
		}
		else               					// three end pt stencils and one zero
		{
			if (pz == 0)
			{
				dpdz = part_dpdz[1] + e1 * cp[  0];
				dqdz = part_dqdz[1] + e1 * cp[320];
			}
			else if (pz == 1)
			{
				dpdz =              - e1 * cp[  0] + e1 * cp[    80];
				dqdz =              - e1 * cp[320] + e1 * cp[320+80];
			}
			else if (pz == 2)
			{
				dpdz =                             - e1 * cp[    80] + e1 * cp[    160];
				dqdz =                             - e1 * cp[320+80] + e1 * cp[320+160];
			}
			else if (pz == 3)
			{
				dpdz = 0.0f;
				dqdz = 0.0f;
			}
		}

		float2 dz_val;
		dz_val.x = dpdz;
		dz_val.y = dqdz;
		dz[tid] = dz_val;
		dz += 4 * stride_z_dz;

		__syncthreads();
	}
}

__device__ 
void Compute_DXED(
	int iDenAng,
	float dpdx,
	float dqdx,
	float dpdy,
	float dqdy,
	float dpdz,
	float dqdz,
        float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
        float Dip_min,
        float Dip_scaler,
        float Azm_min,
        float Azm_scaler,
	float& dxed1_p,
	float& dxed1_q,
	float& dxed2_p,
	float& dxed2_q
	)
{
	float inv_Q, Den, Dip, Azm;
	Decode_DenAng(iDenAng,inv_Q,Den,Dip,Azm,inv_Q_min,inv_Q_scaler,Den_min,Den_scaler,Dip_min,Dip_scaler,Azm_min,Azm_scaler);
	float buoy = 1.0f / Den;

	float cAzm, sAzm, cDip, sDip;
	__sincosf(Dip,&sDip,&cDip);
	__sincosf(Azm,&sAzm,&cAzm);

	// compute dxed1, dxed2, dyed1, dyed2, dzed1 and dzed2
	float temp_p = cAzm * dpdx + sAzm * dpdy;
	float temp2_p = buoy * ( cAzm * dpdy - sAzm * dpdx );
	float temp3_p = buoy * ( cDip * temp_p - sDip * dpdz );
	float temp4_p = buoy * ( sDip * temp_p + cDip * dpdz );

	dxed1_p = -sAzm * temp2_p + cDip * cAzm * temp3_p;
	dxed2_p = sDip * cAzm * temp4_p;

	float temp_q = cAzm * dqdx + sAzm * dqdy;
	float temp2_q = buoy * ( cAzm * dqdy - sAzm * dqdx );
	float temp3_q = buoy * ( cDip * temp_q - sDip * dqdz );
	float temp4_q = buoy * ( sDip * temp_q + cDip * dqdz );

	dxed1_q = -sAzm * temp2_q + cDip * cAzm * temp3_q;
	dxed2_q = sDip * cAzm * temp4_q;
}

__device__ 
void Compute_DYED(
	int iDenAng,
	float dpdx,
	float dqdx,
	float dpdy,
	float dqdy,
	float dpdz,
	float dqdz,
        float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
        float Dip_min,
        float Dip_scaler,
        float Azm_min,
        float Azm_scaler,
	float& dyed1_p,
	float& dyed1_q,
	float& dyed2_p,
	float& dyed2_q
	)
{
	float inv_Q, Den, Dip, Azm;
	Decode_DenAng(iDenAng,inv_Q,Den,Dip,Azm,inv_Q_min,inv_Q_scaler,Den_min,Den_scaler,Dip_min,Dip_scaler,Azm_min,Azm_scaler);
	float buoy = 1.0f / Den;

	float cAzm, sAzm, cDip, sDip;
	__sincosf(Dip,&sDip,&cDip);
	__sincosf(Azm,&sAzm,&cAzm);

	// compute dxed1, dxed2, dyed1, dyed2, dzed1 and dzed2
	float temp_p = cAzm * dpdx + sAzm * dpdy;
	float temp2_p = buoy * ( cAzm * dpdy - sAzm * dpdx );
	float temp3_p = buoy * ( cDip * temp_p - sDip * dpdz );
	float temp4_p = buoy * ( sDip * temp_p + cDip * dpdz );

	dyed1_p = cAzm * temp2_p + cDip * sAzm * temp3_p;
	dyed2_p = sDip * sAzm * temp4_p;

	float temp_q = cAzm * dqdx + sAzm * dqdy;
	float temp2_q = buoy * ( cAzm * dqdy - sAzm * dqdx );
	float temp3_q = buoy * ( cDip * temp_q - sDip * dqdz );
	float temp4_q = buoy * ( sDip * temp_q + cDip * dqdz );

	dyed1_q = cAzm * temp2_q + cDip * sAzm * temp3_q;
	dyed2_q = sDip * sAzm * temp4_q;
}

__device__ 
void Compute_DXED_DYED_DZED(
	int iDenAng,
	float dpdx,
	float dqdx,
	float dpdy,
	float dqdy,
	float dpdz,
	float dqdz,
        float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
        float Dip_min,
        float Dip_scaler,
        float Azm_min,
        float Azm_scaler,
	float& dxed1_p,
	float& dxed1_q,
	float& dxed2_p,
	float& dxed2_q,
	float& dyed1_p,
	float& dyed1_q,
	float& dyed2_p,
	float& dyed2_q,
	float& dzed1_p,
	float& dzed1_q,
	float& dzed2_p,
	float& dzed2_q
	)
{

	float inv_Q, Den, Dip, Azm;
	Decode_DenAng(iDenAng,inv_Q,Den,Dip,Azm,inv_Q_min,inv_Q_scaler,Den_min,Den_scaler,Dip_min,Dip_scaler,Azm_min,Azm_scaler);
	float buoy = 1.0f / Den;

	float cAzm, sAzm, cDip, sDip;
	__sincosf(Dip,&sDip,&cDip);
	__sincosf(Azm,&sAzm,&cAzm);

	// compute dxed1, dxed2, dyed1, dyed2, dzed1 and dzed2
	float temp_p = cAzm * dpdx + sAzm * dpdy;
	float temp2_p = buoy * ( cAzm * dpdy - sAzm * dpdx );
	float temp3_p = buoy * ( cDip * temp_p - sDip * dpdz );
	float temp4_p = buoy * ( sDip * temp_p + cDip * dpdz );

	dxed1_p = -sAzm * temp2_p + cDip * cAzm * temp3_p;
	dyed1_p = cAzm * temp2_p + cDip * sAzm * temp3_p;
	dzed1_p = sDip * temp3_p;

	dxed2_p = sDip * cAzm * temp4_p;
	dyed2_p = sDip * sAzm * temp4_p;
	dzed2_p = cDip * temp4_p;

	float temp_q = cAzm * dqdx + sAzm * dqdy;
	float temp2_q = buoy * ( cAzm * dqdy - sAzm * dqdx );
	float temp3_q = buoy * ( cDip * temp_q - sDip * dqdz );
	float temp4_q = buoy * ( sDip * temp_q + cDip * dqdz );

	dxed1_q = -sAzm * temp2_q + cDip * cAzm * temp3_q;
	dyed1_q = cAzm * temp2_q + cDip * sAzm * temp3_q;
	dzed1_q = sDip * temp3_q;

	dxed2_q = sDip * cAzm * temp4_q;
	dyed2_q = sDip * sAzm * temp4_q;
	dzed2_q = cDip * temp4_q;
}

__device__ 
void Compute_VTI_DXED(
	int iDenAng,
	float dpdx,
	float dqdx,
	float dpdy,
	float dqdy,
	float dpdz,
	float dqdz,
        float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
	float& dxed1_p,
	float& dxed1_q
	)
{

	float inv_Q, Den;
	Decode_DenAng(iDenAng,inv_Q,Den,inv_Q_min,inv_Q_scaler,Den_min,Den_scaler);
	float buoy = 1.0f / Den;

	// compute dxed1, dxed2, dyed1, dyed2, dzed1 and dzed2
	dxed1_p = buoy * dpdx;
	dxed1_q = buoy * dqdx;
}

__device__ 
void Compute_VTI_DYED(
	int iDenAng,
	float dpdx,
	float dqdx,
	float dpdy,
	float dqdy,
	float dpdz,
	float dqdz,
        float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
	float& dyed1_p,
	float& dyed1_q
	)
{

	float inv_Q, Den;
	Decode_DenAng(iDenAng,inv_Q,Den,inv_Q_min,inv_Q_scaler,Den_min,Den_scaler);
	float buoy = 1.0f / Den;

	// compute dxed1, dxed2, dyed1, dyed2, dzed1 and dzed2
	dyed1_p = buoy * dpdy;
	dyed1_q = buoy * dqdy;
}

__device__ 
void Compute_VTI_DXED_DYED_DZED(
	int iDenAng,
	float dpdx,
	float dqdx,
	float dpdy,
	float dqdy,
	float dpdz,
	float dqdz,
        float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
	float& dxed1_p,
	float& dxed1_q,
	float& dyed1_p,
	float& dyed1_q,
	float& dzed2_p,
	float& dzed2_q
	)
{

	float inv_Q, Den;
	Decode_DenAng(iDenAng,inv_Q,Den,inv_Q_min,inv_Q_scaler,Den_min,Den_scaler);
	float buoy = 1.0f / Den;

	// compute dxed1, dxed2, dyed1, dyed2, dzed1 and dzed2
	dxed1_p = buoy * dpdx;
	dyed1_p = buoy * dpdy;
	dzed2_p = buoy * dpdz;

	dxed1_q = buoy * dqdx;
	dyed1_q = buoy * dqdy;
	dzed2_q = buoy * dqdz;
}

__device__ 
float Compute_DX_From_Shared(
	float* v,
	float a1h,
	float a2h,
	float a3h,
	float a4h,
	float a5h
	)
{
	// stencil is anti-symmetric, but there is no saving in exploiting this on architectures with MAD instructions.
	return 
		a5h * v[ 5] + 
		a4h * v[ 4] + 
		a3h * v[ 3] + 
		a2h * v[ 2] + 
		a1h * v[ 1] -
		a1h * v[ 0] -
		a2h * v[-1] - 
		a3h * v[-2] - 
		a4h * v[-3] -
		a5h * v[-4];
}

__device__ 
void Compute_VTI_DYED1and2_Part_V4_V5_YHalo(
	float* dx_dy_dz_DenAng,
	float* d_dx,
        float* d_dy,
        float* d_dz,
        const int* d_em0,
        const int* d_em1,
        const int* d_em2,
        const int* d_em3,
        int dimz,
        const float dt,
        const float a1,
        const float a2,
        const float a3,
        const float a4,
        const float a5,
        const float e1,
        const float a1h,
        const float a2h,
        const float a3h,
        const float a4h,
        const float a5h,
        float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
        float C44C33_min,
        float C44C33_scaler,
        float Vel_min,
        float Vel_scaler,
        float Del_min,
        float Del_scaler,
        float Eps_min,
        float Eps_scaler
	)
{
	int stride_z_dy = 64;
	int stride_y_dy = stride_z_dy * dimz;

	int stride_z_dz = 48;
	int stride_y_dz = stride_z_dz * dimz;

	int stride_z_em = 16;
	int stride_y_em = stride_z_em * dimz;

	int thr_y = blockIdx.y + 4;

	//
	// block shape is 32 by 4 by 1 (x by y by z)
	// algorithm loops over z.
	//

#ifdef NAN_DESU_KA
	int fu2 = 0;
#endif

	for (int iZ = 0;  iZ < dimz;  iZ+=4)
	{
		if (threadIdx.y == 0)
		{
			float2* dx = (float2*)d_dx + thr_y * stride_y_dy + iZ * stride_z_dy;
			const int* em = d_em0 + thr_y * 2 * stride_y_em + iZ * stride_z_em;

			float2 v = dx[threadIdx.x];
			dx_dy_dz_DenAng[    threadIdx.x] = v.x;
			dx_dy_dz_DenAng[192+threadIdx.x] = v.y;

			v = dx[32+threadIdx.x];
			dx_dy_dz_DenAng[ 32+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[224+threadIdx.x] = v.y;

			v = dx[64+threadIdx.x];
			dx_dy_dz_DenAng[ 64+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[256+threadIdx.x] = v.y;

			v = dx[96+threadIdx.x];
			dx_dy_dz_DenAng[ 96+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[288+threadIdx.x] = v.y;

			v = dx[128+threadIdx.x];
			dx_dy_dz_DenAng[128+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[320+threadIdx.x] = v.y;

			v = dx[160+threadIdx.x];
			dx_dy_dz_DenAng[160+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[352+threadIdx.x] = v.y;

			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		else if (threadIdx.y == 1)
		{
			float2* dx = (float2*)d_dy + thr_y * stride_y_dy + iZ * stride_z_dy;
			const int* em = d_em1 + thr_y * 2 * stride_y_em + iZ * stride_z_em;

			float2 v = dx[threadIdx.x];
			dx_dy_dz_DenAng[384+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[576+threadIdx.x] = v.y;

			v = dx[32+threadIdx.x];
			dx_dy_dz_DenAng[416+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[608+threadIdx.x] = v.y;

			v = dx[64+threadIdx.x];
			dx_dy_dz_DenAng[448+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[640+threadIdx.x] = v.y;

			v = dx[96+threadIdx.x];
			dx_dy_dz_DenAng[480+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[672+threadIdx.x] = v.y;

			v = dx[128+threadIdx.x];
			dx_dy_dz_DenAng[512+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[704+threadIdx.x] = v.y;

			v = dx[160+threadIdx.x];
			dx_dy_dz_DenAng[544+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[736+threadIdx.x] = v.y;

			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		else if (threadIdx.y == 2)
		{
			float2* dx = (float2*)d_dz + thr_y * stride_y_dz + iZ * stride_z_dz;
			const int* em = d_em2 + thr_y * 2 * stride_y_em + iZ * stride_z_em;

			float2 v = dx[threadIdx.x];
			dx_dy_dz_DenAng[ 768+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[ 960+threadIdx.x] = v.y;

			v = dx[32+threadIdx.x];
			dx_dy_dz_DenAng[ 800+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[ 992+threadIdx.x] = v.y;

			v = dx[64+threadIdx.x];
			dx_dy_dz_DenAng[ 832+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1024+threadIdx.x] = v.y;

			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		else if (threadIdx.y == 3)
		{
			float2* dx = (float2*)d_dz + thr_y * stride_y_dz + iZ * stride_z_dz;
			const int* em = d_em3 + thr_y * 2 * stride_y_em + iZ * stride_z_em;

			float2 v = dx[96+threadIdx.x];
			dx_dy_dz_DenAng[ 864+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1056+threadIdx.x] = v.y;

			v = dx[128+threadIdx.x];
			dx_dy_dz_DenAng[ 896+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1088+threadIdx.x] = v.y;

			v = dx[160+threadIdx.x];
			dx_dy_dz_DenAng[ 928+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1120+threadIdx.x] = v.y;

			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		__syncthreads();

		// dx is in dx_dy_dz_DenAng[  0, 384>, [0,96> contains 1st line, [96,192> contains 2nd line etc.
		// dy is in dx_dy_dz_DenAng[384, 768>
		// dz is in dx_dy_dz_DenAng[768,1152>

		// compute dxed1, dxed2, dyed1, dyed2, dzed1 and dzed2 for center points first.

		float dyed1_p, dyed1_q;
		Compute_VTI_DYED(
				((int*)dx_dy_dz_DenAng)[1152+threadIdx.y*64+threadIdx.x+16],
				dx_dy_dz_DenAng[    threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[192+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[384+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[576+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[768+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[960+threadIdx.y*48+threadIdx.x+8],
				inv_Q_min,inv_Q_scaler,
				Den_min,Den_scaler,
				dyed1_p, dyed1_q
			    );

		// write out dyed1 and dyed2
		float2 f2_dyed1;
		f2_dyed1.x = dyed1_p;
		f2_dyed1.y = dyed1_q;
		float2* dx = ((float2*)d_dx) + thr_y * stride_y_dy + iZ * stride_z_dy;
		dx[threadIdx.y*64+threadIdx.x] = f2_dyed1;
	
		__syncthreads();
	}
}

__device__
void Compute_VTI_DYED1and2_Part_V4_V5_Main(
	float* dx_dy_dz_DenAng,
	float* d_dx,
	float* d_dy,
	float* d_dz,
	const int* d_em0,
	const int* d_em1,
	const int* d_em2,
	const int* d_em3,
	int dimz,
	const float dt,
        const float a1,
        const float a2,
        const float a3,
        const float a4,
        const float a5,
        const float e1,
        const float a1h,
        const float a2h,
        const float a3h,
        const float a4h,
        const float a5h,
        float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
        float C44C33_min,
        float C44C33_scaler,
        float Vel_min,
        float Vel_scaler,
        float Del_min,
        float Del_scaler,
        float Eps_min,
        float Eps_scaler
	)
{
	const int stride_z_dy = 64;
	const int stride_y_dy = stride_z_dy * dimz;

	const int stride_z_dz = 48;
	const int stride_y_dz = stride_z_dz * dimz;

	const int stride_z_em = 16;
	const int stride_y_em = stride_z_em * dimz;

	const int thr_y = blockIdx.y + 4;

	//
	// block shape is 32 by 4 by 1 (x by y by z)
	// algorithm loops over z.
	//

	float part_ddzed2dz_p[2] = {0.0f, 0.0f};
	float part_ddzed2dz_q[2] = {0.0f, 0.0f};

	float next_d_dxed1_dx_p=0.0f, next_d_dxed1_dx_q=0.0f;
	for (int iZ = 0;  iZ < dimz;  iZ+=4)
	{
		if (threadIdx.y == 0)
		{
			float2* dx = (float2*)d_dx + thr_y * stride_y_dy + iZ * stride_z_dy;
			const int* em = d_em0 + thr_y * 2 * stride_y_em + iZ * stride_z_em;
			
			float2 v = dx[threadIdx.x];
			dx_dy_dz_DenAng[    threadIdx.x] = v.x;
			dx_dy_dz_DenAng[192+threadIdx.x] = v.y;
			
			v = dx[32+threadIdx.x];
			dx_dy_dz_DenAng[ 32+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[224+threadIdx.x] = v.y;

			v = dx[64+threadIdx.x];
			dx_dy_dz_DenAng[ 64+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[256+threadIdx.x] = v.y;

			v = dx[96+threadIdx.x];
			dx_dy_dz_DenAng[ 96+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[288+threadIdx.x] = v.y;

			v = dx[128+threadIdx.x];
			dx_dy_dz_DenAng[128+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[320+threadIdx.x] = v.y;
			
			v = dx[160+threadIdx.x];
			dx_dy_dz_DenAng[160+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[352+threadIdx.x] = v.y;
		
			// NB!!!!
			// Without typecast, input value from em ptr will be converted to float before storage.
			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		else if (threadIdx.y == 1)
		{
			float2* dx = (float2*)d_dy + thr_y * stride_y_dy + iZ * stride_z_dy;
			const int* em = d_em1 + thr_y * 2 * stride_y_em + iZ * stride_z_em;

			float2 v = dx[threadIdx.x];
			dx_dy_dz_DenAng[384+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[576+threadIdx.x] = v.y;
			
			v = dx[32+threadIdx.x];
			dx_dy_dz_DenAng[416+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[608+threadIdx.x] = v.y;

			v = dx[64+threadIdx.x];
			dx_dy_dz_DenAng[448+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[640+threadIdx.x] = v.y;

			v = dx[96+threadIdx.x];
			dx_dy_dz_DenAng[480+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[672+threadIdx.x] = v.y;

			v = dx[128+threadIdx.x];
			dx_dy_dz_DenAng[512+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[704+threadIdx.x] = v.y;
			
			v = dx[160+threadIdx.x];
			dx_dy_dz_DenAng[544+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[736+threadIdx.x] = v.y;
		
			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		else if (threadIdx.y == 2)
		{
			float2* dx = (float2*)d_dz + thr_y * stride_y_dz + iZ * stride_z_dz;
			const int* em = d_em2 + thr_y * 2 * stride_y_em + iZ * stride_z_em;

			float2 v = dx[threadIdx.x];
			dx_dy_dz_DenAng[ 768+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[ 960+threadIdx.x] = v.y;
			
			v = dx[32+threadIdx.x];
			dx_dy_dz_DenAng[ 800+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[ 992+threadIdx.x] = v.y;

			v = dx[64+threadIdx.x];
			dx_dy_dz_DenAng[ 832+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1024+threadIdx.x] = v.y;

			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		else if (threadIdx.y == 3)
		{
			float2* dx = (float2*)d_dz + thr_y * stride_y_dz + iZ * stride_z_dz;
			const int* em = d_em3 + thr_y * 2 * stride_y_em + iZ * stride_z_em;

			float2 v = dx[96+threadIdx.x];
			dx_dy_dz_DenAng[ 864+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1056+threadIdx.x] = v.y;

			v = dx[128+threadIdx.x];
			dx_dy_dz_DenAng[ 896+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1088+threadIdx.x] = v.y;
			
			v = dx[160+threadIdx.x];
			dx_dy_dz_DenAng[ 928+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1120+threadIdx.x] = v.y;
			
			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		__syncthreads();

		// dx is in dx_dy_dz_DenAng[  0, 384>, [0,96> contains 1st line, [96,192> contains 2nd line etc.
		// dy is in dx_dy_dz_DenAng[384, 768>
		// dz is in dx_dy_dz_DenAng[768,1152>

		// compute dxed1, dxed2, dyed1, dyed2, dzed1 and dzed2 for center points first.
		float dyed1_p, dyed1_q;
		float dzed2_p, dzed2_q;
		Compute_VTI_DXED_DYED_DZED(
				((int*)dx_dy_dz_DenAng)[1152+threadIdx.y*64+threadIdx.x+16],
				dx_dy_dz_DenAng[    threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[192+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[384+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[576+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[768+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[960+threadIdx.y*48+threadIdx.x+8],
				inv_Q_min,inv_Q_scaler,
				Den_min,Den_scaler,
				dx_dy_dz_DenAng[    threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[192+threadIdx.y*48+threadIdx.x+8],
				dyed1_p, dyed1_q, 
				dzed2_p, dzed2_q
				);

		// write out dyed1 and dyed2
		float2 f2_dyed1;
		f2_dyed1.x = dyed1_p;
		f2_dyed1.y = dyed1_q;
		float2* dx = ((float2*)d_dx) + thr_y * stride_y_dy + iZ * stride_z_dy;
		dx[threadIdx.y*64+threadIdx.x] = f2_dyed1;

		// compute dxed1 and dxed2 for halo cells
		if (threadIdx.y == 0)
		{
			// left halos
			int thr_x = threadIdx.x & 7;
			int thr_y = threadIdx.x & 24;
			Compute_VTI_DXED(
					((int*)dx_dy_dz_DenAng)[1152+thr_y*8+thr_x+8],
					dx_dy_dz_DenAng[    thr_y*6+thr_x],
					dx_dy_dz_DenAng[192+thr_y*6+thr_x],
					dx_dy_dz_DenAng[384+thr_y*6+thr_x],
					dx_dy_dz_DenAng[576+thr_y*6+thr_x],
					dx_dy_dz_DenAng[768+thr_y*6+thr_x],
					dx_dy_dz_DenAng[960+thr_y*6+thr_x],
					inv_Q_min,inv_Q_scaler,
					Den_min,Den_scaler,
					dx_dy_dz_DenAng[    thr_y*6+thr_x],
					dx_dy_dz_DenAng[192+thr_y*6+thr_x]
				    );
		}
		else if (threadIdx.y == 1)
		{
			// right halos
			int thr_x = threadIdx.x & 7;
			int thr_y = threadIdx.x & 24;
			Compute_VTI_DXED(
					((int*)dx_dy_dz_DenAng)[1152+thr_y*8+thr_x+48],
					dx_dy_dz_DenAng[    thr_y*6+thr_x+40],
					dx_dy_dz_DenAng[192+thr_y*6+thr_x+40],
					dx_dy_dz_DenAng[384+thr_y*6+thr_x+40],
					dx_dy_dz_DenAng[576+thr_y*6+thr_x+40],
					dx_dy_dz_DenAng[768+thr_y*6+thr_x+40],
					dx_dy_dz_DenAng[960+thr_y*6+thr_x+40],
					inv_Q_min,inv_Q_scaler,
					Den_min,Den_scaler,
					dx_dy_dz_DenAng[    thr_y*6+thr_x+40],
					dx_dy_dz_DenAng[192+thr_y*6+thr_x+40]
				    );
		}
		__syncthreads();

		// delay d(dxed1)dx and d(dxed2)dx by 4 Zzz's
		float d_dxed1_dx_p = next_d_dxed1_dx_p;
		float d_dxed1_dx_q = next_d_dxed1_dx_q; 

		// compute d(dxed1)dx and d(dxed2)dx
		next_d_dxed1_dx_p = Compute_DX_From_Shared(dx_dy_dz_DenAng +       threadIdx.y*48 + threadIdx.x + 7,a1h,a2h,a3h,a4h,a5h);
		next_d_dxed1_dx_q = Compute_DX_From_Shared(dx_dy_dz_DenAng + 192 + threadIdx.y*48 + threadIdx.x + 7,a1h,a2h,a3h,a4h,a5h);

		__syncthreads();

		// compute d(dzed1)dz and d(dzed2)dz
		dx_dy_dz_DenAng[256+threadIdx.y*32+threadIdx.x] = dzed2_p;
		dx_dy_dz_DenAng[384+threadIdx.y*32+threadIdx.x] = dzed2_q;

		__syncthreads();

		float* dzed1 = dx_dy_dz_DenAng + threadIdx.x;
		float* share_part_ddp = dx_dy_dz_DenAng + 1408 + threadIdx.x;

		float d_dzed2_dz_p, d_dzed2_dz_q;
		if (iZ > 4 && iZ < dimz-8)
		{
			// regular case
			if (threadIdx.y == 0)
			{
				d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed2_p    * a5;
				part_ddzed2dz_p[0] = part_ddzed2dz_p[1] + dzed2_p    * a1 + dzed1[256+32] * a2 + dzed1[256+64] * a3 + dzed1[256+96] * a4;
				part_ddzed2dz_p[1] = share_part_ddp[64] - dzed2_p    * a4 - dzed1[256+32] * a3 - dzed1[256+64] * a2 - dzed1[256+96] * a1;
				share_part_ddp[64] =                                                                                - dzed1[256+96] * a5;

				d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed2_q    * a5;
				part_ddzed2dz_q[0] = part_ddzed2dz_q[1] + dzed2_q    * a1 + dzed1[384+32] * a2 + dzed1[384+64] * a3 + dzed1[384+96] * a4;
				part_ddzed2dz_q[1] = share_part_ddp[96] - dzed2_q    * a4 - dzed1[384+32] * a3 - dzed1[384+64] * a2 - dzed1[384+96] * a1;
				share_part_ddp[96] =                                                                                - dzed1[384+96] * a5;
			}
			else if (threadIdx.y == 1)
			{
				d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a4 + dzed2_p       * a5;
				part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a1 + dzed2_p       * a1 + dzed1[256+64] * a2 + dzed1[256+96] * a3;
				part_ddzed2dz_p[1] =                    - dzed1[256] * a5 - dzed2_p       * a4 - dzed1[256+64] * a3 - dzed1[256+96] * a2;

				d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a4 + dzed2_q       * a5;
				part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a1 + dzed2_q       * a1 + dzed1[384+64] * a2 + dzed1[384+96] * a3;
				part_ddzed2dz_q[1] =                    - dzed1[384] * a5 - dzed2_q       * a4 - dzed1[384+64] * a3 - dzed1[384+96] * a2;
			}
			else if (threadIdx.y == 2)
			{
				d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a3 + dzed1[256+32] * a4 + dzed2_p       * a5;
				part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a2 - dzed1[256+32] * a1 + dzed2_p       * a1 + dzed1[256+96] * a2;
				part_ddzed2dz_p[1] =                                      - dzed1[256+32] * a5 - dzed2_p       * a4 - dzed1[256+96] * a3;

				d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a3 + dzed1[384+32] * a4 + dzed2_q       * a5;
				part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a2 - dzed1[384+32] * a1 + dzed2_q       * a1 + dzed1[384+96] * a2;
				part_ddzed2dz_q[1] =                                      - dzed1[384+32] * a5 - dzed2_q       * a4 - dzed1[384+96] * a3;
			}
			else if (threadIdx.y == 3)
			{
				d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a2 + dzed1[256+32] * a3 + dzed1[256+64] * a4 + dzed2_p       * a5;
				part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a3 - dzed1[256+32] * a2 - dzed1[256+64] * a1 + dzed2_p       * a1;
				part_ddzed2dz_p[1] =                                                           - dzed1[256+64] * a5 - dzed2_p       * a4;

				d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a2 + dzed1[384+32] * a3 + dzed1[384+64] * a4 + dzed2_q       * a5;
				part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a3 - dzed1[384+32] * a2 - dzed1[384+64] * a1 + dzed2_q       * a1;
				part_ddzed2dz_q[1] =                                                           - dzed1[384+64] * a5 - dzed2_q       * a4;
			}
		}
		else
		{
			if (iZ == 0)
			{
				if (threadIdx.y == 0)
				{
					part_ddzed2dz_p[1] = dzed2_p    * (-a4 - a5) - dzed1[256+32] *   a3       - dzed1[256+64] *  a2       - dzed1[256+96] *  a1;
					share_part_ddp[64] =                                                                                  - dzed1[256+96] *  a5;

					part_ddzed2dz_q[1] = dzed2_q    * (-a4 - a5) - dzed1[384+32] *   a3       - dzed1[384+64] *  a2       - dzed1[384+96] *  a1;
					share_part_ddp[96] =                                                                                  - dzed1[384+96] *  a5;
				}
				else if (threadIdx.y == 1)
				{
					part_ddzed2dz_p[0] = dzed1[256] * (-a1 - a2) + dzed2_p       * ( a1 - a3) + dzed1[256+64] * (a2 - a4) + dzed1[256+96] * (a3 - a5);
					part_ddzed2dz_p[1] = -dzed1[256] *  a5       - dzed2_p       *   a4       - dzed1[256+64] *  a3       - dzed1[256+96] *  a2;

					part_ddzed2dz_q[0] = dzed1[384] * (-a1 - a2) + dzed2_q       * ( a1 - a3) + dzed1[384+64] * (a2 - a4) + dzed1[384+96] * (a3 - a5);
					part_ddzed2dz_q[1] = -dzed1[384] *  a5       - dzed2_q       *   a4       - dzed1[384+64] *  a3       - dzed1[384+96] *  a2;
				}
				else if (threadIdx.y == 2)
				{
					part_ddzed2dz_p[0] = dzed1[256] * (-a2 - a3) + dzed1[256+32] * (-a1 - a4) + dzed2_p       * (a1 - a5) + dzed1[256+96] *  a2;
					part_ddzed2dz_p[1] =                           dzed1[256+32] *  -a5       - dzed2_p       *  a4       - dzed1[256+96] *  a3;

					part_ddzed2dz_q[0] = dzed1[384] * (-a2 - a3) + dzed1[384+32] * (-a1 - a4) + dzed2_q       * (a1 - a5) + dzed1[384+96] *  a2;
					part_ddzed2dz_q[1] =                           dzed1[384+32] *  -a5       - dzed2_q       *  a4       - dzed1[384+96] *  a3;

				}
				else if (threadIdx.y == 3)
				{
					part_ddzed2dz_p[0] = dzed1[256] * (-a3 - a4) + dzed1[256+32] * (-a2 - a5) - dzed1[256+64] *  a1       + dzed2_p       *  a1;
					part_ddzed2dz_p[1] =                                                        dzed1[256+64] * -a5       - dzed2_p       *  a4;

					part_ddzed2dz_q[0] = dzed1[384] * (-a3 - a4) + dzed1[384+32] * (-a2 - a5) - dzed1[384+64] *  a1       + dzed2_q       *  a1;
					part_ddzed2dz_q[1] =                                                        dzed1[384+64] * -a5       - dzed2_q       *  a4;
				}
			}
			else if (iZ == 4)
			{
				if (threadIdx.y == 0)
				{
					d_dzed2_dz_p = 0.0f;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] + dzed2_p    * a1 + dzed1[256+32] * a2 + dzed1[256+64] * a3 + dzed1[256+96] * a4;
					part_ddzed2dz_p[1] = share_part_ddp[64] - dzed2_p    * a4 - dzed1[256+32] * a3 - dzed1[256+64] * a2 - dzed1[256+96] * a1;
					share_part_ddp[64] =                                                                                - dzed1[256+96] * a5;

					d_dzed2_dz_q = 0.0f;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] + dzed2_q    * a1 + dzed1[384+32] * a2 + dzed1[384+64] * a3 + dzed1[384+96] * a4;
					part_ddzed2dz_q[1] = share_part_ddp[96] - dzed2_q    * a4 - dzed1[384+32] * a3 - dzed1[384+64] * a2 - dzed1[384+96] * a1;
					share_part_ddp[96] =                                                                                - dzed1[384+96] * a5;
				}
				else if (threadIdx.y == 1)
				{
					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a1 + dzed2_p       * a1 + dzed1[256+64] * a2 + dzed1[256+96] * a3;
					part_ddzed2dz_p[1] =                    - dzed1[256] * a5 - dzed2_p       * a4 - dzed1[256+64] * a3 - dzed1[256+96] * a2;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a1 + dzed2_q       * a1 + dzed1[384+64] * a2 + dzed1[384+96] * a3;
					part_ddzed2dz_q[1] =                    - dzed1[384] * a5 - dzed2_q       * a4 - dzed1[384+64] * a3 - dzed1[384+96] * a2;
				}
				else if (threadIdx.y == 2)
				{
					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a3 + dzed1[256+32] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a2 - dzed1[256+32] * a1 + dzed2_p       * a1 + dzed1[256+96] * a2;
					part_ddzed2dz_p[1] =                                      - dzed1[256+32] * a5 - dzed2_p       * a4 - dzed1[256+96] * a3;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a3 + dzed1[384+32] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a2 - dzed1[384+32] * a1 + dzed2_q       * a1 + dzed1[384+96] * a2;
					part_ddzed2dz_q[1] =                                      - dzed1[384+32] * a5 - dzed2_q       * a4 - dzed1[384+96] * a3;
				}
				else if (threadIdx.y == 3)
				{
					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a2 + dzed1[256+32] * a3 + dzed1[256+64] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a3 - dzed1[256+32] * a2 - dzed1[256+64] * a1 + dzed2_p       * a1;
					part_ddzed2dz_p[1] =                                                           - dzed1[256+64] * a5 - dzed2_p       * a4;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a2 + dzed1[384+32] * a3 + dzed1[384+64] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a3 - dzed1[384+32] * a2 - dzed1[384+64] * a1 + dzed2_q       * a1;
					part_ddzed2dz_q[1] =                                                           - dzed1[384+64] * a5 - dzed2_q       * a4;
				}
			}
			else if (iZ < dimz-4)
			{
				if (threadIdx.y == 0)
				{
					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed2_p    * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] + dzed2_p    * a1 + dzed1[256+32] * a2 + dzed1[256+64] * a3 + dzed1[256+96] * a4;
					part_ddzed2dz_p[1] =                                                                                - dzed1[256+96] * e1;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed2_q    * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] + dzed2_q    * a1 + dzed1[384+32] * a2 + dzed1[384+64] * a3 + dzed1[384+96] * a4;
					part_ddzed2dz_q[1] =                                                                                - dzed1[384+96] * e1;
				}
				else if (threadIdx.y == 1)
				{
					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a1 + dzed2_p       * a1 + dzed1[256+64] * a2 + dzed1[256+96] * a3;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a1 + dzed2_q       * a1 + dzed1[384+64] * a2 + dzed1[384+96] * a3;
				}
				else if (threadIdx.y == 2)
				{
					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a3 + dzed1[256+32] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a2 - dzed1[256+32] * a1 + dzed2_p       * a1 + dzed1[256+96] * a2;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a3 + dzed1[384+32] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a2 - dzed1[384+32] * a1 + dzed2_q       * a1 + dzed1[384+96] * a2;
				}
				else if (threadIdx.y == 3)
				{
					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a2 + dzed1[256+32] * a3 + dzed1[256+64] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a3 - dzed1[256+32] * a2 - dzed1[256+64] * a1 + dzed2_p       * a1;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a2 + dzed1[384+32] * a3 + dzed1[384+64] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a3 - dzed1[384+32] * a2 - dzed1[384+64] * a1 + dzed2_q       * a1;
				}
			}
			else if (iZ < dimz)
			{
				if (threadIdx.y == 0)
				{
					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed2_p    * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] + dzed2_p    * e1;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed2_q    * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] + dzed2_q    * e1;
				}
				else if (threadIdx.y == 1)
				{
					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] =                    - dzed1[256] * e1 + dzed2_p       * e1;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] =                    - dzed1[384] * e1 + dzed2_q       * e1;
				}
				else if (threadIdx.y == 2)
				{
					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a3 + dzed1[256+32] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] =                                      - dzed1[256+32] * e1 + dzed2_p       * e1;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a3 + dzed1[384+32] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] =                                      - dzed1[384+32] * e1 + dzed2_q       * e1;
				}
				else if (threadIdx.y == 3)
				{
					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a2 + dzed1[256+32] * a3 + dzed1[256+64] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] =                                                           - dzed1[256+64] * e1 + dzed2_p       * e1;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a2 + dzed1[384+32] * a3 + dzed1[384+64] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] =                                                           - dzed1[384+64] * e1 + dzed2_q       * e1;
				}
			}
			else
			{
				d_dzed2_dz_p = part_ddzed2dz_p[0];
				d_dzed2_dz_q = part_ddzed2dz_q[0];
			}
		}

		__syncthreads();

		if (iZ > 0)
		{
			// complete wave equation with X and Z derivatives.
			float V4_p = d_dzed2_dz_p;
			float V4_q = d_dzed2_dz_q;

			float V5_p = d_dxed1_dx_p;
			float V5_q = d_dxed1_dx_q;

			// write out V4 and V5 for dx and dz
			float2* dy = ((float2*)d_dy) + thr_y * stride_y_dy + (iZ-4) * stride_z_dy;

			float2 f2_V4;
			f2_V4.x = V4_p;
			f2_V4.y = V4_q;
			dy[threadIdx.y*64+threadIdx.x] = f2_V4;

			float2 f2_V5;
			f2_V5.x = V5_p;
			f2_V5.y = V5_q;
			dy[threadIdx.y*64+threadIdx.x+32] = f2_V5;
		}
			
		// this is necessary because next iteration will overwrite the part of shared memory that contains DZED1 and DZED2.
		__syncthreads();
	}
}

__global__
#if __CUDA_ARCH__ >= 300
__launch_bounds__(1024)
#elif __CUDA_ARCH__ >= 130
__launch_bounds__(640)
#endif
void Compute_VTI_DYED1and2_Part_V4_V5(
	float* d_dx,
	float* d_dy,
	float* d_dz,
	const int* d_em0,
	const int* d_em1,
	const int* d_em2,
	const int* d_em3,
	int dimy,
	int dimz,
	const float dt,
        const float a1,
        const float a2,
        const float a3,
        const float a4,
        const float a5,
        const float e1,
        const float a1h,
        const float a2h,
        const float a3h,
        const float a4h,
        const float a5h,
        float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
        float C44C33_min,
        float C44C33_scaler,
        float Vel_min,
        float Vel_scaler,
        float Del_min,
        float Del_scaler,
        float Eps_min,
        float Eps_scaler
	)
{
	__shared__ float dx_dy_dz_DenAng[1536];
	/*
	int nn = blockDim.x * blockDim.y * blockDim.z;
	for (int i = 0;  i < 1536;  i+=nn)
	{
		int idx = i+threadIdx.x+threadIdx.y*blockDim.x+threadIdx.z*blockDim.y*blockDim.x;
		if (idx < 1536) dx_dy_dz_DenAng[idx] = 0.0f;
	}
	__syncthreads();
	*/

	if (blockIdx.y < 5 || blockIdx.y >= dimy-13)
	{
		Compute_VTI_DYED1and2_Part_V4_V5_YHalo(
				dx_dy_dz_DenAng,
				d_dx,d_dy,d_dz,d_em0,d_em1,d_em2,d_em3,dimz,
				dt,a1,a2,a3,a4,a5,e1,a1h,a2h,a3h,a4h,a5h,
				inv_Q_min,inv_Q_scaler,Den_min,Den_scaler,C44C33_min,C44C33_scaler,Vel_min,Vel_scaler,Del_min,Del_scaler,Eps_min,Eps_scaler);
	}
	else
	{
		Compute_VTI_DYED1and2_Part_V4_V5_Main(
				dx_dy_dz_DenAng,
				d_dx,d_dy,d_dz,d_em0,d_em1,d_em2,d_em3,dimz,
				dt,a1,a2,a3,a4,a5,e1,a1h,a2h,a3h,a4h,a5h,
				inv_Q_min,inv_Q_scaler,Den_min,Den_scaler,C44C33_min,C44C33_scaler,Vel_min,Vel_scaler,Del_min,Del_scaler,Eps_min,Eps_scaler);
	}
}

__device__ 
void Compute_DYED1and2_Part_V4_V5_YHalo(
	float* dx_dy_dz_DenAng,
	float* d_dx,
        float* d_dy,
        float* d_dz,
        const int* d_em0,
        const int* d_em1,
        const int* d_em2,
        const int* d_em3,
        int dimz,
        const float dt,
        const float a1,
        const float a2,
        const float a3,
        const float a4,
        const float a5,
        const float e1,
        const float a1h,
        const float a2h,
        const float a3h,
        const float a4h,
        const float a5h,
        float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
        float Dip_min,
        float Dip_scaler,
        float Azm_min,
        float Azm_scaler,
        float C44C33_min,
        float C44C33_scaler,
        float Vel_min,
        float Vel_scaler,
        float Del_min,
        float Del_scaler,
        float Eps_min,
        float Eps_scaler
	)
{
	int stride_z_dy = 64;
	int stride_y_dy = stride_z_dy * dimz;

	int stride_z_dz = 48;
	int stride_y_dz = stride_z_dz * dimz;

	int stride_z_em = 16;
	int stride_y_em = stride_z_em * dimz;

	int thr_y = blockIdx.y + 4;

	//
	// block shape is 32 by 4 by 1 (x by y by z)
	// algorithm loops over z.
	//

#ifdef NAN_DESU_KA
	int fu2 = 0;
#endif

	for (int iZ = 0;  iZ < dimz;  iZ+=4)
	{
		if (threadIdx.y == 0)
		{
			float2* dx = (float2*)d_dx + thr_y * stride_y_dy + iZ * stride_z_dy;
			const int* em = d_em0 + thr_y * 2 * stride_y_em + iZ * stride_z_em;

			float2 v = dx[threadIdx.x];
			dx_dy_dz_DenAng[    threadIdx.x] = v.x;
			dx_dy_dz_DenAng[192+threadIdx.x] = v.y;

			v = dx[32+threadIdx.x];
			dx_dy_dz_DenAng[ 32+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[224+threadIdx.x] = v.y;

			v = dx[64+threadIdx.x];
			dx_dy_dz_DenAng[ 64+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[256+threadIdx.x] = v.y;

			v = dx[96+threadIdx.x];
			dx_dy_dz_DenAng[ 96+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[288+threadIdx.x] = v.y;

			v = dx[128+threadIdx.x];
			dx_dy_dz_DenAng[128+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[320+threadIdx.x] = v.y;

			v = dx[160+threadIdx.x];
			dx_dy_dz_DenAng[160+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[352+threadIdx.x] = v.y;

			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		else if (threadIdx.y == 1)
		{
			float2* dx = (float2*)d_dy + thr_y * stride_y_dy + iZ * stride_z_dy;
			const int* em = d_em1 + thr_y * 2 * stride_y_em + iZ * stride_z_em;

			float2 v = dx[threadIdx.x];
			dx_dy_dz_DenAng[384+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[576+threadIdx.x] = v.y;

			v = dx[32+threadIdx.x];
			dx_dy_dz_DenAng[416+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[608+threadIdx.x] = v.y;

			v = dx[64+threadIdx.x];
			dx_dy_dz_DenAng[448+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[640+threadIdx.x] = v.y;

			v = dx[96+threadIdx.x];
			dx_dy_dz_DenAng[480+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[672+threadIdx.x] = v.y;

			v = dx[128+threadIdx.x];
			dx_dy_dz_DenAng[512+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[704+threadIdx.x] = v.y;

			v = dx[160+threadIdx.x];
			dx_dy_dz_DenAng[544+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[736+threadIdx.x] = v.y;

			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		else if (threadIdx.y == 2)
		{
			float2* dx = (float2*)d_dz + thr_y * stride_y_dz + iZ * stride_z_dz;
			const int* em = d_em2 + thr_y * 2 * stride_y_em + iZ * stride_z_em;

			float2 v = dx[threadIdx.x];
			dx_dy_dz_DenAng[ 768+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[ 960+threadIdx.x] = v.y;

			v = dx[32+threadIdx.x];
			dx_dy_dz_DenAng[ 800+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[ 992+threadIdx.x] = v.y;

			v = dx[64+threadIdx.x];
			dx_dy_dz_DenAng[ 832+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1024+threadIdx.x] = v.y;

			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		else if (threadIdx.y == 3)
		{
			float2* dx = (float2*)d_dz + thr_y * stride_y_dz + iZ * stride_z_dz;
			const int* em = d_em3 + thr_y * 2 * stride_y_em + iZ * stride_z_em;

			float2 v = dx[96+threadIdx.x];
			dx_dy_dz_DenAng[ 864+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1056+threadIdx.x] = v.y;

			v = dx[128+threadIdx.x];
			dx_dy_dz_DenAng[ 896+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1088+threadIdx.x] = v.y;

			v = dx[160+threadIdx.x];
			dx_dy_dz_DenAng[ 928+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1120+threadIdx.x] = v.y;

			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		__syncthreads();

		// dx is in dx_dy_dz_DenAng[  0, 384>, [0,96> contains 1st line, [96,192> contains 2nd line etc.
		// dy is in dx_dy_dz_DenAng[384, 768>
		// dz is in dx_dy_dz_DenAng[768,1152>

		// compute dxed1, dxed2, dyed1, dyed2, dzed1 and dzed2 for center points first.

		float dyed1_p, dyed1_q, dyed2_p, dyed2_q;
		Compute_DYED(
				((int*)dx_dy_dz_DenAng)[1152+threadIdx.y*64+threadIdx.x+16],
				dx_dy_dz_DenAng[    threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[192+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[384+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[576+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[768+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[960+threadIdx.y*48+threadIdx.x+8],
				inv_Q_min,inv_Q_scaler,
				Den_min,Den_scaler,
				Dip_min,Dip_scaler,
				Azm_min,Azm_scaler,
				dyed1_p, dyed1_q, dyed2_p, dyed2_q
			    );

		// write out dyed1 and dyed2
		float2 f2_dyed1;
		f2_dyed1.x = dyed1_p;
		f2_dyed1.y = dyed1_q;
		float2* dx = ((float2*)d_dx) + thr_y * stride_y_dy + iZ * stride_z_dy;
		dx[threadIdx.y*64+threadIdx.x] = f2_dyed1;

		float2 f2_dyed2;
		f2_dyed2.x = dyed2_p;
		f2_dyed2.y = dyed2_q;
		dx[threadIdx.y*64+threadIdx.x+32] = f2_dyed2;

#ifdef NAN_DESU_KA
		if (!fu2 && (dyed1_p != 0.0f || dyed2_p != 0.0f || dyed1_q != 0.0f || dyed2_q != 0.0f))
		{
			fu2 = 1;
			printf("Compute_DYED1andDYED2...halo ; threadIdx=[%d,%d,%d] blockIdx=[%d,%d] :: dyed1_p=%e, dyed2_p=%e, dyed1_q=%e, dyed2_q=%e\n",threadIdx.x,threadIdx.y,threadIdx.z,blockIdx.x,blockIdx.y,dyed1_p,dyed2_p,dyed1_q,dyed2_q);
		}
#endif
	
		__syncthreads();
	}
}

__device__
void Compute_DYED1and2_Part_V4_V5_Main(
	float* dx_dy_dz_DenAng,
	float* d_dx,
	float* d_dy,
	float* d_dz,
	const int* d_em0,
	const int* d_em1,
	const int* d_em2,
	const int* d_em3,
	int dimz,
	const float dt,
        const float a1,
        const float a2,
        const float a3,
        const float a4,
        const float a5,
        const float e1,
        const float a1h,
        const float a2h,
        const float a3h,
        const float a4h,
        const float a5h,
        float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
        float Dip_min,
        float Dip_scaler,
        float Azm_min,
        float Azm_scaler,
        float C44C33_min,
        float C44C33_scaler,
        float Vel_min,
        float Vel_scaler,
        float Del_min,
        float Del_scaler,
        float Eps_min,
        float Eps_scaler
	)
{
	const int stride_z_dy = 64;
	const int stride_y_dy = stride_z_dy * dimz;

	const int stride_z_dz = 48;
	const int stride_y_dz = stride_z_dz * dimz;

	const int stride_z_em = 16;
	const int stride_y_em = stride_z_em * dimz;

	const int thr_y = blockIdx.y + 4;

	//
	// block shape is 32 by 4 by 1 (x by y by z)
	// algorithm loops over z.
	//

	float part_ddzed1dz_p[2] = {0.0f, 0.0f};
	float part_ddzed1dz_q[2] = {0.0f, 0.0f};
	float part_ddzed2dz_p[2] = {0.0f, 0.0f};
	float part_ddzed2dz_q[2] = {0.0f, 0.0f};

#ifdef NAN_DESU_KA
	int fu1 = 0, fu2 = 0, fu3 = 0;
#endif

	float next_d_dxed1_dx_p=0.0f, next_d_dxed1_dx_q=0.0f, next_d_dxed2_dx_p=0.0f, next_d_dxed2_dx_q=0.0f;
	for (int iZ = 0;  iZ < dimz;  iZ+=4)
	{
		if (threadIdx.y == 0)
		{
			float2* dx = (float2*)d_dx + thr_y * stride_y_dy + iZ * stride_z_dy;
			const int* em = d_em0 + thr_y * 2 * stride_y_em + iZ * stride_z_em;
			
			float2 v = dx[threadIdx.x];
			dx_dy_dz_DenAng[    threadIdx.x] = v.x;
			dx_dy_dz_DenAng[192+threadIdx.x] = v.y;
			
			v = dx[32+threadIdx.x];
			dx_dy_dz_DenAng[ 32+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[224+threadIdx.x] = v.y;

			v = dx[64+threadIdx.x];
			dx_dy_dz_DenAng[ 64+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[256+threadIdx.x] = v.y;

			v = dx[96+threadIdx.x];
			dx_dy_dz_DenAng[ 96+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[288+threadIdx.x] = v.y;

			v = dx[128+threadIdx.x];
			dx_dy_dz_DenAng[128+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[320+threadIdx.x] = v.y;
			
			v = dx[160+threadIdx.x];
			dx_dy_dz_DenAng[160+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[352+threadIdx.x] = v.y;
		
			// NB!!!!
			// Without typecast, input value from em ptr will be converted to float before storage.
			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		else if (threadIdx.y == 1)
		{
			float2* dx = (float2*)d_dy + thr_y * stride_y_dy + iZ * stride_z_dy;
			const int* em = d_em1 + thr_y * 2 * stride_y_em + iZ * stride_z_em;

			float2 v = dx[threadIdx.x];
			dx_dy_dz_DenAng[384+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[576+threadIdx.x] = v.y;
			
			v = dx[32+threadIdx.x];
			dx_dy_dz_DenAng[416+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[608+threadIdx.x] = v.y;

			v = dx[64+threadIdx.x];
			dx_dy_dz_DenAng[448+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[640+threadIdx.x] = v.y;

			v = dx[96+threadIdx.x];
			dx_dy_dz_DenAng[480+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[672+threadIdx.x] = v.y;

			v = dx[128+threadIdx.x];
			dx_dy_dz_DenAng[512+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[704+threadIdx.x] = v.y;
			
			v = dx[160+threadIdx.x];
			dx_dy_dz_DenAng[544+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[736+threadIdx.x] = v.y;
		
			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		else if (threadIdx.y == 2)
		{
			float2* dx = (float2*)d_dz + thr_y * stride_y_dz + iZ * stride_z_dz;
			const int* em = d_em2 + thr_y * 2 * stride_y_em + iZ * stride_z_em;

			float2 v = dx[threadIdx.x];
			dx_dy_dz_DenAng[ 768+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[ 960+threadIdx.x] = v.y;
			
			v = dx[32+threadIdx.x];
			dx_dy_dz_DenAng[ 800+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[ 992+threadIdx.x] = v.y;

			v = dx[64+threadIdx.x];
			dx_dy_dz_DenAng[ 832+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1024+threadIdx.x] = v.y;

			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		else if (threadIdx.y == 3)
		{
			float2* dx = (float2*)d_dz + thr_y * stride_y_dz + iZ * stride_z_dz;
			const int* em = d_em3 + thr_y * 2 * stride_y_em + iZ * stride_z_em;

			float2 v = dx[96+threadIdx.x];
			dx_dy_dz_DenAng[ 864+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1056+threadIdx.x] = v.y;

			v = dx[128+threadIdx.x];
			dx_dy_dz_DenAng[ 896+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1088+threadIdx.x] = v.y;
			
			v = dx[160+threadIdx.x];
			dx_dy_dz_DenAng[ 928+threadIdx.x] = v.x;
			dx_dy_dz_DenAng[1120+threadIdx.x] = v.y;
			
			((int*)dx_dy_dz_DenAng)[1152+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[   threadIdx.x];
			((int*)dx_dy_dz_DenAng)[1280+(threadIdx.x&15)+((threadIdx.x&16)*4)+threadIdx.y*16] = em[32+threadIdx.x];
		}
		__syncthreads();

		// dx is in dx_dy_dz_DenAng[  0, 384>, [0,96> contains 1st line, [96,192> contains 2nd line etc.
		// dy is in dx_dy_dz_DenAng[384, 768>
		// dz is in dx_dy_dz_DenAng[768,1152>

		// compute dxed1, dxed2, dyed1, dyed2, dzed1 and dzed2 for center points first.
		float dyed1_p, dyed1_q, dyed2_p, dyed2_q;
		float dzed1_p, dzed1_q, dzed2_p, dzed2_q;
		Compute_DXED_DYED_DZED(
				((int*)dx_dy_dz_DenAng)[1152+threadIdx.y*64+threadIdx.x+16],
				dx_dy_dz_DenAng[    threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[192+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[384+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[576+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[768+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[960+threadIdx.y*48+threadIdx.x+8],
				inv_Q_min,inv_Q_scaler,
				Den_min,Den_scaler,
				Dip_min,Dip_scaler,
				Azm_min,Azm_scaler,
				dx_dy_dz_DenAng[    threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[192+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[384+threadIdx.y*48+threadIdx.x+8],
				dx_dy_dz_DenAng[576+threadIdx.y*48+threadIdx.x+8],
				dyed1_p, dyed1_q, dyed2_p, dyed2_q,
				dzed1_p, dzed1_q, dzed2_p, dzed2_q
				);

		// write out dyed1 and dyed2
		float2 f2_dyed1;
		f2_dyed1.x = dyed1_p;
		f2_dyed1.y = dyed1_q;
		float2* dx = ((float2*)d_dx) + thr_y * stride_y_dy + iZ * stride_z_dy;
		dx[threadIdx.y*64+threadIdx.x] = f2_dyed1;

		float2 f2_dyed2;
		f2_dyed2.x = dyed2_p;
		f2_dyed2.y = dyed2_q;
		dx[threadIdx.y*64+threadIdx.x+32] = f2_dyed2;

#ifdef NAN_DESU_KA
		if (!fu2 && (dyed1_p != 0.0f || dyed2_p != 0.0f || dyed1_q != 0.0f || dyed2_q != 0.0f))
		{
			fu2 = 1;
			printf("Compute_DYED1andDYED2... ; threadIdx=[%d,%d,%d] blockIdx=[%d,%d] :: dyed1_p=%e, dyed2_p=%e, dyed1_q=%e, dyed2_q=%e\n",threadIdx.x,threadIdx.y,threadIdx.z,blockIdx.x,blockIdx.y,dyed1_p,dyed2_p,dyed1_q,dyed2_q);
		}
#endif

		// compute dxed1 and dxed2 for halo cells
		if (threadIdx.y == 0)
		{
			// left halos
			int thr_x = threadIdx.x & 7;
			int thr_y = threadIdx.x & 24;
			Compute_DXED(
					((int*)dx_dy_dz_DenAng)[1152+thr_y*8+thr_x+8],
					dx_dy_dz_DenAng[    thr_y*6+thr_x],
					dx_dy_dz_DenAng[192+thr_y*6+thr_x],
					dx_dy_dz_DenAng[384+thr_y*6+thr_x],
					dx_dy_dz_DenAng[576+thr_y*6+thr_x],
					dx_dy_dz_DenAng[768+thr_y*6+thr_x],
					dx_dy_dz_DenAng[960+thr_y*6+thr_x],
					inv_Q_min,inv_Q_scaler,
					Den_min,Den_scaler,
					Dip_min,Dip_scaler,
					Azm_min,Azm_scaler,
					dx_dy_dz_DenAng[    thr_y*6+thr_x],
					dx_dy_dz_DenAng[192+thr_y*6+thr_x],
					dx_dy_dz_DenAng[384+thr_y*6+thr_x],
					dx_dy_dz_DenAng[576+thr_y*6+thr_x]
				    );
		}
		else if (threadIdx.y == 1)
		{
			// right halos
			int thr_x = threadIdx.x & 7;
			int thr_y = threadIdx.x & 24;
			Compute_DXED(
					((int*)dx_dy_dz_DenAng)[1152+thr_y*8+thr_x+48],
					dx_dy_dz_DenAng[    thr_y*6+thr_x+40],
					dx_dy_dz_DenAng[192+thr_y*6+thr_x+40],
					dx_dy_dz_DenAng[384+thr_y*6+thr_x+40],
					dx_dy_dz_DenAng[576+thr_y*6+thr_x+40],
					dx_dy_dz_DenAng[768+thr_y*6+thr_x+40],
					dx_dy_dz_DenAng[960+thr_y*6+thr_x+40],
					inv_Q_min,inv_Q_scaler,
					Den_min,Den_scaler,
					Dip_min,Dip_scaler,
					Azm_min,Azm_scaler,
					dx_dy_dz_DenAng[    thr_y*6+thr_x+40],
					dx_dy_dz_DenAng[192+thr_y*6+thr_x+40],
					dx_dy_dz_DenAng[384+thr_y*6+thr_x+40],
					dx_dy_dz_DenAng[576+thr_y*6+thr_x+40]
				    );
		}
		__syncthreads();

		// delay d(dxed1)dx and d(dxed2)dx by 4 Zzz's
		float d_dxed1_dx_p = next_d_dxed1_dx_p;
		float d_dxed1_dx_q = next_d_dxed1_dx_q; 
		float d_dxed2_dx_p = next_d_dxed2_dx_p;
		float d_dxed2_dx_q = next_d_dxed2_dx_q;

		// compute d(dxed1)dx and d(dxed2)dx
		next_d_dxed1_dx_p = Compute_DX_From_Shared(dx_dy_dz_DenAng +       threadIdx.y*48 + threadIdx.x + 7,a1h,a2h,a3h,a4h,a5h);
		next_d_dxed1_dx_q = Compute_DX_From_Shared(dx_dy_dz_DenAng + 192 + threadIdx.y*48 + threadIdx.x + 7,a1h,a2h,a3h,a4h,a5h);
		next_d_dxed2_dx_p = Compute_DX_From_Shared(dx_dy_dz_DenAng + 384 + threadIdx.y*48 + threadIdx.x + 7,a1h,a2h,a3h,a4h,a5h);
		next_d_dxed2_dx_q = Compute_DX_From_Shared(dx_dy_dz_DenAng + 576 + threadIdx.y*48 + threadIdx.x + 7,a1h,a2h,a3h,a4h,a5h);

		__syncthreads();

		// compute d(dzed1)dz and d(dzed2)dz
		dx_dy_dz_DenAng[    threadIdx.y*32+threadIdx.x] = dzed1_p;
		dx_dy_dz_DenAng[128+threadIdx.y*32+threadIdx.x] = dzed1_q;
		dx_dy_dz_DenAng[256+threadIdx.y*32+threadIdx.x] = dzed2_p;
		dx_dy_dz_DenAng[384+threadIdx.y*32+threadIdx.x] = dzed2_q;

		__syncthreads();

		float* dzed1 = dx_dy_dz_DenAng + threadIdx.x;
		float* share_part_ddp = dx_dy_dz_DenAng + 1408 + threadIdx.x;

		float d_dzed1_dz_p, d_dzed1_dz_q, d_dzed2_dz_p, d_dzed2_dz_q;
		if (iZ > 4 && iZ < dimz-8)
		{
			// regular case
			if (threadIdx.y == 0)
			{
				d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1_p    * a5;
				part_ddzed1dz_p[0] = part_ddzed1dz_p[1] + dzed1_p    * a1 + dzed1[    32] * a2 + dzed1[    64] * a3 + dzed1[    96] * a4;
				part_ddzed1dz_p[1] = share_part_ddp[ 0] - dzed1_p    * a4 - dzed1[    32] * a3 - dzed1[    64] * a2 - dzed1[    96] * a1;
				share_part_ddp[ 0] =                                                                                - dzed1[    96] * a5;

				d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1_q    * a5;
				part_ddzed1dz_q[0] = part_ddzed1dz_q[1] + dzed1_q    * a1 + dzed1[128+32] * a2 + dzed1[128+64] * a3 + dzed1[128+96] * a4;
				part_ddzed1dz_q[1] = share_part_ddp[32] - dzed1_q    * a4 - dzed1[128+32] * a3 - dzed1[128+64] * a2 - dzed1[128+96] * a1;
				share_part_ddp[32] =                                                                                - dzed1[128+96] * a5;

				d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed2_p    * a5;
				part_ddzed2dz_p[0] = part_ddzed2dz_p[1] + dzed2_p    * a1 + dzed1[256+32] * a2 + dzed1[256+64] * a3 + dzed1[256+96] * a4;
				part_ddzed2dz_p[1] = share_part_ddp[64] - dzed2_p    * a4 - dzed1[256+32] * a3 - dzed1[256+64] * a2 - dzed1[256+96] * a1;
				share_part_ddp[64] =                                                                                - dzed1[256+96] * a5;

				d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed2_q    * a5;
				part_ddzed2dz_q[0] = part_ddzed2dz_q[1] + dzed2_q    * a1 + dzed1[384+32] * a2 + dzed1[384+64] * a3 + dzed1[384+96] * a4;
				part_ddzed2dz_q[1] = share_part_ddp[96] - dzed2_q    * a4 - dzed1[384+32] * a3 - dzed1[384+64] * a2 - dzed1[384+96] * a1;
				share_part_ddp[96] =                                                                                - dzed1[384+96] * a5;
			}
			else if (threadIdx.y == 1)
			{
				d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1[  0] * a4 + dzed1_p       * a5;
				part_ddzed1dz_p[0] = part_ddzed1dz_p[1] - dzed1[  0] * a1 + dzed1_p       * a1 + dzed1[    64] * a2 + dzed1[    96] * a3;
				part_ddzed1dz_p[1] =                    - dzed1[  0] * a5 - dzed1_p       * a4 - dzed1[    64] * a3 - dzed1[    96] * a2;

				d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1[128] * a4 + dzed1_q       * a5;
				part_ddzed1dz_q[0] = part_ddzed1dz_q[1] - dzed1[128] * a1 + dzed1_q       * a1 + dzed1[128+64] * a2 + dzed1[128+96] * a3;
				part_ddzed1dz_q[1] =                    - dzed1[128] * a5 - dzed1_q       * a4 - dzed1[128+64] * a3 - dzed1[128+96] * a2;

				d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a4 + dzed2_p       * a5;
				part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a1 + dzed2_p       * a1 + dzed1[256+64] * a2 + dzed1[256+96] * a3;
				part_ddzed2dz_p[1] =                    - dzed1[256] * a5 - dzed2_p       * a4 - dzed1[256+64] * a3 - dzed1[256+96] * a2;

				d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a4 + dzed2_q       * a5;
				part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a1 + dzed2_q       * a1 + dzed1[384+64] * a2 + dzed1[384+96] * a3;
				part_ddzed2dz_q[1] =                    - dzed1[384] * a5 - dzed2_q       * a4 - dzed1[384+64] * a3 - dzed1[384+96] * a2;
			}
			else if (threadIdx.y == 2)
			{
				d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1[  0] * a3 + dzed1[    32] * a4 + dzed1_p       * a5;
				part_ddzed1dz_p[0] = part_ddzed1dz_p[1] - dzed1[  0] * a2 - dzed1[    32] * a1 + dzed1_p       * a1 + dzed1[    96] * a2;
				part_ddzed1dz_p[1] =                                      - dzed1[    32] * a5 - dzed1_p       * a4 - dzed1[    96] * a3;

				d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1[128] * a3 + dzed1[128+32] * a4 + dzed1_q       * a5;
				part_ddzed1dz_q[0] = part_ddzed1dz_q[1] - dzed1[128] * a2 - dzed1[128+32] * a1 + dzed1_q       * a1 + dzed1[128+96] * a2;
				part_ddzed1dz_q[1] =                                      - dzed1[128+32] * a5 - dzed1_q       * a4 - dzed1[128+96] * a3;

				d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a3 + dzed1[256+32] * a4 + dzed2_p       * a5;
				part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a2 - dzed1[256+32] * a1 + dzed2_p       * a1 + dzed1[256+96] * a2;
				part_ddzed2dz_p[1] =                                      - dzed1[256+32] * a5 - dzed2_p       * a4 - dzed1[256+96] * a3;

				d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a3 + dzed1[384+32] * a4 + dzed2_q       * a5;
				part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a2 - dzed1[384+32] * a1 + dzed2_q       * a1 + dzed1[384+96] * a2;
				part_ddzed2dz_q[1] =                                      - dzed1[384+32] * a5 - dzed2_q       * a4 - dzed1[384+96] * a3;
			}
			else if (threadIdx.y == 3)
			{
				d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1[  0] * a2 + dzed1[    32] * a3 + dzed1[    64] * a4 + dzed1_p       * a5;
				part_ddzed1dz_p[0] = part_ddzed1dz_p[1] - dzed1[  0] * a3 - dzed1[    32] * a2 - dzed1[    64] * a1 + dzed1_p       * a1;
				part_ddzed1dz_p[1] =                                                           - dzed1[    64] * a5 - dzed1_p       * a4;

				d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1[128] * a2 + dzed1[128+32] * a3 + dzed1[128+64] * a4 + dzed1_q       * a5;
				part_ddzed1dz_q[0] = part_ddzed1dz_q[1] - dzed1[128] * a3 - dzed1[128+32] * a2 - dzed1[128+64] * a1 + dzed1_q       * a1;
				part_ddzed1dz_q[1] =                                                           - dzed1[128+64] * a5 - dzed1_q       * a4;

				d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a2 + dzed1[256+32] * a3 + dzed1[256+64] * a4 + dzed2_p       * a5;
				part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a3 - dzed1[256+32] * a2 - dzed1[256+64] * a1 + dzed2_p       * a1;
				part_ddzed2dz_p[1] =                                                           - dzed1[256+64] * a5 - dzed2_p       * a4;

				d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a2 + dzed1[384+32] * a3 + dzed1[384+64] * a4 + dzed2_q       * a5;
				part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a3 - dzed1[384+32] * a2 - dzed1[384+64] * a1 + dzed2_q       * a1;
				part_ddzed2dz_q[1] =                                                           - dzed1[384+64] * a5 - dzed2_q       * a4;
			}
		}
		else
		{
			if (iZ == 0)
			{
				if (threadIdx.y == 0)
				{
					part_ddzed1dz_p[1] = dzed1_p    * (-a4 - a5) - dzed1[    32] *   a3       - dzed1[    64] *  a2       - dzed1[    96] *  a1;
					share_part_ddp[ 0] =                                                                                  - dzed1[    96] *  a5;

					part_ddzed1dz_q[1] = dzed1_q    * (-a4 - a5) - dzed1[128+32] *   a3       - dzed1[128+64] *  a2       - dzed1[128+96] *  a1;
					share_part_ddp[32] =                                                                                  - dzed1[128+96] *  a5;

					part_ddzed2dz_p[1] = dzed2_p    * (-a4 - a5) - dzed1[256+32] *   a3       - dzed1[256+64] *  a2       - dzed1[256+96] *  a1;
					share_part_ddp[64] =                                                                                  - dzed1[256+96] *  a5;

					part_ddzed2dz_q[1] = dzed2_q    * (-a4 - a5) - dzed1[384+32] *   a3       - dzed1[384+64] *  a2       - dzed1[384+96] *  a1;
					share_part_ddp[96] =                                                                                  - dzed1[384+96] *  a5;
				}
				else if (threadIdx.y == 1)
				{
					part_ddzed1dz_p[0] = dzed1[  0] * (-a1 - a2) + dzed1_p       * ( a1 - a3) + dzed1[    64] * (a2 - a4) + dzed1[    96] * (a3 - a5);
					part_ddzed1dz_p[1] = -dzed1[  0] *  a5       - dzed1_p       *   a4       - dzed1[    64] *  a3       - dzed1[    96] *  a2;

					part_ddzed1dz_q[0] = dzed1[128] * (-a1 - a2) + dzed1_q       * ( a1 - a3) + dzed1[128+64] * (a2 - a4) + dzed1[128+96] * (a3 - a5);
					part_ddzed1dz_q[1] = -dzed1[128] *  a5       - dzed1_q       *   a4       - dzed1[128+64] *  a3       - dzed1[128+96] *  a2;

					part_ddzed2dz_p[0] = dzed1[256] * (-a1 - a2) + dzed2_p       * ( a1 - a3) + dzed1[256+64] * (a2 - a4) + dzed1[256+96] * (a3 - a5);
					part_ddzed2dz_p[1] = -dzed1[256] *  a5       - dzed2_p       *   a4       - dzed1[256+64] *  a3       - dzed1[256+96] *  a2;

					part_ddzed2dz_q[0] = dzed1[384] * (-a1 - a2) + dzed2_q       * ( a1 - a3) + dzed1[384+64] * (a2 - a4) + dzed1[384+96] * (a3 - a5);
					part_ddzed2dz_q[1] = -dzed1[384] *  a5       - dzed2_q       *   a4       - dzed1[384+64] *  a3       - dzed1[384+96] *  a2;
				}
				else if (threadIdx.y == 2)
				{
					part_ddzed1dz_p[0] = dzed1[  0] * (-a2 - a3) + dzed1[    32] * (-a1 - a4) + dzed1_p       * (a1 - a5) + dzed1[    96] *  a2;
					part_ddzed1dz_p[1] =                           dzed1[    32] *  -a5       - dzed1_p       *  a4       - dzed1[    96] *  a3;

					part_ddzed1dz_q[0] = dzed1[128] * (-a2 - a3) + dzed1[128+32] * (-a1 - a4) + dzed1_q       * (a1 - a5) + dzed1[128+96] *  a2;
					part_ddzed1dz_q[1] =                           dzed1[128+32] *  -a5       - dzed1_q       *  a4       - dzed1[128+96] *  a3;

					part_ddzed2dz_p[0] = dzed1[256] * (-a2 - a3) + dzed1[256+32] * (-a1 - a4) + dzed2_p       * (a1 - a5) + dzed1[256+96] *  a2;
					part_ddzed2dz_p[1] =                           dzed1[256+32] *  -a5       - dzed2_p       *  a4       - dzed1[256+96] *  a3;

					part_ddzed2dz_q[0] = dzed1[384] * (-a2 - a3) + dzed1[384+32] * (-a1 - a4) + dzed2_q       * (a1 - a5) + dzed1[384+96] *  a2;
					part_ddzed2dz_q[1] =                           dzed1[384+32] *  -a5       - dzed2_q       *  a4       - dzed1[384+96] *  a3;

				}
				else if (threadIdx.y == 3)
				{
					part_ddzed1dz_p[0] = dzed1[  0] * (-a3 - a4) + dzed1[    32] * (-a2 - a5) - dzed1[    64] *  a1       + dzed1_p       *  a1;
					part_ddzed1dz_p[1] =                                                        dzed1[    64] * -a5       - dzed1_p       *  a4;

					part_ddzed1dz_q[0] = dzed1[128] * (-a3 - a4) + dzed1[128+32] * (-a2 - a5) - dzed1[128+64] *  a1       + dzed1_q       *  a1;
					part_ddzed1dz_q[1] =                                                        dzed1[128+64] * -a5       - dzed1_q       *  a4;

					part_ddzed2dz_p[0] = dzed1[256] * (-a3 - a4) + dzed1[256+32] * (-a2 - a5) - dzed1[256+64] *  a1       + dzed2_p       *  a1;
					part_ddzed2dz_p[1] =                                                        dzed1[256+64] * -a5       - dzed2_p       *  a4;

					part_ddzed2dz_q[0] = dzed1[384] * (-a3 - a4) + dzed1[384+32] * (-a2 - a5) - dzed1[384+64] *  a1       + dzed2_q       *  a1;
					part_ddzed2dz_q[1] =                                                        dzed1[384+64] * -a5       - dzed2_q       *  a4;
				}
			}
			else if (iZ == 4)
			{
				if (threadIdx.y == 0)
				{
					d_dzed1_dz_p = 0.0f;
					part_ddzed1dz_p[0] = part_ddzed1dz_p[1] + dzed1_p    * a1 + dzed1[    32] * a2 + dzed1[    64] * a3 + dzed1[    96] * a4;
					part_ddzed1dz_p[1] = share_part_ddp[ 0] - dzed1_p    * a4 - dzed1[    32] * a3 - dzed1[    64] * a2 - dzed1[    96] * a1;
					share_part_ddp[ 0] =                                                                                - dzed1[    96] * a5;

					d_dzed1_dz_q = 0.0f;
					part_ddzed1dz_q[0] = part_ddzed1dz_q[1] + dzed1_q    * a1 + dzed1[128+32] * a2 + dzed1[128+64] * a3 + dzed1[128+96] * a4;
					part_ddzed1dz_q[1] = share_part_ddp[32] - dzed1_q    * a4 - dzed1[128+32] * a3 - dzed1[128+64] * a2 - dzed1[128+96] * a1;
					share_part_ddp[32] =                                                                                - dzed1[128+96] * a5;

					d_dzed2_dz_p = 0.0f;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] + dzed2_p    * a1 + dzed1[256+32] * a2 + dzed1[256+64] * a3 + dzed1[256+96] * a4;
					part_ddzed2dz_p[1] = share_part_ddp[64] - dzed2_p    * a4 - dzed1[256+32] * a3 - dzed1[256+64] * a2 - dzed1[256+96] * a1;
					share_part_ddp[64] =                                                                                - dzed1[256+96] * a5;

					d_dzed2_dz_q = 0.0f;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] + dzed2_q    * a1 + dzed1[384+32] * a2 + dzed1[384+64] * a3 + dzed1[384+96] * a4;
					part_ddzed2dz_q[1] = share_part_ddp[96] - dzed2_q    * a4 - dzed1[384+32] * a3 - dzed1[384+64] * a2 - dzed1[384+96] * a1;
					share_part_ddp[96] =                                                                                - dzed1[384+96] * a5;
				}
				else if (threadIdx.y == 1)
				{
					d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1[  0] * a4 + dzed1_p       * a5;
					part_ddzed1dz_p[0] = part_ddzed1dz_p[1] - dzed1[  0] * a1 + dzed1_p       * a1 + dzed1[    64] * a2 + dzed1[    96] * a3;
					part_ddzed1dz_p[1] =                    - dzed1[  0] * a5 - dzed1_p       * a4 - dzed1[    64] * a3 - dzed1[    96] * a2;

					d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1[128] * a4 + dzed1_q       * a5;
					part_ddzed1dz_q[0] = part_ddzed1dz_q[1] - dzed1[128] * a1 + dzed1_q       * a1 + dzed1[128+64] * a2 + dzed1[128+96] * a3;
					part_ddzed1dz_q[1] =                    - dzed1[128] * a5 - dzed1_q       * a4 - dzed1[128+64] * a3 - dzed1[128+96] * a2;

					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a1 + dzed2_p       * a1 + dzed1[256+64] * a2 + dzed1[256+96] * a3;
					part_ddzed2dz_p[1] =                    - dzed1[256] * a5 - dzed2_p       * a4 - dzed1[256+64] * a3 - dzed1[256+96] * a2;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a1 + dzed2_q       * a1 + dzed1[384+64] * a2 + dzed1[384+96] * a3;
					part_ddzed2dz_q[1] =                    - dzed1[384] * a5 - dzed2_q       * a4 - dzed1[384+64] * a3 - dzed1[384+96] * a2;
				}
				else if (threadIdx.y == 2)
				{
					d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1[  0] * a3 + dzed1[    32] * a4 + dzed1_p       * a5;
					part_ddzed1dz_p[0] = part_ddzed1dz_p[1] - dzed1[  0] * a2 - dzed1[    32] * a1 + dzed1_p       * a1 + dzed1[    96] * a2;
					part_ddzed1dz_p[1] =                                      - dzed1[    32] * a5 - dzed1_p       * a4 - dzed1[    96] * a3;

					d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1[128] * a3 + dzed1[128+32] * a4 + dzed1_q       * a5;
					part_ddzed1dz_q[0] = part_ddzed1dz_q[1] - dzed1[128] * a2 - dzed1[128+32] * a1 + dzed1_q       * a1 + dzed1[128+96] * a2;
					part_ddzed1dz_q[1] =                                      - dzed1[128+32] * a5 - dzed1_q       * a4 - dzed1[128+96] * a3;

					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a3 + dzed1[256+32] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a2 - dzed1[256+32] * a1 + dzed2_p       * a1 + dzed1[256+96] * a2;
					part_ddzed2dz_p[1] =                                      - dzed1[256+32] * a5 - dzed2_p       * a4 - dzed1[256+96] * a3;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a3 + dzed1[384+32] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a2 - dzed1[384+32] * a1 + dzed2_q       * a1 + dzed1[384+96] * a2;
					part_ddzed2dz_q[1] =                                      - dzed1[384+32] * a5 - dzed2_q       * a4 - dzed1[384+96] * a3;
				}
				else if (threadIdx.y == 3)
				{
					d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1[  0] * a2 + dzed1[    32] * a3 + dzed1[    64] * a4 + dzed1_p       * a5;
					part_ddzed1dz_p[0] = part_ddzed1dz_p[1] - dzed1[  0] * a3 - dzed1[    32] * a2 - dzed1[    64] * a1 + dzed1_p       * a1;
					part_ddzed1dz_p[1] =                                                           - dzed1[    64] * a5 - dzed1_p       * a4;

					d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1[128] * a2 + dzed1[128+32] * a3 + dzed1[128+64] * a4 + dzed1_q       * a5;
					part_ddzed1dz_q[0] = part_ddzed1dz_q[1] - dzed1[128] * a3 - dzed1[128+32] * a2 - dzed1[128+64] * a1 + dzed1_q       * a1;
					part_ddzed1dz_q[1] =                                                           - dzed1[128+64] * a5 - dzed1_q       * a4;

					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a2 + dzed1[256+32] * a3 + dzed1[256+64] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a3 - dzed1[256+32] * a2 - dzed1[256+64] * a1 + dzed2_p       * a1;
					part_ddzed2dz_p[1] =                                                           - dzed1[256+64] * a5 - dzed2_p       * a4;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a2 + dzed1[384+32] * a3 + dzed1[384+64] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a3 - dzed1[384+32] * a2 - dzed1[384+64] * a1 + dzed2_q       * a1;
					part_ddzed2dz_q[1] =                                                           - dzed1[384+64] * a5 - dzed2_q       * a4;
				}
			}
			else if (iZ < dimz-4)
			{
				if (threadIdx.y == 0)
				{
					d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1_p    * a5;
					part_ddzed1dz_p[0] = part_ddzed1dz_p[1] + dzed1_p    * a1 + dzed1[    32] * a2 + dzed1[    64] * a3 + dzed1[    96] * a4;
					part_ddzed1dz_p[1] =                                                                                - dzed1[    96] * e1;

					d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1_q    * a5;
					part_ddzed1dz_q[0] = part_ddzed1dz_q[1] + dzed1_q    * a1 + dzed1[128+32] * a2 + dzed1[128+64] * a3 + dzed1[128+96] * a4;
					part_ddzed1dz_q[1] =                                                                                - dzed1[128+96] * e1;

					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed2_p    * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] + dzed2_p    * a1 + dzed1[256+32] * a2 + dzed1[256+64] * a3 + dzed1[256+96] * a4;
					part_ddzed2dz_p[1] =                                                                                - dzed1[256+96] * e1;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed2_q    * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] + dzed2_q    * a1 + dzed1[384+32] * a2 + dzed1[384+64] * a3 + dzed1[384+96] * a4;
					part_ddzed2dz_q[1] =                                                                                - dzed1[384+96] * e1;
				}
				else if (threadIdx.y == 1)
				{
					d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1[  0] * a4 + dzed1_p       * a5;
					part_ddzed1dz_p[0] = part_ddzed1dz_p[1] - dzed1[  0] * a1 + dzed1_p       * a1 + dzed1[    64] * a2 + dzed1[    96] * a3;

					d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1[128] * a4 + dzed1_q       * a5;
					part_ddzed1dz_q[0] = part_ddzed1dz_q[1] - dzed1[128] * a1 + dzed1_q       * a1 + dzed1[128+64] * a2 + dzed1[128+96] * a3;

					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a1 + dzed2_p       * a1 + dzed1[256+64] * a2 + dzed1[256+96] * a3;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a1 + dzed2_q       * a1 + dzed1[384+64] * a2 + dzed1[384+96] * a3;
				}
				else if (threadIdx.y == 2)
				{
					d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1[  0] * a3 + dzed1[    32] * a4 + dzed1_p       * a5;
					part_ddzed1dz_p[0] = part_ddzed1dz_p[1] - dzed1[  0] * a2 - dzed1[    32] * a1 + dzed1_p       * a1 + dzed1[    96] * a2;

					d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1[128] * a3 + dzed1[128+32] * a4 + dzed1_q       * a5;
					part_ddzed1dz_q[0] = part_ddzed1dz_q[1] - dzed1[128] * a2 - dzed1[128+32] * a1 + dzed1_q       * a1 + dzed1[128+96] * a2;

					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a3 + dzed1[256+32] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a2 - dzed1[256+32] * a1 + dzed2_p       * a1 + dzed1[256+96] * a2;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a3 + dzed1[384+32] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a2 - dzed1[384+32] * a1 + dzed2_q       * a1 + dzed1[384+96] * a2;
				}
				else if (threadIdx.y == 3)
				{
					d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1[  0] * a2 + dzed1[    32] * a3 + dzed1[    64] * a4 + dzed1_p       * a5;
					part_ddzed1dz_p[0] = part_ddzed1dz_p[1] - dzed1[  0] * a3 - dzed1[    32] * a2 - dzed1[    64] * a1 + dzed1_p       * a1;

					d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1[128] * a2 + dzed1[128+32] * a3 + dzed1[128+64] * a4 + dzed1_q       * a5;
					part_ddzed1dz_q[0] = part_ddzed1dz_q[1] - dzed1[128] * a3 - dzed1[128+32] * a2 - dzed1[128+64] * a1 + dzed1_q       * a1;

					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a2 + dzed1[256+32] * a3 + dzed1[256+64] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] - dzed1[256] * a3 - dzed1[256+32] * a2 - dzed1[256+64] * a1 + dzed2_p       * a1;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a2 + dzed1[384+32] * a3 + dzed1[384+64] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] - dzed1[384] * a3 - dzed1[384+32] * a2 - dzed1[384+64] * a1 + dzed2_q       * a1;
				}
			}
			else if (iZ < dimz)
			{
				if (threadIdx.y == 0)
				{
					d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1_p    * a5;
					part_ddzed1dz_p[0] = part_ddzed1dz_p[1] + dzed1_p    * e1;

					d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1_q    * a5;
					part_ddzed1dz_q[0] = part_ddzed1dz_q[1] + dzed1_q    * e1;

					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed2_p    * a5;
					part_ddzed2dz_p[0] = part_ddzed2dz_p[1] + dzed2_p    * e1;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed2_q    * a5;
					part_ddzed2dz_q[0] = part_ddzed2dz_q[1] + dzed2_q    * e1;
				}
				else if (threadIdx.y == 1)
				{
					d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1[  0] * a4 + dzed1_p       * a5;
					part_ddzed1dz_p[0] =                    - dzed1[  0] * e1 + dzed1_p       * e1;

					d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1[128] * a4 + dzed1_q       * a5;
					part_ddzed1dz_q[0] =                    - dzed1[128] * e1 + dzed1_q       * e1;

					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] =                    - dzed1[256] * e1 + dzed2_p       * e1;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] =                    - dzed1[384] * e1 + dzed2_q       * e1;
				}
				else if (threadIdx.y == 2)
				{
					d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1[  0] * a3 + dzed1[    32] * a4 + dzed1_p       * a5;
					part_ddzed1dz_p[0] =                                      - dzed1[    32] * e1 + dzed1_p       * e1;

					d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1[128] * a3 + dzed1[128+32] * a4 + dzed1_q       * a5;
					part_ddzed1dz_q[0] =                                      - dzed1[128+32] * e1 + dzed1_q       * e1;

					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a3 + dzed1[256+32] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] =                                      - dzed1[256+32] * e1 + dzed2_p       * e1;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a3 + dzed1[384+32] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] =                                      - dzed1[384+32] * e1 + dzed2_q       * e1;
				}
				else if (threadIdx.y == 3)
				{
					d_dzed1_dz_p       = part_ddzed1dz_p[0] + dzed1[  0] * a2 + dzed1[    32] * a3 + dzed1[    64] * a4 + dzed1_p       * a5;
					part_ddzed1dz_p[0] =                                                           - dzed1[    64] * e1 + dzed1_p       * e1;

					d_dzed1_dz_q       = part_ddzed1dz_q[0] + dzed1[128] * a2 + dzed1[128+32] * a3 + dzed1[128+64] * a4 + dzed1_q       * a5;
					part_ddzed1dz_q[0] =                                                           - dzed1[128+64] * e1 + dzed1_q       * e1;

					d_dzed2_dz_p       = part_ddzed2dz_p[0] + dzed1[256] * a2 + dzed1[256+32] * a3 + dzed1[256+64] * a4 + dzed2_p       * a5;
					part_ddzed2dz_p[0] =                                                           - dzed1[256+64] * e1 + dzed2_p       * e1;

					d_dzed2_dz_q       = part_ddzed2dz_q[0] + dzed1[384] * a2 + dzed1[384+32] * a3 + dzed1[384+64] * a4 + dzed2_q       * a5;
					part_ddzed2dz_q[0] =                                                           - dzed1[384+64] * e1 + dzed2_q       * e1;
				}
			}
			else
			{
				d_dzed1_dz_p = part_ddzed1dz_p[0];
				d_dzed1_dz_q = part_ddzed1dz_q[0];
				d_dzed2_dz_p = part_ddzed2dz_p[0];
				d_dzed2_dz_q = part_ddzed2dz_q[0];
			}
		}

		__syncthreads();

		if (iZ > 0)
		{
			// complete wave equation with X and Z derivatives.
			float V4_p = d_dxed2_dx_p + d_dzed2_dz_p;
			float V4_q = d_dxed2_dx_q + d_dzed2_dz_q;

			float V5_p = d_dxed1_dx_p - d_dzed1_dz_p;
			float V5_q = d_dxed1_dx_q - d_dzed1_dz_q;

#ifdef NAN_DESU_KA
			if (!fu1 && (d_dxed1_dx_p != 0.0f || d_dxed1_dx_q != 0.0f || d_dzed1_dz_p != 0.0f || d_dzed1_dz_q != 0.0f || d_dxed2_dx_p != 0.0f || d_dxed2_dx_q != 0.0f || d_dzed2_dz_p != 0.0f || d_dzed2_dz_q != 0.0f))
			{
				fu1 = 1;
				printf("Compute_DYED1andDYED2... ; threadIdx=[%d,%d,%d] blockIdx=[%d,%d] :: d_dxed1_dx_p=%e, d_dxed1_dx_q=%e, d_dzed1_dz_p=%e, d_dzed1_dz_q=%e, d_dxed2_dx_p=%e, d_dxed2_dx_q=%e, d_dzed2_dz_p=%e, d_dzed2_dz_q=%e\n",threadIdx.x,threadIdx.y,threadIdx.z,blockIdx.x,blockIdx.y,d_dxed1_dx_p,d_dxed1_dx_q,d_dzed1_dz_p,d_dzed1_dz_q,d_dxed2_dx_p,d_dxed2_dx_q,d_dzed2_dz_p,d_dzed2_dz_q);
			}
#endif

			// write out V4 and V5 for dx and dz
			float2* dy = ((float2*)d_dy) + thr_y * stride_y_dy + (iZ-4) * stride_z_dy;

			float2 f2_V4;
			f2_V4.x = V4_p;
			f2_V4.y = V4_q;
			dy[threadIdx.y*64+threadIdx.x] = f2_V4;

			float2 f2_V5;
			f2_V5.x = V5_p;
			f2_V5.y = V5_q;
			dy[threadIdx.y*64+threadIdx.x+32] = f2_V5;
		}
			
		// this is necessary because next iteration will overwrite the part of shared memory that contains DZED1 and DZED2.
		__syncthreads();
	}
}

__global__
#if __CUDA_ARCH__ >= 300
__launch_bounds__(1024)
#elif __CUDA_ARCH__ >= 130
__launch_bounds__(640)
#endif
void Compute_DYED1and2_Part_V4_V5(
	float* d_dx,
	float* d_dy,
	float* d_dz,
	const int* d_em0,
	const int* d_em1,
	const int* d_em2,
	const int* d_em3,
	int dimy,
	int dimz,
	const float dt,
        const float a1,
        const float a2,
        const float a3,
        const float a4,
        const float a5,
        const float e1,
        const float a1h,
        const float a2h,
        const float a3h,
        const float a4h,
        const float a5h,
        float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
        float Dip_min,
        float Dip_scaler,
        float Azm_min,
        float Azm_scaler,
        float C44C33_min,
        float C44C33_scaler,
        float Vel_min,
        float Vel_scaler,
        float Del_min,
        float Del_scaler,
        float Eps_min,
        float Eps_scaler
	)
{
	__shared__ float dx_dy_dz_DenAng[1536];
	/*
	int nn = blockDim.x * blockDim.y * blockDim.z;
	for (int i = 0;  i < 1536;  i+=nn)
	{
		int idx = i+threadIdx.x+threadIdx.y*blockDim.x+threadIdx.z*blockDim.y*blockDim.x;
		if (idx < 1536) dx_dy_dz_DenAng[idx] = 0.0f;
	}
	__syncthreads();
	*/

	if (blockIdx.y < 5 || blockIdx.y >= dimy-13)
	{
		Compute_DYED1and2_Part_V4_V5_YHalo(
				dx_dy_dz_DenAng,
				d_dx,d_dy,d_dz,d_em0,d_em1,d_em2,d_em3,dimz,
				dt,a1,a2,a3,a4,a5,e1,a1h,a2h,a3h,a4h,a5h,
				inv_Q_min,inv_Q_scaler,Den_min,Den_scaler,Dip_min,Dip_scaler,Azm_min,Azm_scaler,C44C33_min,C44C33_scaler,Vel_min,Vel_scaler,Del_min,Del_scaler,Eps_min,Eps_scaler);
	}
	else
	{
		Compute_DYED1and2_Part_V4_V5_Main(
				dx_dy_dz_DenAng,
				d_dx,d_dy,d_dz,d_em0,d_em1,d_em2,d_em3,dimz,
				dt,a1,a2,a3,a4,a5,e1,a1h,a2h,a3h,a4h,a5h,
				inv_Q_min,inv_Q_scaler,Den_min,Den_scaler,Dip_min,Dip_scaler,Azm_min,Azm_scaler,C44C33_min,C44C33_scaler,Vel_min,Vel_scaler,Del_min,Del_scaler,Eps_min,Eps_scaler);
	}
}

/*

STAGE 3 ONLY NEEDS TO LOOP OVER blockIdx.y=[4,dimz-5]

*/

//
// Block shape is 32 by 4 by 1.
// Grid shape is 1 by dimz by 1.
// Loop over Y.
// Compute d(dyed1)dy, d(dyed2)dy, compute deltas.
// Read deltas from d(dxed*)dx and d(dzed*)dz.
// Read p-q and previous p-q.
// Compute next p-q.
// Write next p-q.
// 
__global__
#if __CUDA_ARCH__ >= 300
__launch_bounds__(1152)
#elif __CUDA_ARCH__ >= 130
__launch_bounds__(896)
#endif
void Compute_VTI_Next_PQ_T2(
	const float* spg_x,		// ptr adjusted for threadIdx.x == 0
	const float* spg_y,		// ptr adjusted for threadIdx.y == 0
	const float* spg_z,
	const float* d_dyed1and2,  	// dyed1 and dyed2 with no yoffset
	const float* d_V4_V5,  		// V4 and V5 for dx and dz with no offset
	const float* d_pq1_yoff9,	// current pq with yoffset=9
	const float* d_pq2_yoff9,	// 
	const float* d_prev_pq1_yoff9,	// previous pq with yoffset=9
	const float* d_prev_pq2_yoff9,	//
	float* d_next_pq1,		// in: deltas for ddx and ddz. out: next pq
	float* d_next_pq2,		// NB! Apply no yoffset.
	const int* d_em1_yoff9,		// Earth model with yoffset=9
	const int* d_em2_yoff9,		//
	int dimy,
	int dimz,
	const float dt,
        const float a1h,
        const float a2h,
        const float a3h,
        const float a4h,
        const float a5h,
        float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
        float C44C33_min,
        float C44C33_scaler,
        float Vel_min,
        float Vel_scaler,
        float Del_min,
        float Del_scaler,
        float Eps_min,
        float Eps_scaler
	)
{
	const int stride_z_dy = 64;
	const int stride_y_dy = stride_z_dy * dimz;

	int thr_z = blockIdx.y;

	float spgxz1 = spg_x[(threadIdx.x>>1)] * spg_z[thr_z];
	float spgxz2 = spg_x[16+(threadIdx.x>>1)] * spg_z[thr_z];
	
	int rb_idx = 256 + threadIdx.y * 64;

	__shared__ float rb[896];

 	float2* dyed1and2 = ((float2*)d_dyed1and2) + (thr_z>>2) * 4 * stride_z_dy + (thr_z&3) * 64 + (threadIdx.y + 4) * stride_y_dy;
	float2* V4_V5 = ((float2*)d_V4_V5) + (thr_z>>2) * 4 * stride_z_dy + (thr_z&3) * 64 + threadIdx.y * stride_y_dy;

	float part_d_dyed1_dy_p[2] = {0.0f, 0.0f};
	float part_d_dyed1_dy_q[2] = {0.0f, 0.0f};

	// load y=-5.
	if (threadIdx.y == 0)
	{
		float2 v = dyed1and2[threadIdx.x];
		rb[192   +threadIdx.x] = v.x;
		rb[192+32+threadIdx.x] = v.y;
	}
	dyed1and2 += stride_y_dy;
	V4_V5 += stride_y_dy;
	//__syncthreads();  // not needed

	for (int iY = -2;  iY < (dimy+3)/4;  ++iY)
	{
		int thr_y = iY * 4 + threadIdx.y;
		if (thr_y >= dimy) break;

		/*
		if ((iY&31) == 0)
		{
			// load 128 spg_y values
			int iY4 = iY * 4;
			int thr_idx = threadIdx.x + threadIdx.y * 32;
			if (iY4 + thr_idx < dimy)
			{
				rb[1280+thr_idx] = spg_y[iY4+thr_idx];
			}
		}
		*/

		// compute d(dyed1)dy and d(dyed2)dy
		float2 v = dyed1and2[threadIdx.x];
		rb[   rb_idx+threadIdx.x] = v.x;
		rb[32+rb_idx+threadIdx.x] = v.y;
		rb_idx = (rb_idx + 256) & 511;
		dyed1and2 += stride_y_dy * 4;

		__syncthreads();

		int rb_comp_idx = (rb_idx + 64) & 511;  // rb_comp_idx = rb_idx - 7 ...

		int rb_comp_idx2 = (rb_comp_idx +  64) & 511;
		int rb_comp_idx3 = (rb_comp_idx + 128) & 511;
		int rb_comp_idx4 = (rb_comp_idx + 192) & 511;

		float s0 = rb[rb_comp_idx  + threadIdx.x];
		float s1 = rb[rb_comp_idx2 + threadIdx.x];
		float s2 = rb[rb_comp_idx3 + threadIdx.x];
		float s3 = rb[rb_comp_idx4 + threadIdx.x];

		float d_dyed1_dy_p   = part_d_dyed1_dy_p[0] + s3 * a5h + s2 * a4h + s1 * a3h + s0 * a2h;
		part_d_dyed1_dy_p[0] = part_d_dyed1_dy_p[1] + s3 * a1h - s2 * a1h - s1 * a2h - s0 * a3h;
		part_d_dyed1_dy_p[1] =                      - s3 * a4h - s2 * a5h;

		s0 = rb[32 + rb_comp_idx  + threadIdx.x];
		s1 = rb[32 + rb_comp_idx2 + threadIdx.x];
		s2 = rb[32 + rb_comp_idx3 + threadIdx.x];
		s3 = rb[32 + rb_comp_idx4 + threadIdx.x];

		float d_dyed1_dy_q   = part_d_dyed1_dy_q[0] + s3 * a5h + s2 * a4h + s1 * a3h + s0 * a2h;
		part_d_dyed1_dy_q[0] = part_d_dyed1_dy_q[1] + s3 * a1h - s2 * a1h - s1 * a2h - s0 * a3h;
		part_d_dyed1_dy_q[1] =                      - s3 * a4h - s2 * a5h;

		// compute next p-q after lead-in
		if (iY >= 0)
		{
			//int thr_y = iY * 4 + threadIdx.y;

			// complete wave equation with Y derivatives.
			float2 f2_V4 = V4_V5[threadIdx.x];
			float V4_p = f2_V4.x;
			float V4_q = f2_V4.y;

			float2 f2_V5 = V4_V5[threadIdx.x+32];
			float V5_p = d_dyed1_dy_p + f2_V5.x;
			float V5_q = d_dyed1_dy_q + f2_V5.y;

			// ..load earth model for the affected cells. 
			int stride_z_em = 16;
			int stride_y_em = stride_z_em * dimz;

			int idx3 = thr_z * stride_z_em + (((thr_y) * 2) + (threadIdx.x >> 4)) * stride_y_em + (threadIdx.x & 15);
			int lidx3 = threadIdx.y * 32 + (threadIdx.x >> 4) * 128 + (threadIdx.x & 15);

			((int*)rb)[512+lidx3] = d_em1_yoff9[idx3];
			((int*)rb)[512+lidx3+16] = d_em2_yoff9[idx3];
			__syncthreads();
			// rb+1024[  0..127] <= DenAng [0,0] to DenAng [31,3]
			// rb+1024[128..255] <= VelAnis[0,0] to VelAnis[31,3]i

			int idx = threadIdx.x + threadIdx.y * 32;
			float inv_Q, Den, C44C33, Vel, Del, Eps;
			Decode_DenAng_And_VelAnis(((int*)rb)[512+idx],((int*)rb)[512+128+idx],inv_Q,Den,C44C33,Vel,Del,Eps,inv_Q_min,inv_Q_scaler,Den_min,Den_scaler,C44C33_min,C44C33_scaler,Vel_min,Vel_scaler,Del_min,Del_scaler,Eps_min,Eps_scaler);

			float Vp2 = Vel * Vel * dt * dt;
			float C33 = Den * Vp2;
			float C44 = C33 * C44C33;
			float C33mC44 = C33 - C44;
			float C13pC44 = __fsqrt_rn((2.0f * Del * C33 + C33mC44) * C33mC44);
			float C66 = (1.0f + 2.0f * Eps) * C33;

			double delta_p = (double)(C66 * V5_p) + (double)(C44 * V4_p) + (double)(C13pC44 * V4_q);
			double delta_q = (double)(C44 * V5_q) + (double)(C33 * V4_q) + (double)(C13pC44 * V5_p);

			// shuffle deltas through shared memory
			rb[512+idx    ] = (float)delta_p;
			rb[512+idx+128] = (float)delta_q;
			rb[512+idx+256] = inv_Q;
			__syncthreads();
			float delta1 = rb[512+(threadIdx.x&1)*128+(threadIdx.x>>1)+threadIdx.y*32];
			float delta2 = rb[512+16+(threadIdx.x&1)*128+(threadIdx.x>>1)+threadIdx.y*32];
			float inv_Q_1 = rb[768+(threadIdx.x>>1)+threadIdx.y*32];
			float inv_Q_2 = rb[768+16+(threadIdx.x>>1)+threadIdx.y*32];

			const int stride_z = 32;
			const int stride_y = stride_z * dimz;

			int idx4 = thr_y*stride_y + thr_z * stride_z + threadIdx.x;

			//float spgxyz = spgxz * rb[1280+(thr_y&127)];
			float spgxyz1 = spgxz1 * spg_y[thr_y];
			float spgxyz2 = spgxz2 * spg_y[thr_y];

			// add prev, curr pq
			float curr_pq1 = d_pq1_yoff9[idx4];
			float curr_pq2 = d_pq2_yoff9[idx4];
			float prev_pq1 = d_prev_pq1_yoff9[idx4];
			float prev_pq2 = d_prev_pq2_yoff9[idx4];
	
			//float npq1 = delta1 + 2.0f * curr_pq1;
			//npq1 = npq1 - spgxyz1 * prev_pq1;
			//npq1 = inv_Q * spgxyz1 * npq1;
			double npq1 = (double)(2.0f * curr_pq1) - (double)(spgxyz1 * prev_pq1);
			npq1 = npq1 + (double)delta1;
			npq1 = (double)(inv_Q_1 * spgxyz1) * npq1;

			//float npq2 = delta2 + 2.0f * curr_pq2;
			//npq2 = npq2 - spgxyz2 * prev_pq2;
			//npq2 = inv_Q * spgxyz2 * npq2;
			double npq2 = (double)(2.0f * curr_pq2) - (double)(spgxyz2 * prev_pq2);
                        npq2 = npq2 + (double)delta2;
                        npq2 = (double)(inv_Q_2 * spgxyz2) * npq2;

			//float npq1 = 2.0f * curr_pq1 - spgxyz * prev_pq1;
			//float npq2 = 2.0f * curr_pq2 - spgxyz * prev_pq2;

			// sponge and Q
			//npq1 = inv_Q * spgxyz * ( npq1 + delta1 );
			//npq2 = inv_Q * spgxyz * ( npq2 + delta2 );
#ifdef NAN_DESU_KA
			if (!fu1 && (npq1 != 0.0f || npq2 != 0.0f))
			{
				fu1 = 1;
				printf("threadIdx=[%d,%d,%d] blockIdx=[%d,%d] :: npq1=%e, npq2=%e, curr_pq1=%e, curr_pq2=%e, prev_pq1=%e, prev_pq2=%e, xy_delta1=%e, xy_delta2=%e, d_dyed1_dy_p=%e, d_dyed1_dy_q=%e, d_dyed2_dy_p=%e, d_dyed2_dy_q=%e\n",threadIdx.x,threadIdx.y,threadIdx.z,blockIdx.x,blockIdx.y,npq1,npq2,curr_pq1,curr_pq2,prev_pq1,prev_pq2,xy_delta1,xy_delta2,d_dyed1_dy_p,d_dyed1_dy_q,d_dyed2_dy_p,d_dyed2_dy_q);
			}
#endif

			// write out results
			d_next_pq1[idx4] = (float)npq1;
			d_next_pq2[idx4] = (float)npq2;
		}
		V4_V5 += stride_y_dy * 4;
	}
}

//
// Block shape is 32 by 4 by 1.
// Grid shape is 1 by dimz by 1.
// Loop over Y.
// Compute d(dyed1)dy, d(dyed2)dy, compute deltas.
// Read deltas from d(dxed*)dx and d(dzed*)dz.
// Read p-q and previous p-q.
// Compute next p-q.
// Write next p-q.
// 
__global__
#if __CUDA_ARCH__ >= 300
__launch_bounds__(1152)
#elif __CUDA_ARCH__ >= 130
__launch_bounds__(896)
#endif
void Compute_Next_PQ_T2(
	const float* spg_x,		// ptr adjusted for threadIdx.x == 0
	const float* spg_y,		// ptr adjusted for threadIdx.y == 0
	const float* spg_z,
	const float* d_dyed1and2,  	// dyed1 and dyed2 with no yoffset
	const float* d_V4_V5,  		// V4 and V5 for dx and dz with no offset
	const float* d_pq1_yoff9,	// current pq with yoffset=9
	const float* d_pq2_yoff9,	// 
	const float* d_prev_pq1_yoff9,	// previous pq with yoffset=9
	const float* d_prev_pq2_yoff9,	//
	float* d_next_pq1,		// in: deltas for ddx and ddz. out: next pq
	float* d_next_pq2,		// NB! Apply no yoffset.
	const int* d_em1_yoff9,		// Earth model with yoffset=9
	const int* d_em2_yoff9,		//
	int dimy,
	int dimz,
	const float dt,
        const float a1h,
        const float a2h,
        const float a3h,
        const float a4h,
        const float a5h,
        float inv_Q_min,
        float inv_Q_scaler,
        float Den_min,
        float Den_scaler,
        float Dip_min,
        float Dip_scaler,
        float Azm_min,
        float Azm_scaler,
        float C44C33_min,
        float C44C33_scaler,
        float Vel_min,
        float Vel_scaler,
        float Del_min,
        float Del_scaler,
        float Eps_min,
        float Eps_scaler
	)
{
	const int stride_z_dy = 64;
	const int stride_y_dy = stride_z_dy * dimz;

	int thr_z = blockIdx.y;

	float spgxz1 = spg_x[(threadIdx.x>>1)] * spg_z[thr_z];
	float spgxz2 = spg_x[16+(threadIdx.x>>1)] * spg_z[thr_z];
	
	int rb_idx = 512 + threadIdx.y * 128;

	__shared__ float rb[1408];

 	float2* dyed1and2 = ((float2*)d_dyed1and2) + (thr_z>>2) * 4 * stride_z_dy + (thr_z&3) * 64 + (threadIdx.y + 4) * stride_y_dy;
	float2* V4_V5 = ((float2*)d_V4_V5) + (thr_z>>2) * 4 * stride_z_dy + (thr_z&3) * 64 + threadIdx.y * stride_y_dy;

	float part_d_dyed1_dy_p[2] = {0.0f, 0.0f};
	float part_d_dyed1_dy_q[2] = {0.0f, 0.0f};
	float part_d_dyed2_dy_p[2] = {0.0f, 0.0f};
	float part_d_dyed2_dy_q[2] = {0.0f, 0.0f};

	// load y=-5.
	if (threadIdx.y == 0)
	{
		float2 v = dyed1and2[threadIdx.x];
		rb[384   +threadIdx.x] = v.x;
		rb[384+32+threadIdx.x] = v.y;
		v = dyed1and2[threadIdx.x+32];
		rb[384+64+threadIdx.x] = v.x;
		rb[384+96+threadIdx.x] = v.y;
	}
	dyed1and2 += stride_y_dy;
	V4_V5 += stride_y_dy;
	//__syncthreads();  // not needed

#ifdef NAN_DESU_KA
	int fu1 = 0;
#endif
	
	for (int iY = -2;  iY < (dimy+3)/4;  ++iY)
	{
		int thr_y = iY * 4 + threadIdx.y;
		if (thr_y >= dimy) break;

		/*
		if ((iY&31) == 0)
		{
			// load 128 spg_y values
			int iY4 = iY * 4;
			int thr_idx = threadIdx.x + threadIdx.y * 32;
			if (iY4 + thr_idx < dimy)
			{
				rb[1280+thr_idx] = spg_y[iY4+thr_idx];
			}
		}
		*/

		// compute d(dyed1)dy and d(dyed2)dy
		float2 v = dyed1and2[threadIdx.x];
		rb[   rb_idx+threadIdx.x] = v.x;
		rb[32+rb_idx+threadIdx.x] = v.y;
		v = dyed1and2[threadIdx.x+32];
		rb[64+rb_idx+threadIdx.x] = v.x;
		rb[96+rb_idx+threadIdx.x] = v.y;
		rb_idx = (rb_idx + 512) & 1023;
		dyed1and2 += stride_y_dy * 4;

		__syncthreads();

		int rb_comp_idx = (rb_idx + 128) & 1023;  // rb_comp_idx = rb_idx - 7 ...

		int rb_comp_idx2 = (rb_comp_idx + 128) & 1023;
		int rb_comp_idx3 = (rb_comp_idx + 256) & 1023;
		int rb_comp_idx4 = (rb_comp_idx + 384) & 1023;

		float s0 = rb[rb_comp_idx  + threadIdx.x];
		float s1 = rb[rb_comp_idx2 + threadIdx.x];
		float s2 = rb[rb_comp_idx3 + threadIdx.x];
		float s3 = rb[rb_comp_idx4 + threadIdx.x];

		float d_dyed1_dy_p   = part_d_dyed1_dy_p[0] + s3 * a5h + s2 * a4h + s1 * a3h + s0 * a2h;
		part_d_dyed1_dy_p[0] = part_d_dyed1_dy_p[1] + s3 * a1h - s2 * a1h - s1 * a2h - s0 * a3h;
		part_d_dyed1_dy_p[1] =                      - s3 * a4h - s2 * a5h;

		s0 = rb[32 + rb_comp_idx  + threadIdx.x];
		s1 = rb[32 + rb_comp_idx2 + threadIdx.x];
		s2 = rb[32 + rb_comp_idx3 + threadIdx.x];
		s3 = rb[32 + rb_comp_idx4 + threadIdx.x];

		float d_dyed1_dy_q   = part_d_dyed1_dy_q[0] + s3 * a5h + s2 * a4h + s1 * a3h + s0 * a2h;
		part_d_dyed1_dy_q[0] = part_d_dyed1_dy_q[1] + s3 * a1h - s2 * a1h - s1 * a2h - s0 * a3h;
		part_d_dyed1_dy_q[1] =                      - s3 * a4h - s2 * a5h;

		s0 = rb[64 + rb_comp_idx  + threadIdx.x];
		s1 = rb[64 + rb_comp_idx2 + threadIdx.x];
		s2 = rb[64 + rb_comp_idx3 + threadIdx.x];
		s3 = rb[64 + rb_comp_idx4 + threadIdx.x];

		float d_dyed2_dy_p   = part_d_dyed2_dy_p[0] + s3 * a5h + s2 * a4h + s1 * a3h + s0 * a2h;
		part_d_dyed2_dy_p[0] = part_d_dyed2_dy_p[1] + s3 * a1h - s2 * a1h - s1 * a2h - s0 * a3h;
		part_d_dyed2_dy_p[1] =                      - s3 * a4h - s2 * a5h;

		s0 = rb[96 + rb_comp_idx  + threadIdx.x];
		s1 = rb[96 + rb_comp_idx2 + threadIdx.x];
		s2 = rb[96 + rb_comp_idx3 + threadIdx.x];
		s3 = rb[96 + rb_comp_idx4 + threadIdx.x];

		float d_dyed2_dy_q   = part_d_dyed2_dy_q[0] + s3 * a5h + s2 * a4h + s1 * a3h + s0 * a2h;
		part_d_dyed2_dy_q[0] = part_d_dyed2_dy_q[1] + s3 * a1h - s2 * a1h - s1 * a2h - s0 * a3h;
		part_d_dyed2_dy_q[1] =                      - s3 * a4h - s2 * a5h;

		// compute next p-q after lead-in
		if (iY >= 0)
		{
			//int thr_y = iY * 4 + threadIdx.y;

			// complete wave equation with Y derivatives.
			float2 f2_V4 = V4_V5[threadIdx.x];
			float V4_p = d_dyed2_dy_p + f2_V4.x;
			float V4_q = d_dyed2_dy_q + f2_V4.y;

			float2 f2_V5 = V4_V5[threadIdx.x+32];
			float V5_p = d_dyed1_dy_p + f2_V5.x;
			float V5_q = d_dyed1_dy_q + f2_V5.y;

			// ..load earth model for the affected cells. 
			int stride_z_em = 16;
			int stride_y_em = stride_z_em * dimz;

			int idx3 = thr_z * stride_z_em + (((thr_y) * 2) + (threadIdx.x >> 4)) * stride_y_em + (threadIdx.x & 15);
			int lidx3 = threadIdx.y * 32 + (threadIdx.x >> 4) * 128 + (threadIdx.x & 15);

			((int*)rb)[1024+lidx3] = d_em1_yoff9[idx3];
			((int*)rb)[1024+lidx3+16] = d_em2_yoff9[idx3];
			__syncthreads();
			// rb+1024[  0..127] <= DenAng [0,0] to DenAng [31,3]
			// rb+1024[128..255] <= VelAnis[0,0] to VelAnis[31,3]i

			int idx = threadIdx.x + threadIdx.y * 32;
			float inv_Q, Den, C44C33, Vel, Del, Eps;
			Decode_DenAng_And_VelAnis(((int*)rb)[1024+idx],((int*)rb)[1024+128+idx],inv_Q,Den,C44C33,Vel,Del,Eps,inv_Q_min,inv_Q_scaler,Den_min,Den_scaler,C44C33_min,C44C33_scaler,Vel_min,Vel_scaler,Del_min,Del_scaler,Eps_min,Eps_scaler);

			float Vp2 = Vel * Vel * dt * dt;
			float C33 = Den * Vp2;
			float C44 = C33 * C44C33;
			float C33mC44 = C33 - C44;
			float C13pC44 = __fsqrt_rn((2.0f * Del * C33 + C33mC44) * C33mC44);
			float C66 = (1.0f + 2.0f * Eps) * C33;

			double delta_p = (double)(C66 * V5_p) + (double)(C44 * V4_p) + (double)(C13pC44 * V4_q);
			double delta_q = (double)(C44 * V5_q) + (double)(C33 * V4_q) + (double)(C13pC44 * V5_p);

			// shuffle deltas through shared memory
			rb[1024+idx    ] = (float)delta_p;
			rb[1024+idx+128] = (float)delta_q;
			rb[1024+idx+256] = inv_Q;
			__syncthreads();
			float delta1 = rb[1024+(threadIdx.x&1)*128+(threadIdx.x>>1)+threadIdx.y*32];
			float delta2 = rb[1024+16+(threadIdx.x&1)*128+(threadIdx.x>>1)+threadIdx.y*32];
			float inv_Q_1 = rb[1280+(threadIdx.x>>1)+threadIdx.y*32];
			float inv_Q_2 = rb[1280+16+(threadIdx.x>>1)+threadIdx.y*32];

			const int stride_z = 32;
			const int stride_y = stride_z * dimz;

			int idx4 = thr_y*stride_y + thr_z * stride_z + threadIdx.x;

			//float spgxyz = spgxz * rb[1280+(thr_y&127)];
			float spgxyz1 = spgxz1 * spg_y[thr_y];
			float spgxyz2 = spgxz2 * spg_y[thr_y];

			// add prev, curr pq
			float curr_pq1 = d_pq1_yoff9[idx4];
			float curr_pq2 = d_pq2_yoff9[idx4];
			float prev_pq1 = d_prev_pq1_yoff9[idx4];
			float prev_pq2 = d_prev_pq2_yoff9[idx4];
	
			//float npq1 = delta1 + 2.0f * curr_pq1;
			//npq1 = npq1 - spgxyz1 * prev_pq1;
			//npq1 = inv_Q * spgxyz1 * npq1;
			double npq1 = (double)(2.0f * curr_pq1) - (double)(spgxyz1 * prev_pq1);
			npq1 = npq1 + (double)delta1;
			npq1 = (double)(inv_Q_1 * spgxyz1) * npq1;

			//float npq2 = delta2 + 2.0f * curr_pq2;
			//npq2 = npq2 - spgxyz2 * prev_pq2;
			//npq2 = inv_Q * spgxyz2 * npq2;
			double npq2 = (double)(2.0f * curr_pq2) - (double)(spgxyz2 * prev_pq2);
                        npq2 = npq2 + (double)delta2;
                        npq2 = (double)(inv_Q_2 * spgxyz2) * npq2;

			// TMJ 09/05/13
			// invQ should be applied before shuffle, or should be shuffled.

			//float npq1 = 2.0f * curr_pq1 - spgxyz * prev_pq1;
			//float npq2 = 2.0f * curr_pq2 - spgxyz * prev_pq2;

			// sponge and Q
			//npq1 = inv_Q * spgxyz * ( npq1 + delta1 );
			//npq2 = inv_Q * spgxyz * ( npq2 + delta2 );
#ifdef NAN_DESU_KA
			if (!fu1 && (npq1 != 0.0f || npq2 != 0.0f))
			{
				fu1 = 1;
				printf("threadIdx=[%d,%d,%d] blockIdx=[%d,%d] :: npq1=%e, npq2=%e, curr_pq1=%e, curr_pq2=%e, prev_pq1=%e, prev_pq2=%e, xy_delta1=%e, xy_delta2=%e, d_dyed1_dy_p=%e, d_dyed1_dy_q=%e, d_dyed2_dy_p=%e, d_dyed2_dy_q=%e\n",threadIdx.x,threadIdx.y,threadIdx.z,blockIdx.x,blockIdx.y,npq1,npq2,curr_pq1,curr_pq2,prev_pq1,prev_pq2,xy_delta1,xy_delta2,d_dyed1_dy_p,d_dyed1_dy_q,d_dyed2_dy_p,d_dyed2_dy_q);
			}
#endif

			// write out results
			d_next_pq1[idx4] = (float)npq1;
			d_next_pq2[idx4] = (float)npq2;
		}
		V4_V5 += stride_y_dy * 4;
	}
}

//
// Inject source at the computed src_idx and (optionally) ghost index. 
// Only one thread is needed since the source is always injected at a single point.
// As you can imagine, this kernel finishes almost before it starts.
//
__global__ 
void Inject_Source(
	float* d_pq1,
	int src_idx1,
	int ghost_idx1,
	float src_val1,
	float* d_pq2,
	int src_idx2,
	int ghost_idx2,
	float src_val2
	)
{
	if (threadIdx.x == 0)
	{
		if (src_val1 != 0.0f)
		{
			d_pq1[src_idx1  ] += src_val1;
			d_pq1[src_idx1+1] += src_val1;
			if (ghost_idx1 >= 0)
			{
				d_pq1[ghost_idx1  ] -= src_val1;
				d_pq1[ghost_idx1+1] -= src_val1;
			}
		}
		if (src_val2 != 0.0f)
		{
			d_pq2[src_idx2  ] += src_val2;
			d_pq2[src_idx2+1] += src_val2;
			if (ghost_idx2 >= 0)
			{
				d_pq2[ghost_idx2  ] -= src_val2;
				d_pq2[ghost_idx2+1] -= src_val2;
			}
		}
	}
}

__global__
void Extract_Receiver_Values(
	float tfrac,
	int* d_rcx_1, 
	int* d_rcy_1, 
	int* d_rcz_1, 
	float* d_out_1, 
	int recnum_1,
	const float* d_pq,  // pq and rs are aligned in Y so that same index in both will return same cell
	const float* d_rs,
	const float* d_spg_x,
	const float* d_spg_y,
	const float* d_spg_z,
	int x0,
	int y0,
	int dimx,
	int dimy,
	int dimz,
	int recghost_Flag,
	int z_freesurface
	)
{
	long nthreads = blockDim.y * blockDim.x;
	long threadnum = threadIdx.y * blockDim.x + threadIdx.x;
	long rcidx = blockIdx.x * nthreads + threadnum;
	
	if (rcidx < recnum_1)
	{
		long iX = d_rcx_1[rcidx] - x0;
		long iY = d_rcy_1[rcidx] - y0;
		long iZ = d_rcz_1[rcidx];

		float spgfac = d_spg_x[iX] * d_spg_y[iY] * d_spg_z[iZ];

		long stride_z = dimx;
		long stride_y = stride_z * dimz;
		long pqidx = iZ * stride_z + iY * stride_y + iX;	
	
		float2 vpq = ((float2*)d_pq)[pqidx];
		float p = vpq.x;
		float q = vpq.y;

		float2 vrs = ((float2*)d_rs)[pqidx];
		float r = vrs.x;
		float s = vrs.y;

		float val = tfrac * spgfac * ( 2.0f * p + q ) + (1.0f - tfrac) * ( 2.0f * r + s);
		//printf("stripe %d :: iX=%ld, iY=%ld, iZ=%ld, spgfac=%f, pqidx=%ld, rcidx=%ld, val=%e\n",x0/16,iX,iY,iZ,spgfac,pqidx,rcidx,val);

		if (recghost_Flag)
		{
			long iZ_ghost = 2 * z_freesurface - iZ;
			long pqidx_ghost = iZ_ghost * stride_z + iY * stride_y + iX;

			float2 vpq_ghost = ((float2*)d_pq)[pqidx_ghost];
			float p_ghost = vpq_ghost.x;
			float q_ghost = vpq_ghost.y;

			float2 vrs_ghost = ((float2*)d_rs)[pqidx_ghost];
			float r_ghost = vrs_ghost.x;
			float s_ghost = vrs_ghost.y;

			float val_ghost = tfrac * spgfac * ( 2.0f * p_ghost + q_ghost ) + (1.0f - tfrac) * ( 2.0f * r_ghost + s_ghost);

			val = val - val_ghost;
		}

		d_out_1[rcidx] = val;
	}
}

void 
Host_Inject_Source(
	cudaStream_t compute_stream,
	float src_val_1,
	int stripe_xsrc_1,
	int stripe_ysrc_1,
	int stripe_zsrc_1,
	int stripe_zsrcghost_1,
	float* d_curr_pq1,
	float src_val_2,
	int stripe_xsrc_2,
	int stripe_ysrc_2,
	int stripe_zsrc_2,
	int stripe_zsrcghost_2,
	float* d_curr_pq2,
	int stride_y,
	int stride_z,
	int curr_yoff9_idx
	)
{
	if (src_val_1 != 0.0f || src_val_2 != 0.0f)
	{
		int src_idx_1 = stripe_xsrc_1 * 2 + stripe_ysrc_1 * stride_y + stripe_zsrc_1 * stride_z;
		int ghost_idx_1 = stripe_zsrcghost_1 < 0 ? -1 : stripe_xsrc_1 * 2 + stripe_ysrc_1 * stride_y + stripe_zsrcghost_1 * stride_z;

		int src_idx_2 = stripe_xsrc_2 * 2 + stripe_ysrc_2 * stride_y + stripe_zsrc_2 * stride_z;
		int ghost_idx_2 = stripe_zsrcghost_2 < 0 ? -1 : stripe_xsrc_2 * 2 + stripe_ysrc_2 * stride_y + stripe_zsrcghost_2 * stride_z;
		
		dim3 blockShape0(32,1,1);
		dim3 gridShape0(1,1);
		Inject_Source<<<gridShape0,blockShape0,0,compute_stream>>>(
			d_curr_pq1+curr_yoff9_idx,src_idx_1,ghost_idx_1,src_val_1,
			d_curr_pq2+curr_yoff9_idx,src_idx_2,ghost_idx_2,src_val_2
			);
	}
}

/*
** Advance P-Q wavefields one timestep.
*/
void 
TTIDenQ_T2_GPU_Timestep(
	int Kernel_Type,
	cudaStream_t compute_stream,
	const float* d_prev_pq1,
	const float* d_prev_pq2,
	const float* d_curr_pq0,
	const float* d_curr_pq1,
	const float* d_curr_pq2,
	const float* d_curr_pq3,
	float* d_next_pq1,
	float* d_next_pq2,
	const int* d_em0,
	const int* d_em1,
	const int* d_em2,
	const int* d_em3,
	const float* d_spg_x,
	const float* d_spg_y,
	const float* d_spg_z,
	const int prev_y_offset,
	const int curr_y_offset,
	const int next_y_offset,
	const int em_y_offset,
	float* d_dx,
	float* d_dy,
	float* d_dz,
	int x0_1,
	int x0_2,
	int y0,
	int dimy,
	int dimz,
	const float dt,
	const float a1h,
	const float a2h,
	const float a3h,
	const float a4h,
	const float a5h,
	const float a1z,
	const float a2z,
	const float a3z,
	const float a4z,
	const float a5z,
	const float e1z,
	const float inv_Q_min,
	const float inv_Q_scaler,
	const float Den_min,
	const float Den_scaler,
	const float Dip_min,
	const float Dip_scaler,
	const float Azm_min,
	const float Azm_scaler,
	const float C44C33_min,
	const float C44C33_scaler,
	const float Vel_min,
	const float Vel_scaler,
	const float Del_min,
	const float Del_scaler,
	const float Eps_min,
	const float Eps_scaler,
	int stripe_xsrc_1, 
	int stripe_ysrc_1,
	int stripe_zsrc_1,
	int stripe_zsrcghost_1,
	float src_val_1,
	int stripe_xsrc_2, 
	int stripe_ysrc_2,
	int stripe_zsrc_2,
	int stripe_zsrcghost_2,
	float src_val_2,
	float tfrac1,
	int* d_rcx_1,
	int* d_rcy_1,
	int* d_rcz_1,
	float* d_out_1,
	int recnum_1,
	float tfrac2,
	int* d_rcx_2,
	int* d_rcy_2,
	int* d_rcz_2,
	float* d_out_2,
	int recnum_2,
	int recghost_Flag,
	int z_freesurface
	)
{
#ifdef NAN_DESU_KA
	printf("prev_y_offset = %d, curr_y_offset = %d, next_y_offset = %d, em_y_offset = %d\n",prev_y_offset,curr_y_offset,next_y_offset,em_y_offset);
	printf("dimy = %d, dimz = %d\n",dimy,dimz);
	printf("dt = %e, a1h = %e, a2h = %e, a3h = %e, a4h = %e, a5h = %e, a1z = %e, a2z = %e, a3z = %e, a4z = %e, a5z = %e, e1z = %e\n",dt,a1h,a2h,a3h,a4h,a5h,a1z,a2z,a3z,a4z,a5z,e1z);
	printf("inv_Q_min = %e, inv_Q_scaler = %e\n",inv_Q_min,inv_Q_scaler);
	printf("Den_min = %e, Den_scaler = %e\n",Den_min,Den_scaler);
	printf("Dip_min = %e, Dip_scaler = %e\n",Dip_min,Dip_scaler);
	printf("Azm_min = %e, Azm_scaler = %e\n",Azm_min,Azm_scaler);
	printf("C44C33_min = %e, C44C33_scaler = %e\n",C44C33_min,C44C33_scaler);
	printf("Vel_min = %e, Vel_scaler = %e\n",Vel_min,Vel_scaler);
	printf("Del_min = %e, Del_scaler = %e\n",Del_min,Del_scaler);
	printf("Eps_min = %e, Eps_scaler = %e\n",Eps_min,Eps_scaler);
	printf("stripe_xsrc = %d, stripe_ysrc = %d, stripe_zsrc = %d, stripe_zsrcghost = %d, src_val = %e\n",stripe_xsrc,stripe_ysrc,stripe_zsrc,stripe_zsrcghost);
#endif
	const int stride_z = 32;
	const int stride_y = stride_z * dimz;

	const int stride_z_em = 16;
	const int stride_wf_em = stride_z_em * dimz;
	const int stride_y_em = 2 * stride_wf_em;

	int em_idx = em_y_offset * stride_y_em;
	int curr_idx = curr_y_offset * stride_y;
        int next_idx = next_y_offset * stride_y;  // NB! This index points to first output sample and is thus equal to next_yoff9_idx.

	int em_yoff9_idx = (em_y_offset + 9) * stride_y_em;
	int curr_yoff9_idx = (curr_y_offset + 9) * stride_y;
       	int prev_yoff9_idx = (prev_y_offset + 9) * stride_y;

	Host_Inject_Source(
		compute_stream,
		src_val_1,stripe_xsrc_1,stripe_ysrc_1,stripe_zsrc_1,stripe_zsrcghost_1,(float*)d_curr_pq1,
		src_val_2,stripe_xsrc_2,stripe_ysrc_2,stripe_zsrc_2,stripe_zsrcghost_2,(float*)d_curr_pq2,
		stride_y,stride_z,curr_yoff9_idx
		);

	dim3 blockShape1(32,6,1);
	dim3 gridShape1(1,dimz/4);

	dim3 blockShape2(32,6,1);
        dim3 gridShape2(1,dimy-9);

	Compute_DX_DY_Main<<<gridShape1,blockShape1,0,compute_stream>>>(
			d_curr_pq0+curr_idx,		d_curr_pq1+curr_idx,		d_curr_pq2+curr_idx,		d_curr_pq3+curr_idx,
			d_dx,				d_dy,
			dimy,dimz,
			a1h,a2h,a3h,a4h,a5h);

	Compute_DZ_Main<<<gridShape2,blockShape2,0,compute_stream>>>(
			d_curr_pq0+curr_idx,		d_curr_pq1+curr_idx,		d_curr_pq2+curr_idx,		d_curr_pq3+curr_idx,
			d_dz,
			dimy,dimz,
			a1z,a2z,a3z,a4z,a5z,e1z);

	int Is_TTI = Kernel_Type == 2 ? 1 : 0;
	
	dim3 blockShape3(32,4,1);
	dim3 gridShape3(1,dimy-9);
	if (Is_TTI)
	{
		Compute_DYED1and2_Part_V4_V5<<<gridShape3,blockShape3,0,compute_stream>>>(
				d_dx,				d_dy,				d_dz,
				d_em0+em_idx,			d_em1+em_idx,			d_em2+em_idx,			d_em3+em_idx,
				dimy,dimz,
				dt,
				a1z,a2z,a3z,a4z,a5z,e1z,
				a1h,a2h,a3h,a4h,a5h,
				inv_Q_min,inv_Q_scaler,Den_min,Den_scaler,Dip_min,Dip_scaler,Azm_min,Azm_scaler,C44C33_min,C44C33_scaler,Vel_min,Vel_scaler,Del_min,Del_scaler,Eps_min,Eps_scaler);
	}
	else
	{
		Compute_VTI_DYED1and2_Part_V4_V5<<<gridShape3,blockShape3,0,compute_stream>>>(
				d_dx,				d_dy,				d_dz,
				d_em0+em_idx,			d_em1+em_idx,			d_em2+em_idx,			d_em3+em_idx,
				dimy,dimz,
				dt,
				a1z,a2z,a3z,a4z,a5z,e1z,
				a1h,a2h,a3h,a4h,a5h,
				inv_Q_min,inv_Q_scaler,Den_min,Den_scaler,C44C33_min,C44C33_scaler,Vel_min,Vel_scaler,Del_min,Del_scaler,Eps_min,Eps_scaler);
	}

	dim3 blockShape4(32,4,1);
	dim3 gridShape4(1,dimz);
	if (Is_TTI)
	{
		Compute_Next_PQ_T2<<<gridShape4,blockShape4,0,compute_stream>>>(
				d_spg_x + x0_1,			d_spg_y + y0,			d_spg_z,
				d_dx, 				d_dy,
				d_curr_pq1+curr_yoff9_idx, 	d_curr_pq2+curr_yoff9_idx,
				d_prev_pq1+prev_yoff9_idx, 	d_prev_pq2+prev_yoff9_idx,
				d_next_pq1+next_idx, 		d_next_pq2+next_idx,
				d_em1+em_yoff9_idx, 		d_em2+em_yoff9_idx,
				dimy-18,dimz,
				dt,
				a1h,a2h,a3h,a4h,a5h,
				inv_Q_min,inv_Q_scaler,Den_min,Den_scaler,Dip_min,Dip_scaler,Azm_min,Azm_scaler,C44C33_min,C44C33_scaler,Vel_min,Vel_scaler,Del_min,Del_scaler,Eps_min,Eps_scaler
				);
	}
	else
	{
		Compute_VTI_Next_PQ_T2<<<gridShape4,blockShape4,0,compute_stream>>>(
				d_spg_x + x0_1,			d_spg_y + y0,			d_spg_z,
				d_dx, 				d_dy,
				d_curr_pq1+curr_yoff9_idx, 	d_curr_pq2+curr_yoff9_idx,
				d_prev_pq1+prev_yoff9_idx, 	d_prev_pq2+prev_yoff9_idx,
				d_next_pq1+next_idx, 		d_next_pq2+next_idx,
				d_em1+em_yoff9_idx, 		d_em2+em_yoff9_idx,
				dimy-18,dimz,
				dt,
				a1h,a2h,a3h,a4h,a5h,
				inv_Q_min,inv_Q_scaler,Den_min,Den_scaler,C44C33_min,C44C33_scaler,Vel_min,Vel_scaler,Del_min,Del_scaler,Eps_min,Eps_scaler
				);
	}

	if (recnum_1 > 0 || recnum_2 > 0)
	{
		// temporarily remove source from previous timestep
		Host_Inject_Source(
				compute_stream,
				-src_val_1,stripe_xsrc_1,stripe_ysrc_1,stripe_zsrc_1,stripe_zsrcghost_1,(float*)d_curr_pq1,
				-src_val_2,stripe_xsrc_2,stripe_ysrc_2,stripe_zsrc_2,stripe_zsrcghost_2,(float*)d_curr_pq2,
				stride_y,stride_z,curr_yoff9_idx
				);

		if (recnum_1 > 0)
		{
			// output receiver values for first stripe
			//printf("GPU kernel to extract receiver values x=%d, y=%d\n",x0_1,y0);
			dim3 blockShape5(32,4,1);
			dim3 gridShape5((recnum_1+127)/128,1);
			Extract_Receiver_Values<<<gridShape5,blockShape5,0,compute_stream>>>(
					tfrac1,
					d_rcx_1, 			d_rcy_1, 			d_rcz_1, 			d_out_1, 		recnum_1,
					d_curr_pq1 + curr_yoff9_idx,	d_next_pq1 + next_idx,
					d_spg_x + x0_1,			d_spg_y + y0,			d_spg_z,
					x0_1,				y0,
					16,				dimy,				dimz,
					recghost_Flag,			z_freesurface
					);
		}
		if (recnum_2 > 0)
		{
			// output receiver values for second stripe
			//printf("GPU kernel to extract receiver values x=%d, y=%d\n",x0_2,y0);
			dim3 blockShape6(32,4,1);
			dim3 gridShape6((recnum_2+127)/128,1);
			Extract_Receiver_Values<<<gridShape6,blockShape6,0,compute_stream>>>(
					tfrac2,
					d_rcx_2, 			d_rcy_2, 			d_rcz_2, 			d_out_2, 		recnum_2,
					d_curr_pq2 + curr_yoff9_idx,	d_next_pq2 + next_idx,
					d_spg_x + x0_2,			d_spg_y + y0,			d_spg_z,
					x0_2,				y0,
					16,				dimy,				dimz,
					recghost_Flag,			z_freesurface
					);
		}

		// add back source term
		Host_Inject_Source(
				compute_stream,
				src_val_1,stripe_xsrc_1,stripe_ysrc_1,stripe_zsrc_1,stripe_zsrcghost_1,(float*)d_curr_pq1,
				src_val_2,stripe_xsrc_2,stripe_ysrc_2,stripe_zsrc_2,stripe_zsrcghost_2,(float*)d_curr_pq2,
				stride_y,stride_z,curr_yoff9_idx
				);
	}
}

