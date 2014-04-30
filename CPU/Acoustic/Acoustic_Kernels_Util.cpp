#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <omp.h>
#include <time.h>

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

#define OSXSAVEFlag (1UL<<27)
#define AVXFlag     ((1UL<<28)|OSXSAVEFlag)
#define FMAFlag     ((1UL<<12)|AVXFlag|OSXSAVEFlag)
#define CLMULFlag   ((1UL<< 1)|AVXFlag|OSXSAVEFlag)
#define VAESFlag    ((1UL<<25)|AVXFlag|OSXSAVEFlag)

#define U32 unsigned int
#define __cpuid(level, a, b, c, d)                      \
__asm__ ("cpuid\n\t"                                        \
        : "=a" (a), "=b" (b), "=c" (c), "=d" (d)     \
        : "0" (level))

/*
inline void cpuid(U32 level, int& a, int& b, int& c, int& d)
{
        __asm__ ("cpuid\n\t"
                : "=a" (a), "=b" (b), "=c" (c), "=d" (d)
                : "0" (level));
}
*/

inline bool SimdDetectFeature(U32 idFeature)
{
    int EAX, EBX, ECX, EDX;
    __cpuid(1, EAX, EBX, ECX, EDX);
    if((ECX & idFeature) != idFeature)
        return false;
    return true;
}

inline bool SupportsAVX()
{
        return SimdDetectFeature(AVXFlag);
}

// dip and azm are in degrees.
void Compute_Earth_Model(int iso, float dip, float azm)
{
	lutVp2  = new float[1];
	lutEps = new float[1];
	lutDel = new float[1];
	lutDen = new float[1];
	lutBuoy = new float[1];
	lutsDip = new float[1];
	lutcDip = new float[1];
	lutsAzm = new float[1];
	lutcAzm = new float[1];
	lutc44c33 = new float[2];

	float deg2rad = 0.01745329251994329509f;

	lutVp2[0] = 3.439728e+01f;
	lutDen[0] = 2.0f;
	lutBuoy[0] = 5e-1f;
	lutc44c33[1] = 0.0f;

	if (iso)
	{
		lutEps[0] = 3.4e-1f;
		lutDel[0] = 1e-1f;
		//lutEps[0] = 0.0f;
		//lutDel[0] = 0.0f;
		lutsDip[0] = sin(dip*deg2rad);
		lutcDip[0] = cos(dip*deg2rad);
		lutsAzm[0] = sin(azm*deg2rad);
		lutcAzm[0] = cos(azm*deg2rad);
		//lutsDip[0] = 0.0f;
		//lutcDip[0] = 1.0f;
		//lutsAzm[0] = 0.0f;
		//lutcAzm[0] = 1.0f;
		lutc44c33[0] = 0.0f;
	}
	else
	{
		lutEps[0] = 3.4e-1f;
		lutDel[0] = 1e-1f;
		lutsDip[0] = sin(dip*deg2rad);
		lutcDip[0] = cos(dip*deg2rad);
		lutsAzm[0] = sin(azm*deg2rad);
		lutcAzm[0] = cos(azm*deg2rad);
		lutc44c33[0] = 1e-2f;
	}
}

void Compute_Stencils(
	int OTflag,
	float dh,
	float dz
	)
{
	//
	// 10pt stencils on staggered grid
	//

        float a1, a2, a3, a4, a5;
        float b1, b2, b3, b4;
        float c1, c2, c3;
        float d1, d2;
        float e1;

        if(OTflag==4)
        {
                a1 =  1.250394163714f;  a2 = -0.119656543874f;  a3 = 0.031206223579f;
                a4 = -0.009128136972f;  a5 =  0.001882183398f;
        }
        else
        {
                a1 =  1.248489029341f;  a2 = -0.120133754290f;  a3 = 0.031688119039f;
                a4 = -0.008048796917f;  a5 =  0.001090357653f;
        }
        b1 =  1.231650129521f;  b2 = -0.103861125624f;  b3 = 0.020166542235f;
        b4 = -0.002985637689f;
        c1 =  1.199634495725f;  c2 = -0.080370339530f;  c3 = 0.008295304573f;
        d1 =  1.134389630713f;  d2 = -0.044796543571f;
        e1 =  0.5f;

        A1h = a1/dh;  A2h = a2/dh;  A3h = a3/dh;  A4h = a4/dh;  A5h = a5/dh;
        B1h = b1/dh;  B2h = b2/dh;  B3h = b3/dh;  B4h = b4/dh;
        C1h = c1/dh;  C2h = c2/dh;  C3h = c3/dh;
        D1h = d1/dh;  D2h = d2/dh;
        E1h = e1/dh;

	/*
	printf("A1h=%e\n",A1h);
	printf("A2h=%e\n",A2h);
	printf("A3h=%e\n",A3h);
	printf("A4h=%e\n",A4h);
	printf("A5h=%e\n",A5h);
	*/

        A1 = a1/dz;  A2 = a2/dz;  A3 = a3/dz;  A4 = a4/dz;  A5 = a5/dz;
        B1 = b1/dz;  B2 = b2/dz;  B3 = b3/dz;  B4 = b4/dz;
        C1 = c1/dz;  C2 = c2/dz;  C3 = c3/dz;
        D1 = d1/dz;  D2 = d2/dz;
        E1 = e1/dz;

	/*
	printf("A1=%e\n",A1);
	printf("A2=%e\n",A2);
	printf("A3=%e\n",A3);
	printf("A4=%e\n",A4);
	printf("A5=%e\n",A5);
	*/

	A5h_A4h  = _mm_setr_ps( A5h,  A5h,  A4h,  A4h);
	A4h_A3h  = _mm_setr_ps( A4h,  A4h,  A3h,  A3h);
	A3h_A2h  = _mm_setr_ps( A3h,  A3h,  A2h,  A2h);
	A2h_A1h  = _mm_setr_ps( A2h,  A2h,  A1h,  A1h);
	A5h_mA5h = _mm_setr_ps( A5h,  A5h, -A5h, -A5h);
	mA1h_A1h = _mm_setr_ps(-A1h, -A1h,  A1h,  A1h);

	A5h_A4h_A3h_A2h = _mm_setr_ps(A5h, A4h, A3h, A2h);
	A1h_A1h_A1h_A1h = _mm_set1_ps(A1h);

        A5_A4_A3_A2 = _mm_setr_ps(A5, A4, A3, A2);
	A1_A1_A1_A1 = _mm_set1_ps(A1);

	//
	// 11pt central difference stencils
	// 

	A1h_1st_11pt =  0.833333333333f / dh;
	A2h_1st_11pt = -0.238095238100f / dh;
	A3h_1st_11pt =  5.952380952e-2f / dh;
	A4h_1st_11pt = -9.920634921e-3f / dh;
	A5h_1st_11pt =  7.936507937e-4f / dh;

	zzz_A5h_A4h_A3h_1st_11pt = _mm_setr_ps(0.0f, A5h_1st_11pt, A4h_1st_11pt, A3h_1st_11pt);
	A2h_A1h_zzz_zzz_1st_11pt = _mm_setr_ps(A2h_1st_11pt, A1h_1st_11pt, 0.0f, 0.0f);

	A0h_2nd_11pt = -3.1487087031f / (dh*dh);
	A1h_2nd_11pt =  1.8617697535f / (dh*dh);
	A2h_2nd_11pt = -0.3705141600f / (dh*dh);
	A3h_2nd_11pt =  0.1068406382f / (dh*dh);
	A4h_2nd_11pt = -0.0281145606f / (dh*dh);
	A5h_2nd_11pt =  0.0043726804f / (dh*dh);

	zzz_A5h_A4h_A3h_2nd_11pt = _mm_setr_ps(0.0f, A5h_2nd_11pt, A4h_2nd_11pt, A3h_2nd_11pt);
	A2h_A1h_A0h_zzz_2nd_11pt = _mm_setr_ps(A2h_2nd_11pt, A1h_2nd_11pt, A0h_2nd_11pt, 0.0f);

	A1_1st_11pt =  0.833333333333f / dz;
	A2_1st_11pt = -0.238095238100f / dz;
	A3_1st_11pt =  5.952380952e-2f / dz;
	A4_1st_11pt = -9.920634921e-3f / dz;
	A5_1st_11pt =  7.936507937e-4f / dz;

	zz_A5_A4_A3_1st_11pt = _mm_setr_ps(0.0f, A5_1st_11pt, A4_1st_11pt, A3_1st_11pt);
	A2_A1_zz_zz_1st_11pt = _mm_setr_ps(A2_1st_11pt, A1_1st_11pt, 0.0f, 0.0f);

	A0_2nd_11pt = -3.1487087031f / (dz*dz);
	A1_2nd_11pt =  1.8617697535f / (dz*dz);
	A2_2nd_11pt = -0.3705141600f / (dz*dz);
	A3_2nd_11pt =  0.1068406382f / (dz*dz);
	A4_2nd_11pt = -0.0281145606f / (dz*dz);
	A5_2nd_11pt =  0.0043726804f / (dz*dz);

	zz_A5_A4_A3_2nd_11pt = _mm_setr_ps(0.0f, A5_2nd_11pt, A4_2nd_11pt, A3_2nd_11pt);
	A2_A1_A0_zz_2nd_11pt = _mm_setr_ps(A2_2nd_11pt, A1_2nd_11pt, A0_2nd_11pt, 0.0f);
}

/*
** num_x_zeropad	If > 0, right side of x sponge is shifted this many cells to the left. The remanining cells on the right side are zeroed.
**			This will happen if actual dimx is not a multiple of 4, in which case num_x_zeropad equals dimx - actual_dimx.
*/
void Compute_Sponges(
	float spongecoeff_x,
	float spongecoeff_y,
	float spongecoeff_z_lo,
	float spongecoeff_z_hi,
	int spongewidth_x,
	int spongewidth_y,
	int spongewidth_z_lo,
	int spongewidth_z_hi,
	int absorbz0,
	int dimx,
	int num_x_zeropad,
	int dimy,
	int dimz,
	float* spgx,
	float* spgy,
	float* spgz
	)
{
	int xstart = dimx - spongewidth_x;
	spgx[0] = spgx[dimx-1-num_x_zeropad] = 0.0f;
	for (int i = 1;  i < spongewidth_x;  ++i)
	{
		spgx[i] = 1.0 - spongecoeff_x * (spongewidth_x - i) * (spongewidth_x - i);
		//printf("spgx[%d] = %f\n",i,spgx[i]);
	}
	for (int i = spongewidth_x;  i < xstart;  ++i)
	{
		spgx[i] = 1.0f;
		//printf("spgx[%d] = %f\n",i,spgx[i]);
	}
	for (int i = xstart;  i < dimx-1;  ++i)
	{
		spgx[i-num_x_zeropad] = 1.0f - spongecoeff_x * (i - xstart + 1) * (i - xstart + 1);
		//printf("spgx[%d] = %f\n",i,spgx[i]);
	}
	if (num_x_zeropad > 0)
	{
		for (int i = dimx-1-num_x_zeropad;  i < dimx;  ++i)
		{
			spgx[i] = 0.0f;
		}
	}

	int ystart = dimy - spongewidth_y;
	spgy[0] = spgy[dimy-1] = 0.0f;
	for (int i = 1;  i < spongewidth_y;  ++i)
	{
		spgy[i] = 1.0 - spongecoeff_y * (spongewidth_y - i) * (spongewidth_y - i);
		//printf("spgy[%d] = %f\n",i,spgy[i]);
	}
	for (int i = spongewidth_y;  i < ystart;  ++i)
	{
		spgy[i] = 1.0f;
		//printf("spgy[%d] = %f\n",i,spgy[i]);
	}
	for (int i = ystart;  i < dimy-1;  ++i)
	{
		spgy[i] = 1.0f - spongecoeff_y * (i - ystart + 1) * (i - ystart + 1);
		//printf("spgy[%d] = %f\n",i,spgy[i]);
	}

	int zstart = dimz - spongewidth_z_hi;
	spgz[0] = absorbz0 ? 0.0f : 1.0f;
	spgz[dimz-1] = 0.0f;
	for (int i = 1;  i < spongewidth_z_lo;  ++i)
	{
		spgz[i] = absorbz0 ? 1.0 - spongecoeff_z_lo * (spongewidth_z_lo - i) * (spongewidth_z_lo - i) : 1.0f;
		//printf("spgz[%d] = %f\n",i,spgz[i]);
	}
	for (int i = spongewidth_z_lo;  i < zstart;  ++i)
	{
		spgz[i] = 1.0f;
		//printf("spgz[%d] = %f\n",i,spgz[i]);
	}
	for (int i = zstart;  i < dimz-1;  ++i)
	{
		spgz[i] = 1.0f - spongecoeff_z_hi * (i - zstart + 1) * (i - zstart + 1);
		//printf("spgz[%d] = %f\n",i,spgz[i]);
	}
}

void Precompute_Anisotropy_Parameters()
{
	for (int i = 0;  i <= EPSMASK;  ++i)
	{
		lutEps[i] = 1.0f + 2.0f * lutEps[i];
	}
	for (int i = 0;  i <= DELMASK;  ++i)
	{
		lutDel[i] = 2.0f * lutDel[i];
	}
}

void Set_P(
        __m128* pq,
        int dimx,
        int dimy,
        int xh,
        int yh,
	int zh,
        int iX,
        int iY,
        int iZ,
        float val
        )
{
        float* f = (float*)pq;
        int dimxh = dimx + 2*xh;
        int dimyh = dimy + 2*yh;
	unsigned long idx = ((unsigned long)(iZ+zh) * (unsigned long)dimyh * (unsigned long)dimxh + (unsigned long)(iY+yh) * (unsigned long)dimxh + (unsigned long)(iX+xh)) << 1;
        f[idx] = val;
}

void Set_Q(
        __m128* pq,
        int dimx,
        int dimy,
        int xh,
        int yh,
	int zh,
        int iX,
        int iY,
        int iZ,
        float val
        )
{
        float* f = (float*)pq;
        int dimxh = dimx + 2*xh;
        int dimyh = dimy + 2*yh;
	unsigned long idx = ((unsigned long)(iZ+zh) * (unsigned long)dimyh * (unsigned long)dimxh + (unsigned long)(iY+yh) * (unsigned long)dimxh + (unsigned long)(iX+xh)) << 1;
        f[idx+1] = val;
}

float Get_P(
	int KernelType,
        __m128* pq,
        int dimx,
        int dimy,
        int xh,
        int yh,
	int zh,
        int iX,
        int iY,
        int iZ
        )
{
        float* f = (float*)pq;
        int dimxh = dimx + 2*xh;
        int dimyh = dimy + 2*yh;
	unsigned long idx = ((unsigned long)(iZ+zh) * (unsigned long)dimyh * (unsigned long)dimxh + (unsigned long)(iY+yh) * (unsigned long)dimxh + (unsigned long)(iX+xh));
	if (KernelType > 0) idx = idx << 1;
        return f[idx];
}

float Get_Q(
        __m128* pq,
        int dimx,
        int dimy,
        int xh,
        int yh,
	int zh,
        int iX,
        int iY,
        int iZ
        )
{
        float* f = (float*)pq;
        int dimxh = dimx + 2*xh;
        int dimyh = dimy + 2*yh;
	unsigned long idx = ((unsigned long)(iZ+zh) * (unsigned long)dimyh * (unsigned long)dimxh + (unsigned long)(iY+yh) * (unsigned long)dimxh + (unsigned long)(iX+xh)) << 1;
        return f[idx+1];
}

char* Clean_Filename(int logLevel, char* filename)
{
	char* fname = new char[1024];
	if (filename[0] == '~')
	{
		char* HOME = getenv("HOME");
		sprintf(fname, "%s%s", HOME, filename+1);
		if (logLevel >= 6)
		{
			printf("Clean_Filename :: %s => %s\n", filename, fname);
			fflush(stdout);
		}
	}
	else
	{
		sprintf(fname, "%s", filename);
	}
	return fname;
}

FILE* Clean_FOPEN(int logLevel, char* filename, char* mode)
{
	char* fname = Clean_Filename(logLevel,filename);
	FILE* fp = fopen(fname, mode);
	delete [] fname;
	return fp;
}

void Dump_P(
	int KernelType,
        char* wfName,
        __m128* pq,
        int dimx,
        int dimy,
        int dimz,
        int xh,
        int yh,
        int zh
        )
{
        float* f = (float*)pq;
	printf("\nSTART - P -  %s :: %s :: %s - P - START\n",wfName,wfName,wfName);
        for (int iZ = 0;  iZ < dimz;  ++iZ)
        {
                for (int iY = 0;  iY < dimy;  ++iY)
                {
                        for (int iX = 0;  iX < dimx;  ++iX)
                        {
                                float val = Get_P(KernelType,pq,dimx,dimy,xh,yh,zh,iX,iY,iZ);
                                if ( val != 0.0f)
                                {
                                        printf("%s[%3d,%3d,%3d] = %e\n",wfName,iX,iY,iZ,val);
                                        //if (iZ == 286)
                                        //{
                                        //      printf("%d,%d,%d,%e\n",iX,iY,iZ,val);
                                        //}
                                }
                        }
                }
        }
	printf("END - P - %s :: %s :: %s - P - END\n\n",wfName,wfName,wfName);
}

void Dump_Q(
        char* wfName,
        __m128* pq,
        int dimx,
        int dimy,
        int dimz,
        int xh,
        int yh,
        int zh
        )
{
        float* f = (float*)pq;
	printf("\nSTART - Q -  %s :: %s :: %s - Q - START\n",wfName,wfName,wfName);
        for (int iZ = 0;  iZ < dimz;  ++iZ)
        {
                for (int iY = 0;  iY < dimy;  ++iY)
                {
                        for (int iX = 0;  iX < dimx;  ++iX)
                        {
                                float val = Get_Q(pq,dimx,dimy,xh,yh,zh,iX,iY,iZ);
                                if ( val != 0.0f)
                                {
                                        printf("%s[%3d,%3d,%3d] = %e\n",wfName,iX,iY,iZ,val);
                                        //if (iZ == 286)
                                        //{
                                        //      printf("%d,%d,%d,%e\n",iX,iY,iZ,val);
                                        //}
                                }
                        }
                }
        }
	printf("END - Q -  %s :: %s :: %s - Q - END\n\n",wfName,wfName,wfName);
}

void Write_XY_Slice_P(
	int logLevel,
	int KernelType,
        char* filename,
        __m128* pq,
        int dimx,
        int dimy,
        int dimz,
        int xh,
        int yh,
        int zh,
        int iZ
        )
{
        FILE* fp = Clean_FOPEN(logLevel,filename,"w");
        if (fp != 0L)
        {
                float* f = (float*)pq;
                for (int iY = 0;  iY < dimy;  ++iY)
                {
                        for (int iX = 0;  iX < dimx;  ++iX)
                        {
                                float val = Get_P(KernelType,pq,dimx,dimy,xh,yh,zh,iX,iY,iZ);
                                fprintf(fp, "%d %d %e\n", iX, iY, val);
                        }
                        fprintf(fp,"\n");
                }
                fclose(fp);
        }
}

void Write_XZ_Slice_P(
	int logLevel,
	int KernelType,
        char* filename,
        __m128* pq,
        int dimx,
        int dimy,
        int dimz,
        int xh,
        int yh,
        int zh,
        int iY
        )
{
        FILE* fp = Clean_FOPEN(logLevel,filename,"w");
        if (fp != 0L)
        {
                float* f = (float*)pq;
                for (int iZ = 0;  iZ < dimz;  ++iZ)
                {
                        for (int iX = 0;  iX < dimx;  ++iX)
                        {
                                float val = Get_P(KernelType,pq,dimx,dimy,xh,yh,zh,iX,iY,iZ);
                                fprintf(fp, "%d %d %e\n", iX, iZ, val);
                        }
                        fprintf(fp,"\n");
                }
                fclose(fp);
        }
}

void Write_XZ_Slice_EM(
        int logLevel,
        char* filename,
        int* PadDenAng,
        int* PadVelAnis,
        int KernelType,
        int dimx,
        int dimy,
        int dimz,
        int xh,
        int yh,
        int zh,
        int iY
        )
{
        char str[256];

        sprintf(str, "%s_Vp2_y=%d", filename, iY);
        FILE* fp_Vp2 = Clean_FOPEN(logLevel,str, "w");

        sprintf(str, "%s_Den_y=%d", filename, iY);
        FILE* fp_Den = Clean_FOPEN(logLevel,str, "w");

        sprintf(str, "%s_Buoy_y=%d", filename, iY);
        FILE* fp_Buoy = Clean_FOPEN(logLevel,str, "w");

        sprintf(str, "%s_Qatten_y=%d", filename, iY);
        FILE* fp_Qatten = Clean_FOPEN(logLevel,str, "w");

        FILE* fp_Eps = 0L;
        FILE* fp_Del = 0L;
        FILE* fp_C44C33 = 0L;
        if (KernelType > 0)
        {
                sprintf(str, "%s_Eps_y=%d", filename, iY);
                fp_Eps = Clean_FOPEN(logLevel,str, "w");

                sprintf(str, "%s_Del_y=%d", filename, iY);
                fp_Del = Clean_FOPEN(logLevel,str, "w");

                sprintf(str, "%s_C44C33_y=%d", filename, iY);
                fp_C44C33 = Clean_FOPEN(logLevel,str, "w");
        }

        FILE* fp_cAzm = 0L;
        FILE* fp_sAzm = 0L;
        FILE* fp_cDip = 0L; 
        FILE* fp_sDip = 0L;
        if (KernelType > 1)
        {
                sprintf(str, "%s_cAzm_y=%d", filename, iY);
                fp_cAzm = Clean_FOPEN(logLevel,str, "w");

                sprintf(str, "%s_sAzm_y=%d", filename, iY);
                fp_sAzm = Clean_FOPEN(logLevel,str, "w");

                sprintf(str, "%s_cDip_y=%d", filename, iY);
                fp_cDip = Clean_FOPEN(logLevel,str, "w");

                sprintf(str, "%s_sDip_y=%d", filename, iY);
                fp_sDip = Clean_FOPEN(logLevel,str, "w");
        }

        int stride_y_em = dimx + 2*xh;
        int stride_z_em = (dimy + 2*yh) * stride_y_em;

        for (int iZ = 0;  iZ < dimz;  ++iZ)
        {
                for (int iX = 0;  iX < dimx;  iX+=4)
                {
                        unsigned long EM_Offset = (unsigned long)(iZ+zh) * (unsigned long)stride_z_em + (unsigned long)(iY+yh) * (unsigned long)stride_y_em + (unsigned long)(iX+xh);
                        int *cDenAng = PadDenAng + EM_Offset;
                        int *cVelAnis = PadVelAnis + EM_Offset;

                        __m128 Buoy, cDip, sDip, cAzm, sAzm;
                        if (KernelType > 1)
                        {
                                TTIDenQ_Get_EM1(cDenAng,Buoy,cDip,sDip,cAzm,sAzm);
                        }
                        else if (KernelType > 0)
                        {
                                VTIDenQ_Get_EM1(cDenAng,Buoy);
                        }
			else
			{
				ISODenQ_Get_EM1(cVelAnis,Buoy);
			}

                        __m128 Den, Qatten, Vp2, C44C33, Del, Eps;
			if (KernelType > 0)
			{
				TTIDenQ_Get_EM2_Raw(cVelAnis,cDenAng, Den, Qatten, Vp2, C44C33, Del, Eps);
			}
			else
			{
				ISODenQ_Get_EM2_Raw(cVelAnis, Den, Qatten, Vp2);
			}

                        float* p = (float*)&Vp2;
                        fprintf(fp_Vp2, "%d %d %e\n",iX  ,iZ,p[0]);
                        fprintf(fp_Vp2, "%d %d %e\n",iX+1,iZ,p[1]);
                        fprintf(fp_Vp2, "%d %d %e\n",iX+2,iZ,p[2]);
                        fprintf(fp_Vp2, "%d %d %e\n",iX+3,iZ,p[3]);

                        p = (float*)&Den;
                        fprintf(fp_Den, "%d %d %e\n",iX  ,iZ,p[0]);
                        fprintf(fp_Den, "%d %d %e\n",iX+1,iZ,p[1]);
                        fprintf(fp_Den, "%d %d %e\n",iX+2,iZ,p[2]);
                        fprintf(fp_Den, "%d %d %e\n",iX+3,iZ,p[3]);

                        p = (float*)&Qatten;
                        fprintf(fp_Qatten, "%d %d %e\n",iX  ,iZ,p[0]);
                        fprintf(fp_Qatten, "%d %d %e\n",iX+1,iZ,p[1]);
                        fprintf(fp_Qatten, "%d %d %e\n",iX+2,iZ,p[2]);
                        fprintf(fp_Qatten, "%d %d %e\n",iX+3,iZ,p[3]);

                        p = (float*)&Buoy;
                        fprintf(fp_Buoy, "%d %d %e\n",iX  ,iZ,p[0]);
                        fprintf(fp_Buoy, "%d %d %e\n",iX+1,iZ,p[1]);
                        fprintf(fp_Buoy, "%d %d %e\n",iX+2,iZ,p[2]);
                        fprintf(fp_Buoy, "%d %d %e\n",iX+3,iZ,p[3]);

                        if (KernelType > 0)
                        {
                                p = (float*)&Eps;
                                fprintf(fp_Eps, "%d %d %e\n",iX  ,iZ,p[0]);
                                fprintf(fp_Eps, "%d %d %e\n",iX+1,iZ,p[1]);
                                fprintf(fp_Eps, "%d %d %e\n",iX+2,iZ,p[2]);
                                fprintf(fp_Eps, "%d %d %e\n",iX+3,iZ,p[3]);

                                p = (float*)&Del;
                                fprintf(fp_Del, "%d %d %e\n",iX  ,iZ,p[0]);
                                fprintf(fp_Del, "%d %d %e\n",iX+1,iZ,p[1]);
                                fprintf(fp_Del, "%d %d %e\n",iX+2,iZ,p[2]);
                                fprintf(fp_Del, "%d %d %e\n",iX+3,iZ,p[3]);

                                p = (float*)&C44C33;
                                fprintf(fp_C44C33, "%d %d %e\n",iX  ,iZ,p[0]);
                                fprintf(fp_C44C33, "%d %d %e\n",iX+1,iZ,p[1]);
                                fprintf(fp_C44C33, "%d %d %e\n",iX+2,iZ,p[2]);
                                fprintf(fp_C44C33, "%d %d %e\n",iX+3,iZ,p[3]);
                        }

                        if (KernelType > 1)
                        {
                                p = (float*)&cAzm;
                                fprintf(fp_cAzm, "%d %d %e\n",iX  ,iZ,p[0]);
                                fprintf(fp_cAzm, "%d %d %e\n",iX+1,iZ,p[1]);
                                fprintf(fp_cAzm, "%d %d %e\n",iX+2,iZ,p[2]);
                                fprintf(fp_cAzm, "%d %d %e\n",iX+3,iZ,p[3]);

                                p = (float*)&sAzm;
                                fprintf(fp_sAzm, "%d %d %e\n",iX  ,iZ,p[0]);
                                fprintf(fp_sAzm, "%d %d %e\n",iX+1,iZ,p[1]);
                                fprintf(fp_sAzm, "%d %d %e\n",iX+2,iZ,p[2]);
                                fprintf(fp_sAzm, "%d %d %e\n",iX+3,iZ,p[3]);

                                p = (float*)&cDip;
                                fprintf(fp_cDip, "%d %d %e\n",iX  ,iZ,p[0]);
                                fprintf(fp_cDip, "%d %d %e\n",iX+1,iZ,p[1]);
                                fprintf(fp_cDip, "%d %d %e\n",iX+2,iZ,p[2]);
                                fprintf(fp_cDip, "%d %d %e\n",iX+3,iZ,p[3]);

                                p = (float*)&sDip;
                                fprintf(fp_sDip, "%d %d %e\n",iX  ,iZ,p[0]);
                                fprintf(fp_sDip, "%d %d %e\n",iX+1,iZ,p[1]);
                                fprintf(fp_sDip, "%d %d %e\n",iX+2,iZ,p[2]);
                                fprintf(fp_sDip, "%d %d %e\n",iX+3,iZ,p[3]);
                        }
                }
                fprintf(fp_Vp2, "\n");
                fprintf(fp_Den, "\n");
                fprintf(fp_Qatten, "\n");
                fprintf(fp_Buoy, "\n");
                if (KernelType > 0)
                {
                        fprintf(fp_Eps, "\n");
                        fprintf(fp_Del, "\n");
                        fprintf(fp_C44C33, "\n");
                }
                if (KernelType > 1)
                {
                        fprintf(fp_cAzm, "\n");
                        fprintf(fp_sAzm, "\n");
                        fprintf(fp_cDip, "\n");
                        fprintf(fp_sDip, "\n");
                }
        }

        if (fp_Vp2 != 0L) fclose(fp_Vp2);
        if (fp_Eps != 0L) fclose(fp_Eps);
        if (fp_Del != 0L) fclose(fp_Del);
        if (fp_C44C33 != 0L) fclose(fp_C44C33);
        if (fp_Den != 0L) fclose(fp_Den);
        if (fp_Qatten != 0L) fclose(fp_Qatten);
        if (fp_Buoy != 0L) fclose(fp_Buoy);
        if (fp_cAzm != 0L) fclose(fp_cAzm);
        if (fp_sAzm != 0L) fclose(fp_sAzm);
        if (fp_cDip != 0L) fclose(fp_cDip);
        if (fp_sDip != 0L) fclose(fp_sDip);
}

void Range_P(
	int KernelType,
        __m128* pq,
        int dimx,
        int dimy,
        int dimz,
        int xh,
        int yh,
        int zh
        )
{
        float* f = (float*)pq;
        int min_X = dimx;
        int max_X = -1;
        int min_Y = dimy;
        int max_Y = -1;
        int min_Z = dimz;
        int max_Z = -1;
        float min_abs_val = 1e20f;
        float max_abs_val = 0.0f;
        int dimxh = dimx + 2*xh;
        int dimyh = dimy + 2*yh;
        int dimzh = dimz + 2*zh;
        for (int iZ = 0;  iZ < dimz;  ++iZ)
        {
                for (int iY = 0;  iY < dimy;  ++iY)
                {
                        for (int iX = 0;  iX < dimx;  ++iX)
                        {
                                float val = Get_P(KernelType,pq,dimx,dimy,xh,yh,zh,iX,iY,iZ);
                                if ( val != 0.0f)
                                {
                                        if (iX < min_X) min_X = iX;
                                        if (iX > max_X) max_X = iX;
                                        if (iY < min_Y) min_Y = iY;
                                        if (iY > max_Y) max_Y = iY;
                                        if (iZ < min_Z) min_Z = iZ;
                                        if (iZ > max_Z) max_Z = iZ;
                                        float abs_val = val < 0.0f ? -val : val;
                                        if (abs_val < min_abs_val) min_abs_val = abs_val;
                                        if (abs_val > max_abs_val) max_abs_val = abs_val;
                                }
                        }
                }
        }
        if (min_abs_val <= max_abs_val)
        {
                printf("x=[%d,%d] - y=[%d,%d] - z=[%d,%d] - abs_val=[%e,%e]\n",min_X,max_X,min_Y,max_Y,min_Z,max_Z,min_abs_val,max_abs_val);
        }
}

void TTIDenQ_Compute_Performance(
	int support_AVX,			// 0->SSE, 1->AVX
	double mcells,				// size of propagated volume : 1.0 -> one million cells
	double elapsed_time,			// elapsed time in seconds
	int num_threads,			// number of openmp threads
	int bsX,				// block size
	int bsY,
	int KernelType,				// 0->ISO, 1->VTI, 2->TTI
	int stage,				// 0 -> 2nd order in time, 1 -> 4th order, 1st call, 2 -> 4th order, 2nd call, 3 -> 4th order, both calls.
	double& overcomp,			// Overcompute factor caused by halos. 1.0 -> no overcompute.
	double& flops_per_cell,			// Average floating point operations per cell. Usually not an integer due to overcomp.
	double& gflops,				// Achieved average GFLOPS.
	double& eff_freq,			// Effective frequency is a measure of efficiency. Effective frequency is real frequency multiplied by efficiency.
	double& net_comp_throughput		// Net computational throughput in million cells per second.
	)
{
	overcomp = (double)((bsX+12)*(bsY+9)) / (double)(bsX*bsY);
	flops_per_cell = 0.0;
	if (stage == 0)
	{
		// 2nd order in time
		if (KernelType == 0 || KernelType == 1)
		{
			flops_per_cell  = (overcomp - 1.0) * (15.0 + 1.0) + (43.0 + 3.0);  	// DXED1, DYED1, DZED2
			flops_per_cell += 43.0;						 	// d(DXED1)dx, d(DYED1)dy, d(DZED2)dz
			flops_per_cell += 1.0;							// d(DXED1)dx + d(DYED1)dy
			if (KernelType == 1)
			{
				flops_per_cell *= 2.0;							// P and Q
				flops_per_cell += 29.0;							// Wave equation same for VTI and TTI
			}
			else
			{
				flops_per_cell += 9.0;
			}
		}
		else if (KernelType == 2)
		{
			flops_per_cell += 2.0 * (overcomp * (43.0 + 31.0) + 86.0) + 29.0;
		}
	}
	else
	{
		if (KernelType == 0 || KernelType == 1)
		{
			if ((stage & 1) == 1)
			{
				// 4th order in time - first call
				flops_per_cell += (overcomp - 1.0) * (15.0 + 1.0) + (43.0 + 3.0);  	// DXED1, DYED1, DZED2
				flops_per_cell += 43.0;						 	// d(DXED1)dx, d(DYED1)dy, d(DZED2)dz
	                        flops_per_cell += 1.0;                                                  // d(DXED1)dx + d(DYED1)dy
				if (KernelType == 1)
				{
					flops_per_cell *= 2.0;							// P and Q
					flops_per_cell += 18.0;							// Wave equation same for VTI and TTI
				}
				else
				{
					flops_per_cell += 9.0;
				}
			}
			if ((stage & 2) == 2)
			{
				// 4th order in time - second call
				flops_per_cell += (overcomp - 1.0) * (15.0 + 1.0) + (43.0 + 3.0);  	// DXED1, DYED1, DZED2
				flops_per_cell += 43.0;						 	// d(DXED1)dx, d(DYED1)dy, d(DZED2)dz
				flops_per_cell += 1.0;                                                  // d(DXED1)dx + d(DYED1)dy
				if (KernelType == 1)
				{
					flops_per_cell *= 2.0;							// P and Q
					flops_per_cell += 33.0;							// Wave equation same for VTI and TTI
				}
				else
				{
					flops_per_cell += 11.0;
				}
			}
		}
		else if (KernelType == 2)
		{
			if ((stage & 1) == 1)
			{
				// 4th order in time - first call
				flops_per_cell += 2.0 * (overcomp * (43.0 + 31.0) + 86.0) + 18.0;
			}
			if ((stage & 2) == 2)
			{
				// 4th order in time - second call
				flops_per_cell += 2.0 * (overcomp * (43.0 + 31.0) + 86.0) + 33.0;
			}
		}
	}
	gflops = mcells * flops_per_cell / (1e3 * elapsed_time);

	eff_freq = gflops / ((support_AVX ? 8.0 : 4.0) * (double)num_threads);
	net_comp_throughput = mcells / elapsed_time;
}

float Cmp_Spg_Coeff(int iX, int dimx, int absorb0, float spongecoeff, int spongewidth)
{
	int inv_iX = dimx - iX - 1;
	if (iX < spongewidth)
	{
		if (absorb0)
		{
			if (iX == 0)
			{
				return 0.0f;
			}
			else
			{
				float term = (float)(spongewidth - iX);
				return 1.0f - spongecoeff*term*term;
			}
		}
		else
		{
			return 1.0f;
		}
	}
	else if (inv_iX < spongewidth)
	{
		if (iX == 0)
		{
			return 0.0f;
		}
		else
		{
			float term = (float)(spongewidth - inv_iX);
			return 1.0f - spongecoeff*term*term;
		}
	}
	else
	{
		return 1.0f;
	}
}

void ABCsponge(
	int KernelType,
        __m128* pq,
	float* spgx,
	float* spgy,
	float* spgz,
        int dimx,
        int dimy,
        int dimz,
        int xh,
        int yh,
        int zh,
        int x0,
        int x1,
        int y0,
        int y1,
        int z0,
        int z1
        )
{
	unsigned long dimxh = dimx + 2 * xh;
	unsigned long dimyh = dimy + 2 * yh;
#pragma omp parallel for schedule(dynamic)
	for (int iZ = 0;  iZ < dimz;  ++iZ)
	{
		float spg_z = spgz[iZ];
		for (int iY = 0;  iY < dimy; ++iY)
		{
			float spg_y = spgy[iY];
			for (int iX = 0;  iX < dimx;  ++iX)
			{
				float spg_x = spgx[iX];
				float spg = spg_x * spg_y * spg_z;
				if (KernelType > 0)
				{
					unsigned long idx = ( ( (unsigned long)(iZ+zh) * (unsigned long)dimyh + (unsigned long)(iY+yh) ) * (unsigned long)dimxh + (unsigned long)(iX+xh) ) << 1;
					float* fpq = (float*)pq;
					fpq[idx] *= spg;
					fpq[idx+1] *= spg;
				}
				else
				{
					unsigned long idx = ( ( (unsigned long)(iZ+zh) * (unsigned long)dimyh + (unsigned long)(iY+yh) ) * (unsigned long)dimxh + (unsigned long)(iX+xh) );
					float* fpq = (float*)pq;
					fpq[idx] *= spg;
				}
			}
		}
	}
}

