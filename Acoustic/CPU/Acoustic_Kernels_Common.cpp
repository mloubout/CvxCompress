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

inline void ISODenQ_Get_EM1(
        int* pVelDen,
        __m128& Buoy
        )
{
        unsigned int velden1 = pVelDen[0];
        unsigned int velden2 = pVelDen[1];
        unsigned int velden3 = pVelDen[2];
        unsigned int velden4 = pVelDen[3];

        unsigned int indexDen1 = (velden1 >> SHIFTDen) & DENMASK;
        unsigned int indexDen2 = (velden2 >> SHIFTDen) & DENMASK;
        unsigned int indexDen3 = (velden3 >> SHIFTDen) & DENMASK;
        unsigned int indexDen4 = (velden4 >> SHIFTDen) & DENMASK;

        Buoy = _mm_setr_ps(lutBuoy[indexDen1], lutBuoy[indexDen2], lutBuoy[indexDen3], lutBuoy[indexDen4]);
}

inline void ISODenQ_Get_EM2_Raw(
	int* pVelDen,
	__m128& Den,
	__m128& Qatten,
	__m128& Vp2
	)
{
	/*
        unsigned int velden1 = pVelDen[0];
        unsigned int velden2 = pVelDen[1];
        unsigned int velden3 = pVelDen[2];
        unsigned int velden4 = pVelDen[3];

        unsigned int indexDen1 = (velden1 >> SHIFTDen) & DENMASK;
        unsigned int indexDen2 = (velden2 >> SHIFTDen) & DENMASK;
        unsigned int indexDen3 = (velden3 >> SHIFTDen) & DENMASK;
        unsigned int indexDen4 = (velden4 >> SHIFTDen) & DENMASK;

	unsigned int indexQ1 = (velden1 >> SHIFTQ) & QMASK;
	unsigned int indexQ2 = (velden2 >> SHIFTQ) & QMASK;
	unsigned int indexQ3 = (velden3 >> SHIFTQ) & QMASK;
	unsigned int indexQ4 = (velden4 >> SHIFTQ) & QMASK;

        unsigned int indexVpdt2_1  =  velden1 & VELMASK;
        unsigned int indexVpdt2_2  =  velden2 & VELMASK;
        unsigned int indexVpdt2_3  =  velden3 & VELMASK;
        unsigned int indexVpdt2_4  =  velden4 & VELMASK;

        Den = _mm_setr_ps(lutDen[indexDen1], lutDen[indexDen2], lutDen[indexDen3], lutDen[indexDen4]);
	Qatten = _mm_setr_ps(lutQ[indexQ1], lutQ[indexQ2], lutQ[indexQ3], lutQ[indexQ4]);  
        Vp2 = _mm_setr_ps(lutVp2[indexVpdt2_1], lutVp2[indexVpdt2_2], lutVp2[indexVpdt2_3], lutVp2[indexVpdt2_4]);
	*/
	
	__m128i velden = *((__m128i*)pVelDen);
	
	__m128i iDen = _mm_and_si128(_mm_srli_epi32(velden, SHIFTDen), _mm_set1_epi32(DENMASK));
	Den = _mm_add_ps(_mm_mul_ps(_mm_cvtepi32_ps(iDen), _mm_set1_ps(_mm_Denbinsize)), _mm_set1_ps(_mm_Denmin));
	
	__m128i iQ = _mm_and_si128(_mm_srli_epi32(velden, SHIFTQ), _mm_set1_epi32(QMASK));
	Qatten = _mm_add_ps(_mm_mul_ps(_mm_cvtepi32_ps(iQ), _mm_set1_ps(_mm_Qbinsize)), _mm_set1_ps(_mm_Qmin));

	__m128i iVel = _mm_and_si128(velden, _mm_set1_epi32(VELMASK));
	__m128 Vel = _mm_add_ps(_mm_mul_ps(_mm_cvtepi32_ps(iVel), _mm_set1_ps(_mm_Velbinsize)), _mm_set1_ps(_mm_Velmin));
	Vp2 = _mm_mul_ps(_mm_mul_ps(Vel, Vel), _mm_set1_ps(_mm_dt2));
}

inline void ISODenQ_Get_EM2(
        int* pVelDen,
	__m128& Qatten,
        __m128& C33
        )
{
	__m128 Den, Vp2;
	ISODenQ_Get_EM2_Raw(pVelDen,Den,Qatten,Vp2);

        C33 = _mm_mul_ps(Den, Vp2);
}

inline void VTIDenQ_Get_EM1(
        int* pDenAng,
        __m128& Buoy
        )
{
        unsigned int denang1 = pDenAng[0];
        unsigned int denang2 = pDenAng[1];
        unsigned int denang3 = pDenAng[2];
        unsigned int denang4 = pDenAng[3];

        unsigned int indexDen1 = (denang1 >> SHIFTDen) & DENMASK;
        unsigned int indexDen2 = (denang2 >> SHIFTDen) & DENMASK;
        unsigned int indexDen3 = (denang3 >> SHIFTDen) & DENMASK;
        unsigned int indexDen4 = (denang4 >> SHIFTDen) & DENMASK;

        Buoy = _mm_setr_ps(lutBuoy[indexDen1], lutBuoy[indexDen2], lutBuoy[indexDen3], lutBuoy[indexDen4]);
}

inline void TTIDenQ_Get_EM1(
        int* pDenAng,
        __m128& Buoy,
        __m128& cDip,
        __m128& sDip,
        __m128& cAzm,
        __m128& sAzm
        )
{
        unsigned int denang1 = pDenAng[0];
        unsigned int denang2 = pDenAng[1];
        unsigned int denang3 = pDenAng[2];
        unsigned int denang4 = pDenAng[3];

        unsigned int indexDip1 =  denang1              & DIPMASK;
        unsigned int indexAzm1 = (denang1 >> SHIFTAzm) & AZMMASK;
        unsigned int indexDen1 = (denang1 >> SHIFTDen) & DENMASK;

        unsigned int indexDip2 =  denang2              & DIPMASK;
        unsigned int indexAzm2 = (denang2 >> SHIFTAzm) & AZMMASK;
        unsigned int indexDen2 = (denang2 >> SHIFTDen) & DENMASK;

        unsigned int indexDip3 =  denang3              & DIPMASK;
        unsigned int indexAzm3 = (denang3 >> SHIFTAzm) & AZMMASK;
        unsigned int indexDen3 = (denang3 >> SHIFTDen) & DENMASK;

        unsigned int indexDip4 =  denang4              & DIPMASK;
        unsigned int indexAzm4 = (denang4 >> SHIFTAzm) & AZMMASK;
        unsigned int indexDen4 = (denang4 >> SHIFTDen) & DENMASK;

        Buoy = _mm_setr_ps(lutBuoy[indexDen1], lutBuoy[indexDen2], lutBuoy[indexDen3], lutBuoy[indexDen4]);
        cDip = _mm_setr_ps(lutcDip[indexDip1], lutcDip[indexDip2], lutcDip[indexDip3], lutcDip[indexDip4]);
        sDip = _mm_setr_ps(lutsDip[indexDip1], lutsDip[indexDip2], lutsDip[indexDip3], lutsDip[indexDip4]);
        cAzm = _mm_setr_ps(lutcAzm[indexAzm1], lutcAzm[indexAzm2], lutcAzm[indexAzm3], lutcAzm[indexAzm4]);
        sAzm = _mm_setr_ps(lutsAzm[indexAzm1], lutsAzm[indexAzm2], lutsAzm[indexAzm3], lutsAzm[indexAzm4]);
}

inline void TTIDenQ_Get_EM2_Raw(
	int* pVelAnis,
        int* pDenAng,
	__m128& Den,
	__m128& Qatten,
	__m128& Vp2,
	__m128& C44C33,
	__m128& Del,
	__m128& Eps
	)
{
	/*
        unsigned int denang1 = pDenAng[0];
        unsigned int denang2 = pDenAng[1];
        unsigned int denang3 = pDenAng[2];
        unsigned int denang4 = pDenAng[3];

        unsigned int indexDen1 = (denang1 >> SHIFTDen) & DENMASK;
        unsigned int indexDen2 = (denang2 >> SHIFTDen) & DENMASK;
        unsigned int indexDen3 = (denang3 >> SHIFTDen) & DENMASK;
        unsigned int indexDen4 = (denang4 >> SHIFTDen) & DENMASK;

	unsigned int indexQ1 = (denang1 >> SHIFTQ) & QMASK;
	unsigned int indexQ2 = (denang2 >> SHIFTQ) & QMASK;
	unsigned int indexQ3 = (denang3 >> SHIFTQ) & QMASK;
	unsigned int indexQ4 = (denang4 >> SHIFTQ) & QMASK;

        Den = _mm_setr_ps(lutDen[indexDen1], lutDen[indexDen2], lutDen[indexDen3], lutDen[indexDen4]);
	Qatten = _mm_setr_ps(lutQ[indexQ1], lutQ[indexQ2], lutQ[indexQ3], lutQ[indexQ4]);  

        unsigned int velanis1 = pVelAnis[0];
        unsigned int velanis2 = pVelAnis[1];
        unsigned int velanis3 = pVelAnis[2];
        unsigned int velanis4 = pVelAnis[3];

        unsigned int indexVpdt2_1  =  velanis1                 & VELMASK;
        unsigned int indexEps_1    = (velanis1 >> SHIFTEps   ) & EPSMASK;
        unsigned int indexDel_1    = (velanis1 >> SHIFTDel   ) & DELMASK;
        unsigned int indexC44C33_1 = (velanis1 >> SHIFTC44C33) & C44C33MASK;

        unsigned int indexVpdt2_2  =  velanis2                 & VELMASK;
        unsigned int indexEps_2    = (velanis2 >> SHIFTEps   ) & EPSMASK;
        unsigned int indexDel_2    = (velanis2 >> SHIFTDel   ) & DELMASK;
        unsigned int indexC44C33_2 = (velanis2 >> SHIFTC44C33) & C44C33MASK;

        unsigned int indexVpdt2_3  =  velanis3                 & VELMASK;
        unsigned int indexEps_3    = (velanis3 >> SHIFTEps   ) & EPSMASK;
        unsigned int indexDel_3    = (velanis3 >> SHIFTDel   ) & DELMASK;
        unsigned int indexC44C33_3 = (velanis3 >> SHIFTC44C33) & C44C33MASK;

        unsigned int indexVpdt2_4  =  velanis4                 & VELMASK;
        unsigned int indexEps_4    = (velanis4 >> SHIFTEps   ) & EPSMASK;
        unsigned int indexDel_4    = (velanis4 >> SHIFTDel   ) & DELMASK;
        unsigned int indexC44C33_4 = (velanis4 >> SHIFTC44C33) & C44C33MASK;

        Vp2 = _mm_setr_ps(lutVp2[indexVpdt2_1], lutVp2[indexVpdt2_2], lutVp2[indexVpdt2_3], lutVp2[indexVpdt2_4]);
        C44C33 = _mm_setr_ps(lutc44c33[indexC44C33_1], lutc44c33[indexC44C33_2], lutc44c33[indexC44C33_3], lutc44c33[indexC44C33_4]);
        Del = _mm_setr_ps(lutDel[indexDel_1], lutDel[indexDel_2], lutDel[indexDel_3], lutDel[indexDel_4]);
        Eps = _mm_setr_ps(lutEps[indexEps_1], lutEps[indexEps_2], lutEps[indexEps_3], lutEps[indexEps_4]);
	*/
	
	__m128i denang = *((__m128i*)pDenAng);
	
	__m128i iDen = _mm_and_si128(_mm_srli_epi32(denang, SHIFTDen), _mm_set1_epi32(DENMASK));
	Den = _mm_add_ps(_mm_mul_ps(_mm_cvtepi32_ps(iDen), _mm_load1_ps(&_mm_Denbinsize)), _mm_load1_ps(&_mm_Denmin));
	
	__m128i iQ = _mm_and_si128(_mm_srli_epi32(denang, SHIFTQ), _mm_set1_epi32(QMASK));
	Qatten = _mm_add_ps(_mm_mul_ps(_mm_cvtepi32_ps(iQ), _mm_load1_ps(&_mm_Qbinsize)), _mm_load1_ps(&_mm_Qmin));

	__m128i velanis = *((__m128i*)pVelAnis);

	__m128i iVel = _mm_and_si128(velanis, _mm_set1_epi32(VELMASK));
	__m128 Vel = _mm_add_ps(_mm_mul_ps(_mm_cvtepi32_ps(iVel), _mm_load1_ps(&_mm_Velbinsize)), _mm_load1_ps(&_mm_Velmin));
	Vp2 = _mm_mul_ps(_mm_mul_ps(Vel, Vel), _mm_load1_ps(&_mm_dt2));

	__m128i iEps = _mm_and_si128(_mm_srli_epi32(velanis, SHIFTEps), _mm_set1_epi32(EPSMASK));
	Eps = _mm_add_ps(_mm_mul_ps(_mm_cvtepi32_ps(iEps), _mm_load1_ps(&_mm_Epsbinsize)), _mm_load1_ps(&_mm_Epsmin));
	Eps = _mm_add_ps(_mm_set1_ps(1.0f), _mm_add_ps(Eps, Eps));  // ...from Precompute_Anisotropy_Parameters

	__m128i iDel = _mm_and_si128(_mm_srli_epi32(velanis, SHIFTDel), _mm_set1_epi32(DELMASK));
	Del = _mm_add_ps(_mm_mul_ps(_mm_cvtepi32_ps(iDel), _mm_load1_ps(&_mm_Delbinsize)), _mm_load1_ps(&_mm_Delmin));
	Del = _mm_add_ps(Del, Del);  // ...from Precompute_Anisotropy_Parameters

	__m128i iC44C33 = _mm_and_si128(_mm_srli_epi32(velanis, SHIFTC44C33), _mm_set1_epi32(C44C33MASK));
	C44C33 = _mm_add_ps(_mm_mul_ps(_mm_cvtepi32_ps(iC44C33), _mm_load1_ps(&_mm_C44C33binsize)), _mm_load1_ps(&_mm_C44C33min));
}

inline void TTIDenQ_Get_EM2(
        int* pVelAnis,
        int* pDenAng,
	__m128& Qatten,
        __m128& C66C44_01,
        __m128& C66C44_23,
        __m128& C44C33_01,
        __m128& C44C33_23,
        __m128& C55_01,
        __m128& C55_23
        )
{
	__m128 Den, Vp2, C44C33, Del, Eps;
	TTIDenQ_Get_EM2_Raw(pVelAnis,pDenAng,Den,Qatten,Vp2,C44C33,Del,Eps);

        __m128 C33 = _mm_mul_ps(Den, Vp2);  // C33 == Bulk modulus
        __m128 C44 = _mm_mul_ps(C33, C44C33);
        __m128 C33mC44 = _mm_sub_ps(C33, C44);

	// value in variable Del is actually 2.0 * Del
        __m128 C13pC44 = _mm_mul_ps(Del, C33);
        C13pC44 = _mm_add_ps(C13pC44, C33mC44);
        C13pC44 = _mm_mul_ps(C33mC44, C13pC44);
        C13pC44 = _mm_sqrt_ps(C13pC44);  // NB reciprocal square root approximation cannot be used here. Proven experimentally to be too inaccurate.

	// value in variable Eps is actually 1.0 + 2.0 * Eps
        __m128 C66 = _mm_mul_ps(Eps, C33);

        C55_01 = _mm_shuffle_ps(C13pC44, C13pC44, 0x50);
        C55_23 = _mm_shuffle_ps(C13pC44, C13pC44, 0xFA);

        __m128 C66C44_lo = _mm_shuffle_ps(C66, C44, 0x44);
        __m128 C66C44_hi = _mm_shuffle_ps(C66, C44, 0xEE);

        C66C44_01 = _mm_shuffle_ps(C66C44_lo, C66C44_lo, 0xD8);
        C66C44_23 = _mm_shuffle_ps(C66C44_hi, C66C44_hi, 0xD8);

        __m128 C44C33_lo = _mm_shuffle_ps(C44, C33, 0x44);
        __m128 C44C33_hi = _mm_shuffle_ps(C44, C33, 0xEE);

        C44C33_01 = _mm_shuffle_ps(C44C33_lo, C44C33_lo, 0xD8);
        C44C33_23 = _mm_shuffle_ps(C44C33_hi, C44C33_hi, 0xD8);

	/*
	float _c33, _c44, _c55, _c66;
	_c33 = *(((float*)&C33));
	_c44 = *(((float*)&C44));
	_c55 = *(((float*)&C13pC44));
	_c66 = *(((float*)&C66));

	static int been_there_2 = 0;
	if (!been_there_2)
	{
		been_there_2 = 1;
		printf("C33 = %f, C44 = %f, C55 = %f, C66 = %f\n",_c33,_c44,_c55,_c66);
	}
	*/
}

void ISODenQ_Shift_Queues(
	int QMaxShift,
	int QLen,
	__m128* Q,
	__m128*& V4,
	int& QShift
	)
{
	if (QShift > 0)
	{
		V4 -= 1;
		--QShift;
	}
	else
	{
		int QOff = QMaxShift;
		for (int iQ = QLen-1;  iQ >= 0;  --iQ)
		{
			Q[iQ+QOff] = Q[iQ];
		}
		QShift = QMaxShift - 1;
		V4 = Q + QShift;
	}
}

void VTIDenQ_Shift_Queues(
	int QMaxShift,
	int QLen,
	__m128* Q,
	__m128*& V4,
	int& QShift
	)
{
	if (QShift > 0)
	{
		V4 -= 2;
		--QShift;
	}
	else
	{
		int QOff = QMaxShift << 1;
		for (int iQ = QLen-1;  iQ >= 0;  --iQ)
		{
			Q[iQ+QOff] = Q[iQ];
		}
		QShift = QMaxShift - 1;
		V4 = Q + QShift*2;
	}
}

void TTIDenQ_Shift_Queues(
	int QMaxShift,
	int QLen,
	__m128* Q,
	__m128*& V4,
	int& QShift
	)
{
	if (QShift > 0)
	{
		V4 -= 4;
		--QShift;
	}
	else
	{
		int QOff = QMaxShift << 2;
		for (int iQ = QLen-1;  iQ >= 0;  --iQ)
		{
			Q[iQ+QOff] = Q[iQ];
		}
		QShift = QMaxShift - 1;
		V4 = Q + QShift*4;
	}
}

static unsigned long _BufSize = 0;
static __m128* _Buf = 0L;

void ISODenQ_Comp_Buf_Size(
	int bsX,
	int bsY,
	int dimz,
	unsigned long& stripe_m128,
	unsigned long& dyed_buf_m128,
	unsigned long& spg_buf_m128,
	int& QMaxShift,
	unsigned long& QLen,
	unsigned long& QSize,
	unsigned long& BufSize_m128,
	unsigned long& NetBufSize_m128
	)
{
	unsigned long bsX_2 = bsX >> 1;
	unsigned long bsX_4 = bsX >> 2;
	stripe_m128 = (bsX_2+6);
	dyed_buf_m128 = bsX_2 * (bsY+9);
	spg_buf_m128 = bsX_4 + bsY + 11;
	QMaxShift = (dimz+14);
	QLen = bsX_4 * (bsY+1) * 17;
	QSize = QLen + (unsigned long)(QMaxShift);
	BufSize_m128 = stripe_m128 + dyed_buf_m128 + spg_buf_m128 + QSize;
	// round up to whole number of cache pages. add extra cache page for good measure.
	// performance drops into the toilet if two threads keep updating the same cache page.
	BufSize_m128 = (((BufSize_m128 + 3) >> 2) << 2) + 4;
	NetBufSize_m128 = BufSize_m128 - (QSize - QLen);  // Portion of buffer that should remain in cache all the time.
}

void VTIDenQ_Comp_Buf_Size(
	int bsX,
	int bsY,
	int dimz,
	unsigned long& stripe_m128,
	unsigned long& dyed_buf_m128,
	unsigned long& spg_buf_m128,
	int& QMaxShift,
	unsigned long& QLen,
	unsigned long& QSize,
	unsigned long& BufSize_m128,
	unsigned long& NetBufSize_m128
	)
{
	unsigned long bsX_2 = bsX >> 1;
	unsigned long bsX_4 = bsX >> 2;
	stripe_m128 = 2 * (bsX_2+6);
	dyed_buf_m128 = 2 * bsX_2 * (bsY+9);
	spg_buf_m128 = bsX_4 + bsY + 11;
	QMaxShift = (dimz+14);
	QLen = bsX_2 * (bsY+1) * 17;
	QSize = QLen + (unsigned long)(QMaxShift<<1);
	BufSize_m128 = stripe_m128 + dyed_buf_m128 + spg_buf_m128 + QSize;
	// round up to whole number of cache pages. add extra cache page for good measure.
	// performance drops into the toilet if two threads keep updating the same cache page.
	BufSize_m128 = (((BufSize_m128 + 3) >> 2) << 2) + 4;
	NetBufSize_m128 = BufSize_m128 - (QSize - QLen);  // Portion of buffer that should remain in cache all the time.
}

void TTIDenQ_Comp_Buf_Size(
	int bsX,
	int bsY,
	int dimz,
	unsigned long& stripe_m128,
	unsigned long& dyed_buf_m128,
	unsigned long& spg_buf_m128,
	int& QMaxShift,
	unsigned long& QLen,
	unsigned long& QSize,
	unsigned long& BufSize_m128,
	unsigned long& NetBufSize_m128
	)
{
	unsigned long bsX_2 = bsX >> 1;
	unsigned long bsX_4 = bsX >> 2;
	stripe_m128 = 2 * (bsX_2+6) * 2;
	dyed_buf_m128 = 2 * bsX_2 * (bsY+9) * 2;
	spg_buf_m128 = bsX_4 + bsY + 11;
	QMaxShift = (dimz+14);
	QLen = bsX_2 * (bsY+1) * 34;
	QSize = QLen + (unsigned long)(QMaxShift<<2);
	BufSize_m128 = stripe_m128 + dyed_buf_m128 + spg_buf_m128 + QSize;
	// round up to whole number of cache pages. add extra cache page for good measure.
	// performance drops into the toilet if two threads keep updating the same cache page.
	BufSize_m128 = (((BufSize_m128 + 3) >> 2) << 2) + 4;
	NetBufSize_m128 = BufSize_m128 - (QSize - QLen);  // Portion of buffer that should remain in cache all the time.
}

void ISODenQ_Comp_Cache_Usage(
	int bsX,
	int bsY,
	unsigned long& Cache_Usage_Per_Core_KB
	)
{
	// estimate size of compute buffer in __m128 words (16 bytes).
	unsigned long stripe_m128, dyed_buf_m128, spg_buf_m128, QLen, QSize, BufSize_m128, NetBufSize_m128;
        int QMaxShift;
	ISODenQ_Comp_Buf_Size(bsX, bsY, 128, stripe_m128, dyed_buf_m128, spg_buf_m128, QMaxShift, QLen, QSize, BufSize_m128, NetBufSize_m128);
	
	// estimate space used by input data
	// ..input wavefield
	int npages_X = (bsX+12+9+7) >> 4;
	int XY_plane_size = npages_X * (bsY+18) * 3;
	int npages_XZ = (bsX+5+4+7) >> 4;
	int Z_size = npages_XZ * (bsY+9) * 9;
	// ..earth model
	int npages_EM1 = (((bsX+12) + 15) >> 4) * (bsY+9);
	int npages_EM2 = ((bsX+7) >> 4) * bsY;
	// ,.pq, rs in wave equation
	int npages_WE = ((bsX+3) >> 3) * bsY * 2;
	// ..total
	int npages = XY_plane_size + Z_size + npages_EM1 + npages_EM2 + npages_WE;
	//printf("XY_plane_size = %d, Z_size = %d, npages_EM1 = %d, npages_EM2 = %d, npages_WE = %d\n",XY_plane_size,Z_size,npages_EM1,npages_EM2,npages_WE);

	// total up
	Cache_Usage_Per_Core_KB = ((NetBufSize_m128 * 16) + (unsigned long)(npages * 64)) / 1024;
}

void VTIDenQ_Comp_Cache_Usage(
	int bsX,
	int bsY,
	unsigned long& Cache_Usage_Per_Core_KB
	)
{
	// estimate size of compute buffer in __m128 words (16 bytes).
	unsigned long stripe_m128, dyed_buf_m128, spg_buf_m128, QLen, QSize, BufSize_m128, NetBufSize_m128;
        int QMaxShift;
	VTIDenQ_Comp_Buf_Size(bsX, bsY, 128, stripe_m128, dyed_buf_m128, spg_buf_m128, QMaxShift, QLen, QSize, BufSize_m128, NetBufSize_m128);
	
	// estimate space used by input data
	// ..input wavefield
	int npages_X = (bsX+12+9+7) >> 3;
	int XY_plane_size = npages_X * (bsY+18) * 3;
	int npages_XZ = (bsX+5+4+7) >> 3;
	int Z_size = npages_XZ * (bsY+9) * 9;
	// ..earth model
	int npages_EM1 = (((bsX+12) + 15) >> 4) * (bsY+9);
	int npages_EM2 = ((bsX+7) >> 3) * bsY;
	// ,.pq, rs in wave equation
	int npages_WE = ((bsX+3) >> 2) * bsY * 2;
	// ..total
	int npages = XY_plane_size + Z_size + npages_EM1 + npages_EM2 + npages_WE;
	//printf("XY_plane_size = %d, Z_size = %d, npages_EM1 = %d, npages_EM2 = %d, npages_WE = %d\n",XY_plane_size,Z_size,npages_EM1,npages_EM2,npages_WE);

	// total up
	Cache_Usage_Per_Core_KB = ((NetBufSize_m128 * 16) + (unsigned long)(npages * 64)) / 1024;
}

void TTIDenQ_Comp_Cache_Usage(
	int bsX,
	int bsY,
	unsigned long& Cache_Usage_Per_Core_KB
	)
{
	// estimate size of compute buffer in __m128 words (16 bytes).
	unsigned long stripe_m128, dyed_buf_m128, spg_buf_m128, QLen, QSize, BufSize_m128, NetBufSize_m128;
        int QMaxShift;
	TTIDenQ_Comp_Buf_Size(bsX, bsY, 128, stripe_m128, dyed_buf_m128, spg_buf_m128, QMaxShift, QLen, QSize, BufSize_m128, NetBufSize_m128);
	
	// estimate space used by input data
	// ..input wavefield
	int npages_X = (bsX+12+9+7) >> 3;
	int XY_plane_size = npages_X * (bsY+18) * 3;
	int npages_XZ = (bsX+5+4+7) >> 3;
	int Z_size = npages_XZ * (bsY+9) * 9;
	// ..earth model
	int npages_EM1 = (((bsX+12) + 15) >> 4) * (bsY+9);
	int npages_EM2 = ((bsX+7) >> 3) * bsY;
	// ,.pq, rs in wave equation
	int npages_WE = ((bsX+3) >> 2) * bsY * 2;
	// ..total
	int npages = XY_plane_size + Z_size + npages_EM1 + npages_EM2 + npages_WE;
	//printf("XY_plane_size = %d, Z_size = %d, npages_EM1 = %d, npages_EM2 = %d, npages_WE = %d\n",XY_plane_size,Z_size,npages_EM1,npages_EM2,npages_WE);

	// total up
	Cache_Usage_Per_Core_KB = ((NetBufSize_m128 * 16) + (unsigned long)(npages * 64)) / 1024;
}
