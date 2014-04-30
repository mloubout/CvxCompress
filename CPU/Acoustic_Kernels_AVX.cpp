#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <immintrin.h>
#include <omp.h>
#include <time.h>

extern __m128 A5h_A4h;
extern __m128 A4h_A3h;
extern __m128 A3h_A2h;
extern __m128 A2h_A1h;
extern __m128 A5h_mA5h;
extern __m128 mA1h_A1h;

extern __m128 A5h_A4h_A3h_A2h;
extern __m128 A1h_A1h_A1h_A1h;

extern __m128 A5_A4_A3_A2;
extern __m128 A1_A1_A1_A1;

extern float A1h, A2h, A3h, A4h, A5h;
extern float B1h, B2h, B3h, B4h;
extern float C1h, C2h, C3h;
extern float D1h, D2h;
extern float E1h;

extern float A1, A2, A3, A4, A5;
extern float B1, B2, B3, B4;
extern float C1, C2, C3;
extern float D1, D2;
extern float E1;

extern float *lutVp2, *lutEps, *lutDel, *lutDen, *lutBuoy, *lutsDip, *lutcDip, *lutsAzm, *lutcAzm, *lutQ, *lutc44c33;

extern float _mm_Velmin;
extern float _mm_Velbinsize; 
extern float _mm_dt2; 
extern float _mm_Epsmin; 
extern float _mm_Epsbinsize; 
extern float _mm_Delmin; 
extern float _mm_Delbinsize; 
extern float _mm_Denmin; 
extern float _mm_Denbinsize; 
extern float _mm_Qmin; 
extern float _mm_Qbinsize; 
extern float _mm_C44C33min; 
extern float _mm_C44C33binsize; 

namespace {

#include "Acoustic_Kernels_Common.cpp"

inline __m256 conc(__m128 a, __m128 b)
{
	__m256 c = _mm256_castps128_ps256(a);
	c = _mm256_insertf128_ps(c, b, 1);
	return c;
}

// TO-DO: This is an actual AVX2 instruction. Delete when AVX2 is available.
inline __m256 _mm256_alignr_epi8(__m256 a, __m256 b, int n)
{
	__m128 v0 = _mm256_extractf128_ps(a, 1);
	__m128 v1 = _mm256_castps256_ps128(a);
	__m128 w0 = _mm256_extractf128_ps(b, 1);
	__m128 w1 = _mm256_castps256_ps128(b);
	__m256 yy = _mm256_castps128_ps256(_mm_castsi128_ps(_mm_alignr_epi8(_mm_castps_si128(v0), _mm_castps_si128(v1), 8)));
	yy = _mm256_insertf128_ps(yy, _mm_castsi128_ps(_mm_alignr_epi8(_mm_castps_si128(w0), _mm_castps_si128(w1), 8)), 1);
	return yy;
}

inline void AVX_Cmp_DDX(
	__m256* pq,
	__m256& ddx_0123
	)
{
	__m256 v0 = pq[-1];
	__m256 v1 = pq[ 0];
	__m256 v2 = pq[ 1];
	__m256 v3 = pq[ 2];

	__m256 m3 = _mm256_permute2f128_ps(v0, v1, 33);
	__m256 m4 = _mm256_alignr_epi8(v0, m3, 8);
	ddx_0123 = _mm256_mul_ps(_mm256_sub_ps(v2, m4), _mm256_broadcast_ss(&A4h));

	__m256 p2 = _mm256_permute2f128_ps(v1, v2, 33);

	__m256 p1 = _mm256_alignr_epi8(v1, p2, 8);
	ddx_0123 = _mm256_add_ps(ddx_0123, _mm256_mul_ps(_mm256_sub_ps(p1, v1), _mm256_broadcast_ss(&A1h)));

	__m256 m2 = _mm256_permute2f128_ps(m4, p1, 33);
	ddx_0123 = _mm256_add_ps(ddx_0123, _mm256_mul_ps(_mm256_sub_ps(p2, m2), _mm256_broadcast_ss(&A2h)));

	__m256 p6 = _mm256_permute2f128_ps(v2, v3, 33);
	__m256 p5 = _mm256_alignr_epi8(v2, p6, 8);
	ddx_0123 = _mm256_add_ps(ddx_0123, _mm256_mul_ps(_mm256_sub_ps(p5, v0), _mm256_broadcast_ss(&A5h)));

	__m256 p3 = _mm256_permute2f128_ps(p1, p5, 33);
	ddx_0123 = _mm256_add_ps(ddx_0123, _mm256_mul_ps(_mm256_sub_ps(p3, m3), _mm256_broadcast_ss(&A3h)));
}

inline void AVX_Cmp_DDX_2X(
	__m256* pq,
	__m256& ddx_x0_0123,
	__m256& ddx_x1_0123
	)
{
	__m256 v0 = pq[-1];
	__m256 v1 = pq[ 0];
	__m256 v2 = pq[ 1];
	__m256 v3 = pq[ 2];
	__m256 v4 = pq[ 3];

	__m256 m3 = _mm256_permute2f128_ps(v0, v1, 33);
	__m256 m4 = _mm256_alignr_epi8(v0, m3, 8);
	__m256 c4h = _mm256_broadcast_ss(&A4h);
	ddx_x0_0123 = _mm256_mul_ps(_mm256_sub_ps(v2, m4), c4h);
	__m256 p2 = _mm256_permute2f128_ps(v1, v2, 33);
	__m256 p1 = _mm256_alignr_epi8(v1, p2, 8);
	ddx_x1_0123 = _mm256_mul_ps(_mm256_sub_ps(v3, p1), c4h);

	__m256 c1h = _mm256_broadcast_ss(&A1h);
	ddx_x0_0123 = _mm256_add_ps(ddx_x0_0123, _mm256_mul_ps(_mm256_sub_ps(p1, v1), c1h));
	__m256 p6 = _mm256_permute2f128_ps(v2, v3, 33);
	__m256 p5 = _mm256_alignr_epi8(v2, p6, 8);
	ddx_x1_0123 = _mm256_add_ps(ddx_x1_0123, _mm256_mul_ps(_mm256_sub_ps(p5, v2), c1h));

	__m256 m2 = _mm256_permute2f128_ps(m4, p1, 33);
	__m256 c2h = _mm256_broadcast_ss(&A2h);
	ddx_x0_0123 = _mm256_add_ps(ddx_x0_0123, _mm256_mul_ps(_mm256_sub_ps(p2, m2), c2h));
	__m256 p3 = _mm256_permute2f128_ps(p1, p5, 33);
	ddx_x1_0123 = _mm256_add_ps(ddx_x1_0123, _mm256_mul_ps(_mm256_sub_ps(p6, p3), c2h));

	__m256 c5h = _mm256_broadcast_ss(&A5h);
	ddx_x0_0123 = _mm256_add_ps(ddx_x0_0123, _mm256_mul_ps(_mm256_sub_ps(p5, v0), c5h));
	__m256 p6_x1 = _mm256_permute2f128_ps(v3, v4, 33);
	__m256 p5_x1 = _mm256_alignr_epi8(v3, p6_x1, 8);
	ddx_x1_0123 = _mm256_add_ps(ddx_x1_0123, _mm256_mul_ps(_mm256_sub_ps(p5_x1, v1), c5h));

	__m256 c3h = _mm256_broadcast_ss(&A3h);
	ddx_x0_0123 = _mm256_add_ps(ddx_x0_0123, _mm256_mul_ps(_mm256_sub_ps(p3, m3), c3h));
	__m256 p3_x1 = _mm256_permute2f128_ps(p5, p5_x1, 33);
	ddx_x1_0123 = _mm256_add_ps(ddx_x1_0123, _mm256_mul_ps(_mm256_sub_ps(p3_x1, p2), c3h));
}

inline void AVX_Cmp_DDX_1x_stride4(
	__m256* pq,
	__m256& ddx_0123
	)
{
	__m256 v0 = pq[-4];
	__m256 v1 = pq[-2];
	__m256 v2 = pq[ 0];
	__m256 v3 = pq[ 2];

	__m256 m2 = _mm256_permute2f128_ps(v1, v2, 33);
	__m256 p3 = _mm256_permute2f128_ps(v2, v3, 33);
	__m256 m3 = _mm256_alignr_epi8(v1, m2, 8);
	ddx_0123 = _mm256_mul_ps(_mm256_sub_ps(p3, m3), _mm256_broadcast_ss(&A3h));

	__m256 p2 = _mm256_alignr_epi8(v2, p3, 8);
	ddx_0123 = _mm256_add_ps(ddx_0123, _mm256_mul_ps(_mm256_sub_ps(p2, m2), _mm256_broadcast_ss(&A2h)));

	__m256 p4 = _mm256_alignr_epi8(p3, v3, 8);
	ddx_0123 = _mm256_add_ps(ddx_0123, _mm256_mul_ps(_mm256_sub_ps(p4, v1), _mm256_broadcast_ss(&A4h)));

	__m256 m1 = _mm256_permute2f128_ps(m3, p2, 33);
	ddx_0123 = _mm256_add_ps(ddx_0123, _mm256_mul_ps(_mm256_sub_ps(v2, m1), _mm256_broadcast_ss(&A1h)));

	__m256 m6 = _mm256_permute2f128_ps(v0, v1, 33);
	__m256 m5 = _mm256_alignr_epi8(m6, v1, 8);
	ddx_0123 = _mm256_add_ps(ddx_0123, _mm256_mul_ps(_mm256_sub_ps(v3, m5), _mm256_broadcast_ss(&A5h)));
}

inline void AVX_Cmp_DDX_1x_stride8(
	__m256* pq,
	__m256& ddx_0123
	)
{
	__m256 v0 = pq[-8];
	__m256 v1 = pq[-4];
	__m256 v2 = pq[ 0];
	__m256 v3 = pq[ 4];

	__m256 m2 = _mm256_permute2f128_ps(v1, v2, 33);
	__m256 p3 = _mm256_permute2f128_ps(v2, v3, 33);
	__m256 m3 = _mm256_alignr_epi8(v1, m2, 8);
	ddx_0123 = _mm256_mul_ps(_mm256_sub_ps(p3, m3), _mm256_broadcast_ss(&A3h));

	__m256 p2 = _mm256_alignr_epi8(v2, p3, 8);
	ddx_0123 = _mm256_add_ps(ddx_0123, _mm256_mul_ps(_mm256_sub_ps(p2, m2), _mm256_broadcast_ss(&A2h)));

	__m256 p4 = _mm256_alignr_epi8(p3, v3, 8);
	ddx_0123 = _mm256_add_ps(ddx_0123, _mm256_mul_ps(_mm256_sub_ps(p4, v1), _mm256_broadcast_ss(&A4h)));

	__m256 m1 = _mm256_permute2f128_ps(m3, p2, 33);
	ddx_0123 = _mm256_add_ps(ddx_0123, _mm256_mul_ps(_mm256_sub_ps(v2, m1), _mm256_broadcast_ss(&A1h)));

	__m256 m6 = _mm256_permute2f128_ps(v0, v1, 33);
	__m256 m5 = _mm256_alignr_epi8(m6, v1, 8);
	ddx_0123 = _mm256_add_ps(ddx_0123, _mm256_mul_ps(_mm256_sub_ps(v3, m5), _mm256_broadcast_ss(&A5h)));
}

inline void AVX_Cmp_DDY(
	__m256* pq,
	int stride_y_m256,
	__m256& ddy_0123
	)
{
	/*
	__m256 fact = _mm256_broadcast_ss(&A5h);
	ddy_0123 = _mm256_mul_ps(fact, pq[5*stride_y_m256]);
	ddy_0123 = _mm256_sub_ps(ddy_0123, _mm256_mul_ps(fact, pq[-4*stride_y_m256]));
	fact = _mm256_broadcast_ss(&A4h);
	ddy_0123 = _mm256_add_ps(ddy_0123, _mm256_mul_ps(fact, pq[4*stride_y_m256]));
	ddy_0123 = _mm256_sub_ps(ddy_0123, _mm256_mul_ps(fact, pq[-3*stride_y_m256]));
	fact = _mm256_broadcast_ss(&A3h);
	ddy_0123 = _mm256_add_ps(ddy_0123, _mm256_mul_ps(fact, pq[3*stride_y_m256]));
        ddy_0123 = _mm256_sub_ps(ddy_0123, _mm256_mul_ps(fact, pq[-2*stride_y_m256]));
	fact = _mm256_broadcast_ss(&A2h);
	ddy_0123 = _mm256_add_ps(ddy_0123, _mm256_mul_ps(fact, pq[2*stride_y_m256]));
        ddy_0123 = _mm256_sub_ps(ddy_0123, _mm256_mul_ps(fact, pq[-stride_y_m256]));
	fact = _mm256_broadcast_ss(&A1h);
	ddy_0123 = _mm256_add_ps(ddy_0123, _mm256_mul_ps(fact, pq[stride_y_m256]));
        ddy_0123 = _mm256_sub_ps(ddy_0123, _mm256_mul_ps(fact, pq[0]));
	*/
	
	ddy_0123 = _mm256_mul_ps(_mm256_broadcast_ss(&A5h), _mm256_sub_ps(pq[5*stride_y_m256], pq[-4*stride_y_m256]));
	ddy_0123 = _mm256_add_ps(ddy_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A4h), _mm256_sub_ps(pq[4*stride_y_m256], pq[-3*stride_y_m256])));
	ddy_0123 = _mm256_add_ps(ddy_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A3h), _mm256_sub_ps(pq[3*stride_y_m256], pq[-2*stride_y_m256])));
	ddy_0123 = _mm256_add_ps(ddy_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A2h), _mm256_sub_ps(pq[2*stride_y_m256], pq[-1*stride_y_m256])));
	ddy_0123 = _mm256_add_ps(ddy_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A1h), _mm256_sub_ps(pq[stride_y_m256], pq[0])));
}

inline void AVX_Cmp_DDZ(
        __m256* pq,
        int stride_z_m256,
        __m256& ddz_0123
        )
{
        ddz_0123 = _mm256_mul_ps(_mm256_broadcast_ss(&A5), _mm256_sub_ps(pq[5*stride_z_m256], pq[-4*stride_z_m256]));
        ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A4), _mm256_sub_ps(pq[4*stride_z_m256], pq[-3*stride_z_m256])));
        ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A3), _mm256_sub_ps(pq[3*stride_z_m256], pq[-2*stride_z_m256])));
        ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A2), _mm256_sub_ps(pq[2*stride_z_m256], pq[-1*stride_z_m256])));
        ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A1), _mm256_sub_ps(pq[stride_z_m256], pq[0])));
}

inline void AVX_Cmp_DDZ_2Z(
        __m256* pq,
        int stride_z_m256,
        __m256& ddz_z0_0123,
	__m256& ddz_z1_0123
        )
{
	__m256 m5 = _mm256_broadcast_ss(&A5);
        ddz_z0_0123 = _mm256_mul_ps(m5, _mm256_sub_ps(pq[5*stride_z_m256], pq[-4*stride_z_m256]));
        ddz_z1_0123 = _mm256_mul_ps(m5, _mm256_sub_ps(pq[6*stride_z_m256], pq[-3*stride_z_m256]));

	__m256 m4 = _mm256_broadcast_ss(&A4);
        ddz_z0_0123 = _mm256_add_ps(ddz_z0_0123, _mm256_mul_ps(m4, _mm256_sub_ps(pq[4*stride_z_m256], pq[-3*stride_z_m256])));
        ddz_z1_0123 = _mm256_add_ps(ddz_z1_0123, _mm256_mul_ps(m4, _mm256_sub_ps(pq[5*stride_z_m256], pq[-2*stride_z_m256])));

	__m256 m3 = _mm256_broadcast_ss(&A3);
        ddz_z0_0123 = _mm256_add_ps(ddz_z0_0123, _mm256_mul_ps(m3, _mm256_sub_ps(pq[3*stride_z_m256], pq[-2*stride_z_m256])));
        ddz_z1_0123 = _mm256_add_ps(ddz_z1_0123, _mm256_mul_ps(m3, _mm256_sub_ps(pq[4*stride_z_m256], pq[-1*stride_z_m256])));

	__m256 m2 = _mm256_broadcast_ss(&A2);
        ddz_z0_0123 = _mm256_add_ps(ddz_z0_0123, _mm256_mul_ps(m2, _mm256_sub_ps(pq[2*stride_z_m256], pq[-1*stride_z_m256])));
        ddz_z1_0123 = _mm256_add_ps(ddz_z1_0123, _mm256_mul_ps(m2, _mm256_sub_ps(pq[3*stride_z_m256], pq[               0])));

	__m256 m1 = _mm256_broadcast_ss(&A1);
        ddz_z0_0123 = _mm256_add_ps(ddz_z0_0123, _mm256_mul_ps(m1, _mm256_sub_ps(pq[  stride_z_m256], pq[            0])));
        ddz_z1_0123 = _mm256_add_ps(ddz_z1_0123, _mm256_mul_ps(m1, _mm256_sub_ps(pq[2*stride_z_m256], pq[stride_z_m256])));
}

inline void AVX_Cmp_DDZ_EE(
        __m256* pq,
        int stride_z_m256,
        __m256& ddz_0123
        )
{
	__m256 E1E1E1E1 = _mm256_broadcast_ss(&E1);
	ddz_0123 = _mm256_mul_ps(E1E1E1E1, _mm256_sub_ps(pq[stride_z_m256], pq[0]));
}

inline void AVX_Cmp_DDZ_Z0(
        __m256* pq,
        int stride_z_m256,
        __m256& ddz_0123
        )
{
	ddz_0123 = _mm256_mul_ps(_mm256_broadcast_ss(&A5), _mm256_add_ps(pq[5*stride_z_m256], pq[4*stride_z_m256]));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A4), _mm256_add_ps(pq[4*stride_z_m256], pq[3*stride_z_m256])));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A3), _mm256_add_ps(pq[3*stride_z_m256], pq[2*stride_z_m256])));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A2), _mm256_add_ps(pq[2*stride_z_m256], pq[stride_z_m256])));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A1), _mm256_sub_ps(pq[stride_z_m256], pq[0])));
}

inline void AVX_Cmp_DDZ_Z1(
        __m256* pq,
        int stride_z_m256,
        __m256& ddz_0123
        )
{
	ddz_0123 = _mm256_mul_ps(_mm256_broadcast_ss(&A5), _mm256_add_ps(pq[5*stride_z_m256], pq[2*stride_z_m256]));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A4), _mm256_add_ps(pq[4*stride_z_m256], pq[stride_z_m256])));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A3), _mm256_add_ps(pq[3*stride_z_m256], pq[0])));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A2), _mm256_sub_ps(pq[2*stride_z_m256], pq[-stride_z_m256])));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A1), _mm256_sub_ps(pq[stride_z_m256], pq[0])));
}

inline void AVX_Cmp_DDZ_Z2(
        __m256* pq,
        int stride_z_m256,
        __m256& ddz_0123
        )
{
	ddz_0123 = _mm256_mul_ps(_mm256_broadcast_ss(&A5), _mm256_add_ps(pq[5*stride_z_m256], pq[0]));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A4), _mm256_add_ps(pq[4*stride_z_m256], pq[-1*stride_z_m256])));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A3), _mm256_sub_ps(pq[3*stride_z_m256], pq[-2*stride_z_m256])));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A2), _mm256_sub_ps(pq[2*stride_z_m256], pq[-stride_z_m256])));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A1), _mm256_sub_ps(pq[stride_z_m256], pq[0])));
}

inline void AVX_Cmp_DDZ_Z3(
        __m256* pq,
        int stride_z_m256,
        __m256& ddz_0123
        )
{
	ddz_0123 = _mm256_mul_ps(_mm256_broadcast_ss(&A5), _mm256_add_ps(pq[5*stride_z_m256], pq[-2*stride_z_m256]));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A4), _mm256_sub_ps(pq[4*stride_z_m256], pq[-3*stride_z_m256])));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A3), _mm256_sub_ps(pq[3*stride_z_m256], pq[-2*stride_z_m256])));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A2), _mm256_sub_ps(pq[2*stride_z_m256], pq[-stride_z_m256])));
	ddz_0123 = _mm256_add_ps(ddz_0123, _mm256_mul_ps(_mm256_broadcast_ss(&A1), _mm256_sub_ps(pq[stride_z_m256], pq[0])));
}

inline void AVX_TTIDenQ_Cmp_DXED_DYED_DZED(
	__m256& ddx_0123,
	__m256& ddy_0123,
	__m256& ddz_0123,
	__m128& Boy,
	__m128& cDip,
	__m128& sDip,
	__m128& cAzm,
	__m128& sAzm,
        __m256& dxed1_0123,
        __m256& dxed2_0123,
        __m256& dyed1_0123,
        __m256& dyed2_0123,
	__m256& dzed1_0123,
	__m256& dzed2_0123
	)
{
        __m128 Boy_01 = _mm_shuffle_ps(Boy, Boy, 0x50);
        __m128 Boy_23 = _mm_shuffle_ps(Boy, Boy, 0xFA);
	__m256 Boy_0123 = conc(Boy_01, Boy_23);

        __m128 cAzm_01 = _mm_shuffle_ps(cAzm, cAzm, 0x50);
        __m128 cAzm_23 = _mm_shuffle_ps(cAzm, cAzm, 0xFA);
	__m256 cAzm_0123 = conc(cAzm_01, cAzm_23);
	
        __m128 sAzm_01 = _mm_shuffle_ps(sAzm, sAzm, 0x50);
        __m128 sAzm_23 = _mm_shuffle_ps(sAzm, sAzm, 0xFA);
	__m256 sAzm_0123 = conc(sAzm_01, sAzm_23);

        __m128 cDip_01 = _mm_shuffle_ps(cDip, cDip, 0x50);
        __m128 cDip_23 = _mm_shuffle_ps(cDip, cDip, 0xFA);
	__m256 cDip_0123 = conc(cDip_01, cDip_23);
	
        __m128 sDip_01 = _mm_shuffle_ps(sDip, sDip, 0x50);
        __m128 sDip_23 = _mm_shuffle_ps(sDip, sDip, 0xFA);
	__m256 sDip_0123 = conc(sDip_01, sDip_23);

	__m256 temp_0123 = _mm256_add_ps(_mm256_mul_ps(cAzm_0123, ddx_0123), _mm256_mul_ps(sAzm_0123, ddy_0123));
	__m256 temp2_0123 = _mm256_mul_ps(Boy_0123, _mm256_sub_ps(_mm256_mul_ps(cAzm_0123, ddy_0123), _mm256_mul_ps(sAzm_0123, ddx_0123)));
	__m256 temp3_0123 = _mm256_mul_ps(Boy_0123, _mm256_sub_ps(_mm256_mul_ps(cDip_0123, temp_0123), _mm256_mul_ps(sDip_0123, ddz_0123)));
	__m256 temp4_0123 = _mm256_mul_ps(Boy_0123, _mm256_add_ps(_mm256_mul_ps(sDip_0123, temp_0123), _mm256_mul_ps(cDip_0123, ddz_0123)));
	
	dxed1_0123 = _mm256_sub_ps(_mm256_mul_ps(_mm256_mul_ps(cDip_0123, cAzm_0123), temp3_0123), _mm256_mul_ps(sAzm_0123, temp2_0123));
	dyed1_0123 = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(cDip_0123, sAzm_0123), temp3_0123), _mm256_mul_ps(cAzm_0123, temp2_0123));
	dzed1_0123 = _mm256_mul_ps(sDip_0123, temp3_0123);

	dxed2_0123 = _mm256_mul_ps(_mm256_mul_ps(sDip_0123, cAzm_0123), temp4_0123);
	dyed2_0123 = _mm256_mul_ps(_mm256_mul_ps(sDip_0123, sAzm_0123), temp4_0123);
	dzed2_0123 = _mm256_mul_ps(cDip_0123, temp4_0123);
}

inline void AVX_TTIDenQ_Cmp_DXED(
	__m256& ddx,
	__m256& ddy,
	__m256& ddz,
	__m128& Boy,
	__m128& cDip,
	__m128& sDip,
	__m128& cAzm,
	__m128& sAzm,
        __m256& dxed1,
        __m256& dxed2
	)
{
        __m128 Boy_01 = _mm_shuffle_ps(Boy, Boy, 0x50);
        __m128 Boy_23 = _mm_shuffle_ps(Boy, Boy, 0xFA);
	__m256 Boy_0123 = conc(Boy_01, Boy_23);

	__m128 cAzm_01 = _mm_shuffle_ps(cAzm, cAzm, 0x50);
	__m128 cAzm_23 = _mm_shuffle_ps(cAzm, cAzm, 0xFA);
	__m256 cAzm_0123 = conc(cAzm_01, cAzm_23);

	__m128 sAzm_01 = _mm_shuffle_ps(sAzm, sAzm, 0x50);
	__m128 sAzm_23 = _mm_shuffle_ps(sAzm, sAzm, 0xFA);
	__m256 sAzm_0123 = conc(sAzm_01, sAzm_23);

	__m256 temp_0123 = _mm256_add_ps(_mm256_mul_ps(cAzm_0123, ddx), _mm256_mul_ps(sAzm_0123, ddy));
	__m256 temp2_0123 = _mm256_mul_ps(Boy_0123, _mm256_sub_ps(_mm256_mul_ps(cAzm_0123, ddy), _mm256_mul_ps(sAzm_0123, ddx)));

	__m128 cDip_01 = _mm_shuffle_ps(cDip, cDip, 0x50);
        __m128 cDip_23 = _mm_shuffle_ps(cDip, cDip, 0xFA);
	__m256 cDip_0123 = conc(cDip_01, cDip_23);

        __m128 sDip_01 = _mm_shuffle_ps(sDip, sDip, 0x50);
        __m128 sDip_23 = _mm_shuffle_ps(sDip, sDip, 0xFA);
	__m256 sDip_0123 = conc(sDip_01, sDip_23);

	__m256 temp3_0123 = _mm256_mul_ps(Boy_0123, _mm256_sub_ps(_mm256_mul_ps(cDip_0123, temp_0123), _mm256_mul_ps(sDip_0123, ddz)));
	__m256 temp4_0123 = _mm256_mul_ps(Boy_0123, _mm256_add_ps(_mm256_mul_ps(sDip_0123, temp_0123), _mm256_mul_ps(cDip_0123, ddz)));
	
	dxed1 = _mm256_sub_ps(_mm256_mul_ps(_mm256_mul_ps(cDip_0123, cAzm_0123), temp3_0123), _mm256_mul_ps(sAzm_0123, temp2_0123));
	dxed2 = _mm256_mul_ps(_mm256_mul_ps(sDip_0123, cAzm_0123), temp4_0123);
}

inline void AVX_TTIDenQ_Cmp_DYED(
	__m256& ddx,
	__m256& ddy,
	__m256& ddz,
	__m128& Boy,
	__m128& cDip,
	__m128& sDip,
	__m128& cAzm,
	__m128& sAzm,
        __m256& dyed1,
        __m256& dyed2
	)
{
        __m128 Boy_01 = _mm_shuffle_ps(Boy, Boy, 0x50);
        __m128 Boy_23 = _mm_shuffle_ps(Boy, Boy, 0xFA);
	__m256 Boy_0123 = conc(Boy_01, Boy_23);

	__m128 cAzm_01 = _mm_shuffle_ps(cAzm, cAzm, 0x50);
        __m128 cAzm_23 = _mm_shuffle_ps(cAzm, cAzm, 0xFA);
	__m256 cAzm_0123 = conc(cAzm_01, cAzm_23);

	__m128 sAzm_01 = _mm_shuffle_ps(sAzm, sAzm, 0x50);
        __m128 sAzm_23 = _mm_shuffle_ps(sAzm, sAzm, 0xFA);
	__m256 sAzm_0123 = conc(sAzm_01, sAzm_23);

	__m256 temp_0123 = _mm256_add_ps(_mm256_mul_ps(cAzm_0123, ddx), _mm256_mul_ps(sAzm_0123, ddy));
	__m256 temp2_0123 = _mm256_mul_ps(Boy_0123, _mm256_sub_ps(_mm256_mul_ps(cAzm_0123, ddy), _mm256_mul_ps(sAzm_0123, ddx)));

	__m128 cDip_01 = _mm_shuffle_ps(cDip, cDip, 0x50);
        __m128 cDip_23 = _mm_shuffle_ps(cDip, cDip, 0xFA);
	__m256 cDip_0123 = conc(cDip_01, cDip_23);

        __m128 sDip_01 = _mm_shuffle_ps(sDip, sDip, 0x50);
        __m128 sDip_23 = _mm_shuffle_ps(sDip, sDip, 0xFA);
	__m256 sDip_0123 = conc(sDip_01, sDip_23);

	__m256 temp3_0123 = _mm256_mul_ps(Boy_0123, _mm256_sub_ps(_mm256_mul_ps(cDip_0123, temp_0123), _mm256_mul_ps(sDip_0123, ddz)));
	__m256 temp4_0123 = _mm256_mul_ps(Boy_0123, _mm256_add_ps(_mm256_mul_ps(sDip_0123, temp_0123), _mm256_mul_ps(cDip_0123, ddz)));
	
	dyed1 = _mm256_add_ps(_mm256_mul_ps(_mm256_mul_ps(cDip_0123, sAzm_0123), temp3_0123), _mm256_mul_ps(cAzm_0123, temp2_0123));
	dyed2 = _mm256_mul_ps(_mm256_mul_ps(sDip_0123, sAzm_0123), temp4_0123);
}

inline void AVX_TTIDenQ_Cmp_Wave_Equation(
	int stage,
	__m256* pq1,
	__m256* rs1,
	__m256* Apq1,
	int* VelAnis1,
	int* DenAng1,
	__m128& spgzyx,
	__m256& V4,
	__m256& V5
	)
{
	__m128 Qatten, C66C44_01, C66C44_23, C44C33_01, C44C33_23, C55_01, C55_23;
	TTIDenQ_Get_EM2(VelAnis1,DenAng1,Qatten,C66C44_01,C66C44_23,C44C33_01,C44C33_23,C55_01,C55_23);
	//
	// r = 2 * p - r + C33 * (1+2*eps) * V5_p + C44 * V4_p + C13pC44 * V4_q
	// s = 2 * q - s + C33 *             V4_q + C44 * V5_q   C13pC44 * V5_p
	//
	// C55 = C13pC44
	// C66 = C33 * (1+2*eps)
	// 
	// r = 2 * p - r + C66 * V5_p + C44 * V4_p + C55 * V4_q
	// s = 2 * q - s + C44 * V5_q + C33 * V4_q + C55 * V5_p
	//
	// Outputs from Get_EM2 :: C66_C44, C44_C33, C55
	//
	// r = 2 * p - r + C55 * V4_q + C44 * V4_p + C66 * V5_p
	// s = 2 * q - s + C55 * V5_q + C44 * V5_p + C33 * V4_p
	//
	//
	// ISO:
	//
	// eps = del = C44C33 = 0
	// C44 = 0
	// C13pC344 = C33
	// V4 = ddx2 + ddy2 
	// V5 = ddz2
	// p(...) == q...) initially
	//
	// r = 2 * p - r + C33 * V5_p + C33 * V4_q
	// s = 2 * q - s + C33 * V4_q + C33 * V5_p
	//
	// r = 2 * p - r + C33 * ( d(buoy * dpdx)dx + d(buoy * dpdy)dy + d(buoy * dpdz)dz )
	// s = r
	// 
	// p and q are equal at time = 0, so r and s will be equal, which means p and q will never be different.
	// only need to propagate p
	// 

	// debugging
	//rs1[0] = rs1[1] = _mm_set1_ps(0.0f);
	//pq1[0] = pq1[1] = _mm_set1_ps(1.0f);
	//V4_01 = V4_23 = _mm_set1_ps(1.8741e-2f);
	//V5_01 = V5_23 = _mm_set1_ps(7.4895e-2f);

	__m256 C66C44 = conc(C66C44_01, C66C44_23);
	__m256 C44C33 = conc(C44C33_01, C44C33_23);
	__m256 C55 = conc(C55_01, C55_23);

	__m256 rs = _mm256_add_ps(_mm256_mul_ps(V5, C66C44), _mm256_mul_ps(V4, C44C33));
	__m256 V4qV5p = _mm256_shuffle_ps(V4, V5, 0x8D);
	V4qV5p = _mm256_shuffle_ps(V4qV5p, V4qV5p, 0xD8);
	rs = _mm256_add_ps(rs, _mm256_mul_ps(V4qV5p, C55));

	if (stage == 0)
	{
		// 2nd order in time
		__m128 spgzyx_01 = _mm_shuffle_ps(spgzyx, spgzyx, 0x50);
		__m128 spgzyx_23 = _mm_shuffle_ps(spgzyx, spgzyx, 0xFA);
		__m256 spgzyx = conc(spgzyx_01, spgzyx_23);

		__m128 Qatten_01 = _mm_shuffle_ps(Qatten, Qatten, 0x50);
		__m128 Qatten_23 = _mm_shuffle_ps(Qatten, Qatten, 0xFA);
		__m256 Qatten_0123 = conc(Qatten_01, Qatten_23);

		rs = _mm256_add_ps(rs, _mm256_add_ps(pq1[0], pq1[0]));
                rs = _mm256_sub_ps(rs, _mm256_mul_ps(spgzyx, rs1[0]));
                rs1[0] = _mm256_mul_ps(_mm256_mul_ps(rs, spgzyx), Qatten_0123);
	}
	else if (stage == 1)
	{
		// 4th order in time - 1st call
		_mm256_stream_ps((float*)Apq1, rs);
	}
	else
	{
		// 4th order in time - 2nd and final call
		__m128 spgzyx_01 = _mm_shuffle_ps(spgzyx, spgzyx, 0x50);
		__m128 spgzyx_23 = _mm_shuffle_ps(spgzyx, spgzyx, 0xFA);
		__m256 spgzyx = conc(spgzyx_01, spgzyx_23);

                __m128 Qatten_01 = _mm_shuffle_ps(Qatten, Qatten, 0x50);
                __m128 Qatten_23 = _mm_shuffle_ps(Qatten, Qatten, 0xFA);
                __m256 Qatten_0123 = conc(Qatten_01, Qatten_23);

		rs = _mm256_mul_ps(rs, _mm256_set1_ps(1.0f/12.0f));
		rs = _mm256_add_ps(rs, Apq1[0]);
                rs = _mm256_add_ps(rs, _mm256_add_ps(pq1[0], pq1[0]));
                rs = _mm256_sub_ps(rs, _mm256_mul_ps(spgzyx, rs1[0]));
		rs1[0] = _mm256_mul_ps(_mm256_mul_ps(rs, spgzyx), Qatten_0123);
	}
}

inline void AVX_VTIDenQ_Cmp_DXED(
	__m256& ddx,
	__m128& Boy,
        __m256& dxed1
	)
{
        __m128 Boy_01 = _mm_shuffle_ps(Boy, Boy, 0x50);
        __m128 Boy_23 = _mm_shuffle_ps(Boy, Boy, 0xFA);
	__m256 Boy_0123 = conc(Boy_01, Boy_23);

	dxed1 = _mm256_mul_ps(Boy_0123, ddx);
}

inline void AVX_VTIDenQ_Cmp_DYED(
	__m256& ddy,
	__m128& Boy,
        __m256& dyed1
	)
{
        __m128 Boy_01 = _mm_shuffle_ps(Boy, Boy, 0x50);
        __m128 Boy_23 = _mm_shuffle_ps(Boy, Boy, 0xFA);
	__m256 Boy_0123 = conc(Boy_01, Boy_23);

	dyed1 = _mm256_mul_ps(Boy_0123, ddy);
}

inline void AVX_VTIDenQ_Cmp_DXED_DYED_DZED(
	__m256& ddx,
	__m256& ddy,
	__m256& ddz,
	__m128& Boy,
        __m256& dxed1,
        __m256& dyed1,
	__m256& dzed2
	)
{
        __m128 Boy_01 = _mm_shuffle_ps(Boy, Boy, 0x50);
        __m128 Boy_23 = _mm_shuffle_ps(Boy, Boy, 0xFA);
	__m256 Boy_0123 = conc(Boy_01, Boy_23);

	dxed1 = _mm256_mul_ps(Boy_0123, ddx);
	dyed1 = _mm256_mul_ps(Boy_0123, ddy);
	dzed2 = _mm256_mul_ps(Boy_0123, ddz);
}

// 
// Fill up the vertical queues prior to looping over Z.
// This routine actually processes one Z more than it needs to.
// This is faster (proven experimentally) than adding extra logic to prevent it.
//
void AVX_VTIDenQ_Process_Patch_Leadin(
	int logLevel,
	__m256* pq,		// pq for stage 0 and 1, Apq for stage 2
	int* DenAng,
	int* VelAnis,
	__m256* stripe,
	__m256* dyed_buf,
	__m256* V4,
	int iX0,
	int iY0,
	int iZ,
	int iXN_halo,
	int iXN,
	int iYN_halo,
	int iYN,
	int stride_y_m256,
	int stride_z_m256,
	int stride_y_em,
	int stride_z_em
	)
{
	int abs_iZ0 = (iZ >= 0) ? iZ : -(iZ+1);
	int abs_iZ1 = ((iZ+1) >= 0) ? (iZ+1) : -(iZ+2);

	__m256* pq0 = pq + (unsigned long)abs_iZ0 * (unsigned long)stride_z_m256 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m256 + (unsigned long)((iX0-4) >> 1);
	int* DenAng0 = DenAng + (unsigned long)abs_iZ0 * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

	int abs_stride_z_m256 = (abs_iZ1 - abs_iZ0) * stride_z_m256;
	int abs_stride_z_em = (abs_iZ1 - abs_iZ0) * stride_z_em;
	if (logLevel >= 6)
	{
		printf("iZ = %d, abs_iZ0 = %d, abs_iZ1 = %d, abs_stride_z_m256 = %d, abs_stride_em = %d\n",iZ,abs_iZ0,abs_iZ1,abs_stride_z_m256,abs_stride_z_em);
		fflush(stdout);
	}

	__m256* V4_0 = V4;
	__m256* V4_1 = V4;

	__m256* dyed_buf_0 = dyed_buf;
	__m256* dyed_buf_1 = dyed_buf + 8 * iXN;

	for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
	{
		pq0 += 2;
		DenAng0 += 8;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m256 ddy;
			AVX_Cmp_DDY(pq0,stride_y_m256,ddy);

			__m128 Boy;
			VTIDenQ_Get_EM1(DenAng0, Boy);

			AVX_VTIDenQ_Cmp_DYED(
					ddy,  
					Boy,
					dyed_buf_0[0]
					);

			AVX_Cmp_DDY(pq0+abs_stride_z_m256,stride_y_m256,ddy);

			VTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy);

			AVX_VTIDenQ_Cmp_DYED(
					ddy,
					Boy,
					dyed_buf_0[1]
					);
			dyed_buf_0 += 2;

			pq0 += 1;
			DenAng0 += 4;
		}

		pq0 += stride_y_m256 - (iXN+2);
		DenAng0 += stride_y_em - (iXN+2)*4;
	}

	for (int iY = (iYN_halo-5);  iY > 4;  --iY)
	{
		__m256* stripe0 = stripe;

		for (int iX = 2;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m256 ddx;
			AVX_Cmp_DDX(pq0,ddx);

			__m128 Boy;
			VTIDenQ_Get_EM1(DenAng0, Boy);

			AVX_VTIDenQ_Cmp_DXED(
					ddx,  
					Boy,
					stripe0[ 0]
					);

			AVX_Cmp_DDX(pq0+abs_stride_z_m256,ddx);

			VTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy);

			AVX_VTIDenQ_Cmp_DXED(
					ddx,  
					Boy,
					stripe0[ 1]
					);
			stripe0 += 2;

			pq0 += 1;
			DenAng0 += 4;
		}

		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m256 ddx, ddy, ddz;
			AVX_Cmp_DDX(pq0,ddx);
			AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
			if (abs_iZ0 == 0)
				AVX_Cmp_DDZ_Z0(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 1)
				AVX_Cmp_DDZ_Z1(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 2)
				AVX_Cmp_DDZ_Z2(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 3)
				AVX_Cmp_DDZ_Z3(pq0,stride_z_m256,ddz);
			else
				AVX_Cmp_DDZ(pq0,stride_z_m256,ddz);

			__m128 Boy;
			VTIDenQ_Get_EM1(DenAng0, Boy);

			AVX_VTIDenQ_Cmp_DXED_DYED_DZED(
					ddx, ddy, ddz, 
					Boy,
					stripe0[ 0], 
					dyed_buf_0[0],
					V4_0[7]
					);

			AVX_Cmp_DDX(pq0+abs_stride_z_m256,ddx);
			AVX_Cmp_DDY(pq0+abs_stride_z_m256,stride_y_m256,ddy);
			if (abs_iZ1 == 0)
				AVX_Cmp_DDZ_Z0(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 1)
				AVX_Cmp_DDZ_Z1(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 2)
				AVX_Cmp_DDZ_Z2(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 3)
				AVX_Cmp_DDZ_Z3(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else
				AVX_Cmp_DDZ(pq0+abs_stride_z_m256,stride_z_m256,ddz);

			VTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy);

			AVX_VTIDenQ_Cmp_DXED_DYED_DZED(
					ddx, ddy, ddz, 
					Boy,
					stripe0[ 1],
					dyed_buf_0[1],
					V4_0[6]
					);
			stripe0 += 2;	
			dyed_buf_0 += 2;
			V4_0 += 17;

			pq0 += 1;
			DenAng0 += 4;
		}
		V4_0 -= 17 * iXN;

		for (int iX = 1;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m256 ddx;
			AVX_Cmp_DDX(pq0,ddx);

			__m128 Boy;
			VTIDenQ_Get_EM1(DenAng0, Boy);

			AVX_VTIDenQ_Cmp_DXED(
					ddx,   
					Boy,
					stripe0[ 0]
					);

			AVX_Cmp_DDX(pq0+abs_stride_z_m256,ddx);

			VTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy);

			AVX_VTIDenQ_Cmp_DXED(
					ddx,  
					Boy,
					stripe0[ 1]
					);
			stripe0 += 2;

			pq0 += 1;
			DenAng0 += 4;
		}

		pq0 += stride_y_m256 - iXN_halo;
		DenAng0 += stride_y_em - iXN_halo*4;

		stripe0 = stripe + 4;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			if (iY <= (iYN_halo-9))
			{
				// dydyed1, dydyed2
				__m256 dydyed1;
				AVX_Cmp_DDY(dyed_buf_1,2*iXN,dydyed1);
				V4_1[1] = _mm256_add_ps(V4_1[1], dydyed1);

				AVX_Cmp_DDY(dyed_buf_1+1,2*iXN,dydyed1);
				V4_1[0] = _mm256_add_ps(V4_1[0], dydyed1);

				dyed_buf_1 += 2;
				V4_1 += 17;
			}

			// dxdxed1, dxdxed2
			__m256 dxdxed1;
			AVX_Cmp_DDX_1x_stride4(stripe0,dxdxed1);
			V4_0[1] = dxdxed1;

			AVX_Cmp_DDX_1x_stride4(stripe0+1,dxdxed1);
			V4_0[0] = dxdxed1;

			V4_0 += 17;

			stripe0 += 2;	
		}
	}

	for (int iY = 4;  iY > 0;  --iY)
	{
		pq0 += 2;
		DenAng0 += 8;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m256 ddy;
			AVX_Cmp_DDY(pq0,stride_y_m256,ddy);

			__m128 Boy;
			VTIDenQ_Get_EM1(DenAng0, Boy);

			AVX_VTIDenQ_Cmp_DYED(
					ddy,  
					Boy,
					dyed_buf_0[0]
					);

			AVX_Cmp_DDY(pq0+abs_stride_z_m256,stride_y_m256,ddy);

			VTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy);

			AVX_VTIDenQ_Cmp_DYED(
					ddy,  
					Boy,
					dyed_buf_0[1]
					);
			dyed_buf_0 += 2;

			pq0 += 1;
			DenAng0 += 4;

			// dydyed1, dydyed2
			__m256 dydyed1;
			AVX_Cmp_DDY(dyed_buf_1,2*iXN,dydyed1);
			V4_1[1] = _mm256_add_ps(V4_1[1], dydyed1);

			AVX_Cmp_DDY(dyed_buf_1+1,2*iXN,dydyed1);
			V4_1[0] = _mm256_add_ps(V4_1[0], dydyed1);

			dyed_buf_1 += 2;
			V4_1 += 17;
		}

		pq0 += stride_y_m256 - (iXN+2);
		DenAng0 += stride_y_em - (iXN+2)*4;
	}
}

//
// Valid from [1,dimz-1>
//
void AVX_VTIDenQ_Process_Patch(
	int logLevel,
	int stage,
	__m256* pq,
	__m256* rs,
	__m256* Apq,			// pass pq ptr if 2nd order in time
	int* DenAng,
	int* VelAnis,
	__m128* _mm_spgx,
	__m128* _mm_spgy,
	__m128 _mm_spgz0,
	__m128 _mm_spgz1,
	__m256* stripe,
	__m256* dyed_buf,
	__m256* V4,
	int iX0,
	int iY0,
	int iZ,
	int dimz,
	int iXN_halo,
	int iXN,
	int iYN_halo,
	int iYN,
	int stride_y_m256,
	int stride_z_m256,
	int stride_y_em,
	int stride_z_em
#ifdef TMJ_TIMING
	,unsigned long& thrcnt2
#endif
	)
{
	int dz_edge_diff = dimz - iZ - 1;  // dimz-1 -> 0, dimz-2 -> 1 etc.

	__m256* pq1 = pq + (unsigned long)iZ * (unsigned long)stride_z_m256 + (unsigned long)iY0 * (unsigned long)stride_y_m256 + (unsigned long)(iX0 >> 1);
	__m256* rs1 = rs + (unsigned long)iZ * (unsigned long)stride_z_m256 + (unsigned long)iY0 * (unsigned long)stride_y_m256 + (unsigned long)(iX0 >> 1);
	__m256* Apq1 = Apq + (unsigned long)iZ * (unsigned long)stride_z_m256 + (unsigned long)iY0 * (unsigned long)stride_y_m256 + (unsigned long)(iX0 >> 1);
	int* VelAnis1 = VelAnis + (unsigned long)iZ * (unsigned long)stride_z_em + (unsigned long)iY0 * (unsigned long)stride_y_em + (unsigned long)(iX0*2);
	int* DenAng1 = DenAng + (unsigned long)iZ * (unsigned long)stride_z_em + (unsigned long)iY0 * (unsigned long)stride_y_em + (unsigned long)(iX0*2);

	if (logLevel >= 6)
	{
		printf("Process_Patch :: iZ = %d\n",iZ);
		fflush(stdout);
	}

	// process next tile ("patch")
	if (dz_edge_diff > 5)
	{
	 	__m256* pq0 = ((stage == 0 || stage == 1) ? pq : Apq) + (unsigned long)(iZ+4) * (unsigned long)stride_z_m256 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m256 + (unsigned long)((iX0-4) >> 1);
		int* DenAng0 = DenAng + (unsigned long)(iZ+4) * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

		__m256* dyed_buf_0 = dyed_buf;
		__m256* dyed_buf_1 = dyed_buf + 8 * iXN;

		__m256* V4_0 = V4;
		__m256* V4_1 = V4;

		for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
		{
			pq0 += 2;
			DenAng0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m256), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m256+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m256), _MM_HINT_T0);

				__m256 ddy;
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				AVX_VTIDenQ_Cmp_DYED(
						ddy,  
						Boy,
						dyed_buf_0[0]
						);

				AVX_Cmp_DDY(pq0+stride_z_m256,stride_y_m256,ddy);

				VTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy);

				AVX_VTIDenQ_Cmp_DYED(
						ddy,  
						Boy,
                                                dyed_buf_0[1]
						);

				dyed_buf_0 += 2;
				pq0 += 1;
				DenAng0 += 4;
			}
			pq0 += stride_y_m256 - (iXN+2);
			DenAng0 += stride_y_em - (iXN+2)*4;
		}

		for (int iY = (iYN_halo-5);  iY > 4;  --iY)
		{
			__m256* stripe0 = stripe;

			for (int iX = 2;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m256), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m256+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m256), _MM_HINT_T0);

				__m256 ddx;
				AVX_Cmp_DDX(pq0,ddx);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				AVX_VTIDenQ_Cmp_DXED(
						ddx,  
						Boy,
						stripe0[ 0]
						);

				AVX_Cmp_DDX(pq0+stride_z_m256,ddx);

				VTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy);

				AVX_VTIDenQ_Cmp_DXED(
						ddx,  
						Boy,
						stripe0[ 1]
						);
				stripe0 += 2;

				pq0 += 1;
				DenAng0 += 4;
			}

			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m256), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m256+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m256), _MM_HINT_T0);

				__m256 ddx, ddy, ddz;
				AVX_Cmp_DDX(pq0,ddx);
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
				if (dz_edge_diff > 9)
				{
					AVX_Cmp_DDZ(pq0,stride_z_m256,ddz);
				}
				else
				{
					AVX_Cmp_DDZ_EE(pq0,stride_z_m256,ddz);
				}

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				AVX_VTIDenQ_Cmp_DXED_DYED_DZED(
						ddx, ddy, ddz, 
						Boy,
						stripe0[ 0],
						dyed_buf_0[0],
						V4_0[7]
						);

				AVX_Cmp_DDX(pq0+stride_z_m256,ddx);
				AVX_Cmp_DDY(pq0+stride_z_m256,stride_y_m256,ddy);
				if (dz_edge_diff > 10)
				{
					AVX_Cmp_DDZ(pq0+stride_z_m256,stride_z_m256,ddz);
				}
				else
				{
					AVX_Cmp_DDZ_EE(pq0+stride_z_m256,stride_z_m256,ddz);
				}

				VTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy);

				AVX_VTIDenQ_Cmp_DXED_DYED_DZED(
						ddx, ddy, ddz, 
						Boy,
						stripe0[ 1],
                                                dyed_buf_0[1],
						V4_0[6]
						);

				stripe0 += 2;
				V4_0 += 17;

				dyed_buf_0 += 2;
				pq0 += 1;
				DenAng0 += 4;
			}
			V4_0 -= 17 * iXN;

			for (int iX = 1;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m256), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m256+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m256), _MM_HINT_T0);

				__m256 ddx;
				AVX_Cmp_DDX(pq0,ddx);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				AVX_VTIDenQ_Cmp_DXED(
						ddx,  
						Boy,
						stripe0[ 0]
						);

				AVX_Cmp_DDX(pq0+stride_z_m256,ddx);

				VTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy);

				AVX_VTIDenQ_Cmp_DXED(
						ddx,  
						Boy,
						stripe0[ 1]
						);
				stripe0 += 2;

				pq0 += 1;
				DenAng0 += 4;
			}

			pq0 += stride_y_m256 - iXN_halo;
			DenAng0 += stride_y_em - iXN_halo*4;

			stripe0 = stripe + 4;
			__m128 _mm_spgzy0 = _mm_mul_ps(_mm_spgz0, _mm_spgy[iY]);
			__m128 _mm_spgzy1 = _mm_mul_ps(_mm_spgz1, _mm_spgy[iY]);
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng1), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng1+stride_z_em), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelAnis1), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelAnis1+stride_z_em), _MM_HINT_NTA);
				if (stage == 0)
				{
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m256), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m256), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m256), _MM_HINT_T0);
				}

				if (iY <= (iYN_halo-9))
				{	
					// dydyed1, dydyed2
					__m256 dydyed1;
					AVX_Cmp_DDY(dyed_buf_1,2*iXN,dydyed1);
					V4_1[1] = _mm256_add_ps(V4_1[1], dydyed1);

					AVX_Cmp_DDY(dyed_buf_1+1,2*iXN,dydyed1);
					V4_1[0] = _mm256_add_ps(V4_1[0], dydyed1);

					dyed_buf_1 += 2;
					V4_1 += 17;
				}

				// dxdxed1, dxdxed2
				__m256 dxdxed1;
				AVX_Cmp_DDX_1x_stride4(stripe0,dxdxed1);
				V4_0[1] = dxdxed1;

				AVX_Cmp_DDX_1x_stride4(stripe0+1,dxdxed1);
				V4_0[0] = dxdxed1;

				// dzdzed1, dzdzed2
				__m256 V5 = V4_0[5];

				__m256 dzdzed2;
				AVX_Cmp_DDZ(V4_0+12,-1,dzdzed2);
				__m256 V4 = dzdzed2;

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				AVX_TTIDenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelAnis1,DenAng1,spgzyx0,V4,V5
						);

				V5 = V4_0[ 4];

				AVX_Cmp_DDZ(V4_0+11,-1,dzdzed2);
				V4 = dzdzed2;

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				AVX_TTIDenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m256,rs1+stride_z_m256,Apq1+stride_z_m256,VelAnis1+stride_z_em,DenAng1+stride_z_em,spgzyx1,V4,V5
						);

				V4_0 += 17;

				pq1 += 1;
				rs1 += 1;
				Apq1 += 1;
				VelAnis1 += 4;
				DenAng1 += 4;
				stripe0 += 2;	
#ifdef TMJ_TIMING
				thrcnt2+=2;
#endif
			}
			pq1 += stride_y_m256 - iXN;
			rs1 += stride_y_m256 - iXN;
			Apq1 += stride_y_m256 - iXN;
			VelAnis1 += stride_y_em - iXN*4;
			DenAng1 += stride_y_em - iXN*4;
		}

		for (int iY = 4;  iY > 0;  --iY)
		{
			pq0 += 2;
			DenAng0 += 8;	
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m256), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m256+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m256), _MM_HINT_T0);

				__m256 ddy;
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				AVX_VTIDenQ_Cmp_DYED(
						ddy,  
						Boy,
						dyed_buf_0[0]
						);

				AVX_Cmp_DDY(pq0+stride_z_m256,stride_y_m256,ddy);

				VTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy);

				AVX_VTIDenQ_Cmp_DYED(
						ddy,  
						Boy,
                                                dyed_buf_0[1]
						);

				dyed_buf_0 += 2;
				pq0 += 1;
				DenAng0 += 4;

				// dydyed1, dydyed2
				__m256 dydyed1;
				AVX_Cmp_DDY(dyed_buf_1,2*iXN,dydyed1);
				V4_1[1] = _mm256_add_ps(V4_1[1], dydyed1);

				AVX_Cmp_DDY(dyed_buf_1+1,2*iXN,dydyed1);
				V4_1[0] = _mm256_add_ps(V4_1[0], dydyed1);

				dyed_buf_1 += 2;
				V4_1 += 17;
			}
	
			pq0 += stride_y_m256 - (iXN+2);
			DenAng0 += stride_y_em - (iXN+2)*4;
		}
	}
	else if (dz_edge_diff == 5)
	{
	 	__m256* pq0 = ((stage == 0 || stage == 1) ? pq : Apq) + (unsigned long)(iZ+4) * (unsigned long)stride_z_m256 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m256 + (unsigned long)((iX0-4) >> 1);
		int* DenAng0 = DenAng + (unsigned long)(iZ+4) * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

		__m256* dyed_buf_0 = dyed_buf;
		__m256* dyed_buf_1 = dyed_buf + 8 * iXN;

		__m256* V4_0 = V4;
		__m256* V4_1 = V4;

		for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
		{
			pq0 += 1;
			DenAng0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m256 ddy;
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				AVX_VTIDenQ_Cmp_DYED(
						ddy,  
						Boy,
						dyed_buf_0[0]
						);
				dyed_buf_0 += 2;

				pq0 += 1;
				DenAng0 += 4;
			}

			pq0 += stride_y_m256 - (iXN+2);
			DenAng0 += stride_y_em - (iXN+2)*4;
		}

		for (int iY = (iYN_halo-5);  iY > 4;  --iY)
		{
			__m256* stripe0 = stripe;

			for (int iX = 2;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m256 ddx;
				AVX_Cmp_DDX(pq0,ddx);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				AVX_VTIDenQ_Cmp_DXED(
						ddx,  
						Boy,
						stripe0[ 0]
						);
				stripe0 += 2;

				pq0 += 1;
				DenAng0 += 4;
			}

			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m256 ddx, ddy, ddz;
				AVX_Cmp_DDX(pq0,ddx);
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
				AVX_Cmp_DDZ_EE(pq0,stride_z_m256,ddz);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				AVX_VTIDenQ_Cmp_DXED_DYED_DZED(
						ddx, ddy, ddz, 
						Boy,
						stripe0[ 0],
						dyed_buf_0[0],
						V4_0[7]
						);
				stripe0 += 2;
				dyed_buf_0 += 2;
				V4_0 += 17;

				pq0 += 1;
				DenAng0 += 4;
			}
			V4_0 -= 17 * iXN;

			for (int iX = 1;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m256 ddx;
				AVX_Cmp_DDX(pq0,ddx);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				AVX_VTIDenQ_Cmp_DXED(
						ddx,  
						Boy,
						stripe0[ 0]
						);
				stripe0 += 2;

				pq0 += 1;
				DenAng0 += 4;
			}

			pq0 += stride_y_m256 - iXN_halo;
			DenAng0 += stride_y_em - iXN_halo*4;

			stripe0 = stripe + 4;
			__m128 _mm_spgzy0 = _mm_mul_ps(_mm_spgz0, _mm_spgy[iY]);
			__m128 _mm_spgzy1 = _mm_mul_ps(_mm_spgz1, _mm_spgy[iY]);
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng1), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng1+stride_z_em), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelAnis1), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelAnis1+stride_z_em), _MM_HINT_NTA);
				if (stage == 0)
				{
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m256), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m256), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m256), _MM_HINT_T0);
				}

				if (iY <= (iYN_halo-9))
				{
					// dydyed1, dydyed2
					__m256 dydyed1;
					AVX_Cmp_DDY(dyed_buf_1,2*iXN,dydyed1);
					V4_1[1] = _mm256_add_ps(V4_1[1], dydyed1);

					dyed_buf_1 += 2;
					V4_1 += 17;
				}

				// dxdxed1, dxdxed2
				__m256 dxdxed1;
				AVX_Cmp_DDX_1x_stride4(stripe0,dxdxed1);
				V4_0[1] = dxdxed1;

				// dzdzed1, dzdzed2
				__m256 V5 = V4_0[5];

				__m256 dzdzed2;
				AVX_Cmp_DDZ(V4_0+12,-1,dzdzed2);
				__m256 V4 = dzdzed2;

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				AVX_TTIDenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelAnis1,DenAng1,spgzyx0,V4,V5
						);

				V5 = V4_0[ 4];

				AVX_Cmp_DDZ_EE(V4_0+11,-1,dzdzed2);
				V4 = dzdzed2;

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				AVX_TTIDenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m256,rs1+stride_z_m256,Apq1+stride_z_m256,VelAnis1+stride_z_em,DenAng1+stride_z_em,spgzyx1,V4,V5
						);

				V4_0 += 17;

				pq1 += 1;
				rs1 += 1;
				Apq1 += 1;
				VelAnis1 += 4;
				DenAng1 += 4;
				stripe0 += 2;	
#ifdef TMJ_TIMING
				thrcnt2+=2;
#endif
			}
			pq1 += stride_y_m256 - iXN;
			rs1 += stride_y_m256 - iXN;
			Apq1 += stride_y_m256 - iXN;
			VelAnis1 += stride_y_em - iXN*4;
			DenAng1 += stride_y_em - iXN*4;
		}

		for (int iY = 4;  iY > 0;  --iY)
		{
			pq0 += 2;
			DenAng0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m256 ddy;
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				AVX_VTIDenQ_Cmp_DYED(
						ddy,  
						Boy,
						dyed_buf_0[0]
						);
				dyed_buf_0 += 2;

				pq0 += 1;
				DenAng0 += 4;

				// dydyed1, dydyed2
				__m256 dydyed1;
				AVX_Cmp_DDY(dyed_buf_1,2*iXN,dydyed1);
				V4_1[1] = _mm256_add_ps(V4_1[1], dydyed1);

				dyed_buf_1 += 2;
				V4_1 += 17;
			}

			pq0 += stride_y_m256 - (iXN+2);
			DenAng0 += stride_y_em - (iXN+2)*4;
		}
	}
	else // no more tiles to process, finish rest of outputs from queue
	{
		__m256* V4_0 = V4;
		for (int iY = iYN;  iY > 0;  --iY)
		{
			__m128 _mm_spgzy0 = _mm_mul_ps(_mm_spgz0, _mm_spgy[iY+4]);
			__m128 _mm_spgzy1 = _mm_mul_ps(_mm_spgz1, _mm_spgy[iY+4]);
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng1), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng1+stride_z_em), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelAnis1), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelAnis1+stride_z_em), _MM_HINT_NTA);
				if (stage == 0)
				{
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m256), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m256), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m256), _MM_HINT_T0);
				}

				// dzdzed1, dzdzed2
				// use shorter two point stencil
				__m256 V5 = V4_0[5];

				__m256 dzdzed2;
				AVX_Cmp_DDZ_EE(V4_0+12,-1,dzdzed2);
				__m256 V4 = dzdzed2;

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				AVX_TTIDenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelAnis1,DenAng1,spgzyx0,V4,V5
						);

				V5 = V4_0[ 4];

				AVX_Cmp_DDZ_EE(V4_0+11,-1,dzdzed2);
				V4 = dzdzed2;

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				AVX_TTIDenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m256,rs1+stride_z_m256,Apq1+stride_z_m256,VelAnis1+stride_z_em,DenAng1+stride_z_em,spgzyx1,V4,V5
						);

				V4_0 += 17;

				pq1 += 1;
				rs1 += 1;
				Apq1 += 1;
				VelAnis1 += 4;
				DenAng1 += 4;
#ifdef TMJ_TIMING
				thrcnt2 += 2;
#endif
			}
			pq1 += stride_y_m256 - iXN;
			rs1 += stride_y_m256 - iXN;
			Apq1 += stride_y_m256 - iXN;
			VelAnis1 += stride_y_em - iXN*4;
			DenAng1 += stride_y_em - iXN*4;
		}
	}
}

// 
// Fill up the vertical queues prior to looping over Z.
// This routine actually processes one Z more than it needs to.
// This is faster (proven experimentally) than adding extra logic to prevent it.
//
void AVX_TTIDenQ_Process_Patch_Leadin(
	int logLevel,
	__m256* pq,		// pq for stage 0 and 1, Apq for stage 2
	int* DenAng,
	int* VelAnis,
	__m256* stripe,
	__m256* dyed_buf,
	__m256* V4,
	int iX0,
	int iY0,
	int iZ,
	int iXN_halo,
	int iXN,
	int iYN_halo,
	int iYN,
	int stride_y_m256,
	int stride_z_m256,
	int stride_y_em,
	int stride_z_em
	)
{
	int abs_iZ0 = (iZ >= 0) ? iZ : -(iZ+1);
	int abs_iZ1 = ((iZ+1) >= 0) ? (iZ+1) : -(iZ+2);

	__m256* pq0 = pq + (unsigned long)abs_iZ0 * (unsigned long)stride_z_m256 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m256 + (unsigned long)((iX0-4) >> 1);
	int* DenAng0 = DenAng + (unsigned long)abs_iZ0 * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

	int abs_stride_z_m256 = (abs_iZ1 - abs_iZ0) * stride_z_m256;
	int abs_stride_z_em = (abs_iZ1 - abs_iZ0) * stride_z_em;
	if (logLevel >= 6)
	{
		printf("iZ = %d, abs_iZ0 = %d, abs_iZ1 = %d, abs_stride_z_m256 = %d, abs_stride_em = %d\n",iZ,abs_iZ0,abs_iZ1,abs_stride_z_m256,abs_stride_z_em);
		fflush(stdout);
	}

	__m256* V4_0 = V4;
	__m256* V4_1 = V4;

	__m256* dyed_buf_0 = dyed_buf;
	__m256* dyed_buf_1 = dyed_buf + 16 * iXN;

	for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
	{
		pq0 += 2;
		DenAng0 += 8;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m256 ddx, ddy, ddz;
			AVX_Cmp_DDX(pq0,ddx);
			AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
			if (abs_iZ0 == 0)
				AVX_Cmp_DDZ_Z0(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 1)
				AVX_Cmp_DDZ_Z1(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 2)
				AVX_Cmp_DDZ_Z2(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 3)
				AVX_Cmp_DDZ_Z3(pq0,stride_z_m256,ddz);
			else
				AVX_Cmp_DDZ(pq0,stride_z_m256,ddz);

			__m128 Boy, cDip, sDip, cAzm, sAzm;
			TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

			AVX_TTIDenQ_Cmp_DYED(
					ddx, ddy, ddz, 
					Boy, cDip, sDip, cAzm, sAzm,
					dyed_buf_0[0], dyed_buf_0[1]
					);

			AVX_Cmp_DDX(pq0+abs_stride_z_m256,ddx);
			AVX_Cmp_DDY(pq0+abs_stride_z_m256,stride_y_m256,ddy);
			if (abs_iZ1 == 0)
				AVX_Cmp_DDZ_Z0(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 1)
				AVX_Cmp_DDZ_Z1(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 2)
				AVX_Cmp_DDZ_Z2(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 3)
				AVX_Cmp_DDZ_Z3(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else
				AVX_Cmp_DDZ(pq0+abs_stride_z_m256,stride_z_m256,ddz);

			TTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

			AVX_TTIDenQ_Cmp_DYED(
					ddx, ddy, ddz, 
					Boy, cDip, sDip, cAzm, sAzm,
					dyed_buf_0[2], dyed_buf_0[3]
					);
			dyed_buf_0 += 4;

			pq0 += 1;
			DenAng0 += 4;
		}

		pq0 += stride_y_m256 - (iXN+2);
		DenAng0 += stride_y_em - (iXN+2)*4;
	}

	for (int iY = (iYN_halo-5);  iY > 4;  --iY)
	{
		__m256* stripe0 = stripe;

		for (int iX = 2;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m256 ddx, ddy, ddz;
			AVX_Cmp_DDX(pq0,ddx);
			AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
			if (abs_iZ0 == 0)
				AVX_Cmp_DDZ_Z0(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 1)
				AVX_Cmp_DDZ_Z1(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 2)
				AVX_Cmp_DDZ_Z2(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 3)
				AVX_Cmp_DDZ_Z3(pq0,stride_z_m256,ddz);
			else
				AVX_Cmp_DDZ(pq0,stride_z_m256,ddz);

			__m128 Boy, cDip, sDip, cAzm, sAzm;
			TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

			AVX_TTIDenQ_Cmp_DXED(
					ddx, ddy, ddz, 
					Boy, cDip, sDip, cAzm, sAzm,
					stripe0[ 0], stripe0[ 1]
					);

			AVX_Cmp_DDX(pq0+abs_stride_z_m256,ddx);
			AVX_Cmp_DDY(pq0+abs_stride_z_m256,stride_y_m256,ddy);
			if (abs_iZ1 == 0)
				AVX_Cmp_DDZ_Z0(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 1)
				AVX_Cmp_DDZ_Z1(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 2)
				AVX_Cmp_DDZ_Z2(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 3)
				AVX_Cmp_DDZ_Z3(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else
				AVX_Cmp_DDZ(pq0+abs_stride_z_m256,stride_z_m256,ddz);

			TTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

			AVX_TTIDenQ_Cmp_DXED(
					ddx, ddy, ddz, 
					Boy, cDip, sDip, cAzm, sAzm,
					stripe0[ 2], stripe0[ 3]
					);
			stripe0 += 4;

			pq0 += 1;
			DenAng0 += 4;
		}

		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m256 ddx, ddy, ddz;
			AVX_Cmp_DDX(pq0,ddx);
			AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
			if (abs_iZ0 == 0)
				AVX_Cmp_DDZ_Z0(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 1)
				AVX_Cmp_DDZ_Z1(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 2)
				AVX_Cmp_DDZ_Z2(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 3)
				AVX_Cmp_DDZ_Z3(pq0,stride_z_m256,ddz);
			else
				AVX_Cmp_DDZ(pq0,stride_z_m256,ddz);

			__m128 Boy, cDip, sDip, cAzm, sAzm;
			TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

			AVX_TTIDenQ_Cmp_DXED_DYED_DZED(
					ddx, ddy, ddz, 
					Boy, cDip, sDip, cAzm, sAzm,
					stripe0[ 0], stripe0[ 1],
					dyed_buf_0[0], dyed_buf_0[1],
					V4_0[14], V4_0[15]
					);

			AVX_Cmp_DDX(pq0+abs_stride_z_m256,ddx);
			AVX_Cmp_DDY(pq0+abs_stride_z_m256,stride_y_m256,ddy);
			if (abs_iZ1 == 0)
				AVX_Cmp_DDZ_Z0(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 1)
				AVX_Cmp_DDZ_Z1(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 2)
				AVX_Cmp_DDZ_Z2(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 3)
				AVX_Cmp_DDZ_Z3(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else
				AVX_Cmp_DDZ(pq0+abs_stride_z_m256,stride_z_m256,ddz);

			TTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

			AVX_TTIDenQ_Cmp_DXED_DYED_DZED(
					ddx, ddy, ddz, 
					Boy, cDip, sDip, cAzm, sAzm,
					stripe0[ 2], stripe0[ 3],
					dyed_buf_0[2], dyed_buf_0[3],
					V4_0[12], V4_0[13]
					);
			stripe0 += 4;	
			dyed_buf_0 += 4;
			V4_0 += 34;

			pq0 += 1;
			DenAng0 += 4;
		}
		V4_0 -= 34 * iXN;

		for (int iX = 1;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m256 ddx, ddy, ddz;
			AVX_Cmp_DDX(pq0,ddx);
			AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
			if (abs_iZ0 == 0)
				AVX_Cmp_DDZ_Z0(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 1)
				AVX_Cmp_DDZ_Z1(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 2)
				AVX_Cmp_DDZ_Z2(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 3)
				AVX_Cmp_DDZ_Z3(pq0,stride_z_m256,ddz);
			else
				AVX_Cmp_DDZ(pq0,stride_z_m256,ddz);

			__m128 Boy, cDip, sDip, cAzm, sAzm;
			TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

			AVX_TTIDenQ_Cmp_DXED(
					ddx, ddy, ddz, 
					Boy, cDip, sDip, cAzm, sAzm,
					stripe0[ 0], stripe0[ 1]
					);

			AVX_Cmp_DDX(pq0+abs_stride_z_m256,ddx);
			AVX_Cmp_DDY(pq0+abs_stride_z_m256,stride_y_m256,ddy);
			if (abs_iZ1 == 0)
				AVX_Cmp_DDZ_Z0(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 1)
				AVX_Cmp_DDZ_Z1(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 2)
				AVX_Cmp_DDZ_Z2(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 3)
				AVX_Cmp_DDZ_Z3(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else
				AVX_Cmp_DDZ(pq0+abs_stride_z_m256,stride_z_m256,ddz);

			TTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

			AVX_TTIDenQ_Cmp_DXED(
					ddx, ddy, ddz, 
					Boy, cDip, sDip, cAzm, sAzm,
					stripe0[ 2], stripe0[ 3]
					);
			stripe0 += 4;

			pq0 += 1;
			DenAng0 += 4;
		}

		pq0 += stride_y_m256 - iXN_halo;
		DenAng0 += stride_y_em - iXN_halo*4;

		stripe0 = stripe + 8;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			if (iY <= (iYN_halo-9))
			{
				// dydyed1, dydyed2
				__m256 dydyed1;
				AVX_Cmp_DDY(dyed_buf_1,4*iXN,dydyed1);
				V4_1[2] = _mm256_add_ps(V4_1[2], dydyed1);

				__m256 dydyed2;
				AVX_Cmp_DDY(dyed_buf_1+1,4*iXN,dydyed2);
				V4_1[3] = _mm256_add_ps(V4_1[3], dydyed2);

				AVX_Cmp_DDY(dyed_buf_1+2,4*iXN,dydyed1);
				V4_1[0] = _mm256_add_ps(V4_1[0], dydyed1);

				AVX_Cmp_DDY(dyed_buf_1+3,4*iXN,dydyed2);
				V4_1[1] = _mm256_add_ps(V4_1[1], dydyed2);
				dyed_buf_1 += 4;
				V4_1 += 34;
			}

			// dxdxed1, dxdxed2
			__m256 dxdxed1;
			AVX_Cmp_DDX_1x_stride8(stripe0,dxdxed1);
			V4_0[2] = dxdxed1;

			__m256 dxdxed2;
			AVX_Cmp_DDX_1x_stride8(stripe0+1,dxdxed2);
			V4_0[3] = dxdxed2;

			AVX_Cmp_DDX_1x_stride8(stripe0+2,dxdxed1);
			V4_0[0] = dxdxed1;

			AVX_Cmp_DDX_1x_stride8(stripe0+3,dxdxed2);
			V4_0[1] = dxdxed2;

			V4_0 += 34;

			stripe0 += 4;	
		}
	}

	for (int iY = 4;  iY > 0;  --iY)
	{
		pq0 += 2;
		DenAng0 += 8;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m256 ddx, ddy, ddz;
			AVX_Cmp_DDX(pq0,ddx);
			AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
			if (abs_iZ0 == 0)
				AVX_Cmp_DDZ_Z0(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 1)
				AVX_Cmp_DDZ_Z1(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 2)
				AVX_Cmp_DDZ_Z2(pq0,stride_z_m256,ddz);
			else if (abs_iZ0 == 3)
				AVX_Cmp_DDZ_Z3(pq0,stride_z_m256,ddz);
			else
				AVX_Cmp_DDZ(pq0,stride_z_m256,ddz);

			__m128 Boy, cDip, sDip, cAzm, sAzm;
			TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

			AVX_TTIDenQ_Cmp_DYED(
					ddx, ddy, ddz, 
					Boy, cDip, sDip, cAzm, sAzm,
					dyed_buf_0[0], dyed_buf_0[1]
					);

			AVX_Cmp_DDX(pq0+abs_stride_z_m256,ddx);
			AVX_Cmp_DDY(pq0+abs_stride_z_m256,stride_y_m256,ddy);
			if (abs_iZ1 == 0)
				AVX_Cmp_DDZ_Z0(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 1)
				AVX_Cmp_DDZ_Z1(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 2)
				AVX_Cmp_DDZ_Z2(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else if (abs_iZ1 == 3)
				AVX_Cmp_DDZ_Z3(pq0+abs_stride_z_m256,stride_z_m256,ddz);
			else
				AVX_Cmp_DDZ(pq0+abs_stride_z_m256,stride_z_m256,ddz);

			TTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

			AVX_TTIDenQ_Cmp_DYED(
					ddx, ddy, ddz, 
					Boy, cDip, sDip, cAzm, sAzm,
					dyed_buf_0[2], dyed_buf_0[3]
					);
			dyed_buf_0 += 4;

			pq0 += 1;
			DenAng0 += 4;

			// dydyed1, dydyed2
			__m256 dydyed1;
			AVX_Cmp_DDY(dyed_buf_1,4*iXN,dydyed1);
			V4_1[2] = _mm256_add_ps(V4_1[2], dydyed1);

			__m256 dydyed2;
			AVX_Cmp_DDY(dyed_buf_1+1,4*iXN,dydyed2);
			V4_1[3] = _mm256_add_ps(V4_1[3], dydyed2);

			AVX_Cmp_DDY(dyed_buf_1+2,4*iXN,dydyed1);
			V4_1[0] = _mm256_add_ps(V4_1[0], dydyed1);

			AVX_Cmp_DDY(dyed_buf_1+3,4*iXN,dydyed2);
			V4_1[1] = _mm256_add_ps(V4_1[1], dydyed2);
			dyed_buf_1 += 4;
			V4_1 += 34;
		}

		pq0 += stride_y_m256 - (iXN+2);
		DenAng0 += stride_y_em - (iXN+2)*4;
	}
}

//
// Valid from [1,dimz-1>
//
void AVX_TTIDenQ_Process_Patch(
	int logLevel,
	int stage,
	restrict __m256* pq,
	restrict __m256* rs,
	restrict __m256* Apq,			// pass pq ptr if 2nd order in time
	int* DenAng,
	int* VelAnis,
	restrict __m128* _mm_spgx,
	restrict __m128* _mm_spgy,
	__m128 _mm_spgz0,
	__m128 _mm_spgz1,
	restrict __m256* stripe,
	restrict __m256* dyed_buf,
	restrict __m256* V4,
	int iX0,
	int iY0,
	int iZ,
	int dimz,
	int iXN_halo,
	int iXN,
	int iYN_halo,
	int iYN,
	int stride_y_m256,
	int stride_z_m256,
	int stride_y_em,
	int stride_z_em
#ifdef TMJ_TIMING
	,unsigned long& thrcnt2
#endif
	)
{
	int dz_edge_diff = dimz - iZ - 1;  // dimz-1 -> 0, dimz-2 -> 1 etc.

	restrict __m256* pq1 = pq + (unsigned long)iZ * (unsigned long)stride_z_m256 + (unsigned long)iY0 * (unsigned long)stride_y_m256 + (unsigned long)(iX0 >> 1);
	restrict __m256* rs1 = rs + (unsigned long)iZ * (unsigned long)stride_z_m256 + (unsigned long)iY0 * (unsigned long)stride_y_m256 + (unsigned long)(iX0 >> 1);
	restrict __m256* Apq1 = Apq + (unsigned long)iZ * (unsigned long)stride_z_m256 + (unsigned long)iY0 * (unsigned long)stride_y_m256 + (unsigned long)(iX0 >> 1);
	int* VelAnis1 = VelAnis + (unsigned long)iZ * (unsigned long)stride_z_em + (unsigned long)iY0 * (unsigned long)stride_y_em + (unsigned long)(iX0*2);
	int* DenAng1 = DenAng + (unsigned long)iZ * (unsigned long)stride_z_em + (unsigned long)iY0 * (unsigned long)stride_y_em + (unsigned long)(iX0*2);

	if (logLevel >= 6)
	{
		printf("Process_Patch :: iZ = %d\n",iZ);
		fflush(stdout);
	}

	// process next tile ("patch")
	if (dz_edge_diff > 5)
	{
	 	restrict __m256* pq0 = ((stage == 0 || stage == 1) ? pq : Apq) + (unsigned long)(iZ+4) * (unsigned long)stride_z_m256 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m256 + (unsigned long)((iX0-4) >> 1);
		int* DenAng0 = DenAng + (unsigned long)(iZ+4) * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

		restrict __m256* dyed_buf_0 = dyed_buf;
		restrict __m256* dyed_buf_1 = dyed_buf + 16 * iXN;

		restrict __m256* V4_0 = V4;
		restrict __m256* V4_1 = V4;

		for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
		{
			pq0 += 2;
			DenAng0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m256), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m256+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m256), _MM_HINT_T0);

				__m256 ddx, ddy, ddz_z0, ddz_z1;
				AVX_Cmp_DDX(pq0,ddx);
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
				if (dz_edge_diff > 10)
				{
					AVX_Cmp_DDZ_2Z(pq0,stride_z_m256,ddz_z0,ddz_z1);
				}
				else if (dz_edge_diff > 9)
				{
					AVX_Cmp_DDZ(pq0,stride_z_m256,ddz_z0);
					AVX_Cmp_DDZ_EE(pq0+stride_z_m256,stride_z_m256,ddz_z1);
				}
				else
				{
					AVX_Cmp_DDZ_EE(pq0,stride_z_m256,ddz_z0);
					AVX_Cmp_DDZ_EE(pq0+stride_z_m256,stride_z_m256,ddz_z1);
				}

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DYED(
						ddx, ddy, ddz_z0, 
						Boy, cDip, sDip, cAzm, sAzm,
						dyed_buf_0[0], dyed_buf_0[1]
						);

				AVX_Cmp_DDX(pq0+stride_z_m256,ddx);
				AVX_Cmp_DDY(pq0+stride_z_m256,stride_y_m256,ddy);

				TTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DYED(
						ddx, ddy, ddz_z1, 
						Boy, cDip, sDip, cAzm, sAzm,
                                                dyed_buf_0[2], dyed_buf_0[3]
						);

				dyed_buf_0 += 4;
				pq0 += 1;
				DenAng0 += 4;
			}
			pq0 += stride_y_m256 - (iXN+2);
			DenAng0 += stride_y_em - (iXN+2)*4;
		}

		for (int iY = (iYN_halo-5);  iY > 4;  --iY)
		{
			__m256* stripe0 = stripe;

			for (int iX = 2;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m256), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m256+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m256), _MM_HINT_T0);

				__m256 ddx, ddy, ddz_z0, ddz_z1;
				AVX_Cmp_DDX(pq0,ddx);
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
				if (dz_edge_diff > 10)
				{
					AVX_Cmp_DDZ_2Z(pq0,stride_z_m256,ddz_z0,ddz_z1);
				}
				else if (dz_edge_diff > 9)
				{
					AVX_Cmp_DDZ(pq0,stride_z_m256,ddz_z0);
					AVX_Cmp_DDZ_EE(pq0+stride_z_m256,stride_z_m256,ddz_z1);
				}
				else
				{
					AVX_Cmp_DDZ_EE(pq0,stride_z_m256,ddz_z0);
					AVX_Cmp_DDZ_EE(pq0+stride_z_m256,stride_z_m256,ddz_z1);
				}

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DXED(
						ddx, ddy, ddz_z0, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 0], stripe0[ 1]
						);

				AVX_Cmp_DDX(pq0+stride_z_m256,ddx);
				AVX_Cmp_DDY(pq0+stride_z_m256,stride_y_m256,ddy);

				TTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DXED(
						ddx, ddy, ddz_z1, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 2], stripe0[ 3]
						);
				stripe0 += 4;

				pq0 += 1;
				DenAng0 += 4;
			}

			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m256), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m256+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m256), _MM_HINT_T0);

				__m256 ddx, ddy, ddz_z0, ddz_z1;
				AVX_Cmp_DDX(pq0,ddx);
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
				if (dz_edge_diff > 10)
				{
					AVX_Cmp_DDZ_2Z(pq0,stride_z_m256,ddz_z0,ddz_z1);
				}
				else if (dz_edge_diff > 9)
				{
					AVX_Cmp_DDZ(pq0,stride_z_m256,ddz_z0);
					AVX_Cmp_DDZ_EE(pq0+stride_z_m256,stride_z_m256,ddz_z1);
				}
				else
				{
					AVX_Cmp_DDZ_EE(pq0,stride_z_m256,ddz_z0);
					AVX_Cmp_DDZ_EE(pq0+stride_z_m256,stride_z_m256,ddz_z1);
				}

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DXED_DYED_DZED(
						ddx, ddy, ddz_z0, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 0], stripe0[ 1],
						dyed_buf_0[0], dyed_buf_0[1],
						V4_0[14], V4_0[15]
						);

				AVX_Cmp_DDX(pq0+stride_z_m256,ddx);
				AVX_Cmp_DDY(pq0+stride_z_m256,stride_y_m256,ddy);

				TTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DXED_DYED_DZED(
						ddx, ddy, ddz_z1, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 2], stripe0[ 3],
                                                dyed_buf_0[2], dyed_buf_0[3],
						V4_0[12], V4_0[13]
						);
				stripe0 += 4;
				V4_0 += 34;

				dyed_buf_0 += 4;
				pq0 += 1;
				DenAng0 += 4;
			}
			/*
			for (int iX = iXN;  iX > 0;  iX-=2)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m256), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m256+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m256), _MM_HINT_T0);

				__m256 ddz_x0_z0, ddz_x1_z0, ddz_x0_z1, ddz_x1_z1;
				if (dz_edge_diff > 10)
                                {
                                        AVX_Cmp_DDZ_2Z(pq0,stride_z_m256,ddz_x0_z0,ddz_x0_z1);
                                        AVX_Cmp_DDZ_2Z(pq0+1,stride_z_m256,ddz_x1_z0,ddz_x1_z1);
                                }
                                else if (dz_edge_diff > 9)
                                {
                                        AVX_Cmp_DDZ(pq0,stride_z_m256,ddz_x0_z0);
                                        AVX_Cmp_DDZ(pq0+1,stride_z_m256,ddz_x1_z0);
                                        AVX_Cmp_DDZ_EE(pq0+stride_z_m256,stride_z_m256,ddz_x0_z1);
                                        AVX_Cmp_DDZ_EE(pq0+1+stride_z_m256,stride_z_m256,ddz_x1_z1);
                                }
                                else
                                {
                                        AVX_Cmp_DDZ_EE(pq0,stride_z_m256,ddz_x0_z0);
                                        AVX_Cmp_DDZ_EE(pq0+1,stride_z_m256,ddz_x1_z0);
                                        AVX_Cmp_DDZ_EE(pq0+stride_z_m256,stride_z_m256,ddz_x0_z1);
                                        AVX_Cmp_DDZ_EE(pq0+1+stride_z_m256,stride_z_m256,ddz_x1_z1);
                                }

				__m256 ddx_x0_z0, ddx_x1_z0, ddx_x0_z1, ddx_x1_z1;
				AVX_Cmp_DDX_2X(pq0,ddx_x0_z0,ddx_x1_z0);
				AVX_Cmp_DDX_2X(pq0+stride_z_m256,ddx_x0_z1,ddx_x1_z1);

				// x0, z0
				__m256 ddy;
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);

				__m128 Boy, cDip, sDip, cAzm, sAzm;
                                TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

                                AVX_TTIDenQ_Cmp_DXED_DYED_DZED(
                                                ddx_x0_z0, ddy, ddz_x0_z0,
                                                Boy, cDip, sDip, cAzm, sAzm,
                                                stripe0[ 0], stripe0[ 1],
                                                dyed_buf_0[0], dyed_buf_0[1],
                                                V4_0[14], V4_0[15]
                                                );

				// x0, z1
				AVX_Cmp_DDY(pq0+stride_z_m256,stride_y_m256,ddy);

				TTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DXED_DYED_DZED(
                                                ddx_x0_z1, ddy, ddz_x0_z1,
                                                Boy, cDip, sDip, cAzm, sAzm,
                                                stripe0[ 2], stripe0[ 3],
                                                dyed_buf_0[2], dyed_buf_0[3],
                                                V4_0[12], V4_0[13]
                                                );

				// x1, z0
				AVX_Cmp_DDY(pq0+1,stride_y_m256,ddy);

                                TTIDenQ_Get_EM1(DenAng0+4, Boy, cDip, sDip, cAzm, sAzm);

                                AVX_TTIDenQ_Cmp_DXED_DYED_DZED(
                                                ddx_x1_z0, ddy, ddz_x1_z0,
                                                Boy, cDip, sDip, cAzm, sAzm,
                                                stripe0[ 4], stripe0[ 5],
                                                dyed_buf_0[4], dyed_buf_0[5],
                                                V4_0[14+34], V4_0[15+34]
                                                );

				// x1, z1
				AVX_Cmp_DDY(pq0+1+stride_z_m256,stride_y_m256,ddy);

				TTIDenQ_Get_EM1(DenAng0+4+stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DXED_DYED_DZED(
                                                ddx_x1_z1, ddy, ddz_x1_z1,
                                                Boy, cDip, sDip, cAzm, sAzm,
                                                stripe0[ 2+4], stripe0[ 3+4],
                                                dyed_buf_0[2+4], dyed_buf_0[3+4],
                                                V4_0[12+34], V4_0[13+34]
                                                );

				stripe0 += 8;
				V4_0 += 68;

				dyed_buf_0 += 8;
				pq0 += 2;
				DenAng0 += 8;
			}
			*/
			V4_0 -= 34 * iXN;

			for (int iX = 1;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m256), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m256+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m256), _MM_HINT_T0);

				__m256 ddx, ddy, ddz_z0, ddz_z1;
				AVX_Cmp_DDX(pq0,ddx);
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
				if (dz_edge_diff > 10)
				{
					AVX_Cmp_DDZ_2Z(pq0,stride_z_m256,ddz_z0,ddz_z1);
				}
				else if (dz_edge_diff > 9)
				{
					AVX_Cmp_DDZ(pq0,stride_z_m256,ddz_z0);
					AVX_Cmp_DDZ_EE(pq0+stride_z_m256,stride_z_m256,ddz_z1);
				}
				else
				{
					AVX_Cmp_DDZ_EE(pq0,stride_z_m256,ddz_z0);
					AVX_Cmp_DDZ_EE(pq0+stride_z_m256,stride_z_m256,ddz_z1);
				}

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DXED(
						ddx, ddy, ddz_z0, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 0], stripe0[ 1]
						);

				AVX_Cmp_DDX(pq0+stride_z_m256,ddx);
				AVX_Cmp_DDY(pq0+stride_z_m256,stride_y_m256,ddy);

				TTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DXED(
						ddx, ddy, ddz_z1, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 2], stripe0[ 3]
						);
				stripe0 += 4;

				pq0 += 1;
				DenAng0 += 4;
			}

			pq0 += stride_y_m256 - iXN_halo;
			DenAng0 += stride_y_em - iXN_halo*4;

			stripe0 = stripe + 8;
			__m128 _mm_spgzy0 = _mm_mul_ps(_mm_spgz0, _mm_spgy[iY]);
			__m128 _mm_spgzy1 = _mm_mul_ps(_mm_spgz1, _mm_spgy[iY]);
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng1), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng1+stride_z_em), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelAnis1), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelAnis1+stride_z_em), _MM_HINT_NTA);
				if (stage == 0)
				{
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m256), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m256), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m256), _MM_HINT_T0);
				}

				if (iY <= (iYN_halo-9))
				{	
					// dydyed1, dydyed2
					__m256 dydyed1;
					AVX_Cmp_DDY(dyed_buf_1,4*iXN,dydyed1);
					V4_1[2] = _mm256_add_ps(V4_1[2], dydyed1);

					__m256 dydyed2;
					AVX_Cmp_DDY(dyed_buf_1+1,4*iXN,dydyed2);
					V4_1[3] = _mm256_add_ps(V4_1[3], dydyed2);

					AVX_Cmp_DDY(dyed_buf_1+2,4*iXN,dydyed1);
					V4_1[0] = _mm256_add_ps(V4_1[0], dydyed1);

					AVX_Cmp_DDY(dyed_buf_1+3,4*iXN,dydyed2);
					V4_1[1] = _mm256_add_ps(V4_1[1], dydyed2);
					dyed_buf_1 += 4;
					V4_1 += 34;
				}

				// dxdxed1, dxdxed2
				__m256 dxdxed1;
				AVX_Cmp_DDX_1x_stride8(stripe0,dxdxed1);
				V4_0[2] = dxdxed1;

				__m256 dxdxed2;
				AVX_Cmp_DDX_1x_stride8(stripe0+1,dxdxed2);
				V4_0[3] = dxdxed2;

				AVX_Cmp_DDX_1x_stride8(stripe0+2,dxdxed1);
				V4_0[0] = dxdxed1;

				AVX_Cmp_DDX_1x_stride8(stripe0+3,dxdxed2);
				V4_0[1] = dxdxed2;

				// dzdzed1, dzdzed2
				__m256 dzdzed1_z0, dzdzed1_z1;
				AVX_Cmp_DDZ_2Z(V4_0+24,-2,dzdzed1_z0, dzdzed1_z1);
				__m256 V5 = _mm256_sub_ps(V4_0[10], dzdzed1_z0);

				__m256 dzdzed2_z0, dzdzed2_z1;
				AVX_Cmp_DDZ_2Z(V4_0+25,-2,dzdzed2_z0, dzdzed2_z1);
				__m256 V4 = _mm256_add_ps(V4_0[11], dzdzed2_z0);

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				AVX_TTIDenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelAnis1,DenAng1,spgzyx0,V4,V5
						);

				V5 = _mm256_sub_ps(V4_0[ 8], dzdzed1_z1);

				V4 = _mm256_add_ps(V4_0[ 9], dzdzed2_z1);

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				AVX_TTIDenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m256,rs1+stride_z_m256,Apq1+stride_z_m256,VelAnis1+stride_z_em,DenAng1+stride_z_em,spgzyx1,V4,V5
						);

				V4_0 += 34;

				pq1 += 1;
				rs1 += 1;
				Apq1 += 1;
				VelAnis1 += 4;
				DenAng1 += 4;
				stripe0 += 4;	
#ifdef TMJ_TIMING
				thrcnt2+=2;
#endif
			}
			pq1 += stride_y_m256 - iXN;
			rs1 += stride_y_m256 - iXN;
			Apq1 += stride_y_m256 - iXN;
			VelAnis1 += stride_y_em - iXN*4;
			DenAng1 += stride_y_em - iXN*4;
		}

		for (int iY = 4;  iY > 0;  --iY)
		{
			pq0 += 2;
			DenAng0 += 8;	
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m256), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m256+5*stride_y_m256), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m256), _MM_HINT_T0);

				__m256 ddx, ddy, ddz_z0, ddz_z1;
				AVX_Cmp_DDX(pq0,ddx);
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
				if (dz_edge_diff > 10)
				{
					AVX_Cmp_DDZ_2Z(pq0,stride_z_m256,ddz_z0,ddz_z1);
				}
				else if (dz_edge_diff > 9)
				{
					AVX_Cmp_DDZ(pq0,stride_z_m256,ddz_z0);
					AVX_Cmp_DDZ_EE(pq0+stride_z_m256,stride_z_m256,ddz_z1);
				}
				else
				{
					AVX_Cmp_DDZ_EE(pq0,stride_z_m256,ddz_z0);
					AVX_Cmp_DDZ_EE(pq0+stride_z_m256,stride_z_m256,ddz_z1);
				}

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DYED(
						ddx, ddy, ddz_z0, 
						Boy, cDip, sDip, cAzm, sAzm,
						dyed_buf_0[0], dyed_buf_0[1]
						);

				AVX_Cmp_DDX(pq0+stride_z_m256,ddx);
				AVX_Cmp_DDY(pq0+stride_z_m256,stride_y_m256,ddy);

				TTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DYED(
						ddx, ddy, ddz_z1, 
						Boy, cDip, sDip, cAzm, sAzm,
                                                dyed_buf_0[2], dyed_buf_0[3]
						);

				dyed_buf_0 += 4;
				pq0 += 1;
				DenAng0 += 4;

				// dydyed1, dydyed2
				__m256 dydyed1;
				AVX_Cmp_DDY(dyed_buf_1,4*iXN,dydyed1);
				V4_1[2] = _mm256_add_ps(V4_1[2], dydyed1);

				__m256 dydyed2;
				AVX_Cmp_DDY(dyed_buf_1+1,4*iXN,dydyed2);
				V4_1[3] = _mm256_add_ps(V4_1[3], dydyed2);

				AVX_Cmp_DDY(dyed_buf_1+2,4*iXN,dydyed1);
				V4_1[0] = _mm256_add_ps(V4_1[0], dydyed1);

				AVX_Cmp_DDY(dyed_buf_1+3,4*iXN,dydyed2);
				V4_1[1] = _mm256_add_ps(V4_1[1], dydyed2);
				dyed_buf_1 += 4;
				V4_1 += 34;
			}
	
			pq0 += stride_y_m256 - (iXN+2);
			DenAng0 += stride_y_em - (iXN+2)*4;
		}
	}
	else if (dz_edge_diff == 5)
	{
	 	__m256* pq0 = ((stage == 0 || stage == 1) ? pq : Apq) + (unsigned long)(iZ+4) * (unsigned long)stride_z_m256 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m256 + (unsigned long)((iX0-4) >> 1);
		int* DenAng0 = DenAng + (unsigned long)(iZ+4) * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

		__m256* dyed_buf_0 = dyed_buf;
		__m256* dyed_buf_1 = dyed_buf + 16 * iXN;

		__m256* V4_0 = V4;
		__m256* V4_1 = V4;

		for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
		{
			pq0 += 1;
			DenAng0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m256 ddx, ddy, ddz;
				AVX_Cmp_DDX(pq0,ddx);
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
				AVX_Cmp_DDZ_EE(pq0,stride_z_m256,ddz);

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DYED(
						ddx, ddy, ddz, 
						Boy, cDip, sDip, cAzm, sAzm,
						dyed_buf_0[0], dyed_buf_0[1]
						);
				dyed_buf_0 += 4;

				pq0 += 1;
				DenAng0 += 4;
			}

			pq0 += stride_y_m256 - (iXN+2);
			DenAng0 += stride_y_em - (iXN+2)*4;
		}

		for (int iY = (iYN_halo-5);  iY > 4;  --iY)
		{
			__m256* stripe0 = stripe;

			for (int iX = 2;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m256 ddx, ddy, ddz;
				AVX_Cmp_DDX(pq0,ddx);
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
				AVX_Cmp_DDZ_EE(pq0,stride_z_m256,ddz);

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DXED(
						ddx, ddy, ddz, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 0], stripe0[ 1]
						);
				stripe0 += 4;

				pq0 += 1;
				DenAng0 += 4;
			}

			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m256 ddx, ddy, ddz;
				AVX_Cmp_DDX(pq0,ddx);
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
				AVX_Cmp_DDZ_EE(pq0,stride_z_m256,ddz);

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DXED_DYED_DZED(
						ddx, ddy, ddz, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 0], stripe0[ 1], 
						dyed_buf_0[0], dyed_buf_0[1],
						V4_0[14], V4_0[15]
						);
				stripe0 += 4;
				dyed_buf_0 += 4;
				V4_0 += 34;

				pq0 += 1;
				DenAng0 += 4;
			}
			V4_0 -= 34 * iXN;

			for (int iX = 1;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m256 ddx, ddy, ddz;
				AVX_Cmp_DDX(pq0,ddx);
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
				AVX_Cmp_DDZ_EE(pq0,stride_z_m256,ddz);

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DXED(
						ddx, ddy, ddz, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 0], stripe0[ 1]
						);
				stripe0 += 4;

				pq0 += 1;
				DenAng0 += 4;
			}

			pq0 += stride_y_m256 - iXN_halo;
			DenAng0 += stride_y_em - iXN_halo*4;

			stripe0 = stripe + 8;
			__m128 _mm_spgzy0 = _mm_mul_ps(_mm_spgz0, _mm_spgy[iY]);
			__m128 _mm_spgzy1 = _mm_mul_ps(_mm_spgz1, _mm_spgy[iY]);
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng1), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng1+stride_z_em), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelAnis1), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelAnis1+stride_z_em), _MM_HINT_NTA);
				if (stage == 0)
				{
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m256), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m256), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m256), _MM_HINT_T0);
				}

				if (iY <= (iYN_halo-9))
				{
					// dydyed1, dydyed2
					__m256 dydyed1;
					AVX_Cmp_DDY(dyed_buf_1,4*iXN,dydyed1);
					V4_1[2] = _mm256_add_ps(V4_1[2], dydyed1);

					__m256 dydyed2;
					AVX_Cmp_DDY(dyed_buf_1+1,4*iXN,dydyed2);
					V4_1[3] = _mm256_add_ps(V4_1[3], dydyed2);
					dyed_buf_1 += 4;
					V4_1 += 34;
				}

				// dxdxed1, dxdxed2
				__m256 dxdxed1;
				AVX_Cmp_DDX_1x_stride8(stripe0,dxdxed1);
				V4_0[2] = dxdxed1;

				__m256 dxdxed2;
				AVX_Cmp_DDX_1x_stride8(stripe0+1,dxdxed2);
				V4_0[3] = dxdxed2;

				// dzdzed1, dzdzed2
				__m256 dzdzed1;
				AVX_Cmp_DDZ(V4_0+24,-2,dzdzed1);
				__m256 V5 = _mm256_sub_ps(V4_0[10], dzdzed1);

				__m256 dzdzed2;
				AVX_Cmp_DDZ(V4_0+25,-2,dzdzed2);
				__m256 V4 = _mm256_add_ps(V4_0[11], dzdzed2);

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				AVX_TTIDenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelAnis1,DenAng1,spgzyx0,V4,V5
						);

				AVX_Cmp_DDZ_EE(V4_0+22,-2,dzdzed1);
				V5 = _mm256_sub_ps(V4_0[ 8], dzdzed1);

				AVX_Cmp_DDZ_EE(V4_0+23,-2,dzdzed2);
				V4 = _mm256_add_ps(V4_0[ 9], dzdzed2);

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				AVX_TTIDenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m256,rs1+stride_z_m256,Apq1+stride_z_m256,VelAnis1+stride_z_em,DenAng1+stride_z_em,spgzyx1,V4,V5
						);

				V4_0 += 34;

				pq1 += 1;
				rs1 += 1;
				Apq1 += 1;
				VelAnis1 += 4;
				DenAng1 += 4;
				stripe0 += 4;	
#ifdef TMJ_TIMING
				thrcnt2+=2;
#endif
			}
			pq1 += stride_y_m256 - iXN;
			rs1 += stride_y_m256 - iXN;
			Apq1 += stride_y_m256 - iXN;
			VelAnis1 += stride_y_em - iXN*4;
			DenAng1 += stride_y_em - iXN*4;
		}

		for (int iY = 4;  iY > 0;  --iY)
		{
			pq0 += 2;
			DenAng0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m256 ddx, ddy, ddz;
				AVX_Cmp_DDX(pq0,ddx);
				AVX_Cmp_DDY(pq0,stride_y_m256,ddy);
				AVX_Cmp_DDZ_EE(pq0,stride_z_m256,ddz);

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				AVX_TTIDenQ_Cmp_DYED(
						ddx, ddy, ddz, 
						Boy, cDip, sDip, cAzm, sAzm,
						dyed_buf_0[0], dyed_buf_0[1]
						);
				dyed_buf_0 += 4;

				pq0 += 1;
				DenAng0 += 4;

				// dydyed1, dydyed2
				__m256 dydyed1;
				AVX_Cmp_DDY(dyed_buf_1,4*iXN,dydyed1);
				V4_1[2] = _mm256_add_ps(V4_1[2], dydyed1);

				__m256 dydyed2;
				AVX_Cmp_DDY(dyed_buf_1+1,4*iXN,dydyed2);
				V4_1[3] = _mm256_add_ps(V4_1[3], dydyed2);
				dyed_buf_1 += 4;
				V4_1 += 34;
			}

			pq0 += stride_y_m256 - (iXN+2);
			DenAng0 += stride_y_em - (iXN+2)*4;
		}
	}
	else // no more tiles to process, finish rest of outputs from queue
	{
		__m256* V4_0 = V4;
		for (int iY = iYN;  iY > 0;  --iY)
		{
			__m128 _mm_spgzy0 = _mm_mul_ps(_mm_spgz0, _mm_spgy[iY+4]);
			__m128 _mm_spgzy1 = _mm_mul_ps(_mm_spgz1, _mm_spgy[iY+4]);
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng1), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng1+stride_z_em), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelAnis1), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelAnis1+stride_z_em), _MM_HINT_NTA);
				if (stage == 0)
				{
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m256), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m256), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m256), _MM_HINT_T0);
				}

				// dzdzed1, dzdzed2
				// use shorter two point stencil
				__m256 dzdzed1;
				AVX_Cmp_DDZ_EE(V4_0+24,-2,dzdzed1);
				__m256 V5 = _mm256_sub_ps(V4_0[10], dzdzed1);

				__m256 dzdzed2;
				AVX_Cmp_DDZ_EE(V4_0+25,-2,dzdzed2);
				__m256 V4 = _mm256_add_ps(V4_0[11], dzdzed2);

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				AVX_TTIDenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelAnis1,DenAng1,spgzyx0,V4,V5
						);

				AVX_Cmp_DDZ_EE(V4_0+22,-2,dzdzed1);
				V5 = _mm256_sub_ps(V4_0[ 8], dzdzed1);

				AVX_Cmp_DDZ_EE(V4_0+23,-2,dzdzed2);
				V4 = _mm256_add_ps(V4_0[ 9], dzdzed2);

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				AVX_TTIDenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m256,rs1+stride_z_m256,Apq1+stride_z_m256,VelAnis1+stride_z_em,DenAng1+stride_z_em,spgzyx1,V4,V5
						);

				V4_0 += 34;

				pq1 += 1;
				rs1 += 1;
				Apq1 += 1;
				VelAnis1 += 4;
				DenAng1 += 4;
#ifdef TMJ_TIMING
				thrcnt2 += 2;
#endif
			}
			pq1 += stride_y_m256 - iXN;
			rs1 += stride_y_m256 - iXN;
			Apq1 += stride_y_m256 - iXN;
			VelAnis1 += stride_y_em - iXN*4;
			DenAng1 += stride_y_em - iXN*4;
		}
	}
}

} // end of anonymous namespace

void AVX_VTIDenQ_TimeStep(
	int logLevel,
	int stage,
	__m128* pq,
	__m128* rs,
	__m128* Apq,
	int* VelAnis,
        int* DenAng,
	float* spgx,
	float* spgy,
	float* spgz,
        int dimx,
        int dimy,
        int dimz,
	int xh,
	int yh,
	int zh,
        int bsX,
        int bsY,
	int x0,
	int x1,
	int y0,
	int y1,
	int z0,
	int z1
        )
{
	// adjust x0 so that it is always a multiple of MIN_BSX
	int xx0 = x0 >> 2;
	int xx1 = x1 >> 2;

	int dimxh = dimx + 2*xh;
	int dimyh = dimy + 2*yh;
	int dimzh = dimz + 2*zh;

        int stride_y_m128 = dimxh / 2;
        int stride_z_m128 = stride_y_m128 * dimyh;
        int dimx_2 = dimx / 2;

	int stride_y_em = dimxh;
	int stride_z_em = stride_y_em * dimyh;

	int num_threads = 0;
#pragma omp parallel
        {
                num_threads = omp_get_num_threads();
        }

	//
	// tile XY face
	//

	int bsXX = ((bsX + 7) >> 3) << 3;
	if (bsXX > bsX)
	{
		printf("WARNING: Block Size X must be a multiple of 8. Changed from %d to %d\n",bsX,bsXX);
		bsX = bsXX;
	}

	int nnY = y1 - y0 + 1;
	int nnX = xx1 - xx0 + 1;

	int bsX_2 = bsX >> 1;
	int bsX_4 = bsX >> 2;
	int tileX = (nnX+bsX_4-1) / bsX_4;
	int tileY = (nnY+bsY-1) / bsY;
	int tileNN = tileX * tileY;

	if (logLevel >= 5)
	{
		printf("Number of tiles is %d\n",tileNN);
		fflush(stdout);
	}

	//
	// thread queues
	//

	unsigned long stripe_m128, dyed_buf_m128, spg_buf_m128, QLen, QSize, BufSize_m128, NetBufSize_m128;
        int QMaxShift;
	VTIDenQ_Comp_Buf_Size(bsX, bsY, 128, stripe_m128, dyed_buf_m128, spg_buf_m128, QMaxShift, QLen, QSize, BufSize_m128, NetBufSize_m128);

	if (BufSize_m128*(unsigned long)num_threads > _BufSize)
	{
		unsigned long NewBufSize = (unsigned long)num_threads * BufSize_m128 * 16;
		if (_Buf != 0L) free(_Buf);
		posix_memalign((void**)&_Buf, 64, NewBufSize);
		_BufSize = BufSize_m128 * (unsigned long)num_threads;
	}

	if (logLevel >= 5)
	{
		printf("Cached buffer size is %.2f KB\n", (double)((unsigned long)num_threads * NetBufSize_m128 * 16) / 1024.0);
		fflush(stdout);
	}

#ifdef TMJ_TIMING
	unsigned long thrcnt2 = 0;
	struct timespec ts0, ts1;
	clock_gettime(CLOCK_REALTIME, &ts0);
	double wall_time_0 = (double)ts0.tv_sec + (double)ts0.tv_nsec * 1e-9;
#endif

	// skip halos
	__m128* pqnh = pq + zh * stride_z_m128 + yh * stride_y_m128 + (xh >> 1);
	__m128* rsnh = rs + zh * stride_z_m128 + yh * stride_y_m128 + (xh >> 1);
	__m128* Apqnh = Apq + zh * stride_z_m128 + yh * stride_y_m128 + (xh >> 1);
	
	int* DenAng_nh = DenAng + zh * stride_z_em + yh * stride_y_em + xh;
	int* VelAnis_nh = VelAnis + zh * stride_z_em + yh * stride_y_em + xh;

#ifdef TMJ_TIMING
#pragma omp parallel for schedule(dynamic) reduction(+:thrcnt2)
#else
#pragma omp parallel for schedule(dynamic)
#endif
	for (int i = 0;  i < tileNN;  ++i)
	{
		int iTileX = i % tileX;
		int iTileY = i / tileX;

		// skip corner tiles because we don't know how to handle these yet.
		//if (iTileX == 0 || iTileX == (tileX-1) || iTileY == 0 || iTileY == (tileY-1)) continue; 

		int iY0 = y0 + iTileY * bsY;
		int iY1 = iY0 + bsY - 1;
		if (iY1 > y1) iY1 = y1;
		int iYN = iY1 - iY0 + 1;
		int iYN_halo = iYN + 9;

		int iX0 = xx0 + iTileX * bsX_4;
		int iX1 = iX0 + bsX_4 - 1;
		if (iX1 > xx1) iX1 = xx1;

		int iXN = iX1 - iX0 + 1;
		int iXN_halo = iXN + 3;

		int dyed_stride_m128 = iXN * 4;

		if (logLevel >= 5 && (i & 15) == 1)
		{
			printf("tile %d :: x=[%d,%d] iXN=%d y=[%d,%d] iYN=%d\n",i,iX0*4,iX1*4+3,iXN,iY0,iY1,iYN);
			fflush(stdout);
		}

		// get thread buffer
		int thread_num;
		thread_num = omp_get_thread_num();
		__m128* thrBuf = _Buf + (unsigned long)thread_num * BufSize_m128;
		__m128* stripe = thrBuf;
		__m128* dyed_buf = stripe + stripe_m128;
		__m128* spg_buf = dyed_buf + dyed_buf_m128;
		__m128* Q = spg_buf + spg_buf_m128;
		int QShift = QMaxShift-1;
		__m128* V4 = Q + (QShift<<1);

		// copy the X sponges into properly aligned array
		// also flip the order so it matches the iX values during iterations
		__m128* _mm_spgx = spg_buf;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			int x0 = (iX0<<2) + ((iXN-iX) << 2);
			_mm_spgx[iX] = _mm_setr_ps(spgx[x0],spgx[x0+1],spgx[x0+2],spgx[x0+3]);
		}
		_mm_spgx[0] = _mm_set1_ps(0.0f);

		__m128* _mm_spgy = _mm_spgx + iXN + 1;
		for (int iY = iYN_halo;  iY > 0;  --iY)
		{
			if (iY >= 5 &&  iY <= (iYN_halo-5))
			{
				int yy = iYN_halo - iY + iY0 - 5;
				_mm_spgy[iY] = _mm_set1_ps(spgy[yy]);
			}
			else
			{
				_mm_spgy[iY] = _mm_set1_ps(0.0f);
			}
		}
		_mm_spgy[0] = _mm_set1_ps(0.0f);

		//
		// Do Z leadin. Process_Patch_Leadin correctly handles negative Zzz's.
		//
		for (int iZ = z0-5;  iZ < z0+4;  iZ+=2)
		{
			VTIDenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);
			VTIDenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);
			AVX_VTIDenQ_Process_Patch_Leadin(
				logLevel,
				(__m256*)((stage == 0 || stage == 1) ? pqnh : Apqnh),DenAng_nh,VelAnis_nh,
				(__m256*)stripe,(__m256*)dyed_buf,(__m256*)V4,
				iX0*2,iY0,iZ,iXN_halo,iXN,iYN_halo,iYN,
				stride_y_m128/2,stride_z_m128/2,stride_y_em,stride_z_em);
		}
		VTIDenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);	

		//
		// Internal Z cells and lead-out
		//
		for (int iZ = z0;  iZ <= z1;  iZ+=2)
		{
			if (logLevel >= 6)
			{
				printf("thread %d - tile %d : iZ = %d\n",thread_num,i,iZ);
				fflush(stdout);
			}

			__m128 _mm_spgz0 = _mm_set1_ps(spgz[iZ]);
			__m128 _mm_spgz1 = _mm_set1_ps(spgz[iZ+1]);
			AVX_VTIDenQ_Process_Patch(
				logLevel,
				stage,(__m256*)pqnh,(__m256*)rsnh,(__m256*)Apqnh,DenAng_nh,VelAnis_nh,_mm_spgx,_mm_spgy,_mm_spgz0,_mm_spgz1,
				(__m256*)stripe,(__m256*)dyed_buf,(__m256*)V4,
				iX0*2,iY0,iZ,dimz,iXN_halo,iXN,iYN_halo,iYN,stride_y_m128/2,stride_z_m128/2,stride_y_em,stride_z_em
#ifdef TMJ_TIMING
				,thrcnt2
#endif
				);
			VTIDenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);
			VTIDenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);
		}
	}
#ifdef TMJ_TIMING
	clock_gettime(CLOCK_REALTIME, &ts1);
	double wall_time_1 = (double)ts1.tv_sec + (double)ts1.tv_nsec * 1e-9;
	double elapsed_time_1 = wall_time_1 - wall_time_0;

	double overcomp, flops_per_cell, gflops, eff_freq, net_comp_throughput;
	TTIDenQ_Compute_Performance((double)thrcnt2*4e-6,elapsed_time_1,bsX,bsY,stage,overcomp,flops_per_cell,gflops,eff_freq,net_comp_throughput);

	printf("elapsed_time = %.2f, cnt2 = %ld, overcomp = %.2f, Net Comp Throughput = %.0f MCells/s, GFLOPS = %.2f, Effective Clock Frequency = %.2f\n", elapsed_time_1, thrcnt2, overcomp, net_comp_throughput, gflops, eff_freq);
#endif

	//if (Buf != 0L) free((void*)Buf);
}

void AVX_TTIDenQ_TimeStep(
	int logLevel,
	int stage,
	__m128* pq,
	__m128* rs,
	__m128* Apq,
	int* VelAnis,
        int* DenAng,
	float* spgx,
	float* spgy,
	float* spgz,
        int dimx,
        int dimy,
        int dimz,
	int xh,
	int yh,
	int zh,
        int bsX,
        int bsY,
	int x0,
	int x1,
	int y0,
	int y1,
	int z0,
	int z1
        )
{
	if (logLevel >= 5)
	{
		printf("AVX_TTIDenQ_TimeStep\n");
		fflush(stdout);
	}

	// adjust x0 so that it is always a multiple of MIN_BSX
	int xx0 = x0 >> 2;
	int xx1 = x1 >> 2;

	int dimxh = dimx + 2*xh;
	int dimyh = dimy + 2*yh;
	int dimzh = dimz + 2*zh;

        int stride_y_m128 = dimxh / 2;
        int stride_z_m128 = stride_y_m128 * dimyh;
        int dimx_2 = dimx / 2;

	int stride_y_em = dimxh;
	int stride_z_em = stride_y_em * dimyh;

	int num_threads = 0;
#pragma omp parallel
        {
                num_threads = omp_get_num_threads();
        }

	//
	// tile XY face
	//

	int bsXX = ((bsX + 7) >> 3) << 3;
	if (bsXX > bsX)
	{
		printf("WARNING: Block Size X must be a multiple of 8. Changed from %d to %d\n",bsX,bsXX);
		bsX = bsXX;
	}

	int nnY = y1 - y0 + 1;
	int nnX = xx1 - xx0 + 1;

	int bsX_2 = bsX >> 1;
	int bsX_4 = bsX >> 2;
	int tileX = (nnX+bsX_4-1) / bsX_4;
	int tileY = (nnY+bsY-1) / bsY;
	int tileNN = tileX * tileY;

	if (logLevel >= 5)
	{
		printf("Number of tiles is %d\n",tileNN);
		fflush(stdout);
	}

	//
	// thread queues
	//

	unsigned long stripe_m128, dyed_buf_m128, spg_buf_m128, QLen, QSize, BufSize_m128, NetBufSize_m128;
        int QMaxShift;
	TTIDenQ_Comp_Buf_Size(bsX, bsY, 128, stripe_m128, dyed_buf_m128, spg_buf_m128, QMaxShift, QLen, QSize, BufSize_m128, NetBufSize_m128);

	if (BufSize_m128*(unsigned long)num_threads > _BufSize)
	{
		unsigned long NewBufSize = (unsigned long)num_threads * BufSize_m128 * 16;
		if (_Buf != 0L) free(_Buf);
		posix_memalign((void**)&_Buf, 64, NewBufSize);
		_BufSize = BufSize_m128 * (unsigned long)num_threads;
	}

	if (logLevel >= 5)
	{
		printf("Cached buffer size is %.2f KB\n", (double)((unsigned long)num_threads * NetBufSize_m128 * 16) / 1024.0);
		fflush(stdout);
	}

#ifdef TMJ_TIMING
	unsigned long thrcnt2 = 0;
	struct timespec ts0, ts1;
	clock_gettime(CLOCK_REALTIME, &ts0);
	double wall_time_0 = (double)ts0.tv_sec + (double)ts0.tv_nsec * 1e-9;
#endif

	// skip halos
	__m128* pqnh = pq + zh * stride_z_m128 + yh * stride_y_m128 + (xh >> 1);
	__m128* rsnh = rs + zh * stride_z_m128 + yh * stride_y_m128 + (xh >> 1);
	__m128* Apqnh = Apq + zh * stride_z_m128 + yh * stride_y_m128 + (xh >> 1);
	
	int* DenAng_nh = DenAng + zh * stride_z_em + yh * stride_y_em + xh;
	int* VelAnis_nh = VelAnis + zh * stride_z_em + yh * stride_y_em + xh;

#ifdef TMJ_TIMING
#pragma omp parallel for schedule(dynamic) reduction(+:thrcnt2) 
#else
#pragma omp parallel for schedule(dynamic) 
#endif
	for (int i = 0;  i < tileNN;  ++i)
	{
		int iTileX = i % tileX;
		int iTileY = i / tileX;

		int iY0 = y0 + iTileY * bsY;
		int iY1 = iY0 + bsY - 1;
		if (iY1 > y1) iY1 = y1;
		int iYN = iY1 - iY0 + 1;
		int iYN_halo = iYN + 9;

		int iX0 = xx0 + iTileX * bsX_4;
		int iX1 = iX0 + bsX_4 - 1;
		if (iX1 > xx1) iX1 = xx1;

		int iXN = iX1 - iX0 + 1;
		int iXN_halo = iXN + 3;

		int dyed_stride_m128 = iXN * 4;

		if (logLevel >= 5 && (i & 15) == 1)
		{
			printf("tile %d :: x=[%d,%d] iXN=%d y=[%d,%d] iYN=%d\n",i,iX0*4,iX1*4+3,iXN,iY0,iY1,iYN);
			fflush(stdout);
		}

		// get thread buffer
		int thread_num;
		thread_num = omp_get_thread_num();
		__m128* thrBuf = _Buf + (unsigned long)thread_num * BufSize_m128;
		__m128* stripe = thrBuf;
		__m128* dyed_buf = stripe + stripe_m128;
		__m128* spg_buf = dyed_buf + dyed_buf_m128;
		__m128* Q = spg_buf + spg_buf_m128;
		int QShift = QMaxShift-1;
		__m128* V4 = Q + (QShift<<2);

		// copy the X sponges into properly aligned array
		// also flip the order so it matches the iX values during iterations
		__m128* _mm_spgx = spg_buf;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			int x0 = (iX0<<2) + ((iXN-iX) << 2);
			_mm_spgx[iX] = _mm_setr_ps(spgx[x0],spgx[x0+1],spgx[x0+2],spgx[x0+3]);
		}
		_mm_spgx[0] = _mm_set1_ps(0.0f);

		__m128* _mm_spgy = _mm_spgx + iXN + 1;
		for (int iY = iYN_halo;  iY > 0;  --iY)
		{
			if (iY >= 5 &&  iY <= (iYN_halo-5))
			{
				int yy = iYN_halo - iY + iY0 - 5;
				_mm_spgy[iY] = _mm_set1_ps(spgy[yy]);
			}
			else
			{
				_mm_spgy[iY] = _mm_set1_ps(0.0f);
			}
		}
		_mm_spgy[0] = _mm_set1_ps(0.0f);

		//
		// Do Z leadin. Process_Patch_Leadin correctly handles negative Zzz's.
		//
		for (int iZ = z0-5;  iZ < z0+4;  iZ+=2)
		{
			TTIDenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);
			TTIDenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);
			AVX_TTIDenQ_Process_Patch_Leadin(
				logLevel,
				(__m256*)((stage == 0 || stage == 1) ? pqnh : Apqnh),DenAng_nh,VelAnis_nh,
				(__m256*)stripe,(__m256*)dyed_buf,(__m256*)V4,
				iX0*2,iY0,iZ,iXN_halo,iXN,iYN_halo,iYN,
				stride_y_m128/2,stride_z_m128/2,stride_y_em,stride_z_em);
		}
		TTIDenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);	

		//
		// Internal Z cells and lead-out
		//
		for (int iZ = z0;  iZ <= z1;  iZ+=2)
		{
			if (logLevel >= 6)
			{
				printf("thread %d - tile %d : iZ = %d\n",thread_num,i,iZ);
				fflush(stdout);
			}

			__m128 _mm_spgz0 = _mm_set1_ps(spgz[iZ]);
			__m128 _mm_spgz1 = _mm_set1_ps(spgz[iZ+1]);
			AVX_TTIDenQ_Process_Patch(
				logLevel,
				stage,(__m256*)pqnh,(__m256*)rsnh,(__m256*)Apqnh,DenAng_nh,VelAnis_nh,_mm_spgx,_mm_spgy,_mm_spgz0,_mm_spgz1,
				(__m256*)stripe,(__m256*)dyed_buf,(__m256*)V4,
				iX0*2,iY0,iZ,dimz,iXN_halo,iXN,iYN_halo,iYN,stride_y_m128/2,stride_z_m128/2,stride_y_em,stride_z_em
#ifdef TMJ_TIMING
				,thrcnt2
#endif
				);
			TTIDenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);
			TTIDenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);
		}
	}
#ifdef TMJ_TIMING
	clock_gettime(CLOCK_REALTIME, &ts1);
	double wall_time_1 = (double)ts1.tv_sec + (double)ts1.tv_nsec * 1e-9;
	double elapsed_time_1 = wall_time_1 - wall_time_0;

	double overcomp, flops_per_cell, gflops, eff_freq, net_comp_throughput;
	TTIDenQ_Compute_Performance((double)thrcnt2*4e-6,elapsed_time_1,bsX,bsY,stage,overcomp,flops_per_cell,gflops,eff_freq,net_comp_throughput);

	printf("elapsed_time = %.2f, cnt2 = %ld, overcomp = %.2f, Net Comp Throughput = %.0f MCells/s, GFLOPS = %.2f, Effective Clock Frequency = %.2f\n", elapsed_time_1, thrcnt2, overcomp, net_comp_throughput, gflops, eff_freq);
#endif

	//if (Buf != 0L) free((void*)Buf);
}
