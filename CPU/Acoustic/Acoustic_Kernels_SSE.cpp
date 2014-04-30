#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <xmmintrin.h>
#include <pmmintrin.h>
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

extern __m128 A2h_A1h_zzz_zzz_1st_11pt;
extern __m128 zzz_A5h_A4h_A3h_1st_11pt;
extern __m128 A2h_A1h_A0h_zzz_2nd_11pt;
extern __m128 zzz_A5h_A4h_A3h_2nd_11pt;

extern float A1h_1st_11pt, A2h_1st_11pt, A3h_1st_11pt, A4h_1st_11pt, A5h_1st_11pt;
extern float A0h_2nd_11pt, A1h_2nd_11pt, A2h_2nd_11pt, A3h_2nd_11pt, A4h_2nd_11pt, A5h_2nd_11pt;

extern __m128 A2_A1_zz_zz_1st_11pt;
extern __m128 zz_A5_A4_A3_1st_11pt;
extern __m128 A2_A1_A0_zz_2nd_11pt;
extern __m128 zz_A5_A4_A3_2nd_11pt;

extern float A1_1st_11pt, A2_1st_11pt, A3_1st_11pt, A4_1st_11pt, A5_1st_11pt;
extern float A0_2nd_11pt, A1_2nd_11pt, A2_2nd_11pt, A3_2nd_11pt, A4_2nd_11pt, A5_2nd_11pt;

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

/*

V0		V1		V2		V3
-1		0		1		2
m5 m4 m3 m2	m1 p1 p2 p3	p4 p5 p6 p7	p8 p9 pA pB
m6 m5 m4 m3	m2 m1 p1 p2	p3 p4 p5 p6	p7 p8 p9 pA
m7 m6 m5 m4	m3 m2 m1 p1	p2 p3 p4 p5	p6 p7 p8 p9
m8 m7 m6 m5	m4 m3 m2 m1	p1 p2 p3 p4	p5 p6 p7 p8

*/
inline void ISO_Cmp_DDX(
	__m128* pq,
	__m128& ddx_0123
	)
{
	__m128 v0 = pq[-1];
	__m128 v1 = pq[ 0];
	__m128 v2 = pq[ 1];
	__m128 v3 = pq[ 2];

	__m128 m5 = v0;
	__m128 p5 = _mm_alignr_epi8(v3, v2, 4);
	ddx_0123 = _mm_mul_ps(_mm_sub_ps(p5, m5), _mm_load1_ps(&A5h));
	
	__m128 m1 = v1;
	__m128 p1 = _mm_alignr_epi8(v2, v1, 4);
	ddx_0123 = _mm_add_ps(ddx_0123, _mm_mul_ps(_mm_sub_ps(p1, m1), _mm_load1_ps(&A1h)));

	__m128 p4 = v2;
	__m128 m4 = _mm_alignr_epi8(v1, v0, 4);
	ddx_0123 = _mm_add_ps(ddx_0123, _mm_mul_ps(_mm_sub_ps(p4, m4), _mm_load1_ps(&A4h)));

	__m128 m2 = _mm_alignr_epi8(v1, v0, 12);
	__m128 p2 = _mm_alignr_epi8(v2, v1, 8);
	ddx_0123 = _mm_add_ps(ddx_0123, _mm_mul_ps(_mm_sub_ps(p2, m2), _mm_load1_ps(&A2h)));

	__m128 m3 = _mm_alignr_epi8(v1, v0, 8);
	__m128 p3 = _mm_alignr_epi8(v2, v1, 12);
	ddx_0123 = _mm_add_ps(ddx_0123, _mm_mul_ps(_mm_sub_ps(p3, m3), _mm_load1_ps(&A3h)));
}

/*

v0		v1		v2		v3
-2		-1		0		1
m8 m7 m6 m5 	m4 m3 m2 m1	p1 p2 p3 p4	p5 p6 p7 p8
m9 m8 m7 m6	m5 m4 m3 m2	m1 p1 p2 p3	p4 p5 p6 p7
mA m9 m8 m7	m6 m5 m4 m3	m2 m1 p1 p2	p3 p4 p5 p6
mB mA m9 m8	m7 m6 m5 m4	m3 m2 m1 p1	p2 p3 p4 p5

*/
inline void ISO_Cmp_DDX_1x_stride2(
	__m128* pq,
	__m128& ddx_0123
	)
{
	__m128 v0 = pq[-4];
	__m128 v1 = pq[-2];
	__m128 v2 = pq[ 0];
	__m128 v3 = pq[ 2];

	__m128 p1 = v2;
	__m128 m1 = _mm_alignr_epi8(v2, v1, 12);
	ddx_0123 = _mm_mul_ps(_mm_sub_ps(p1, m1), _mm_load1_ps(&A1h));

	__m128 p5 = v3;
	__m128 m5 = _mm_alignr_epi8(v1, v0, 12);
	ddx_0123 = _mm_add_ps(ddx_0123, _mm_mul_ps(_mm_sub_ps(p5, m5), _mm_load1_ps(&A5h)));

	__m128 m4 = v1;
	__m128 p4 = _mm_alignr_epi8(v3, v2, 12);
	ddx_0123 = _mm_add_ps(ddx_0123, _mm_mul_ps(_mm_sub_ps(p4, m4), _mm_load1_ps(&A4h)));

	__m128 m3 = _mm_alignr_epi8(v2, v1, 4);
	__m128 p3 = _mm_alignr_epi8(v3, v2, 8);
	ddx_0123 = _mm_add_ps(ddx_0123, _mm_mul_ps(_mm_sub_ps(p3, m3), _mm_load1_ps(&A3h)));

	__m128 m2 = _mm_alignr_epi8(v2, v1, 8);
	__m128 p2 = _mm_alignr_epi8(v3, v2, 4);
	ddx_0123 = _mm_add_ps(ddx_0123, _mm_mul_ps(_mm_sub_ps(p2, m2), _mm_load1_ps(&A2h)));
}

inline void Cmp_DDX(
        __m128* pq,                     // points to pq[x,y,z]
        __m128& ddx_01,
        __m128& ddx_23
        )
{
        __m128 v0 = pq[-2];
        __m128 v1 = pq[-1];
        __m128 v2 = pq[ 0];
        __m128 v3 = pq[ 1];
        __m128 v4 = pq[ 2];
        __m128 r4 = _mm_shuffle_ps(v4, v4, 78);

        __m128 A5hA4hA3hA2h = A5h_A4h_A3h_A2h;
        __m128 acc0 = _mm_mul_ps(_mm_sub_ps(r4, v0), _mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 0x50)); // A5h_A4h
        __m128 acc1 = _mm_mul_ps(_mm_sub_ps(r4, v1), _mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 0xA5)); // A4h_A3h
        __m128 acc2 = _mm_mul_ps(_mm_sub_ps(r4, v2), _mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 0xFA)); // A3h_A2h
        __m128 acc3 = _mm_mul_ps(_mm_sub_ps(r4, v3), A2h_A1h);

        __m128 r3 = _mm_shuffle_ps(v3, v3, 78);
        acc0 = _mm_add_ps(acc0, _mm_mul_ps(_mm_sub_ps(r3, v1), _mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 0xFA)));      // A3h_A2h
        acc1 = _mm_add_ps(acc1, _mm_mul_ps(_mm_sub_ps(r3, v2), A2h_A1h));
        __m128 v5 = pq[ 3];
        __m128 r5 = _mm_shuffle_ps(v5, v5, 78);
        acc2 = _mm_add_ps(acc2, _mm_mul_ps(_mm_sub_ps(r5, v1), _mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 0x50)));      // A5h_A4h
        acc3 = _mm_add_ps(acc3, _mm_mul_ps(_mm_sub_ps(r5, v2), _mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 0xA5)));      // A4h_A3h

        acc0 = _mm_add_ps(acc0, _mm_mul_ps(v2, mA1h_A1h));
        acc1 = _mm_add_ps(acc1, _mm_mul_ps(_mm_shuffle_ps(v5, v0, 228), A5h_mA5h));
        acc2 = _mm_add_ps(acc2, _mm_mul_ps(v3, mA1h_A1h));
        __m128 v6 = pq[ 4];
        acc3 = _mm_add_ps(acc3, _mm_mul_ps(_mm_shuffle_ps(v6, v1, 228), A5h_mA5h));

        ddx_01 = _mm_hadd_ps(_mm_shuffle_ps(acc0, acc0, 216), _mm_shuffle_ps(acc1, acc1, 216));
        ddx_23 = _mm_hadd_ps(_mm_shuffle_ps(acc2, acc2, 216), _mm_shuffle_ps(acc3, acc3, 216));
}

// Verified 06/20/13 TMJ
inline void Cmp_11pt_DDX(
	__m128* pq,
	__m128& ddx_01,
	__m128& ddx_23
	)
{
	//  v0         v1         v2         v3         v4         v5         v6         v7
	// -06 -05    -04 -03    -02 -01    +00 +01    +02 +03    +04 +05    +06 +07    +08 +09
	//      m5     m4  m3     m2  m1     zz  p1     p2  p3     p4  p5
        //             m5  m4     m3  m2     m1  zz     p1  p2     p3  p4     p5
	//                 m5     m4  m3     m2  m1     zz  p1     p2  p3     p4  p5
	//                        m5  m4     m3  m2     m1  zz     p1  p2     p3  p4     p5

	__m128 v7 = pq[ 4];
	__m128 A5_zz_1st = _mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 0x05);
	__m128 acc3_1st = _mm_mul_ps(v7, A5_zz_1st);

	__m128 v6 = pq[ 3];
	__m128 acc1_1st = _mm_mul_ps(v6, A5_zz_1st);
	__m128 A4_A5_1st = _mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 0x5A);
	__m128 acc2_1st = _mm_mul_ps(v6, A4_A5_1st);
	__m128 A3_A4_1st = _mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 0xAF);
	acc3_1st = _mm_add_ps(acc3_1st, _mm_mul_ps(v6, A3_A4_1st));

	__m128 v5 = pq[ 2];
	__m128 acc0_1st = _mm_mul_ps(v5, A4_A5_1st);
	acc1_1st = _mm_add_ps(acc1_1st, _mm_mul_ps(v5, A3_A4_1st));
	__m128 A2_A3_1st = _mm_shuffle_ps(A2h_A1h_zzz_zzz_1st_11pt ,zzz_A5h_A4h_A3h_1st_11pt, 0xF0);
	acc2_1st = _mm_add_ps(acc2_1st, _mm_mul_ps(v5, A2_A3_1st));
	__m128 A1_A2_1st = _mm_shuffle_ps(A2h_A1h_zzz_zzz_1st_11pt, A2h_A1h_zzz_zzz_1st_11pt, 0x05);
	acc3_1st = _mm_add_ps(acc3_1st, _mm_mul_ps(v5, A1_A2_1st));

	__m128 v4 = pq[ 1];
	acc0_1st = _mm_add_ps(acc0_1st, _mm_mul_ps(v4, A2_A3_1st));
	acc1_1st = _mm_add_ps(acc1_1st, _mm_mul_ps(v4, A1_A2_1st));
	__m128 zz_A1_1st = _mm_shuffle_ps(A2h_A1h_zzz_zzz_1st_11pt, A2h_A1h_zzz_zzz_1st_11pt, 0x5A);
	acc2_1st = _mm_add_ps(acc2_1st, _mm_mul_ps(v4, zz_A1_1st));
	__m128 A1_zz_1st = _mm_shuffle_ps(A2h_A1h_zzz_zzz_1st_11pt, A2h_A1h_zzz_zzz_1st_11pt, 0xA5);
	acc3_1st = _mm_sub_ps(acc3_1st, _mm_mul_ps(v4, A1_zz_1st));

	__m128 v3 = pq[ 0];
	acc0_1st = _mm_add_ps(acc0_1st, _mm_mul_ps(v3, zz_A1_1st));
	acc1_1st = _mm_sub_ps(acc1_1st, _mm_mul_ps(v3, A1_zz_1st));
	__m128 A2_A1_1st = _mm_shuffle_ps(A2h_A1h_zzz_zzz_1st_11pt, A2h_A1h_zzz_zzz_1st_11pt, 0x50);
	acc2_1st = _mm_sub_ps(acc2_1st, _mm_mul_ps(v3, A2_A1_1st));
	__m128 A3_A2_1st = _mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, A2h_A1h_zzz_zzz_1st_11pt, 0x0F);
	acc3_1st = _mm_sub_ps(acc3_1st, _mm_mul_ps(v3, A3_A2_1st));

	__m128 v2 = pq[-1];
	acc0_1st = _mm_sub_ps(acc0_1st, _mm_mul_ps(v2, A2_A1_1st));
	acc1_1st = _mm_sub_ps(acc1_1st, _mm_mul_ps(v2, A3_A2_1st));
	__m128 A4_A3_1st = _mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 0xFA);
	acc2_1st = _mm_sub_ps(acc2_1st, _mm_mul_ps(v2, A4_A3_1st));
	__m128 A5_A4_1st = _mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 0xA5);
	acc3_1st = _mm_sub_ps(acc3_1st, _mm_mul_ps(v2, A5_A4_1st));

	__m128 v1 = pq[-2];
	acc0_1st = _mm_sub_ps(acc0_1st, _mm_mul_ps(v1, A4_A3_1st));
	acc1_1st = _mm_sub_ps(acc1_1st, _mm_mul_ps(v1, A5_A4_1st));
	__m128 zz_A5_1st = _mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 0x50);
	acc2_1st = _mm_sub_ps(acc2_1st, _mm_mul_ps(v1, zz_A5_1st));

	__m128 v0 = pq[-3];
	acc0_1st = _mm_sub_ps(acc0_1st, _mm_mul_ps(v0, zz_A5_1st));

        ddx_01 = _mm_hadd_ps(_mm_shuffle_ps(acc0_1st, acc0_1st, 216), _mm_shuffle_ps(acc1_1st, acc1_1st, 216));
        ddx_23 = _mm_hadd_ps(_mm_shuffle_ps(acc2_1st, acc2_1st, 216), _mm_shuffle_ps(acc3_1st, acc3_1st, 216));
}

// Verified 06/20/13 TMJ
inline void Cmp_11pt_DDX_DDX2(
	__m128* pq,
	__m128& ddx_01,
	__m128& ddx_23,
	__m128& ddx2_01,
	__m128& ddx2_23
	)
{
	//  v0         v1         v2         v3         v4         v5         v6         v7
	// -06 -05    -04 -03    -02 -01    +00 +01    +02 +03    +04 +05    +06 +07    +08 +09
	//      m5     m4  m3     m2  m1     zz  p1     p2  p3     p4  p5
        //             m5  m4     m3  m2     m1  zz     p1  p2     p3  p4     p5
	//                 m5     m4  m3     m2  m1     zz  p1     p2  p3     p4  p5
	//                        m5  m4     m3  m2     m1  zz     p1  p2     p3  p4     p5

	__m128 v7 = pq[ 4];
	__m128 A5_zz_1st = _mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 0x05);
	__m128 acc3_1st = _mm_mul_ps(v7, A5_zz_1st);
	__m128 A5_zz_2nd = _mm_shuffle_ps(zzz_A5h_A4h_A3h_2nd_11pt, zzz_A5h_A4h_A3h_2nd_11pt, 0x05);
	__m128 acc3_2nd = _mm_mul_ps(v7, A5_zz_2nd);

	__m128 v6 = pq[ 3];
	__m128 acc1_1st = _mm_mul_ps(v6, A5_zz_1st);
	__m128 A4_A5_1st = _mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 0x5A);
	__m128 acc2_1st = _mm_mul_ps(v6, A4_A5_1st);
	__m128 A3_A4_1st = _mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 0xAF);
	acc3_1st = _mm_add_ps(acc3_1st, _mm_mul_ps(v6, A3_A4_1st));
	__m128 acc1_2nd = _mm_mul_ps(v6, A5_zz_2nd);
	__m128 A4_A5_2nd = _mm_shuffle_ps(zzz_A5h_A4h_A3h_2nd_11pt, zzz_A5h_A4h_A3h_2nd_11pt, 0x5A);
	__m128 acc2_2nd = _mm_mul_ps(v6, A4_A5_2nd);
	__m128 A3_A4_2nd = _mm_shuffle_ps(zzz_A5h_A4h_A3h_2nd_11pt, zzz_A5h_A4h_A3h_2nd_11pt, 0xAF);
	acc3_2nd = _mm_add_ps(acc3_2nd, _mm_mul_ps(v6, A3_A4_2nd));

	__m128 v5 = pq[ 2];
	__m128 acc0_1st = _mm_mul_ps(v5, A4_A5_1st);
	acc1_1st = _mm_add_ps(acc1_1st, _mm_mul_ps(v5, A3_A4_1st));
	__m128 A2_A3_1st = _mm_shuffle_ps(A2h_A1h_zzz_zzz_1st_11pt ,zzz_A5h_A4h_A3h_1st_11pt, 0xF0);
	acc2_1st = _mm_add_ps(acc2_1st, _mm_mul_ps(v5, A2_A3_1st));
	__m128 A1_A2_1st = _mm_shuffle_ps(A2h_A1h_zzz_zzz_1st_11pt, A2h_A1h_zzz_zzz_1st_11pt, 0x05);
	acc3_1st = _mm_add_ps(acc3_1st, _mm_mul_ps(v5, A1_A2_1st));
	__m128 acc0_2nd = _mm_mul_ps(v5, A4_A5_2nd);
	acc1_2nd = _mm_add_ps(acc1_2nd, _mm_mul_ps(v5, A3_A4_2nd));
	__m128 A2_A3_2nd = _mm_shuffle_ps(A2h_A1h_A0h_zzz_2nd_11pt ,zzz_A5h_A4h_A3h_2nd_11pt, 0xF0);
	acc2_2nd = _mm_add_ps(acc2_2nd, _mm_mul_ps(v5, A2_A3_2nd));
	__m128 A1_A2_2nd = _mm_shuffle_ps(A2h_A1h_A0h_zzz_2nd_11pt, A2h_A1h_A0h_zzz_2nd_11pt, 0x05);
	acc3_2nd = _mm_add_ps(acc3_2nd, _mm_mul_ps(v5, A1_A2_2nd));

	__m128 v4 = pq[ 1];
	acc0_1st = _mm_add_ps(acc0_1st, _mm_mul_ps(v4, A2_A3_1st));
	acc1_1st = _mm_add_ps(acc1_1st, _mm_mul_ps(v4, A1_A2_1st));
	__m128 zz_A1_1st = _mm_shuffle_ps(A2h_A1h_zzz_zzz_1st_11pt, A2h_A1h_zzz_zzz_1st_11pt, 0x5A);
	acc2_1st = _mm_add_ps(acc2_1st, _mm_mul_ps(v4, zz_A1_1st));
	__m128 A1_zz_1st = _mm_shuffle_ps(A2h_A1h_zzz_zzz_1st_11pt, A2h_A1h_zzz_zzz_1st_11pt, 0xA5);
	acc3_1st = _mm_sub_ps(acc3_1st, _mm_mul_ps(v4, A1_zz_1st));
	acc0_2nd = _mm_add_ps(acc0_2nd, _mm_mul_ps(v4, A2_A3_2nd));
	acc1_2nd = _mm_add_ps(acc1_2nd, _mm_mul_ps(v4, A1_A2_2nd));
	__m128 zz_A1_2nd = _mm_shuffle_ps(A2h_A1h_A0h_zzz_2nd_11pt, A2h_A1h_A0h_zzz_2nd_11pt, 0x5A);
	acc2_2nd = _mm_add_ps(acc2_2nd, _mm_mul_ps(v4, zz_A1_2nd));
	__m128 A1_zz_2nd = _mm_shuffle_ps(A2h_A1h_A0h_zzz_2nd_11pt, A2h_A1h_A0h_zzz_2nd_11pt, 0xA5);
	acc3_2nd = _mm_add_ps(acc3_2nd, _mm_mul_ps(v4, A1_zz_2nd));

	__m128 v3 = pq[ 0];
	acc0_1st = _mm_add_ps(acc0_1st, _mm_mul_ps(v3, zz_A1_1st));
	acc1_1st = _mm_sub_ps(acc1_1st, _mm_mul_ps(v3, A1_zz_1st));
	__m128 A2_A1_1st = _mm_shuffle_ps(A2h_A1h_zzz_zzz_1st_11pt, A2h_A1h_zzz_zzz_1st_11pt, 0x50);
	acc2_1st = _mm_sub_ps(acc2_1st, _mm_mul_ps(v3, A2_A1_1st));
	__m128 A3_A2_1st = _mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, A2h_A1h_zzz_zzz_1st_11pt, 0x0F);
	acc3_1st = _mm_sub_ps(acc3_1st, _mm_mul_ps(v3, A3_A2_1st));
	acc0_2nd = _mm_add_ps(acc0_2nd, _mm_mul_ps(v3, zz_A1_2nd));
	acc1_2nd = _mm_add_ps(acc1_2nd, _mm_mul_ps(v3, A1_zz_2nd));
	__m128 A2_A1_2nd = _mm_shuffle_ps(A2h_A1h_A0h_zzz_2nd_11pt, A2h_A1h_A0h_zzz_2nd_11pt, 0x50);
	acc2_2nd = _mm_add_ps(acc2_2nd, _mm_mul_ps(v3, A2_A1_2nd));
	__m128 A3_A2_2nd = _mm_shuffle_ps(zzz_A5h_A4h_A3h_2nd_11pt, A2h_A1h_A0h_zzz_2nd_11pt, 0x0F);
	acc3_2nd = _mm_add_ps(acc3_2nd, _mm_mul_ps(v3, A3_A2_2nd));

	__m128 v2 = pq[-1];
	acc0_1st = _mm_sub_ps(acc0_1st, _mm_mul_ps(v2, A2_A1_1st));
	acc1_1st = _mm_sub_ps(acc1_1st, _mm_mul_ps(v2, A3_A2_1st));
	__m128 A4_A3_1st = _mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 0xFA);
	acc2_1st = _mm_sub_ps(acc2_1st, _mm_mul_ps(v2, A4_A3_1st));
	__m128 A5_A4_1st = _mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 0xA5);
	acc3_1st = _mm_sub_ps(acc3_1st, _mm_mul_ps(v2, A5_A4_1st));
	acc0_2nd = _mm_add_ps(acc0_2nd, _mm_mul_ps(v2, A2_A1_2nd));
	acc1_2nd = _mm_add_ps(acc1_2nd, _mm_mul_ps(v2, A3_A2_2nd));
	__m128 A4_A3_2nd = _mm_shuffle_ps(zzz_A5h_A4h_A3h_2nd_11pt, zzz_A5h_A4h_A3h_2nd_11pt, 0xFA);
	acc2_2nd = _mm_add_ps(acc2_2nd, _mm_mul_ps(v2, A4_A3_2nd));
	__m128 A5_A4_2nd = _mm_shuffle_ps(zzz_A5h_A4h_A3h_2nd_11pt, zzz_A5h_A4h_A3h_2nd_11pt, 0xA5);
	acc3_2nd = _mm_add_ps(acc3_2nd, _mm_mul_ps(v2, A5_A4_2nd));

	__m128 v1 = pq[-2];
	acc0_1st = _mm_sub_ps(acc0_1st, _mm_mul_ps(v1, A4_A3_1st));
	acc1_1st = _mm_sub_ps(acc1_1st, _mm_mul_ps(v1, A5_A4_1st));
	__m128 zz_A5_1st = _mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 0x50);
	acc2_1st = _mm_sub_ps(acc2_1st, _mm_mul_ps(v1, zz_A5_1st));
	acc0_2nd = _mm_add_ps(acc0_2nd, _mm_mul_ps(v1, A4_A3_2nd));
	acc1_2nd = _mm_add_ps(acc1_2nd, _mm_mul_ps(v1, A5_A4_2nd));
	__m128 zz_A5_2nd = _mm_shuffle_ps(zzz_A5h_A4h_A3h_2nd_11pt, zzz_A5h_A4h_A3h_2nd_11pt, 0x50);
	acc2_2nd = _mm_add_ps(acc2_2nd, _mm_mul_ps(v1, zz_A5_2nd));

	__m128 v0 = pq[-3];
	acc0_1st = _mm_sub_ps(acc0_1st, _mm_mul_ps(v0, zz_A5_1st));
	acc0_2nd = _mm_add_ps(acc0_2nd, _mm_mul_ps(v0, zz_A5_2nd));

        ddx_01  = _mm_hadd_ps(_mm_shuffle_ps(acc0_1st, acc0_1st, 216), _mm_shuffle_ps(acc1_1st, acc1_1st, 216));
        ddx_23  = _mm_hadd_ps(_mm_shuffle_ps(acc2_1st, acc2_1st, 216), _mm_shuffle_ps(acc3_1st, acc3_1st, 216));
        ddx2_01 = _mm_hadd_ps(_mm_shuffle_ps(acc0_2nd, acc0_2nd, 216), _mm_shuffle_ps(acc1_2nd, acc1_2nd, 216));
        ddx2_23 = _mm_hadd_ps(_mm_shuffle_ps(acc2_2nd, acc2_2nd, 216), _mm_shuffle_ps(acc3_2nd, acc3_2nd, 216));
}

inline void Cmp_11pt_DDY(
	__m128* pq,
	int stride_y_m128,
	__m128& ddy_01,
	__m128& ddy_23
	)
{
	ddy_01 = _mm_mul_ps(_mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 85), _mm_sub_ps(pq[5*stride_y_m128], pq[-5*stride_y_m128]));
	ddy_01 = _mm_add_ps(ddy_01, _mm_mul_ps(_mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 170), _mm_sub_ps(pq[4*stride_y_m128], pq[-4*stride_y_m128])));
	ddy_01 = _mm_add_ps(ddy_01, _mm_mul_ps(_mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 255), _mm_sub_ps(pq[3*stride_y_m128], pq[-3*stride_y_m128])));
	ddy_01 = _mm_add_ps(ddy_01, _mm_mul_ps(_mm_shuffle_ps(A2h_A1h_zzz_zzz_1st_11pt, A2h_A1h_zzz_zzz_1st_11pt, 0), _mm_sub_ps(pq[2*stride_y_m128], pq[-2*stride_y_m128])));
	ddy_01 = _mm_add_ps(ddy_01, _mm_mul_ps(_mm_shuffle_ps(A2h_A1h_zzz_zzz_1st_11pt, A2h_A1h_zzz_zzz_1st_11pt, 85), _mm_sub_ps(pq[1*stride_y_m128], pq[-stride_y_m128])));

	ddy_23 = _mm_mul_ps(_mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 85), _mm_sub_ps(pq[5*stride_y_m128+1], pq[-5*stride_y_m128+1]));
	ddy_23 = _mm_add_ps(ddy_23, _mm_mul_ps(_mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 170), _mm_sub_ps(pq[4*stride_y_m128+1], pq[-4*stride_y_m128+1])));
	ddy_23 = _mm_add_ps(ddy_23, _mm_mul_ps(_mm_shuffle_ps(zzz_A5h_A4h_A3h_1st_11pt, zzz_A5h_A4h_A3h_1st_11pt, 255), _mm_sub_ps(pq[3*stride_y_m128+1], pq[-3*stride_y_m128+1])));
	ddy_23 = _mm_add_ps(ddy_23, _mm_mul_ps(_mm_shuffle_ps(A2h_A1h_zzz_zzz_1st_11pt, A2h_A1h_zzz_zzz_1st_11pt, 0), _mm_sub_ps(pq[2*stride_y_m128+1], pq[-2*stride_y_m128+1])));
	ddy_23 = _mm_add_ps(ddy_23, _mm_mul_ps(_mm_shuffle_ps(A2h_A1h_zzz_zzz_1st_11pt, A2h_A1h_zzz_zzz_1st_11pt, 85), _mm_sub_ps(pq[stride_y_m128+1], pq[-stride_y_m128+1])));
}

inline void Cmp_11pt_DDY2(
	__m128* pq,
	int stride_y_m128,
	__m128& ddy2_01,
	__m128& ddy2_23
	)
{
	ddy2_01 = _mm_mul_ps(_mm_shuffle_ps(zzz_A5h_A4h_A3h_2nd_11pt, zzz_A5h_A4h_A3h_2nd_11pt, 85), _mm_add_ps(pq[5*stride_y_m128], pq[-5*stride_y_m128]));
	ddy2_01 = _mm_add_ps(ddy2_01, _mm_mul_ps(_mm_shuffle_ps(zzz_A5h_A4h_A3h_2nd_11pt, zzz_A5h_A4h_A3h_2nd_11pt, 170), _mm_add_ps(pq[4*stride_y_m128], pq[-4*stride_y_m128])));
	ddy2_01 = _mm_add_ps(ddy2_01, _mm_mul_ps(_mm_shuffle_ps(zzz_A5h_A4h_A3h_2nd_11pt, zzz_A5h_A4h_A3h_2nd_11pt, 255), _mm_add_ps(pq[3*stride_y_m128], pq[-3*stride_y_m128])));
	ddy2_01 = _mm_add_ps(ddy2_01, _mm_mul_ps(_mm_shuffle_ps(A2h_A1h_A0h_zzz_2nd_11pt, A2h_A1h_A0h_zzz_2nd_11pt, 0), _mm_add_ps(pq[2*stride_y_m128], pq[-2*stride_y_m128])));
	ddy2_01 = _mm_add_ps(ddy2_01, _mm_mul_ps(_mm_shuffle_ps(A2h_A1h_A0h_zzz_2nd_11pt, A2h_A1h_A0h_zzz_2nd_11pt, 85), _mm_add_ps(pq[1*stride_y_m128], pq[-stride_y_m128])));
	ddy2_01 = _mm_add_ps(ddy2_01, _mm_mul_ps(_mm_shuffle_ps(A2h_A1h_A0h_zzz_2nd_11pt, A2h_A1h_A0h_zzz_2nd_11pt, 170), *pq));

	ddy2_23 = _mm_mul_ps(_mm_shuffle_ps(zzz_A5h_A4h_A3h_2nd_11pt, zzz_A5h_A4h_A3h_2nd_11pt, 85), _mm_add_ps(pq[5*stride_y_m128+1], pq[-5*stride_y_m128+1]));
	ddy2_23 = _mm_add_ps(ddy2_23, _mm_mul_ps(_mm_shuffle_ps(zzz_A5h_A4h_A3h_2nd_11pt, zzz_A5h_A4h_A3h_2nd_11pt, 170), _mm_add_ps(pq[4*stride_y_m128+1], pq[-4*stride_y_m128+1])));
	ddy2_23 = _mm_add_ps(ddy2_23, _mm_mul_ps(_mm_shuffle_ps(zzz_A5h_A4h_A3h_2nd_11pt, zzz_A5h_A4h_A3h_2nd_11pt, 255), _mm_add_ps(pq[3*stride_y_m128+1], pq[-3*stride_y_m128+1])));
	ddy2_23 = _mm_add_ps(ddy2_23, _mm_mul_ps(_mm_shuffle_ps(A2h_A1h_A0h_zzz_2nd_11pt, A2h_A1h_A0h_zzz_2nd_11pt, 0), _mm_add_ps(pq[2*stride_y_m128+1], pq[-2*stride_y_m128+1])));
	ddy2_23 = _mm_add_ps(ddy2_23, _mm_mul_ps(_mm_shuffle_ps(A2h_A1h_A0h_zzz_2nd_11pt, A2h_A1h_A0h_zzz_2nd_11pt, 85), _mm_add_ps(pq[stride_y_m128+1], pq[-stride_y_m128+1])));
	ddy2_23 = _mm_add_ps(ddy2_23, _mm_mul_ps(_mm_shuffle_ps(A2h_A1h_A0h_zzz_2nd_11pt, A2h_A1h_A0h_zzz_2nd_11pt, 170), pq[1]));
}

inline void Cmp_11pt_DDZ(
	__m128* pq,
	int stride_z_m128,
	__m128& ddz_01,
	__m128& ddz_23
	)
{
	ddz_01 = _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_1st_11pt, zz_A5_A4_A3_1st_11pt, 85), _mm_sub_ps(pq[5*stride_z_m128], pq[-5*stride_z_m128]));
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_1st_11pt, zz_A5_A4_A3_1st_11pt, 170), _mm_sub_ps(pq[4*stride_z_m128], pq[-4*stride_z_m128])));
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_1st_11pt, zz_A5_A4_A3_1st_11pt, 255), _mm_sub_ps(pq[3*stride_z_m128], pq[-3*stride_z_m128])));
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A2_A1_zz_zz_1st_11pt, A2_A1_zz_zz_1st_11pt, 0), _mm_sub_ps(pq[2*stride_z_m128], pq[-2*stride_z_m128])));
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A2_A1_zz_zz_1st_11pt, A2_A1_zz_zz_1st_11pt, 85), _mm_sub_ps(pq[1*stride_z_m128], pq[-stride_z_m128])));

	ddz_23 = _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_1st_11pt, zz_A5_A4_A3_1st_11pt, 85), _mm_sub_ps(pq[5*stride_z_m128+1], pq[-5*stride_z_m128+1]));
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_1st_11pt, zz_A5_A4_A3_1st_11pt, 170), _mm_sub_ps(pq[4*stride_z_m128+1], pq[-4*stride_z_m128+1])));
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_1st_11pt, zz_A5_A4_A3_1st_11pt, 255), _mm_sub_ps(pq[3*stride_z_m128+1], pq[-3*stride_z_m128+1])));
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A2_A1_zz_zz_1st_11pt, A2_A1_zz_zz_1st_11pt, 0), _mm_sub_ps(pq[2*stride_z_m128+1], pq[-2*stride_z_m128+1])));
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A2_A1_zz_zz_1st_11pt, A2_A1_zz_zz_1st_11pt, 85), _mm_sub_ps(pq[stride_z_m128+1], pq[-stride_z_m128+1])));
}

inline void Cmp_11pt_DDZ_DDZ2(
	__m128* pq,
	int stride_z_m128,
	__m128& ddz_01,
	__m128& ddz_23,
	__m128& ddz2_01,
	__m128& ddz2_23
	)
{
	ddz_01  = _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_1st_11pt, zz_A5_A4_A3_1st_11pt, 85), _mm_sub_ps(pq[5*stride_z_m128], pq[-5*stride_z_m128]));
	ddz2_01 = _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_2nd_11pt, zz_A5_A4_A3_2nd_11pt, 85), _mm_add_ps(pq[5*stride_z_m128], pq[-5*stride_z_m128]));
	ddz_01  = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_1st_11pt, zz_A5_A4_A3_1st_11pt, 170), _mm_sub_ps(pq[4*stride_z_m128], pq[-4*stride_z_m128])));
	ddz2_01 = _mm_add_ps(ddz2_01, _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_2nd_11pt, zz_A5_A4_A3_2nd_11pt, 170), _mm_add_ps(pq[4*stride_z_m128], pq[-4*stride_z_m128])));
	ddz_01  = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_1st_11pt, zz_A5_A4_A3_1st_11pt, 255), _mm_sub_ps(pq[3*stride_z_m128], pq[-3*stride_z_m128])));
	ddz2_01 = _mm_add_ps(ddz2_01, _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_2nd_11pt, zz_A5_A4_A3_2nd_11pt, 255), _mm_add_ps(pq[3*stride_z_m128], pq[-3*stride_z_m128])));
	ddz_01  = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A2_A1_zz_zz_1st_11pt, A2_A1_zz_zz_1st_11pt, 0), _mm_sub_ps(pq[2*stride_z_m128], pq[-2*stride_z_m128])));
	ddz2_01 = _mm_add_ps(ddz2_01, _mm_mul_ps(_mm_shuffle_ps(A2_A1_A0_zz_2nd_11pt, A2_A1_A0_zz_2nd_11pt, 0), _mm_add_ps(pq[2*stride_z_m128], pq[-2*stride_z_m128])));
	ddz_01  = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A2_A1_zz_zz_1st_11pt, A2_A1_zz_zz_1st_11pt, 85), _mm_sub_ps(pq[1*stride_z_m128], pq[-stride_z_m128])));
	ddz2_01 = _mm_add_ps(ddz2_01, _mm_mul_ps(_mm_shuffle_ps(A2_A1_A0_zz_2nd_11pt, A2_A1_A0_zz_2nd_11pt, 85), _mm_add_ps(pq[1*stride_z_m128], pq[-stride_z_m128])));
	ddz2_01 = _mm_add_ps(ddz2_01, _mm_mul_ps(_mm_shuffle_ps(A2_A1_A0_zz_2nd_11pt, A2_A1_A0_zz_2nd_11pt, 170), *pq));

	ddz_23  = _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_1st_11pt, zz_A5_A4_A3_1st_11pt, 85), _mm_sub_ps(pq[5*stride_z_m128+1], pq[-5*stride_z_m128+1]));
	ddz2_23 = _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_2nd_11pt, zz_A5_A4_A3_2nd_11pt, 85), _mm_add_ps(pq[5*stride_z_m128+1], pq[-5*stride_z_m128+1]));
	ddz_23  = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_1st_11pt, zz_A5_A4_A3_1st_11pt, 170), _mm_sub_ps(pq[4*stride_z_m128+1], pq[-4*stride_z_m128+1])));
	ddz2_23 = _mm_add_ps(ddz2_23, _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_2nd_11pt, zz_A5_A4_A3_2nd_11pt, 170), _mm_add_ps(pq[4*stride_z_m128+1], pq[-4*stride_z_m128+1])));
	ddz_23  = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_1st_11pt, zz_A5_A4_A3_1st_11pt, 255), _mm_sub_ps(pq[3*stride_z_m128+1], pq[-3*stride_z_m128+1])));
	ddz2_23 = _mm_add_ps(ddz2_23, _mm_mul_ps(_mm_shuffle_ps(zz_A5_A4_A3_2nd_11pt, zz_A5_A4_A3_2nd_11pt, 255), _mm_add_ps(pq[3*stride_z_m128+1], pq[-3*stride_z_m128+1])));
	ddz_23  = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A2_A1_zz_zz_1st_11pt, A2_A1_zz_zz_1st_11pt, 0), _mm_sub_ps(pq[2*stride_z_m128+1], pq[-2*stride_z_m128+1])));
	ddz2_23 = _mm_add_ps(ddz2_23, _mm_mul_ps(_mm_shuffle_ps(A2_A1_A0_zz_2nd_11pt, A2_A1_A0_zz_2nd_11pt, 0), _mm_add_ps(pq[2*stride_z_m128+1], pq[-2*stride_z_m128+1])));
	ddz_23  = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A2_A1_zz_zz_1st_11pt, A2_A1_zz_zz_1st_11pt, 85), _mm_sub_ps(pq[stride_z_m128+1], pq[-stride_z_m128+1])));
	ddz2_23 = _mm_add_ps(ddz2_23, _mm_mul_ps(_mm_shuffle_ps(A2_A1_A0_zz_2nd_11pt, A2_A1_A0_zz_2nd_11pt, 85), _mm_add_ps(pq[stride_z_m128+1], pq[-stride_z_m128+1])));
	ddz2_23 = _mm_add_ps(ddz2_23, _mm_mul_ps(_mm_shuffle_ps(A2_A1_A0_zz_2nd_11pt, A2_A1_A0_zz_2nd_11pt, 170), pq[1]));
}

inline void Cmp_DDX_1x_stride4(
        __m128* pq,
        __m128& ddx_01,
        __m128& ddx_23
        )
{
        __m128 v0 = pq[ -7];
        __m128 v1 = pq[ -4];
        __m128 v2 = pq[ -3];
        __m128 v3 = pq[  0];
        __m128 v4 = pq[  1];
        __m128 v5 = pq[  4];
        __m128 v6 = pq[  5];
          
        __m128 r4 = _mm_shuffle_ps(v4, v4, 78);
          
        __m128 acc0 = _mm_mul_ps(_mm_sub_ps(r4, v1), A4h_A3h);
        __m128 acc1 = _mm_mul_ps(_mm_sub_ps(r4, v2), A3h_A2h);
        __m128 acc2 = _mm_mul_ps(_mm_sub_ps(r4, v3), A2h_A1h);
        __m128 acc3 = _mm_mul_ps(v4, mA1h_A1h);

        acc0 = _mm_add_ps(acc0, _mm_mul_ps(_mm_shuffle_ps(v5, v0, 228), A5h_mA5h));
        __m128 r5 = _mm_shuffle_ps(v5, v5, 78);
        acc1 = _mm_add_ps(acc1, _mm_mul_ps(_mm_sub_ps(r5, v1), A5h_A4h));
        acc2 = _mm_add_ps(acc2, _mm_mul_ps(_mm_sub_ps(r5, v2), A4h_A3h));
        acc3 = _mm_add_ps(acc3, _mm_mul_ps(_mm_sub_ps(r5, v3), A3h_A2h));

        __m128 r3 = _mm_shuffle_ps(v3, v3, 78); 
        acc0 = _mm_add_ps(acc0, _mm_mul_ps(_mm_sub_ps(r3, v2), A2h_A1h));
        acc1 = _mm_add_ps(acc1, _mm_mul_ps(v3, mA1h_A1h));
        acc2 = _mm_add_ps(acc2, _mm_mul_ps(_mm_shuffle_ps(v6, v1, 228), A5h_mA5h));
        __m128 r6 = _mm_shuffle_ps(v6, v6, 78);
        acc3 = _mm_add_ps(acc3, _mm_mul_ps(_mm_sub_ps(r6, v2), A5h_A4h));

        ddx_01 = _mm_hadd_ps(_mm_shuffle_ps(acc0, acc0, 216), _mm_shuffle_ps(acc1, acc1, 216));
        ddx_23 = _mm_hadd_ps(_mm_shuffle_ps(acc2, acc2, 216), _mm_shuffle_ps(acc3, acc3, 216));
}

inline void Cmp_DDX_1x_stride8(
        __m128* pq,
        __m128& ddx_01,
        __m128& ddx_23
        )
{
        __m128 v0 = pq[-15];
        __m128 v1 = pq[ -8];
        __m128 v2 = pq[ -7];
        __m128 v3 = pq[  0];
        __m128 v4 = pq[  1];
        __m128 v5 = pq[  8];
        __m128 v6 = pq[  9];

        __m128 r4 = _mm_shuffle_ps(v4, v4, 78);

        __m128 acc0 = _mm_mul_ps(_mm_sub_ps(r4, v1), A4h_A3h);
        __m128 acc1 = _mm_mul_ps(_mm_sub_ps(r4, v2), A3h_A2h);
        __m128 acc2 = _mm_mul_ps(_mm_sub_ps(r4, v3), A2h_A1h);
        __m128 acc3 = _mm_mul_ps(v4, mA1h_A1h);

        acc0 = _mm_add_ps(acc0, _mm_mul_ps(_mm_shuffle_ps(v5, v0, 228), A5h_mA5h));
        __m128 r5 = _mm_shuffle_ps(v5, v5, 78);
        acc1 = _mm_add_ps(acc1, _mm_mul_ps(_mm_sub_ps(r5, v1), A5h_A4h));
        acc2 = _mm_add_ps(acc2, _mm_mul_ps(_mm_sub_ps(r5, v2), A4h_A3h));
        acc3 = _mm_add_ps(acc3, _mm_mul_ps(_mm_sub_ps(r5, v3), A3h_A2h));

        __m128 r3 = _mm_shuffle_ps(v3, v3, 78);
        acc0 = _mm_add_ps(acc0, _mm_mul_ps(_mm_sub_ps(r3, v2), A2h_A1h));
        acc1 = _mm_add_ps(acc1, _mm_mul_ps(v3, mA1h_A1h));
        acc2 = _mm_add_ps(acc2, _mm_mul_ps(_mm_shuffle_ps(v6, v1, 228), A5h_mA5h));
        __m128 r6 = _mm_shuffle_ps(v6, v6, 78);
        acc3 = _mm_add_ps(acc3, _mm_mul_ps(_mm_sub_ps(r6, v2), A5h_A4h));

        ddx_01 = _mm_hadd_ps(_mm_shuffle_ps(acc0, acc0, 216), _mm_shuffle_ps(acc1, acc1, 216));
        ddx_23 = _mm_hadd_ps(_mm_shuffle_ps(acc2, acc2, 216), _mm_shuffle_ps(acc3, acc3, 216));
}

inline void ISO_Cmp_DDY(
	__m128* pq,
	int stride_y_m128,
	__m128& ddy_0123
	)
{
	__m128 A5hA4hA3hA2h = A5h_A4h_A3h_A2h;
	
	ddy_0123 = _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h,   0), _mm_sub_ps(pq[5*stride_y_m128], pq[-4*stride_y_m128]));				// A5h
	ddy_0123 = _mm_add_ps(ddy_0123, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h,  85), _mm_sub_ps(pq[4*stride_y_m128], pq[-3*stride_y_m128])));        // A4h
	ddy_0123 = _mm_add_ps(ddy_0123, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 170), _mm_sub_ps(pq[3*stride_y_m128], pq[-2*stride_y_m128])));        // A3h
	ddy_0123 = _mm_add_ps(ddy_0123, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 255), _mm_sub_ps(pq[2*stride_y_m128], pq[-stride_y_m128])));          // A2h
	ddy_0123 = _mm_add_ps(ddy_0123, _mm_mul_ps(A1h_A1h_A1h_A1h, _mm_sub_ps(pq[stride_y_m128], pq[0])));								// A1h
}

inline void Cmp_DDY(
	__m128* pq,
	int stride_y_m128,
	__m128& ddy_01,
	__m128& ddy_23
	)
{
	__m128 A5hA4hA3hA2h = A5h_A4h_A3h_A2h;
	
	/*
	ddy_01 = _mm_mul_ps(A1h_A1h_A1h_A1h, _mm_sub_ps(pq[stride_y_m128], pq[0]));	// A1h
	ddy_01 = _mm_add_ps(ddy_01, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 255), _mm_sub_ps(pq[2*stride_y_m128], pq[-stride_y_m128])));		// A2h
        ddy_01 = _mm_add_ps(ddy_01, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 170), _mm_sub_ps(pq[3*stride_y_m128], pq[-2*stride_y_m128])));	// A3h
        ddy_01 = _mm_add_ps(ddy_01, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h,  85), _mm_sub_ps(pq[4*stride_y_m128], pq[-3*stride_y_m128])));	// A4h
        ddy_01 = _mm_add_ps(ddy_01, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h,   0), _mm_sub_ps(pq[5*stride_y_m128], pq[-4*stride_y_m128])));	// A5h

        ddy_23 = _mm_mul_ps(A1h_A1h_A1h_A1h, _mm_sub_ps(pq[1+stride_y_m128], pq[1]));    // A1h
        ddy_23 = _mm_add_ps(ddy_23, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 255), _mm_sub_ps(pq[1+2*stride_y_m128], pq[1-stride_y_m128]))); 	// A2h
        ddy_23 = _mm_add_ps(ddy_23, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 170), _mm_sub_ps(pq[1+3*stride_y_m128], pq[1-2*stride_y_m128]))); 	// A3h
        ddy_23 = _mm_add_ps(ddy_23, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h,  85), _mm_sub_ps(pq[1+4*stride_y_m128], pq[1-3*stride_y_m128]))); 	// A4h
        ddy_23 = _mm_add_ps(ddy_23, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h,   0), _mm_sub_ps(pq[1+5*stride_y_m128], pq[1-4*stride_y_m128]))); 	// A5h
	*/

	ddy_01 = _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h,   0), _mm_sub_ps(pq[5*stride_y_m128], pq[-4*stride_y_m128]));				// A5h
	ddy_01 = _mm_add_ps(ddy_01, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h,  85), _mm_sub_ps(pq[4*stride_y_m128], pq[-3*stride_y_m128])));        // A4h
	ddy_01 = _mm_add_ps(ddy_01, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 170), _mm_sub_ps(pq[3*stride_y_m128], pq[-2*stride_y_m128])));        // A3h
	ddy_01 = _mm_add_ps(ddy_01, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 255), _mm_sub_ps(pq[2*stride_y_m128], pq[-stride_y_m128])));          // A2h
	ddy_01 = _mm_add_ps(ddy_01, _mm_mul_ps(A1h_A1h_A1h_A1h, _mm_sub_ps(pq[stride_y_m128], pq[0])));								// A1h

	ddy_23 = _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h,   0), _mm_sub_ps(pq[1+5*stride_y_m128], pq[1-4*stride_y_m128]));				// A5h
	ddy_23 = _mm_add_ps(ddy_23, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h,  85), _mm_sub_ps(pq[1+4*stride_y_m128], pq[1-3*stride_y_m128])));     // A4h
	ddy_23 = _mm_add_ps(ddy_23, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 170), _mm_sub_ps(pq[1+3*stride_y_m128], pq[1-2*stride_y_m128])));     // A3h
	ddy_23 = _mm_add_ps(ddy_23, _mm_mul_ps(_mm_shuffle_ps(A5hA4hA3hA2h, A5hA4hA3hA2h, 255), _mm_sub_ps(pq[1+2*stride_y_m128], pq[1-stride_y_m128])));       // A2h
	ddy_23 = _mm_add_ps(ddy_23, _mm_mul_ps(A1h_A1h_A1h_A1h, _mm_sub_ps(pq[1+stride_y_m128], pq[1])));							// A1h
}

inline void ISO_Cmp_DDZ(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_0123
        )
{
        __m128 A5A4A3A2 = A5_A4_A3_A2;

	ddz_0123 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_sub_ps(pq[5*stride_z_m128], pq[-4*stride_z_m128]));				// A5
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_sub_ps(pq[4*stride_z_m128], pq[-3*stride_z_m128])));        // A4
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[3*stride_z_m128], pq[-2*stride_z_m128])));        // A3
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0])));							// A1
}

inline void ISO_Cmp_DDZ_EE(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_0123
        )
{
	__m128 E1E1E1E1 = _mm_set1_ps(E1);
	ddz_0123 = _mm_mul_ps(E1E1E1E1, _mm_sub_ps(pq[stride_z_m128], pq[0]));
}

inline void ISO_Cmp_DDZ_Z0(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_0123
        )
{
        __m128 A5A4A3A2 = _mm_load_ps((float*)&A5_A4_A3_A2);

	ddz_0123 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], pq[4*stride_z_m128]));				// A5
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[4*stride_z_m128], pq[3*stride_z_m128])));        	// A4
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_add_ps(pq[3*stride_z_m128], pq[2*stride_z_m128])));        	// A3
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_add_ps(pq[2*stride_z_m128], pq[stride_z_m128])));          	// A2
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0])));							// A1
}

inline void ISO_Cmp_DDZ_Z1(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_0123
        )
{
        __m128 A5A4A3A2 = _mm_load_ps((float*)&A5_A4_A3_A2);

	ddz_0123 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], pq[2*stride_z_m128]));				// A5
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[4*stride_z_m128], pq[stride_z_m128])));           // A4
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_add_ps(pq[3*stride_z_m128], pq[0])));                       // A3
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0])));							// A1
}

inline void ISO_Cmp_DDZ_Z2(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_0123
        )
{
        __m128 A5A4A3A2 = _mm_load_ps((float*)&A5_A4_A3_A2);

	ddz_0123 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], pq[0]));						// A5
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[4*stride_z_m128], pq[-1*stride_z_m128])));        // A4
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[3*stride_z_m128], pq[-2*stride_z_m128])));        // A3
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0])));							// A1
}

inline void ISO_Cmp_DDZ_Z3(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_0123
        )
{
        __m128 A5A4A3A2 = _mm_load_ps((float*)&A5_A4_A3_A2);

	ddz_0123 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], pq[-2*stride_z_m128]));				// A5
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_sub_ps(pq[4*stride_z_m128], pq[-3*stride_z_m128])));        // A4
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[3*stride_z_m128], pq[-2*stride_z_m128])));        // A3
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0])));							// A1
}

/*
inline void ISO_Cmp_DDZ_Z0(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_0123
        )
{
        __m128 A5A4A3A2 = _mm_load_ps((float*)&A5_A4_A3_A2);

	ddz_0123 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], _mm_mul_ps(_mm_set1_ps(0.8f), pq[4*stride_z_m128])));				// A5
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[4*stride_z_m128], _mm_mul_ps(_mm_set1_ps(0.8f), pq[3*stride_z_m128]))));        	// A4
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_add_ps(pq[3*stride_z_m128], _mm_mul_ps(_mm_set1_ps(0.8f), pq[2*stride_z_m128]))));        	// A3
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_add_ps(pq[2*stride_z_m128], _mm_mul_ps(_mm_set1_ps(0.8f), pq[stride_z_m128]))));          	// A2
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0])));							// A1
}

inline void ISO_Cmp_DDZ_Z1(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_0123
        )
{
        __m128 A5A4A3A2 = _mm_load_ps((float*)&A5_A4_A3_A2);

	ddz_0123 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], _mm_mul_ps(_mm_set1_ps(0.8f), pq[2*stride_z_m128])));				// A5
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[4*stride_z_m128], _mm_mul_ps(_mm_set1_ps(0.8f), pq[stride_z_m128]))));           // A4
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_add_ps(pq[3*stride_z_m128], _mm_mul_ps(_mm_set1_ps(0.8f), pq[0]))));                       // A3
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0])));							// A1
}

inline void ISO_Cmp_DDZ_Z2(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_0123
        )
{
        __m128 A5A4A3A2 = _mm_load_ps((float*)&A5_A4_A3_A2);

	ddz_0123 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], _mm_mul_ps(_mm_set1_ps(0.8f), pq[0])));						// A5
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[4*stride_z_m128], _mm_mul_ps(_mm_set1_ps(0.8f), pq[-1*stride_z_m128]))));        // A4
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[3*stride_z_m128], pq[-2*stride_z_m128])));        // A3
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0])));							// A1
}

inline void ISO_Cmp_DDZ_Z3(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_0123
        )
{
        __m128 A5A4A3A2 = _mm_load_ps((float*)&A5_A4_A3_A2);

	ddz_0123 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], _mm_mul_ps(_mm_set1_ps(0.8f), pq[-2*stride_z_m128])));				// A5
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_sub_ps(pq[4*stride_z_m128], pq[-3*stride_z_m128])));        // A4
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[3*stride_z_m128], pq[-2*stride_z_m128])));        // A3
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
	ddz_0123 = _mm_add_ps(ddz_0123, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0])));							// A1
}
*/

inline void Cmp_DDZ(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_01,
        __m128& ddz_23
        )
{
        __m128 A5A4A3A2 = A5_A4_A3_A2;

	/*
        ddz_01 = _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0]));     // A1
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[3*stride_z_m128], pq[-2*stride_z_m128])));        // A3
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_sub_ps(pq[4*stride_z_m128], pq[-3*stride_z_m128])));        // A4
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_sub_ps(pq[5*stride_z_m128], pq[-4*stride_z_m128])));        // A5

        ddz_23 = _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[1+stride_z_m128], pq[1]));    // A1
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[1+2*stride_z_m128], pq[1-stride_z_m128])));       // A2
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[1+3*stride_z_m128], pq[1-2*stride_z_m128])));     // A3
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_sub_ps(pq[1+4*stride_z_m128], pq[1-3*stride_z_m128])));     // A4
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_sub_ps(pq[1+5*stride_z_m128], pq[1-4*stride_z_m128])));     // A5  
	*/

	ddz_01 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_sub_ps(pq[5*stride_z_m128], pq[-4*stride_z_m128]));				// A5
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_sub_ps(pq[4*stride_z_m128], pq[-3*stride_z_m128])));        // A4
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[3*stride_z_m128], pq[-2*stride_z_m128])));        // A3
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0])));							// A1
	
	ddz_23 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_sub_ps(pq[1+5*stride_z_m128], pq[1-4*stride_z_m128]));				// A5  
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_sub_ps(pq[1+4*stride_z_m128], pq[1-3*stride_z_m128])));     // A4
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[1+3*stride_z_m128], pq[1-2*stride_z_m128])));     // A3
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[1+2*stride_z_m128], pq[1-stride_z_m128])));       // A2
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[1+stride_z_m128], pq[1])));							// A1
}

inline void Cmp_DDZ_2Z(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_z0_01,
        __m128& ddz_z0_23,
        __m128& ddz_z1_01,
        __m128& ddz_z1_23
        )
{
        __m128 A5A4A3A2 = A5_A4_A3_A2;

	__m128 mA5 = _mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0);
	ddz_z1_01 = _mm_mul_ps(mA5, _mm_sub_ps(pq[  6*stride_z_m128], pq[ -3*stride_z_m128]));				// A5
	ddz_z1_23 = _mm_mul_ps(mA5, _mm_sub_ps(pq[1+6*stride_z_m128], pq[1-3*stride_z_m128]));				// A5  
	ddz_z0_01 = _mm_mul_ps(mA5, _mm_sub_ps(pq[  5*stride_z_m128], pq[ -4*stride_z_m128]));				// A5
	ddz_z0_23 = _mm_mul_ps(mA5, _mm_sub_ps(pq[1+5*stride_z_m128], pq[1-4*stride_z_m128]));				// A5  

	__m128 mA4 = _mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85);
	ddz_z1_01 = _mm_add_ps(ddz_z1_01, _mm_mul_ps(mA4, _mm_sub_ps(pq[  5*stride_z_m128], pq[ -2*stride_z_m128])));        // A4
	ddz_z1_23 = _mm_add_ps(ddz_z1_23, _mm_mul_ps(mA4, _mm_sub_ps(pq[1+5*stride_z_m128], pq[1-2*stride_z_m128])));     // A4
	ddz_z0_01 = _mm_add_ps(ddz_z0_01, _mm_mul_ps(mA4, _mm_sub_ps(pq[  4*stride_z_m128], pq[ -3*stride_z_m128])));        // A4
	ddz_z0_23 = _mm_add_ps(ddz_z0_23, _mm_mul_ps(mA4, _mm_sub_ps(pq[1+4*stride_z_m128], pq[1-3*stride_z_m128])));     // A4

	__m128 mA3 = _mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170);
	ddz_z1_01 = _mm_add_ps(ddz_z1_01, _mm_mul_ps(mA3, _mm_sub_ps(pq[  4*stride_z_m128], pq[ -1*stride_z_m128])));        // A3
	ddz_z1_23 = _mm_add_ps(ddz_z1_23, _mm_mul_ps(mA3, _mm_sub_ps(pq[1+4*stride_z_m128], pq[1-1*stride_z_m128])));     // A3
	ddz_z0_01 = _mm_add_ps(ddz_z0_01, _mm_mul_ps(mA3, _mm_sub_ps(pq[  3*stride_z_m128], pq[ -2*stride_z_m128])));        // A3
	ddz_z0_23 = _mm_add_ps(ddz_z0_23, _mm_mul_ps(mA3, _mm_sub_ps(pq[1+3*stride_z_m128], pq[1-2*stride_z_m128])));     // A3

	__m128 mA2 = _mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255);
	ddz_z1_01 = _mm_add_ps(ddz_z1_01, _mm_mul_ps(mA2, _mm_sub_ps(pq[  3*stride_z_m128], pq[0])));          // A2
	ddz_z1_23 = _mm_add_ps(ddz_z1_23, _mm_mul_ps(mA2, _mm_sub_ps(pq[1+3*stride_z_m128], pq[1])));       // A2
	ddz_z0_01 = _mm_add_ps(ddz_z0_01, _mm_mul_ps(mA2, _mm_sub_ps(pq[  2*stride_z_m128], pq[-stride_z_m128])));          // A2
	ddz_z0_23 = _mm_add_ps(ddz_z0_23, _mm_mul_ps(mA2, _mm_sub_ps(pq[1+2*stride_z_m128], pq[1-stride_z_m128])));       // A2

	ddz_z1_01 = _mm_add_ps(ddz_z1_01, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[  2*stride_z_m128], pq[stride_z_m128])));							// A1
	ddz_z1_23 = _mm_add_ps(ddz_z1_23, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[1+2*stride_z_m128], pq[1+stride_z_m128])));							// A1
	ddz_z0_01 = _mm_add_ps(ddz_z0_01, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[    stride_z_m128], pq[0])));							// A1
	ddz_z0_23 = _mm_add_ps(ddz_z0_23, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[1+  stride_z_m128], pq[1])));							// A1
}

inline void Cmp_DDZ_EE(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_01,
        __m128& ddz_23
        )
{
	__m128 E1E1E1E1 = _mm_set1_ps(E1);
	ddz_01 = _mm_mul_ps(E1E1E1E1, _mm_sub_ps(pq[stride_z_m128], pq[0]));
	ddz_23 = _mm_mul_ps(E1E1E1E1, _mm_sub_ps(pq[1+stride_z_m128], pq[1]));
}

inline void Cmp_DDZ_Z0(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_01,
        __m128& ddz_23
        )
{
        __m128 A5A4A3A2 = _mm_load_ps((float*)&A5_A4_A3_A2);

	/*
        ddz_01 = _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0]));     // A1
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_add_ps(pq[2*stride_z_m128], pq[stride_z_m128])));          // A2
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_add_ps(pq[3*stride_z_m128], pq[2*stride_z_m128])));        // A3
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[4*stride_z_m128], pq[3*stride_z_m128])));        // A4
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], pq[4*stride_z_m128])));        // A5

        ddz_23 = _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[1+stride_z_m128], pq[1]));    // A1
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_add_ps(pq[1+2*stride_z_m128], pq[1+stride_z_m128])));       // A2
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_add_ps(pq[1+3*stride_z_m128], pq[1+2*stride_z_m128])));     // A3
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[1+4*stride_z_m128], pq[1+3*stride_z_m128])));     // A4
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[1+5*stride_z_m128], pq[1+4*stride_z_m128])));     // A5  
	*/

	ddz_01 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], pq[4*stride_z_m128]));				// A5
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[4*stride_z_m128], pq[3*stride_z_m128])));        	// A4
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_add_ps(pq[3*stride_z_m128], pq[2*stride_z_m128])));        	// A3
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_add_ps(pq[2*stride_z_m128], pq[stride_z_m128])));          	// A2
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0])));							// A1

	ddz_23 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[1+5*stride_z_m128], pq[1+4*stride_z_m128]));				// A5  
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[1+4*stride_z_m128], pq[1+3*stride_z_m128])));     // A4
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_add_ps(pq[1+3*stride_z_m128], pq[1+2*stride_z_m128])));     // A3
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_add_ps(pq[1+2*stride_z_m128], pq[1+stride_z_m128])));       // A2
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[1+stride_z_m128], pq[1])));							// A1
}

inline void Cmp_DDZ_Z1(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_01,
        __m128& ddz_23
        )
{
        __m128 A5A4A3A2 = _mm_load_ps((float*)&A5_A4_A3_A2);

	/*
        ddz_01 = _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0]));     									// A1
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_add_ps(pq[3*stride_z_m128], pq[0])));        		// A3
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[4*stride_z_m128], pq[stride_z_m128])));        	// A4
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], pq[2*stride_z_m128])));        	// A5

        ddz_23 = _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[1+stride_z_m128], pq[1]));    									// A1
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[1+2*stride_z_m128], pq[1-stride_z_m128])));       // A2
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_add_ps(pq[1+3*stride_z_m128], pq[1])));     		// A3
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[1+4*stride_z_m128], pq[1+stride_z_m128])));     	// A4
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[1+5*stride_z_m128], pq[1+2*stride_z_m128])));     // A5  
	*/

	ddz_01 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], pq[2*stride_z_m128]));				// A5
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[4*stride_z_m128], pq[stride_z_m128])));           // A4
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_add_ps(pq[3*stride_z_m128], pq[0])));                       // A3
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0])));							// A1
	
	ddz_23 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[1+5*stride_z_m128], pq[1+2*stride_z_m128]));				// A5  
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[1+4*stride_z_m128], pq[1+stride_z_m128])));       // A4
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_add_ps(pq[1+3*stride_z_m128], pq[1])));                     // A3
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[1+2*stride_z_m128], pq[1-stride_z_m128])));       // A2
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[1+stride_z_m128], pq[1])));							// A1
}

inline void Cmp_DDZ_Z2(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_01,
        __m128& ddz_23
        )
{
        __m128 A5A4A3A2 = _mm_load_ps((float*)&A5_A4_A3_A2);

	/*
        ddz_01 = _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0]));     									// A1
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[3*stride_z_m128], pq[-2*stride_z_m128])));        // A3
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[4*stride_z_m128], pq[-1*stride_z_m128])));        // A4
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], pq[0])));        		// A5

        ddz_23 = _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[1+stride_z_m128], pq[1]));    									// A1
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[1+2*stride_z_m128], pq[1-stride_z_m128])));       // A2
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[1+3*stride_z_m128], pq[1-2*stride_z_m128])));     // A3
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[1+4*stride_z_m128], pq[1-1*stride_z_m128])));     // A4
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[1+5*stride_z_m128], pq[1])));     		// A5  
	*/

	ddz_01 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], pq[0]));						// A5
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[4*stride_z_m128], pq[-1*stride_z_m128])));        // A4
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[3*stride_z_m128], pq[-2*stride_z_m128])));        // A3
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0])));							// A1
	
	ddz_23 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[1+5*stride_z_m128], pq[1]));						// A5  
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_add_ps(pq[1+4*stride_z_m128], pq[1-1*stride_z_m128])));     // A4
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[1+3*stride_z_m128], pq[1-2*stride_z_m128])));     // A3
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[1+2*stride_z_m128], pq[1-stride_z_m128])));       // A2
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[1+stride_z_m128], pq[1])));							// A1
}

inline void Cmp_DDZ_Z3(
        __m128* pq,
        int stride_z_m128,
        __m128& ddz_01,
        __m128& ddz_23
        )
{
        __m128 A5A4A3A2 = _mm_load_ps((float*)&A5_A4_A3_A2);

	/*
        ddz_01 = _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0]));     									// A1
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[3*stride_z_m128], pq[-2*stride_z_m128])));        // A3
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_sub_ps(pq[4*stride_z_m128], pq[-3*stride_z_m128])));        // A4
        ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], pq[-2*stride_z_m128])));        // A5

        ddz_23 = _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[1+stride_z_m128], pq[1]));    									// A1
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[1+2*stride_z_m128], pq[1-stride_z_m128])));       // A2
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[1+3*stride_z_m128], pq[1-2*stride_z_m128])));     // A3
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_sub_ps(pq[1+4*stride_z_m128], pq[1-3*stride_z_m128])));     // A4
        ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[1+5*stride_z_m128], pq[1-2*stride_z_m128])));     // A5  
	*/

	ddz_01 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[5*stride_z_m128], pq[-2*stride_z_m128]));				// A5
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_sub_ps(pq[4*stride_z_m128], pq[-3*stride_z_m128])));        // A4
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[3*stride_z_m128], pq[-2*stride_z_m128])));        // A3
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[2*stride_z_m128], pq[-stride_z_m128])));          // A2
	ddz_01 = _mm_add_ps(ddz_01, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[stride_z_m128], pq[0])));							// A1
	
	ddz_23 = _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,   0), _mm_add_ps(pq[1+5*stride_z_m128], pq[1-2*stride_z_m128]));				// A5  
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2,  85), _mm_sub_ps(pq[1+4*stride_z_m128], pq[1-3*stride_z_m128])));     // A4
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 170), _mm_sub_ps(pq[1+3*stride_z_m128], pq[1-2*stride_z_m128])));     // A3
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(_mm_shuffle_ps(A5A4A3A2, A5A4A3A2, 255), _mm_sub_ps(pq[1+2*stride_z_m128], pq[1-stride_z_m128])));       // A2
	ddz_23 = _mm_add_ps(ddz_23, _mm_mul_ps(A1_A1_A1_A1, _mm_sub_ps(pq[1+stride_z_m128], pq[1])));							// A1
}

inline void ISODenQ_Cmp_DXED(
	__m128& ddx_0123,
	__m128& Boy,
	__m128& dxed1_0123
	)
{
	dxed1_0123 = _mm_mul_ps(ddx_0123, Boy);
}

inline void ISODenQ_Cmp_DYED(
	__m128& ddy_0123,
	__m128& Boy,
	__m128& dyed1_0123
	)
{
	dyed1_0123 = _mm_mul_ps(ddy_0123, Boy);
}

inline void ISODenQ_Cmp_DXED_DYED_DZED(
	__m128& ddx_0123,
	__m128& ddy_0123,
	__m128& ddz_0123,
	__m128& Boy,
	__m128& dxed1_0123,
	__m128& dyed1_0123,
	__m128& dzed1_0123
	)
{
	dxed1_0123 = _mm_mul_ps(ddx_0123, Boy);
	dyed1_0123 = _mm_mul_ps(ddy_0123, Boy);
	dzed1_0123 = _mm_mul_ps(ddz_0123, Boy);
}

inline void VTIDenQ_Cmp_DXED(
	__m128& ddx_01,
	__m128& ddx_23,
	__m128& Boy,
        __m128& dxed1_01,
        __m128& dxed1_23
	)
{
	/*
	temp = cAzm * ddx + sAzm * ddy
	temp2 = buoy * ( cAzm * ddy - sAzm * ddx )
	temp3 = buoy * ( cDip * temp - sDip * ddz )
	temp4 = buoy * ( sdip * temp + cDip * ddz )

	dxed1 = -sAzm * temp2 + cDip * cAzm * temp3
	dyed1 = cAzm * temp2 + cDip * sAzm * temp3
	dzed1 = -sDip * temp3

	dxed2 = sDip * cAzm * temp4
	dyed2 = sDip * sAzm * temp4
	dzed2 = cDip * temp4

	VTI:

	cDip = cAzm = 1
	sDip = sAzm = 0
	
	temp  = ddx
	temp2 = buoy * ddy
	temp3 = buoy * ddx
	temp4 = buoy * ddz

	dxed1 = buoy * ddx
	dyed1 = buoy * ddy
	dzed2 = buoy * ddz

	=> VTI requires 3 stencils instead of 9 : d_dx2, d_dy2 and d_dz2.


	****

	temp = ddx
	temp2 = buoy * ddy
	temp3 = buoy * temp = buoy * ddx
	temp4 = buoy * ddz

	dxed1 = temp3 = buoy * ddx
	dyed1 = temp2 = buoy * ddy
	dzed1 = 0

	V4 = d(buoy * ddx)dx + d(buoy * ddy)dy
	V5 = d(buoy * ddz)dz

	dxed2 = 0
	dyed2 = 0
	dzed2 = temp4 = buoy * ddz

	*/

	//ddx_01 = ddx_23 = _mm_set1_ps(1.0f);
	//ddy_01 = ddy_23 = _mm_set1_ps(2.0f);
	//ddz_01 = ddz_23 = _mm_set1_ps(3.0f);

        __m128 Boy_01 = _mm_shuffle_ps(Boy, Boy, 0x50);
	dxed1_01 = _mm_mul_ps(Boy_01, ddx_01);

	// 23

        __m128 Boy_23 = _mm_shuffle_ps(Boy, Boy, 0xFA);
	dxed1_23 = _mm_mul_ps(Boy_23, ddx_23);
}

inline void VTIDenQ_Cmp_DYED(
	__m128& ddy_01,
	__m128& ddy_23,
	__m128& Boy,
        __m128& dyed1_01,
        __m128& dyed1_23
	)
{
	/*
	temp = cAzm * ddx + sAzm * ddy
	temp2 = buoy * ( cAzm * ddy - sAzm * ddx )
	temp3 = buoy * ( cDip * temp - sDip * ddz )
	temp4 = buoy * ( sdip * temp + cDip * ddz )

	dxed1 = -sAzm * temp2 + cDip * cAzm * temp3
	dyed1 = cAzm * temp2 + cDip * sAzm * temp3
	dzed1 = -sDip * temp3

	dxed2 = sDip * cAzm * temp4
	dyed2 = sDip * sAzm * temp4
	dzed2 = cDip * temp4

	VTI:

	cDip = cAzm = 1
	sDip = sAzm = 0
	
	temp  = ddx
	temp2 = buoy * ddy
	temp3 = buoy * ddx
	temp4 = buoy * ddz

	dxed1 = buoy * ddx
	dyed1 = buoy * ddy
	dzed2 = buoy * ddz

	=> VTI requires 3 stencils instead of 9 : d_dx2, d_dy2 and d_dz2.


	****

	temp = ddx
	temp2 = buoy * ddy
	temp3 = buoy * temp = buoy * ddx
	temp4 = buoy * ddz

	dxed1 = temp3 = buoy * ddx
	dyed1 = temp2 = buoy * ddy
	dzed1 = 0

	V4 = d(buoy * ddx)dx + d(buoy * ddy)dy
	V5 = d(buoy * ddz)dz

	dxed2 = 0
	dyed2 = 0
	dzed2 = temp4 = buoy * ddz

	*/

	//ddx_01 = ddx_23 = _mm_set1_ps(1.0f);
	//ddy_01 = ddy_23 = _mm_set1_ps(2.0f);
	//ddz_01 = ddz_23 = _mm_set1_ps(3.0f);

        __m128 Boy_01 = _mm_shuffle_ps(Boy, Boy, 0x50);
	dyed1_01 = _mm_mul_ps(Boy_01, ddy_01);

	// 23

        __m128 Boy_23 = _mm_shuffle_ps(Boy, Boy, 0xFA);
	dyed1_23 = _mm_mul_ps(Boy_23, ddy_23);
}

inline void VTIDenQ_Cmp_DXED_DYED_DZED(
	__m128& ddx_01,
	__m128& ddx_23,
	__m128& ddy_01,
	__m128& ddy_23,
	__m128& ddz_01,
	__m128& ddz_23,
	__m128& Boy,
        __m128& dxed1_01,
        __m128& dxed1_23,
        __m128& dyed1_01,
        __m128& dyed1_23,
	__m128& dzed2_01,
	__m128& dzed2_23
	)
{
	/*
	temp = cAzm * ddx + sAzm * ddy
	temp2 = buoy * ( cAzm * ddy - sAzm * ddx )
	temp3 = buoy * ( cDip * temp - sDip * ddz )
	temp4 = buoy * ( sdip * temp + cDip * ddz )

	dxed1 = -sAzm * temp2 + cDip * cAzm * temp3
	dyed1 = cAzm * temp2 + cDip * sAzm * temp3
	dzed1 = -sDip * temp3

	dxed2 = sDip * cAzm * temp4
	dyed2 = sDip * sAzm * temp4
	dzed2 = cDip * temp4

	VTI:

	cDip = cAzm = 1
	sDip = sAzm = 0
	
	temp  = ddx
	temp2 = buoy * ddy
	temp3 = buoy * ddx
	temp4 = buoy * ddz

	dxed1 = buoy * ddx
	dyed1 = buoy * ddy
	dzed2 = buoy * ddz

	=> VTI requires 3 stencils instead of 9 : d_dx2, d_dy2 and d_dz2.


	****

	temp = ddx
	temp2 = buoy * ddy
	temp3 = buoy * temp = buoy * ddx
	temp4 = buoy * ddz

	dxed1 = temp3 = buoy * ddx
	dyed1 = temp2 = buoy * ddy
	dzed1 = 0

	V4 = d(buoy * ddx)dx + d(buoy * ddy)dy
	V5 = d(buoy * ddz)dz

	dxed2 = 0
	dyed2 = 0
	dzed2 = temp4 = buoy * ddz

	*/

	//ddx_01 = ddx_23 = _mm_set1_ps(1.0f);
	//ddy_01 = ddy_23 = _mm_set1_ps(2.0f);
	//ddz_01 = ddz_23 = _mm_set1_ps(3.0f);

        __m128 Boy_01 = _mm_shuffle_ps(Boy, Boy, 0x50);
	dxed1_01 = _mm_mul_ps(Boy_01, ddx_01);
	dyed1_01 = _mm_mul_ps(Boy_01, ddy_01);
	dzed2_01 = _mm_mul_ps(Boy_01, ddz_01);

	// 23

        __m128 Boy_23 = _mm_shuffle_ps(Boy, Boy, 0xFA);
	dxed1_23 = _mm_mul_ps(Boy_23, ddx_23);
	dyed1_23 = _mm_mul_ps(Boy_23, ddy_23);
	dzed2_23 = _mm_mul_ps(Boy_23, ddz_23);
}

inline void TTIDenQ_Cmp_DXED(
	__m128& ddx_01,
	__m128& ddx_23,
	__m128& ddy_01,
	__m128& ddy_23,
	__m128& ddz_01,
	__m128& ddz_23,
	__m128& Boy,
	__m128& cDip,
	__m128& sDip,
	__m128& cAzm,
	__m128& sAzm,
        __m128& dxed1_01,
        __m128& dxed1_23,
        __m128& dxed2_01,
        __m128& dxed2_23
	)
{
	/*
	temp = cAzm * ddx + sAzm * ddy
	temp2 = buoy * ( cAzm * ddy - sAzm * ddx )
	temp3 = buoy * ( cDip * temp - sDip * ddz )
	temp4 = buoy * ( sdip * temp + cDip * ddz )

	dxed1 = -sAzm * temp2 + cDip * cAzm * temp3
	dyed1 = cAzm * temp2 + cDip * sAzm * temp3
	dzed1 = -sDip * temp3

	dxed2 = sDip * cAzm * temp4
	dyed2 = sDip * sAzm * temp4
	dzed2 = cDip * temp4

	VTI:

	cDip = cAzm = 1
	sDip = sAzm = 0
	
	temp  = ddx
	temp2 = buoy * ddy
	temp3 = buoy * ddx
	temp4 = buoy * ddz

	dxed1 = buoy * ddx
	dyed1 = buoy * ddy
	dzed2 = buoy * ddz

	=> VTI requires 3 stencils instead of 9 : d_dx2, d_dy2 and d_dz2.


	****

	temp = ddx
	temp2 = buoy * ddy
	temp3 = buoy * temp = buoy * ddx
	temp4 = buoy * ddz

	dxed1 = temp3 = buoy * ddx
	dyed1 = temp2 = buoy * ddy
	dzed1 = 0

	V4 = d(buoy * ddx)dx + d(buoy * ddy)dy
	V5 = d(buoy * ddz)dz

	dxed2 = 0
	dyed2 = 0
	dzed2 = temp4 = buoy * ddz

	*/

	//ddx_01 = ddx_23 = _mm_set1_ps(1.0f);
	//ddy_01 = ddy_23 = _mm_set1_ps(2.0f);
	//ddz_01 = ddz_23 = _mm_set1_ps(3.0f);

        __m128 Boy_01 = _mm_shuffle_ps(Boy, Boy, 0x50);
	__m128 cAzm_01 = _mm_shuffle_ps(cAzm, cAzm, 0x50);
	__m128 sAzm_01 = _mm_shuffle_ps(sAzm, sAzm, 0x50);
	__m128 temp_01 = _mm_add_ps(_mm_mul_ps(cAzm_01, ddx_01), _mm_mul_ps(sAzm_01, ddy_01));
	__m128 temp2_01 = _mm_mul_ps(Boy_01, _mm_sub_ps(_mm_mul_ps(cAzm_01, ddy_01), _mm_mul_ps(sAzm_01, ddx_01)));
	__m128 cDip_01 = _mm_shuffle_ps(cDip, cDip, 0x50);
        __m128 sDip_01 = _mm_shuffle_ps(sDip, sDip, 0x50);
	__m128 temp3_01 = _mm_mul_ps(Boy_01, _mm_sub_ps(_mm_mul_ps(cDip_01, temp_01), _mm_mul_ps(sDip_01, ddz_01)));
	__m128 temp4_01 = _mm_mul_ps(Boy_01, _mm_add_ps(_mm_mul_ps(sDip_01, temp_01), _mm_mul_ps(cDip_01, ddz_01)));
	
	dxed1_01 = _mm_sub_ps(_mm_mul_ps(_mm_mul_ps(cDip_01, cAzm_01), temp3_01), _mm_mul_ps(sAzm_01, temp2_01));
	dxed2_01 = _mm_mul_ps(_mm_mul_ps(sDip_01, cAzm_01), temp4_01);

	// 23

        __m128 Boy_23 = _mm_shuffle_ps(Boy, Boy, 0xFA);
	__m128 cAzm_23 = _mm_shuffle_ps(cAzm, cAzm, 0xFA);
	__m128 sAzm_23 = _mm_shuffle_ps(sAzm, sAzm, 0xFA);
	__m128 temp_23 = _mm_add_ps(_mm_mul_ps(cAzm_23, ddx_23), _mm_mul_ps(sAzm_23, ddy_23));
	__m128 temp2_23 = _mm_mul_ps(Boy_23, _mm_sub_ps(_mm_mul_ps(cAzm_23, ddy_23), _mm_mul_ps(sAzm_23, ddx_23)));
	__m128 cDip_23 = _mm_shuffle_ps(cDip, cDip, 0xFA);
        __m128 sDip_23 = _mm_shuffle_ps(sDip, sDip, 0xFA);
	__m128 temp3_23 = _mm_mul_ps(Boy_23, _mm_sub_ps(_mm_mul_ps(cDip_23, temp_23), _mm_mul_ps(sDip_23, ddz_23)));
	__m128 temp4_23 = _mm_mul_ps(Boy_23, _mm_add_ps(_mm_mul_ps(sDip_23, temp_23), _mm_mul_ps(cDip_23, ddz_23)));
	
	dxed1_23 = _mm_sub_ps(_mm_mul_ps(_mm_mul_ps(cDip_23, cAzm_23), temp3_23), _mm_mul_ps(sAzm_23, temp2_23));
	dxed2_23 = _mm_mul_ps(_mm_mul_ps(sDip_23, cAzm_23), temp4_23);
}

inline void TTIDenQ_Cmp_DYED(
	__m128& ddx_01,
	__m128& ddx_23,
	__m128& ddy_01,
	__m128& ddy_23,
	__m128& ddz_01,
	__m128& ddz_23,
	__m128& Boy,
	__m128& cDip,
	__m128& sDip,
	__m128& cAzm,
	__m128& sAzm,
        __m128& dyed1_01,
        __m128& dyed1_23,
        __m128& dyed2_01,
        __m128& dyed2_23
	)
{
	/*
	temp = cAzm * ddx + sAzm * ddy
	temp2 = buoy * ( cAzm * ddy - sAzm * ddx )
	temp3 = buoy * ( cDip * temp - sDip * ddz )
	temp4 = buoy * ( sdip * temp + cDip * ddz )

	dxed1 = -sAzm * temp2 + cDip * cAzm * temp3
	dyed1 = cAzm * temp2 + cDip * sAzm * temp3
	dzed1 = -sDip * temp3

	dxed2 = sDip * cAzm * temp4
	dyed2 = sDip * sAzm * temp4
	dzed2 = cDip * temp4

	VTI:

	cDip = cAzm = 1
	sDip = sAzm = 0
	
	temp  = ddx
	temp2 = buoy * ddy
	temp3 = buoy * ddx
	temp4 = buoy * ddz

	dxed1 = buoy * ddx
	dyed1 = buoy * ddy
	dzed2 = buoy * ddz

	=> VTI requires 3 stencils instead of 9 : d_dx2, d_dy2 and d_dz2.


	****

	temp = ddx
	temp2 = buoy * ddy
	temp3 = buoy * temp = buoy * ddx
	temp4 = buoy * ddz

	dxed1 = temp3 = buoy * ddx
	dyed1 = temp2 = buoy * ddy
	dzed1 = 0

	V4 = d(buoy * ddx)dx + d(buoy * ddy)dy
	V5 = d(buoy * ddz)dz

	dxed2 = 0
	dyed2 = 0
	dzed2 = temp4 = buoy * ddz

	*/

	//ddx_01 = ddx_23 = _mm_set1_ps(1.0f);
	//ddy_01 = ddy_23 = _mm_set1_ps(2.0f);
	//ddz_01 = ddz_23 = _mm_set1_ps(3.0f);

        __m128 Boy_01 = _mm_shuffle_ps(Boy, Boy, 0x50);
	__m128 cAzm_01 = _mm_shuffle_ps(cAzm, cAzm, 0x50);
	__m128 sAzm_01 = _mm_shuffle_ps(sAzm, sAzm, 0x50);
	__m128 temp_01 = _mm_add_ps(_mm_mul_ps(cAzm_01, ddx_01), _mm_mul_ps(sAzm_01, ddy_01));
	__m128 temp2_01 = _mm_mul_ps(Boy_01, _mm_sub_ps(_mm_mul_ps(cAzm_01, ddy_01), _mm_mul_ps(sAzm_01, ddx_01)));
	__m128 cDip_01 = _mm_shuffle_ps(cDip, cDip, 0x50);
        __m128 sDip_01 = _mm_shuffle_ps(sDip, sDip, 0x50);
	__m128 temp3_01 = _mm_mul_ps(Boy_01, _mm_sub_ps(_mm_mul_ps(cDip_01, temp_01), _mm_mul_ps(sDip_01, ddz_01)));
	__m128 temp4_01 = _mm_mul_ps(Boy_01, _mm_add_ps(_mm_mul_ps(sDip_01, temp_01), _mm_mul_ps(cDip_01, ddz_01)));
	
	dyed1_01 = _mm_add_ps(_mm_mul_ps(_mm_mul_ps(cDip_01, sAzm_01), temp3_01), _mm_mul_ps(cAzm_01, temp2_01));
	dyed2_01 = _mm_mul_ps(_mm_mul_ps(sDip_01, sAzm_01), temp4_01);

	// 23

        __m128 Boy_23 = _mm_shuffle_ps(Boy, Boy, 0xFA);
	__m128 cAzm_23 = _mm_shuffle_ps(cAzm, cAzm, 0xFA);
	__m128 sAzm_23 = _mm_shuffle_ps(sAzm, sAzm, 0xFA);
	__m128 temp_23 = _mm_add_ps(_mm_mul_ps(cAzm_23, ddx_23), _mm_mul_ps(sAzm_23, ddy_23));
	__m128 temp2_23 = _mm_mul_ps(Boy_23, _mm_sub_ps(_mm_mul_ps(cAzm_23, ddy_23), _mm_mul_ps(sAzm_23, ddx_23)));
	__m128 cDip_23 = _mm_shuffle_ps(cDip, cDip, 0xFA);
        __m128 sDip_23 = _mm_shuffle_ps(sDip, sDip, 0xFA);
	__m128 temp3_23 = _mm_mul_ps(Boy_23, _mm_sub_ps(_mm_mul_ps(cDip_23, temp_23), _mm_mul_ps(sDip_23, ddz_23)));
	__m128 temp4_23 = _mm_mul_ps(Boy_23, _mm_add_ps(_mm_mul_ps(sDip_23, temp_23), _mm_mul_ps(cDip_23, ddz_23)));
	
	dyed1_23 = _mm_add_ps(_mm_mul_ps(_mm_mul_ps(cDip_23, sAzm_23), temp3_23), _mm_mul_ps(cAzm_23, temp2_23));
	dyed2_23 = _mm_mul_ps(_mm_mul_ps(sDip_23, sAzm_23), temp4_23);
}

inline void TTIDenQ_Cmp_DXED_DYED_DZED(
	__m128& ddx_01,
	__m128& ddx_23,
	__m128& ddy_01,
	__m128& ddy_23,
	__m128& ddz_01,
	__m128& ddz_23,
	__m128& Boy,
	__m128& cDip,
	__m128& sDip,
	__m128& cAzm,
	__m128& sAzm,
        __m128& dxed1_01,
        __m128& dxed1_23,
        __m128& dxed2_01,
        __m128& dxed2_23,
        __m128& dyed1_01,
        __m128& dyed1_23,
        __m128& dyed2_01,
        __m128& dyed2_23,
	__m128& dzed1_01,
	__m128& dzed1_23,
	__m128& dzed2_01,
	__m128& dzed2_23
	)
{
	/*
	temp = cAzm * ddx + sAzm * ddy
	temp2 = buoy * ( cAzm * ddy - sAzm * ddx )
	temp3 = buoy * ( cDip * temp - sDip * ddz )
	temp4 = buoy * ( sdip * temp + cDip * ddz )

	dxed1 = -sAzm * temp2 + cDip * cAzm * temp3
	dyed1 = cAzm * temp2 + cDip * sAzm * temp3
	dzed1 = sDip * temp3

	dxed2 = sDip * cAzm * temp4
	dyed2 = sDip * sAzm * temp4
	dzed2 = cDip * temp4

	VTI:

	cDip = cAzm = 1
	sDip = sAzm = 0
	
	temp  = ddx
	temp2 = buoy * ddy
	temp3 = buoy * ddx
	temp4 = buoy * ddz

	dxed1 = buoy * ddx
	dyed1 = buoy * ddy
	dzed2 = buoy * ddz

	=> VTI requires 3 stencils instead of 9 : d_dx2, d_dy2 and d_dz2.


	****

	temp = ddx
	temp2 = buoy * ddy
	temp3 = buoy * temp = buoy * ddx
	temp4 = buoy * ddz

	dxed1 = temp3 = buoy * ddx
	dyed1 = temp2 = buoy * ddy
	dzed1 = 0

	V4 = d(buoy * ddx)dx + d(buoy * ddy)dy
	V5 = d(buoy * ddz)dz

	dxed2 = 0
	dyed2 = 0
	dzed2 = temp4 = buoy * ddz

	*/

	/*
	ddx_01 = ddx_23 = _mm_set1_ps(1.0f);
	ddy_01 = ddy_23 = _mm_set1_ps(1.0f);
	ddz_01 = ddz_23 = _mm_set1_ps(1.0f);
	*/

        __m128 Boy_01 = _mm_shuffle_ps(Boy, Boy, 0x50);
	__m128 cAzm_01 = _mm_shuffle_ps(cAzm, cAzm, 0x50);
	__m128 sAzm_01 = _mm_shuffle_ps(sAzm, sAzm, 0x50);
	__m128 temp_01 = _mm_add_ps(_mm_mul_ps(cAzm_01, ddx_01), _mm_mul_ps(sAzm_01, ddy_01));
	__m128 temp2_01 = _mm_mul_ps(Boy_01, _mm_sub_ps(_mm_mul_ps(cAzm_01, ddy_01), _mm_mul_ps(sAzm_01, ddx_01)));
	__m128 cDip_01 = _mm_shuffle_ps(cDip, cDip, 0x50);
        __m128 sDip_01 = _mm_shuffle_ps(sDip, sDip, 0x50);
	__m128 temp3_01 = _mm_mul_ps(Boy_01, _mm_sub_ps(_mm_mul_ps(cDip_01, temp_01), _mm_mul_ps(sDip_01, ddz_01)));
	__m128 temp4_01 = _mm_mul_ps(Boy_01, _mm_add_ps(_mm_mul_ps(sDip_01, temp_01), _mm_mul_ps(cDip_01, ddz_01)));
	
	dxed1_01 = _mm_sub_ps(_mm_mul_ps(_mm_mul_ps(cDip_01, cAzm_01), temp3_01), _mm_mul_ps(sAzm_01, temp2_01));
	dyed1_01 = _mm_add_ps(_mm_mul_ps(_mm_mul_ps(cDip_01, sAzm_01), temp3_01), _mm_mul_ps(cAzm_01, temp2_01));
	dzed1_01 = _mm_mul_ps(sDip_01, temp3_01);

	dxed2_01 = _mm_mul_ps(_mm_mul_ps(sDip_01, cAzm_01), temp4_01);
	dyed2_01 = _mm_mul_ps(_mm_mul_ps(sDip_01, sAzm_01), temp4_01);
	dzed2_01 = _mm_mul_ps(cDip_01, temp4_01);

	// 23

        __m128 Boy_23 = _mm_shuffle_ps(Boy, Boy, 0xFA);
	__m128 cAzm_23 = _mm_shuffle_ps(cAzm, cAzm, 0xFA);
	__m128 sAzm_23 = _mm_shuffle_ps(sAzm, sAzm, 0xFA);
	__m128 temp_23 = _mm_add_ps(_mm_mul_ps(cAzm_23, ddx_23), _mm_mul_ps(sAzm_23, ddy_23));
	__m128 temp2_23 = _mm_mul_ps(Boy_23, _mm_sub_ps(_mm_mul_ps(cAzm_23, ddy_23), _mm_mul_ps(sAzm_23, ddx_23)));
	__m128 cDip_23 = _mm_shuffle_ps(cDip, cDip, 0xFA);
        __m128 sDip_23 = _mm_shuffle_ps(sDip, sDip, 0xFA);
	__m128 temp3_23 = _mm_mul_ps(Boy_23, _mm_sub_ps(_mm_mul_ps(cDip_23, temp_23), _mm_mul_ps(sDip_23, ddz_23)));
	__m128 temp4_23 = _mm_mul_ps(Boy_23, _mm_add_ps(_mm_mul_ps(sDip_23, temp_23), _mm_mul_ps(cDip_23, ddz_23)));
	
	dxed1_23 = _mm_sub_ps(_mm_mul_ps(_mm_mul_ps(cDip_23, cAzm_23), temp3_23), _mm_mul_ps(sAzm_23, temp2_23));
	dyed1_23 = _mm_add_ps(_mm_mul_ps(_mm_mul_ps(cDip_23, sAzm_23), temp3_23), _mm_mul_ps(cAzm_23, temp2_23));
	dzed1_23 = _mm_mul_ps(sDip_23, temp3_23);

	dxed2_23 = _mm_mul_ps(_mm_mul_ps(sDip_23, cAzm_23), temp4_23);
	dyed2_23 = _mm_mul_ps(_mm_mul_ps(sDip_23, sAzm_23), temp4_23);
	dzed2_23 = _mm_mul_ps(cDip_23, temp4_23);

	/*
	float* dxed1 = (float*)&dxed1_01;
	float* dyed1 = (float*)&dyed1_01;
	float* dzed1 = (float*)&dzed1_01;

	float* dxed2 = (float*)&dxed2_01;
	float* dyed2 = (float*)&dyed2_01;
	float* dzed2 = (float*)&dzed2_01;

	if (*dxed1 > 1e-3f && *dyed1 > 1e-3f)
	{

		printf("dxed1 = %e\n",*dxed1);
		printf("dyed1 = %e\n",*dyed1);
		printf("dzed1 = %e\n",*dzed1);

		printf("dxed2 = %e\n",*dxed2);
		printf("dyed2 = %e\n",*dyed2);
		printf("dzed2 = %e\n",*dzed2);

		printf("\n");

		exit(0);
	}
	*/
}

inline void ISODenQ_Cmp_Wave_Equation(
	int stage,
	__m128* pq1,
	__m128* rs1,
	__m128* Apq1,
	int* VelDen1,
	__m128& spgzyx,
	__m128& V4_0123,
	__m128& V5_0123
	)
{
	__m128 Qatten, C33;
	ISODenQ_Get_EM2(VelDen1, Qatten, C33);
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

	__m128 rs_0123 = _mm_mul_ps(C33, _mm_add_ps(V5_0123, V4_0123));

	if (stage == 0)
	{
		// 2nd order in time
		rs1[0] = _mm_mul_ps(_mm_mul_ps(_mm_add_ps(_mm_sub_ps(_mm_add_ps(pq1[0], pq1[0]), _mm_mul_ps(spgzyx, rs1[0])), rs_0123), spgzyx), Qatten);

		//rs_0123 = _mm_add_ps(rs_0123, _mm_add_ps(pq1[0], pq1[0]));
		//rs_0123 = _mm_sub_ps(rs_0123, _mm_mul_ps(spgzyx, rs1[0]));
		//rs1[0] = _mm_mul_ps(rs_0123, spgzyx);
	}
	else if (stage == 1)
	{
		// 4th order in time - 1st call
		_mm_stream_ps((float*)(Apq1), rs_0123);
		//Apq1[0] = rs_01;
		//Apq1[1] = rs_23;
	}
	else
	{
		// 4th order in time - 2nd and final call
		rs_0123 = _mm_mul_ps(rs_0123, _mm_set1_ps(1.0f/12.0f));
		rs_0123 = _mm_add_ps(rs_0123, Apq1[0]);
                rs_0123 = _mm_add_ps(rs_0123, _mm_add_ps(pq1[0], pq1[0]));
                rs_0123 = _mm_sub_ps(rs_0123, _mm_mul_ps(spgzyx, rs1[0]));
		rs1[0] = _mm_mul_ps(_mm_mul_ps(rs_0123, spgzyx), Qatten);
	}
}

inline void TTIDenQ_Cmp_Wave_Equation(
	int stage,
	__m128* pq1,
	__m128* rs1,
	__m128* Apq1,
	int* VelAnis1,
	int* DenAng1,
	__m128& spgzyx,
	__m128& V4_01,
	__m128& V4_23,
	__m128& V5_01,
	__m128& V5_23
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
	/*
	rs1[0] = rs1[1] = _mm_set1_ps(0.0f);
	pq1[0] = pq1[1] = _mm_set1_ps(0.0f);
	V4_01 = V4_23 = _mm_set1_ps(1.0f);
	V5_01 = V5_23 = _mm_set1_ps(1.0f);
	*/

	__m128 rs_01 = _mm_add_ps(_mm_mul_ps(V5_01, C66C44_01), _mm_mul_ps(V4_01, C44C33_01));
	__m128 V4qV5p_01 = _mm_shuffle_ps(V4_01, V5_01, 0x8D);
	V4qV5p_01 = _mm_shuffle_ps(V4qV5p_01, V4qV5p_01, 0xD8);
	rs_01 = _mm_add_ps(rs_01, _mm_mul_ps(V4qV5p_01, C55_01));

	__m128 rs_23 = _mm_add_ps(_mm_mul_ps(V5_23, C66C44_23), _mm_mul_ps(V4_23, C44C33_23));
	__m128 V4qV5p_23 = _mm_shuffle_ps(V4_23, V5_23, 0x8D);
	V4qV5p_23 = _mm_shuffle_ps(V4qV5p_23, V4qV5p_23, 0xD8);
	rs_23 = _mm_add_ps(rs_23, _mm_mul_ps(V4qV5p_23, C55_23));

	/*
	static int been_there_done_that = 0;
	if (!been_there_done_that)
	{
		been_there_done_that = 1;
		float _r, _s;
		_r = ((float*)&rs_01)[0];
		_s = ((float*)&rs_01)[1];
		printf("r = %f, s = %f\n",_r,_s);
	}
	*/

	if (stage == 0)
	{
		// 2nd order in time
		__m128 spgzyx_01 = _mm_shuffle_ps(spgzyx, spgzyx, 0x50);
		rs_01 = _mm_add_ps(rs_01, _mm_add_ps(pq1[0], pq1[0]));
		rs_01 = _mm_sub_ps(rs_01, _mm_mul_ps(spgzyx_01, rs1[0]));
		__m128 Qatten_01 = _mm_shuffle_ps(Qatten, Qatten, 0x50);
		rs1[0] = _mm_mul_ps(_mm_mul_ps(rs_01, spgzyx_01), Qatten_01);

		__m128 spgzyx_23 = _mm_shuffle_ps(spgzyx, spgzyx, 0xFA);
		rs_23 = _mm_add_ps(rs_23, _mm_add_ps(pq1[1], pq1[1]));
		rs_23 = _mm_sub_ps(rs_23, _mm_mul_ps(spgzyx_23, rs1[1]));
		__m128 Qatten_23 = _mm_shuffle_ps(Qatten, Qatten, 0xFA);
		rs1[1] = _mm_mul_ps(_mm_mul_ps(rs_23, spgzyx_23), Qatten_23);
	}
	else if (stage == 1)
	{
		// 4th order in time - 1st call
		_mm_stream_ps((float*)(Apq1), rs_01);
		_mm_stream_ps((float*)(Apq1+1), rs_23);
		//Apq1[0] = rs_01;
		//Apq1[1] = rs_23;
	}
	else
	{
		// 4th order in time - 2nd and final call
		__m128 spgzyx_01 = _mm_shuffle_ps(spgzyx, spgzyx, 0x50);
		rs_01 = _mm_mul_ps(rs_01, _mm_set1_ps(1.0f/12.0f));
		rs_01 = _mm_add_ps(rs_01, Apq1[0]);
                rs_01 = _mm_add_ps(rs_01, _mm_add_ps(pq1[0], pq1[0]));
                rs_01 = _mm_sub_ps(rs_01, _mm_mul_ps(spgzyx_01, rs1[0]));
		__m128 Qatten_01 = _mm_shuffle_ps(Qatten, Qatten, 0x50);
		rs1[0] = _mm_mul_ps(_mm_mul_ps(rs_01, spgzyx_01), Qatten_01);

		__m128 spgzyx_23 = _mm_shuffle_ps(spgzyx, spgzyx, 0xFA);
		rs_23 = _mm_mul_ps(rs_23, _mm_set1_ps(1.0f/12.0f));
		rs_23 = _mm_add_ps(rs_23, Apq1[1]);
                rs_23 = _mm_add_ps(rs_23, _mm_add_ps(pq1[1], pq1[1]));
                rs_23 = _mm_sub_ps(rs_23, _mm_mul_ps(spgzyx_23, rs1[1]));
		__m128 Qatten_23 = _mm_shuffle_ps(Qatten, Qatten, 0xFA);
		rs1[1] = _mm_mul_ps(_mm_mul_ps(rs_23, spgzyx_23), Qatten_23);
	}
}

void ISODenQ_Process_Patch_Leadin(
	int logLevel,
	__m128* pq,		// pq for stage 0 and 1, Apq for stage 2
	int* VelDen,
	__m128* stripe,
	__m128* dyed_buf,
	__m128* V4,
	int iX0,
	int iY0,
	int iZ,
	int iXN_halo,
	int iXN,
	int iYN_halo,
	int iYN,
	int stride_y_m128,
	int stride_z_m128,
	int stride_y_em,
	int stride_z_em
	)
{
	//float neg_iZ_fact = (iZ < 0) ? 0.8f : 1.0f;
	int abs_iZ0 = (iZ >= 0) ? iZ : -(iZ+1);
	int abs_iZ1 = ((iZ+1) >= 0) ? (iZ+1) : -(iZ+2);

	__m128* pq0 = pq + (unsigned long)abs_iZ0 * (unsigned long)stride_z_m128 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m128 + (unsigned long)((iX0-4) >> 1);
	int* VelDen0 = VelDen + (unsigned long)abs_iZ0 * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

	int abs_stride_z_m128 = (abs_iZ1 - abs_iZ0) * stride_z_m128;
	int abs_stride_z_em = (abs_iZ1 - abs_iZ0) * stride_z_em;
	if (logLevel >= 6)
	{
		printf("iZ = %d, abs_iZ0 = %d, abs_iZ1 = %d, abs_stride_z_m128 = %d, abs_stride_em = %d\n",iZ,abs_iZ0,abs_iZ1,abs_stride_z_m128,abs_stride_z_em);
		fflush(stdout);
	}

	__m128* V4_0 = V4;
	__m128* V4_1 = V4;

	__m128* dyed_buf_0 = dyed_buf;
	__m128* dyed_buf_1 = dyed_buf + 8 * iXN;

	for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
	{
		pq0 += 2;
		VelDen0 += 8;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
			_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

			__m128 ddy_0123;
			ISO_Cmp_DDY(pq0,stride_y_m128,ddy_0123);

			__m128 Boy;
			ISODenQ_Get_EM1(VelDen0, Boy);

			ISODenQ_Cmp_DYED(
					ddy_0123,  
					Boy,
					dyed_buf_0[0]
					);

			ISO_Cmp_DDY(pq0+abs_stride_z_m128,stride_y_m128,ddy_0123);

			ISODenQ_Get_EM1(VelDen0+abs_stride_z_em, Boy);

			ISODenQ_Cmp_DYED(
					ddy_0123,  
					Boy,
					dyed_buf_0[1]
					);
			dyed_buf_0 += 2;

			pq0 += 1;
			VelDen0 += 4;
		}

		pq0 += stride_y_m128 - (iXN+2);
		VelDen0 += stride_y_em - (iXN+2)*4;
	}

	for (int iY = (iYN_halo-5);  iY > 4;  --iY)
	{
		__m128* stripe0 = stripe;

		for (int iX = 2;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
			_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

			__m128 ddx_0123;
			ISO_Cmp_DDX(pq0,ddx_0123);

			__m128 Boy;
			ISODenQ_Get_EM1(VelDen0, Boy);

			ISODenQ_Cmp_DXED(
					ddx_0123,  
					Boy,
					stripe0[ 0]
					);

			ISO_Cmp_DDX(pq0+abs_stride_z_m128,ddx_0123);

			ISODenQ_Get_EM1(VelDen0+abs_stride_z_em, Boy);

			ISODenQ_Cmp_DXED(
					ddx_0123,  
					Boy,
					stripe0[ 1]
					);
			stripe0 += 2;

			pq0 += 1;
			VelDen0 += 4;
		}

		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
			_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

			__m128 ddx_0123, ddy_0123, ddz_0123;
			ISO_Cmp_DDX(pq0,ddx_0123);
			ISO_Cmp_DDY(pq0,stride_y_m128,ddy_0123);
			if (abs_iZ0 == 0)
				ISO_Cmp_DDZ_Z0(pq0,stride_z_m128,ddz_0123);
			else if (abs_iZ0 == 1)
				ISO_Cmp_DDZ_Z1(pq0,stride_z_m128,ddz_0123);
			else if (abs_iZ0 == 2)
				ISO_Cmp_DDZ_Z2(pq0,stride_z_m128,ddz_0123);
			else if (abs_iZ0 == 3)
				ISO_Cmp_DDZ_Z3(pq0,stride_z_m128,ddz_0123);
			else
				ISO_Cmp_DDZ(pq0,stride_z_m128,ddz_0123);

			__m128 Boy;
			ISODenQ_Get_EM1(VelDen0, Boy);

			ISODenQ_Cmp_DXED_DYED_DZED(
					ddx_0123, ddy_0123, ddz_0123, 
					Boy,
					stripe0[ 0],
					dyed_buf_0[0],
					V4_0[7]
					);
			//V4_0[7] = _mm_mul_ps(_mm_set1_ps(neg_iZ_fact), V4_0[7]);

			ISO_Cmp_DDX(pq0+abs_stride_z_m128,ddx_0123);
			ISO_Cmp_DDY(pq0+abs_stride_z_m128,stride_y_m128,ddy_0123);
			if (abs_iZ1 == 0)
				ISO_Cmp_DDZ_Z0(pq0+abs_stride_z_m128,stride_z_m128,ddz_0123);
			else if (abs_iZ1 == 1)
				ISO_Cmp_DDZ_Z1(pq0+abs_stride_z_m128,stride_z_m128,ddz_0123);
			else if (abs_iZ1 == 2)
				ISO_Cmp_DDZ_Z2(pq0+abs_stride_z_m128,stride_z_m128,ddz_0123);
			else if (abs_iZ1 == 3)
				ISO_Cmp_DDZ_Z3(pq0+abs_stride_z_m128,stride_z_m128,ddz_0123);
			else
				ISO_Cmp_DDZ(pq0+abs_stride_z_m128,stride_z_m128,ddz_0123);

			ISODenQ_Get_EM1(VelDen0+abs_stride_z_em, Boy);

			ISODenQ_Cmp_DXED_DYED_DZED(
					ddx_0123, ddy_0123, ddz_0123, 
					Boy,
					stripe0[ 1],
					dyed_buf_0[1],
					V4_0[6]
					);
			//V4_0[6] = _mm_mul_ps(_mm_set1_ps(neg_iZ_fact), V4_0[6]);

			stripe0 += 2;	
			dyed_buf_0 += 2;
			V4_0 += 17;

			pq0 += 1;
			VelDen0 += 4;
		}
		V4_0 -= 17 * iXN;

		for (int iX = 1;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
			_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

			__m128 ddx_0123;
			ISO_Cmp_DDX(pq0,ddx_0123);

			__m128 Boy;
			ISODenQ_Get_EM1(VelDen0, Boy);

			ISODenQ_Cmp_DXED(
					ddx_0123,   
					Boy,
					stripe0[ 0]
					);

			ISO_Cmp_DDX(pq0+abs_stride_z_m128,ddx_0123);

			ISODenQ_Get_EM1(VelDen0+abs_stride_z_em, Boy);

			ISODenQ_Cmp_DXED(
					ddx_0123,  
					Boy,
					stripe0[ 1]
					);
			stripe0 += 2;

			pq0 += 1;
			VelDen0 += 4;
		}

		pq0 += stride_y_m128 - iXN_halo;
		VelDen0 += stride_y_em - iXN_halo*4;

		stripe0 = stripe + 4;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			if (iY <= (iYN_halo-9))
			{
				// dydyed1, dydyed2
				__m128 dydyed1_0123;
				ISO_Cmp_DDY(dyed_buf_1,2*iXN,dydyed1_0123);
				V4_1[1] = _mm_add_ps(V4_1[1], dydyed1_0123);

				ISO_Cmp_DDY(dyed_buf_1+1,2*iXN,dydyed1_0123);
				V4_1[0] = _mm_add_ps(V4_1[0], dydyed1_0123);

				dyed_buf_1 += 2;
				V4_1 += 17;
			}

			// dxdxed1, dxdxed2
			__m128 dxdxed1_0123;
			ISO_Cmp_DDX_1x_stride2(stripe0,dxdxed1_0123);
			V4_0[1] = dxdxed1_0123;

			ISO_Cmp_DDX_1x_stride2(stripe0+1,dxdxed1_0123);
			V4_0[0] = dxdxed1_0123;

			V4_0 += 17;

			stripe0 += 2;	
		}
	}

	for (int iY = 4;  iY > 0;  --iY)
	{
		pq0 += 2;
		VelDen0 += 8;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
			_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

			__m128 ddy_0123;
			ISO_Cmp_DDY(pq0,stride_y_m128,ddy_0123);

			__m128 Boy;
			ISODenQ_Get_EM1(VelDen0, Boy);

			ISODenQ_Cmp_DYED(
					ddy_0123,  
					Boy,
					dyed_buf_0[0]
					);

			ISO_Cmp_DDY(pq0+abs_stride_z_m128,stride_y_m128,ddy_0123);

			ISODenQ_Get_EM1(VelDen0+abs_stride_z_em, Boy);

			ISODenQ_Cmp_DYED(
					ddy_0123,  
					Boy,
					dyed_buf_0[1]
					);
			dyed_buf_0 += 2;

			pq0 += 1;
			VelDen0 += 4;

			// dydyed1, dydyed2
			__m128 dydyed1_0123;
			ISO_Cmp_DDY(dyed_buf_1,2*iXN,dydyed1_0123);
			V4_1[1] = _mm_add_ps(V4_1[1], dydyed1_0123);

			ISO_Cmp_DDY(dyed_buf_1+1,2*iXN,dydyed1_0123);
			V4_1[0] = _mm_add_ps(V4_1[0], dydyed1_0123);

			dyed_buf_1 += 2;
			V4_1 += 17;
		}

		pq0 += stride_y_m128 - (iXN+2);
		VelDen0 += stride_y_em - (iXN+2)*4;
	}
}

//
// Valid from [1,dimz-1>
//
void ISODenQ_Process_Patch(
	int logLevel,
	int stage,
	__m128* pq,
	__m128* rs,
	__m128* Apq,			// pass pq ptr if 2nd order in time
	int* VelDen,
	__m128* _mm_spgx,
	__m128* _mm_spgy,
	__m128 _mm_spgz0,
	__m128 _mm_spgz1,
	__m128* stripe,
	__m128* dyed_buf,
	__m128* V4,
	int iX0,
	int iY0,
	int iZ,
	int dimz,
	int iXN_halo,
	int iXN,
	int iYN_halo,
	int iYN,
	int stride_y_m128,
	int stride_z_m128,
	int stride_y_em,
	int stride_z_em
#ifdef TMJ_TIMING
	,unsigned long& thrcnt2
#endif
	)
{
	int dz_edge_diff = dimz - iZ - 1;  // dimz-1 -> 0, dimz-2 -> 1 etc.

	__m128* pq1 = pq + (unsigned long)iZ * (unsigned long)stride_z_m128 + (unsigned long)iY0 * (unsigned long)stride_y_m128 + (unsigned long)(iX0 >> 1);
	__m128* rs1 = rs + (unsigned long)iZ * (unsigned long)stride_z_m128 + (unsigned long)iY0 * (unsigned long)stride_y_m128 + (unsigned long)(iX0 >> 1);
	__m128* Apq1 = Apq + (unsigned long)iZ * (unsigned long)stride_z_m128 + (unsigned long)iY0 * (unsigned long)stride_y_m128 + (unsigned long)(iX0 >> 1);
	int* VelDen1 = VelDen + (unsigned long)iZ * (unsigned long)stride_z_em + (unsigned long)iY0 * (unsigned long)stride_y_em + (unsigned long)(iX0*2);

	if (logLevel >= 6)
	{
		printf("Process_Patch :: iZ = %d\n",iZ);
		fflush(stdout);
	}

	// process next tile ("patch")
	if (dz_edge_diff > 5)
	{
	 	__m128* pq0 = ((stage == 0 || stage == 1) ? pq : Apq) + (unsigned long)(iZ+4) * (unsigned long)stride_z_m128 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m128 + (unsigned long)((iX0-4) >> 1);
		int* VelDen0 = VelDen + (unsigned long)(iZ+4) * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

		__m128* dyed_buf_0 = dyed_buf;
		__m128* dyed_buf_1 = dyed_buf + 8 * iXN;

		__m128* V4_0 = V4;
		__m128* V4_1 = V4;

		for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
		{
			pq0 += 2;
			VelDen0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddy_0123;
				ISO_Cmp_DDY(pq0,stride_y_m128,ddy_0123);

				__m128 Boy;
				ISODenQ_Get_EM1(VelDen0, Boy);

				ISODenQ_Cmp_DYED(
						ddy_0123,  
						Boy,
						dyed_buf_0[0]
						);

				ISO_Cmp_DDY(pq0+stride_z_m128,stride_y_m128,ddy_0123);

				ISODenQ_Get_EM1(VelDen0+stride_z_em, Boy);

				ISODenQ_Cmp_DYED(
						ddy_0123,  
						Boy,
                                                dyed_buf_0[1]
						);

				dyed_buf_0 += 2;
				pq0 += 1;
				VelDen0 += 4;
			}
			pq0 += stride_y_m128 - (iXN+2);
			VelDen0 += stride_y_em - (iXN+2)*4;
		}

		for (int iY = (iYN_halo-5);  iY > 4;  --iY)
		{
			__m128* stripe0 = stripe;

			for (int iX = 2;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddx_0123;
				ISO_Cmp_DDX(pq0,ddx_0123);

				__m128 Boy;
				ISODenQ_Get_EM1(VelDen0, Boy);

				ISODenQ_Cmp_DXED(
						ddx_0123,  
						Boy,
						stripe0[ 0]
						);

				ISO_Cmp_DDX(pq0+stride_z_m128,ddx_0123);

				ISODenQ_Get_EM1(VelDen0+stride_z_em, Boy);

				ISODenQ_Cmp_DXED(
						ddx_0123,  
						Boy,
						stripe0[ 1]
						);
				stripe0 += 2;

				pq0 += 1;
				VelDen0 += 4;
			}

			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddx_0123, ddy_0123, ddz_0123;
				ISO_Cmp_DDX(pq0,ddx_0123);
				ISO_Cmp_DDY(pq0,stride_y_m128,ddy_0123);
				if (dz_edge_diff > 9)
				{
					ISO_Cmp_DDZ(pq0,stride_z_m128,ddz_0123);
				}
				else
				{
					ISO_Cmp_DDZ_EE(pq0,stride_z_m128,ddz_0123);
				}

				__m128 Boy;
				ISODenQ_Get_EM1(VelDen0, Boy);

				ISODenQ_Cmp_DXED_DYED_DZED(
						ddx_0123, ddy_0123, ddz_0123, 
						Boy,
						stripe0[ 0],
						dyed_buf_0[0],
						V4_0[7]
						);

				ISO_Cmp_DDX(pq0+stride_z_m128,ddx_0123);
				ISO_Cmp_DDY(pq0+stride_z_m128,stride_y_m128,ddy_0123);
				if (dz_edge_diff > 10)
				{
					ISO_Cmp_DDZ(pq0+stride_z_m128,stride_z_m128,ddz_0123);
				}
				else
				{
					ISO_Cmp_DDZ_EE(pq0+stride_z_m128,stride_z_m128,ddz_0123);
				}

				ISODenQ_Get_EM1(VelDen0+stride_z_em, Boy);

				ISODenQ_Cmp_DXED_DYED_DZED(
						ddx_0123, ddy_0123, ddz_0123, 
						Boy,
						stripe0[ 1],
                                                dyed_buf_0[1],
						V4_0[6]
						);
				stripe0 += 2;
				V4_0 += 17;

				dyed_buf_0 += 2;
				pq0 += 1;
				VelDen0 += 4;
			}
			V4_0 -= 17 * iXN;

			for (int iX = 1;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddx_0123;
				ISO_Cmp_DDX(pq0,ddx_0123);

				__m128 Boy;
				ISODenQ_Get_EM1(VelDen0, Boy);

				ISODenQ_Cmp_DXED(
						ddx_0123,  
						Boy,
						stripe0[ 0]
						);

				ISO_Cmp_DDX(pq0+stride_z_m128,ddx_0123);

				ISODenQ_Get_EM1(VelDen0+stride_z_em, Boy);

				ISODenQ_Cmp_DXED(
						ddx_0123,  
						Boy,
						stripe0[ 1]
						);
				stripe0 += 2;

				pq0 += 1;
				VelDen0 += 4;
			}

			pq0 += stride_y_m128 - iXN_halo;
			VelDen0 += stride_y_em - iXN_halo*4;

			stripe0 = stripe + 4;
			__m128 _mm_spgzy0 = _mm_mul_ps(_mm_spgz0, _mm_spgy[iY]);
			__m128 _mm_spgzy1 = _mm_mul_ps(_mm_spgz1, _mm_spgy[iY]);
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(VelDen1), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelDen1+stride_z_em), _MM_HINT_NTA);
				if (stage == 0)
				{
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m128), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}

				if (iY <= (iYN_halo-9))
				{	
					// dydyed1, dydyed2
					__m128 dydyed1_0123;
					ISO_Cmp_DDY(dyed_buf_1,2*iXN,dydyed1_0123);
					V4_1[1] = _mm_add_ps(V4_1[1], dydyed1_0123);

					ISO_Cmp_DDY(dyed_buf_1+1,2*iXN,dydyed1_0123);
					V4_1[0] = _mm_add_ps(V4_1[0], dydyed1_0123);

					dyed_buf_1 += 2;
					V4_1 += 17;
				}

				// dxdxed1, dxdxed2
				__m128 dxdxed1_0123;
				ISO_Cmp_DDX_1x_stride2(stripe0,dxdxed1_0123);
				V4_0[1] = dxdxed1_0123;

				ISO_Cmp_DDX_1x_stride2(stripe0+1,dxdxed1_0123);
				V4_0[0] = dxdxed1_0123;

				// dzdzed1, dzdzed2
				__m128 V5_0123 = V4_0[5];

				__m128 dzdzed2_0123;
				ISO_Cmp_DDZ(V4_0+12,-1,dzdzed2_0123);
				__m128 V4_0123 = dzdzed2_0123;

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				ISODenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelDen1,spgzyx0,V4_0123,V5_0123
						);

				V5_0123 = V4_0[ 4];

				ISO_Cmp_DDZ(V4_0+11,-1,dzdzed2_0123);
				V4_0123 = dzdzed2_0123;

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				ISODenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m128,rs1+stride_z_m128,Apq1+stride_z_m128,VelDen1+stride_z_em,spgzyx1,V4_0123,V5_0123
						);

				V4_0 += 17;

				pq1 += 1;
				rs1 += 1;
				Apq1 += 1;
				VelDen1 += 4;
				stripe0 += 2;	
#ifdef TMJ_TIMING
				thrcnt2+=2;
#endif
			}
			pq1 += stride_y_m128 - iXN;
			rs1 += stride_y_m128 - iXN;
			Apq1 += stride_y_m128 - iXN;
			VelDen1 += stride_y_em - iXN*4;
		}

		for (int iY = 4;  iY > 0;  --iY)
		{
			pq0 += 2;
			VelDen0 += 8;	
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddy_0123;
				ISO_Cmp_DDY(pq0,stride_y_m128,ddy_0123);

				__m128 Boy;
				ISODenQ_Get_EM1(VelDen0, Boy);

				ISODenQ_Cmp_DYED(
						ddy_0123,  
						Boy,
						dyed_buf_0[0]
						);

				ISO_Cmp_DDY(pq0+stride_z_m128,stride_y_m128,ddy_0123);

				ISODenQ_Get_EM1(VelDen0+stride_z_em, Boy);

				ISODenQ_Cmp_DYED(
						ddy_0123,  
						Boy,
                                                dyed_buf_0[1]
						);

				dyed_buf_0 += 2;
				pq0 += 1;
				VelDen0 += 4;

				// dydyed1, dydyed2
				__m128 dydyed1_0123;
				ISO_Cmp_DDY(dyed_buf_1,2*iXN,dydyed1_0123);
				V4_1[1] = _mm_add_ps(V4_1[1], dydyed1_0123);

				ISO_Cmp_DDY(dyed_buf_1+1,2*iXN,dydyed1_0123);
				V4_1[0] = _mm_add_ps(V4_1[0], dydyed1_0123);

				dyed_buf_1 += 2;
				V4_1 += 17;
			}
	
			pq0 += stride_y_m128 - (iXN+2);
			VelDen0 += stride_y_em - (iXN+2)*4;
		}
	}
	else if (dz_edge_diff == 5)
	{
	 	__m128* pq0 = ((stage == 0 || stage == 1) ? pq : Apq) + (unsigned long)(iZ+4) * (unsigned long)stride_z_m128 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m128 + (unsigned long)((iX0-4) >> 1);
		int* VelDen0 = VelDen + (unsigned long)(iZ+4) * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

		__m128* dyed_buf_0 = dyed_buf;
		__m128* dyed_buf_1 = dyed_buf + 8 * iXN;

		__m128* V4_0 = V4;
		__m128* V4_1 = V4;

		for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
		{
			pq0 += 2;
			VelDen0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

				__m128 ddy_0123;
				ISO_Cmp_DDY(pq0,stride_y_m128,ddy_0123);

				__m128 Boy;
				ISODenQ_Get_EM1(VelDen0, Boy);

				ISODenQ_Cmp_DYED(
						ddy_0123,  
						Boy,
						dyed_buf_0[0]
						);
				dyed_buf_0 += 2;

				pq0 += 1;
				VelDen0 += 4;
			}

			pq0 += stride_y_m128 - (iXN+2);
			VelDen0 += stride_y_em - (iXN+2)*4;
		}

		for (int iY = (iYN_halo-5);  iY > 4;  --iY)
		{
			__m128* stripe0 = stripe;

			for (int iX = 2;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

				__m128 ddx_0123;
				ISO_Cmp_DDX(pq0,ddx_0123);

				__m128 Boy;
				ISODenQ_Get_EM1(VelDen0, Boy);

				ISODenQ_Cmp_DXED(
						ddx_0123,  
						Boy,
						stripe0[ 0]
						);
				stripe0 += 2;

				pq0 += 1;
				VelDen0 += 4;
			}

			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

				__m128 ddx_0123, ddy_0123, ddz_0123;
				ISO_Cmp_DDX(pq0,ddx_0123);
				ISO_Cmp_DDY(pq0,stride_y_m128,ddy_0123);
				ISO_Cmp_DDZ_EE(pq0,stride_z_m128,ddz_0123);

				__m128 Boy;
				ISODenQ_Get_EM1(VelDen0, Boy);

				ISODenQ_Cmp_DXED_DYED_DZED(
						ddx_0123, ddy_0123, ddz_0123, 
						Boy,
						stripe0[ 0],
						dyed_buf_0[0],
						V4_0[7]
						);
				stripe0 += 2;
				dyed_buf_0 += 2;
				V4_0 += 17;

				pq0 += 1;
				VelDen0 += 4;
			}
			V4_0 -= 17 * iXN;

			for (int iX = 1;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

				__m128 ddx_0123;
				ISO_Cmp_DDX(pq0,ddx_0123);

				__m128 Boy;
				ISODenQ_Get_EM1(VelDen0, Boy);

				ISODenQ_Cmp_DXED(
						ddx_0123,  
						Boy,
						stripe0[ 0]
						);
				stripe0 += 2;

				pq0 += 1;
				VelDen0 += 4;
			}

			pq0 += stride_y_m128 - iXN_halo;
			VelDen0 += stride_y_em - iXN_halo*4;

			stripe0 = stripe + 4;
			__m128 _mm_spgzy0 = _mm_mul_ps(_mm_spgz0, _mm_spgy[iY]);
			__m128 _mm_spgzy1 = _mm_mul_ps(_mm_spgz1, _mm_spgy[iY]);
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(VelDen1), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelDen1+stride_z_em), _MM_HINT_NTA);
				if (stage == 0)
				{
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m128), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}

				if (iY <= (iYN_halo-9))
				{
					// dydyed1, dydyed2
					__m128 dydyed1_0123;
					ISO_Cmp_DDY(dyed_buf_1,2*iXN,dydyed1_0123);
					V4_1[1] = _mm_add_ps(V4_1[1], dydyed1_0123);

					dyed_buf_1 += 2;
					V4_1 += 17;
				}

				// dxdxed1, dxdxed2
				__m128 dxdxed1_0123;
				ISO_Cmp_DDX_1x_stride2(stripe0,dxdxed1_0123);
				V4_0[1] = dxdxed1_0123;

				// dzdzed1, dzdzed2
				__m128 V5_0123 = V4_0[5];

				__m128 dzdzed2_0123;
				ISO_Cmp_DDZ(V4_0+12,-1,dzdzed2_0123);
				__m128 V4_0123 = dzdzed2_0123;

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				ISODenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelDen1,spgzyx0,V4_0123,V5_0123
						);

				V5_0123 = V4_0[ 4];

				ISO_Cmp_DDZ_EE(V4_0+11,-1,dzdzed2_0123);
				V4_0123 = dzdzed2_0123;

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				ISODenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m128,rs1+stride_z_m128,Apq1+stride_z_m128,VelDen1+stride_z_em,spgzyx1,V4_0123,V5_0123
						);

				V4_0 += 17;

				pq1 += 1;
				rs1 += 1;
				Apq1 += 1;
				VelDen1 += 4;
				stripe0 += 2;	
#ifdef TMJ_TIMING
				thrcnt2+=2;
#endif
			}
			pq1 += stride_y_m128 - iXN;
			rs1 += stride_y_m128 - iXN;
			Apq1 += stride_y_m128 - iXN;
			VelDen1 += stride_y_em - iXN*4;
		}

		for (int iY = 4;  iY > 0;  --iY)
		{
			pq0 += 2;
			VelDen0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(VelDen0), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelDen0+stride_z_em), _MM_HINT_NTA);

				__m128 ddy_0123;
				ISO_Cmp_DDY(pq0,stride_y_m128,ddy_0123);

				__m128 Boy;
				ISODenQ_Get_EM1(VelDen0, Boy);

				ISODenQ_Cmp_DYED(
						ddy_0123,  
						Boy,
						dyed_buf_0[0]
						);
				dyed_buf_0 += 2;

				pq0 += 1;
				VelDen0 += 4;

				// dydyed1, dydyed2
				__m128 dydyed1_0123;
				ISO_Cmp_DDY(dyed_buf_1,2*iXN,dydyed1_0123);
				V4_1[1] = _mm_add_ps(V4_1[1], dydyed1_0123);

				dyed_buf_1 += 2;
				V4_1 += 17;
			}

			pq0 += stride_y_m128 - (iXN+2);
			VelDen0 += stride_y_em - (iXN+2)*4;
		}
	}
	else // no more tiles to process, finish rest of outputs from queue
	{
		__m128* V4_0 = V4;
		for (int iY = iYN;  iY > 0;  --iY)
		{
			__m128 _mm_spgzy0 = _mm_mul_ps(_mm_spgz0, _mm_spgy[iY+4]);
			__m128 _mm_spgzy1 = _mm_mul_ps(_mm_spgz1, _mm_spgy[iY+4]);
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(VelDen1), _MM_HINT_NTA);
				_mm_prefetch((char*)(VelDen1+stride_z_em), _MM_HINT_NTA);
				if (stage == 0)
				{
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m128), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}

				// dzdzed1, dzdzed2
				// use shorter two point stencil
				__m128 V5_0123 = V4_0[5];

				__m128 dzdzed2_0123;
				ISO_Cmp_DDZ_EE(V4_0+12,-1,dzdzed2_0123);
				__m128 V4_0123 = dzdzed2_0123;

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				ISODenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelDen1,spgzyx0,V4_0123,V5_0123
						);

				V5_0123 = V4_0[ 4];

				ISO_Cmp_DDZ_EE(V4_0+11,-1,dzdzed2_0123);
				V4_0123 = dzdzed2_0123;

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				ISODenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m128,rs1+stride_z_m128,Apq1+stride_z_m128,VelDen1+stride_z_em,spgzyx1,V4_0123,V5_0123
						);

				V4_0 += 17;

				pq1 += 1;
				rs1 += 1;
				Apq1 += 1;
				VelDen1 += 4;
#ifdef TMJ_TIMING
				thrcnt2 += 2;
#endif
			}
			pq1 += stride_y_m128 - iXN;
			rs1 += stride_y_m128 - iXN;
			Apq1 += stride_y_m128 - iXN;
			VelDen1 += stride_y_em - iXN*4;
		}
	}
}

// 
// Fill up the vertical queues prior to looping over Z.
// This routine actually processes one Z more than it needs to.
// This is faster (proven experimentally) than adding extra logic to prevent it.
//
void VTIDenQ_Process_Patch_Leadin(
	int logLevel,
	__m128* pq,		// pq for stage 0 and 1, Apq for stage 2
	int* DenAng,
	int* VelAnis,
	__m128* stripe,
	__m128* dyed_buf,
	__m128* V4,
	int iX0,
	int iY0,
	int iZ,
	int iXN_halo,
	int iXN,
	int iYN_halo,
	int iYN,
	int stride_y_m128,
	int stride_z_m128,
	int stride_y_em,
	int stride_z_em
	)
{
	int abs_iZ0 = (iZ >= 0) ? iZ : -(iZ+1);
	int abs_iZ1 = ((iZ+1) >= 0) ? (iZ+1) : -(iZ+2);

	__m128* pq0 = pq + (unsigned long)abs_iZ0 * (unsigned long)stride_z_m128 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m128 + (unsigned long)(iX0-4);
	int* DenAng0 = DenAng + (unsigned long)abs_iZ0 * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

	int abs_stride_z_m128 = (abs_iZ1 - abs_iZ0) * stride_z_m128;
	int abs_stride_z_em = (abs_iZ1 - abs_iZ0) * stride_z_em;
	if (logLevel >= 6)
	{
		printf("iZ = %d, abs_iZ0 = %d, abs_iZ1 = %d, abs_stride_z_m128 = %d, abs_stride_em = %d\n",iZ,abs_iZ0,abs_iZ1,abs_stride_z_m128,abs_stride_z_em);
		fflush(stdout);
	}

	__m128* V4_0 = V4;
	__m128* V4_1 = V4;

	__m128* dyed_buf_0 = dyed_buf;
	__m128* dyed_buf_1 = dyed_buf + 16 * iXN;

	for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
	{
		pq0 += 4;
		DenAng0 += 8;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m128 ddy_01, ddy_23;
			Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);

			__m128 Boy;
			VTIDenQ_Get_EM1(DenAng0, Boy);

			VTIDenQ_Cmp_DYED(
					ddy_01, ddy_23,  
					Boy,
					dyed_buf_0[0], dyed_buf_0[1]
					);

			Cmp_DDY(pq0+abs_stride_z_m128,stride_y_m128,ddy_01,ddy_23);

			VTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy);

			VTIDenQ_Cmp_DYED(
					ddy_01, ddy_23,  
					Boy,
					dyed_buf_0[2], dyed_buf_0[3]
					);
			dyed_buf_0 += 4;

			pq0 += 2;
			DenAng0 += 4;
		}

		pq0 += stride_y_m128 - (iXN+2)*2;
		DenAng0 += stride_y_em - (iXN+2)*4;
	}

	for (int iY = (iYN_halo-5);  iY > 4;  --iY)
	{
		__m128* stripe0 = stripe;

		for (int iX = 2;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m128 ddx_01, ddx_23;
			Cmp_DDX(pq0,ddx_01,ddx_23);

			__m128 Boy;
			VTIDenQ_Get_EM1(DenAng0, Boy);

			VTIDenQ_Cmp_DXED(
					ddx_01, ddx_23,  
					Boy,
					stripe0[ 0], stripe0[ 1]
					);

			Cmp_DDX(pq0+abs_stride_z_m128,ddx_01,ddx_23);

			VTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy);

			VTIDenQ_Cmp_DXED(
					ddx_01, ddx_23,  
					Boy,
					stripe0[ 2], stripe0[ 3]
					);
			stripe0 += 4;

			pq0 += 2;
			DenAng0 += 4;
		}

		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23;
			Cmp_DDX(pq0,ddx_01,ddx_23);
			Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
			if (abs_iZ0 == 0)
				Cmp_DDZ_Z0(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 1)
				Cmp_DDZ_Z1(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 2)
				Cmp_DDZ_Z2(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 3)
				Cmp_DDZ_Z3(pq0,stride_z_m128,ddz_01,ddz_23);
			else
				Cmp_DDZ(pq0,stride_z_m128,ddz_01,ddz_23);

			__m128 Boy;
			VTIDenQ_Get_EM1(DenAng0, Boy);

			VTIDenQ_Cmp_DXED_DYED_DZED(
					ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
					Boy,
					stripe0[ 0], stripe0[ 1],
					dyed_buf_0[0], dyed_buf_0[1],
					V4_0[14], V4_0[15]
					);

			Cmp_DDX(pq0+abs_stride_z_m128,ddx_01,ddx_23);
			Cmp_DDY(pq0+abs_stride_z_m128,stride_y_m128,ddy_01,ddy_23);
			if (abs_iZ1 == 0)
				Cmp_DDZ_Z0(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 1)
				Cmp_DDZ_Z1(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 2)
				Cmp_DDZ_Z2(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 3)
				Cmp_DDZ_Z3(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else
				Cmp_DDZ(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);

			VTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy);

			VTIDenQ_Cmp_DXED_DYED_DZED(
					ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
					Boy,
					stripe0[ 2], stripe0[ 3],
					dyed_buf_0[2], dyed_buf_0[3],
					V4_0[12], V4_0[13]
					);
			stripe0 += 4;	
			dyed_buf_0 += 4;
			V4_0 += 34;

			pq0 += 2;
			DenAng0 += 4;
		}
		V4_0 -= 34 * iXN;

		for (int iX = 1;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m128 ddx_01, ddx_23;
			Cmp_DDX(pq0,ddx_01,ddx_23);

			__m128 Boy;
			VTIDenQ_Get_EM1(DenAng0, Boy);

			VTIDenQ_Cmp_DXED(
					ddx_01, ddx_23,   
					Boy,
					stripe0[ 0], stripe0[ 1]
					);

			Cmp_DDX(pq0+abs_stride_z_m128,ddx_01,ddx_23);

			VTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy);

			VTIDenQ_Cmp_DXED(
					ddx_01, ddx_23,  
					Boy,
					stripe0[ 2], stripe0[ 3]
					);
			stripe0 += 4;

			pq0 += 2;
			DenAng0 += 4;
		}

		pq0 += stride_y_m128 - iXN_halo*2;
		DenAng0 += stride_y_em - iXN_halo*4;

		stripe0 = stripe + 8;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			if (iY <= (iYN_halo-9))
			{
				// dydyed1, dydyed2
				__m128 dydyed1_01, dydyed1_23;
				Cmp_DDY(dyed_buf_1,4*iXN,dydyed1_01,dydyed1_23);
				V4_1[2] = _mm_add_ps(V4_1[2], dydyed1_01);
				V4_1[3] = _mm_add_ps(V4_1[3], dydyed1_23);

				Cmp_DDY(dyed_buf_1+2,4*iXN,dydyed1_01,dydyed1_23);
				V4_1[0] = _mm_add_ps(V4_1[0], dydyed1_01);
				V4_1[1] = _mm_add_ps(V4_1[1], dydyed1_23);

				dyed_buf_1 += 4;
				V4_1 += 34;
			}

			// dxdxed1, dxdxed2
			__m128 dxdxed1_01, dxdxed1_23;
			Cmp_DDX_1x_stride4(stripe0,dxdxed1_01,dxdxed1_23);
			V4_0[2] = dxdxed1_01;
			V4_0[3] = dxdxed1_23;

			Cmp_DDX_1x_stride4(stripe0+2,dxdxed1_01,dxdxed1_23);
			V4_0[0] = dxdxed1_01;
			V4_0[1] = dxdxed1_23;

			V4_0 += 34;

			stripe0 += 4;	
		}
	}

	for (int iY = 4;  iY > 0;  --iY)
	{
		pq0 += 4;
		DenAng0 += 8;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m128 ddy_01, ddy_23;
			Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);

			__m128 Boy;
			VTIDenQ_Get_EM1(DenAng0, Boy);

			VTIDenQ_Cmp_DYED(
					ddy_01, ddy_23,  
					Boy,
					dyed_buf_0[0], dyed_buf_0[1]
					);

			Cmp_DDY(pq0+abs_stride_z_m128,stride_y_m128,ddy_01,ddy_23);

			VTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy);

			VTIDenQ_Cmp_DYED(
					ddy_01, ddy_23,  
					Boy,
					dyed_buf_0[2], dyed_buf_0[3]
					);
			dyed_buf_0 += 4;

			pq0 += 2;
			DenAng0 += 4;

			// dydyed1, dydyed2
			__m128 dydyed1_01, dydyed1_23;
			Cmp_DDY(dyed_buf_1,4*iXN,dydyed1_01,dydyed1_23);
			V4_1[2] = _mm_add_ps(V4_1[2], dydyed1_01);
			V4_1[3] = _mm_add_ps(V4_1[3], dydyed1_23);

			Cmp_DDY(dyed_buf_1+2,4*iXN,dydyed1_01,dydyed1_23);
			V4_1[0] = _mm_add_ps(V4_1[0], dydyed1_01);
			V4_1[1] = _mm_add_ps(V4_1[1], dydyed1_23);

			dyed_buf_1 += 4;
			V4_1 += 34;
		}

		pq0 += stride_y_m128 - (iXN+2)*2;
		DenAng0 += stride_y_em - (iXN+2)*4;
	}
}

//
// Valid from [1,dimz-1>
//
void VTIDenQ_Process_Patch(
	int logLevel,
	int stage,
	__m128* pq,
	__m128* rs,
	__m128* Apq,			// pass pq ptr if 2nd order in time
	int* DenAng,
	int* VelAnis,
	__m128* _mm_spgx,
	__m128* _mm_spgy,
	__m128 _mm_spgz0,
	__m128 _mm_spgz1,
	__m128* stripe,
	__m128* dyed_buf,
	__m128* V4,
	int iX0,
	int iY0,
	int iZ,
	int dimz,
	int iXN_halo,
	int iXN,
	int iYN_halo,
	int iYN,
	int stride_y_m128,
	int stride_z_m128,
	int stride_y_em,
	int stride_z_em
#ifdef TMJ_TIMING
	,unsigned long& thrcnt2
#endif
	)
{
	int dz_edge_diff = dimz - iZ - 1;  // dimz-1 -> 0, dimz-2 -> 1 etc.

	__m128* pq1 = pq + (unsigned long)iZ * (unsigned long)stride_z_m128 + (unsigned long)iY0 * (unsigned long)stride_y_m128 + (unsigned long)iX0;
	__m128* rs1 = rs + (unsigned long)iZ * (unsigned long)stride_z_m128 + (unsigned long)iY0 * (unsigned long)stride_y_m128 + (unsigned long)iX0;
	__m128* Apq1 = Apq + (unsigned long)iZ * (unsigned long)stride_z_m128 + (unsigned long)iY0 * (unsigned long)stride_y_m128 + (unsigned long)iX0;
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
	 	__m128* pq0 = ((stage == 0 || stage == 1) ? pq : Apq) + (unsigned long)(iZ+4) * (unsigned long)stride_z_m128 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m128 + (unsigned long)(iX0-4);
		int* DenAng0 = DenAng + (unsigned long)(iZ+4) * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

		__m128* dyed_buf_0 = dyed_buf;
		__m128* dyed_buf_1 = dyed_buf + 16 * iXN;

		__m128* V4_0 = V4;
		__m128* V4_1 = V4;

		for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
		{
			pq0 += 4;
			DenAng0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddy_01, ddy_23;
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				VTIDenQ_Cmp_DYED(
						ddy_01, ddy_23,  
						Boy,
						dyed_buf_0[0], dyed_buf_0[1]
						);

				Cmp_DDY(pq0+stride_z_m128,stride_y_m128,ddy_01,ddy_23);

				VTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy);

				VTIDenQ_Cmp_DYED(
						ddy_01, ddy_23,  
						Boy,
                                                dyed_buf_0[2], dyed_buf_0[3]
						);

				dyed_buf_0 += 4;
				pq0 += 2;
				DenAng0 += 4;
			}
			pq0 += stride_y_m128 - (iXN+2)*2;
			DenAng0 += stride_y_em - (iXN+2)*4;
		}

		for (int iY = (iYN_halo-5);  iY > 4;  --iY)
		{
			__m128* stripe0 = stripe;

			for (int iX = 2;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddx_01, ddx_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				VTIDenQ_Cmp_DXED(
						ddx_01, ddx_23,  
						Boy,
						stripe0[ 0], stripe0[ 1]
						);

				Cmp_DDX(pq0+stride_z_m128,ddx_01,ddx_23);

				VTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy);

				VTIDenQ_Cmp_DXED(
						ddx_01, ddx_23,  
						Boy,
						stripe0[ 2], stripe0[ 3]
						);
				stripe0 += 4;

				pq0 += 2;
				DenAng0 += 4;
			}

			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
				if (dz_edge_diff > 9)
				{
					Cmp_DDZ(pq0,stride_z_m128,ddz_01,ddz_23);
				}
				else
				{
					Cmp_DDZ_EE(pq0,stride_z_m128,ddz_01,ddz_23);
				}

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				VTIDenQ_Cmp_DXED_DYED_DZED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
						Boy,
						stripe0[ 0], stripe0[ 1],
						dyed_buf_0[0], dyed_buf_0[1],
						V4_0[14], V4_0[15]
						);

				Cmp_DDX(pq0+stride_z_m128,ddx_01,ddx_23);
				Cmp_DDY(pq0+stride_z_m128,stride_y_m128,ddy_01,ddy_23);
				if (dz_edge_diff > 10)
				{
					Cmp_DDZ(pq0+stride_z_m128,stride_z_m128,ddz_01,ddz_23);
				}
				else
				{
					Cmp_DDZ_EE(pq0+stride_z_m128,stride_z_m128,ddz_01,ddz_23);
				}

				VTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy);

				VTIDenQ_Cmp_DXED_DYED_DZED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
						Boy,
						stripe0[ 2], stripe0[ 3],
                                                dyed_buf_0[2], dyed_buf_0[3],
						V4_0[12], V4_0[13]
						);
				stripe0 += 4;
				V4_0 += 34;

				dyed_buf_0 += 4;
				pq0 += 2;
				DenAng0 += 4;
			}
			V4_0 -= 34 * iXN;

			for (int iX = 1;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddx_01, ddx_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				VTIDenQ_Cmp_DXED(
						ddx_01, ddx_23,  
						Boy,
						stripe0[ 0], stripe0[ 1]
						);

				Cmp_DDX(pq0+stride_z_m128,ddx_01,ddx_23);

				VTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy);

				VTIDenQ_Cmp_DXED(
						ddx_01, ddx_23,  
						Boy,
						stripe0[ 2], stripe0[ 3]
						);
				stripe0 += 4;

				pq0 += 2;
				DenAng0 += 4;
			}

			pq0 += stride_y_m128 - iXN_halo*2;
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
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m128), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}

				if (iY <= (iYN_halo-9))
				{	
					// dydyed1, dydyed2
					__m128 dydyed1_01, dydyed1_23;
					Cmp_DDY(dyed_buf_1,4*iXN,dydyed1_01,dydyed1_23);
					V4_1[2] = _mm_add_ps(V4_1[2], dydyed1_01);
					V4_1[3] = _mm_add_ps(V4_1[3], dydyed1_23);

					Cmp_DDY(dyed_buf_1+2,4*iXN,dydyed1_01,dydyed1_23);
					V4_1[0] = _mm_add_ps(V4_1[0], dydyed1_01);
					V4_1[1] = _mm_add_ps(V4_1[1], dydyed1_23);

					dyed_buf_1 += 4;
					V4_1 += 34;
				}

				// dxdxed1, dxdxed2
				__m128 dxdxed1_01, dxdxed1_23;
				Cmp_DDX_1x_stride4(stripe0,dxdxed1_01,dxdxed1_23);
				V4_0[2] = dxdxed1_01;
				V4_0[3] = dxdxed1_23;

				Cmp_DDX_1x_stride4(stripe0+2,dxdxed1_01,dxdxed1_23);
				V4_0[0] = dxdxed1_01;
				V4_0[1] = dxdxed1_23;

				// dzdzed1, dzdzed2
				__m128 V5_01 = V4_0[10];
				__m128 V5_23 = V4_0[11];

				__m128 dzdzed2_01, dzdzed2_23;
				Cmp_DDZ(V4_0+24,-2,dzdzed2_01,dzdzed2_23);
				__m128 V4_01 = dzdzed2_01;
				__m128 V4_23 = dzdzed2_23;

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				TTIDenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelAnis1,DenAng1,spgzyx0,V4_01,V4_23,V5_01,V5_23
						);

				V5_01 = V4_0[ 8];
				V5_23 = V4_0[ 9];

				Cmp_DDZ(V4_0+22,-2,dzdzed2_01,dzdzed2_23);
				V4_01 = dzdzed2_01;
				V4_23 = dzdzed2_23;

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				TTIDenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m128,rs1+stride_z_m128,Apq1+stride_z_m128,VelAnis1+stride_z_em,DenAng1+stride_z_em,spgzyx1,V4_01,V4_23,V5_01,V5_23
						);

				V4_0 += 34;

				pq1 += 2;
				rs1 += 2;
				Apq1 += 2;
				VelAnis1 += 4;
				DenAng1 += 4;
				stripe0 += 4;	
#ifdef TMJ_TIMING
				thrcnt2+=2;
#endif
			}
			pq1 += stride_y_m128 - iXN*2;
			rs1 += stride_y_m128 - iXN*2;
			Apq1 += stride_y_m128 - iXN*2;
			VelAnis1 += stride_y_em - iXN*4;
			DenAng1 += stride_y_em - iXN*4;
		}

		for (int iY = 4;  iY > 0;  --iY)
		{
			pq0 += 4;
			DenAng0 += 8;	
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddy_01, ddy_23;
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				VTIDenQ_Cmp_DYED(
						ddy_01, ddy_23,  
						Boy,
						dyed_buf_0[0], dyed_buf_0[1]
						);

				Cmp_DDY(pq0+stride_z_m128,stride_y_m128,ddy_01,ddy_23);

				VTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy);

				VTIDenQ_Cmp_DYED(
						ddy_01, ddy_23,  
						Boy,
                                                dyed_buf_0[2], dyed_buf_0[3]
						);

				dyed_buf_0 += 4;
				pq0 += 2;
				DenAng0 += 4;

				// dydyed1, dydyed2
				__m128 dydyed1_01, dydyed1_23;
				Cmp_DDY(dyed_buf_1,4*iXN,dydyed1_01,dydyed1_23);
				V4_1[2] = _mm_add_ps(V4_1[2], dydyed1_01);
				V4_1[3] = _mm_add_ps(V4_1[3], dydyed1_23);

				Cmp_DDY(dyed_buf_1+2,4*iXN,dydyed1_01,dydyed1_23);
				V4_1[0] = _mm_add_ps(V4_1[0], dydyed1_01);
				V4_1[1] = _mm_add_ps(V4_1[1], dydyed1_23);

				dyed_buf_1 += 4;
				V4_1 += 34;
			}
	
			pq0 += stride_y_m128 - (iXN+2)*2;
			DenAng0 += stride_y_em - (iXN+2)*4;
		}
	}
	else if (dz_edge_diff == 5)
	{
	 	__m128* pq0 = ((stage == 0 || stage == 1) ? pq : Apq) + (unsigned long)(iZ+4) * (unsigned long)stride_z_m128 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m128 + (unsigned long)(iX0-4);
		int* DenAng0 = DenAng + (unsigned long)(iZ+4) * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

		__m128* dyed_buf_0 = dyed_buf;
		__m128* dyed_buf_1 = dyed_buf + 16 * iXN;

		__m128* V4_0 = V4;
		__m128* V4_1 = V4;

		for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
		{
			pq0 += 4;
			DenAng0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m128 ddy_01, ddy_23;
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				VTIDenQ_Cmp_DYED(
						ddy_01, ddy_23,  
						Boy,
						dyed_buf_0[0], dyed_buf_0[1]
						);
				dyed_buf_0 += 4;

				pq0 += 2;
				DenAng0 += 4;
			}

			pq0 += stride_y_m128 - (iXN+2)*2;
			DenAng0 += stride_y_em - (iXN+2)*4;
		}

		for (int iY = (iYN_halo-5);  iY > 4;  --iY)
		{
			__m128* stripe0 = stripe;

			for (int iX = 2;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m128 ddx_01, ddx_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				VTIDenQ_Cmp_DXED(
						ddx_01, ddx_23,  
						Boy,
						stripe0[ 0], stripe0[ 1]
						);
				stripe0 += 4;

				pq0 += 2;
				DenAng0 += 4;
			}

			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
				Cmp_DDZ_EE(pq0,stride_z_m128,ddz_01,ddz_23);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				VTIDenQ_Cmp_DXED_DYED_DZED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
						Boy,
						stripe0[ 0], stripe0[ 1],
						dyed_buf_0[0], dyed_buf_0[1],
						V4_0[14], V4_0[15]
						);
				stripe0 += 4;
				dyed_buf_0 += 4;
				V4_0 += 34;

				pq0 += 2;
				DenAng0 += 4;
			}
			V4_0 -= 34 * iXN;

			for (int iX = 1;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m128 ddx_01, ddx_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				VTIDenQ_Cmp_DXED(
						ddx_01, ddx_23,  
						Boy,
						stripe0[ 0], stripe0[ 1]
						);
				stripe0 += 4;

				pq0 += 2;
				DenAng0 += 4;
			}

			pq0 += stride_y_m128 - iXN_halo*2;
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
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m128), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}

				if (iY <= (iYN_halo-9))
				{
					// dydyed1, dydyed2
					__m128 dydyed1_01, dydyed1_23;
					Cmp_DDY(dyed_buf_1,4*iXN,dydyed1_01,dydyed1_23);
					V4_1[2] = _mm_add_ps(V4_1[2], dydyed1_01);
					V4_1[3] = _mm_add_ps(V4_1[3], dydyed1_23);

					dyed_buf_1 += 4;
					V4_1 += 34;
				}

				// dxdxed1, dxdxed2
				__m128 dxdxed1_01, dxdxed1_23;
				Cmp_DDX_1x_stride4(stripe0,dxdxed1_01,dxdxed1_23);
				V4_0[2] = dxdxed1_01;
				V4_0[3] = dxdxed1_23;

				// dzdzed1, dzdzed2
				__m128 V5_01 = V4_0[10];
				__m128 V5_23 = V4_0[11];

				__m128 dzdzed2_01, dzdzed2_23;
				Cmp_DDZ(V4_0+24,-2,dzdzed2_01,dzdzed2_23);
				__m128 V4_01 = dzdzed2_01;
				__m128 V4_23 = dzdzed2_23;

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				TTIDenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelAnis1,DenAng1,spgzyx0,V4_01,V4_23,V5_01,V5_23
						);

				V5_01 = V4_0[ 8];
				V5_23 = V4_0[ 9];

				Cmp_DDZ_EE(V4_0+22,-2,dzdzed2_01,dzdzed2_23);
				V4_01 = dzdzed2_01;
				V4_23 = dzdzed2_23;

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				TTIDenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m128,rs1+stride_z_m128,Apq1+stride_z_m128,VelAnis1+stride_z_em,DenAng1+stride_z_em,spgzyx1,V4_01,V4_23,V5_01,V5_23
						);

				V4_0 += 34;

				pq1 += 2;
				rs1 += 2;
				Apq1 += 2;
				VelAnis1 += 4;
				DenAng1 += 4;
				stripe0 += 4;	
#ifdef TMJ_TIMING
				thrcnt2+=2;
#endif
			}
			pq1 += stride_y_m128 - iXN*2;
			rs1 += stride_y_m128 - iXN*2;
			Apq1 += stride_y_m128 - iXN*2;
			VelAnis1 += stride_y_em - iXN*4;
			DenAng1 += stride_y_em - iXN*4;
		}

		for (int iY = 4;  iY > 0;  --iY)
		{
			pq0 += 4;
			DenAng0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m128 ddy_01, ddy_23;
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);

				__m128 Boy;
				VTIDenQ_Get_EM1(DenAng0, Boy);

				VTIDenQ_Cmp_DYED(
						ddy_01, ddy_23,  
						Boy,
						dyed_buf_0[0], dyed_buf_0[1]
						);
				dyed_buf_0 += 4;

				pq0 += 2;
				DenAng0 += 4;

				// dydyed1, dydyed2
				__m128 dydyed1_01, dydyed1_23;
				Cmp_DDY(dyed_buf_1,4*iXN,dydyed1_01,dydyed1_23);
				V4_1[2] = _mm_add_ps(V4_1[2], dydyed1_01);
				V4_1[3] = _mm_add_ps(V4_1[3], dydyed1_23);

				dyed_buf_1 += 4;
				V4_1 += 34;
			}

			pq0 += stride_y_m128 - (iXN+2)*2;
			DenAng0 += stride_y_em - (iXN+2)*4;
		}
	}
	else // no more tiles to process, finish rest of outputs from queue
	{
		__m128* V4_0 = V4;
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
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m128), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}

				// dzdzed1, dzdzed2
				// use shorter two point stencil
				__m128 V5_01 = V4_0[10];
				__m128 V5_23 = V4_0[11];

				__m128 dzdzed2_01, dzdzed2_23;
				Cmp_DDZ_EE(V4_0+24,-2,dzdzed2_01,dzdzed2_23);
				__m128 V4_01 = dzdzed2_01;
				__m128 V4_23 = dzdzed2_23;

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				TTIDenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelAnis1,DenAng1,spgzyx0,V4_01,V4_23,V5_01,V5_23
						);

				V5_01 = V4_0[ 8];
				V5_23 = V4_0[ 9];

				Cmp_DDZ_EE(V4_0+22,-2,dzdzed2_01,dzdzed2_23);
				V4_01 = dzdzed2_01;
				V4_23 = dzdzed2_23;

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				TTIDenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m128,rs1+stride_z_m128,Apq1+stride_z_m128,VelAnis1+stride_z_em,DenAng1+stride_z_em,spgzyx1,V4_01,V4_23,V5_01,V5_23
						);

				V4_0 += 34;

				pq1 += 2;
				rs1 += 2;
				Apq1 += 2;
				VelAnis1 += 4;
				DenAng1 += 4;
#ifdef TMJ_TIMING
				thrcnt2 += 2;
#endif
			}
			pq1 += stride_y_m128 - iXN*2;
			rs1 += stride_y_m128 - iXN*2;
			Apq1 += stride_y_m128 - iXN*2;
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
void TTIDenQ_Process_Patch_Leadin(
	int logLevel,
	__m128* pq,		// pq for stage 0 and 1, Apq for stage 2
	int* DenAng,
	int* VelAnis,
	__m128* stripe,
	__m128* dyed_buf,
	__m128* V4,
	int iX0,
	int iY0,
	int iZ,
	int iXN_halo,
	int iXN,
	int iYN_halo,
	int iYN,
	int stride_y_m128,
	int stride_z_m128,
	int stride_y_em,
	int stride_z_em
	)
{
	int abs_iZ0 = (iZ >= 0) ? iZ : -(iZ+1);
	int abs_iZ1 = ((iZ+1) >= 0) ? (iZ+1) : -(iZ+2);

	__m128* pq0 = pq + (unsigned long)abs_iZ0 * (unsigned long)stride_z_m128 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m128 + (unsigned long)(iX0-4);
	int* DenAng0 = DenAng + (unsigned long)abs_iZ0 * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

	int abs_stride_z_m128 = (abs_iZ1 - abs_iZ0) * stride_z_m128;
	int abs_stride_z_em = (abs_iZ1 - abs_iZ0) * stride_z_em;
	if (logLevel >= 6)
	{
		printf("iZ = %d, abs_iZ0 = %d, abs_iZ1 = %d, abs_stride_z_m128 = %d, abs_stride_em = %d\n",iZ,abs_iZ0,abs_iZ1,abs_stride_z_m128,abs_stride_z_em);
		fflush(stdout);
	}

	__m128* V4_0 = V4;
	__m128* V4_1 = V4;

	__m128* dyed_buf_0 = dyed_buf;
	__m128* dyed_buf_1 = dyed_buf + 32 * iXN;

	for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
	{
		pq0 += 4;
		DenAng0 += 8;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23;
			Cmp_DDX(pq0,ddx_01,ddx_23);
			Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
			if (abs_iZ0 == 0)
				Cmp_DDZ_Z0(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 1)
				Cmp_DDZ_Z1(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 2)
				Cmp_DDZ_Z2(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 3)
				Cmp_DDZ_Z3(pq0,stride_z_m128,ddz_01,ddz_23);
			else
				Cmp_DDZ(pq0,stride_z_m128,ddz_01,ddz_23);

			__m128 Boy, cDip, sDip, cAzm, sAzm;
			TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

			TTIDenQ_Cmp_DYED(
					ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
					Boy, cDip, sDip, cAzm, sAzm,
					dyed_buf_0[0], dyed_buf_0[1], dyed_buf_0[2], dyed_buf_0[3]
					);

			Cmp_DDX(pq0+abs_stride_z_m128,ddx_01,ddx_23);
			Cmp_DDY(pq0+abs_stride_z_m128,stride_y_m128,ddy_01,ddy_23);
			if (abs_iZ1 == 0)
				Cmp_DDZ_Z0(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 1)
				Cmp_DDZ_Z1(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 2)
				Cmp_DDZ_Z2(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 3)
				Cmp_DDZ_Z3(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else
				Cmp_DDZ(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);

			TTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

			TTIDenQ_Cmp_DYED(
					ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
					Boy, cDip, sDip, cAzm, sAzm,
					dyed_buf_0[4], dyed_buf_0[5], dyed_buf_0[6], dyed_buf_0[7]
					);
			dyed_buf_0 += 8;

			pq0 += 2;
			DenAng0 += 4;
		}

		pq0 += stride_y_m128 - (iXN+2)*2;
		DenAng0 += stride_y_em - (iXN+2)*4;
	}

	for (int iY = (iYN_halo-5);  iY > 4;  --iY)
	{
		__m128* stripe0 = stripe;

		for (int iX = 2;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23;
			Cmp_DDX(pq0,ddx_01,ddx_23);
			Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
			if (abs_iZ0 == 0)
				Cmp_DDZ_Z0(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 1)
				Cmp_DDZ_Z1(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 2)
				Cmp_DDZ_Z2(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 3)
				Cmp_DDZ_Z3(pq0,stride_z_m128,ddz_01,ddz_23);
			else
				Cmp_DDZ(pq0,stride_z_m128,ddz_01,ddz_23);

			__m128 Boy, cDip, sDip, cAzm, sAzm;
			TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

			TTIDenQ_Cmp_DXED(
					ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
					Boy, cDip, sDip, cAzm, sAzm,
					stripe0[ 0], stripe0[ 1], stripe0[ 2], stripe0[ 3]
					);

			Cmp_DDX(pq0+abs_stride_z_m128,ddx_01,ddx_23);
			Cmp_DDY(pq0+abs_stride_z_m128,stride_y_m128,ddy_01,ddy_23);
			if (abs_iZ1 == 0)
				Cmp_DDZ_Z0(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 1)
				Cmp_DDZ_Z1(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 2)
				Cmp_DDZ_Z2(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 3)
				Cmp_DDZ_Z3(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else
				Cmp_DDZ(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);

			TTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

			TTIDenQ_Cmp_DXED(
					ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
					Boy, cDip, sDip, cAzm, sAzm,
					stripe0[ 4], stripe0[ 5], stripe0[ 6], stripe0[ 7]
					);
			stripe0 += 8;

			pq0 += 2;
			DenAng0 += 4;
		}

		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23;
			Cmp_DDX(pq0,ddx_01,ddx_23);
			Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
			if (abs_iZ0 == 0)
				Cmp_DDZ_Z0(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 1)
				Cmp_DDZ_Z1(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 2)
				Cmp_DDZ_Z2(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 3)
				Cmp_DDZ_Z3(pq0,stride_z_m128,ddz_01,ddz_23);
			else
				Cmp_DDZ(pq0,stride_z_m128,ddz_01,ddz_23);

			__m128 Boy, cDip, sDip, cAzm, sAzm;
			TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

			TTIDenQ_Cmp_DXED_DYED_DZED(
					ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
					Boy, cDip, sDip, cAzm, sAzm,
					stripe0[ 0], stripe0[ 1], stripe0[ 2], stripe0[ 3],
					dyed_buf_0[0], dyed_buf_0[1], dyed_buf_0[2], dyed_buf_0[3],
					V4_0[28], V4_0[29], V4_0[30], V4_0[31]
					);

			Cmp_DDX(pq0+abs_stride_z_m128,ddx_01,ddx_23);
			Cmp_DDY(pq0+abs_stride_z_m128,stride_y_m128,ddy_01,ddy_23);
			if (abs_iZ1 == 0)
				Cmp_DDZ_Z0(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 1)
				Cmp_DDZ_Z1(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 2)
				Cmp_DDZ_Z2(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 3)
				Cmp_DDZ_Z3(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else
				Cmp_DDZ(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);

			TTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

			TTIDenQ_Cmp_DXED_DYED_DZED(
					ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
					Boy, cDip, sDip, cAzm, sAzm,
					stripe0[ 4], stripe0[ 5], stripe0[ 6], stripe0[ 7],
					dyed_buf_0[4], dyed_buf_0[5], dyed_buf_0[6], dyed_buf_0[7],
					V4_0[24], V4_0[25], V4_0[26], V4_0[27]
					);
			stripe0 += 8;	
			dyed_buf_0 += 8;
			V4_0 += 68;

			pq0 += 2;
			DenAng0 += 4;
		}
		V4_0 -= 68 * iXN;

		for (int iX = 1;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23;
			Cmp_DDX(pq0,ddx_01,ddx_23);
			Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
			if (abs_iZ0 == 0)
				Cmp_DDZ_Z0(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 1)
				Cmp_DDZ_Z1(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 2)
				Cmp_DDZ_Z2(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 3)
				Cmp_DDZ_Z3(pq0,stride_z_m128,ddz_01,ddz_23);
			else
				Cmp_DDZ(pq0,stride_z_m128,ddz_01,ddz_23);

			__m128 Boy, cDip, sDip, cAzm, sAzm;
			TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

			TTIDenQ_Cmp_DXED(
					ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
					Boy, cDip, sDip, cAzm, sAzm,
					stripe0[ 0], stripe0[ 1], stripe0[ 2], stripe0[ 3]
					);

			Cmp_DDX(pq0+abs_stride_z_m128,ddx_01,ddx_23);
			Cmp_DDY(pq0+abs_stride_z_m128,stride_y_m128,ddy_01,ddy_23);
			if (abs_iZ1 == 0)
				Cmp_DDZ_Z0(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 1)
				Cmp_DDZ_Z1(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 2)
				Cmp_DDZ_Z2(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 3)
				Cmp_DDZ_Z3(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else
				Cmp_DDZ(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);

			TTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

			TTIDenQ_Cmp_DXED(
					ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
					Boy, cDip, sDip, cAzm, sAzm,
					stripe0[ 4], stripe0[ 5], stripe0[ 6], stripe0[ 7]
					);
			stripe0 += 8;

			pq0 += 2;
			DenAng0 += 4;
		}

		pq0 += stride_y_m128 - iXN_halo*2;
		DenAng0 += stride_y_em - iXN_halo*4;

		stripe0 = stripe + 16;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			if (iY <= (iYN_halo-9))
			{
				// dydyed1, dydyed2
				__m128 dydyed1_01, dydyed1_23;
				Cmp_DDY(dyed_buf_1,8*iXN,dydyed1_01,dydyed1_23);
				V4_1[4] = _mm_add_ps(V4_1[4], dydyed1_01);
				V4_1[5] = _mm_add_ps(V4_1[5], dydyed1_23);

				__m128 dydyed2_01, dydyed2_23;
				Cmp_DDY(dyed_buf_1+2,8*iXN,dydyed2_01,dydyed2_23);
				V4_1[6] = _mm_add_ps(V4_1[6], dydyed2_01);
				V4_1[7] = _mm_add_ps(V4_1[7], dydyed2_23);

				Cmp_DDY(dyed_buf_1+4,8*iXN,dydyed1_01,dydyed1_23);
				V4_1[0] = _mm_add_ps(V4_1[0], dydyed1_01);
				V4_1[1] = _mm_add_ps(V4_1[1], dydyed1_23);

				Cmp_DDY(dyed_buf_1+6,8*iXN,dydyed2_01,dydyed2_23);
				V4_1[2] = _mm_add_ps(V4_1[2], dydyed2_01);
				V4_1[3] = _mm_add_ps(V4_1[3], dydyed2_23);
				dyed_buf_1 += 8;
				V4_1 += 68;
			}

			// dxdxed1, dxdxed2
			__m128 dxdxed1_01, dxdxed1_23;
			Cmp_DDX_1x_stride8(stripe0,dxdxed1_01,dxdxed1_23);
			V4_0[4] = dxdxed1_01;
			V4_0[5] = dxdxed1_23;

			__m128 dxdxed2_01, dxdxed2_23;
			Cmp_DDX_1x_stride8(stripe0+2,dxdxed2_01,dxdxed2_23);
			V4_0[6] = dxdxed2_01;
			V4_0[7] = dxdxed2_23;

			Cmp_DDX_1x_stride8(stripe0+4,dxdxed1_01,dxdxed1_23);
			V4_0[0] = dxdxed1_01;
			V4_0[1] = dxdxed1_23;

			Cmp_DDX_1x_stride8(stripe0+6,dxdxed2_01,dxdxed2_23);
			V4_0[2] = dxdxed2_01;
			V4_0[3] = dxdxed2_23;

			V4_0 += 68;

			stripe0 += 8;	
		}
	}

	for (int iY = 4;  iY > 0;  --iY)
	{
		pq0 += 4;
		DenAng0 += 8;
		for (int iX = iXN;  iX > 0;  --iX)
		{
			_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
			_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

			__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23;
			Cmp_DDX(pq0,ddx_01,ddx_23);
			Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
			if (abs_iZ0 == 0)
				Cmp_DDZ_Z0(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 1)
				Cmp_DDZ_Z1(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 2)
				Cmp_DDZ_Z2(pq0,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ0 == 3)
				Cmp_DDZ_Z3(pq0,stride_z_m128,ddz_01,ddz_23);
			else
				Cmp_DDZ(pq0,stride_z_m128,ddz_01,ddz_23);

			__m128 Boy, cDip, sDip, cAzm, sAzm;
			TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

			TTIDenQ_Cmp_DYED(
					ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
					Boy, cDip, sDip, cAzm, sAzm,
					dyed_buf_0[0], dyed_buf_0[1], dyed_buf_0[2], dyed_buf_0[3]
					);

			Cmp_DDX(pq0+abs_stride_z_m128,ddx_01,ddx_23);
			Cmp_DDY(pq0+abs_stride_z_m128,stride_y_m128,ddy_01,ddy_23);
			if (abs_iZ1 == 0)
				Cmp_DDZ_Z0(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 1)
				Cmp_DDZ_Z1(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 2)
				Cmp_DDZ_Z2(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else if (abs_iZ1 == 3)
				Cmp_DDZ_Z3(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);
			else
				Cmp_DDZ(pq0+abs_stride_z_m128,stride_z_m128,ddz_01,ddz_23);

			TTIDenQ_Get_EM1(DenAng0+abs_stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

			TTIDenQ_Cmp_DYED(
					ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
					Boy, cDip, sDip, cAzm, sAzm,
					dyed_buf_0[4], dyed_buf_0[5], dyed_buf_0[6], dyed_buf_0[7]
					);
			dyed_buf_0 += 8;

			pq0 += 2;
			DenAng0 += 4;

			// dydyed1, dydyed2
			__m128 dydyed1_01, dydyed1_23;
			Cmp_DDY(dyed_buf_1,8*iXN,dydyed1_01,dydyed1_23);
			V4_1[4] = _mm_add_ps(V4_1[4], dydyed1_01);
			V4_1[5] = _mm_add_ps(V4_1[5], dydyed1_23);

			__m128 dydyed2_01, dydyed2_23;
			Cmp_DDY(dyed_buf_1+2,8*iXN,dydyed2_01,dydyed2_23);
			V4_1[6] = _mm_add_ps(V4_1[6], dydyed2_01);
			V4_1[7] = _mm_add_ps(V4_1[7], dydyed2_23);

			Cmp_DDY(dyed_buf_1+4,8*iXN,dydyed1_01,dydyed1_23);
			V4_1[0] = _mm_add_ps(V4_1[0], dydyed1_01);
			V4_1[1] = _mm_add_ps(V4_1[1], dydyed1_23);

			Cmp_DDY(dyed_buf_1+6,8*iXN,dydyed2_01,dydyed2_23);
			V4_1[2] = _mm_add_ps(V4_1[2], dydyed2_01);
			V4_1[3] = _mm_add_ps(V4_1[3], dydyed2_23);
			dyed_buf_1 += 8;
			V4_1 += 68;
		}

		pq0 += stride_y_m128 - (iXN+2)*2;
		DenAng0 += stride_y_em - (iXN+2)*4;
	}
}

//
// Valid from [1,dimz-1>
//
void TTIDenQ_Process_Patch(
	int logLevel,
	int stage,
	__m128* pq,
	__m128* rs,
	__m128* Apq,			// pass pq ptr if 2nd order in time
	int* DenAng,
	int* VelAnis,
	__m128* _mm_spgx,
	__m128* _mm_spgy,
	__m128 _mm_spgz0,
	__m128 _mm_spgz1,
	__m128* stripe,
	__m128* dyed_buf,
	__m128* V4,
	int iX0,
	int iY0,
	int iZ,
	int dimz,
	int iXN_halo,
	int iXN,
	int iYN_halo,
	int iYN,
	int stride_y_m128,
	int stride_z_m128,
	int stride_y_em,
	int stride_z_em
#ifdef TMJ_TIMING
	,unsigned long& thrcnt2
#endif
	)
{
	int dz_edge_diff = dimz - iZ - 1;  // dimz-1 -> 0, dimz-2 -> 1 etc.

	__m128* pq1 = pq + (unsigned long)iZ * (unsigned long)stride_z_m128 + (unsigned long)iY0 * (unsigned long)stride_y_m128 + (unsigned long)iX0;
	__m128* rs1 = rs + (unsigned long)iZ * (unsigned long)stride_z_m128 + (unsigned long)iY0 * (unsigned long)stride_y_m128 + (unsigned long)iX0;
	__m128* Apq1 = Apq + (unsigned long)iZ * (unsigned long)stride_z_m128 + (unsigned long)iY0 * (unsigned long)stride_y_m128 + (unsigned long)iX0;
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
	 	__m128* pq0 = ((stage == 0 || stage == 1) ? pq : Apq) + (unsigned long)(iZ+4) * (unsigned long)stride_z_m128 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m128 + (unsigned long)(iX0-4);
		int* DenAng0 = DenAng + (unsigned long)(iZ+4) * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

		__m128* dyed_buf_0 = dyed_buf;
		__m128* dyed_buf_1 = dyed_buf + 32 * iXN;

		__m128* V4_0 = V4;
		__m128* V4_1 = V4;

		for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
		{
			pq0 += 4;
			DenAng0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_z0_01, ddz_z0_23, ddz_z1_01, ddz_z1_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
				if (dz_edge_diff > 10)
				{
					Cmp_DDZ_2Z(pq0,stride_z_m128,ddz_z0_01,ddz_z0_23,ddz_z1_01,ddz_z1_23);
				}
				else if (dz_edge_diff > 9)
				{
					Cmp_DDZ(pq0,stride_z_m128,ddz_z0_01,ddz_z0_23);
					Cmp_DDZ_EE(pq0+stride_z_m128,stride_z_m128,ddz_z1_01,ddz_z1_23);
				}
				else
				{
					Cmp_DDZ_EE(pq0,stride_z_m128,ddz_z0_01,ddz_z0_23);
					Cmp_DDZ_EE(pq0+stride_z_m128,stride_z_m128,ddz_z1_01,ddz_z1_23);
				}

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DYED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_z0_01, ddz_z0_23, 
						Boy, cDip, sDip, cAzm, sAzm,
						dyed_buf_0[0], dyed_buf_0[1], dyed_buf_0[2], dyed_buf_0[3]
						);

				Cmp_DDX(pq0+stride_z_m128,ddx_01,ddx_23);
				Cmp_DDY(pq0+stride_z_m128,stride_y_m128,ddy_01,ddy_23);

				TTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DYED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_z1_01, ddz_z1_23, 
						Boy, cDip, sDip, cAzm, sAzm,
                                                dyed_buf_0[4], dyed_buf_0[5], dyed_buf_0[6], dyed_buf_0[7]
						);

				dyed_buf_0 += 8;
				pq0 += 2;
				DenAng0 += 4;
			}
			pq0 += stride_y_m128 - (iXN+2)*2;
			DenAng0 += stride_y_em - (iXN+2)*4;
		}

		for (int iY = (iYN_halo-5);  iY > 4;  --iY)
		{
			__m128* stripe0 = stripe;

			for (int iX = 2;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_z0_01, ddz_z0_23, ddz_z1_01, ddz_z1_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
				if (dz_edge_diff > 10)
				{
					Cmp_DDZ_2Z(pq0,stride_z_m128,ddz_z0_01,ddz_z0_23,ddz_z1_01,ddz_z1_23);
				}
				else if (dz_edge_diff > 9)
				{
					Cmp_DDZ(pq0,stride_z_m128,ddz_z0_01,ddz_z0_23);
					Cmp_DDZ_EE(pq0+stride_z_m128,stride_z_m128,ddz_z1_01,ddz_z1_23);
				}
				else
				{
					Cmp_DDZ_EE(pq0,stride_z_m128,ddz_z0_01,ddz_z0_23);
					Cmp_DDZ_EE(pq0+stride_z_m128,stride_z_m128,ddz_z1_01,ddz_z1_23);
				}

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DXED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_z0_01, ddz_z0_23, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 0], stripe0[ 1], stripe0[ 2], stripe0[ 3]
						);

				Cmp_DDX(pq0+stride_z_m128,ddx_01,ddx_23);
				Cmp_DDY(pq0+stride_z_m128,stride_y_m128,ddy_01,ddy_23);

				TTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DXED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_z1_01, ddz_z1_23, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 4], stripe0[ 5], stripe0[ 6], stripe0[ 7]
						);
				stripe0 += 8;

				pq0 += 2;
				DenAng0 += 4;
			}

			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_z0_01, ddz_z0_23, ddz_z1_01, ddz_z1_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
				if (dz_edge_diff > 10)
				{
					Cmp_DDZ_2Z(pq0,stride_z_m128,ddz_z0_01,ddz_z0_23,ddz_z1_01,ddz_z1_23);
				}
				else if (dz_edge_diff > 9)
				{
					Cmp_DDZ(pq0,stride_z_m128,ddz_z0_01,ddz_z0_23);
					Cmp_DDZ_EE(pq0+stride_z_m128,stride_z_m128,ddz_z1_01,ddz_z1_23);
				}
				else
				{
					Cmp_DDZ_EE(pq0,stride_z_m128,ddz_z0_01,ddz_z0_23);
					Cmp_DDZ_EE(pq0+stride_z_m128,stride_z_m128,ddz_z1_01,ddz_z1_23);
				}

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DXED_DYED_DZED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_z0_01, ddz_z0_23, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 0], stripe0[ 1], stripe0[ 2], stripe0[ 3],
						dyed_buf_0[0], dyed_buf_0[1], dyed_buf_0[2], dyed_buf_0[3],
						V4_0[28], V4_0[29], V4_0[30], V4_0[31]
						);

				Cmp_DDX(pq0+stride_z_m128,ddx_01,ddx_23);
				Cmp_DDY(pq0+stride_z_m128,stride_y_m128,ddy_01,ddy_23);

				TTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DXED_DYED_DZED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_z1_01, ddz_z1_23, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 4], stripe0[ 5], stripe0[ 6], stripe0[ 7],
                                                dyed_buf_0[4], dyed_buf_0[5], dyed_buf_0[6], dyed_buf_0[7],
						V4_0[24], V4_0[25], V4_0[26], V4_0[27]
						);
				stripe0 += 8;
				V4_0 += 68;

				dyed_buf_0 += 8;
				pq0 += 2;
				DenAng0 += 4;
			}
			V4_0 -= 68 * iXN;

			for (int iX = 1;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_z0_01, ddz_z0_23, ddz_z1_01, ddz_z1_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
				if (dz_edge_diff > 10)
				{
					Cmp_DDZ_2Z(pq0,stride_z_m128,ddz_z0_01,ddz_z0_23,ddz_z1_01,ddz_z1_23);
				}
				else if (dz_edge_diff > 9)
				{
					Cmp_DDZ(pq0,stride_z_m128,ddz_z0_01,ddz_z0_23);
					Cmp_DDZ_EE(pq0+stride_z_m128,stride_z_m128,ddz_z1_01,ddz_z1_23);
				}
				else
				{
					Cmp_DDZ_EE(pq0,stride_z_m128,ddz_z0_01,ddz_z0_23);
					Cmp_DDZ_EE(pq0+stride_z_m128,stride_z_m128,ddz_z1_01,ddz_z1_23);
				}

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DXED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_z0_01, ddz_z0_23, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 0], stripe0[ 1], stripe0[ 2], stripe0[ 3]
						);

				Cmp_DDX(pq0+stride_z_m128,ddx_01,ddx_23);
				Cmp_DDY(pq0+stride_z_m128,stride_y_m128,ddy_01,ddy_23);

				TTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DXED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_z1_01, ddz_z1_23, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 4], stripe0[ 5], stripe0[ 6], stripe0[ 7]
						);
				stripe0 += 8;

				pq0 += 2;
				DenAng0 += 4;
			}

			pq0 += stride_y_m128 - iXN_halo*2;
			DenAng0 += stride_y_em - iXN_halo*4;

			stripe0 = stripe + 16;
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
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m128), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}

				if (iY <= (iYN_halo-9))
				{	
					// dydyed1, dydyed2
					__m128 dydyed1_01, dydyed1_23;
					Cmp_DDY(dyed_buf_1,8*iXN,dydyed1_01,dydyed1_23);
					V4_1[4] = _mm_add_ps(V4_1[4], dydyed1_01);
					V4_1[5] = _mm_add_ps(V4_1[5], dydyed1_23);

					__m128 dydyed2_01, dydyed2_23;
					Cmp_DDY(dyed_buf_1+2,8*iXN,dydyed2_01,dydyed2_23);
					V4_1[6] = _mm_add_ps(V4_1[6], dydyed2_01);
					V4_1[7] = _mm_add_ps(V4_1[7], dydyed2_23);

					Cmp_DDY(dyed_buf_1+4,8*iXN,dydyed1_01,dydyed1_23);
					V4_1[0] = _mm_add_ps(V4_1[0], dydyed1_01);
					V4_1[1] = _mm_add_ps(V4_1[1], dydyed1_23);

					Cmp_DDY(dyed_buf_1+6,8*iXN,dydyed2_01,dydyed2_23);
					V4_1[2] = _mm_add_ps(V4_1[2], dydyed2_01);
					V4_1[3] = _mm_add_ps(V4_1[3], dydyed2_23);
					dyed_buf_1 += 8;
					V4_1 += 68;
				}

				// dxdxed1, dxdxed2
				__m128 dxdxed1_01, dxdxed1_23;
				Cmp_DDX_1x_stride8(stripe0,dxdxed1_01,dxdxed1_23);
				V4_0[4] = dxdxed1_01;
				V4_0[5] = dxdxed1_23;

				__m128 dxdxed2_01, dxdxed2_23;
				Cmp_DDX_1x_stride8(stripe0+2,dxdxed2_01,dxdxed2_23);
				V4_0[6] = dxdxed2_01;
				V4_0[7] = dxdxed2_23;

				Cmp_DDX_1x_stride8(stripe0+4,dxdxed1_01,dxdxed1_23);
				V4_0[0] = dxdxed1_01;
				V4_0[1] = dxdxed1_23;

				Cmp_DDX_1x_stride8(stripe0+6,dxdxed2_01,dxdxed2_23);
				V4_0[2] = dxdxed2_01;
				V4_0[3] = dxdxed2_23;

				// dzdzed1, dzdzed2
				__m128 dzdzed1_z0_01, dzdzed1_z0_23, dzdzed1_z1_01, dzdzed1_z1_23;
				__m128 dzdzed2_z0_01, dzdzed2_z0_23, dzdzed2_z1_01, dzdzed2_z1_23;
				Cmp_DDZ_2Z(V4_0+48,-4,dzdzed1_z0_01, dzdzed1_z0_23, dzdzed1_z1_01, dzdzed1_z1_23);
				Cmp_DDZ_2Z(V4_0+50,-4,dzdzed2_z0_01, dzdzed2_z0_23, dzdzed2_z1_01, dzdzed2_z1_23);

				__m128 V5_01 = _mm_sub_ps(V4_0[20], dzdzed1_z0_01);
                                __m128 V5_23 = _mm_sub_ps(V4_0[21], dzdzed1_z0_23);
	
				__m128 V4_01 = _mm_add_ps(V4_0[22], dzdzed2_z0_01);
                                __m128 V4_23 = _mm_add_ps(V4_0[23], dzdzed2_z0_23);

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				TTIDenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelAnis1,DenAng1,spgzyx0,V4_01,V4_23,V5_01,V5_23
						);

				V5_01 = _mm_sub_ps(V4_0[16], dzdzed1_z1_01);
                                V5_23 = _mm_sub_ps(V4_0[17], dzdzed1_z1_23);

                                V4_01 = _mm_add_ps(V4_0[18], dzdzed2_z1_01);
                                V4_23 = _mm_add_ps(V4_0[19], dzdzed2_z1_23);

                                __m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
                                TTIDenQ_Cmp_Wave_Equation(
                                                stage,pq1+stride_z_m128,rs1+stride_z_m128,Apq1+stride_z_m128,VelAnis1+stride_z_em,DenAng1+stride_z_em,spgzyx1,V4_01,V4_23,V5_01,V5_23
                                                );

				/*
				__m128 dzdzed1_01, dzdzed1_23;
				Cmp_DDZ(V4_0+48,-4,dzdzed1_01,dzdzed1_23);
				__m128 V5_01 = _mm_sub_ps(V4_0[20], dzdzed1_01);
				__m128 V5_23 = _mm_sub_ps(V4_0[21], dzdzed1_23);

				__m128 dzdzed2_01, dzdzed2_23;
				Cmp_DDZ(V4_0+50,-4,dzdzed2_01,dzdzed2_23);
				__m128 V4_01 = _mm_add_ps(V4_0[22], dzdzed2_01);
				__m128 V4_23 = _mm_add_ps(V4_0[23], dzdzed2_23);

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				TTIDenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelAnis1,DenAng1,spgzyx0,V4_01,V4_23,V5_01,V5_23
						);

				Cmp_DDZ(V4_0+44,-4,dzdzed1_01,dzdzed1_23);
				V5_01 = _mm_sub_ps(V4_0[16], dzdzed1_01);
				V5_23 = _mm_sub_ps(V4_0[17], dzdzed1_23);

				Cmp_DDZ(V4_0+46,-4,dzdzed2_01,dzdzed2_23);
				V4_01 = _mm_add_ps(V4_0[18], dzdzed2_01);
				V4_23 = _mm_add_ps(V4_0[19], dzdzed2_23);

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				TTIDenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m128,rs1+stride_z_m128,Apq1+stride_z_m128,VelAnis1+stride_z_em,DenAng1+stride_z_em,spgzyx1,V4_01,V4_23,V5_01,V5_23
						);
				*/

				V4_0 += 68;

				pq1 += 2;
				rs1 += 2;
				Apq1 += 2;
				VelAnis1 += 4;
				DenAng1 += 4;
				stripe0 += 8;	
#ifdef TMJ_TIMING
				thrcnt2+=2;
#endif
			}
			pq1 += stride_y_m128 - iXN*2;
			rs1 += stride_y_m128 - iXN*2;
			Apq1 += stride_y_m128 - iXN*2;
			VelAnis1 += stride_y_em - iXN*4;
			DenAng1 += stride_y_em - iXN*4;
		}

		for (int iY = 4;  iY > 0;  --iY)
		{
			pq0 += 4;
			DenAng0 += 8;	
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				_mm_prefetch((char*)(pq0+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+5*stride_z_m128), _MM_HINT_T0);

				_mm_prefetch((char*)(pq0+stride_z_m128+5*stride_y_m128), _MM_HINT_T0);
				_mm_prefetch((char*)(pq0+6*stride_z_m128), _MM_HINT_T0);

				__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
				if (dz_edge_diff > 9)
				{
					Cmp_DDZ(pq0,stride_z_m128,ddz_01,ddz_23);
				}
				else
				{
					Cmp_DDZ_EE(pq0,stride_z_m128,ddz_01,ddz_23);
				}

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DYED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
						Boy, cDip, sDip, cAzm, sAzm,
						dyed_buf_0[0], dyed_buf_0[1], dyed_buf_0[2], dyed_buf_0[3]
						);

				Cmp_DDX(pq0+stride_z_m128,ddx_01,ddx_23);
				Cmp_DDY(pq0+stride_z_m128,stride_y_m128,ddy_01,ddy_23);
				if (dz_edge_diff > 10)
				{
					Cmp_DDZ(pq0+stride_z_m128,stride_z_m128,ddz_01,ddz_23);
				}
				else
				{
					Cmp_DDZ_EE(pq0+stride_z_m128,stride_z_m128,ddz_01,ddz_23);
				}

				TTIDenQ_Get_EM1(DenAng0+stride_z_em, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DYED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
						Boy, cDip, sDip, cAzm, sAzm,
                                                dyed_buf_0[4], dyed_buf_0[5], dyed_buf_0[6], dyed_buf_0[7]
						);

				dyed_buf_0 += 8;
				pq0 += 2;
				DenAng0 += 4;

				// dydyed1, dydyed2
				__m128 dydyed1_01, dydyed1_23;
				Cmp_DDY(dyed_buf_1,8*iXN,dydyed1_01,dydyed1_23);
				V4_1[4] = _mm_add_ps(V4_1[4], dydyed1_01);
				V4_1[5] = _mm_add_ps(V4_1[5], dydyed1_23);

				__m128 dydyed2_01, dydyed2_23;
				Cmp_DDY(dyed_buf_1+2,8*iXN,dydyed2_01,dydyed2_23);
				V4_1[6] = _mm_add_ps(V4_1[6], dydyed2_01);
				V4_1[7] = _mm_add_ps(V4_1[7], dydyed2_23);

				Cmp_DDY(dyed_buf_1+4,8*iXN,dydyed1_01,dydyed1_23);
				V4_1[0] = _mm_add_ps(V4_1[0], dydyed1_01);
				V4_1[1] = _mm_add_ps(V4_1[1], dydyed1_23);

				Cmp_DDY(dyed_buf_1+6,8*iXN,dydyed2_01,dydyed2_23);
				V4_1[2] = _mm_add_ps(V4_1[2], dydyed2_01);
				V4_1[3] = _mm_add_ps(V4_1[3], dydyed2_23);
				dyed_buf_1 += 8;
				V4_1 += 68;
			}
	
			pq0 += stride_y_m128 - (iXN+2)*2;
			DenAng0 += stride_y_em - (iXN+2)*4;
		}
	}
	else if (dz_edge_diff == 5)
	{
	 	__m128* pq0 = ((stage == 0 || stage == 1) ? pq : Apq) + (unsigned long)(iZ+4) * (unsigned long)stride_z_m128 + (unsigned long)(iY0-5) * (unsigned long)stride_y_m128 + (unsigned long)(iX0-4);
		int* DenAng0 = DenAng + (unsigned long)(iZ+4) * (unsigned long)stride_z_em + (unsigned long)(iY0-5) * (unsigned long)stride_y_em + (unsigned long)((iX0-4)*2);

		__m128* dyed_buf_0 = dyed_buf;
		__m128* dyed_buf_1 = dyed_buf + 32 * iXN;

		__m128* V4_0 = V4;
		__m128* V4_1 = V4;

		for (int iY = iYN_halo;  iY > (iYN_halo-5);  --iY)
		{
			pq0 += 4;
			DenAng0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
				Cmp_DDZ_EE(pq0,stride_z_m128,ddz_01,ddz_23);

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DYED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
						Boy, cDip, sDip, cAzm, sAzm,
						dyed_buf_0[0], dyed_buf_0[1], dyed_buf_0[2], dyed_buf_0[3]
						);
				dyed_buf_0 += 8;

				pq0 += 2;
				DenAng0 += 4;
			}

			pq0 += stride_y_m128 - (iXN+2)*2;
			DenAng0 += stride_y_em - (iXN+2)*4;
		}

		for (int iY = (iYN_halo-5);  iY > 4;  --iY)
		{
			__m128* stripe0 = stripe;

			for (int iX = 2;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
				Cmp_DDZ_EE(pq0,stride_z_m128,ddz_01,ddz_23);

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DXED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 0], stripe0[ 1], stripe0[ 2], stripe0[ 3]
						);
				stripe0 += 8;

				pq0 += 2;
				DenAng0 += 4;
			}

			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
				Cmp_DDZ_EE(pq0,stride_z_m128,ddz_01,ddz_23);

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DXED_DYED_DZED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 0], stripe0[ 1], stripe0[ 2], stripe0[ 3],
						dyed_buf_0[0], dyed_buf_0[1], dyed_buf_0[2], dyed_buf_0[3],
						V4_0[28], V4_0[29], V4_0[30], V4_0[31]
						);
				stripe0 += 8;
				dyed_buf_0 += 8;
				V4_0 += 68;

				pq0 += 2;
				DenAng0 += 4;
			}
			V4_0 -= 68 * iXN;

			for (int iX = 1;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
				Cmp_DDZ_EE(pq0,stride_z_m128,ddz_01,ddz_23);

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DXED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
						Boy, cDip, sDip, cAzm, sAzm,
						stripe0[ 0], stripe0[ 1], stripe0[ 2], stripe0[ 3]
						);
				stripe0 += 8;

				pq0 += 2;
				DenAng0 += 4;
			}

			pq0 += stride_y_m128 - iXN_halo*2;
			DenAng0 += stride_y_em - iXN_halo*4;

			stripe0 = stripe + 16;
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
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m128), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}

				if (iY <= (iYN_halo-9))
				{
					// dydyed1, dydyed2
					__m128 dydyed1_01, dydyed1_23;
					Cmp_DDY(dyed_buf_1,8*iXN,dydyed1_01,dydyed1_23);
					V4_1[4] = _mm_add_ps(V4_1[4], dydyed1_01);
					V4_1[5] = _mm_add_ps(V4_1[5], dydyed1_23);

					__m128 dydyed2_01, dydyed2_23;
					Cmp_DDY(dyed_buf_1+2,8*iXN,dydyed2_01,dydyed2_23);
					V4_1[6] = _mm_add_ps(V4_1[6], dydyed2_01);
					V4_1[7] = _mm_add_ps(V4_1[7], dydyed2_23);
					dyed_buf_1 += 8;
					V4_1 += 68;
				}

				// dxdxed1, dxdxed2
				__m128 dxdxed1_01, dxdxed1_23;
				Cmp_DDX_1x_stride8(stripe0,dxdxed1_01,dxdxed1_23);
				V4_0[4] = dxdxed1_01;
				V4_0[5] = dxdxed1_23;

				__m128 dxdxed2_01, dxdxed2_23;
				Cmp_DDX_1x_stride8(stripe0+2,dxdxed2_01,dxdxed2_23);
				V4_0[6] = dxdxed2_01;
				V4_0[7] = dxdxed2_23;

				// dzdzed1, dzdzed2
				__m128 dzdzed1_01, dzdzed1_23;
				Cmp_DDZ(V4_0+48,-4,dzdzed1_01,dzdzed1_23);
				__m128 V5_01 = _mm_sub_ps(V4_0[20], dzdzed1_01);
				__m128 V5_23 = _mm_sub_ps(V4_0[21], dzdzed1_23);

				__m128 dzdzed2_01, dzdzed2_23;
				Cmp_DDZ(V4_0+50,-4,dzdzed2_01,dzdzed2_23);
				__m128 V4_01 = _mm_add_ps(V4_0[22], dzdzed2_01);
				__m128 V4_23 = _mm_add_ps(V4_0[23], dzdzed2_23);

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				TTIDenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelAnis1,DenAng1,spgzyx0,V4_01,V4_23,V5_01,V5_23
						);

				Cmp_DDZ_EE(V4_0+44,-4,dzdzed1_01,dzdzed1_23);
				V5_01 = _mm_sub_ps(V4_0[16], dzdzed1_01);
				V5_23 = _mm_sub_ps(V4_0[17], dzdzed1_23);

				Cmp_DDZ_EE(V4_0+46,-4,dzdzed2_01,dzdzed2_23);
				V4_01 = _mm_add_ps(V4_0[18], dzdzed2_01);
				V4_23 = _mm_add_ps(V4_0[19], dzdzed2_23);

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				TTIDenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m128,rs1+stride_z_m128,Apq1+stride_z_m128,VelAnis1+stride_z_em,DenAng1+stride_z_em,spgzyx1,V4_01,V4_23,V5_01,V5_23
						);

				V4_0 += 68;

				pq1 += 2;
				rs1 += 2;
				Apq1 += 2;
				VelAnis1 += 4;
				DenAng1 += 4;
				stripe0 += 8;	
#ifdef TMJ_TIMING
				thrcnt2+=2;
#endif
			}
			pq1 += stride_y_m128 - iXN*2;
			rs1 += stride_y_m128 - iXN*2;
			Apq1 += stride_y_m128 - iXN*2;
			VelAnis1 += stride_y_em - iXN*4;
			DenAng1 += stride_y_em - iXN*4;
		}

		for (int iY = 4;  iY > 0;  --iY)
		{
			pq0 += 4;
			DenAng0 += 8;
			for (int iX = iXN;  iX > 0;  --iX)
			{
				_mm_prefetch((char*)(DenAng0), _MM_HINT_NTA);
				_mm_prefetch((char*)(DenAng0+stride_z_em), _MM_HINT_NTA);

				__m128 ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23;
				Cmp_DDX(pq0,ddx_01,ddx_23);
				Cmp_DDY(pq0,stride_y_m128,ddy_01,ddy_23);
				Cmp_DDZ_EE(pq0,stride_z_m128,ddz_01,ddz_23);

				__m128 Boy, cDip, sDip, cAzm, sAzm;
				TTIDenQ_Get_EM1(DenAng0, Boy, cDip, sDip, cAzm, sAzm);

				TTIDenQ_Cmp_DYED(
						ddx_01, ddx_23, ddy_01, ddy_23, ddz_01, ddz_23, 
						Boy, cDip, sDip, cAzm, sAzm,
						dyed_buf_0[0], dyed_buf_0[1], dyed_buf_0[2], dyed_buf_0[3]
						);
				dyed_buf_0 += 8;

				pq0 += 2;
				DenAng0 += 4;

				// dydyed1, dydyed2
				__m128 dydyed1_01, dydyed1_23;
				Cmp_DDY(dyed_buf_1,8*iXN,dydyed1_01,dydyed1_23);
				V4_1[4] = _mm_add_ps(V4_1[4], dydyed1_01);
				V4_1[5] = _mm_add_ps(V4_1[5], dydyed1_23);

				__m128 dydyed2_01, dydyed2_23;
				Cmp_DDY(dyed_buf_1+2,8*iXN,dydyed2_01,dydyed2_23);
				V4_1[6] = _mm_add_ps(V4_1[6], dydyed2_01);
				V4_1[7] = _mm_add_ps(V4_1[7], dydyed2_23);
				dyed_buf_1 += 8;
				V4_1 += 68;
			}

			pq0 += stride_y_m128 - (iXN+2)*2;
			DenAng0 += stride_y_em - (iXN+2)*4;
		}
	}
	else // no more tiles to process, finish rest of outputs from queue
	{
		__m128* V4_0 = V4;
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
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}
				else if (stage == 2)
				{
					_mm_prefetch((char*)pq1, _MM_HINT_NTA);
					_mm_prefetch((char*)(pq1+stride_z_m128), _MM_HINT_NTA);
					_mm_prefetch((char*)rs1, _MM_HINT_T0);
					_mm_prefetch((char*)(rs1+stride_z_m128), _MM_HINT_T0);
				}

				// dzdzed1, dzdzed2
				// use shorter two point stencil
				__m128 dzdzed1_01, dzdzed1_23;
				Cmp_DDZ_EE(V4_0+48,-4,dzdzed1_01,dzdzed1_23);
				__m128 V5_01 = _mm_sub_ps(V4_0[20], dzdzed1_01);
				__m128 V5_23 = _mm_sub_ps(V4_0[21], dzdzed1_23);

				__m128 dzdzed2_01, dzdzed2_23;
				Cmp_DDZ_EE(V4_0+50,-4,dzdzed2_01,dzdzed2_23);
				__m128 V4_01 = _mm_add_ps(V4_0[22], dzdzed2_01);
				__m128 V4_23 = _mm_add_ps(V4_0[23], dzdzed2_23);

				__m128 spgzyx0 = _mm_mul_ps(_mm_spgzy0, _mm_spgx[iX]);
				TTIDenQ_Cmp_Wave_Equation(
						stage,pq1,rs1,Apq1,VelAnis1,DenAng1,spgzyx0,V4_01,V4_23,V5_01,V5_23
						);

				Cmp_DDZ_EE(V4_0+44,-4,dzdzed1_01,dzdzed1_23);
				V5_01 = _mm_sub_ps(V4_0[16], dzdzed1_01);
				V5_23 = _mm_sub_ps(V4_0[17], dzdzed1_23);

				Cmp_DDZ_EE(V4_0+46,-4,dzdzed2_01,dzdzed2_23);
				V4_01 = _mm_add_ps(V4_0[18], dzdzed2_01);
				V4_23 = _mm_add_ps(V4_0[19], dzdzed2_23);

				__m128 spgzyx1 = _mm_mul_ps(_mm_spgzy1, _mm_spgx[iX]);
				TTIDenQ_Cmp_Wave_Equation(
						stage,pq1+stride_z_m128,rs1+stride_z_m128,Apq1+stride_z_m128,VelAnis1+stride_z_em,DenAng1+stride_z_em,spgzyx1,V4_01,V4_23,V5_01,V5_23
						);

				V4_0 += 68;

				pq1 += 2;
				rs1 += 2;
				Apq1 += 2;
				VelAnis1 += 4;
				DenAng1 += 4;
#ifdef TMJ_TIMING
				thrcnt2 += 2;
#endif
			}
			pq1 += stride_y_m128 - iXN*2;
			rs1 += stride_y_m128 - iXN*2;
			Apq1 += stride_y_m128 - iXN*2;
			VelAnis1 += stride_y_em - iXN*4;
			DenAng1 += stride_y_em - iXN*4;
		}
	}
}

} // end of anonymous namespace

void ISODenQ_TimeStep(
	int logLevel,
	int stage,
	__m128* pq,
	__m128* rs,
	__m128* Apq,
	int* VelDen,
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
	/*
	__m128* p = new __m128[64];
	__m128* o = new __m128[64];
	for (int i = 0;  i < 64;  ++i)
	{
		p[i] = o[i] = _mm_set1_ps(0.0f);
	}
	p[32] = _mm_setr_ps(1.0f, 0.0f, 0.0f, 0.0f);

	for (int i = 8;  i < 56;  ++i)
	{
		ISO_Cmp_DDX_1x_stride2(p+i,o[i]);
	}

	float* fo = (float*)o;
	for (int i = 0;  i < 256;  ++i)
	{
		float val = fo[i];
		if (val != 0.0f)
		{
			printf("fo[%d] = %e\n",i,val);
		}
	}
	
	exit(0);
	*/

	// adjust x0 so that it is always a multiple of MIN_BSX
	int xx0 = x0 >> 2;
	int xx1 = x1 >> 2;

	int dimxh = dimx + 2*xh;
	int dimyh = dimy + 2*yh;
	int dimzh = dimz + 2*zh;

        int stride_y_m128 = dimxh / 4;
        int stride_z_m128 = stride_y_m128 * dimyh;

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
	ISODenQ_Comp_Buf_Size(bsX, bsY, 128, stripe_m128, dyed_buf_m128, spg_buf_m128, QMaxShift, QLen, QSize, BufSize_m128, NetBufSize_m128);

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
	__m128* pqnh = pq + zh * stride_z_m128 + yh * stride_y_m128 + (xh >> 2);
	__m128* rsnh = rs + zh * stride_z_m128 + yh * stride_y_m128 + (xh >> 2);
	__m128* Apqnh = Apq + zh * stride_z_m128 + yh * stride_y_m128 + (xh >> 2);
	
	int* VelDen_nh = VelDen + zh * stride_z_em + yh * stride_y_em + xh;

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
		__m128* V4 = Q + QShift;

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
			ISODenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);
			ISODenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);
			ISODenQ_Process_Patch_Leadin(
				logLevel,
				((stage == 0 || stage == 1) ? pqnh : Apqnh),VelDen_nh,
				stripe,dyed_buf,V4,
				iX0*2,iY0,iZ,iXN_halo,iXN,iYN_halo,iYN,
				stride_y_m128,stride_z_m128,stride_y_em,stride_z_em);
		}
		ISODenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);	

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
			ISODenQ_Process_Patch(
				logLevel,
				stage,pqnh,rsnh,Apqnh,VelDen_nh,_mm_spgx,_mm_spgy,_mm_spgz0,_mm_spgz1,
				stripe,dyed_buf,V4,
				iX0*2,iY0,iZ,dimz,iXN_halo,iXN,iYN_halo,iYN,stride_y_m128,stride_z_m128,stride_y_em,stride_z_em
#ifdef TMJ_TIMING
				,thrcnt2
#endif
				);
			ISODenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);
			ISODenQ_Shift_Queues(QMaxShift,QLen,Q,V4,QShift);
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

void VTIDenQ_TimeStep(
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
			VTIDenQ_Process_Patch_Leadin(
				logLevel,
				((stage == 0 || stage == 1) ? pqnh : Apqnh),DenAng_nh,VelAnis_nh,
				stripe,dyed_buf,V4,
				iX0*2,iY0,iZ,iXN_halo,iXN,iYN_halo,iYN,
				stride_y_m128,stride_z_m128,stride_y_em,stride_z_em);
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
			VTIDenQ_Process_Patch(
				logLevel,
				stage,pqnh,rsnh,Apqnh,DenAng_nh,VelAnis_nh,_mm_spgx,_mm_spgy,_mm_spgz0,_mm_spgz1,
				stripe,dyed_buf,V4,
				iX0*2,iY0,iZ,dimz,iXN_halo,iXN,iYN_halo,iYN,stride_y_m128,stride_z_m128,stride_y_em,stride_z_em
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

void TTIDenQ_TimeStep(
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
		printf("TTIDenQ_TimeStep\n");
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
			TTIDenQ_Process_Patch_Leadin(
				logLevel,
				((stage == 0 || stage == 1) ? pqnh : Apqnh),DenAng_nh,VelAnis_nh,
				stripe,dyed_buf,V4,
				iX0*2,iY0,iZ,iXN_halo,iXN,iYN_halo,iYN,
				stride_y_m128,stride_z_m128,stride_y_em,stride_z_em);
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
			TTIDenQ_Process_Patch(
				logLevel,
				stage,pqnh,rsnh,Apqnh,DenAng_nh,VelAnis_nh,_mm_spgx,_mm_spgy,_mm_spgz0,_mm_spgz1,
				stripe,dyed_buf,V4,
				iX0*2,iY0,iZ,dimz,iXN_halo,iXN,iYN_halo,iYN,stride_y_m128,stride_z_m128,stride_y_em,stride_z_em
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

void Call_Cmp_11pt_DDX(
	__m128* pq,
	__m128& ddx_01,
	__m128& ddx_23
	)
{
	Cmp_11pt_DDX(pq,ddx_01,ddx_23);
}

void Call_Cmp_11pt_DDX_DDX2(
	__m128* pq,
	__m128& ddx_01,
	__m128& ddx_23,
	__m128& ddx2_01,
	__m128& ddx2_23
	)
{
	Cmp_11pt_DDX_DDX2(pq,ddx_01,ddx_23,ddx2_01,ddx2_23);
}

