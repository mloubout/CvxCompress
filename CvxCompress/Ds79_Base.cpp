/*!
 * Don't edit this code, it was automatically generated.
 * Base functions for wavelet transforms of length 2 to 256.
 */
#include <immintrin.h>

/*
 * Define coefficients for Antonini 7-9 tap filter.
 */
#define al0  8.526986790094000e-001f
#define al1  3.774028556126500e-001f
#define al2 -1.106244044184200e-001f
#define al3 -2.384946501938001e-002f
#define al4  3.782845550699501e-002f

#define ah0  7.884856164056601e-001f
#define ah1 -4.180922732222101e-001f
#define ah2 -4.068941760955800e-002f
#define ah3  6.453888262893799e-002f

#define _mm_al0 _mm256_set1_ps(al0)
#define _mm_al1 _mm256_set1_ps(al1)
#define _mm_al2 _mm256_set1_ps(al2)
#define _mm_al3 _mm256_set1_ps(al3)
#define _mm_al4 _mm256_set1_ps(al4)

#define _mm_ah0 _mm256_set1_ps(ah0)
#define _mm_ah1 _mm256_set1_ps(ah1)
#define _mm_ah2 _mm256_set1_ps(ah2)
#define _mm_ah3 _mm256_set1_ps(ah3)

#ifdef __AVX2__

static inline void _Ds79_AVX_2(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v0 = data[0*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v0,v0));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v1,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v0,v0),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v0,acc1);
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v0));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v0,v0),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v1,acc1);
	data[1*stride] = acc1;
}

static inline void _Ds79_AVX_4(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v2 = data[2*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v2));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v3,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v2,v2),acc1);
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v1),acc1);
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_al0, v0,acc1);
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v2));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v0,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v1,acc1);
	data[2*stride] = acc1;

	// lower band :: ix=1
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v0));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v1,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v0,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v2,acc1);
	data[1*stride] = acc1;

	// upper band :: ix=1
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v0));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v2,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v3,acc1);
	data[3*stride] = acc1;
}

static inline void _Ds79_AVX_8(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v4 = data[4*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v4));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v3,v3),acc1);
	__m256 v2 = data[2*stride];
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v2,v2),acc1);
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v1),acc1);
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_al0, v0,acc1);
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v4));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v0,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v1,acc1);
	data[4*stride] = acc1;

	// lower band :: ix=1
	__m256 v6 = data[6*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v6));
	__m256 v5 = data[5*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v1,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v0,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v2,acc1);
	data[1*stride] = acc1;

	// upper band :: ix=1
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v6));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v2,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v3,acc1);
	data[5*stride] = acc1;

	// lower band :: ix=2
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v0,v6));
	__m256 v7 = data[7*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v1,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v2,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v3,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v4,acc1);
	data[2*stride] = acc1;

	// upper band :: ix=2
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v6));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v3,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v4,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v5,acc1);
	data[6*stride] = acc1;

	// lower band :: ix=3
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v4));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v3,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v4,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v5,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v6,acc1);
	data[3*stride] = acc1;

	// upper band :: ix=3
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v4,v4));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v5,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v6,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v7,acc1);
	data[7*stride] = acc1;
}

static inline void _Ds79_AVX_16(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v4 = data[4*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v4));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v3,v3),acc1);
	__m256 v2 = data[2*stride];
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v2,v2),acc1);
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v1),acc1);
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_al0, v0,acc1);
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v4));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v0,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v1,acc1);
	__m256 v8 = data[8*stride];
	data[8*stride] = acc1;

	// lower band :: ix=1
	__m256 v6 = data[6*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v6));
	__m256 v5 = data[5*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v1,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v0,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v2,acc1);
	data[1*stride] = acc1;

	// upper band :: ix=1
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v6));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v2,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v3,acc1);
	__m256 v9 = data[9*stride];
	data[9*stride] = acc1;

	// lower band :: ix=2
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v0,v8));
	__m256 v7 = data[7*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v1,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v2,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v3,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v4,acc1);
	data[2*stride] = acc1;

	// upper band :: ix=2
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v8));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v3,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v4,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v5,acc1);
	__m256 v10 = data[10*stride];
	data[10*stride] = acc1;

	// lower band :: ix=3
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v10));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v3,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v4,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v5,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v6,acc1);
	data[3*stride] = acc1;

	// upper band :: ix=3
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v4,v10));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v5,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v6,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v7,acc1);
	__m256 v11 = data[11*stride];
	data[11*stride] = acc1;

	// lower band :: ix=4
	__m256 v12 = data[12*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v12));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v5,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v6,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v7,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v8,acc1);
	data[4*stride] = acc1;

	// upper band :: ix=4
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v6,v12));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v7,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v8,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v9,acc1);
	data[12*stride] = acc1;

	// lower band :: ix=5
	__m256 v14 = data[14*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v6,v14));
	__m256 v13 = data[13*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v7,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v8,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v9,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v10,acc1);
	data[5*stride] = acc1;

	// upper band :: ix=5
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v8,v14));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v9,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v10,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v11,acc1);
	data[13*stride] = acc1;

	// lower band :: ix=6
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v8,v14));
	__m256 v15 = data[15*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v9,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v10,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v11,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v12,acc1);
	data[6*stride] = acc1;

	// upper band :: ix=6
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v10,v14));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v11,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v12,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v13,acc1);
	data[14*stride] = acc1;

	// lower band :: ix=7
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v10,v12));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v11,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v12,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v13,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v14,acc1);
	data[7*stride] = acc1;

	// upper band :: ix=7
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v12,v12));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v13,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v14,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v15,acc1);
	data[15*stride] = acc1;
}

static inline void _Ds79_AVX_32(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v4 = data[4*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v4));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v3,v3),acc1);
	__m256 v2 = data[2*stride];
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v2,v2),acc1);
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v1),acc1);
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_al0, v0,acc1);
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v4));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v0,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v1,acc1);
	__m256 v16 = data[16*stride];
	data[16*stride] = acc1;

	// lower band :: ix=1
	__m256 v6 = data[6*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v6));
	__m256 v5 = data[5*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v1,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v0,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v2,acc1);
	data[1*stride] = acc1;

	// upper band :: ix=1
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v6));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v2,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v3,acc1);
	__m256 v17 = data[17*stride];
	data[17*stride] = acc1;

	// lower band :: ix=2
	__m256 v8 = data[8*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v0,v8));
	__m256 v7 = data[7*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v1,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v2,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v3,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v4,acc1);
	data[2*stride] = acc1;

	// upper band :: ix=2
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v8));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v3,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v4,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v5,acc1);
	__m256 v18 = data[18*stride];
	data[18*stride] = acc1;

	// lower band :: ix=3
	__m256 v10 = data[10*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v10));
	__m256 v9 = data[9*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v3,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v4,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v5,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v6,acc1);
	data[3*stride] = acc1;

	// upper band :: ix=3
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v4,v10));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v5,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v6,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v7,acc1);
	__m256 v19 = data[19*stride];
	data[19*stride] = acc1;

	// lower band :: ix=4
	__m256 v12 = data[12*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v12));
	__m256 v11 = data[11*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v5,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v6,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v7,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v8,acc1);
	data[4*stride] = acc1;

	// upper band :: ix=4
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v6,v12));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v7,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v8,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v9,acc1);
	__m256 v20 = data[20*stride];
	data[20*stride] = acc1;

	// lower band :: ix=5
	__m256 v14 = data[14*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v6,v14));
	__m256 v13 = data[13*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v7,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v8,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v9,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v10,acc1);
	data[5*stride] = acc1;

	// upper band :: ix=5
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v8,v14));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v9,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v10,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v11,acc1);
	__m256 v21 = data[21*stride];
	data[21*stride] = acc1;

	// lower band :: ix=6
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v8,v16));
	__m256 v15 = data[15*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v9,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v10,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v11,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v12,acc1);
	data[6*stride] = acc1;

	// upper band :: ix=6
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v10,v16));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v11,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v12,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v13,acc1);
	__m256 v22 = data[22*stride];
	data[22*stride] = acc1;

	// lower band :: ix=7
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v10,v18));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v11,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v12,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v13,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v14,acc1);
	data[7*stride] = acc1;

	// upper band :: ix=7
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v12,v18));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v13,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v14,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v15,acc1);
	__m256 v23 = data[23*stride];
	data[23*stride] = acc1;

	// lower band :: ix=8
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v12,v20));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v13,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v14,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v15,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v16,acc1);
	data[8*stride] = acc1;

	// upper band :: ix=8
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v14,v20));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v15,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v16,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v17,acc1);
	__m256 v24 = data[24*stride];
	data[24*stride] = acc1;

	// lower band :: ix=9
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v14,v22));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v15,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v16,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v17,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v18,acc1);
	data[9*stride] = acc1;

	// upper band :: ix=9
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v16,v22));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v17,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v18,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v19,acc1);
	__m256 v25 = data[25*stride];
	data[25*stride] = acc1;

	// lower band :: ix=10
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v16,v24));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v17,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v18,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v19,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v20,acc1);
	data[10*stride] = acc1;

	// upper band :: ix=10
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v18,v24));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v19,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v20,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v21,acc1);
	__m256 v26 = data[26*stride];
	data[26*stride] = acc1;

	// lower band :: ix=11
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v18,v26));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v19,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v20,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v21,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v22,acc1);
	data[11*stride] = acc1;

	// upper band :: ix=11
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v20,v26));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v21,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v22,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v23,acc1);
	__m256 v27 = data[27*stride];
	data[27*stride] = acc1;

	// lower band :: ix=12
	__m256 v28 = data[28*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v20,v28));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v21,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v22,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v23,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v24,acc1);
	data[12*stride] = acc1;

	// upper band :: ix=12
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v22,v28));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v23,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v24,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v25,acc1);
	data[28*stride] = acc1;

	// lower band :: ix=13
	__m256 v30 = data[30*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v22,v30));
	__m256 v29 = data[29*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v23,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v24,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v25,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v26,acc1);
	data[13*stride] = acc1;

	// upper band :: ix=13
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v24,v30));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v25,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v26,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v27,acc1);
	data[29*stride] = acc1;

	// lower band :: ix=14
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v24,v30));
	__m256 v31 = data[31*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v25,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v26,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v27,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v28,acc1);
	data[14*stride] = acc1;

	// upper band :: ix=14
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v26,v30));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v27,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v28,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v29,acc1);
	data[30*stride] = acc1;

	// lower band :: ix=15
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v26,v28));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v27,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v28,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v29,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v30,acc1);
	data[15*stride] = acc1;

	// upper band :: ix=15
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v28,v28));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v29,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v30,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v31,acc1);
	data[31*stride] = acc1;
}

static void _Ds79_AVX_64(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v4 = data[4*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v4));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v3,v3),acc1);
	__m256 v2 = data[2*stride];
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v2,v2),acc1);
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v1),acc1);
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_al0, v0,acc1);
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v4));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v0,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v1,acc1);
	__m256 v32 = data[32*stride];
	data[32*stride] = acc1;

	// lower band :: ix=1
	__m256 v6 = data[6*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v6));
	__m256 v5 = data[5*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v1,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v0,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v2,acc1);
	data[1*stride] = acc1;

	// upper band :: ix=1
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v6));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v2,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v3,acc1);
	__m256 v33 = data[33*stride];
	data[33*stride] = acc1;

	// lower band :: ix=2
	__m256 v8 = data[8*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v0,v8));
	__m256 v7 = data[7*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v1,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v2,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v3,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v4,acc1);
	data[2*stride] = acc1;

	// upper band :: ix=2
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v8));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v3,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v4,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v5,acc1);
	__m256 v34 = data[34*stride];
	data[34*stride] = acc1;

	// lower band :: ix=3
	__m256 v10 = data[10*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v10));
	__m256 v9 = data[9*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v3,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v4,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v5,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v6,acc1);
	data[3*stride] = acc1;

	// upper band :: ix=3
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v4,v10));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v5,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v6,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v7,acc1);
	__m256 v35 = data[35*stride];
	data[35*stride] = acc1;

	// lower band :: ix=4
	__m256 v12 = data[12*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v12));
	__m256 v11 = data[11*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v5,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v6,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v7,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v8,acc1);
	data[4*stride] = acc1;

	// upper band :: ix=4
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v6,v12));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v7,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v8,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v9,acc1);
	__m256 v36 = data[36*stride];
	data[36*stride] = acc1;

	// lower band :: ix=5
	__m256 v14 = data[14*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v6,v14));
	__m256 v13 = data[13*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v7,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v8,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v9,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v10,acc1);
	data[5*stride] = acc1;

	// upper band :: ix=5
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v8,v14));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v9,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v10,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v11,acc1);
	__m256 v37 = data[37*stride];
	data[37*stride] = acc1;

	// lower band :: ix=6
	__m256 v16 = data[16*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v8,v16));
	__m256 v15 = data[15*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v9,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v10,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v11,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v12,acc1);
	data[6*stride] = acc1;

	// upper band :: ix=6
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v10,v16));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v11,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v12,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v13,acc1);
	__m256 v38 = data[38*stride];
	data[38*stride] = acc1;

	// lower band :: ix=7
	__m256 v18 = data[18*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v10,v18));
	__m256 v17 = data[17*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v11,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v12,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v13,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v14,acc1);
	data[7*stride] = acc1;

	// upper band :: ix=7
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v12,v18));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v13,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v14,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v15,acc1);
	__m256 v39 = data[39*stride];
	data[39*stride] = acc1;

	// lower band :: ix=8
	__m256 v20 = data[20*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v12,v20));
	__m256 v19 = data[19*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v13,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v14,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v15,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v16,acc1);
	data[8*stride] = acc1;

	// upper band :: ix=8
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v14,v20));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v15,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v16,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v17,acc1);
	__m256 v40 = data[40*stride];
	data[40*stride] = acc1;

	// lower band :: ix=9
	__m256 v22 = data[22*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v14,v22));
	__m256 v21 = data[21*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v15,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v16,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v17,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v18,acc1);
	data[9*stride] = acc1;

	// upper band :: ix=9
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v16,v22));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v17,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v18,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v19,acc1);
	__m256 v41 = data[41*stride];
	data[41*stride] = acc1;

	// lower band :: ix=10
	__m256 v24 = data[24*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v16,v24));
	__m256 v23 = data[23*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v17,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v18,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v19,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v20,acc1);
	data[10*stride] = acc1;

	// upper band :: ix=10
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v18,v24));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v19,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v20,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v21,acc1);
	__m256 v42 = data[42*stride];
	data[42*stride] = acc1;

	// lower band :: ix=11
	__m256 v26 = data[26*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v18,v26));
	__m256 v25 = data[25*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v19,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v20,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v21,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v22,acc1);
	data[11*stride] = acc1;

	// upper band :: ix=11
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v20,v26));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v21,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v22,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v23,acc1);
	__m256 v43 = data[43*stride];
	data[43*stride] = acc1;

	// lower band :: ix=12
	__m256 v28 = data[28*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v20,v28));
	__m256 v27 = data[27*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v21,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v22,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v23,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v24,acc1);
	data[12*stride] = acc1;

	// upper band :: ix=12
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v22,v28));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v23,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v24,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v25,acc1);
	__m256 v44 = data[44*stride];
	data[44*stride] = acc1;

	// lower band :: ix=13
	__m256 v30 = data[30*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v22,v30));
	__m256 v29 = data[29*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v23,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v24,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v25,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v26,acc1);
	data[13*stride] = acc1;

	// upper band :: ix=13
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v24,v30));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v25,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v26,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v27,acc1);
	__m256 v45 = data[45*stride];
	data[45*stride] = acc1;

	// lower band :: ix=14
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v24,v32));
	__m256 v31 = data[31*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v25,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v26,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v27,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v28,acc1);
	data[14*stride] = acc1;

	// upper band :: ix=14
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v26,v32));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v27,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v28,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v29,acc1);
	__m256 v46 = data[46*stride];
	data[46*stride] = acc1;

	// lower band :: ix=15
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v26,v34));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v27,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v28,v32),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v29,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v30,acc1);
	data[15*stride] = acc1;

	// upper band :: ix=15
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v28,v34));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v29,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v30,v32),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v31,acc1);
	__m256 v47 = data[47*stride];
	data[47*stride] = acc1;

	// lower band :: ix=16
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v28,v36));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v29,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v30,v34),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v31,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v32,acc1);
	data[16*stride] = acc1;

	// upper band :: ix=16
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v30,v36));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v31,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v32,v34),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v33,acc1);
	__m256 v48 = data[48*stride];
	data[48*stride] = acc1;

	// lower band :: ix=17
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v30,v38));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v31,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v32,v36),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v33,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v34,acc1);
	data[17*stride] = acc1;

	// upper band :: ix=17
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v32,v38));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v33,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v34,v36),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v35,acc1);
	__m256 v49 = data[49*stride];
	data[49*stride] = acc1;

	// lower band :: ix=18
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v32,v40));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v33,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v34,v38),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v35,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v36,acc1);
	data[18*stride] = acc1;

	// upper band :: ix=18
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v34,v40));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v35,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v36,v38),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v37,acc1);
	__m256 v50 = data[50*stride];
	data[50*stride] = acc1;

	// lower band :: ix=19
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v34,v42));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v35,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v36,v40),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v37,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v38,acc1);
	data[19*stride] = acc1;

	// upper band :: ix=19
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v36,v42));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v37,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v38,v40),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v39,acc1);
	__m256 v51 = data[51*stride];
	data[51*stride] = acc1;

	// lower band :: ix=20
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v36,v44));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v37,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v38,v42),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v39,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v40,acc1);
	data[20*stride] = acc1;

	// upper band :: ix=20
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v38,v44));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v39,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v40,v42),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v41,acc1);
	__m256 v52 = data[52*stride];
	data[52*stride] = acc1;

	// lower band :: ix=21
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v38,v46));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v39,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v40,v44),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v41,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v42,acc1);
	data[21*stride] = acc1;

	// upper band :: ix=21
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v40,v46));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v41,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v42,v44),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v43,acc1);
	__m256 v53 = data[53*stride];
	data[53*stride] = acc1;

	// lower band :: ix=22
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v40,v48));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v41,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v42,v46),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v43,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v44,acc1);
	data[22*stride] = acc1;

	// upper band :: ix=22
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v42,v48));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v43,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v44,v46),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v45,acc1);
	__m256 v54 = data[54*stride];
	data[54*stride] = acc1;

	// lower band :: ix=23
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v42,v50));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v43,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v44,v48),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v45,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v46,acc1);
	data[23*stride] = acc1;

	// upper band :: ix=23
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v44,v50));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v45,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v46,v48),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v47,acc1);
	__m256 v55 = data[55*stride];
	data[55*stride] = acc1;

	// lower band :: ix=24
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v44,v52));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v45,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v46,v50),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v47,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v48,acc1);
	data[24*stride] = acc1;

	// upper band :: ix=24
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v46,v52));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v47,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v48,v50),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v49,acc1);
	__m256 v56 = data[56*stride];
	data[56*stride] = acc1;

	// lower band :: ix=25
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v46,v54));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v47,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v48,v52),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v49,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v50,acc1);
	data[25*stride] = acc1;

	// upper band :: ix=25
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v48,v54));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v49,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v50,v52),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v51,acc1);
	__m256 v57 = data[57*stride];
	data[57*stride] = acc1;

	// lower band :: ix=26
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v48,v56));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v49,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v50,v54),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v51,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v52,acc1);
	data[26*stride] = acc1;

	// upper band :: ix=26
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v50,v56));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v51,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v52,v54),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v53,acc1);
	__m256 v58 = data[58*stride];
	data[58*stride] = acc1;

	// lower band :: ix=27
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v50,v58));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v51,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v52,v56),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v53,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v54,acc1);
	data[27*stride] = acc1;

	// upper band :: ix=27
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v52,v58));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v53,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v54,v56),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v55,acc1);
	__m256 v59 = data[59*stride];
	data[59*stride] = acc1;

	// lower band :: ix=28
	__m256 v60 = data[60*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v52,v60));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v53,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v54,v58),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v55,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v56,acc1);
	data[28*stride] = acc1;

	// upper band :: ix=28
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v54,v60));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v55,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v56,v58),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v57,acc1);
	data[60*stride] = acc1;

	// lower band :: ix=29
	__m256 v62 = data[62*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v54,v62));
	__m256 v61 = data[61*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v55,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v56,v60),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v57,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v58,acc1);
	data[29*stride] = acc1;

	// upper band :: ix=29
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v56,v62));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v57,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v58,v60),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v59,acc1);
	data[61*stride] = acc1;

	// lower band :: ix=30
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v56,v62));
	__m256 v63 = data[63*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v57,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v58,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v59,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v60,acc1);
	data[30*stride] = acc1;

	// upper band :: ix=30
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v58,v62));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v59,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v60,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v61,acc1);
	data[62*stride] = acc1;

	// lower band :: ix=31
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v58,v60));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v59,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v60,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v61,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v62,acc1);
	data[31*stride] = acc1;

	// upper band :: ix=31
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v60,v60));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v61,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v62,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v63,acc1);
	data[63*stride] = acc1;
}

static void _Ds79_AVX_128(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v4 = data[4*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v4));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v3,v3),acc1);
	__m256 v2 = data[2*stride];
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v2,v2),acc1);
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v1),acc1);
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_al0, v0,acc1);
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v4));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v0,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v1,acc1);
	__m256 v64 = data[64*stride];
	data[64*stride] = acc1;

	// lower band :: ix=1
	__m256 v6 = data[6*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v6));
	__m256 v5 = data[5*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v1,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v0,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v2,acc1);
	data[1*stride] = acc1;

	// upper band :: ix=1
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v6));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v2,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v3,acc1);
	__m256 v65 = data[65*stride];
	data[65*stride] = acc1;

	// lower band :: ix=2
	__m256 v8 = data[8*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v0,v8));
	__m256 v7 = data[7*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v1,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v2,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v3,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v4,acc1);
	data[2*stride] = acc1;

	// upper band :: ix=2
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v8));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v3,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v4,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v5,acc1);
	__m256 v66 = data[66*stride];
	data[66*stride] = acc1;

	// lower band :: ix=3
	__m256 v10 = data[10*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v10));
	__m256 v9 = data[9*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v3,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v4,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v5,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v6,acc1);
	data[3*stride] = acc1;

	// upper band :: ix=3
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v4,v10));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v5,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v6,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v7,acc1);
	__m256 v67 = data[67*stride];
	data[67*stride] = acc1;

	// lower band :: ix=4
	__m256 v12 = data[12*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v12));
	__m256 v11 = data[11*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v5,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v6,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v7,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v8,acc1);
	data[4*stride] = acc1;

	// upper band :: ix=4
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v6,v12));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v7,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v8,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v9,acc1);
	__m256 v68 = data[68*stride];
	data[68*stride] = acc1;

	// lower band :: ix=5
	__m256 v14 = data[14*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v6,v14));
	__m256 v13 = data[13*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v7,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v8,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v9,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v10,acc1);
	data[5*stride] = acc1;

	// upper band :: ix=5
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v8,v14));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v9,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v10,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v11,acc1);
	__m256 v69 = data[69*stride];
	data[69*stride] = acc1;

	// lower band :: ix=6
	__m256 v16 = data[16*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v8,v16));
	__m256 v15 = data[15*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v9,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v10,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v11,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v12,acc1);
	data[6*stride] = acc1;

	// upper band :: ix=6
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v10,v16));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v11,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v12,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v13,acc1);
	__m256 v70 = data[70*stride];
	data[70*stride] = acc1;

	// lower band :: ix=7
	__m256 v18 = data[18*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v10,v18));
	__m256 v17 = data[17*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v11,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v12,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v13,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v14,acc1);
	data[7*stride] = acc1;

	// upper band :: ix=7
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v12,v18));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v13,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v14,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v15,acc1);
	__m256 v71 = data[71*stride];
	data[71*stride] = acc1;

	// lower band :: ix=8
	__m256 v20 = data[20*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v12,v20));
	__m256 v19 = data[19*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v13,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v14,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v15,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v16,acc1);
	data[8*stride] = acc1;

	// upper band :: ix=8
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v14,v20));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v15,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v16,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v17,acc1);
	__m256 v72 = data[72*stride];
	data[72*stride] = acc1;

	// lower band :: ix=9
	__m256 v22 = data[22*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v14,v22));
	__m256 v21 = data[21*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v15,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v16,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v17,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v18,acc1);
	data[9*stride] = acc1;

	// upper band :: ix=9
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v16,v22));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v17,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v18,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v19,acc1);
	__m256 v73 = data[73*stride];
	data[73*stride] = acc1;

	// lower band :: ix=10
	__m256 v24 = data[24*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v16,v24));
	__m256 v23 = data[23*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v17,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v18,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v19,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v20,acc1);
	data[10*stride] = acc1;

	// upper band :: ix=10
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v18,v24));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v19,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v20,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v21,acc1);
	__m256 v74 = data[74*stride];
	data[74*stride] = acc1;

	// lower band :: ix=11
	__m256 v26 = data[26*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v18,v26));
	__m256 v25 = data[25*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v19,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v20,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v21,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v22,acc1);
	data[11*stride] = acc1;

	// upper band :: ix=11
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v20,v26));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v21,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v22,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v23,acc1);
	__m256 v75 = data[75*stride];
	data[75*stride] = acc1;

	// lower band :: ix=12
	__m256 v28 = data[28*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v20,v28));
	__m256 v27 = data[27*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v21,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v22,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v23,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v24,acc1);
	data[12*stride] = acc1;

	// upper band :: ix=12
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v22,v28));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v23,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v24,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v25,acc1);
	__m256 v76 = data[76*stride];
	data[76*stride] = acc1;

	// lower band :: ix=13
	__m256 v30 = data[30*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v22,v30));
	__m256 v29 = data[29*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v23,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v24,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v25,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v26,acc1);
	data[13*stride] = acc1;

	// upper band :: ix=13
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v24,v30));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v25,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v26,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v27,acc1);
	__m256 v77 = data[77*stride];
	data[77*stride] = acc1;

	// lower band :: ix=14
	__m256 v32 = data[32*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v24,v32));
	__m256 v31 = data[31*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v25,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v26,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v27,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v28,acc1);
	data[14*stride] = acc1;

	// upper band :: ix=14
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v26,v32));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v27,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v28,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v29,acc1);
	__m256 v78 = data[78*stride];
	data[78*stride] = acc1;

	// lower band :: ix=15
	__m256 v34 = data[34*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v26,v34));
	__m256 v33 = data[33*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v27,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v28,v32),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v29,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v30,acc1);
	data[15*stride] = acc1;

	// upper band :: ix=15
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v28,v34));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v29,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v30,v32),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v31,acc1);
	__m256 v79 = data[79*stride];
	data[79*stride] = acc1;

	// lower band :: ix=16
	__m256 v36 = data[36*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v28,v36));
	__m256 v35 = data[35*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v29,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v30,v34),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v31,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v32,acc1);
	data[16*stride] = acc1;

	// upper band :: ix=16
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v30,v36));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v31,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v32,v34),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v33,acc1);
	__m256 v80 = data[80*stride];
	data[80*stride] = acc1;

	// lower band :: ix=17
	__m256 v38 = data[38*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v30,v38));
	__m256 v37 = data[37*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v31,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v32,v36),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v33,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v34,acc1);
	data[17*stride] = acc1;

	// upper band :: ix=17
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v32,v38));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v33,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v34,v36),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v35,acc1);
	__m256 v81 = data[81*stride];
	data[81*stride] = acc1;

	// lower band :: ix=18
	__m256 v40 = data[40*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v32,v40));
	__m256 v39 = data[39*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v33,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v34,v38),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v35,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v36,acc1);
	data[18*stride] = acc1;

	// upper band :: ix=18
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v34,v40));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v35,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v36,v38),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v37,acc1);
	__m256 v82 = data[82*stride];
	data[82*stride] = acc1;

	// lower band :: ix=19
	__m256 v42 = data[42*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v34,v42));
	__m256 v41 = data[41*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v35,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v36,v40),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v37,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v38,acc1);
	data[19*stride] = acc1;

	// upper band :: ix=19
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v36,v42));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v37,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v38,v40),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v39,acc1);
	__m256 v83 = data[83*stride];
	data[83*stride] = acc1;

	// lower band :: ix=20
	__m256 v44 = data[44*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v36,v44));
	__m256 v43 = data[43*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v37,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v38,v42),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v39,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v40,acc1);
	data[20*stride] = acc1;

	// upper band :: ix=20
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v38,v44));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v39,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v40,v42),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v41,acc1);
	__m256 v84 = data[84*stride];
	data[84*stride] = acc1;

	// lower band :: ix=21
	__m256 v46 = data[46*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v38,v46));
	__m256 v45 = data[45*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v39,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v40,v44),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v41,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v42,acc1);
	data[21*stride] = acc1;

	// upper band :: ix=21
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v40,v46));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v41,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v42,v44),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v43,acc1);
	__m256 v85 = data[85*stride];
	data[85*stride] = acc1;

	// lower band :: ix=22
	__m256 v48 = data[48*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v40,v48));
	__m256 v47 = data[47*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v41,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v42,v46),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v43,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v44,acc1);
	data[22*stride] = acc1;

	// upper band :: ix=22
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v42,v48));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v43,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v44,v46),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v45,acc1);
	__m256 v86 = data[86*stride];
	data[86*stride] = acc1;

	// lower band :: ix=23
	__m256 v50 = data[50*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v42,v50));
	__m256 v49 = data[49*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v43,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v44,v48),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v45,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v46,acc1);
	data[23*stride] = acc1;

	// upper band :: ix=23
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v44,v50));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v45,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v46,v48),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v47,acc1);
	__m256 v87 = data[87*stride];
	data[87*stride] = acc1;

	// lower band :: ix=24
	__m256 v52 = data[52*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v44,v52));
	__m256 v51 = data[51*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v45,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v46,v50),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v47,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v48,acc1);
	data[24*stride] = acc1;

	// upper band :: ix=24
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v46,v52));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v47,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v48,v50),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v49,acc1);
	__m256 v88 = data[88*stride];
	data[88*stride] = acc1;

	// lower band :: ix=25
	__m256 v54 = data[54*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v46,v54));
	__m256 v53 = data[53*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v47,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v48,v52),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v49,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v50,acc1);
	data[25*stride] = acc1;

	// upper band :: ix=25
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v48,v54));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v49,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v50,v52),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v51,acc1);
	__m256 v89 = data[89*stride];
	data[89*stride] = acc1;

	// lower band :: ix=26
	__m256 v56 = data[56*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v48,v56));
	__m256 v55 = data[55*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v49,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v50,v54),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v51,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v52,acc1);
	data[26*stride] = acc1;

	// upper band :: ix=26
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v50,v56));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v51,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v52,v54),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v53,acc1);
	__m256 v90 = data[90*stride];
	data[90*stride] = acc1;

	// lower band :: ix=27
	__m256 v58 = data[58*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v50,v58));
	__m256 v57 = data[57*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v51,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v52,v56),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v53,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v54,acc1);
	data[27*stride] = acc1;

	// upper band :: ix=27
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v52,v58));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v53,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v54,v56),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v55,acc1);
	__m256 v91 = data[91*stride];
	data[91*stride] = acc1;

	// lower band :: ix=28
	__m256 v60 = data[60*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v52,v60));
	__m256 v59 = data[59*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v53,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v54,v58),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v55,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v56,acc1);
	data[28*stride] = acc1;

	// upper band :: ix=28
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v54,v60));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v55,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v56,v58),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v57,acc1);
	__m256 v92 = data[92*stride];
	data[92*stride] = acc1;

	// lower band :: ix=29
	__m256 v62 = data[62*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v54,v62));
	__m256 v61 = data[61*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v55,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v56,v60),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v57,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v58,acc1);
	data[29*stride] = acc1;

	// upper band :: ix=29
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v56,v62));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v57,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v58,v60),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v59,acc1);
	__m256 v93 = data[93*stride];
	data[93*stride] = acc1;

	// lower band :: ix=30
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v56,v64));
	__m256 v63 = data[63*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v57,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v58,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v59,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v60,acc1);
	data[30*stride] = acc1;

	// upper band :: ix=30
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v58,v64));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v59,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v60,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v61,acc1);
	__m256 v94 = data[94*stride];
	data[94*stride] = acc1;

	// lower band :: ix=31
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v58,v66));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v59,v65),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v60,v64),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v61,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v62,acc1);
	data[31*stride] = acc1;

	// upper band :: ix=31
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v60,v66));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v61,v65),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v62,v64),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v63,acc1);
	__m256 v95 = data[95*stride];
	data[95*stride] = acc1;

	// lower band :: ix=32
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v60,v68));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v61,v67),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v62,v66),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v63,v65),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v64,acc1);
	data[32*stride] = acc1;

	// upper band :: ix=32
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v62,v68));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v63,v67),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v64,v66),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v65,acc1);
	__m256 v96 = data[96*stride];
	data[96*stride] = acc1;

	// lower band :: ix=33
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v62,v70));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v63,v69),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v64,v68),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v65,v67),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v66,acc1);
	data[33*stride] = acc1;

	// upper band :: ix=33
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v64,v70));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v65,v69),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v66,v68),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v67,acc1);
	__m256 v97 = data[97*stride];
	data[97*stride] = acc1;

	// lower band :: ix=34
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v64,v72));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v65,v71),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v66,v70),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v67,v69),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v68,acc1);
	data[34*stride] = acc1;

	// upper band :: ix=34
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v66,v72));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v67,v71),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v68,v70),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v69,acc1);
	__m256 v98 = data[98*stride];
	data[98*stride] = acc1;

	// lower band :: ix=35
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v66,v74));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v67,v73),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v68,v72),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v69,v71),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v70,acc1);
	data[35*stride] = acc1;

	// upper band :: ix=35
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v68,v74));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v69,v73),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v70,v72),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v71,acc1);
	__m256 v99 = data[99*stride];
	data[99*stride] = acc1;

	// lower band :: ix=36
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v68,v76));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v69,v75),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v70,v74),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v71,v73),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v72,acc1);
	data[36*stride] = acc1;

	// upper band :: ix=36
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v70,v76));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v71,v75),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v72,v74),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v73,acc1);
	__m256 v100 = data[100*stride];
	data[100*stride] = acc1;

	// lower band :: ix=37
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v70,v78));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v71,v77),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v72,v76),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v73,v75),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v74,acc1);
	data[37*stride] = acc1;

	// upper band :: ix=37
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v72,v78));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v73,v77),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v74,v76),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v75,acc1);
	__m256 v101 = data[101*stride];
	data[101*stride] = acc1;

	// lower band :: ix=38
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v72,v80));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v73,v79),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v74,v78),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v75,v77),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v76,acc1);
	data[38*stride] = acc1;

	// upper band :: ix=38
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v74,v80));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v75,v79),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v76,v78),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v77,acc1);
	__m256 v102 = data[102*stride];
	data[102*stride] = acc1;

	// lower band :: ix=39
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v74,v82));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v75,v81),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v76,v80),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v77,v79),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v78,acc1);
	data[39*stride] = acc1;

	// upper band :: ix=39
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v76,v82));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v77,v81),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v78,v80),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v79,acc1);
	__m256 v103 = data[103*stride];
	data[103*stride] = acc1;

	// lower band :: ix=40
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v76,v84));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v77,v83),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v78,v82),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v79,v81),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v80,acc1);
	data[40*stride] = acc1;

	// upper band :: ix=40
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v78,v84));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v79,v83),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v80,v82),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v81,acc1);
	__m256 v104 = data[104*stride];
	data[104*stride] = acc1;

	// lower band :: ix=41
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v78,v86));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v79,v85),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v80,v84),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v81,v83),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v82,acc1);
	data[41*stride] = acc1;

	// upper band :: ix=41
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v80,v86));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v81,v85),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v82,v84),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v83,acc1);
	__m256 v105 = data[105*stride];
	data[105*stride] = acc1;

	// lower band :: ix=42
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v80,v88));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v81,v87),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v82,v86),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v83,v85),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v84,acc1);
	data[42*stride] = acc1;

	// upper band :: ix=42
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v82,v88));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v83,v87),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v84,v86),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v85,acc1);
	__m256 v106 = data[106*stride];
	data[106*stride] = acc1;

	// lower band :: ix=43
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v82,v90));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v83,v89),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v84,v88),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v85,v87),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v86,acc1);
	data[43*stride] = acc1;

	// upper band :: ix=43
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v84,v90));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v85,v89),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v86,v88),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v87,acc1);
	__m256 v107 = data[107*stride];
	data[107*stride] = acc1;

	// lower band :: ix=44
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v84,v92));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v85,v91),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v86,v90),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v87,v89),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v88,acc1);
	data[44*stride] = acc1;

	// upper band :: ix=44
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v86,v92));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v87,v91),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v88,v90),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v89,acc1);
	__m256 v108 = data[108*stride];
	data[108*stride] = acc1;

	// lower band :: ix=45
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v86,v94));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v87,v93),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v88,v92),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v89,v91),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v90,acc1);
	data[45*stride] = acc1;

	// upper band :: ix=45
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v88,v94));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v89,v93),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v90,v92),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v91,acc1);
	__m256 v109 = data[109*stride];
	data[109*stride] = acc1;

	// lower band :: ix=46
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v88,v96));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v89,v95),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v90,v94),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v91,v93),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v92,acc1);
	data[46*stride] = acc1;

	// upper band :: ix=46
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v90,v96));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v91,v95),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v92,v94),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v93,acc1);
	__m256 v110 = data[110*stride];
	data[110*stride] = acc1;

	// lower band :: ix=47
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v90,v98));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v91,v97),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v92,v96),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v93,v95),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v94,acc1);
	data[47*stride] = acc1;

	// upper band :: ix=47
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v92,v98));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v93,v97),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v94,v96),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v95,acc1);
	__m256 v111 = data[111*stride];
	data[111*stride] = acc1;

	// lower band :: ix=48
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v92,v100));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v93,v99),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v94,v98),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v95,v97),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v96,acc1);
	data[48*stride] = acc1;

	// upper band :: ix=48
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v94,v100));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v95,v99),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v96,v98),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v97,acc1);
	__m256 v112 = data[112*stride];
	data[112*stride] = acc1;

	// lower band :: ix=49
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v94,v102));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v95,v101),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v96,v100),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v97,v99),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v98,acc1);
	data[49*stride] = acc1;

	// upper band :: ix=49
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v96,v102));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v97,v101),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v98,v100),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v99,acc1);
	__m256 v113 = data[113*stride];
	data[113*stride] = acc1;

	// lower band :: ix=50
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v96,v104));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v97,v103),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v98,v102),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v99,v101),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v100,acc1);
	data[50*stride] = acc1;

	// upper band :: ix=50
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v98,v104));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v99,v103),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v100,v102),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v101,acc1);
	__m256 v114 = data[114*stride];
	data[114*stride] = acc1;

	// lower band :: ix=51
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v98,v106));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v99,v105),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v100,v104),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v101,v103),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v102,acc1);
	data[51*stride] = acc1;

	// upper band :: ix=51
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v100,v106));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v101,v105),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v102,v104),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v103,acc1);
	__m256 v115 = data[115*stride];
	data[115*stride] = acc1;

	// lower band :: ix=52
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v100,v108));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v101,v107),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v102,v106),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v103,v105),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v104,acc1);
	data[52*stride] = acc1;

	// upper band :: ix=52
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v102,v108));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v103,v107),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v104,v106),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v105,acc1);
	__m256 v116 = data[116*stride];
	data[116*stride] = acc1;

	// lower band :: ix=53
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v102,v110));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v103,v109),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v104,v108),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v105,v107),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v106,acc1);
	data[53*stride] = acc1;

	// upper band :: ix=53
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v104,v110));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v105,v109),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v106,v108),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v107,acc1);
	__m256 v117 = data[117*stride];
	data[117*stride] = acc1;

	// lower band :: ix=54
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v104,v112));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v105,v111),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v106,v110),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v107,v109),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v108,acc1);
	data[54*stride] = acc1;

	// upper band :: ix=54
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v106,v112));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v107,v111),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v108,v110),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v109,acc1);
	__m256 v118 = data[118*stride];
	data[118*stride] = acc1;

	// lower band :: ix=55
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v106,v114));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v107,v113),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v108,v112),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v109,v111),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v110,acc1);
	data[55*stride] = acc1;

	// upper band :: ix=55
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v108,v114));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v109,v113),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v110,v112),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v111,acc1);
	__m256 v119 = data[119*stride];
	data[119*stride] = acc1;

	// lower band :: ix=56
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v108,v116));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v109,v115),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v110,v114),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v111,v113),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v112,acc1);
	data[56*stride] = acc1;

	// upper band :: ix=56
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v110,v116));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v111,v115),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v112,v114),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v113,acc1);
	__m256 v120 = data[120*stride];
	data[120*stride] = acc1;

	// lower band :: ix=57
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v110,v118));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v111,v117),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v112,v116),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v113,v115),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v114,acc1);
	data[57*stride] = acc1;

	// upper band :: ix=57
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v112,v118));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v113,v117),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v114,v116),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v115,acc1);
	__m256 v121 = data[121*stride];
	data[121*stride] = acc1;

	// lower band :: ix=58
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v112,v120));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v113,v119),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v114,v118),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v115,v117),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v116,acc1);
	data[58*stride] = acc1;

	// upper band :: ix=58
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v114,v120));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v115,v119),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v116,v118),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v117,acc1);
	__m256 v122 = data[122*stride];
	data[122*stride] = acc1;

	// lower band :: ix=59
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v114,v122));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v115,v121),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v116,v120),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v117,v119),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v118,acc1);
	data[59*stride] = acc1;

	// upper band :: ix=59
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v116,v122));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v117,v121),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v118,v120),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v119,acc1);
	__m256 v123 = data[123*stride];
	data[123*stride] = acc1;

	// lower band :: ix=60
	__m256 v124 = data[124*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v116,v124));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v117,v123),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v118,v122),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v119,v121),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v120,acc1);
	data[60*stride] = acc1;

	// upper band :: ix=60
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v118,v124));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v119,v123),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v120,v122),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v121,acc1);
	data[124*stride] = acc1;

	// lower band :: ix=61
	__m256 v126 = data[126*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v118,v126));
	__m256 v125 = data[125*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v119,v125),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v120,v124),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v121,v123),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v122,acc1);
	data[61*stride] = acc1;

	// upper band :: ix=61
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v120,v126));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v121,v125),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v122,v124),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v123,acc1);
	data[125*stride] = acc1;

	// lower band :: ix=62
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v120,v126));
	__m256 v127 = data[127*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v121,v127),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v122,v126),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v123,v125),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v124,acc1);
	data[62*stride] = acc1;

	// upper band :: ix=62
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v122,v126));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v123,v127),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v124,v126),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v125,acc1);
	data[126*stride] = acc1;

	// lower band :: ix=63
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v122,v124));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v123,v125),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v124,v126),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v125,v127),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v126,acc1);
	data[63*stride] = acc1;

	// upper band :: ix=63
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v124,v124));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v125,v125),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v126,v126),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v127,acc1);
	data[127*stride] = acc1;
}

static void _Ds79_AVX_256(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v4 = data[4*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v4));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v3,v3),acc1);
	__m256 v2 = data[2*stride];
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v2,v2),acc1);
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v1),acc1);
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_al0, v0,acc1);
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v4));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v0,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v1,acc1);
	__m256 v128 = data[128*stride];
	data[128*stride] = acc1;

	// lower band :: ix=1
	__m256 v6 = data[6*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v6));
	__m256 v5 = data[5*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v1,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v0,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v2,acc1);
	data[1*stride] = acc1;

	// upper band :: ix=1
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v6));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v1,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v2,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v3,acc1);
	__m256 v129 = data[129*stride];
	data[129*stride] = acc1;

	// lower band :: ix=2
	__m256 v8 = data[8*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v0,v8));
	__m256 v7 = data[7*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v1,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v2,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v3,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v4,acc1);
	data[2*stride] = acc1;

	// upper band :: ix=2
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v8));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v3,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v4,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v5,acc1);
	__m256 v130 = data[130*stride];
	data[130*stride] = acc1;

	// lower band :: ix=3
	__m256 v10 = data[10*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v10));
	__m256 v9 = data[9*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v3,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v4,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v5,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v6,acc1);
	data[3*stride] = acc1;

	// upper band :: ix=3
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v4,v10));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v5,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v6,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v7,acc1);
	__m256 v131 = data[131*stride];
	data[131*stride] = acc1;

	// lower band :: ix=4
	__m256 v12 = data[12*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v12));
	__m256 v11 = data[11*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v5,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v6,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v7,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v8,acc1);
	data[4*stride] = acc1;

	// upper band :: ix=4
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v6,v12));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v7,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v8,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v9,acc1);
	__m256 v132 = data[132*stride];
	data[132*stride] = acc1;

	// lower band :: ix=5
	__m256 v14 = data[14*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v6,v14));
	__m256 v13 = data[13*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v7,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v8,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v9,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v10,acc1);
	data[5*stride] = acc1;

	// upper band :: ix=5
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v8,v14));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v9,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v10,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v11,acc1);
	__m256 v133 = data[133*stride];
	data[133*stride] = acc1;

	// lower band :: ix=6
	__m256 v16 = data[16*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v8,v16));
	__m256 v15 = data[15*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v9,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v10,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v11,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v12,acc1);
	data[6*stride] = acc1;

	// upper band :: ix=6
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v10,v16));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v11,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v12,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v13,acc1);
	__m256 v134 = data[134*stride];
	data[134*stride] = acc1;

	// lower band :: ix=7
	__m256 v18 = data[18*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v10,v18));
	__m256 v17 = data[17*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v11,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v12,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v13,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v14,acc1);
	data[7*stride] = acc1;

	// upper band :: ix=7
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v12,v18));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v13,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v14,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v15,acc1);
	__m256 v135 = data[135*stride];
	data[135*stride] = acc1;

	// lower band :: ix=8
	__m256 v20 = data[20*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v12,v20));
	__m256 v19 = data[19*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v13,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v14,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v15,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v16,acc1);
	data[8*stride] = acc1;

	// upper band :: ix=8
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v14,v20));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v15,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v16,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v17,acc1);
	__m256 v136 = data[136*stride];
	data[136*stride] = acc1;

	// lower band :: ix=9
	__m256 v22 = data[22*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v14,v22));
	__m256 v21 = data[21*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v15,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v16,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v17,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v18,acc1);
	data[9*stride] = acc1;

	// upper band :: ix=9
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v16,v22));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v17,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v18,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v19,acc1);
	__m256 v137 = data[137*stride];
	data[137*stride] = acc1;

	// lower band :: ix=10
	__m256 v24 = data[24*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v16,v24));
	__m256 v23 = data[23*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v17,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v18,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v19,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v20,acc1);
	data[10*stride] = acc1;

	// upper band :: ix=10
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v18,v24));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v19,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v20,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v21,acc1);
	__m256 v138 = data[138*stride];
	data[138*stride] = acc1;

	// lower band :: ix=11
	__m256 v26 = data[26*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v18,v26));
	__m256 v25 = data[25*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v19,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v20,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v21,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v22,acc1);
	data[11*stride] = acc1;

	// upper band :: ix=11
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v20,v26));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v21,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v22,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v23,acc1);
	__m256 v139 = data[139*stride];
	data[139*stride] = acc1;

	// lower band :: ix=12
	__m256 v28 = data[28*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v20,v28));
	__m256 v27 = data[27*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v21,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v22,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v23,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v24,acc1);
	data[12*stride] = acc1;

	// upper band :: ix=12
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v22,v28));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v23,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v24,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v25,acc1);
	__m256 v140 = data[140*stride];
	data[140*stride] = acc1;

	// lower band :: ix=13
	__m256 v30 = data[30*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v22,v30));
	__m256 v29 = data[29*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v23,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v24,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v25,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v26,acc1);
	data[13*stride] = acc1;

	// upper band :: ix=13
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v24,v30));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v25,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v26,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v27,acc1);
	__m256 v141 = data[141*stride];
	data[141*stride] = acc1;

	// lower band :: ix=14
	__m256 v32 = data[32*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v24,v32));
	__m256 v31 = data[31*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v25,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v26,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v27,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v28,acc1);
	data[14*stride] = acc1;

	// upper band :: ix=14
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v26,v32));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v27,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v28,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v29,acc1);
	__m256 v142 = data[142*stride];
	data[142*stride] = acc1;

	// lower band :: ix=15
	__m256 v34 = data[34*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v26,v34));
	__m256 v33 = data[33*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v27,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v28,v32),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v29,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v30,acc1);
	data[15*stride] = acc1;

	// upper band :: ix=15
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v28,v34));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v29,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v30,v32),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v31,acc1);
	__m256 v143 = data[143*stride];
	data[143*stride] = acc1;

	// lower band :: ix=16
	__m256 v36 = data[36*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v28,v36));
	__m256 v35 = data[35*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v29,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v30,v34),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v31,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v32,acc1);
	data[16*stride] = acc1;

	// upper band :: ix=16
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v30,v36));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v31,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v32,v34),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v33,acc1);
	__m256 v144 = data[144*stride];
	data[144*stride] = acc1;

	// lower band :: ix=17
	__m256 v38 = data[38*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v30,v38));
	__m256 v37 = data[37*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v31,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v32,v36),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v33,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v34,acc1);
	data[17*stride] = acc1;

	// upper band :: ix=17
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v32,v38));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v33,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v34,v36),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v35,acc1);
	__m256 v145 = data[145*stride];
	data[145*stride] = acc1;

	// lower band :: ix=18
	__m256 v40 = data[40*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v32,v40));
	__m256 v39 = data[39*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v33,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v34,v38),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v35,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v36,acc1);
	data[18*stride] = acc1;

	// upper band :: ix=18
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v34,v40));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v35,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v36,v38),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v37,acc1);
	__m256 v146 = data[146*stride];
	data[146*stride] = acc1;

	// lower band :: ix=19
	__m256 v42 = data[42*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v34,v42));
	__m256 v41 = data[41*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v35,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v36,v40),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v37,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v38,acc1);
	data[19*stride] = acc1;

	// upper band :: ix=19
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v36,v42));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v37,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v38,v40),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v39,acc1);
	__m256 v147 = data[147*stride];
	data[147*stride] = acc1;

	// lower band :: ix=20
	__m256 v44 = data[44*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v36,v44));
	__m256 v43 = data[43*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v37,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v38,v42),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v39,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v40,acc1);
	data[20*stride] = acc1;

	// upper band :: ix=20
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v38,v44));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v39,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v40,v42),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v41,acc1);
	__m256 v148 = data[148*stride];
	data[148*stride] = acc1;

	// lower band :: ix=21
	__m256 v46 = data[46*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v38,v46));
	__m256 v45 = data[45*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v39,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v40,v44),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v41,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v42,acc1);
	data[21*stride] = acc1;

	// upper band :: ix=21
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v40,v46));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v41,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v42,v44),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v43,acc1);
	__m256 v149 = data[149*stride];
	data[149*stride] = acc1;

	// lower band :: ix=22
	__m256 v48 = data[48*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v40,v48));
	__m256 v47 = data[47*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v41,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v42,v46),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v43,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v44,acc1);
	data[22*stride] = acc1;

	// upper band :: ix=22
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v42,v48));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v43,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v44,v46),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v45,acc1);
	__m256 v150 = data[150*stride];
	data[150*stride] = acc1;

	// lower band :: ix=23
	__m256 v50 = data[50*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v42,v50));
	__m256 v49 = data[49*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v43,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v44,v48),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v45,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v46,acc1);
	data[23*stride] = acc1;

	// upper band :: ix=23
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v44,v50));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v45,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v46,v48),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v47,acc1);
	__m256 v151 = data[151*stride];
	data[151*stride] = acc1;

	// lower band :: ix=24
	__m256 v52 = data[52*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v44,v52));
	__m256 v51 = data[51*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v45,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v46,v50),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v47,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v48,acc1);
	data[24*stride] = acc1;

	// upper band :: ix=24
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v46,v52));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v47,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v48,v50),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v49,acc1);
	__m256 v152 = data[152*stride];
	data[152*stride] = acc1;

	// lower band :: ix=25
	__m256 v54 = data[54*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v46,v54));
	__m256 v53 = data[53*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v47,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v48,v52),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v49,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v50,acc1);
	data[25*stride] = acc1;

	// upper band :: ix=25
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v48,v54));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v49,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v50,v52),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v51,acc1);
	__m256 v153 = data[153*stride];
	data[153*stride] = acc1;

	// lower band :: ix=26
	__m256 v56 = data[56*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v48,v56));
	__m256 v55 = data[55*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v49,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v50,v54),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v51,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v52,acc1);
	data[26*stride] = acc1;

	// upper band :: ix=26
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v50,v56));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v51,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v52,v54),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v53,acc1);
	__m256 v154 = data[154*stride];
	data[154*stride] = acc1;

	// lower band :: ix=27
	__m256 v58 = data[58*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v50,v58));
	__m256 v57 = data[57*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v51,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v52,v56),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v53,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v54,acc1);
	data[27*stride] = acc1;

	// upper band :: ix=27
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v52,v58));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v53,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v54,v56),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v55,acc1);
	__m256 v155 = data[155*stride];
	data[155*stride] = acc1;

	// lower band :: ix=28
	__m256 v60 = data[60*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v52,v60));
	__m256 v59 = data[59*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v53,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v54,v58),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v55,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v56,acc1);
	data[28*stride] = acc1;

	// upper band :: ix=28
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v54,v60));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v55,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v56,v58),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v57,acc1);
	__m256 v156 = data[156*stride];
	data[156*stride] = acc1;

	// lower band :: ix=29
	__m256 v62 = data[62*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v54,v62));
	__m256 v61 = data[61*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v55,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v56,v60),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v57,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v58,acc1);
	data[29*stride] = acc1;

	// upper band :: ix=29
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v56,v62));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v57,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v58,v60),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v59,acc1);
	__m256 v157 = data[157*stride];
	data[157*stride] = acc1;

	// lower band :: ix=30
	__m256 v64 = data[64*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v56,v64));
	__m256 v63 = data[63*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v57,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v58,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v59,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v60,acc1);
	data[30*stride] = acc1;

	// upper band :: ix=30
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v58,v64));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v59,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v60,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v61,acc1);
	__m256 v158 = data[158*stride];
	data[158*stride] = acc1;

	// lower band :: ix=31
	__m256 v66 = data[66*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v58,v66));
	__m256 v65 = data[65*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v59,v65),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v60,v64),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v61,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v62,acc1);
	data[31*stride] = acc1;

	// upper band :: ix=31
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v60,v66));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v61,v65),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v62,v64),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v63,acc1);
	__m256 v159 = data[159*stride];
	data[159*stride] = acc1;

	// lower band :: ix=32
	__m256 v68 = data[68*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v60,v68));
	__m256 v67 = data[67*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v61,v67),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v62,v66),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v63,v65),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v64,acc1);
	data[32*stride] = acc1;

	// upper band :: ix=32
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v62,v68));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v63,v67),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v64,v66),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v65,acc1);
	__m256 v160 = data[160*stride];
	data[160*stride] = acc1;

	// lower band :: ix=33
	__m256 v70 = data[70*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v62,v70));
	__m256 v69 = data[69*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v63,v69),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v64,v68),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v65,v67),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v66,acc1);
	data[33*stride] = acc1;

	// upper band :: ix=33
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v64,v70));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v65,v69),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v66,v68),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v67,acc1);
	__m256 v161 = data[161*stride];
	data[161*stride] = acc1;

	// lower band :: ix=34
	__m256 v72 = data[72*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v64,v72));
	__m256 v71 = data[71*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v65,v71),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v66,v70),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v67,v69),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v68,acc1);
	data[34*stride] = acc1;

	// upper band :: ix=34
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v66,v72));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v67,v71),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v68,v70),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v69,acc1);
	__m256 v162 = data[162*stride];
	data[162*stride] = acc1;

	// lower band :: ix=35
	__m256 v74 = data[74*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v66,v74));
	__m256 v73 = data[73*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v67,v73),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v68,v72),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v69,v71),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v70,acc1);
	data[35*stride] = acc1;

	// upper band :: ix=35
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v68,v74));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v69,v73),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v70,v72),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v71,acc1);
	__m256 v163 = data[163*stride];
	data[163*stride] = acc1;

	// lower band :: ix=36
	__m256 v76 = data[76*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v68,v76));
	__m256 v75 = data[75*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v69,v75),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v70,v74),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v71,v73),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v72,acc1);
	data[36*stride] = acc1;

	// upper band :: ix=36
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v70,v76));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v71,v75),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v72,v74),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v73,acc1);
	__m256 v164 = data[164*stride];
	data[164*stride] = acc1;

	// lower band :: ix=37
	__m256 v78 = data[78*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v70,v78));
	__m256 v77 = data[77*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v71,v77),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v72,v76),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v73,v75),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v74,acc1);
	data[37*stride] = acc1;

	// upper band :: ix=37
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v72,v78));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v73,v77),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v74,v76),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v75,acc1);
	__m256 v165 = data[165*stride];
	data[165*stride] = acc1;

	// lower band :: ix=38
	__m256 v80 = data[80*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v72,v80));
	__m256 v79 = data[79*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v73,v79),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v74,v78),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v75,v77),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v76,acc1);
	data[38*stride] = acc1;

	// upper band :: ix=38
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v74,v80));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v75,v79),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v76,v78),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v77,acc1);
	__m256 v166 = data[166*stride];
	data[166*stride] = acc1;

	// lower band :: ix=39
	__m256 v82 = data[82*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v74,v82));
	__m256 v81 = data[81*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v75,v81),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v76,v80),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v77,v79),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v78,acc1);
	data[39*stride] = acc1;

	// upper band :: ix=39
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v76,v82));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v77,v81),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v78,v80),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v79,acc1);
	__m256 v167 = data[167*stride];
	data[167*stride] = acc1;

	// lower band :: ix=40
	__m256 v84 = data[84*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v76,v84));
	__m256 v83 = data[83*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v77,v83),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v78,v82),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v79,v81),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v80,acc1);
	data[40*stride] = acc1;

	// upper band :: ix=40
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v78,v84));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v79,v83),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v80,v82),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v81,acc1);
	__m256 v168 = data[168*stride];
	data[168*stride] = acc1;

	// lower band :: ix=41
	__m256 v86 = data[86*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v78,v86));
	__m256 v85 = data[85*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v79,v85),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v80,v84),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v81,v83),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v82,acc1);
	data[41*stride] = acc1;

	// upper band :: ix=41
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v80,v86));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v81,v85),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v82,v84),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v83,acc1);
	__m256 v169 = data[169*stride];
	data[169*stride] = acc1;

	// lower band :: ix=42
	__m256 v88 = data[88*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v80,v88));
	__m256 v87 = data[87*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v81,v87),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v82,v86),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v83,v85),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v84,acc1);
	data[42*stride] = acc1;

	// upper band :: ix=42
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v82,v88));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v83,v87),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v84,v86),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v85,acc1);
	__m256 v170 = data[170*stride];
	data[170*stride] = acc1;

	// lower band :: ix=43
	__m256 v90 = data[90*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v82,v90));
	__m256 v89 = data[89*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v83,v89),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v84,v88),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v85,v87),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v86,acc1);
	data[43*stride] = acc1;

	// upper band :: ix=43
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v84,v90));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v85,v89),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v86,v88),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v87,acc1);
	__m256 v171 = data[171*stride];
	data[171*stride] = acc1;

	// lower band :: ix=44
	__m256 v92 = data[92*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v84,v92));
	__m256 v91 = data[91*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v85,v91),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v86,v90),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v87,v89),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v88,acc1);
	data[44*stride] = acc1;

	// upper band :: ix=44
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v86,v92));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v87,v91),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v88,v90),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v89,acc1);
	__m256 v172 = data[172*stride];
	data[172*stride] = acc1;

	// lower band :: ix=45
	__m256 v94 = data[94*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v86,v94));
	__m256 v93 = data[93*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v87,v93),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v88,v92),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v89,v91),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v90,acc1);
	data[45*stride] = acc1;

	// upper band :: ix=45
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v88,v94));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v89,v93),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v90,v92),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v91,acc1);
	__m256 v173 = data[173*stride];
	data[173*stride] = acc1;

	// lower band :: ix=46
	__m256 v96 = data[96*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v88,v96));
	__m256 v95 = data[95*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v89,v95),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v90,v94),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v91,v93),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v92,acc1);
	data[46*stride] = acc1;

	// upper band :: ix=46
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v90,v96));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v91,v95),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v92,v94),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v93,acc1);
	__m256 v174 = data[174*stride];
	data[174*stride] = acc1;

	// lower band :: ix=47
	__m256 v98 = data[98*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v90,v98));
	__m256 v97 = data[97*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v91,v97),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v92,v96),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v93,v95),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v94,acc1);
	data[47*stride] = acc1;

	// upper band :: ix=47
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v92,v98));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v93,v97),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v94,v96),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v95,acc1);
	__m256 v175 = data[175*stride];
	data[175*stride] = acc1;

	// lower band :: ix=48
	__m256 v100 = data[100*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v92,v100));
	__m256 v99 = data[99*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v93,v99),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v94,v98),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v95,v97),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v96,acc1);
	data[48*stride] = acc1;

	// upper band :: ix=48
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v94,v100));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v95,v99),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v96,v98),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v97,acc1);
	__m256 v176 = data[176*stride];
	data[176*stride] = acc1;

	// lower band :: ix=49
	__m256 v102 = data[102*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v94,v102));
	__m256 v101 = data[101*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v95,v101),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v96,v100),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v97,v99),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v98,acc1);
	data[49*stride] = acc1;

	// upper band :: ix=49
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v96,v102));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v97,v101),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v98,v100),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v99,acc1);
	__m256 v177 = data[177*stride];
	data[177*stride] = acc1;

	// lower band :: ix=50
	__m256 v104 = data[104*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v96,v104));
	__m256 v103 = data[103*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v97,v103),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v98,v102),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v99,v101),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v100,acc1);
	data[50*stride] = acc1;

	// upper band :: ix=50
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v98,v104));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v99,v103),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v100,v102),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v101,acc1);
	__m256 v178 = data[178*stride];
	data[178*stride] = acc1;

	// lower band :: ix=51
	__m256 v106 = data[106*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v98,v106));
	__m256 v105 = data[105*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v99,v105),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v100,v104),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v101,v103),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v102,acc1);
	data[51*stride] = acc1;

	// upper band :: ix=51
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v100,v106));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v101,v105),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v102,v104),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v103,acc1);
	__m256 v179 = data[179*stride];
	data[179*stride] = acc1;

	// lower band :: ix=52
	__m256 v108 = data[108*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v100,v108));
	__m256 v107 = data[107*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v101,v107),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v102,v106),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v103,v105),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v104,acc1);
	data[52*stride] = acc1;

	// upper band :: ix=52
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v102,v108));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v103,v107),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v104,v106),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v105,acc1);
	__m256 v180 = data[180*stride];
	data[180*stride] = acc1;

	// lower band :: ix=53
	__m256 v110 = data[110*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v102,v110));
	__m256 v109 = data[109*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v103,v109),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v104,v108),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v105,v107),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v106,acc1);
	data[53*stride] = acc1;

	// upper band :: ix=53
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v104,v110));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v105,v109),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v106,v108),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v107,acc1);
	__m256 v181 = data[181*stride];
	data[181*stride] = acc1;

	// lower band :: ix=54
	__m256 v112 = data[112*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v104,v112));
	__m256 v111 = data[111*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v105,v111),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v106,v110),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v107,v109),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v108,acc1);
	data[54*stride] = acc1;

	// upper band :: ix=54
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v106,v112));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v107,v111),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v108,v110),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v109,acc1);
	__m256 v182 = data[182*stride];
	data[182*stride] = acc1;

	// lower band :: ix=55
	__m256 v114 = data[114*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v106,v114));
	__m256 v113 = data[113*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v107,v113),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v108,v112),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v109,v111),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v110,acc1);
	data[55*stride] = acc1;

	// upper band :: ix=55
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v108,v114));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v109,v113),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v110,v112),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v111,acc1);
	__m256 v183 = data[183*stride];
	data[183*stride] = acc1;

	// lower band :: ix=56
	__m256 v116 = data[116*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v108,v116));
	__m256 v115 = data[115*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v109,v115),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v110,v114),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v111,v113),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v112,acc1);
	data[56*stride] = acc1;

	// upper band :: ix=56
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v110,v116));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v111,v115),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v112,v114),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v113,acc1);
	__m256 v184 = data[184*stride];
	data[184*stride] = acc1;

	// lower band :: ix=57
	__m256 v118 = data[118*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v110,v118));
	__m256 v117 = data[117*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v111,v117),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v112,v116),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v113,v115),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v114,acc1);
	data[57*stride] = acc1;

	// upper band :: ix=57
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v112,v118));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v113,v117),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v114,v116),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v115,acc1);
	__m256 v185 = data[185*stride];
	data[185*stride] = acc1;

	// lower band :: ix=58
	__m256 v120 = data[120*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v112,v120));
	__m256 v119 = data[119*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v113,v119),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v114,v118),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v115,v117),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v116,acc1);
	data[58*stride] = acc1;

	// upper band :: ix=58
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v114,v120));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v115,v119),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v116,v118),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v117,acc1);
	__m256 v186 = data[186*stride];
	data[186*stride] = acc1;

	// lower band :: ix=59
	__m256 v122 = data[122*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v114,v122));
	__m256 v121 = data[121*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v115,v121),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v116,v120),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v117,v119),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v118,acc1);
	data[59*stride] = acc1;

	// upper band :: ix=59
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v116,v122));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v117,v121),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v118,v120),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v119,acc1);
	__m256 v187 = data[187*stride];
	data[187*stride] = acc1;

	// lower band :: ix=60
	__m256 v124 = data[124*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v116,v124));
	__m256 v123 = data[123*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v117,v123),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v118,v122),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v119,v121),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v120,acc1);
	data[60*stride] = acc1;

	// upper band :: ix=60
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v118,v124));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v119,v123),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v120,v122),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v121,acc1);
	__m256 v188 = data[188*stride];
	data[188*stride] = acc1;

	// lower band :: ix=61
	__m256 v126 = data[126*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v118,v126));
	__m256 v125 = data[125*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v119,v125),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v120,v124),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v121,v123),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v122,acc1);
	data[61*stride] = acc1;

	// upper band :: ix=61
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v120,v126));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v121,v125),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v122,v124),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v123,acc1);
	__m256 v189 = data[189*stride];
	data[189*stride] = acc1;

	// lower band :: ix=62
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v120,v128));
	__m256 v127 = data[127*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v121,v127),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v122,v126),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v123,v125),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v124,acc1);
	data[62*stride] = acc1;

	// upper band :: ix=62
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v122,v128));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v123,v127),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v124,v126),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v125,acc1);
	__m256 v190 = data[190*stride];
	data[190*stride] = acc1;

	// lower band :: ix=63
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v122,v130));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v123,v129),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v124,v128),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v125,v127),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v126,acc1);
	data[63*stride] = acc1;

	// upper band :: ix=63
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v124,v130));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v125,v129),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v126,v128),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v127,acc1);
	__m256 v191 = data[191*stride];
	data[191*stride] = acc1;

	// lower band :: ix=64
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v124,v132));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v125,v131),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v126,v130),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v127,v129),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v128,acc1);
	data[64*stride] = acc1;

	// upper band :: ix=64
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v126,v132));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v127,v131),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v128,v130),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v129,acc1);
	__m256 v192 = data[192*stride];
	data[192*stride] = acc1;

	// lower band :: ix=65
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v126,v134));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v127,v133),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v128,v132),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v129,v131),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v130,acc1);
	data[65*stride] = acc1;

	// upper band :: ix=65
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v128,v134));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v129,v133),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v130,v132),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v131,acc1);
	__m256 v193 = data[193*stride];
	data[193*stride] = acc1;

	// lower band :: ix=66
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v128,v136));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v129,v135),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v130,v134),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v131,v133),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v132,acc1);
	data[66*stride] = acc1;

	// upper band :: ix=66
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v130,v136));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v131,v135),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v132,v134),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v133,acc1);
	__m256 v194 = data[194*stride];
	data[194*stride] = acc1;

	// lower band :: ix=67
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v130,v138));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v131,v137),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v132,v136),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v133,v135),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v134,acc1);
	data[67*stride] = acc1;

	// upper band :: ix=67
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v132,v138));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v133,v137),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v134,v136),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v135,acc1);
	__m256 v195 = data[195*stride];
	data[195*stride] = acc1;

	// lower band :: ix=68
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v132,v140));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v133,v139),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v134,v138),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v135,v137),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v136,acc1);
	data[68*stride] = acc1;

	// upper band :: ix=68
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v134,v140));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v135,v139),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v136,v138),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v137,acc1);
	__m256 v196 = data[196*stride];
	data[196*stride] = acc1;

	// lower band :: ix=69
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v134,v142));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v135,v141),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v136,v140),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v137,v139),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v138,acc1);
	data[69*stride] = acc1;

	// upper band :: ix=69
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v136,v142));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v137,v141),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v138,v140),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v139,acc1);
	__m256 v197 = data[197*stride];
	data[197*stride] = acc1;

	// lower band :: ix=70
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v136,v144));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v137,v143),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v138,v142),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v139,v141),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v140,acc1);
	data[70*stride] = acc1;

	// upper band :: ix=70
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v138,v144));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v139,v143),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v140,v142),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v141,acc1);
	__m256 v198 = data[198*stride];
	data[198*stride] = acc1;

	// lower band :: ix=71
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v138,v146));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v139,v145),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v140,v144),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v141,v143),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v142,acc1);
	data[71*stride] = acc1;

	// upper band :: ix=71
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v140,v146));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v141,v145),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v142,v144),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v143,acc1);
	__m256 v199 = data[199*stride];
	data[199*stride] = acc1;

	// lower band :: ix=72
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v140,v148));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v141,v147),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v142,v146),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v143,v145),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v144,acc1);
	data[72*stride] = acc1;

	// upper band :: ix=72
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v142,v148));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v143,v147),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v144,v146),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v145,acc1);
	__m256 v200 = data[200*stride];
	data[200*stride] = acc1;

	// lower band :: ix=73
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v142,v150));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v143,v149),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v144,v148),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v145,v147),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v146,acc1);
	data[73*stride] = acc1;

	// upper band :: ix=73
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v144,v150));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v145,v149),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v146,v148),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v147,acc1);
	__m256 v201 = data[201*stride];
	data[201*stride] = acc1;

	// lower band :: ix=74
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v144,v152));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v145,v151),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v146,v150),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v147,v149),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v148,acc1);
	data[74*stride] = acc1;

	// upper band :: ix=74
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v146,v152));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v147,v151),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v148,v150),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v149,acc1);
	__m256 v202 = data[202*stride];
	data[202*stride] = acc1;

	// lower band :: ix=75
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v146,v154));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v147,v153),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v148,v152),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v149,v151),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v150,acc1);
	data[75*stride] = acc1;

	// upper band :: ix=75
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v148,v154));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v149,v153),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v150,v152),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v151,acc1);
	__m256 v203 = data[203*stride];
	data[203*stride] = acc1;

	// lower band :: ix=76
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v148,v156));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v149,v155),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v150,v154),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v151,v153),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v152,acc1);
	data[76*stride] = acc1;

	// upper band :: ix=76
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v150,v156));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v151,v155),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v152,v154),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v153,acc1);
	__m256 v204 = data[204*stride];
	data[204*stride] = acc1;

	// lower band :: ix=77
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v150,v158));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v151,v157),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v152,v156),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v153,v155),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v154,acc1);
	data[77*stride] = acc1;

	// upper band :: ix=77
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v152,v158));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v153,v157),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v154,v156),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v155,acc1);
	__m256 v205 = data[205*stride];
	data[205*stride] = acc1;

	// lower band :: ix=78
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v152,v160));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v153,v159),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v154,v158),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v155,v157),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v156,acc1);
	data[78*stride] = acc1;

	// upper band :: ix=78
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v154,v160));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v155,v159),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v156,v158),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v157,acc1);
	__m256 v206 = data[206*stride];
	data[206*stride] = acc1;

	// lower band :: ix=79
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v154,v162));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v155,v161),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v156,v160),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v157,v159),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v158,acc1);
	data[79*stride] = acc1;

	// upper band :: ix=79
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v156,v162));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v157,v161),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v158,v160),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v159,acc1);
	__m256 v207 = data[207*stride];
	data[207*stride] = acc1;

	// lower band :: ix=80
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v156,v164));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v157,v163),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v158,v162),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v159,v161),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v160,acc1);
	data[80*stride] = acc1;

	// upper band :: ix=80
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v158,v164));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v159,v163),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v160,v162),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v161,acc1);
	__m256 v208 = data[208*stride];
	data[208*stride] = acc1;

	// lower band :: ix=81
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v158,v166));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v159,v165),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v160,v164),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v161,v163),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v162,acc1);
	data[81*stride] = acc1;

	// upper band :: ix=81
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v160,v166));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v161,v165),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v162,v164),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v163,acc1);
	__m256 v209 = data[209*stride];
	data[209*stride] = acc1;

	// lower band :: ix=82
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v160,v168));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v161,v167),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v162,v166),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v163,v165),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v164,acc1);
	data[82*stride] = acc1;

	// upper band :: ix=82
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v162,v168));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v163,v167),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v164,v166),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v165,acc1);
	__m256 v210 = data[210*stride];
	data[210*stride] = acc1;

	// lower band :: ix=83
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v162,v170));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v163,v169),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v164,v168),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v165,v167),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v166,acc1);
	data[83*stride] = acc1;

	// upper band :: ix=83
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v164,v170));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v165,v169),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v166,v168),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v167,acc1);
	__m256 v211 = data[211*stride];
	data[211*stride] = acc1;

	// lower band :: ix=84
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v164,v172));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v165,v171),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v166,v170),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v167,v169),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v168,acc1);
	data[84*stride] = acc1;

	// upper band :: ix=84
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v166,v172));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v167,v171),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v168,v170),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v169,acc1);
	__m256 v212 = data[212*stride];
	data[212*stride] = acc1;

	// lower band :: ix=85
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v166,v174));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v167,v173),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v168,v172),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v169,v171),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v170,acc1);
	data[85*stride] = acc1;

	// upper band :: ix=85
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v168,v174));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v169,v173),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v170,v172),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v171,acc1);
	__m256 v213 = data[213*stride];
	data[213*stride] = acc1;

	// lower band :: ix=86
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v168,v176));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v169,v175),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v170,v174),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v171,v173),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v172,acc1);
	data[86*stride] = acc1;

	// upper band :: ix=86
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v170,v176));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v171,v175),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v172,v174),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v173,acc1);
	__m256 v214 = data[214*stride];
	data[214*stride] = acc1;

	// lower band :: ix=87
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v170,v178));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v171,v177),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v172,v176),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v173,v175),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v174,acc1);
	data[87*stride] = acc1;

	// upper band :: ix=87
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v172,v178));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v173,v177),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v174,v176),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v175,acc1);
	__m256 v215 = data[215*stride];
	data[215*stride] = acc1;

	// lower band :: ix=88
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v172,v180));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v173,v179),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v174,v178),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v175,v177),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v176,acc1);
	data[88*stride] = acc1;

	// upper band :: ix=88
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v174,v180));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v175,v179),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v176,v178),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v177,acc1);
	__m256 v216 = data[216*stride];
	data[216*stride] = acc1;

	// lower band :: ix=89
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v174,v182));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v175,v181),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v176,v180),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v177,v179),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v178,acc1);
	data[89*stride] = acc1;

	// upper band :: ix=89
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v176,v182));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v177,v181),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v178,v180),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v179,acc1);
	__m256 v217 = data[217*stride];
	data[217*stride] = acc1;

	// lower band :: ix=90
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v176,v184));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v177,v183),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v178,v182),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v179,v181),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v180,acc1);
	data[90*stride] = acc1;

	// upper band :: ix=90
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v178,v184));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v179,v183),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v180,v182),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v181,acc1);
	__m256 v218 = data[218*stride];
	data[218*stride] = acc1;

	// lower band :: ix=91
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v178,v186));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v179,v185),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v180,v184),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v181,v183),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v182,acc1);
	data[91*stride] = acc1;

	// upper band :: ix=91
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v180,v186));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v181,v185),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v182,v184),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v183,acc1);
	__m256 v219 = data[219*stride];
	data[219*stride] = acc1;

	// lower band :: ix=92
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v180,v188));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v181,v187),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v182,v186),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v183,v185),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v184,acc1);
	data[92*stride] = acc1;

	// upper band :: ix=92
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v182,v188));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v183,v187),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v184,v186),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v185,acc1);
	__m256 v220 = data[220*stride];
	data[220*stride] = acc1;

	// lower band :: ix=93
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v182,v190));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v183,v189),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v184,v188),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v185,v187),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v186,acc1);
	data[93*stride] = acc1;

	// upper band :: ix=93
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v184,v190));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v185,v189),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v186,v188),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v187,acc1);
	__m256 v221 = data[221*stride];
	data[221*stride] = acc1;

	// lower band :: ix=94
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v184,v192));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v185,v191),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v186,v190),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v187,v189),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v188,acc1);
	data[94*stride] = acc1;

	// upper band :: ix=94
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v186,v192));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v187,v191),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v188,v190),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v189,acc1);
	__m256 v222 = data[222*stride];
	data[222*stride] = acc1;

	// lower band :: ix=95
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v186,v194));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v187,v193),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v188,v192),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v189,v191),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v190,acc1);
	data[95*stride] = acc1;

	// upper band :: ix=95
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v188,v194));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v189,v193),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v190,v192),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v191,acc1);
	__m256 v223 = data[223*stride];
	data[223*stride] = acc1;

	// lower band :: ix=96
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v188,v196));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v189,v195),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v190,v194),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v191,v193),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v192,acc1);
	data[96*stride] = acc1;

	// upper band :: ix=96
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v190,v196));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v191,v195),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v192,v194),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v193,acc1);
	__m256 v224 = data[224*stride];
	data[224*stride] = acc1;

	// lower band :: ix=97
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v190,v198));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v191,v197),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v192,v196),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v193,v195),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v194,acc1);
	data[97*stride] = acc1;

	// upper band :: ix=97
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v192,v198));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v193,v197),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v194,v196),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v195,acc1);
	__m256 v225 = data[225*stride];
	data[225*stride] = acc1;

	// lower band :: ix=98
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v192,v200));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v193,v199),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v194,v198),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v195,v197),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v196,acc1);
	data[98*stride] = acc1;

	// upper band :: ix=98
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v194,v200));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v195,v199),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v196,v198),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v197,acc1);
	__m256 v226 = data[226*stride];
	data[226*stride] = acc1;

	// lower band :: ix=99
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v194,v202));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v195,v201),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v196,v200),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v197,v199),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v198,acc1);
	data[99*stride] = acc1;

	// upper band :: ix=99
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v196,v202));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v197,v201),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v198,v200),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v199,acc1);
	__m256 v227 = data[227*stride];
	data[227*stride] = acc1;

	// lower band :: ix=100
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v196,v204));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v197,v203),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v198,v202),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v199,v201),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v200,acc1);
	data[100*stride] = acc1;

	// upper band :: ix=100
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v198,v204));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v199,v203),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v200,v202),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v201,acc1);
	__m256 v228 = data[228*stride];
	data[228*stride] = acc1;

	// lower band :: ix=101
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v198,v206));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v199,v205),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v200,v204),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v201,v203),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v202,acc1);
	data[101*stride] = acc1;

	// upper band :: ix=101
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v200,v206));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v201,v205),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v202,v204),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v203,acc1);
	__m256 v229 = data[229*stride];
	data[229*stride] = acc1;

	// lower band :: ix=102
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v200,v208));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v201,v207),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v202,v206),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v203,v205),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v204,acc1);
	data[102*stride] = acc1;

	// upper band :: ix=102
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v202,v208));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v203,v207),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v204,v206),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v205,acc1);
	__m256 v230 = data[230*stride];
	data[230*stride] = acc1;

	// lower band :: ix=103
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v202,v210));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v203,v209),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v204,v208),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v205,v207),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v206,acc1);
	data[103*stride] = acc1;

	// upper band :: ix=103
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v204,v210));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v205,v209),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v206,v208),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v207,acc1);
	__m256 v231 = data[231*stride];
	data[231*stride] = acc1;

	// lower band :: ix=104
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v204,v212));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v205,v211),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v206,v210),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v207,v209),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v208,acc1);
	data[104*stride] = acc1;

	// upper band :: ix=104
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v206,v212));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v207,v211),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v208,v210),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v209,acc1);
	__m256 v232 = data[232*stride];
	data[232*stride] = acc1;

	// lower band :: ix=105
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v206,v214));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v207,v213),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v208,v212),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v209,v211),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v210,acc1);
	data[105*stride] = acc1;

	// upper band :: ix=105
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v208,v214));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v209,v213),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v210,v212),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v211,acc1);
	__m256 v233 = data[233*stride];
	data[233*stride] = acc1;

	// lower band :: ix=106
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v208,v216));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v209,v215),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v210,v214),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v211,v213),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v212,acc1);
	data[106*stride] = acc1;

	// upper band :: ix=106
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v210,v216));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v211,v215),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v212,v214),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v213,acc1);
	__m256 v234 = data[234*stride];
	data[234*stride] = acc1;

	// lower band :: ix=107
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v210,v218));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v211,v217),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v212,v216),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v213,v215),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v214,acc1);
	data[107*stride] = acc1;

	// upper band :: ix=107
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v212,v218));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v213,v217),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v214,v216),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v215,acc1);
	__m256 v235 = data[235*stride];
	data[235*stride] = acc1;

	// lower band :: ix=108
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v212,v220));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v213,v219),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v214,v218),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v215,v217),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v216,acc1);
	data[108*stride] = acc1;

	// upper band :: ix=108
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v214,v220));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v215,v219),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v216,v218),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v217,acc1);
	__m256 v236 = data[236*stride];
	data[236*stride] = acc1;

	// lower band :: ix=109
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v214,v222));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v215,v221),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v216,v220),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v217,v219),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v218,acc1);
	data[109*stride] = acc1;

	// upper band :: ix=109
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v216,v222));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v217,v221),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v218,v220),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v219,acc1);
	__m256 v237 = data[237*stride];
	data[237*stride] = acc1;

	// lower band :: ix=110
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v216,v224));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v217,v223),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v218,v222),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v219,v221),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v220,acc1);
	data[110*stride] = acc1;

	// upper band :: ix=110
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v218,v224));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v219,v223),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v220,v222),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v221,acc1);
	__m256 v238 = data[238*stride];
	data[238*stride] = acc1;

	// lower band :: ix=111
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v218,v226));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v219,v225),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v220,v224),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v221,v223),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v222,acc1);
	data[111*stride] = acc1;

	// upper band :: ix=111
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v220,v226));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v221,v225),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v222,v224),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v223,acc1);
	__m256 v239 = data[239*stride];
	data[239*stride] = acc1;

	// lower band :: ix=112
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v220,v228));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v221,v227),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v222,v226),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v223,v225),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v224,acc1);
	data[112*stride] = acc1;

	// upper band :: ix=112
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v222,v228));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v223,v227),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v224,v226),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v225,acc1);
	__m256 v240 = data[240*stride];
	data[240*stride] = acc1;

	// lower band :: ix=113
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v222,v230));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v223,v229),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v224,v228),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v225,v227),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v226,acc1);
	data[113*stride] = acc1;

	// upper band :: ix=113
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v224,v230));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v225,v229),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v226,v228),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v227,acc1);
	__m256 v241 = data[241*stride];
	data[241*stride] = acc1;

	// lower band :: ix=114
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v224,v232));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v225,v231),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v226,v230),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v227,v229),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v228,acc1);
	data[114*stride] = acc1;

	// upper band :: ix=114
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v226,v232));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v227,v231),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v228,v230),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v229,acc1);
	__m256 v242 = data[242*stride];
	data[242*stride] = acc1;

	// lower band :: ix=115
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v226,v234));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v227,v233),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v228,v232),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v229,v231),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v230,acc1);
	data[115*stride] = acc1;

	// upper band :: ix=115
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v228,v234));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v229,v233),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v230,v232),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v231,acc1);
	__m256 v243 = data[243*stride];
	data[243*stride] = acc1;

	// lower band :: ix=116
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v228,v236));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v229,v235),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v230,v234),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v231,v233),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v232,acc1);
	data[116*stride] = acc1;

	// upper band :: ix=116
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v230,v236));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v231,v235),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v232,v234),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v233,acc1);
	__m256 v244 = data[244*stride];
	data[244*stride] = acc1;

	// lower band :: ix=117
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v230,v238));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v231,v237),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v232,v236),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v233,v235),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v234,acc1);
	data[117*stride] = acc1;

	// upper band :: ix=117
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v232,v238));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v233,v237),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v234,v236),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v235,acc1);
	__m256 v245 = data[245*stride];
	data[245*stride] = acc1;

	// lower band :: ix=118
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v232,v240));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v233,v239),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v234,v238),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v235,v237),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v236,acc1);
	data[118*stride] = acc1;

	// upper band :: ix=118
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v234,v240));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v235,v239),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v236,v238),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v237,acc1);
	__m256 v246 = data[246*stride];
	data[246*stride] = acc1;

	// lower band :: ix=119
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v234,v242));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v235,v241),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v236,v240),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v237,v239),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v238,acc1);
	data[119*stride] = acc1;

	// upper band :: ix=119
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v236,v242));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v237,v241),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v238,v240),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v239,acc1);
	__m256 v247 = data[247*stride];
	data[247*stride] = acc1;

	// lower band :: ix=120
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v236,v244));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v237,v243),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v238,v242),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v239,v241),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v240,acc1);
	data[120*stride] = acc1;

	// upper band :: ix=120
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v238,v244));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v239,v243),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v240,v242),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v241,acc1);
	__m256 v248 = data[248*stride];
	data[248*stride] = acc1;

	// lower band :: ix=121
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v238,v246));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v239,v245),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v240,v244),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v241,v243),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v242,acc1);
	data[121*stride] = acc1;

	// upper band :: ix=121
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v240,v246));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v241,v245),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v242,v244),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v243,acc1);
	__m256 v249 = data[249*stride];
	data[249*stride] = acc1;

	// lower band :: ix=122
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v240,v248));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v241,v247),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v242,v246),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v243,v245),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v244,acc1);
	data[122*stride] = acc1;

	// upper band :: ix=122
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v242,v248));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v243,v247),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v244,v246),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v245,acc1);
	__m256 v250 = data[250*stride];
	data[250*stride] = acc1;

	// lower band :: ix=123
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v242,v250));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v243,v249),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v244,v248),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v245,v247),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v246,acc1);
	data[123*stride] = acc1;

	// upper band :: ix=123
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v244,v250));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v245,v249),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v246,v248),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v247,acc1);
	__m256 v251 = data[251*stride];
	data[251*stride] = acc1;

	// lower band :: ix=124
	__m256 v252 = data[252*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v244,v252));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v245,v251),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v246,v250),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v247,v249),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v248,acc1);
	data[124*stride] = acc1;

	// upper band :: ix=124
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v246,v252));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v247,v251),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v248,v250),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v249,acc1);
	data[252*stride] = acc1;

	// lower band :: ix=125
	__m256 v254 = data[254*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v246,v254));
	__m256 v253 = data[253*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v247,v253),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v248,v252),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v249,v251),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v250,acc1);
	data[125*stride] = acc1;

	// upper band :: ix=125
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v248,v254));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v249,v253),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v250,v252),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v251,acc1);
	data[253*stride] = acc1;

	// lower band :: ix=126
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v248,v254));
	__m256 v255 = data[255*stride];
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v249,v255),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v250,v254),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v251,v253),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v252,acc1);
	data[126*stride] = acc1;

	// upper band :: ix=126
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v250,v254));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v251,v255),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v252,v254),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v253,acc1);
	data[254*stride] = acc1;

	// lower band :: ix=127
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v250,v252));
	acc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v251,v253),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v252,v254),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v253,v255),acc1);
	acc1 = _mm256_fmadd_ps(_mm_al0, v254,acc1);
	data[127*stride] = acc1;

	// upper band :: ix=127
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v252,v252));
	acc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v253,v253),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v254,v254),acc1);
	acc1 = _mm256_fmadd_ps(_mm_ah0,v255,acc1);
	data[255*stride] = acc1;
}

#else

static inline void _Ds79_AVX_2(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v0 = data[0*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v0,v0));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v1,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v0,v0)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v0));
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v0));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v0,v0)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v1));
	data[1*stride] = acc1;
}

static inline void _Ds79_AVX_4(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v2 = data[2*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v2));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v3,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v2,v2)));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v1)));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v0));
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v2));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v0,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v1));
	data[2*stride] = acc1;

	// lower band :: ix=1
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v0));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v1,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v0,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v2));
	data[1*stride] = acc1;

	// upper band :: ix=1
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v0));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v2,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v3));
	data[3*stride] = acc1;
}

static inline void _Ds79_AVX_8(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v4 = data[4*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v4));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v3,v3)));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v2,v2)));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v1)));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v0));
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v4));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v0,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v1));
	data[4*stride] = acc1;

	// lower band :: ix=1
	__m256 v6 = data[6*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v6));
	__m256 v5 = data[5*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v1,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v0,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v2));
	data[1*stride] = acc1;

	// upper band :: ix=1
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v6));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v2,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v3));
	data[5*stride] = acc1;

	// lower band :: ix=2
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v0,v6));
	__m256 v7 = data[7*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v1,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v2,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v3,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v4));
	data[2*stride] = acc1;

	// upper band :: ix=2
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v6));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v3,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v4,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v5));
	data[6*stride] = acc1;

	// lower band :: ix=3
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v4));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v3,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v4,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v5,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v6));
	data[3*stride] = acc1;

	// upper band :: ix=3
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v4,v4));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v5,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v6,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v7));
	data[7*stride] = acc1;
}

static inline void _Ds79_AVX_16(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v4 = data[4*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v4));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v3,v3)));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v2,v2)));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v1)));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v0));
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v4));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v0,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v1));
	__m256 v8 = data[8*stride];
	data[8*stride] = acc1;

	// lower band :: ix=1
	__m256 v6 = data[6*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v6));
	__m256 v5 = data[5*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v1,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v0,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v2));
	data[1*stride] = acc1;

	// upper band :: ix=1
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v6));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v2,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v3));
	__m256 v9 = data[9*stride];
	data[9*stride] = acc1;

	// lower band :: ix=2
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v0,v8));
	__m256 v7 = data[7*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v1,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v2,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v3,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v4));
	data[2*stride] = acc1;

	// upper band :: ix=2
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v8));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v3,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v4,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v5));
	__m256 v10 = data[10*stride];
	data[10*stride] = acc1;

	// lower band :: ix=3
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v10));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v3,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v4,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v5,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v6));
	data[3*stride] = acc1;

	// upper band :: ix=3
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v4,v10));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v5,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v6,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v7));
	__m256 v11 = data[11*stride];
	data[11*stride] = acc1;

	// lower band :: ix=4
	__m256 v12 = data[12*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v12));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v5,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v6,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v7,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v8));
	data[4*stride] = acc1;

	// upper band :: ix=4
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v6,v12));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v7,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v8,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v9));
	data[12*stride] = acc1;

	// lower band :: ix=5
	__m256 v14 = data[14*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v6,v14));
	__m256 v13 = data[13*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v7,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v8,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v9,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v10));
	data[5*stride] = acc1;

	// upper band :: ix=5
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v8,v14));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v9,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v10,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v11));
	data[13*stride] = acc1;

	// lower band :: ix=6
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v8,v14));
	__m256 v15 = data[15*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v9,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v10,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v11,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v12));
	data[6*stride] = acc1;

	// upper band :: ix=6
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v10,v14));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v11,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v12,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v13));
	data[14*stride] = acc1;

	// lower band :: ix=7
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v10,v12));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v11,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v12,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v13,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v14));
	data[7*stride] = acc1;

	// upper band :: ix=7
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v12,v12));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v13,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v14,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v15));
	data[15*stride] = acc1;
}

static inline void _Ds79_AVX_32(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v4 = data[4*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v4));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v3,v3)));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v2,v2)));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v1)));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v0));
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v4));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v0,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v1));
	__m256 v16 = data[16*stride];
	data[16*stride] = acc1;

	// lower band :: ix=1
	__m256 v6 = data[6*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v6));
	__m256 v5 = data[5*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v1,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v0,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v2));
	data[1*stride] = acc1;

	// upper band :: ix=1
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v6));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v2,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v3));
	__m256 v17 = data[17*stride];
	data[17*stride] = acc1;

	// lower band :: ix=2
	__m256 v8 = data[8*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v0,v8));
	__m256 v7 = data[7*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v1,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v2,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v3,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v4));
	data[2*stride] = acc1;

	// upper band :: ix=2
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v8));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v3,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v4,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v5));
	__m256 v18 = data[18*stride];
	data[18*stride] = acc1;

	// lower band :: ix=3
	__m256 v10 = data[10*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v10));
	__m256 v9 = data[9*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v3,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v4,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v5,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v6));
	data[3*stride] = acc1;

	// upper band :: ix=3
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v4,v10));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v5,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v6,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v7));
	__m256 v19 = data[19*stride];
	data[19*stride] = acc1;

	// lower band :: ix=4
	__m256 v12 = data[12*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v12));
	__m256 v11 = data[11*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v5,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v6,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v7,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v8));
	data[4*stride] = acc1;

	// upper band :: ix=4
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v6,v12));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v7,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v8,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v9));
	__m256 v20 = data[20*stride];
	data[20*stride] = acc1;

	// lower band :: ix=5
	__m256 v14 = data[14*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v6,v14));
	__m256 v13 = data[13*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v7,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v8,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v9,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v10));
	data[5*stride] = acc1;

	// upper band :: ix=5
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v8,v14));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v9,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v10,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v11));
	__m256 v21 = data[21*stride];
	data[21*stride] = acc1;

	// lower band :: ix=6
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v8,v16));
	__m256 v15 = data[15*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v9,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v10,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v11,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v12));
	data[6*stride] = acc1;

	// upper band :: ix=6
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v10,v16));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v11,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v12,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v13));
	__m256 v22 = data[22*stride];
	data[22*stride] = acc1;

	// lower band :: ix=7
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v10,v18));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v11,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v12,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v13,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v14));
	data[7*stride] = acc1;

	// upper band :: ix=7
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v12,v18));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v13,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v14,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v15));
	__m256 v23 = data[23*stride];
	data[23*stride] = acc1;

	// lower band :: ix=8
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v12,v20));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v13,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v14,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v15,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v16));
	data[8*stride] = acc1;

	// upper band :: ix=8
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v14,v20));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v15,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v16,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v17));
	__m256 v24 = data[24*stride];
	data[24*stride] = acc1;

	// lower band :: ix=9
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v14,v22));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v15,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v16,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v17,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v18));
	data[9*stride] = acc1;

	// upper band :: ix=9
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v16,v22));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v17,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v18,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v19));
	__m256 v25 = data[25*stride];
	data[25*stride] = acc1;

	// lower band :: ix=10
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v16,v24));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v17,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v18,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v19,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v20));
	data[10*stride] = acc1;

	// upper band :: ix=10
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v18,v24));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v19,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v20,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v21));
	__m256 v26 = data[26*stride];
	data[26*stride] = acc1;

	// lower band :: ix=11
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v18,v26));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v19,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v20,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v21,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v22));
	data[11*stride] = acc1;

	// upper band :: ix=11
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v20,v26));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v21,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v22,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v23));
	__m256 v27 = data[27*stride];
	data[27*stride] = acc1;

	// lower band :: ix=12
	__m256 v28 = data[28*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v20,v28));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v21,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v22,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v23,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v24));
	data[12*stride] = acc1;

	// upper band :: ix=12
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v22,v28));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v23,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v24,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v25));
	data[28*stride] = acc1;

	// lower band :: ix=13
	__m256 v30 = data[30*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v22,v30));
	__m256 v29 = data[29*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v23,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v24,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v25,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v26));
	data[13*stride] = acc1;

	// upper band :: ix=13
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v24,v30));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v25,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v26,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v27));
	data[29*stride] = acc1;

	// lower band :: ix=14
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v24,v30));
	__m256 v31 = data[31*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v25,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v26,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v27,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v28));
	data[14*stride] = acc1;

	// upper band :: ix=14
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v26,v30));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v27,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v28,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v29));
	data[30*stride] = acc1;

	// lower band :: ix=15
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v26,v28));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v27,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v28,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v29,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v30));
	data[15*stride] = acc1;

	// upper band :: ix=15
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v28,v28));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v29,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v30,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v31));
	data[31*stride] = acc1;
}

static void _Ds79_AVX_64(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v4 = data[4*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v4));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v3,v3)));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v2,v2)));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v1)));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v0));
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v4));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v0,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v1));
	__m256 v32 = data[32*stride];
	data[32*stride] = acc1;

	// lower band :: ix=1
	__m256 v6 = data[6*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v6));
	__m256 v5 = data[5*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v1,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v0,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v2));
	data[1*stride] = acc1;

	// upper band :: ix=1
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v6));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v2,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v3));
	__m256 v33 = data[33*stride];
	data[33*stride] = acc1;

	// lower band :: ix=2
	__m256 v8 = data[8*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v0,v8));
	__m256 v7 = data[7*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v1,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v2,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v3,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v4));
	data[2*stride] = acc1;

	// upper band :: ix=2
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v8));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v3,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v4,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v5));
	__m256 v34 = data[34*stride];
	data[34*stride] = acc1;

	// lower band :: ix=3
	__m256 v10 = data[10*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v10));
	__m256 v9 = data[9*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v3,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v4,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v5,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v6));
	data[3*stride] = acc1;

	// upper band :: ix=3
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v4,v10));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v5,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v6,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v7));
	__m256 v35 = data[35*stride];
	data[35*stride] = acc1;

	// lower band :: ix=4
	__m256 v12 = data[12*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v12));
	__m256 v11 = data[11*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v5,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v6,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v7,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v8));
	data[4*stride] = acc1;

	// upper band :: ix=4
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v6,v12));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v7,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v8,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v9));
	__m256 v36 = data[36*stride];
	data[36*stride] = acc1;

	// lower band :: ix=5
	__m256 v14 = data[14*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v6,v14));
	__m256 v13 = data[13*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v7,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v8,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v9,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v10));
	data[5*stride] = acc1;

	// upper band :: ix=5
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v8,v14));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v9,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v10,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v11));
	__m256 v37 = data[37*stride];
	data[37*stride] = acc1;

	// lower band :: ix=6
	__m256 v16 = data[16*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v8,v16));
	__m256 v15 = data[15*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v9,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v10,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v11,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v12));
	data[6*stride] = acc1;

	// upper band :: ix=6
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v10,v16));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v11,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v12,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v13));
	__m256 v38 = data[38*stride];
	data[38*stride] = acc1;

	// lower band :: ix=7
	__m256 v18 = data[18*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v10,v18));
	__m256 v17 = data[17*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v11,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v12,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v13,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v14));
	data[7*stride] = acc1;

	// upper band :: ix=7
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v12,v18));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v13,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v14,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v15));
	__m256 v39 = data[39*stride];
	data[39*stride] = acc1;

	// lower band :: ix=8
	__m256 v20 = data[20*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v12,v20));
	__m256 v19 = data[19*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v13,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v14,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v15,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v16));
	data[8*stride] = acc1;

	// upper band :: ix=8
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v14,v20));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v15,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v16,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v17));
	__m256 v40 = data[40*stride];
	data[40*stride] = acc1;

	// lower band :: ix=9
	__m256 v22 = data[22*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v14,v22));
	__m256 v21 = data[21*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v15,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v16,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v17,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v18));
	data[9*stride] = acc1;

	// upper band :: ix=9
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v16,v22));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v17,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v18,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v19));
	__m256 v41 = data[41*stride];
	data[41*stride] = acc1;

	// lower band :: ix=10
	__m256 v24 = data[24*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v16,v24));
	__m256 v23 = data[23*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v17,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v18,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v19,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v20));
	data[10*stride] = acc1;

	// upper band :: ix=10
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v18,v24));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v19,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v20,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v21));
	__m256 v42 = data[42*stride];
	data[42*stride] = acc1;

	// lower band :: ix=11
	__m256 v26 = data[26*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v18,v26));
	__m256 v25 = data[25*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v19,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v20,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v21,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v22));
	data[11*stride] = acc1;

	// upper band :: ix=11
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v20,v26));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v21,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v22,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v23));
	__m256 v43 = data[43*stride];
	data[43*stride] = acc1;

	// lower band :: ix=12
	__m256 v28 = data[28*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v20,v28));
	__m256 v27 = data[27*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v21,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v22,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v23,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v24));
	data[12*stride] = acc1;

	// upper band :: ix=12
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v22,v28));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v23,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v24,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v25));
	__m256 v44 = data[44*stride];
	data[44*stride] = acc1;

	// lower band :: ix=13
	__m256 v30 = data[30*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v22,v30));
	__m256 v29 = data[29*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v23,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v24,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v25,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v26));
	data[13*stride] = acc1;

	// upper band :: ix=13
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v24,v30));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v25,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v26,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v27));
	__m256 v45 = data[45*stride];
	data[45*stride] = acc1;

	// lower band :: ix=14
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v24,v32));
	__m256 v31 = data[31*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v25,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v26,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v27,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v28));
	data[14*stride] = acc1;

	// upper band :: ix=14
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v26,v32));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v27,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v28,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v29));
	__m256 v46 = data[46*stride];
	data[46*stride] = acc1;

	// lower band :: ix=15
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v26,v34));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v27,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v28,v32)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v29,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v30));
	data[15*stride] = acc1;

	// upper band :: ix=15
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v28,v34));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v29,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v30,v32)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v31));
	__m256 v47 = data[47*stride];
	data[47*stride] = acc1;

	// lower band :: ix=16
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v28,v36));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v29,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v30,v34)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v31,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v32));
	data[16*stride] = acc1;

	// upper band :: ix=16
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v30,v36));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v31,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v32,v34)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v33));
	__m256 v48 = data[48*stride];
	data[48*stride] = acc1;

	// lower band :: ix=17
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v30,v38));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v31,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v32,v36)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v33,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v34));
	data[17*stride] = acc1;

	// upper band :: ix=17
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v32,v38));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v33,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v34,v36)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v35));
	__m256 v49 = data[49*stride];
	data[49*stride] = acc1;

	// lower band :: ix=18
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v32,v40));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v33,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v34,v38)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v35,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v36));
	data[18*stride] = acc1;

	// upper band :: ix=18
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v34,v40));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v35,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v36,v38)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v37));
	__m256 v50 = data[50*stride];
	data[50*stride] = acc1;

	// lower band :: ix=19
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v34,v42));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v35,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v36,v40)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v37,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v38));
	data[19*stride] = acc1;

	// upper band :: ix=19
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v36,v42));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v37,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v38,v40)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v39));
	__m256 v51 = data[51*stride];
	data[51*stride] = acc1;

	// lower band :: ix=20
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v36,v44));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v37,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v38,v42)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v39,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v40));
	data[20*stride] = acc1;

	// upper band :: ix=20
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v38,v44));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v39,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v40,v42)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v41));
	__m256 v52 = data[52*stride];
	data[52*stride] = acc1;

	// lower band :: ix=21
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v38,v46));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v39,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v40,v44)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v41,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v42));
	data[21*stride] = acc1;

	// upper band :: ix=21
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v40,v46));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v41,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v42,v44)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v43));
	__m256 v53 = data[53*stride];
	data[53*stride] = acc1;

	// lower band :: ix=22
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v40,v48));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v41,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v42,v46)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v43,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v44));
	data[22*stride] = acc1;

	// upper band :: ix=22
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v42,v48));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v43,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v44,v46)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v45));
	__m256 v54 = data[54*stride];
	data[54*stride] = acc1;

	// lower band :: ix=23
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v42,v50));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v43,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v44,v48)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v45,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v46));
	data[23*stride] = acc1;

	// upper band :: ix=23
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v44,v50));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v45,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v46,v48)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v47));
	__m256 v55 = data[55*stride];
	data[55*stride] = acc1;

	// lower band :: ix=24
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v44,v52));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v45,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v46,v50)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v47,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v48));
	data[24*stride] = acc1;

	// upper band :: ix=24
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v46,v52));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v47,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v48,v50)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v49));
	__m256 v56 = data[56*stride];
	data[56*stride] = acc1;

	// lower band :: ix=25
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v46,v54));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v47,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v48,v52)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v49,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v50));
	data[25*stride] = acc1;

	// upper band :: ix=25
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v48,v54));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v49,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v50,v52)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v51));
	__m256 v57 = data[57*stride];
	data[57*stride] = acc1;

	// lower band :: ix=26
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v48,v56));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v49,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v50,v54)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v51,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v52));
	data[26*stride] = acc1;

	// upper band :: ix=26
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v50,v56));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v51,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v52,v54)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v53));
	__m256 v58 = data[58*stride];
	data[58*stride] = acc1;

	// lower band :: ix=27
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v50,v58));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v51,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v52,v56)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v53,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v54));
	data[27*stride] = acc1;

	// upper band :: ix=27
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v52,v58));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v53,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v54,v56)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v55));
	__m256 v59 = data[59*stride];
	data[59*stride] = acc1;

	// lower band :: ix=28
	__m256 v60 = data[60*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v52,v60));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v53,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v54,v58)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v55,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v56));
	data[28*stride] = acc1;

	// upper band :: ix=28
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v54,v60));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v55,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v56,v58)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v57));
	data[60*stride] = acc1;

	// lower band :: ix=29
	__m256 v62 = data[62*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v54,v62));
	__m256 v61 = data[61*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v55,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v56,v60)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v57,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v58));
	data[29*stride] = acc1;

	// upper band :: ix=29
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v56,v62));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v57,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v58,v60)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v59));
	data[61*stride] = acc1;

	// lower band :: ix=30
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v56,v62));
	__m256 v63 = data[63*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v57,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v58,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v59,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v60));
	data[30*stride] = acc1;

	// upper band :: ix=30
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v58,v62));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v59,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v60,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v61));
	data[62*stride] = acc1;

	// lower band :: ix=31
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v58,v60));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v59,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v60,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v61,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v62));
	data[31*stride] = acc1;

	// upper band :: ix=31
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v60,v60));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v61,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v62,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v63));
	data[63*stride] = acc1;
}

static void _Ds79_AVX_128(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v4 = data[4*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v4));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v3,v3)));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v2,v2)));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v1)));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v0));
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v4));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v0,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v1));
	__m256 v64 = data[64*stride];
	data[64*stride] = acc1;

	// lower band :: ix=1
	__m256 v6 = data[6*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v6));
	__m256 v5 = data[5*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v1,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v0,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v2));
	data[1*stride] = acc1;

	// upper band :: ix=1
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v6));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v2,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v3));
	__m256 v65 = data[65*stride];
	data[65*stride] = acc1;

	// lower band :: ix=2
	__m256 v8 = data[8*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v0,v8));
	__m256 v7 = data[7*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v1,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v2,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v3,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v4));
	data[2*stride] = acc1;

	// upper band :: ix=2
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v8));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v3,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v4,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v5));
	__m256 v66 = data[66*stride];
	data[66*stride] = acc1;

	// lower band :: ix=3
	__m256 v10 = data[10*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v10));
	__m256 v9 = data[9*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v3,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v4,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v5,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v6));
	data[3*stride] = acc1;

	// upper band :: ix=3
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v4,v10));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v5,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v6,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v7));
	__m256 v67 = data[67*stride];
	data[67*stride] = acc1;

	// lower band :: ix=4
	__m256 v12 = data[12*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v12));
	__m256 v11 = data[11*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v5,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v6,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v7,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v8));
	data[4*stride] = acc1;

	// upper band :: ix=4
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v6,v12));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v7,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v8,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v9));
	__m256 v68 = data[68*stride];
	data[68*stride] = acc1;

	// lower band :: ix=5
	__m256 v14 = data[14*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v6,v14));
	__m256 v13 = data[13*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v7,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v8,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v9,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v10));
	data[5*stride] = acc1;

	// upper band :: ix=5
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v8,v14));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v9,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v10,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v11));
	__m256 v69 = data[69*stride];
	data[69*stride] = acc1;

	// lower band :: ix=6
	__m256 v16 = data[16*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v8,v16));
	__m256 v15 = data[15*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v9,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v10,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v11,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v12));
	data[6*stride] = acc1;

	// upper band :: ix=6
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v10,v16));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v11,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v12,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v13));
	__m256 v70 = data[70*stride];
	data[70*stride] = acc1;

	// lower band :: ix=7
	__m256 v18 = data[18*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v10,v18));
	__m256 v17 = data[17*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v11,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v12,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v13,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v14));
	data[7*stride] = acc1;

	// upper band :: ix=7
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v12,v18));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v13,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v14,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v15));
	__m256 v71 = data[71*stride];
	data[71*stride] = acc1;

	// lower band :: ix=8
	__m256 v20 = data[20*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v12,v20));
	__m256 v19 = data[19*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v13,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v14,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v15,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v16));
	data[8*stride] = acc1;

	// upper band :: ix=8
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v14,v20));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v15,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v16,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v17));
	__m256 v72 = data[72*stride];
	data[72*stride] = acc1;

	// lower band :: ix=9
	__m256 v22 = data[22*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v14,v22));
	__m256 v21 = data[21*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v15,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v16,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v17,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v18));
	data[9*stride] = acc1;

	// upper band :: ix=9
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v16,v22));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v17,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v18,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v19));
	__m256 v73 = data[73*stride];
	data[73*stride] = acc1;

	// lower band :: ix=10
	__m256 v24 = data[24*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v16,v24));
	__m256 v23 = data[23*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v17,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v18,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v19,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v20));
	data[10*stride] = acc1;

	// upper band :: ix=10
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v18,v24));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v19,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v20,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v21));
	__m256 v74 = data[74*stride];
	data[74*stride] = acc1;

	// lower band :: ix=11
	__m256 v26 = data[26*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v18,v26));
	__m256 v25 = data[25*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v19,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v20,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v21,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v22));
	data[11*stride] = acc1;

	// upper band :: ix=11
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v20,v26));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v21,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v22,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v23));
	__m256 v75 = data[75*stride];
	data[75*stride] = acc1;

	// lower band :: ix=12
	__m256 v28 = data[28*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v20,v28));
	__m256 v27 = data[27*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v21,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v22,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v23,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v24));
	data[12*stride] = acc1;

	// upper band :: ix=12
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v22,v28));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v23,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v24,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v25));
	__m256 v76 = data[76*stride];
	data[76*stride] = acc1;

	// lower band :: ix=13
	__m256 v30 = data[30*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v22,v30));
	__m256 v29 = data[29*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v23,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v24,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v25,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v26));
	data[13*stride] = acc1;

	// upper band :: ix=13
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v24,v30));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v25,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v26,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v27));
	__m256 v77 = data[77*stride];
	data[77*stride] = acc1;

	// lower band :: ix=14
	__m256 v32 = data[32*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v24,v32));
	__m256 v31 = data[31*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v25,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v26,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v27,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v28));
	data[14*stride] = acc1;

	// upper band :: ix=14
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v26,v32));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v27,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v28,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v29));
	__m256 v78 = data[78*stride];
	data[78*stride] = acc1;

	// lower band :: ix=15
	__m256 v34 = data[34*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v26,v34));
	__m256 v33 = data[33*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v27,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v28,v32)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v29,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v30));
	data[15*stride] = acc1;

	// upper band :: ix=15
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v28,v34));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v29,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v30,v32)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v31));
	__m256 v79 = data[79*stride];
	data[79*stride] = acc1;

	// lower band :: ix=16
	__m256 v36 = data[36*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v28,v36));
	__m256 v35 = data[35*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v29,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v30,v34)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v31,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v32));
	data[16*stride] = acc1;

	// upper band :: ix=16
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v30,v36));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v31,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v32,v34)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v33));
	__m256 v80 = data[80*stride];
	data[80*stride] = acc1;

	// lower band :: ix=17
	__m256 v38 = data[38*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v30,v38));
	__m256 v37 = data[37*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v31,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v32,v36)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v33,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v34));
	data[17*stride] = acc1;

	// upper band :: ix=17
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v32,v38));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v33,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v34,v36)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v35));
	__m256 v81 = data[81*stride];
	data[81*stride] = acc1;

	// lower band :: ix=18
	__m256 v40 = data[40*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v32,v40));
	__m256 v39 = data[39*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v33,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v34,v38)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v35,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v36));
	data[18*stride] = acc1;

	// upper band :: ix=18
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v34,v40));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v35,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v36,v38)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v37));
	__m256 v82 = data[82*stride];
	data[82*stride] = acc1;

	// lower band :: ix=19
	__m256 v42 = data[42*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v34,v42));
	__m256 v41 = data[41*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v35,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v36,v40)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v37,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v38));
	data[19*stride] = acc1;

	// upper band :: ix=19
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v36,v42));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v37,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v38,v40)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v39));
	__m256 v83 = data[83*stride];
	data[83*stride] = acc1;

	// lower band :: ix=20
	__m256 v44 = data[44*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v36,v44));
	__m256 v43 = data[43*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v37,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v38,v42)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v39,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v40));
	data[20*stride] = acc1;

	// upper band :: ix=20
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v38,v44));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v39,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v40,v42)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v41));
	__m256 v84 = data[84*stride];
	data[84*stride] = acc1;

	// lower band :: ix=21
	__m256 v46 = data[46*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v38,v46));
	__m256 v45 = data[45*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v39,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v40,v44)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v41,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v42));
	data[21*stride] = acc1;

	// upper band :: ix=21
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v40,v46));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v41,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v42,v44)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v43));
	__m256 v85 = data[85*stride];
	data[85*stride] = acc1;

	// lower band :: ix=22
	__m256 v48 = data[48*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v40,v48));
	__m256 v47 = data[47*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v41,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v42,v46)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v43,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v44));
	data[22*stride] = acc1;

	// upper band :: ix=22
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v42,v48));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v43,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v44,v46)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v45));
	__m256 v86 = data[86*stride];
	data[86*stride] = acc1;

	// lower band :: ix=23
	__m256 v50 = data[50*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v42,v50));
	__m256 v49 = data[49*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v43,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v44,v48)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v45,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v46));
	data[23*stride] = acc1;

	// upper band :: ix=23
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v44,v50));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v45,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v46,v48)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v47));
	__m256 v87 = data[87*stride];
	data[87*stride] = acc1;

	// lower band :: ix=24
	__m256 v52 = data[52*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v44,v52));
	__m256 v51 = data[51*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v45,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v46,v50)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v47,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v48));
	data[24*stride] = acc1;

	// upper band :: ix=24
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v46,v52));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v47,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v48,v50)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v49));
	__m256 v88 = data[88*stride];
	data[88*stride] = acc1;

	// lower band :: ix=25
	__m256 v54 = data[54*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v46,v54));
	__m256 v53 = data[53*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v47,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v48,v52)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v49,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v50));
	data[25*stride] = acc1;

	// upper band :: ix=25
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v48,v54));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v49,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v50,v52)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v51));
	__m256 v89 = data[89*stride];
	data[89*stride] = acc1;

	// lower band :: ix=26
	__m256 v56 = data[56*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v48,v56));
	__m256 v55 = data[55*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v49,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v50,v54)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v51,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v52));
	data[26*stride] = acc1;

	// upper band :: ix=26
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v50,v56));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v51,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v52,v54)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v53));
	__m256 v90 = data[90*stride];
	data[90*stride] = acc1;

	// lower band :: ix=27
	__m256 v58 = data[58*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v50,v58));
	__m256 v57 = data[57*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v51,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v52,v56)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v53,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v54));
	data[27*stride] = acc1;

	// upper band :: ix=27
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v52,v58));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v53,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v54,v56)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v55));
	__m256 v91 = data[91*stride];
	data[91*stride] = acc1;

	// lower band :: ix=28
	__m256 v60 = data[60*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v52,v60));
	__m256 v59 = data[59*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v53,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v54,v58)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v55,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v56));
	data[28*stride] = acc1;

	// upper band :: ix=28
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v54,v60));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v55,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v56,v58)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v57));
	__m256 v92 = data[92*stride];
	data[92*stride] = acc1;

	// lower band :: ix=29
	__m256 v62 = data[62*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v54,v62));
	__m256 v61 = data[61*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v55,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v56,v60)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v57,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v58));
	data[29*stride] = acc1;

	// upper band :: ix=29
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v56,v62));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v57,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v58,v60)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v59));
	__m256 v93 = data[93*stride];
	data[93*stride] = acc1;

	// lower band :: ix=30
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v56,v64));
	__m256 v63 = data[63*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v57,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v58,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v59,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v60));
	data[30*stride] = acc1;

	// upper band :: ix=30
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v58,v64));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v59,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v60,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v61));
	__m256 v94 = data[94*stride];
	data[94*stride] = acc1;

	// lower band :: ix=31
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v58,v66));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v59,v65)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v60,v64)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v61,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v62));
	data[31*stride] = acc1;

	// upper band :: ix=31
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v60,v66));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v61,v65)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v62,v64)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v63));
	__m256 v95 = data[95*stride];
	data[95*stride] = acc1;

	// lower band :: ix=32
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v60,v68));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v61,v67)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v62,v66)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v63,v65)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v64));
	data[32*stride] = acc1;

	// upper band :: ix=32
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v62,v68));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v63,v67)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v64,v66)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v65));
	__m256 v96 = data[96*stride];
	data[96*stride] = acc1;

	// lower band :: ix=33
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v62,v70));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v63,v69)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v64,v68)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v65,v67)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v66));
	data[33*stride] = acc1;

	// upper band :: ix=33
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v64,v70));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v65,v69)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v66,v68)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v67));
	__m256 v97 = data[97*stride];
	data[97*stride] = acc1;

	// lower band :: ix=34
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v64,v72));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v65,v71)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v66,v70)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v67,v69)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v68));
	data[34*stride] = acc1;

	// upper band :: ix=34
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v66,v72));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v67,v71)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v68,v70)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v69));
	__m256 v98 = data[98*stride];
	data[98*stride] = acc1;

	// lower band :: ix=35
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v66,v74));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v67,v73)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v68,v72)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v69,v71)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v70));
	data[35*stride] = acc1;

	// upper band :: ix=35
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v68,v74));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v69,v73)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v70,v72)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v71));
	__m256 v99 = data[99*stride];
	data[99*stride] = acc1;

	// lower band :: ix=36
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v68,v76));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v69,v75)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v70,v74)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v71,v73)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v72));
	data[36*stride] = acc1;

	// upper band :: ix=36
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v70,v76));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v71,v75)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v72,v74)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v73));
	__m256 v100 = data[100*stride];
	data[100*stride] = acc1;

	// lower band :: ix=37
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v70,v78));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v71,v77)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v72,v76)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v73,v75)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v74));
	data[37*stride] = acc1;

	// upper band :: ix=37
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v72,v78));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v73,v77)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v74,v76)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v75));
	__m256 v101 = data[101*stride];
	data[101*stride] = acc1;

	// lower band :: ix=38
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v72,v80));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v73,v79)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v74,v78)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v75,v77)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v76));
	data[38*stride] = acc1;

	// upper band :: ix=38
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v74,v80));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v75,v79)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v76,v78)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v77));
	__m256 v102 = data[102*stride];
	data[102*stride] = acc1;

	// lower band :: ix=39
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v74,v82));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v75,v81)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v76,v80)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v77,v79)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v78));
	data[39*stride] = acc1;

	// upper band :: ix=39
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v76,v82));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v77,v81)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v78,v80)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v79));
	__m256 v103 = data[103*stride];
	data[103*stride] = acc1;

	// lower band :: ix=40
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v76,v84));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v77,v83)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v78,v82)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v79,v81)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v80));
	data[40*stride] = acc1;

	// upper band :: ix=40
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v78,v84));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v79,v83)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v80,v82)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v81));
	__m256 v104 = data[104*stride];
	data[104*stride] = acc1;

	// lower band :: ix=41
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v78,v86));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v79,v85)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v80,v84)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v81,v83)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v82));
	data[41*stride] = acc1;

	// upper band :: ix=41
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v80,v86));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v81,v85)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v82,v84)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v83));
	__m256 v105 = data[105*stride];
	data[105*stride] = acc1;

	// lower band :: ix=42
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v80,v88));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v81,v87)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v82,v86)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v83,v85)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v84));
	data[42*stride] = acc1;

	// upper band :: ix=42
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v82,v88));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v83,v87)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v84,v86)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v85));
	__m256 v106 = data[106*stride];
	data[106*stride] = acc1;

	// lower band :: ix=43
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v82,v90));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v83,v89)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v84,v88)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v85,v87)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v86));
	data[43*stride] = acc1;

	// upper band :: ix=43
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v84,v90));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v85,v89)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v86,v88)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v87));
	__m256 v107 = data[107*stride];
	data[107*stride] = acc1;

	// lower band :: ix=44
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v84,v92));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v85,v91)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v86,v90)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v87,v89)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v88));
	data[44*stride] = acc1;

	// upper band :: ix=44
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v86,v92));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v87,v91)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v88,v90)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v89));
	__m256 v108 = data[108*stride];
	data[108*stride] = acc1;

	// lower band :: ix=45
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v86,v94));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v87,v93)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v88,v92)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v89,v91)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v90));
	data[45*stride] = acc1;

	// upper band :: ix=45
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v88,v94));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v89,v93)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v90,v92)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v91));
	__m256 v109 = data[109*stride];
	data[109*stride] = acc1;

	// lower band :: ix=46
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v88,v96));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v89,v95)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v90,v94)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v91,v93)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v92));
	data[46*stride] = acc1;

	// upper band :: ix=46
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v90,v96));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v91,v95)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v92,v94)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v93));
	__m256 v110 = data[110*stride];
	data[110*stride] = acc1;

	// lower band :: ix=47
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v90,v98));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v91,v97)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v92,v96)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v93,v95)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v94));
	data[47*stride] = acc1;

	// upper band :: ix=47
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v92,v98));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v93,v97)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v94,v96)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v95));
	__m256 v111 = data[111*stride];
	data[111*stride] = acc1;

	// lower band :: ix=48
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v92,v100));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v93,v99)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v94,v98)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v95,v97)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v96));
	data[48*stride] = acc1;

	// upper band :: ix=48
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v94,v100));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v95,v99)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v96,v98)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v97));
	__m256 v112 = data[112*stride];
	data[112*stride] = acc1;

	// lower band :: ix=49
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v94,v102));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v95,v101)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v96,v100)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v97,v99)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v98));
	data[49*stride] = acc1;

	// upper band :: ix=49
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v96,v102));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v97,v101)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v98,v100)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v99));
	__m256 v113 = data[113*stride];
	data[113*stride] = acc1;

	// lower band :: ix=50
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v96,v104));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v97,v103)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v98,v102)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v99,v101)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v100));
	data[50*stride] = acc1;

	// upper band :: ix=50
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v98,v104));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v99,v103)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v100,v102)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v101));
	__m256 v114 = data[114*stride];
	data[114*stride] = acc1;

	// lower band :: ix=51
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v98,v106));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v99,v105)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v100,v104)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v101,v103)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v102));
	data[51*stride] = acc1;

	// upper band :: ix=51
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v100,v106));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v101,v105)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v102,v104)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v103));
	__m256 v115 = data[115*stride];
	data[115*stride] = acc1;

	// lower band :: ix=52
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v100,v108));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v101,v107)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v102,v106)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v103,v105)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v104));
	data[52*stride] = acc1;

	// upper band :: ix=52
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v102,v108));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v103,v107)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v104,v106)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v105));
	__m256 v116 = data[116*stride];
	data[116*stride] = acc1;

	// lower band :: ix=53
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v102,v110));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v103,v109)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v104,v108)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v105,v107)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v106));
	data[53*stride] = acc1;

	// upper band :: ix=53
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v104,v110));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v105,v109)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v106,v108)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v107));
	__m256 v117 = data[117*stride];
	data[117*stride] = acc1;

	// lower band :: ix=54
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v104,v112));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v105,v111)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v106,v110)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v107,v109)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v108));
	data[54*stride] = acc1;

	// upper band :: ix=54
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v106,v112));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v107,v111)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v108,v110)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v109));
	__m256 v118 = data[118*stride];
	data[118*stride] = acc1;

	// lower band :: ix=55
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v106,v114));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v107,v113)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v108,v112)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v109,v111)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v110));
	data[55*stride] = acc1;

	// upper band :: ix=55
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v108,v114));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v109,v113)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v110,v112)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v111));
	__m256 v119 = data[119*stride];
	data[119*stride] = acc1;

	// lower band :: ix=56
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v108,v116));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v109,v115)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v110,v114)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v111,v113)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v112));
	data[56*stride] = acc1;

	// upper band :: ix=56
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v110,v116));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v111,v115)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v112,v114)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v113));
	__m256 v120 = data[120*stride];
	data[120*stride] = acc1;

	// lower band :: ix=57
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v110,v118));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v111,v117)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v112,v116)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v113,v115)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v114));
	data[57*stride] = acc1;

	// upper band :: ix=57
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v112,v118));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v113,v117)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v114,v116)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v115));
	__m256 v121 = data[121*stride];
	data[121*stride] = acc1;

	// lower band :: ix=58
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v112,v120));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v113,v119)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v114,v118)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v115,v117)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v116));
	data[58*stride] = acc1;

	// upper band :: ix=58
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v114,v120));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v115,v119)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v116,v118)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v117));
	__m256 v122 = data[122*stride];
	data[122*stride] = acc1;

	// lower band :: ix=59
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v114,v122));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v115,v121)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v116,v120)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v117,v119)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v118));
	data[59*stride] = acc1;

	// upper band :: ix=59
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v116,v122));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v117,v121)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v118,v120)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v119));
	__m256 v123 = data[123*stride];
	data[123*stride] = acc1;

	// lower band :: ix=60
	__m256 v124 = data[124*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v116,v124));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v117,v123)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v118,v122)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v119,v121)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v120));
	data[60*stride] = acc1;

	// upper band :: ix=60
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v118,v124));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v119,v123)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v120,v122)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v121));
	data[124*stride] = acc1;

	// lower band :: ix=61
	__m256 v126 = data[126*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v118,v126));
	__m256 v125 = data[125*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v119,v125)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v120,v124)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v121,v123)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v122));
	data[61*stride] = acc1;

	// upper band :: ix=61
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v120,v126));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v121,v125)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v122,v124)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v123));
	data[125*stride] = acc1;

	// lower band :: ix=62
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v120,v126));
	__m256 v127 = data[127*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v121,v127)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v122,v126)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v123,v125)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v124));
	data[62*stride] = acc1;

	// upper band :: ix=62
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v122,v126));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v123,v127)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v124,v126)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v125));
	data[126*stride] = acc1;

	// lower band :: ix=63
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v122,v124));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v123,v125)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v124,v126)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v125,v127)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v126));
	data[63*stride] = acc1;

	// upper band :: ix=63
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v124,v124));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v125,v125)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v126,v126)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v127));
	data[127*stride] = acc1;
}

static void _Ds79_AVX_256(__m256* data, int stride)
{
	__m256 acc1;

	// lower band :: ix=0
	__m256 v4 = data[4*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v4));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v3,v3)));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v2,v2)));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v1)));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v0));
	data[0*stride] = acc1;

	// upper band :: ix=0
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v4));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v0,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v1));
	__m256 v128 = data[128*stride];
	data[128*stride] = acc1;

	// lower band :: ix=1
	__m256 v6 = data[6*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v6));
	__m256 v5 = data[5*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v1,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v0,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v2));
	data[1*stride] = acc1;

	// upper band :: ix=1
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v0,v6));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v1,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v2,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v3));
	__m256 v129 = data[129*stride];
	data[129*stride] = acc1;

	// lower band :: ix=2
	__m256 v8 = data[8*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v0,v8));
	__m256 v7 = data[7*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v1,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v2,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v3,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v4));
	data[2*stride] = acc1;

	// upper band :: ix=2
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v2,v8));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v3,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v4,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v5));
	__m256 v130 = data[130*stride];
	data[130*stride] = acc1;

	// lower band :: ix=3
	__m256 v10 = data[10*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v2,v10));
	__m256 v9 = data[9*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v3,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v4,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v5,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v6));
	data[3*stride] = acc1;

	// upper band :: ix=3
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v4,v10));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v5,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v6,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v7));
	__m256 v131 = data[131*stride];
	data[131*stride] = acc1;

	// lower band :: ix=4
	__m256 v12 = data[12*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v4,v12));
	__m256 v11 = data[11*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v5,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v6,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v7,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v8));
	data[4*stride] = acc1;

	// upper band :: ix=4
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v6,v12));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v7,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v8,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v9));
	__m256 v132 = data[132*stride];
	data[132*stride] = acc1;

	// lower band :: ix=5
	__m256 v14 = data[14*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v6,v14));
	__m256 v13 = data[13*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v7,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v8,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v9,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v10));
	data[5*stride] = acc1;

	// upper band :: ix=5
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v8,v14));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v9,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v10,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v11));
	__m256 v133 = data[133*stride];
	data[133*stride] = acc1;

	// lower band :: ix=6
	__m256 v16 = data[16*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v8,v16));
	__m256 v15 = data[15*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v9,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v10,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v11,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v12));
	data[6*stride] = acc1;

	// upper band :: ix=6
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v10,v16));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v11,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v12,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v13));
	__m256 v134 = data[134*stride];
	data[134*stride] = acc1;

	// lower band :: ix=7
	__m256 v18 = data[18*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v10,v18));
	__m256 v17 = data[17*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v11,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v12,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v13,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v14));
	data[7*stride] = acc1;

	// upper band :: ix=7
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v12,v18));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v13,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v14,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v15));
	__m256 v135 = data[135*stride];
	data[135*stride] = acc1;

	// lower band :: ix=8
	__m256 v20 = data[20*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v12,v20));
	__m256 v19 = data[19*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v13,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v14,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v15,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v16));
	data[8*stride] = acc1;

	// upper band :: ix=8
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v14,v20));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v15,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v16,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v17));
	__m256 v136 = data[136*stride];
	data[136*stride] = acc1;

	// lower band :: ix=9
	__m256 v22 = data[22*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v14,v22));
	__m256 v21 = data[21*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v15,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v16,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v17,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v18));
	data[9*stride] = acc1;

	// upper band :: ix=9
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v16,v22));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v17,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v18,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v19));
	__m256 v137 = data[137*stride];
	data[137*stride] = acc1;

	// lower band :: ix=10
	__m256 v24 = data[24*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v16,v24));
	__m256 v23 = data[23*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v17,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v18,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v19,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v20));
	data[10*stride] = acc1;

	// upper band :: ix=10
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v18,v24));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v19,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v20,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v21));
	__m256 v138 = data[138*stride];
	data[138*stride] = acc1;

	// lower band :: ix=11
	__m256 v26 = data[26*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v18,v26));
	__m256 v25 = data[25*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v19,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v20,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v21,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v22));
	data[11*stride] = acc1;

	// upper band :: ix=11
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v20,v26));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v21,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v22,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v23));
	__m256 v139 = data[139*stride];
	data[139*stride] = acc1;

	// lower band :: ix=12
	__m256 v28 = data[28*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v20,v28));
	__m256 v27 = data[27*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v21,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v22,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v23,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v24));
	data[12*stride] = acc1;

	// upper band :: ix=12
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v22,v28));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v23,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v24,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v25));
	__m256 v140 = data[140*stride];
	data[140*stride] = acc1;

	// lower band :: ix=13
	__m256 v30 = data[30*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v22,v30));
	__m256 v29 = data[29*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v23,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v24,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v25,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v26));
	data[13*stride] = acc1;

	// upper band :: ix=13
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v24,v30));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v25,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v26,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v27));
	__m256 v141 = data[141*stride];
	data[141*stride] = acc1;

	// lower band :: ix=14
	__m256 v32 = data[32*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v24,v32));
	__m256 v31 = data[31*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v25,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v26,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v27,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v28));
	data[14*stride] = acc1;

	// upper band :: ix=14
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v26,v32));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v27,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v28,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v29));
	__m256 v142 = data[142*stride];
	data[142*stride] = acc1;

	// lower band :: ix=15
	__m256 v34 = data[34*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v26,v34));
	__m256 v33 = data[33*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v27,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v28,v32)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v29,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v30));
	data[15*stride] = acc1;

	// upper band :: ix=15
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v28,v34));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v29,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v30,v32)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v31));
	__m256 v143 = data[143*stride];
	data[143*stride] = acc1;

	// lower band :: ix=16
	__m256 v36 = data[36*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v28,v36));
	__m256 v35 = data[35*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v29,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v30,v34)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v31,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v32));
	data[16*stride] = acc1;

	// upper band :: ix=16
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v30,v36));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v31,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v32,v34)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v33));
	__m256 v144 = data[144*stride];
	data[144*stride] = acc1;

	// lower band :: ix=17
	__m256 v38 = data[38*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v30,v38));
	__m256 v37 = data[37*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v31,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v32,v36)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v33,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v34));
	data[17*stride] = acc1;

	// upper band :: ix=17
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v32,v38));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v33,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v34,v36)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v35));
	__m256 v145 = data[145*stride];
	data[145*stride] = acc1;

	// lower band :: ix=18
	__m256 v40 = data[40*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v32,v40));
	__m256 v39 = data[39*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v33,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v34,v38)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v35,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v36));
	data[18*stride] = acc1;

	// upper band :: ix=18
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v34,v40));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v35,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v36,v38)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v37));
	__m256 v146 = data[146*stride];
	data[146*stride] = acc1;

	// lower band :: ix=19
	__m256 v42 = data[42*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v34,v42));
	__m256 v41 = data[41*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v35,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v36,v40)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v37,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v38));
	data[19*stride] = acc1;

	// upper band :: ix=19
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v36,v42));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v37,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v38,v40)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v39));
	__m256 v147 = data[147*stride];
	data[147*stride] = acc1;

	// lower band :: ix=20
	__m256 v44 = data[44*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v36,v44));
	__m256 v43 = data[43*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v37,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v38,v42)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v39,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v40));
	data[20*stride] = acc1;

	// upper band :: ix=20
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v38,v44));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v39,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v40,v42)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v41));
	__m256 v148 = data[148*stride];
	data[148*stride] = acc1;

	// lower band :: ix=21
	__m256 v46 = data[46*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v38,v46));
	__m256 v45 = data[45*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v39,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v40,v44)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v41,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v42));
	data[21*stride] = acc1;

	// upper band :: ix=21
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v40,v46));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v41,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v42,v44)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v43));
	__m256 v149 = data[149*stride];
	data[149*stride] = acc1;

	// lower band :: ix=22
	__m256 v48 = data[48*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v40,v48));
	__m256 v47 = data[47*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v41,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v42,v46)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v43,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v44));
	data[22*stride] = acc1;

	// upper band :: ix=22
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v42,v48));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v43,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v44,v46)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v45));
	__m256 v150 = data[150*stride];
	data[150*stride] = acc1;

	// lower band :: ix=23
	__m256 v50 = data[50*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v42,v50));
	__m256 v49 = data[49*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v43,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v44,v48)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v45,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v46));
	data[23*stride] = acc1;

	// upper band :: ix=23
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v44,v50));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v45,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v46,v48)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v47));
	__m256 v151 = data[151*stride];
	data[151*stride] = acc1;

	// lower band :: ix=24
	__m256 v52 = data[52*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v44,v52));
	__m256 v51 = data[51*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v45,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v46,v50)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v47,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v48));
	data[24*stride] = acc1;

	// upper band :: ix=24
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v46,v52));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v47,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v48,v50)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v49));
	__m256 v152 = data[152*stride];
	data[152*stride] = acc1;

	// lower band :: ix=25
	__m256 v54 = data[54*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v46,v54));
	__m256 v53 = data[53*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v47,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v48,v52)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v49,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v50));
	data[25*stride] = acc1;

	// upper band :: ix=25
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v48,v54));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v49,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v50,v52)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v51));
	__m256 v153 = data[153*stride];
	data[153*stride] = acc1;

	// lower band :: ix=26
	__m256 v56 = data[56*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v48,v56));
	__m256 v55 = data[55*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v49,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v50,v54)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v51,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v52));
	data[26*stride] = acc1;

	// upper band :: ix=26
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v50,v56));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v51,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v52,v54)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v53));
	__m256 v154 = data[154*stride];
	data[154*stride] = acc1;

	// lower band :: ix=27
	__m256 v58 = data[58*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v50,v58));
	__m256 v57 = data[57*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v51,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v52,v56)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v53,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v54));
	data[27*stride] = acc1;

	// upper band :: ix=27
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v52,v58));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v53,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v54,v56)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v55));
	__m256 v155 = data[155*stride];
	data[155*stride] = acc1;

	// lower band :: ix=28
	__m256 v60 = data[60*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v52,v60));
	__m256 v59 = data[59*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v53,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v54,v58)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v55,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v56));
	data[28*stride] = acc1;

	// upper band :: ix=28
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v54,v60));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v55,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v56,v58)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v57));
	__m256 v156 = data[156*stride];
	data[156*stride] = acc1;

	// lower band :: ix=29
	__m256 v62 = data[62*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v54,v62));
	__m256 v61 = data[61*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v55,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v56,v60)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v57,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v58));
	data[29*stride] = acc1;

	// upper band :: ix=29
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v56,v62));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v57,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v58,v60)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v59));
	__m256 v157 = data[157*stride];
	data[157*stride] = acc1;

	// lower band :: ix=30
	__m256 v64 = data[64*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v56,v64));
	__m256 v63 = data[63*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v57,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v58,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v59,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v60));
	data[30*stride] = acc1;

	// upper band :: ix=30
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v58,v64));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v59,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v60,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v61));
	__m256 v158 = data[158*stride];
	data[158*stride] = acc1;

	// lower band :: ix=31
	__m256 v66 = data[66*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v58,v66));
	__m256 v65 = data[65*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v59,v65)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v60,v64)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v61,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v62));
	data[31*stride] = acc1;

	// upper band :: ix=31
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v60,v66));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v61,v65)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v62,v64)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v63));
	__m256 v159 = data[159*stride];
	data[159*stride] = acc1;

	// lower band :: ix=32
	__m256 v68 = data[68*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v60,v68));
	__m256 v67 = data[67*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v61,v67)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v62,v66)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v63,v65)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v64));
	data[32*stride] = acc1;

	// upper band :: ix=32
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v62,v68));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v63,v67)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v64,v66)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v65));
	__m256 v160 = data[160*stride];
	data[160*stride] = acc1;

	// lower band :: ix=33
	__m256 v70 = data[70*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v62,v70));
	__m256 v69 = data[69*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v63,v69)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v64,v68)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v65,v67)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v66));
	data[33*stride] = acc1;

	// upper band :: ix=33
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v64,v70));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v65,v69)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v66,v68)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v67));
	__m256 v161 = data[161*stride];
	data[161*stride] = acc1;

	// lower band :: ix=34
	__m256 v72 = data[72*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v64,v72));
	__m256 v71 = data[71*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v65,v71)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v66,v70)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v67,v69)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v68));
	data[34*stride] = acc1;

	// upper band :: ix=34
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v66,v72));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v67,v71)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v68,v70)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v69));
	__m256 v162 = data[162*stride];
	data[162*stride] = acc1;

	// lower band :: ix=35
	__m256 v74 = data[74*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v66,v74));
	__m256 v73 = data[73*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v67,v73)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v68,v72)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v69,v71)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v70));
	data[35*stride] = acc1;

	// upper band :: ix=35
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v68,v74));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v69,v73)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v70,v72)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v71));
	__m256 v163 = data[163*stride];
	data[163*stride] = acc1;

	// lower band :: ix=36
	__m256 v76 = data[76*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v68,v76));
	__m256 v75 = data[75*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v69,v75)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v70,v74)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v71,v73)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v72));
	data[36*stride] = acc1;

	// upper band :: ix=36
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v70,v76));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v71,v75)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v72,v74)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v73));
	__m256 v164 = data[164*stride];
	data[164*stride] = acc1;

	// lower band :: ix=37
	__m256 v78 = data[78*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v70,v78));
	__m256 v77 = data[77*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v71,v77)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v72,v76)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v73,v75)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v74));
	data[37*stride] = acc1;

	// upper band :: ix=37
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v72,v78));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v73,v77)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v74,v76)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v75));
	__m256 v165 = data[165*stride];
	data[165*stride] = acc1;

	// lower band :: ix=38
	__m256 v80 = data[80*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v72,v80));
	__m256 v79 = data[79*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v73,v79)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v74,v78)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v75,v77)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v76));
	data[38*stride] = acc1;

	// upper band :: ix=38
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v74,v80));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v75,v79)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v76,v78)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v77));
	__m256 v166 = data[166*stride];
	data[166*stride] = acc1;

	// lower band :: ix=39
	__m256 v82 = data[82*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v74,v82));
	__m256 v81 = data[81*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v75,v81)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v76,v80)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v77,v79)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v78));
	data[39*stride] = acc1;

	// upper band :: ix=39
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v76,v82));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v77,v81)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v78,v80)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v79));
	__m256 v167 = data[167*stride];
	data[167*stride] = acc1;

	// lower band :: ix=40
	__m256 v84 = data[84*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v76,v84));
	__m256 v83 = data[83*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v77,v83)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v78,v82)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v79,v81)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v80));
	data[40*stride] = acc1;

	// upper band :: ix=40
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v78,v84));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v79,v83)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v80,v82)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v81));
	__m256 v168 = data[168*stride];
	data[168*stride] = acc1;

	// lower band :: ix=41
	__m256 v86 = data[86*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v78,v86));
	__m256 v85 = data[85*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v79,v85)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v80,v84)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v81,v83)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v82));
	data[41*stride] = acc1;

	// upper band :: ix=41
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v80,v86));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v81,v85)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v82,v84)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v83));
	__m256 v169 = data[169*stride];
	data[169*stride] = acc1;

	// lower band :: ix=42
	__m256 v88 = data[88*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v80,v88));
	__m256 v87 = data[87*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v81,v87)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v82,v86)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v83,v85)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v84));
	data[42*stride] = acc1;

	// upper band :: ix=42
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v82,v88));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v83,v87)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v84,v86)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v85));
	__m256 v170 = data[170*stride];
	data[170*stride] = acc1;

	// lower band :: ix=43
	__m256 v90 = data[90*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v82,v90));
	__m256 v89 = data[89*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v83,v89)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v84,v88)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v85,v87)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v86));
	data[43*stride] = acc1;

	// upper band :: ix=43
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v84,v90));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v85,v89)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v86,v88)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v87));
	__m256 v171 = data[171*stride];
	data[171*stride] = acc1;

	// lower band :: ix=44
	__m256 v92 = data[92*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v84,v92));
	__m256 v91 = data[91*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v85,v91)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v86,v90)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v87,v89)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v88));
	data[44*stride] = acc1;

	// upper band :: ix=44
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v86,v92));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v87,v91)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v88,v90)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v89));
	__m256 v172 = data[172*stride];
	data[172*stride] = acc1;

	// lower band :: ix=45
	__m256 v94 = data[94*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v86,v94));
	__m256 v93 = data[93*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v87,v93)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v88,v92)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v89,v91)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v90));
	data[45*stride] = acc1;

	// upper band :: ix=45
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v88,v94));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v89,v93)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v90,v92)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v91));
	__m256 v173 = data[173*stride];
	data[173*stride] = acc1;

	// lower band :: ix=46
	__m256 v96 = data[96*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v88,v96));
	__m256 v95 = data[95*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v89,v95)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v90,v94)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v91,v93)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v92));
	data[46*stride] = acc1;

	// upper band :: ix=46
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v90,v96));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v91,v95)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v92,v94)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v93));
	__m256 v174 = data[174*stride];
	data[174*stride] = acc1;

	// lower band :: ix=47
	__m256 v98 = data[98*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v90,v98));
	__m256 v97 = data[97*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v91,v97)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v92,v96)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v93,v95)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v94));
	data[47*stride] = acc1;

	// upper band :: ix=47
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v92,v98));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v93,v97)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v94,v96)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v95));
	__m256 v175 = data[175*stride];
	data[175*stride] = acc1;

	// lower band :: ix=48
	__m256 v100 = data[100*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v92,v100));
	__m256 v99 = data[99*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v93,v99)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v94,v98)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v95,v97)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v96));
	data[48*stride] = acc1;

	// upper band :: ix=48
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v94,v100));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v95,v99)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v96,v98)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v97));
	__m256 v176 = data[176*stride];
	data[176*stride] = acc1;

	// lower band :: ix=49
	__m256 v102 = data[102*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v94,v102));
	__m256 v101 = data[101*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v95,v101)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v96,v100)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v97,v99)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v98));
	data[49*stride] = acc1;

	// upper band :: ix=49
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v96,v102));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v97,v101)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v98,v100)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v99));
	__m256 v177 = data[177*stride];
	data[177*stride] = acc1;

	// lower band :: ix=50
	__m256 v104 = data[104*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v96,v104));
	__m256 v103 = data[103*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v97,v103)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v98,v102)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v99,v101)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v100));
	data[50*stride] = acc1;

	// upper band :: ix=50
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v98,v104));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v99,v103)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v100,v102)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v101));
	__m256 v178 = data[178*stride];
	data[178*stride] = acc1;

	// lower band :: ix=51
	__m256 v106 = data[106*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v98,v106));
	__m256 v105 = data[105*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v99,v105)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v100,v104)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v101,v103)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v102));
	data[51*stride] = acc1;

	// upper band :: ix=51
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v100,v106));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v101,v105)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v102,v104)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v103));
	__m256 v179 = data[179*stride];
	data[179*stride] = acc1;

	// lower band :: ix=52
	__m256 v108 = data[108*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v100,v108));
	__m256 v107 = data[107*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v101,v107)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v102,v106)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v103,v105)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v104));
	data[52*stride] = acc1;

	// upper band :: ix=52
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v102,v108));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v103,v107)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v104,v106)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v105));
	__m256 v180 = data[180*stride];
	data[180*stride] = acc1;

	// lower band :: ix=53
	__m256 v110 = data[110*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v102,v110));
	__m256 v109 = data[109*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v103,v109)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v104,v108)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v105,v107)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v106));
	data[53*stride] = acc1;

	// upper band :: ix=53
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v104,v110));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v105,v109)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v106,v108)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v107));
	__m256 v181 = data[181*stride];
	data[181*stride] = acc1;

	// lower band :: ix=54
	__m256 v112 = data[112*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v104,v112));
	__m256 v111 = data[111*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v105,v111)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v106,v110)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v107,v109)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v108));
	data[54*stride] = acc1;

	// upper band :: ix=54
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v106,v112));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v107,v111)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v108,v110)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v109));
	__m256 v182 = data[182*stride];
	data[182*stride] = acc1;

	// lower band :: ix=55
	__m256 v114 = data[114*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v106,v114));
	__m256 v113 = data[113*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v107,v113)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v108,v112)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v109,v111)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v110));
	data[55*stride] = acc1;

	// upper band :: ix=55
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v108,v114));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v109,v113)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v110,v112)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v111));
	__m256 v183 = data[183*stride];
	data[183*stride] = acc1;

	// lower band :: ix=56
	__m256 v116 = data[116*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v108,v116));
	__m256 v115 = data[115*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v109,v115)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v110,v114)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v111,v113)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v112));
	data[56*stride] = acc1;

	// upper band :: ix=56
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v110,v116));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v111,v115)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v112,v114)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v113));
	__m256 v184 = data[184*stride];
	data[184*stride] = acc1;

	// lower band :: ix=57
	__m256 v118 = data[118*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v110,v118));
	__m256 v117 = data[117*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v111,v117)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v112,v116)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v113,v115)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v114));
	data[57*stride] = acc1;

	// upper band :: ix=57
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v112,v118));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v113,v117)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v114,v116)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v115));
	__m256 v185 = data[185*stride];
	data[185*stride] = acc1;

	// lower band :: ix=58
	__m256 v120 = data[120*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v112,v120));
	__m256 v119 = data[119*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v113,v119)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v114,v118)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v115,v117)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v116));
	data[58*stride] = acc1;

	// upper band :: ix=58
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v114,v120));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v115,v119)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v116,v118)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v117));
	__m256 v186 = data[186*stride];
	data[186*stride] = acc1;

	// lower band :: ix=59
	__m256 v122 = data[122*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v114,v122));
	__m256 v121 = data[121*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v115,v121)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v116,v120)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v117,v119)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v118));
	data[59*stride] = acc1;

	// upper band :: ix=59
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v116,v122));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v117,v121)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v118,v120)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v119));
	__m256 v187 = data[187*stride];
	data[187*stride] = acc1;

	// lower band :: ix=60
	__m256 v124 = data[124*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v116,v124));
	__m256 v123 = data[123*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v117,v123)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v118,v122)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v119,v121)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v120));
	data[60*stride] = acc1;

	// upper band :: ix=60
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v118,v124));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v119,v123)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v120,v122)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v121));
	__m256 v188 = data[188*stride];
	data[188*stride] = acc1;

	// lower band :: ix=61
	__m256 v126 = data[126*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v118,v126));
	__m256 v125 = data[125*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v119,v125)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v120,v124)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v121,v123)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v122));
	data[61*stride] = acc1;

	// upper band :: ix=61
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v120,v126));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v121,v125)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v122,v124)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v123));
	__m256 v189 = data[189*stride];
	data[189*stride] = acc1;

	// lower band :: ix=62
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v120,v128));
	__m256 v127 = data[127*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v121,v127)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v122,v126)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v123,v125)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v124));
	data[62*stride] = acc1;

	// upper band :: ix=62
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v122,v128));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v123,v127)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v124,v126)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v125));
	__m256 v190 = data[190*stride];
	data[190*stride] = acc1;

	// lower band :: ix=63
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v122,v130));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v123,v129)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v124,v128)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v125,v127)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v126));
	data[63*stride] = acc1;

	// upper band :: ix=63
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v124,v130));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v125,v129)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v126,v128)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v127));
	__m256 v191 = data[191*stride];
	data[191*stride] = acc1;

	// lower band :: ix=64
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v124,v132));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v125,v131)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v126,v130)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v127,v129)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v128));
	data[64*stride] = acc1;

	// upper band :: ix=64
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v126,v132));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v127,v131)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v128,v130)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v129));
	__m256 v192 = data[192*stride];
	data[192*stride] = acc1;

	// lower band :: ix=65
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v126,v134));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v127,v133)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v128,v132)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v129,v131)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v130));
	data[65*stride] = acc1;

	// upper band :: ix=65
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v128,v134));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v129,v133)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v130,v132)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v131));
	__m256 v193 = data[193*stride];
	data[193*stride] = acc1;

	// lower band :: ix=66
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v128,v136));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v129,v135)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v130,v134)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v131,v133)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v132));
	data[66*stride] = acc1;

	// upper band :: ix=66
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v130,v136));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v131,v135)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v132,v134)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v133));
	__m256 v194 = data[194*stride];
	data[194*stride] = acc1;

	// lower band :: ix=67
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v130,v138));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v131,v137)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v132,v136)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v133,v135)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v134));
	data[67*stride] = acc1;

	// upper band :: ix=67
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v132,v138));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v133,v137)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v134,v136)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v135));
	__m256 v195 = data[195*stride];
	data[195*stride] = acc1;

	// lower band :: ix=68
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v132,v140));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v133,v139)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v134,v138)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v135,v137)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v136));
	data[68*stride] = acc1;

	// upper band :: ix=68
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v134,v140));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v135,v139)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v136,v138)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v137));
	__m256 v196 = data[196*stride];
	data[196*stride] = acc1;

	// lower band :: ix=69
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v134,v142));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v135,v141)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v136,v140)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v137,v139)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v138));
	data[69*stride] = acc1;

	// upper band :: ix=69
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v136,v142));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v137,v141)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v138,v140)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v139));
	__m256 v197 = data[197*stride];
	data[197*stride] = acc1;

	// lower band :: ix=70
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v136,v144));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v137,v143)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v138,v142)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v139,v141)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v140));
	data[70*stride] = acc1;

	// upper band :: ix=70
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v138,v144));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v139,v143)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v140,v142)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v141));
	__m256 v198 = data[198*stride];
	data[198*stride] = acc1;

	// lower band :: ix=71
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v138,v146));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v139,v145)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v140,v144)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v141,v143)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v142));
	data[71*stride] = acc1;

	// upper band :: ix=71
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v140,v146));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v141,v145)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v142,v144)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v143));
	__m256 v199 = data[199*stride];
	data[199*stride] = acc1;

	// lower band :: ix=72
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v140,v148));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v141,v147)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v142,v146)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v143,v145)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v144));
	data[72*stride] = acc1;

	// upper band :: ix=72
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v142,v148));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v143,v147)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v144,v146)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v145));
	__m256 v200 = data[200*stride];
	data[200*stride] = acc1;

	// lower band :: ix=73
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v142,v150));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v143,v149)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v144,v148)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v145,v147)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v146));
	data[73*stride] = acc1;

	// upper band :: ix=73
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v144,v150));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v145,v149)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v146,v148)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v147));
	__m256 v201 = data[201*stride];
	data[201*stride] = acc1;

	// lower band :: ix=74
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v144,v152));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v145,v151)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v146,v150)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v147,v149)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v148));
	data[74*stride] = acc1;

	// upper band :: ix=74
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v146,v152));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v147,v151)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v148,v150)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v149));
	__m256 v202 = data[202*stride];
	data[202*stride] = acc1;

	// lower band :: ix=75
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v146,v154));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v147,v153)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v148,v152)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v149,v151)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v150));
	data[75*stride] = acc1;

	// upper band :: ix=75
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v148,v154));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v149,v153)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v150,v152)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v151));
	__m256 v203 = data[203*stride];
	data[203*stride] = acc1;

	// lower band :: ix=76
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v148,v156));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v149,v155)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v150,v154)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v151,v153)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v152));
	data[76*stride] = acc1;

	// upper band :: ix=76
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v150,v156));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v151,v155)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v152,v154)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v153));
	__m256 v204 = data[204*stride];
	data[204*stride] = acc1;

	// lower band :: ix=77
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v150,v158));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v151,v157)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v152,v156)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v153,v155)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v154));
	data[77*stride] = acc1;

	// upper band :: ix=77
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v152,v158));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v153,v157)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v154,v156)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v155));
	__m256 v205 = data[205*stride];
	data[205*stride] = acc1;

	// lower band :: ix=78
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v152,v160));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v153,v159)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v154,v158)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v155,v157)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v156));
	data[78*stride] = acc1;

	// upper band :: ix=78
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v154,v160));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v155,v159)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v156,v158)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v157));
	__m256 v206 = data[206*stride];
	data[206*stride] = acc1;

	// lower band :: ix=79
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v154,v162));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v155,v161)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v156,v160)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v157,v159)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v158));
	data[79*stride] = acc1;

	// upper band :: ix=79
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v156,v162));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v157,v161)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v158,v160)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v159));
	__m256 v207 = data[207*stride];
	data[207*stride] = acc1;

	// lower band :: ix=80
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v156,v164));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v157,v163)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v158,v162)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v159,v161)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v160));
	data[80*stride] = acc1;

	// upper band :: ix=80
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v158,v164));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v159,v163)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v160,v162)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v161));
	__m256 v208 = data[208*stride];
	data[208*stride] = acc1;

	// lower band :: ix=81
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v158,v166));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v159,v165)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v160,v164)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v161,v163)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v162));
	data[81*stride] = acc1;

	// upper band :: ix=81
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v160,v166));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v161,v165)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v162,v164)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v163));
	__m256 v209 = data[209*stride];
	data[209*stride] = acc1;

	// lower band :: ix=82
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v160,v168));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v161,v167)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v162,v166)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v163,v165)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v164));
	data[82*stride] = acc1;

	// upper band :: ix=82
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v162,v168));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v163,v167)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v164,v166)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v165));
	__m256 v210 = data[210*stride];
	data[210*stride] = acc1;

	// lower band :: ix=83
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v162,v170));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v163,v169)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v164,v168)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v165,v167)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v166));
	data[83*stride] = acc1;

	// upper band :: ix=83
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v164,v170));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v165,v169)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v166,v168)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v167));
	__m256 v211 = data[211*stride];
	data[211*stride] = acc1;

	// lower band :: ix=84
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v164,v172));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v165,v171)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v166,v170)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v167,v169)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v168));
	data[84*stride] = acc1;

	// upper band :: ix=84
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v166,v172));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v167,v171)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v168,v170)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v169));
	__m256 v212 = data[212*stride];
	data[212*stride] = acc1;

	// lower band :: ix=85
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v166,v174));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v167,v173)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v168,v172)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v169,v171)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v170));
	data[85*stride] = acc1;

	// upper band :: ix=85
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v168,v174));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v169,v173)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v170,v172)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v171));
	__m256 v213 = data[213*stride];
	data[213*stride] = acc1;

	// lower band :: ix=86
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v168,v176));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v169,v175)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v170,v174)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v171,v173)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v172));
	data[86*stride] = acc1;

	// upper band :: ix=86
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v170,v176));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v171,v175)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v172,v174)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v173));
	__m256 v214 = data[214*stride];
	data[214*stride] = acc1;

	// lower band :: ix=87
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v170,v178));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v171,v177)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v172,v176)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v173,v175)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v174));
	data[87*stride] = acc1;

	// upper band :: ix=87
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v172,v178));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v173,v177)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v174,v176)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v175));
	__m256 v215 = data[215*stride];
	data[215*stride] = acc1;

	// lower band :: ix=88
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v172,v180));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v173,v179)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v174,v178)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v175,v177)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v176));
	data[88*stride] = acc1;

	// upper band :: ix=88
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v174,v180));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v175,v179)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v176,v178)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v177));
	__m256 v216 = data[216*stride];
	data[216*stride] = acc1;

	// lower band :: ix=89
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v174,v182));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v175,v181)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v176,v180)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v177,v179)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v178));
	data[89*stride] = acc1;

	// upper band :: ix=89
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v176,v182));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v177,v181)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v178,v180)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v179));
	__m256 v217 = data[217*stride];
	data[217*stride] = acc1;

	// lower band :: ix=90
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v176,v184));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v177,v183)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v178,v182)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v179,v181)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v180));
	data[90*stride] = acc1;

	// upper band :: ix=90
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v178,v184));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v179,v183)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v180,v182)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v181));
	__m256 v218 = data[218*stride];
	data[218*stride] = acc1;

	// lower band :: ix=91
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v178,v186));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v179,v185)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v180,v184)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v181,v183)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v182));
	data[91*stride] = acc1;

	// upper band :: ix=91
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v180,v186));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v181,v185)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v182,v184)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v183));
	__m256 v219 = data[219*stride];
	data[219*stride] = acc1;

	// lower band :: ix=92
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v180,v188));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v181,v187)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v182,v186)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v183,v185)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v184));
	data[92*stride] = acc1;

	// upper band :: ix=92
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v182,v188));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v183,v187)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v184,v186)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v185));
	__m256 v220 = data[220*stride];
	data[220*stride] = acc1;

	// lower band :: ix=93
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v182,v190));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v183,v189)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v184,v188)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v185,v187)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v186));
	data[93*stride] = acc1;

	// upper band :: ix=93
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v184,v190));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v185,v189)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v186,v188)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v187));
	__m256 v221 = data[221*stride];
	data[221*stride] = acc1;

	// lower band :: ix=94
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v184,v192));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v185,v191)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v186,v190)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v187,v189)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v188));
	data[94*stride] = acc1;

	// upper band :: ix=94
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v186,v192));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v187,v191)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v188,v190)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v189));
	__m256 v222 = data[222*stride];
	data[222*stride] = acc1;

	// lower band :: ix=95
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v186,v194));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v187,v193)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v188,v192)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v189,v191)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v190));
	data[95*stride] = acc1;

	// upper band :: ix=95
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v188,v194));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v189,v193)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v190,v192)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v191));
	__m256 v223 = data[223*stride];
	data[223*stride] = acc1;

	// lower band :: ix=96
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v188,v196));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v189,v195)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v190,v194)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v191,v193)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v192));
	data[96*stride] = acc1;

	// upper band :: ix=96
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v190,v196));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v191,v195)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v192,v194)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v193));
	__m256 v224 = data[224*stride];
	data[224*stride] = acc1;

	// lower band :: ix=97
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v190,v198));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v191,v197)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v192,v196)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v193,v195)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v194));
	data[97*stride] = acc1;

	// upper band :: ix=97
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v192,v198));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v193,v197)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v194,v196)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v195));
	__m256 v225 = data[225*stride];
	data[225*stride] = acc1;

	// lower band :: ix=98
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v192,v200));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v193,v199)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v194,v198)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v195,v197)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v196));
	data[98*stride] = acc1;

	// upper band :: ix=98
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v194,v200));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v195,v199)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v196,v198)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v197));
	__m256 v226 = data[226*stride];
	data[226*stride] = acc1;

	// lower band :: ix=99
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v194,v202));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v195,v201)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v196,v200)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v197,v199)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v198));
	data[99*stride] = acc1;

	// upper band :: ix=99
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v196,v202));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v197,v201)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v198,v200)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v199));
	__m256 v227 = data[227*stride];
	data[227*stride] = acc1;

	// lower band :: ix=100
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v196,v204));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v197,v203)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v198,v202)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v199,v201)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v200));
	data[100*stride] = acc1;

	// upper band :: ix=100
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v198,v204));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v199,v203)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v200,v202)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v201));
	__m256 v228 = data[228*stride];
	data[228*stride] = acc1;

	// lower band :: ix=101
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v198,v206));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v199,v205)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v200,v204)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v201,v203)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v202));
	data[101*stride] = acc1;

	// upper band :: ix=101
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v200,v206));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v201,v205)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v202,v204)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v203));
	__m256 v229 = data[229*stride];
	data[229*stride] = acc1;

	// lower band :: ix=102
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v200,v208));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v201,v207)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v202,v206)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v203,v205)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v204));
	data[102*stride] = acc1;

	// upper band :: ix=102
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v202,v208));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v203,v207)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v204,v206)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v205));
	__m256 v230 = data[230*stride];
	data[230*stride] = acc1;

	// lower band :: ix=103
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v202,v210));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v203,v209)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v204,v208)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v205,v207)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v206));
	data[103*stride] = acc1;

	// upper band :: ix=103
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v204,v210));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v205,v209)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v206,v208)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v207));
	__m256 v231 = data[231*stride];
	data[231*stride] = acc1;

	// lower band :: ix=104
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v204,v212));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v205,v211)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v206,v210)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v207,v209)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v208));
	data[104*stride] = acc1;

	// upper band :: ix=104
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v206,v212));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v207,v211)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v208,v210)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v209));
	__m256 v232 = data[232*stride];
	data[232*stride] = acc1;

	// lower band :: ix=105
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v206,v214));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v207,v213)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v208,v212)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v209,v211)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v210));
	data[105*stride] = acc1;

	// upper band :: ix=105
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v208,v214));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v209,v213)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v210,v212)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v211));
	__m256 v233 = data[233*stride];
	data[233*stride] = acc1;

	// lower band :: ix=106
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v208,v216));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v209,v215)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v210,v214)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v211,v213)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v212));
	data[106*stride] = acc1;

	// upper band :: ix=106
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v210,v216));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v211,v215)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v212,v214)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v213));
	__m256 v234 = data[234*stride];
	data[234*stride] = acc1;

	// lower band :: ix=107
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v210,v218));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v211,v217)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v212,v216)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v213,v215)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v214));
	data[107*stride] = acc1;

	// upper band :: ix=107
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v212,v218));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v213,v217)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v214,v216)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v215));
	__m256 v235 = data[235*stride];
	data[235*stride] = acc1;

	// lower band :: ix=108
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v212,v220));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v213,v219)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v214,v218)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v215,v217)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v216));
	data[108*stride] = acc1;

	// upper band :: ix=108
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v214,v220));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v215,v219)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v216,v218)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v217));
	__m256 v236 = data[236*stride];
	data[236*stride] = acc1;

	// lower band :: ix=109
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v214,v222));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v215,v221)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v216,v220)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v217,v219)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v218));
	data[109*stride] = acc1;

	// upper band :: ix=109
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v216,v222));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v217,v221)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v218,v220)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v219));
	__m256 v237 = data[237*stride];
	data[237*stride] = acc1;

	// lower band :: ix=110
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v216,v224));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v217,v223)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v218,v222)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v219,v221)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v220));
	data[110*stride] = acc1;

	// upper band :: ix=110
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v218,v224));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v219,v223)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v220,v222)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v221));
	__m256 v238 = data[238*stride];
	data[238*stride] = acc1;

	// lower band :: ix=111
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v218,v226));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v219,v225)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v220,v224)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v221,v223)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v222));
	data[111*stride] = acc1;

	// upper band :: ix=111
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v220,v226));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v221,v225)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v222,v224)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v223));
	__m256 v239 = data[239*stride];
	data[239*stride] = acc1;

	// lower band :: ix=112
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v220,v228));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v221,v227)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v222,v226)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v223,v225)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v224));
	data[112*stride] = acc1;

	// upper band :: ix=112
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v222,v228));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v223,v227)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v224,v226)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v225));
	__m256 v240 = data[240*stride];
	data[240*stride] = acc1;

	// lower band :: ix=113
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v222,v230));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v223,v229)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v224,v228)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v225,v227)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v226));
	data[113*stride] = acc1;

	// upper band :: ix=113
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v224,v230));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v225,v229)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v226,v228)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v227));
	__m256 v241 = data[241*stride];
	data[241*stride] = acc1;

	// lower band :: ix=114
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v224,v232));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v225,v231)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v226,v230)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v227,v229)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v228));
	data[114*stride] = acc1;

	// upper band :: ix=114
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v226,v232));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v227,v231)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v228,v230)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v229));
	__m256 v242 = data[242*stride];
	data[242*stride] = acc1;

	// lower band :: ix=115
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v226,v234));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v227,v233)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v228,v232)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v229,v231)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v230));
	data[115*stride] = acc1;

	// upper band :: ix=115
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v228,v234));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v229,v233)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v230,v232)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v231));
	__m256 v243 = data[243*stride];
	data[243*stride] = acc1;

	// lower band :: ix=116
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v228,v236));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v229,v235)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v230,v234)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v231,v233)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v232));
	data[116*stride] = acc1;

	// upper band :: ix=116
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v230,v236));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v231,v235)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v232,v234)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v233));
	__m256 v244 = data[244*stride];
	data[244*stride] = acc1;

	// lower band :: ix=117
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v230,v238));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v231,v237)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v232,v236)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v233,v235)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v234));
	data[117*stride] = acc1;

	// upper band :: ix=117
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v232,v238));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v233,v237)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v234,v236)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v235));
	__m256 v245 = data[245*stride];
	data[245*stride] = acc1;

	// lower band :: ix=118
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v232,v240));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v233,v239)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v234,v238)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v235,v237)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v236));
	data[118*stride] = acc1;

	// upper band :: ix=118
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v234,v240));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v235,v239)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v236,v238)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v237));
	__m256 v246 = data[246*stride];
	data[246*stride] = acc1;

	// lower band :: ix=119
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v234,v242));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v235,v241)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v236,v240)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v237,v239)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v238));
	data[119*stride] = acc1;

	// upper band :: ix=119
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v236,v242));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v237,v241)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v238,v240)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v239));
	__m256 v247 = data[247*stride];
	data[247*stride] = acc1;

	// lower band :: ix=120
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v236,v244));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v237,v243)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v238,v242)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v239,v241)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v240));
	data[120*stride] = acc1;

	// upper band :: ix=120
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v238,v244));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v239,v243)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v240,v242)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v241));
	__m256 v248 = data[248*stride];
	data[248*stride] = acc1;

	// lower band :: ix=121
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v238,v246));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v239,v245)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v240,v244)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v241,v243)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v242));
	data[121*stride] = acc1;

	// upper band :: ix=121
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v240,v246));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v241,v245)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v242,v244)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v243));
	__m256 v249 = data[249*stride];
	data[249*stride] = acc1;

	// lower band :: ix=122
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v240,v248));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v241,v247)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v242,v246)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v243,v245)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v244));
	data[122*stride] = acc1;

	// upper band :: ix=122
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v242,v248));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v243,v247)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v244,v246)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v245));
	__m256 v250 = data[250*stride];
	data[250*stride] = acc1;

	// lower band :: ix=123
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v242,v250));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v243,v249)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v244,v248)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v245,v247)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v246));
	data[123*stride] = acc1;

	// upper band :: ix=123
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v244,v250));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v245,v249)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v246,v248)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v247));
	__m256 v251 = data[251*stride];
	data[251*stride] = acc1;

	// lower band :: ix=124
	__m256 v252 = data[252*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v244,v252));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v245,v251)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v246,v250)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v247,v249)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v248));
	data[124*stride] = acc1;

	// upper band :: ix=124
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v246,v252));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v247,v251)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v248,v250)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v249));
	data[252*stride] = acc1;

	// lower band :: ix=125
	__m256 v254 = data[254*stride];
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v246,v254));
	__m256 v253 = data[253*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v247,v253)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v248,v252)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v249,v251)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v250));
	data[125*stride] = acc1;

	// upper band :: ix=125
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v248,v254));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v249,v253)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v250,v252)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v251));
	data[253*stride] = acc1;

	// lower band :: ix=126
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v248,v254));
	__m256 v255 = data[255*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v249,v255)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v250,v254)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v251,v253)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v252));
	data[126*stride] = acc1;

	// upper band :: ix=126
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v250,v254));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v251,v255)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v252,v254)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v253));
	data[254*stride] = acc1;

	// lower band :: ix=127
	acc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v250,v252));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v251,v253)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v252,v254)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v253,v255)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v254));
	data[127*stride] = acc1;

	// upper band :: ix=127
	acc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v252,v252));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v253,v253)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v254,v254)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v255));
	data[255*stride] = acc1;
}

#endif
