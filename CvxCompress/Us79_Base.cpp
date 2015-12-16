/*!
 * Don't edit this code, it was automatically generated.
 * Base functions for wavelet transforms of length 2 to 256.
 */
#include <immintrin.h>

/*
 * Define coefficients for Antonini 7-9 tap filter.
 */
#define sl0  7.884856164056601e-001f
#define sl1  4.180922732222101e-001f
#define sl2 -4.068941760955800e-002f
#define sl3 -6.453888262893799e-002f
#define sh0  8.526986790094000e-001f
#define sh1 -3.774028556126500e-001f
#define sh2 -1.106244044184200e-001f
#define sh3  2.384946501938001e-002f
#define sh4  3.782845550699501e-002f

#define _mm_sl0 _mm256_set1_ps(sl0)
#define _mm_sl1 _mm256_set1_ps(sl1)
#define _mm_sl2 _mm256_set1_ps(sl2)
#define _mm_sl3 _mm256_set1_ps(sl3)
#define _mm_sh0 _mm256_set1_ps(sh0)

#define _mm_sh1 _mm256_set1_ps(sh1)
#define _mm_sh2 _mm256_set1_ps(sh2)
#define _mm_sh3 _mm256_set1_ps(sh3)
#define _mm_sh4 _mm256_set1_ps(sh4)

#ifdef __AVX2__

static inline void _Us79_AVX_2(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v1 = data[1*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v1,v1));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v0,v0),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v1,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v0,acc1);
	data[0*stride] = acc1;

	// odd samples :: k=1
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v1,v1));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v0,v0),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v1,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v0,v0),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v1,acc1);
	data[1*stride] = acc1;
}

static inline void _Us79_AVX_4(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v3 = data[3*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v3,v3));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v1,v1),acc1);
	__m256 v2 = data[2*stride];
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v2,v2),acc1);
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl0,v0,acc1);
	data[0*stride] = acc1;

	// odd samples :: k=1
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v3,v2));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v1,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v2,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v0,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v2,acc1);
	data[1*stride] = acc1;

	// even samples :: k=2
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v2,v2));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v0,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v2,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v1,acc1);
	data[2*stride] = acc1;

	// odd samples :: k=3
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v2,v2));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v0,v0),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v2,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v1,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v3,acc1);
	data[3*stride] = acc1;
}

static inline void _Us79_AVX_8(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v5 = data[5*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v5,v5));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v1,v1),acc1);
	__m256 v4 = data[4*stride];
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v4,v4),acc1);
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl0,v0,acc1);
	data[0*stride] = acc1;

	// odd samples :: k=1
	__m256 v6 = data[6*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v5,v6));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v1,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v4,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v0,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v4,acc1);
	data[1*stride] = acc1;

	// even samples :: k=2
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v4,v6));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v0,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v4,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v1,acc1);
	data[2*stride] = acc1;

	// odd samples :: k=3
	__m256 v7 = data[7*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v4,v7));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v0,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v4,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v1,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v5,acc1);
	data[3*stride] = acc1;

	// even samples :: k=4
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v4,v7));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v5,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v2,acc1);
	data[4*stride] = acc1;

	// odd samples :: k=5
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v4,v6));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v5,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v2,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v6,acc1);
	data[5*stride] = acc1;

	// even samples :: k=6
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v5,v6));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v2,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v6,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v3,acc1);
	data[6*stride] = acc1;

	// odd samples :: k=7
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v5,v5));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v2,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v6,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v3,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v7,acc1);
	data[7*stride] = acc1;
}

static inline void _Us79_AVX_16(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v9 = data[9*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v9,v9));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v1,v1),acc1);
	__m256 v8 = data[8*stride];
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v8,v8),acc1);
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl0,v0,acc1);
	data[0*stride] = acc1;

	// odd samples :: k=1
	__m256 v10 = data[10*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v9,v10));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v1,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v8,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v0,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v8,acc1);
	data[1*stride] = acc1;

	// even samples :: k=2
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v8,v10));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v0,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v8,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v1,acc1);
	data[2*stride] = acc1;

	// odd samples :: k=3
	__m256 v11 = data[11*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v8,v11));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v0,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v8,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v1,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v9,acc1);
	data[3*stride] = acc1;

	// even samples :: k=4
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v8,v11));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v9,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v2,acc1);
	__m256 v4 = data[4*stride];
	data[4*stride] = acc1;

	// odd samples :: k=5
	__m256 v12 = data[12*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v8,v12));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v1,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v9,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v2,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v10,acc1);
	__m256 v5 = data[5*stride];
	data[5*stride] = acc1;

	// even samples :: k=6
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v9,v12));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v2,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v10,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v3,acc1);
	__m256 v6 = data[6*stride];
	data[6*stride] = acc1;

	// odd samples :: k=7
	__m256 v13 = data[13*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v9,v13));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v2,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v10,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v3,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v11,acc1);
	__m256 v7 = data[7*stride];
	data[7*stride] = acc1;

	// even samples :: k=8
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v10,v13));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v3,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v11,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v4,acc1);
	data[8*stride] = acc1;

	// odd samples :: k=9
	__m256 v14 = data[14*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v10,v14));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v3,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v11,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v4,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v12,acc1);
	data[9*stride] = acc1;

	// even samples :: k=10
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v11,v14));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v4,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v12,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v5,acc1);
	data[10*stride] = acc1;

	// odd samples :: k=11
	__m256 v15 = data[15*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v11,v15));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v4,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v12,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v5,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v13,acc1);
	data[11*stride] = acc1;

	// even samples :: k=12
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v12,v15));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v5,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v13,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v6,acc1);
	data[12*stride] = acc1;

	// odd samples :: k=13
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v12,v14));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v5,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v13,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v6,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v14,acc1);
	data[13*stride] = acc1;

	// even samples :: k=14
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v13,v14));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v6,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v14,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v7,acc1);
	data[14*stride] = acc1;

	// odd samples :: k=15
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v13,v13));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v6,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v14,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v7,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v15,acc1);
	data[15*stride] = acc1;
}

static inline void _Us79_AVX_32(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v17 = data[17*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v17,v17));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v1,v1),acc1);
	__m256 v16 = data[16*stride];
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v16,v16),acc1);
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl0,v0,acc1);
	data[0*stride] = acc1;

	// odd samples :: k=1
	__m256 v18 = data[18*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v17,v18));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v1,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v16,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v0,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v16,acc1);
	data[1*stride] = acc1;

	// even samples :: k=2
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v16,v18));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v0,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v16,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v1,acc1);
	data[2*stride] = acc1;

	// odd samples :: k=3
	__m256 v19 = data[19*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v16,v19));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v0,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v16,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v1,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v17,acc1);
	data[3*stride] = acc1;

	// even samples :: k=4
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v16,v19));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v17,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v2,acc1);
	__m256 v4 = data[4*stride];
	data[4*stride] = acc1;

	// odd samples :: k=5
	__m256 v20 = data[20*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v16,v20));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v1,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v17,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v2,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v18,acc1);
	__m256 v5 = data[5*stride];
	data[5*stride] = acc1;

	// even samples :: k=6
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v17,v20));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v2,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v18,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v3,acc1);
	__m256 v6 = data[6*stride];
	data[6*stride] = acc1;

	// odd samples :: k=7
	__m256 v21 = data[21*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v17,v21));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v2,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v18,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v3,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v19,acc1);
	__m256 v7 = data[7*stride];
	data[7*stride] = acc1;

	// even samples :: k=8
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v18,v21));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v3,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v19,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v4,acc1);
	__m256 v8 = data[8*stride];
	data[8*stride] = acc1;

	// odd samples :: k=9
	__m256 v22 = data[22*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v18,v22));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v3,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v19,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v4,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v20,acc1);
	__m256 v9 = data[9*stride];
	data[9*stride] = acc1;

	// even samples :: k=10
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v19,v22));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v4,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v20,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v5,acc1);
	__m256 v10 = data[10*stride];
	data[10*stride] = acc1;

	// odd samples :: k=11
	__m256 v23 = data[23*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v19,v23));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v4,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v20,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v5,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v21,acc1);
	__m256 v11 = data[11*stride];
	data[11*stride] = acc1;

	// even samples :: k=12
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v20,v23));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v5,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v21,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v6,acc1);
	__m256 v12 = data[12*stride];
	data[12*stride] = acc1;

	// odd samples :: k=13
	__m256 v24 = data[24*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v20,v24));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v5,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v21,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v6,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v22,acc1);
	__m256 v13 = data[13*stride];
	data[13*stride] = acc1;

	// even samples :: k=14
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v21,v24));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v6,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v22,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v7,acc1);
	__m256 v14 = data[14*stride];
	data[14*stride] = acc1;

	// odd samples :: k=15
	__m256 v25 = data[25*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v21,v25));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v6,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v22,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v7,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v23,acc1);
	__m256 v15 = data[15*stride];
	data[15*stride] = acc1;

	// even samples :: k=16
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v22,v25));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v7,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v23,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v8,acc1);
	data[16*stride] = acc1;

	// odd samples :: k=17
	__m256 v26 = data[26*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v22,v26));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v7,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v23,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v8,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v24,acc1);
	data[17*stride] = acc1;

	// even samples :: k=18
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v23,v26));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v8,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v24,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v9,acc1);
	data[18*stride] = acc1;

	// odd samples :: k=19
	__m256 v27 = data[27*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v23,v27));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v8,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v24,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v9,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v25,acc1);
	data[19*stride] = acc1;

	// even samples :: k=20
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v24,v27));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v9,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v25,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v10,acc1);
	data[20*stride] = acc1;

	// odd samples :: k=21
	__m256 v28 = data[28*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v24,v28));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v9,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v25,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v10,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v26,acc1);
	data[21*stride] = acc1;

	// even samples :: k=22
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v25,v28));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v10,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v26,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v11,acc1);
	data[22*stride] = acc1;

	// odd samples :: k=23
	__m256 v29 = data[29*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v25,v29));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v10,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v26,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v11,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v27,acc1);
	data[23*stride] = acc1;

	// even samples :: k=24
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v26,v29));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v11,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v27,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v12,acc1);
	data[24*stride] = acc1;

	// odd samples :: k=25
	__m256 v30 = data[30*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v26,v30));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v11,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v27,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v12,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v28,acc1);
	data[25*stride] = acc1;

	// even samples :: k=26
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v27,v30));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v12,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v28,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v13,acc1);
	data[26*stride] = acc1;

	// odd samples :: k=27
	__m256 v31 = data[31*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v27,v31));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v12,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v28,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v13,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v29,acc1);
	data[27*stride] = acc1;

	// even samples :: k=28
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v28,v31));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v13,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v29,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v14,acc1);
	data[28*stride] = acc1;

	// odd samples :: k=29
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v28,v30));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v13,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v29,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v14,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v30,acc1);
	data[29*stride] = acc1;

	// even samples :: k=30
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v29,v30));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v14,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v30,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v15,acc1);
	data[30*stride] = acc1;

	// odd samples :: k=31
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v29,v29));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v14,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v30,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v15,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v31,acc1);
	data[31*stride] = acc1;
}

static void _Us79_AVX_64(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v33 = data[33*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v33,v33));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v1,v1),acc1);
	__m256 v32 = data[32*stride];
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v32,v32),acc1);
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl0,v0,acc1);
	data[0*stride] = acc1;

	// odd samples :: k=1
	__m256 v34 = data[34*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v33,v34));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v1,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v32,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v0,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v32,acc1);
	data[1*stride] = acc1;

	// even samples :: k=2
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v32,v34));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v0,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v32,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v1,acc1);
	data[2*stride] = acc1;

	// odd samples :: k=3
	__m256 v35 = data[35*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v32,v35));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v0,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v32,v34),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v1,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v33,acc1);
	data[3*stride] = acc1;

	// even samples :: k=4
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v32,v35));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v33,v34),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v2,acc1);
	__m256 v4 = data[4*stride];
	data[4*stride] = acc1;

	// odd samples :: k=5
	__m256 v36 = data[36*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v32,v36));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v1,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v33,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v2,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v34,acc1);
	__m256 v5 = data[5*stride];
	data[5*stride] = acc1;

	// even samples :: k=6
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v33,v36));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v2,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v34,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v3,acc1);
	__m256 v6 = data[6*stride];
	data[6*stride] = acc1;

	// odd samples :: k=7
	__m256 v37 = data[37*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v33,v37));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v2,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v34,v36),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v3,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v35,acc1);
	__m256 v7 = data[7*stride];
	data[7*stride] = acc1;

	// even samples :: k=8
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v34,v37));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v3,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v35,v36),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v4,acc1);
	__m256 v8 = data[8*stride];
	data[8*stride] = acc1;

	// odd samples :: k=9
	__m256 v38 = data[38*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v34,v38));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v3,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v35,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v4,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v36,acc1);
	__m256 v9 = data[9*stride];
	data[9*stride] = acc1;

	// even samples :: k=10
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v35,v38));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v4,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v36,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v5,acc1);
	__m256 v10 = data[10*stride];
	data[10*stride] = acc1;

	// odd samples :: k=11
	__m256 v39 = data[39*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v35,v39));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v4,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v36,v38),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v5,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v37,acc1);
	__m256 v11 = data[11*stride];
	data[11*stride] = acc1;

	// even samples :: k=12
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v36,v39));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v5,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v37,v38),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v6,acc1);
	__m256 v12 = data[12*stride];
	data[12*stride] = acc1;

	// odd samples :: k=13
	__m256 v40 = data[40*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v36,v40));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v5,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v37,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v6,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v38,acc1);
	__m256 v13 = data[13*stride];
	data[13*stride] = acc1;

	// even samples :: k=14
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v37,v40));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v6,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v38,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v7,acc1);
	__m256 v14 = data[14*stride];
	data[14*stride] = acc1;

	// odd samples :: k=15
	__m256 v41 = data[41*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v37,v41));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v6,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v38,v40),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v7,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v39,acc1);
	__m256 v15 = data[15*stride];
	data[15*stride] = acc1;

	// even samples :: k=16
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v38,v41));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v7,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v39,v40),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v8,acc1);
	__m256 v16 = data[16*stride];
	data[16*stride] = acc1;

	// odd samples :: k=17
	__m256 v42 = data[42*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v38,v42));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v7,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v39,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v8,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v40,acc1);
	__m256 v17 = data[17*stride];
	data[17*stride] = acc1;

	// even samples :: k=18
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v39,v42));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v8,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v40,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v9,acc1);
	__m256 v18 = data[18*stride];
	data[18*stride] = acc1;

	// odd samples :: k=19
	__m256 v43 = data[43*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v39,v43));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v8,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v40,v42),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v9,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v41,acc1);
	__m256 v19 = data[19*stride];
	data[19*stride] = acc1;

	// even samples :: k=20
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v40,v43));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v9,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v41,v42),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v10,acc1);
	__m256 v20 = data[20*stride];
	data[20*stride] = acc1;

	// odd samples :: k=21
	__m256 v44 = data[44*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v40,v44));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v9,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v41,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v10,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v42,acc1);
	__m256 v21 = data[21*stride];
	data[21*stride] = acc1;

	// even samples :: k=22
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v41,v44));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v10,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v42,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v11,acc1);
	__m256 v22 = data[22*stride];
	data[22*stride] = acc1;

	// odd samples :: k=23
	__m256 v45 = data[45*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v41,v45));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v10,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v42,v44),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v11,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v43,acc1);
	__m256 v23 = data[23*stride];
	data[23*stride] = acc1;

	// even samples :: k=24
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v42,v45));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v11,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v43,v44),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v12,acc1);
	__m256 v24 = data[24*stride];
	data[24*stride] = acc1;

	// odd samples :: k=25
	__m256 v46 = data[46*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v42,v46));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v11,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v43,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v12,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v44,acc1);
	__m256 v25 = data[25*stride];
	data[25*stride] = acc1;

	// even samples :: k=26
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v43,v46));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v12,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v44,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v13,acc1);
	__m256 v26 = data[26*stride];
	data[26*stride] = acc1;

	// odd samples :: k=27
	__m256 v47 = data[47*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v43,v47));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v12,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v44,v46),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v13,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v45,acc1);
	__m256 v27 = data[27*stride];
	data[27*stride] = acc1;

	// even samples :: k=28
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v44,v47));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v13,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v45,v46),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v14,acc1);
	__m256 v28 = data[28*stride];
	data[28*stride] = acc1;

	// odd samples :: k=29
	__m256 v48 = data[48*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v44,v48));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v13,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v45,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v14,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v46,acc1);
	__m256 v29 = data[29*stride];
	data[29*stride] = acc1;

	// even samples :: k=30
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v45,v48));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v14,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v46,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v15,acc1);
	__m256 v30 = data[30*stride];
	data[30*stride] = acc1;

	// odd samples :: k=31
	__m256 v49 = data[49*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v45,v49));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v14,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v46,v48),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v15,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v47,acc1);
	__m256 v31 = data[31*stride];
	data[31*stride] = acc1;

	// even samples :: k=32
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v46,v49));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v15,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v47,v48),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v16,acc1);
	data[32*stride] = acc1;

	// odd samples :: k=33
	__m256 v50 = data[50*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v46,v50));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v15,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v47,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v16,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v48,acc1);
	data[33*stride] = acc1;

	// even samples :: k=34
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v47,v50));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v16,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v48,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v17,acc1);
	data[34*stride] = acc1;

	// odd samples :: k=35
	__m256 v51 = data[51*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v47,v51));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v16,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v48,v50),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v17,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v49,acc1);
	data[35*stride] = acc1;

	// even samples :: k=36
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v48,v51));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v17,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v49,v50),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v18,acc1);
	data[36*stride] = acc1;

	// odd samples :: k=37
	__m256 v52 = data[52*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v48,v52));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v17,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v49,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v18,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v50,acc1);
	data[37*stride] = acc1;

	// even samples :: k=38
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v49,v52));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v18,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v50,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v19,acc1);
	data[38*stride] = acc1;

	// odd samples :: k=39
	__m256 v53 = data[53*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v49,v53));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v18,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v50,v52),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v19,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v51,acc1);
	data[39*stride] = acc1;

	// even samples :: k=40
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v50,v53));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v19,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v51,v52),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v20,acc1);
	data[40*stride] = acc1;

	// odd samples :: k=41
	__m256 v54 = data[54*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v50,v54));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v19,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v51,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v20,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v52,acc1);
	data[41*stride] = acc1;

	// even samples :: k=42
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v51,v54));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v20,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v52,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v21,acc1);
	data[42*stride] = acc1;

	// odd samples :: k=43
	__m256 v55 = data[55*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v51,v55));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v20,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v52,v54),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v21,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v53,acc1);
	data[43*stride] = acc1;

	// even samples :: k=44
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v52,v55));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v21,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v53,v54),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v22,acc1);
	data[44*stride] = acc1;

	// odd samples :: k=45
	__m256 v56 = data[56*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v52,v56));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v21,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v53,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v22,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v54,acc1);
	data[45*stride] = acc1;

	// even samples :: k=46
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v53,v56));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v22,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v54,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v23,acc1);
	data[46*stride] = acc1;

	// odd samples :: k=47
	__m256 v57 = data[57*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v53,v57));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v22,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v54,v56),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v23,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v55,acc1);
	data[47*stride] = acc1;

	// even samples :: k=48
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v54,v57));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v23,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v55,v56),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v24,acc1);
	data[48*stride] = acc1;

	// odd samples :: k=49
	__m256 v58 = data[58*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v54,v58));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v23,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v55,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v24,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v56,acc1);
	data[49*stride] = acc1;

	// even samples :: k=50
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v55,v58));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v24,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v56,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v25,acc1);
	data[50*stride] = acc1;

	// odd samples :: k=51
	__m256 v59 = data[59*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v55,v59));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v24,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v56,v58),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v25,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v57,acc1);
	data[51*stride] = acc1;

	// even samples :: k=52
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v56,v59));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v25,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v57,v58),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v26,acc1);
	data[52*stride] = acc1;

	// odd samples :: k=53
	__m256 v60 = data[60*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v56,v60));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v25,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v57,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v26,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v58,acc1);
	data[53*stride] = acc1;

	// even samples :: k=54
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v57,v60));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v26,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v58,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v27,acc1);
	data[54*stride] = acc1;

	// odd samples :: k=55
	__m256 v61 = data[61*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v57,v61));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v26,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v58,v60),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v27,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v59,acc1);
	data[55*stride] = acc1;

	// even samples :: k=56
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v58,v61));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v27,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v59,v60),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v28,acc1);
	data[56*stride] = acc1;

	// odd samples :: k=57
	__m256 v62 = data[62*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v58,v62));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v27,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v59,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v28,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v60,acc1);
	data[57*stride] = acc1;

	// even samples :: k=58
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v59,v62));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v28,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v60,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v29,acc1);
	data[58*stride] = acc1;

	// odd samples :: k=59
	__m256 v63 = data[63*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v59,v63));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v28,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v60,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v29,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v61,acc1);
	data[59*stride] = acc1;

	// even samples :: k=60
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v60,v63));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v29,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v61,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v30,acc1);
	data[60*stride] = acc1;

	// odd samples :: k=61
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v60,v62));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v29,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v61,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v30,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v62,acc1);
	data[61*stride] = acc1;

	// even samples :: k=62
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v61,v62));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v30,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v62,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v31,acc1);
	data[62*stride] = acc1;

	// odd samples :: k=63
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v61,v61));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v30,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v62,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v31,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v63,acc1);
	data[63*stride] = acc1;
}

static void _Us79_AVX_128(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v65 = data[65*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v65,v65));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v1,v1),acc1);
	__m256 v64 = data[64*stride];
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v64,v64),acc1);
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl0,v0,acc1);
	data[0*stride] = acc1;

	// odd samples :: k=1
	__m256 v66 = data[66*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v65,v66));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v1,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v64,v65),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v0,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v64,acc1);
	data[1*stride] = acc1;

	// even samples :: k=2
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v64,v66));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v0,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v64,v65),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v1,acc1);
	data[2*stride] = acc1;

	// odd samples :: k=3
	__m256 v67 = data[67*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v64,v67));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v0,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v64,v66),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v1,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v65,acc1);
	data[3*stride] = acc1;

	// even samples :: k=4
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v64,v67));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v65,v66),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v2,acc1);
	__m256 v4 = data[4*stride];
	data[4*stride] = acc1;

	// odd samples :: k=5
	__m256 v68 = data[68*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v64,v68));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v1,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v65,v67),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v2,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v66,acc1);
	__m256 v5 = data[5*stride];
	data[5*stride] = acc1;

	// even samples :: k=6
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v65,v68));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v2,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v66,v67),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v3,acc1);
	__m256 v6 = data[6*stride];
	data[6*stride] = acc1;

	// odd samples :: k=7
	__m256 v69 = data[69*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v65,v69));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v2,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v66,v68),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v3,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v67,acc1);
	__m256 v7 = data[7*stride];
	data[7*stride] = acc1;

	// even samples :: k=8
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v66,v69));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v3,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v67,v68),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v4,acc1);
	__m256 v8 = data[8*stride];
	data[8*stride] = acc1;

	// odd samples :: k=9
	__m256 v70 = data[70*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v66,v70));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v3,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v67,v69),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v4,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v68,acc1);
	__m256 v9 = data[9*stride];
	data[9*stride] = acc1;

	// even samples :: k=10
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v67,v70));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v4,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v68,v69),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v5,acc1);
	__m256 v10 = data[10*stride];
	data[10*stride] = acc1;

	// odd samples :: k=11
	__m256 v71 = data[71*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v67,v71));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v4,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v68,v70),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v5,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v69,acc1);
	__m256 v11 = data[11*stride];
	data[11*stride] = acc1;

	// even samples :: k=12
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v68,v71));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v5,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v69,v70),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v6,acc1);
	__m256 v12 = data[12*stride];
	data[12*stride] = acc1;

	// odd samples :: k=13
	__m256 v72 = data[72*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v68,v72));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v5,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v69,v71),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v6,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v70,acc1);
	__m256 v13 = data[13*stride];
	data[13*stride] = acc1;

	// even samples :: k=14
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v69,v72));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v6,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v70,v71),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v7,acc1);
	__m256 v14 = data[14*stride];
	data[14*stride] = acc1;

	// odd samples :: k=15
	__m256 v73 = data[73*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v69,v73));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v6,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v70,v72),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v7,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v71,acc1);
	__m256 v15 = data[15*stride];
	data[15*stride] = acc1;

	// even samples :: k=16
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v70,v73));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v7,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v71,v72),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v8,acc1);
	__m256 v16 = data[16*stride];
	data[16*stride] = acc1;

	// odd samples :: k=17
	__m256 v74 = data[74*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v70,v74));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v7,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v71,v73),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v8,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v72,acc1);
	__m256 v17 = data[17*stride];
	data[17*stride] = acc1;

	// even samples :: k=18
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v71,v74));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v8,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v72,v73),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v9,acc1);
	__m256 v18 = data[18*stride];
	data[18*stride] = acc1;

	// odd samples :: k=19
	__m256 v75 = data[75*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v71,v75));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v8,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v72,v74),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v9,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v73,acc1);
	__m256 v19 = data[19*stride];
	data[19*stride] = acc1;

	// even samples :: k=20
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v72,v75));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v9,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v73,v74),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v10,acc1);
	__m256 v20 = data[20*stride];
	data[20*stride] = acc1;

	// odd samples :: k=21
	__m256 v76 = data[76*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v72,v76));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v9,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v73,v75),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v10,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v74,acc1);
	__m256 v21 = data[21*stride];
	data[21*stride] = acc1;

	// even samples :: k=22
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v73,v76));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v10,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v74,v75),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v11,acc1);
	__m256 v22 = data[22*stride];
	data[22*stride] = acc1;

	// odd samples :: k=23
	__m256 v77 = data[77*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v73,v77));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v10,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v74,v76),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v11,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v75,acc1);
	__m256 v23 = data[23*stride];
	data[23*stride] = acc1;

	// even samples :: k=24
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v74,v77));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v11,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v75,v76),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v12,acc1);
	__m256 v24 = data[24*stride];
	data[24*stride] = acc1;

	// odd samples :: k=25
	__m256 v78 = data[78*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v74,v78));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v11,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v75,v77),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v12,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v76,acc1);
	__m256 v25 = data[25*stride];
	data[25*stride] = acc1;

	// even samples :: k=26
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v75,v78));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v12,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v76,v77),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v13,acc1);
	__m256 v26 = data[26*stride];
	data[26*stride] = acc1;

	// odd samples :: k=27
	__m256 v79 = data[79*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v75,v79));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v12,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v76,v78),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v13,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v77,acc1);
	__m256 v27 = data[27*stride];
	data[27*stride] = acc1;

	// even samples :: k=28
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v76,v79));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v13,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v77,v78),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v14,acc1);
	__m256 v28 = data[28*stride];
	data[28*stride] = acc1;

	// odd samples :: k=29
	__m256 v80 = data[80*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v76,v80));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v13,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v77,v79),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v14,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v78,acc1);
	__m256 v29 = data[29*stride];
	data[29*stride] = acc1;

	// even samples :: k=30
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v77,v80));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v14,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v78,v79),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v15,acc1);
	__m256 v30 = data[30*stride];
	data[30*stride] = acc1;

	// odd samples :: k=31
	__m256 v81 = data[81*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v77,v81));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v14,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v78,v80),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v15,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v79,acc1);
	__m256 v31 = data[31*stride];
	data[31*stride] = acc1;

	// even samples :: k=32
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v78,v81));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v15,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v79,v80),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v16,acc1);
	__m256 v32 = data[32*stride];
	data[32*stride] = acc1;

	// odd samples :: k=33
	__m256 v82 = data[82*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v78,v82));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v15,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v79,v81),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v16,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v80,acc1);
	__m256 v33 = data[33*stride];
	data[33*stride] = acc1;

	// even samples :: k=34
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v79,v82));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v16,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v80,v81),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v17,acc1);
	__m256 v34 = data[34*stride];
	data[34*stride] = acc1;

	// odd samples :: k=35
	__m256 v83 = data[83*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v79,v83));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v16,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v80,v82),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v17,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v81,acc1);
	__m256 v35 = data[35*stride];
	data[35*stride] = acc1;

	// even samples :: k=36
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v80,v83));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v17,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v81,v82),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v18,acc1);
	__m256 v36 = data[36*stride];
	data[36*stride] = acc1;

	// odd samples :: k=37
	__m256 v84 = data[84*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v80,v84));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v17,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v81,v83),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v18,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v82,acc1);
	__m256 v37 = data[37*stride];
	data[37*stride] = acc1;

	// even samples :: k=38
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v81,v84));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v18,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v82,v83),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v19,acc1);
	__m256 v38 = data[38*stride];
	data[38*stride] = acc1;

	// odd samples :: k=39
	__m256 v85 = data[85*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v81,v85));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v18,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v82,v84),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v19,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v83,acc1);
	__m256 v39 = data[39*stride];
	data[39*stride] = acc1;

	// even samples :: k=40
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v82,v85));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v19,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v83,v84),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v20,acc1);
	__m256 v40 = data[40*stride];
	data[40*stride] = acc1;

	// odd samples :: k=41
	__m256 v86 = data[86*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v82,v86));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v19,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v83,v85),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v20,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v84,acc1);
	__m256 v41 = data[41*stride];
	data[41*stride] = acc1;

	// even samples :: k=42
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v83,v86));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v20,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v84,v85),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v21,acc1);
	__m256 v42 = data[42*stride];
	data[42*stride] = acc1;

	// odd samples :: k=43
	__m256 v87 = data[87*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v83,v87));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v20,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v84,v86),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v21,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v85,acc1);
	__m256 v43 = data[43*stride];
	data[43*stride] = acc1;

	// even samples :: k=44
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v84,v87));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v21,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v85,v86),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v22,acc1);
	__m256 v44 = data[44*stride];
	data[44*stride] = acc1;

	// odd samples :: k=45
	__m256 v88 = data[88*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v84,v88));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v21,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v85,v87),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v22,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v86,acc1);
	__m256 v45 = data[45*stride];
	data[45*stride] = acc1;

	// even samples :: k=46
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v85,v88));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v22,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v86,v87),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v23,acc1);
	__m256 v46 = data[46*stride];
	data[46*stride] = acc1;

	// odd samples :: k=47
	__m256 v89 = data[89*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v85,v89));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v22,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v86,v88),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v23,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v87,acc1);
	__m256 v47 = data[47*stride];
	data[47*stride] = acc1;

	// even samples :: k=48
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v86,v89));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v23,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v87,v88),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v24,acc1);
	__m256 v48 = data[48*stride];
	data[48*stride] = acc1;

	// odd samples :: k=49
	__m256 v90 = data[90*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v86,v90));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v23,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v87,v89),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v24,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v88,acc1);
	__m256 v49 = data[49*stride];
	data[49*stride] = acc1;

	// even samples :: k=50
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v87,v90));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v24,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v88,v89),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v25,acc1);
	__m256 v50 = data[50*stride];
	data[50*stride] = acc1;

	// odd samples :: k=51
	__m256 v91 = data[91*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v87,v91));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v24,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v88,v90),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v25,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v89,acc1);
	__m256 v51 = data[51*stride];
	data[51*stride] = acc1;

	// even samples :: k=52
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v88,v91));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v25,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v89,v90),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v26,acc1);
	__m256 v52 = data[52*stride];
	data[52*stride] = acc1;

	// odd samples :: k=53
	__m256 v92 = data[92*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v88,v92));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v25,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v89,v91),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v26,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v90,acc1);
	__m256 v53 = data[53*stride];
	data[53*stride] = acc1;

	// even samples :: k=54
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v89,v92));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v26,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v90,v91),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v27,acc1);
	__m256 v54 = data[54*stride];
	data[54*stride] = acc1;

	// odd samples :: k=55
	__m256 v93 = data[93*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v89,v93));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v26,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v90,v92),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v27,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v91,acc1);
	__m256 v55 = data[55*stride];
	data[55*stride] = acc1;

	// even samples :: k=56
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v90,v93));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v27,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v91,v92),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v28,acc1);
	__m256 v56 = data[56*stride];
	data[56*stride] = acc1;

	// odd samples :: k=57
	__m256 v94 = data[94*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v90,v94));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v27,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v91,v93),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v28,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v92,acc1);
	__m256 v57 = data[57*stride];
	data[57*stride] = acc1;

	// even samples :: k=58
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v91,v94));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v28,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v92,v93),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v29,acc1);
	__m256 v58 = data[58*stride];
	data[58*stride] = acc1;

	// odd samples :: k=59
	__m256 v95 = data[95*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v91,v95));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v28,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v92,v94),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v29,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v93,acc1);
	__m256 v59 = data[59*stride];
	data[59*stride] = acc1;

	// even samples :: k=60
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v92,v95));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v29,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v93,v94),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v30,acc1);
	__m256 v60 = data[60*stride];
	data[60*stride] = acc1;

	// odd samples :: k=61
	__m256 v96 = data[96*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v92,v96));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v29,v32),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v93,v95),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v30,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v94,acc1);
	__m256 v61 = data[61*stride];
	data[61*stride] = acc1;

	// even samples :: k=62
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v93,v96));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v30,v32),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v94,v95),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v31,acc1);
	__m256 v62 = data[62*stride];
	data[62*stride] = acc1;

	// odd samples :: k=63
	__m256 v97 = data[97*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v93,v97));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v30,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v94,v96),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v31,v32),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v95,acc1);
	__m256 v63 = data[63*stride];
	data[63*stride] = acc1;

	// even samples :: k=64
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v94,v97));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v31,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v95,v96),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v32,acc1);
	data[64*stride] = acc1;

	// odd samples :: k=65
	__m256 v98 = data[98*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v94,v98));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v31,v34),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v95,v97),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v32,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v96,acc1);
	data[65*stride] = acc1;

	// even samples :: k=66
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v95,v98));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v32,v34),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v96,v97),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v33,acc1);
	data[66*stride] = acc1;

	// odd samples :: k=67
	__m256 v99 = data[99*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v95,v99));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v32,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v96,v98),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v33,v34),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v97,acc1);
	data[67*stride] = acc1;

	// even samples :: k=68
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v96,v99));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v33,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v97,v98),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v34,acc1);
	data[68*stride] = acc1;

	// odd samples :: k=69
	__m256 v100 = data[100*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v96,v100));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v33,v36),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v97,v99),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v34,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v98,acc1);
	data[69*stride] = acc1;

	// even samples :: k=70
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v97,v100));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v34,v36),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v98,v99),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v35,acc1);
	data[70*stride] = acc1;

	// odd samples :: k=71
	__m256 v101 = data[101*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v97,v101));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v34,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v98,v100),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v35,v36),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v99,acc1);
	data[71*stride] = acc1;

	// even samples :: k=72
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v98,v101));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v35,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v99,v100),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v36,acc1);
	data[72*stride] = acc1;

	// odd samples :: k=73
	__m256 v102 = data[102*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v98,v102));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v35,v38),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v99,v101),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v36,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v100,acc1);
	data[73*stride] = acc1;

	// even samples :: k=74
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v99,v102));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v36,v38),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v100,v101),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v37,acc1);
	data[74*stride] = acc1;

	// odd samples :: k=75
	__m256 v103 = data[103*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v99,v103));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v36,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v100,v102),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v37,v38),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v101,acc1);
	data[75*stride] = acc1;

	// even samples :: k=76
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v100,v103));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v37,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v101,v102),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v38,acc1);
	data[76*stride] = acc1;

	// odd samples :: k=77
	__m256 v104 = data[104*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v100,v104));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v37,v40),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v101,v103),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v38,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v102,acc1);
	data[77*stride] = acc1;

	// even samples :: k=78
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v101,v104));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v38,v40),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v102,v103),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v39,acc1);
	data[78*stride] = acc1;

	// odd samples :: k=79
	__m256 v105 = data[105*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v101,v105));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v38,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v102,v104),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v39,v40),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v103,acc1);
	data[79*stride] = acc1;

	// even samples :: k=80
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v102,v105));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v39,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v103,v104),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v40,acc1);
	data[80*stride] = acc1;

	// odd samples :: k=81
	__m256 v106 = data[106*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v102,v106));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v39,v42),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v103,v105),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v40,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v104,acc1);
	data[81*stride] = acc1;

	// even samples :: k=82
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v103,v106));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v40,v42),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v104,v105),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v41,acc1);
	data[82*stride] = acc1;

	// odd samples :: k=83
	__m256 v107 = data[107*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v103,v107));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v40,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v104,v106),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v41,v42),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v105,acc1);
	data[83*stride] = acc1;

	// even samples :: k=84
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v104,v107));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v41,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v105,v106),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v42,acc1);
	data[84*stride] = acc1;

	// odd samples :: k=85
	__m256 v108 = data[108*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v104,v108));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v41,v44),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v105,v107),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v42,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v106,acc1);
	data[85*stride] = acc1;

	// even samples :: k=86
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v105,v108));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v42,v44),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v106,v107),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v43,acc1);
	data[86*stride] = acc1;

	// odd samples :: k=87
	__m256 v109 = data[109*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v105,v109));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v42,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v106,v108),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v43,v44),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v107,acc1);
	data[87*stride] = acc1;

	// even samples :: k=88
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v106,v109));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v43,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v107,v108),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v44,acc1);
	data[88*stride] = acc1;

	// odd samples :: k=89
	__m256 v110 = data[110*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v106,v110));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v43,v46),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v107,v109),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v44,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v108,acc1);
	data[89*stride] = acc1;

	// even samples :: k=90
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v107,v110));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v44,v46),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v108,v109),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v45,acc1);
	data[90*stride] = acc1;

	// odd samples :: k=91
	__m256 v111 = data[111*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v107,v111));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v44,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v108,v110),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v45,v46),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v109,acc1);
	data[91*stride] = acc1;

	// even samples :: k=92
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v108,v111));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v45,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v109,v110),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v46,acc1);
	data[92*stride] = acc1;

	// odd samples :: k=93
	__m256 v112 = data[112*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v108,v112));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v45,v48),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v109,v111),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v46,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v110,acc1);
	data[93*stride] = acc1;

	// even samples :: k=94
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v109,v112));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v46,v48),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v110,v111),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v47,acc1);
	data[94*stride] = acc1;

	// odd samples :: k=95
	__m256 v113 = data[113*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v109,v113));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v46,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v110,v112),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v47,v48),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v111,acc1);
	data[95*stride] = acc1;

	// even samples :: k=96
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v110,v113));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v47,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v111,v112),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v48,acc1);
	data[96*stride] = acc1;

	// odd samples :: k=97
	__m256 v114 = data[114*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v110,v114));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v47,v50),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v111,v113),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v48,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v112,acc1);
	data[97*stride] = acc1;

	// even samples :: k=98
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v111,v114));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v48,v50),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v112,v113),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v49,acc1);
	data[98*stride] = acc1;

	// odd samples :: k=99
	__m256 v115 = data[115*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v111,v115));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v48,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v112,v114),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v49,v50),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v113,acc1);
	data[99*stride] = acc1;

	// even samples :: k=100
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v112,v115));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v49,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v113,v114),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v50,acc1);
	data[100*stride] = acc1;

	// odd samples :: k=101
	__m256 v116 = data[116*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v112,v116));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v49,v52),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v113,v115),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v50,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v114,acc1);
	data[101*stride] = acc1;

	// even samples :: k=102
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v113,v116));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v50,v52),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v114,v115),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v51,acc1);
	data[102*stride] = acc1;

	// odd samples :: k=103
	__m256 v117 = data[117*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v113,v117));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v50,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v114,v116),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v51,v52),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v115,acc1);
	data[103*stride] = acc1;

	// even samples :: k=104
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v114,v117));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v51,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v115,v116),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v52,acc1);
	data[104*stride] = acc1;

	// odd samples :: k=105
	__m256 v118 = data[118*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v114,v118));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v51,v54),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v115,v117),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v52,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v116,acc1);
	data[105*stride] = acc1;

	// even samples :: k=106
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v115,v118));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v52,v54),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v116,v117),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v53,acc1);
	data[106*stride] = acc1;

	// odd samples :: k=107
	__m256 v119 = data[119*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v115,v119));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v52,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v116,v118),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v53,v54),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v117,acc1);
	data[107*stride] = acc1;

	// even samples :: k=108
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v116,v119));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v53,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v117,v118),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v54,acc1);
	data[108*stride] = acc1;

	// odd samples :: k=109
	__m256 v120 = data[120*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v116,v120));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v53,v56),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v117,v119),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v54,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v118,acc1);
	data[109*stride] = acc1;

	// even samples :: k=110
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v117,v120));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v54,v56),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v118,v119),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v55,acc1);
	data[110*stride] = acc1;

	// odd samples :: k=111
	__m256 v121 = data[121*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v117,v121));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v54,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v118,v120),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v55,v56),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v119,acc1);
	data[111*stride] = acc1;

	// even samples :: k=112
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v118,v121));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v55,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v119,v120),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v56,acc1);
	data[112*stride] = acc1;

	// odd samples :: k=113
	__m256 v122 = data[122*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v118,v122));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v55,v58),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v119,v121),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v56,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v120,acc1);
	data[113*stride] = acc1;

	// even samples :: k=114
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v119,v122));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v56,v58),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v120,v121),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v57,acc1);
	data[114*stride] = acc1;

	// odd samples :: k=115
	__m256 v123 = data[123*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v119,v123));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v56,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v120,v122),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v57,v58),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v121,acc1);
	data[115*stride] = acc1;

	// even samples :: k=116
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v120,v123));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v57,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v121,v122),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v58,acc1);
	data[116*stride] = acc1;

	// odd samples :: k=117
	__m256 v124 = data[124*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v120,v124));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v57,v60),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v121,v123),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v58,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v122,acc1);
	data[117*stride] = acc1;

	// even samples :: k=118
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v121,v124));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v58,v60),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v122,v123),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v59,acc1);
	data[118*stride] = acc1;

	// odd samples :: k=119
	__m256 v125 = data[125*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v121,v125));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v58,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v122,v124),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v59,v60),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v123,acc1);
	data[119*stride] = acc1;

	// even samples :: k=120
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v122,v125));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v59,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v123,v124),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v60,acc1);
	data[120*stride] = acc1;

	// odd samples :: k=121
	__m256 v126 = data[126*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v122,v126));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v59,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v123,v125),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v60,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v124,acc1);
	data[121*stride] = acc1;

	// even samples :: k=122
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v123,v126));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v60,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v124,v125),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v61,acc1);
	data[122*stride] = acc1;

	// odd samples :: k=123
	__m256 v127 = data[127*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v123,v127));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v60,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v124,v126),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v61,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v125,acc1);
	data[123*stride] = acc1;

	// even samples :: k=124
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v124,v127));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v61,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v125,v126),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v62,acc1);
	data[124*stride] = acc1;

	// odd samples :: k=125
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v124,v126));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v61,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v125,v127),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v62,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v126,acc1);
	data[125*stride] = acc1;

	// even samples :: k=126
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v125,v126));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v62,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v126,v127),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v63,acc1);
	data[126*stride] = acc1;

	// odd samples :: k=127
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v125,v125));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v62,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v126,v126),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v63,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v127,acc1);
	data[127*stride] = acc1;
}

static void _Us79_AVX_256(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v129 = data[129*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v129,v129));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v1,v1),acc1);
	__m256 v128 = data[128*stride];
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v128,v128),acc1);
	__m256 v0 = data[0*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl0,v0,acc1);
	data[0*stride] = acc1;

	// odd samples :: k=1
	__m256 v130 = data[130*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v129,v130));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v1,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v128,v129),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v0,v1),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v128,acc1);
	data[1*stride] = acc1;

	// even samples :: k=2
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v128,v130));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v0,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v128,v129),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v1,acc1);
	data[2*stride] = acc1;

	// odd samples :: k=3
	__m256 v131 = data[131*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v128,v131));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v0,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v128,v130),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v1,v2),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v129,acc1);
	data[3*stride] = acc1;

	// even samples :: k=4
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v128,v131));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v1,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v129,v130),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v2,acc1);
	__m256 v4 = data[4*stride];
	data[4*stride] = acc1;

	// odd samples :: k=5
	__m256 v132 = data[132*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v128,v132));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v1,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v129,v131),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v2,v3),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v130,acc1);
	__m256 v5 = data[5*stride];
	data[5*stride] = acc1;

	// even samples :: k=6
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v129,v132));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v2,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v130,v131),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v3,acc1);
	__m256 v6 = data[6*stride];
	data[6*stride] = acc1;

	// odd samples :: k=7
	__m256 v133 = data[133*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v129,v133));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v2,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v130,v132),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v3,v4),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v131,acc1);
	__m256 v7 = data[7*stride];
	data[7*stride] = acc1;

	// even samples :: k=8
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v130,v133));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v3,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v131,v132),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v4,acc1);
	__m256 v8 = data[8*stride];
	data[8*stride] = acc1;

	// odd samples :: k=9
	__m256 v134 = data[134*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v130,v134));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v3,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v131,v133),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v4,v5),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v132,acc1);
	__m256 v9 = data[9*stride];
	data[9*stride] = acc1;

	// even samples :: k=10
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v131,v134));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v4,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v132,v133),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v5,acc1);
	__m256 v10 = data[10*stride];
	data[10*stride] = acc1;

	// odd samples :: k=11
	__m256 v135 = data[135*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v131,v135));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v4,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v132,v134),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v5,v6),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v133,acc1);
	__m256 v11 = data[11*stride];
	data[11*stride] = acc1;

	// even samples :: k=12
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v132,v135));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v5,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v133,v134),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v6,acc1);
	__m256 v12 = data[12*stride];
	data[12*stride] = acc1;

	// odd samples :: k=13
	__m256 v136 = data[136*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v132,v136));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v5,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v133,v135),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v6,v7),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v134,acc1);
	__m256 v13 = data[13*stride];
	data[13*stride] = acc1;

	// even samples :: k=14
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v133,v136));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v6,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v134,v135),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v7,acc1);
	__m256 v14 = data[14*stride];
	data[14*stride] = acc1;

	// odd samples :: k=15
	__m256 v137 = data[137*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v133,v137));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v6,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v134,v136),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v7,v8),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v135,acc1);
	__m256 v15 = data[15*stride];
	data[15*stride] = acc1;

	// even samples :: k=16
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v134,v137));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v7,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v135,v136),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v8,acc1);
	__m256 v16 = data[16*stride];
	data[16*stride] = acc1;

	// odd samples :: k=17
	__m256 v138 = data[138*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v134,v138));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v7,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v135,v137),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v8,v9),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v136,acc1);
	__m256 v17 = data[17*stride];
	data[17*stride] = acc1;

	// even samples :: k=18
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v135,v138));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v8,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v136,v137),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v9,acc1);
	__m256 v18 = data[18*stride];
	data[18*stride] = acc1;

	// odd samples :: k=19
	__m256 v139 = data[139*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v135,v139));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v8,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v136,v138),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v9,v10),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v137,acc1);
	__m256 v19 = data[19*stride];
	data[19*stride] = acc1;

	// even samples :: k=20
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v136,v139));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v9,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v137,v138),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v10,acc1);
	__m256 v20 = data[20*stride];
	data[20*stride] = acc1;

	// odd samples :: k=21
	__m256 v140 = data[140*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v136,v140));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v9,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v137,v139),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v10,v11),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v138,acc1);
	__m256 v21 = data[21*stride];
	data[21*stride] = acc1;

	// even samples :: k=22
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v137,v140));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v10,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v138,v139),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v11,acc1);
	__m256 v22 = data[22*stride];
	data[22*stride] = acc1;

	// odd samples :: k=23
	__m256 v141 = data[141*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v137,v141));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v10,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v138,v140),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v11,v12),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v139,acc1);
	__m256 v23 = data[23*stride];
	data[23*stride] = acc1;

	// even samples :: k=24
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v138,v141));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v11,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v139,v140),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v12,acc1);
	__m256 v24 = data[24*stride];
	data[24*stride] = acc1;

	// odd samples :: k=25
	__m256 v142 = data[142*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v138,v142));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v11,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v139,v141),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v12,v13),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v140,acc1);
	__m256 v25 = data[25*stride];
	data[25*stride] = acc1;

	// even samples :: k=26
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v139,v142));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v12,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v140,v141),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v13,acc1);
	__m256 v26 = data[26*stride];
	data[26*stride] = acc1;

	// odd samples :: k=27
	__m256 v143 = data[143*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v139,v143));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v12,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v140,v142),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v13,v14),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v141,acc1);
	__m256 v27 = data[27*stride];
	data[27*stride] = acc1;

	// even samples :: k=28
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v140,v143));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v13,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v141,v142),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v14,acc1);
	__m256 v28 = data[28*stride];
	data[28*stride] = acc1;

	// odd samples :: k=29
	__m256 v144 = data[144*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v140,v144));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v13,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v141,v143),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v14,v15),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v142,acc1);
	__m256 v29 = data[29*stride];
	data[29*stride] = acc1;

	// even samples :: k=30
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v141,v144));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v14,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v142,v143),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v15,acc1);
	__m256 v30 = data[30*stride];
	data[30*stride] = acc1;

	// odd samples :: k=31
	__m256 v145 = data[145*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v141,v145));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v14,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v142,v144),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v15,v16),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v143,acc1);
	__m256 v31 = data[31*stride];
	data[31*stride] = acc1;

	// even samples :: k=32
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v142,v145));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v15,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v143,v144),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v16,acc1);
	__m256 v32 = data[32*stride];
	data[32*stride] = acc1;

	// odd samples :: k=33
	__m256 v146 = data[146*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v142,v146));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v15,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v143,v145),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v16,v17),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v144,acc1);
	__m256 v33 = data[33*stride];
	data[33*stride] = acc1;

	// even samples :: k=34
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v143,v146));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v16,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v144,v145),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v17,acc1);
	__m256 v34 = data[34*stride];
	data[34*stride] = acc1;

	// odd samples :: k=35
	__m256 v147 = data[147*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v143,v147));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v16,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v144,v146),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v17,v18),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v145,acc1);
	__m256 v35 = data[35*stride];
	data[35*stride] = acc1;

	// even samples :: k=36
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v144,v147));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v17,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v145,v146),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v18,acc1);
	__m256 v36 = data[36*stride];
	data[36*stride] = acc1;

	// odd samples :: k=37
	__m256 v148 = data[148*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v144,v148));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v17,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v145,v147),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v18,v19),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v146,acc1);
	__m256 v37 = data[37*stride];
	data[37*stride] = acc1;

	// even samples :: k=38
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v145,v148));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v18,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v146,v147),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v19,acc1);
	__m256 v38 = data[38*stride];
	data[38*stride] = acc1;

	// odd samples :: k=39
	__m256 v149 = data[149*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v145,v149));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v18,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v146,v148),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v19,v20),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v147,acc1);
	__m256 v39 = data[39*stride];
	data[39*stride] = acc1;

	// even samples :: k=40
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v146,v149));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v19,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v147,v148),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v20,acc1);
	__m256 v40 = data[40*stride];
	data[40*stride] = acc1;

	// odd samples :: k=41
	__m256 v150 = data[150*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v146,v150));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v19,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v147,v149),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v20,v21),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v148,acc1);
	__m256 v41 = data[41*stride];
	data[41*stride] = acc1;

	// even samples :: k=42
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v147,v150));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v20,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v148,v149),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v21,acc1);
	__m256 v42 = data[42*stride];
	data[42*stride] = acc1;

	// odd samples :: k=43
	__m256 v151 = data[151*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v147,v151));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v20,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v148,v150),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v21,v22),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v149,acc1);
	__m256 v43 = data[43*stride];
	data[43*stride] = acc1;

	// even samples :: k=44
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v148,v151));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v21,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v149,v150),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v22,acc1);
	__m256 v44 = data[44*stride];
	data[44*stride] = acc1;

	// odd samples :: k=45
	__m256 v152 = data[152*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v148,v152));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v21,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v149,v151),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v22,v23),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v150,acc1);
	__m256 v45 = data[45*stride];
	data[45*stride] = acc1;

	// even samples :: k=46
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v149,v152));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v22,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v150,v151),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v23,acc1);
	__m256 v46 = data[46*stride];
	data[46*stride] = acc1;

	// odd samples :: k=47
	__m256 v153 = data[153*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v149,v153));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v22,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v150,v152),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v23,v24),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v151,acc1);
	__m256 v47 = data[47*stride];
	data[47*stride] = acc1;

	// even samples :: k=48
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v150,v153));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v23,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v151,v152),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v24,acc1);
	__m256 v48 = data[48*stride];
	data[48*stride] = acc1;

	// odd samples :: k=49
	__m256 v154 = data[154*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v150,v154));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v23,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v151,v153),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v24,v25),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v152,acc1);
	__m256 v49 = data[49*stride];
	data[49*stride] = acc1;

	// even samples :: k=50
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v151,v154));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v24,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v152,v153),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v25,acc1);
	__m256 v50 = data[50*stride];
	data[50*stride] = acc1;

	// odd samples :: k=51
	__m256 v155 = data[155*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v151,v155));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v24,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v152,v154),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v25,v26),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v153,acc1);
	__m256 v51 = data[51*stride];
	data[51*stride] = acc1;

	// even samples :: k=52
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v152,v155));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v25,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v153,v154),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v26,acc1);
	__m256 v52 = data[52*stride];
	data[52*stride] = acc1;

	// odd samples :: k=53
	__m256 v156 = data[156*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v152,v156));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v25,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v153,v155),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v26,v27),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v154,acc1);
	__m256 v53 = data[53*stride];
	data[53*stride] = acc1;

	// even samples :: k=54
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v153,v156));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v26,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v154,v155),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v27,acc1);
	__m256 v54 = data[54*stride];
	data[54*stride] = acc1;

	// odd samples :: k=55
	__m256 v157 = data[157*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v153,v157));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v26,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v154,v156),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v27,v28),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v155,acc1);
	__m256 v55 = data[55*stride];
	data[55*stride] = acc1;

	// even samples :: k=56
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v154,v157));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v27,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v155,v156),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v28,acc1);
	__m256 v56 = data[56*stride];
	data[56*stride] = acc1;

	// odd samples :: k=57
	__m256 v158 = data[158*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v154,v158));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v27,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v155,v157),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v28,v29),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v156,acc1);
	__m256 v57 = data[57*stride];
	data[57*stride] = acc1;

	// even samples :: k=58
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v155,v158));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v28,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v156,v157),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v29,acc1);
	__m256 v58 = data[58*stride];
	data[58*stride] = acc1;

	// odd samples :: k=59
	__m256 v159 = data[159*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v155,v159));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v28,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v156,v158),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v29,v30),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v157,acc1);
	__m256 v59 = data[59*stride];
	data[59*stride] = acc1;

	// even samples :: k=60
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v156,v159));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v29,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v157,v158),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v30,acc1);
	__m256 v60 = data[60*stride];
	data[60*stride] = acc1;

	// odd samples :: k=61
	__m256 v160 = data[160*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v156,v160));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v29,v32),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v157,v159),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v30,v31),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v158,acc1);
	__m256 v61 = data[61*stride];
	data[61*stride] = acc1;

	// even samples :: k=62
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v157,v160));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v30,v32),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v158,v159),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v31,acc1);
	__m256 v62 = data[62*stride];
	data[62*stride] = acc1;

	// odd samples :: k=63
	__m256 v161 = data[161*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v157,v161));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v30,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v158,v160),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v31,v32),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v159,acc1);
	__m256 v63 = data[63*stride];
	data[63*stride] = acc1;

	// even samples :: k=64
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v158,v161));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v31,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v159,v160),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v32,acc1);
	__m256 v64 = data[64*stride];
	data[64*stride] = acc1;

	// odd samples :: k=65
	__m256 v162 = data[162*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v158,v162));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v31,v34),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v159,v161),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v32,v33),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v160,acc1);
	__m256 v65 = data[65*stride];
	data[65*stride] = acc1;

	// even samples :: k=66
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v159,v162));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v32,v34),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v160,v161),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v33,acc1);
	__m256 v66 = data[66*stride];
	data[66*stride] = acc1;

	// odd samples :: k=67
	__m256 v163 = data[163*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v159,v163));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v32,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v160,v162),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v33,v34),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v161,acc1);
	__m256 v67 = data[67*stride];
	data[67*stride] = acc1;

	// even samples :: k=68
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v160,v163));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v33,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v161,v162),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v34,acc1);
	__m256 v68 = data[68*stride];
	data[68*stride] = acc1;

	// odd samples :: k=69
	__m256 v164 = data[164*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v160,v164));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v33,v36),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v161,v163),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v34,v35),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v162,acc1);
	__m256 v69 = data[69*stride];
	data[69*stride] = acc1;

	// even samples :: k=70
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v161,v164));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v34,v36),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v162,v163),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v35,acc1);
	__m256 v70 = data[70*stride];
	data[70*stride] = acc1;

	// odd samples :: k=71
	__m256 v165 = data[165*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v161,v165));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v34,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v162,v164),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v35,v36),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v163,acc1);
	__m256 v71 = data[71*stride];
	data[71*stride] = acc1;

	// even samples :: k=72
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v162,v165));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v35,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v163,v164),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v36,acc1);
	__m256 v72 = data[72*stride];
	data[72*stride] = acc1;

	// odd samples :: k=73
	__m256 v166 = data[166*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v162,v166));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v35,v38),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v163,v165),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v36,v37),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v164,acc1);
	__m256 v73 = data[73*stride];
	data[73*stride] = acc1;

	// even samples :: k=74
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v163,v166));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v36,v38),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v164,v165),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v37,acc1);
	__m256 v74 = data[74*stride];
	data[74*stride] = acc1;

	// odd samples :: k=75
	__m256 v167 = data[167*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v163,v167));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v36,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v164,v166),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v37,v38),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v165,acc1);
	__m256 v75 = data[75*stride];
	data[75*stride] = acc1;

	// even samples :: k=76
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v164,v167));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v37,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v165,v166),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v38,acc1);
	__m256 v76 = data[76*stride];
	data[76*stride] = acc1;

	// odd samples :: k=77
	__m256 v168 = data[168*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v164,v168));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v37,v40),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v165,v167),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v38,v39),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v166,acc1);
	__m256 v77 = data[77*stride];
	data[77*stride] = acc1;

	// even samples :: k=78
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v165,v168));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v38,v40),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v166,v167),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v39,acc1);
	__m256 v78 = data[78*stride];
	data[78*stride] = acc1;

	// odd samples :: k=79
	__m256 v169 = data[169*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v165,v169));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v38,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v166,v168),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v39,v40),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v167,acc1);
	__m256 v79 = data[79*stride];
	data[79*stride] = acc1;

	// even samples :: k=80
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v166,v169));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v39,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v167,v168),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v40,acc1);
	__m256 v80 = data[80*stride];
	data[80*stride] = acc1;

	// odd samples :: k=81
	__m256 v170 = data[170*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v166,v170));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v39,v42),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v167,v169),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v40,v41),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v168,acc1);
	__m256 v81 = data[81*stride];
	data[81*stride] = acc1;

	// even samples :: k=82
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v167,v170));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v40,v42),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v168,v169),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v41,acc1);
	__m256 v82 = data[82*stride];
	data[82*stride] = acc1;

	// odd samples :: k=83
	__m256 v171 = data[171*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v167,v171));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v40,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v168,v170),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v41,v42),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v169,acc1);
	__m256 v83 = data[83*stride];
	data[83*stride] = acc1;

	// even samples :: k=84
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v168,v171));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v41,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v169,v170),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v42,acc1);
	__m256 v84 = data[84*stride];
	data[84*stride] = acc1;

	// odd samples :: k=85
	__m256 v172 = data[172*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v168,v172));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v41,v44),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v169,v171),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v42,v43),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v170,acc1);
	__m256 v85 = data[85*stride];
	data[85*stride] = acc1;

	// even samples :: k=86
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v169,v172));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v42,v44),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v170,v171),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v43,acc1);
	__m256 v86 = data[86*stride];
	data[86*stride] = acc1;

	// odd samples :: k=87
	__m256 v173 = data[173*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v169,v173));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v42,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v170,v172),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v43,v44),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v171,acc1);
	__m256 v87 = data[87*stride];
	data[87*stride] = acc1;

	// even samples :: k=88
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v170,v173));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v43,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v171,v172),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v44,acc1);
	__m256 v88 = data[88*stride];
	data[88*stride] = acc1;

	// odd samples :: k=89
	__m256 v174 = data[174*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v170,v174));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v43,v46),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v171,v173),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v44,v45),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v172,acc1);
	__m256 v89 = data[89*stride];
	data[89*stride] = acc1;

	// even samples :: k=90
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v171,v174));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v44,v46),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v172,v173),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v45,acc1);
	__m256 v90 = data[90*stride];
	data[90*stride] = acc1;

	// odd samples :: k=91
	__m256 v175 = data[175*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v171,v175));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v44,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v172,v174),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v45,v46),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v173,acc1);
	__m256 v91 = data[91*stride];
	data[91*stride] = acc1;

	// even samples :: k=92
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v172,v175));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v45,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v173,v174),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v46,acc1);
	__m256 v92 = data[92*stride];
	data[92*stride] = acc1;

	// odd samples :: k=93
	__m256 v176 = data[176*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v172,v176));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v45,v48),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v173,v175),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v46,v47),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v174,acc1);
	__m256 v93 = data[93*stride];
	data[93*stride] = acc1;

	// even samples :: k=94
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v173,v176));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v46,v48),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v174,v175),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v47,acc1);
	__m256 v94 = data[94*stride];
	data[94*stride] = acc1;

	// odd samples :: k=95
	__m256 v177 = data[177*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v173,v177));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v46,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v174,v176),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v47,v48),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v175,acc1);
	__m256 v95 = data[95*stride];
	data[95*stride] = acc1;

	// even samples :: k=96
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v174,v177));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v47,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v175,v176),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v48,acc1);
	__m256 v96 = data[96*stride];
	data[96*stride] = acc1;

	// odd samples :: k=97
	__m256 v178 = data[178*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v174,v178));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v47,v50),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v175,v177),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v48,v49),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v176,acc1);
	__m256 v97 = data[97*stride];
	data[97*stride] = acc1;

	// even samples :: k=98
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v175,v178));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v48,v50),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v176,v177),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v49,acc1);
	__m256 v98 = data[98*stride];
	data[98*stride] = acc1;

	// odd samples :: k=99
	__m256 v179 = data[179*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v175,v179));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v48,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v176,v178),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v49,v50),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v177,acc1);
	__m256 v99 = data[99*stride];
	data[99*stride] = acc1;

	// even samples :: k=100
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v176,v179));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v49,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v177,v178),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v50,acc1);
	__m256 v100 = data[100*stride];
	data[100*stride] = acc1;

	// odd samples :: k=101
	__m256 v180 = data[180*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v176,v180));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v49,v52),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v177,v179),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v50,v51),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v178,acc1);
	__m256 v101 = data[101*stride];
	data[101*stride] = acc1;

	// even samples :: k=102
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v177,v180));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v50,v52),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v178,v179),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v51,acc1);
	__m256 v102 = data[102*stride];
	data[102*stride] = acc1;

	// odd samples :: k=103
	__m256 v181 = data[181*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v177,v181));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v50,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v178,v180),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v51,v52),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v179,acc1);
	__m256 v103 = data[103*stride];
	data[103*stride] = acc1;

	// even samples :: k=104
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v178,v181));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v51,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v179,v180),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v52,acc1);
	__m256 v104 = data[104*stride];
	data[104*stride] = acc1;

	// odd samples :: k=105
	__m256 v182 = data[182*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v178,v182));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v51,v54),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v179,v181),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v52,v53),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v180,acc1);
	__m256 v105 = data[105*stride];
	data[105*stride] = acc1;

	// even samples :: k=106
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v179,v182));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v52,v54),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v180,v181),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v53,acc1);
	__m256 v106 = data[106*stride];
	data[106*stride] = acc1;

	// odd samples :: k=107
	__m256 v183 = data[183*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v179,v183));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v52,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v180,v182),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v53,v54),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v181,acc1);
	__m256 v107 = data[107*stride];
	data[107*stride] = acc1;

	// even samples :: k=108
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v180,v183));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v53,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v181,v182),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v54,acc1);
	__m256 v108 = data[108*stride];
	data[108*stride] = acc1;

	// odd samples :: k=109
	__m256 v184 = data[184*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v180,v184));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v53,v56),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v181,v183),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v54,v55),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v182,acc1);
	__m256 v109 = data[109*stride];
	data[109*stride] = acc1;

	// even samples :: k=110
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v181,v184));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v54,v56),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v182,v183),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v55,acc1);
	__m256 v110 = data[110*stride];
	data[110*stride] = acc1;

	// odd samples :: k=111
	__m256 v185 = data[185*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v181,v185));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v54,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v182,v184),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v55,v56),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v183,acc1);
	__m256 v111 = data[111*stride];
	data[111*stride] = acc1;

	// even samples :: k=112
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v182,v185));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v55,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v183,v184),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v56,acc1);
	__m256 v112 = data[112*stride];
	data[112*stride] = acc1;

	// odd samples :: k=113
	__m256 v186 = data[186*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v182,v186));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v55,v58),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v183,v185),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v56,v57),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v184,acc1);
	__m256 v113 = data[113*stride];
	data[113*stride] = acc1;

	// even samples :: k=114
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v183,v186));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v56,v58),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v184,v185),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v57,acc1);
	__m256 v114 = data[114*stride];
	data[114*stride] = acc1;

	// odd samples :: k=115
	__m256 v187 = data[187*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v183,v187));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v56,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v184,v186),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v57,v58),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v185,acc1);
	__m256 v115 = data[115*stride];
	data[115*stride] = acc1;

	// even samples :: k=116
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v184,v187));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v57,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v185,v186),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v58,acc1);
	__m256 v116 = data[116*stride];
	data[116*stride] = acc1;

	// odd samples :: k=117
	__m256 v188 = data[188*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v184,v188));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v57,v60),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v185,v187),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v58,v59),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v186,acc1);
	__m256 v117 = data[117*stride];
	data[117*stride] = acc1;

	// even samples :: k=118
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v185,v188));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v58,v60),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v186,v187),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v59,acc1);
	__m256 v118 = data[118*stride];
	data[118*stride] = acc1;

	// odd samples :: k=119
	__m256 v189 = data[189*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v185,v189));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v58,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v186,v188),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v59,v60),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v187,acc1);
	__m256 v119 = data[119*stride];
	data[119*stride] = acc1;

	// even samples :: k=120
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v186,v189));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v59,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v187,v188),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v60,acc1);
	__m256 v120 = data[120*stride];
	data[120*stride] = acc1;

	// odd samples :: k=121
	__m256 v190 = data[190*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v186,v190));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v59,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v187,v189),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v60,v61),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v188,acc1);
	__m256 v121 = data[121*stride];
	data[121*stride] = acc1;

	// even samples :: k=122
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v187,v190));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v60,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v188,v189),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v61,acc1);
	__m256 v122 = data[122*stride];
	data[122*stride] = acc1;

	// odd samples :: k=123
	__m256 v191 = data[191*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v187,v191));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v60,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v188,v190),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v61,v62),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v189,acc1);
	__m256 v123 = data[123*stride];
	data[123*stride] = acc1;

	// even samples :: k=124
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v188,v191));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v61,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v189,v190),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v62,acc1);
	__m256 v124 = data[124*stride];
	data[124*stride] = acc1;

	// odd samples :: k=125
	__m256 v192 = data[192*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v188,v192));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v61,v64),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v189,v191),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v62,v63),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v190,acc1);
	__m256 v125 = data[125*stride];
	data[125*stride] = acc1;

	// even samples :: k=126
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v189,v192));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v62,v64),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v190,v191),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v63,acc1);
	__m256 v126 = data[126*stride];
	data[126*stride] = acc1;

	// odd samples :: k=127
	__m256 v193 = data[193*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v189,v193));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v62,v65),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v190,v192),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v63,v64),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v191,acc1);
	__m256 v127 = data[127*stride];
	data[127*stride] = acc1;

	// even samples :: k=128
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v190,v193));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v63,v65),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v191,v192),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v64,acc1);
	data[128*stride] = acc1;

	// odd samples :: k=129
	__m256 v194 = data[194*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v190,v194));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v63,v66),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v191,v193),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v64,v65),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v192,acc1);
	data[129*stride] = acc1;

	// even samples :: k=130
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v191,v194));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v64,v66),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v192,v193),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v65,acc1);
	data[130*stride] = acc1;

	// odd samples :: k=131
	__m256 v195 = data[195*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v191,v195));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v64,v67),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v192,v194),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v65,v66),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v193,acc1);
	data[131*stride] = acc1;

	// even samples :: k=132
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v192,v195));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v65,v67),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v193,v194),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v66,acc1);
	data[132*stride] = acc1;

	// odd samples :: k=133
	__m256 v196 = data[196*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v192,v196));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v65,v68),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v193,v195),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v66,v67),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v194,acc1);
	data[133*stride] = acc1;

	// even samples :: k=134
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v193,v196));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v66,v68),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v194,v195),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v67,acc1);
	data[134*stride] = acc1;

	// odd samples :: k=135
	__m256 v197 = data[197*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v193,v197));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v66,v69),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v194,v196),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v67,v68),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v195,acc1);
	data[135*stride] = acc1;

	// even samples :: k=136
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v194,v197));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v67,v69),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v195,v196),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v68,acc1);
	data[136*stride] = acc1;

	// odd samples :: k=137
	__m256 v198 = data[198*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v194,v198));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v67,v70),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v195,v197),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v68,v69),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v196,acc1);
	data[137*stride] = acc1;

	// even samples :: k=138
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v195,v198));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v68,v70),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v196,v197),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v69,acc1);
	data[138*stride] = acc1;

	// odd samples :: k=139
	__m256 v199 = data[199*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v195,v199));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v68,v71),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v196,v198),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v69,v70),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v197,acc1);
	data[139*stride] = acc1;

	// even samples :: k=140
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v196,v199));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v69,v71),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v197,v198),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v70,acc1);
	data[140*stride] = acc1;

	// odd samples :: k=141
	__m256 v200 = data[200*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v196,v200));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v69,v72),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v197,v199),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v70,v71),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v198,acc1);
	data[141*stride] = acc1;

	// even samples :: k=142
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v197,v200));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v70,v72),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v198,v199),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v71,acc1);
	data[142*stride] = acc1;

	// odd samples :: k=143
	__m256 v201 = data[201*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v197,v201));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v70,v73),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v198,v200),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v71,v72),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v199,acc1);
	data[143*stride] = acc1;

	// even samples :: k=144
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v198,v201));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v71,v73),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v199,v200),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v72,acc1);
	data[144*stride] = acc1;

	// odd samples :: k=145
	__m256 v202 = data[202*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v198,v202));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v71,v74),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v199,v201),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v72,v73),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v200,acc1);
	data[145*stride] = acc1;

	// even samples :: k=146
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v199,v202));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v72,v74),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v200,v201),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v73,acc1);
	data[146*stride] = acc1;

	// odd samples :: k=147
	__m256 v203 = data[203*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v199,v203));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v72,v75),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v200,v202),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v73,v74),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v201,acc1);
	data[147*stride] = acc1;

	// even samples :: k=148
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v200,v203));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v73,v75),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v201,v202),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v74,acc1);
	data[148*stride] = acc1;

	// odd samples :: k=149
	__m256 v204 = data[204*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v200,v204));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v73,v76),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v201,v203),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v74,v75),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v202,acc1);
	data[149*stride] = acc1;

	// even samples :: k=150
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v201,v204));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v74,v76),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v202,v203),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v75,acc1);
	data[150*stride] = acc1;

	// odd samples :: k=151
	__m256 v205 = data[205*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v201,v205));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v74,v77),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v202,v204),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v75,v76),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v203,acc1);
	data[151*stride] = acc1;

	// even samples :: k=152
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v202,v205));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v75,v77),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v203,v204),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v76,acc1);
	data[152*stride] = acc1;

	// odd samples :: k=153
	__m256 v206 = data[206*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v202,v206));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v75,v78),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v203,v205),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v76,v77),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v204,acc1);
	data[153*stride] = acc1;

	// even samples :: k=154
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v203,v206));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v76,v78),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v204,v205),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v77,acc1);
	data[154*stride] = acc1;

	// odd samples :: k=155
	__m256 v207 = data[207*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v203,v207));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v76,v79),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v204,v206),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v77,v78),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v205,acc1);
	data[155*stride] = acc1;

	// even samples :: k=156
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v204,v207));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v77,v79),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v205,v206),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v78,acc1);
	data[156*stride] = acc1;

	// odd samples :: k=157
	__m256 v208 = data[208*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v204,v208));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v77,v80),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v205,v207),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v78,v79),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v206,acc1);
	data[157*stride] = acc1;

	// even samples :: k=158
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v205,v208));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v78,v80),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v206,v207),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v79,acc1);
	data[158*stride] = acc1;

	// odd samples :: k=159
	__m256 v209 = data[209*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v205,v209));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v78,v81),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v206,v208),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v79,v80),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v207,acc1);
	data[159*stride] = acc1;

	// even samples :: k=160
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v206,v209));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v79,v81),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v207,v208),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v80,acc1);
	data[160*stride] = acc1;

	// odd samples :: k=161
	__m256 v210 = data[210*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v206,v210));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v79,v82),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v207,v209),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v80,v81),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v208,acc1);
	data[161*stride] = acc1;

	// even samples :: k=162
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v207,v210));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v80,v82),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v208,v209),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v81,acc1);
	data[162*stride] = acc1;

	// odd samples :: k=163
	__m256 v211 = data[211*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v207,v211));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v80,v83),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v208,v210),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v81,v82),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v209,acc1);
	data[163*stride] = acc1;

	// even samples :: k=164
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v208,v211));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v81,v83),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v209,v210),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v82,acc1);
	data[164*stride] = acc1;

	// odd samples :: k=165
	__m256 v212 = data[212*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v208,v212));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v81,v84),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v209,v211),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v82,v83),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v210,acc1);
	data[165*stride] = acc1;

	// even samples :: k=166
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v209,v212));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v82,v84),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v210,v211),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v83,acc1);
	data[166*stride] = acc1;

	// odd samples :: k=167
	__m256 v213 = data[213*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v209,v213));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v82,v85),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v210,v212),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v83,v84),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v211,acc1);
	data[167*stride] = acc1;

	// even samples :: k=168
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v210,v213));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v83,v85),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v211,v212),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v84,acc1);
	data[168*stride] = acc1;

	// odd samples :: k=169
	__m256 v214 = data[214*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v210,v214));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v83,v86),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v211,v213),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v84,v85),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v212,acc1);
	data[169*stride] = acc1;

	// even samples :: k=170
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v211,v214));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v84,v86),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v212,v213),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v85,acc1);
	data[170*stride] = acc1;

	// odd samples :: k=171
	__m256 v215 = data[215*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v211,v215));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v84,v87),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v212,v214),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v85,v86),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v213,acc1);
	data[171*stride] = acc1;

	// even samples :: k=172
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v212,v215));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v85,v87),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v213,v214),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v86,acc1);
	data[172*stride] = acc1;

	// odd samples :: k=173
	__m256 v216 = data[216*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v212,v216));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v85,v88),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v213,v215),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v86,v87),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v214,acc1);
	data[173*stride] = acc1;

	// even samples :: k=174
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v213,v216));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v86,v88),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v214,v215),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v87,acc1);
	data[174*stride] = acc1;

	// odd samples :: k=175
	__m256 v217 = data[217*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v213,v217));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v86,v89),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v214,v216),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v87,v88),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v215,acc1);
	data[175*stride] = acc1;

	// even samples :: k=176
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v214,v217));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v87,v89),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v215,v216),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v88,acc1);
	data[176*stride] = acc1;

	// odd samples :: k=177
	__m256 v218 = data[218*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v214,v218));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v87,v90),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v215,v217),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v88,v89),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v216,acc1);
	data[177*stride] = acc1;

	// even samples :: k=178
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v215,v218));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v88,v90),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v216,v217),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v89,acc1);
	data[178*stride] = acc1;

	// odd samples :: k=179
	__m256 v219 = data[219*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v215,v219));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v88,v91),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v216,v218),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v89,v90),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v217,acc1);
	data[179*stride] = acc1;

	// even samples :: k=180
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v216,v219));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v89,v91),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v217,v218),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v90,acc1);
	data[180*stride] = acc1;

	// odd samples :: k=181
	__m256 v220 = data[220*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v216,v220));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v89,v92),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v217,v219),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v90,v91),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v218,acc1);
	data[181*stride] = acc1;

	// even samples :: k=182
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v217,v220));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v90,v92),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v218,v219),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v91,acc1);
	data[182*stride] = acc1;

	// odd samples :: k=183
	__m256 v221 = data[221*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v217,v221));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v90,v93),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v218,v220),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v91,v92),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v219,acc1);
	data[183*stride] = acc1;

	// even samples :: k=184
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v218,v221));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v91,v93),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v219,v220),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v92,acc1);
	data[184*stride] = acc1;

	// odd samples :: k=185
	__m256 v222 = data[222*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v218,v222));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v91,v94),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v219,v221),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v92,v93),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v220,acc1);
	data[185*stride] = acc1;

	// even samples :: k=186
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v219,v222));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v92,v94),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v220,v221),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v93,acc1);
	data[186*stride] = acc1;

	// odd samples :: k=187
	__m256 v223 = data[223*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v219,v223));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v92,v95),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v220,v222),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v93,v94),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v221,acc1);
	data[187*stride] = acc1;

	// even samples :: k=188
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v220,v223));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v93,v95),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v221,v222),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v94,acc1);
	data[188*stride] = acc1;

	// odd samples :: k=189
	__m256 v224 = data[224*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v220,v224));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v93,v96),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v221,v223),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v94,v95),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v222,acc1);
	data[189*stride] = acc1;

	// even samples :: k=190
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v221,v224));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v94,v96),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v222,v223),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v95,acc1);
	data[190*stride] = acc1;

	// odd samples :: k=191
	__m256 v225 = data[225*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v221,v225));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v94,v97),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v222,v224),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v95,v96),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v223,acc1);
	data[191*stride] = acc1;

	// even samples :: k=192
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v222,v225));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v95,v97),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v223,v224),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v96,acc1);
	data[192*stride] = acc1;

	// odd samples :: k=193
	__m256 v226 = data[226*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v222,v226));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v95,v98),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v223,v225),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v96,v97),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v224,acc1);
	data[193*stride] = acc1;

	// even samples :: k=194
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v223,v226));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v96,v98),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v224,v225),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v97,acc1);
	data[194*stride] = acc1;

	// odd samples :: k=195
	__m256 v227 = data[227*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v223,v227));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v96,v99),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v224,v226),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v97,v98),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v225,acc1);
	data[195*stride] = acc1;

	// even samples :: k=196
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v224,v227));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v97,v99),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v225,v226),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v98,acc1);
	data[196*stride] = acc1;

	// odd samples :: k=197
	__m256 v228 = data[228*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v224,v228));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v97,v100),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v225,v227),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v98,v99),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v226,acc1);
	data[197*stride] = acc1;

	// even samples :: k=198
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v225,v228));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v98,v100),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v226,v227),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v99,acc1);
	data[198*stride] = acc1;

	// odd samples :: k=199
	__m256 v229 = data[229*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v225,v229));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v98,v101),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v226,v228),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v99,v100),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v227,acc1);
	data[199*stride] = acc1;

	// even samples :: k=200
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v226,v229));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v99,v101),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v227,v228),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v100,acc1);
	data[200*stride] = acc1;

	// odd samples :: k=201
	__m256 v230 = data[230*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v226,v230));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v99,v102),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v227,v229),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v100,v101),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v228,acc1);
	data[201*stride] = acc1;

	// even samples :: k=202
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v227,v230));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v100,v102),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v228,v229),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v101,acc1);
	data[202*stride] = acc1;

	// odd samples :: k=203
	__m256 v231 = data[231*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v227,v231));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v100,v103),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v228,v230),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v101,v102),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v229,acc1);
	data[203*stride] = acc1;

	// even samples :: k=204
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v228,v231));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v101,v103),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v229,v230),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v102,acc1);
	data[204*stride] = acc1;

	// odd samples :: k=205
	__m256 v232 = data[232*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v228,v232));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v101,v104),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v229,v231),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v102,v103),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v230,acc1);
	data[205*stride] = acc1;

	// even samples :: k=206
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v229,v232));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v102,v104),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v230,v231),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v103,acc1);
	data[206*stride] = acc1;

	// odd samples :: k=207
	__m256 v233 = data[233*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v229,v233));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v102,v105),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v230,v232),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v103,v104),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v231,acc1);
	data[207*stride] = acc1;

	// even samples :: k=208
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v230,v233));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v103,v105),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v231,v232),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v104,acc1);
	data[208*stride] = acc1;

	// odd samples :: k=209
	__m256 v234 = data[234*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v230,v234));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v103,v106),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v231,v233),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v104,v105),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v232,acc1);
	data[209*stride] = acc1;

	// even samples :: k=210
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v231,v234));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v104,v106),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v232,v233),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v105,acc1);
	data[210*stride] = acc1;

	// odd samples :: k=211
	__m256 v235 = data[235*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v231,v235));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v104,v107),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v232,v234),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v105,v106),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v233,acc1);
	data[211*stride] = acc1;

	// even samples :: k=212
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v232,v235));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v105,v107),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v233,v234),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v106,acc1);
	data[212*stride] = acc1;

	// odd samples :: k=213
	__m256 v236 = data[236*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v232,v236));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v105,v108),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v233,v235),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v106,v107),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v234,acc1);
	data[213*stride] = acc1;

	// even samples :: k=214
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v233,v236));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v106,v108),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v234,v235),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v107,acc1);
	data[214*stride] = acc1;

	// odd samples :: k=215
	__m256 v237 = data[237*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v233,v237));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v106,v109),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v234,v236),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v107,v108),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v235,acc1);
	data[215*stride] = acc1;

	// even samples :: k=216
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v234,v237));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v107,v109),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v235,v236),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v108,acc1);
	data[216*stride] = acc1;

	// odd samples :: k=217
	__m256 v238 = data[238*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v234,v238));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v107,v110),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v235,v237),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v108,v109),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v236,acc1);
	data[217*stride] = acc1;

	// even samples :: k=218
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v235,v238));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v108,v110),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v236,v237),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v109,acc1);
	data[218*stride] = acc1;

	// odd samples :: k=219
	__m256 v239 = data[239*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v235,v239));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v108,v111),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v236,v238),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v109,v110),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v237,acc1);
	data[219*stride] = acc1;

	// even samples :: k=220
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v236,v239));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v109,v111),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v237,v238),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v110,acc1);
	data[220*stride] = acc1;

	// odd samples :: k=221
	__m256 v240 = data[240*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v236,v240));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v109,v112),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v237,v239),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v110,v111),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v238,acc1);
	data[221*stride] = acc1;

	// even samples :: k=222
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v237,v240));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v110,v112),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v238,v239),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v111,acc1);
	data[222*stride] = acc1;

	// odd samples :: k=223
	__m256 v241 = data[241*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v237,v241));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v110,v113),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v238,v240),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v111,v112),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v239,acc1);
	data[223*stride] = acc1;

	// even samples :: k=224
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v238,v241));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v111,v113),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v239,v240),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v112,acc1);
	data[224*stride] = acc1;

	// odd samples :: k=225
	__m256 v242 = data[242*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v238,v242));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v111,v114),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v239,v241),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v112,v113),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v240,acc1);
	data[225*stride] = acc1;

	// even samples :: k=226
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v239,v242));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v112,v114),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v240,v241),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v113,acc1);
	data[226*stride] = acc1;

	// odd samples :: k=227
	__m256 v243 = data[243*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v239,v243));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v112,v115),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v240,v242),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v113,v114),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v241,acc1);
	data[227*stride] = acc1;

	// even samples :: k=228
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v240,v243));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v113,v115),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v241,v242),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v114,acc1);
	data[228*stride] = acc1;

	// odd samples :: k=229
	__m256 v244 = data[244*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v240,v244));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v113,v116),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v241,v243),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v114,v115),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v242,acc1);
	data[229*stride] = acc1;

	// even samples :: k=230
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v241,v244));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v114,v116),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v242,v243),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v115,acc1);
	data[230*stride] = acc1;

	// odd samples :: k=231
	__m256 v245 = data[245*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v241,v245));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v114,v117),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v242,v244),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v115,v116),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v243,acc1);
	data[231*stride] = acc1;

	// even samples :: k=232
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v242,v245));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v115,v117),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v243,v244),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v116,acc1);
	data[232*stride] = acc1;

	// odd samples :: k=233
	__m256 v246 = data[246*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v242,v246));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v115,v118),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v243,v245),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v116,v117),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v244,acc1);
	data[233*stride] = acc1;

	// even samples :: k=234
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v243,v246));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v116,v118),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v244,v245),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v117,acc1);
	data[234*stride] = acc1;

	// odd samples :: k=235
	__m256 v247 = data[247*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v243,v247));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v116,v119),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v244,v246),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v117,v118),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v245,acc1);
	data[235*stride] = acc1;

	// even samples :: k=236
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v244,v247));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v117,v119),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v245,v246),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v118,acc1);
	data[236*stride] = acc1;

	// odd samples :: k=237
	__m256 v248 = data[248*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v244,v248));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v117,v120),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v245,v247),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v118,v119),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v246,acc1);
	data[237*stride] = acc1;

	// even samples :: k=238
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v245,v248));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v118,v120),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v246,v247),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v119,acc1);
	data[238*stride] = acc1;

	// odd samples :: k=239
	__m256 v249 = data[249*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v245,v249));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v118,v121),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v246,v248),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v119,v120),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v247,acc1);
	data[239*stride] = acc1;

	// even samples :: k=240
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v246,v249));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v119,v121),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v247,v248),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v120,acc1);
	data[240*stride] = acc1;

	// odd samples :: k=241
	__m256 v250 = data[250*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v246,v250));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v119,v122),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v247,v249),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v120,v121),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v248,acc1);
	data[241*stride] = acc1;

	// even samples :: k=242
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v247,v250));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v120,v122),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v248,v249),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v121,acc1);
	data[242*stride] = acc1;

	// odd samples :: k=243
	__m256 v251 = data[251*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v247,v251));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v120,v123),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v248,v250),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v121,v122),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v249,acc1);
	data[243*stride] = acc1;

	// even samples :: k=244
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v248,v251));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v121,v123),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v249,v250),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v122,acc1);
	data[244*stride] = acc1;

	// odd samples :: k=245
	__m256 v252 = data[252*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v248,v252));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v121,v124),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v249,v251),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v122,v123),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v250,acc1);
	data[245*stride] = acc1;

	// even samples :: k=246
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v249,v252));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v122,v124),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v250,v251),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v123,acc1);
	data[246*stride] = acc1;

	// odd samples :: k=247
	__m256 v253 = data[253*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v249,v253));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v122,v125),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v250,v252),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v123,v124),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v251,acc1);
	data[247*stride] = acc1;

	// even samples :: k=248
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v250,v253));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v123,v125),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v251,v252),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v124,acc1);
	data[248*stride] = acc1;

	// odd samples :: k=249
	__m256 v254 = data[254*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v250,v254));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v123,v126),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v251,v253),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v124,v125),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v252,acc1);
	data[249*stride] = acc1;

	// even samples :: k=250
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v251,v254));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v124,v126),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v252,v253),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v125,acc1);
	data[250*stride] = acc1;

	// odd samples :: k=251
	__m256 v255 = data[255*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v251,v255));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v124,v127),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v252,v254),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v125,v126),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v253,acc1);
	data[251*stride] = acc1;

	// even samples :: k=252
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v252,v255));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v125,v127),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v253,v254),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v126,acc1);
	data[252*stride] = acc1;

	// odd samples :: k=253
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v252,v254));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v125,v127),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v253,v255),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v126,v127),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v254,acc1);
	data[253*stride] = acc1;

	// even samples :: k=254
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v253,v254));
	acc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v126,v127),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v254,v255),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl0,v127,acc1);
	data[254*stride] = acc1;

	// odd samples :: k=255
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v253,v253));
	acc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v126,v126),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v254,v254),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v127,v127),acc1);
	acc1 = _mm256_fmadd_ps(_mm_sh0,v255,acc1);
	data[255*stride] = acc1;
}

#else

static inline void _Us79_AVX_2(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v1 = data[1*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v1,v1));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v0,v0)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v1,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v0));
	data[0*stride] = acc1;

	// odd samples :: k=1
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v1,v1));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v0,v0)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v1,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v0,v0)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v1));
	data[1*stride] = acc1;
}

static inline void _Us79_AVX_4(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v3 = data[3*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v3,v3));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v1,v1)));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v2,v2)));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v0));
	data[0*stride] = acc1;

	// odd samples :: k=1
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v3,v2));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v1,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v2,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v0,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v2));
	data[1*stride] = acc1;

	// even samples :: k=2
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v2,v2));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v0,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v2,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v1));
	data[2*stride] = acc1;

	// odd samples :: k=3
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v2,v2));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v0,v0)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v2,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v1,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v3));
	data[3*stride] = acc1;
}

static inline void _Us79_AVX_8(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v5 = data[5*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v5,v5));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v1,v1)));
	__m256 v4 = data[4*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v4,v4)));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v0));
	data[0*stride] = acc1;

	// odd samples :: k=1
	__m256 v6 = data[6*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v5,v6));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v1,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v4,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v0,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v4));
	data[1*stride] = acc1;

	// even samples :: k=2
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v4,v6));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v0,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v4,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v1));
	data[2*stride] = acc1;

	// odd samples :: k=3
	__m256 v7 = data[7*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v4,v7));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v0,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v4,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v1,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v5));
	data[3*stride] = acc1;

	// even samples :: k=4
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v4,v7));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v5,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v2));
	data[4*stride] = acc1;

	// odd samples :: k=5
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v4,v6));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v5,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v2,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v6));
	data[5*stride] = acc1;

	// even samples :: k=6
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v5,v6));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v2,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v6,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v3));
	data[6*stride] = acc1;

	// odd samples :: k=7
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v5,v5));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v2,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v6,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v3,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v7));
	data[7*stride] = acc1;
}

static inline void _Us79_AVX_16(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v9 = data[9*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v9,v9));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v1,v1)));
	__m256 v8 = data[8*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v8,v8)));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v0));
	data[0*stride] = acc1;

	// odd samples :: k=1
	__m256 v10 = data[10*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v9,v10));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v1,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v8,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v0,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v8));
	data[1*stride] = acc1;

	// even samples :: k=2
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v8,v10));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v0,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v8,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v1));
	data[2*stride] = acc1;

	// odd samples :: k=3
	__m256 v11 = data[11*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v8,v11));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v0,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v8,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v1,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v9));
	data[3*stride] = acc1;

	// even samples :: k=4
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v8,v11));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v9,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v2));
	__m256 v4 = data[4*stride];
	data[4*stride] = acc1;

	// odd samples :: k=5
	__m256 v12 = data[12*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v8,v12));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v1,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v9,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v2,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v10));
	__m256 v5 = data[5*stride];
	data[5*stride] = acc1;

	// even samples :: k=6
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v9,v12));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v2,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v10,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v3));
	__m256 v6 = data[6*stride];
	data[6*stride] = acc1;

	// odd samples :: k=7
	__m256 v13 = data[13*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v9,v13));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v2,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v10,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v3,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v11));
	__m256 v7 = data[7*stride];
	data[7*stride] = acc1;

	// even samples :: k=8
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v10,v13));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v3,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v11,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v4));
	data[8*stride] = acc1;

	// odd samples :: k=9
	__m256 v14 = data[14*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v10,v14));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v3,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v11,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v4,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v12));
	data[9*stride] = acc1;

	// even samples :: k=10
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v11,v14));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v4,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v12,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v5));
	data[10*stride] = acc1;

	// odd samples :: k=11
	__m256 v15 = data[15*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v11,v15));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v4,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v12,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v5,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v13));
	data[11*stride] = acc1;

	// even samples :: k=12
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v12,v15));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v5,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v13,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v6));
	data[12*stride] = acc1;

	// odd samples :: k=13
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v12,v14));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v5,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v13,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v6,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v14));
	data[13*stride] = acc1;

	// even samples :: k=14
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v13,v14));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v6,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v14,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v7));
	data[14*stride] = acc1;

	// odd samples :: k=15
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v13,v13));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v6,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v14,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v7,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v15));
	data[15*stride] = acc1;
}

static inline void _Us79_AVX_32(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v17 = data[17*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v17,v17));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v1,v1)));
	__m256 v16 = data[16*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v16,v16)));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v0));
	data[0*stride] = acc1;

	// odd samples :: k=1
	__m256 v18 = data[18*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v17,v18));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v1,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v16,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v0,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v16));
	data[1*stride] = acc1;

	// even samples :: k=2
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v16,v18));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v0,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v16,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v1));
	data[2*stride] = acc1;

	// odd samples :: k=3
	__m256 v19 = data[19*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v16,v19));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v0,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v16,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v1,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v17));
	data[3*stride] = acc1;

	// even samples :: k=4
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v16,v19));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v17,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v2));
	__m256 v4 = data[4*stride];
	data[4*stride] = acc1;

	// odd samples :: k=5
	__m256 v20 = data[20*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v16,v20));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v1,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v17,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v2,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v18));
	__m256 v5 = data[5*stride];
	data[5*stride] = acc1;

	// even samples :: k=6
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v17,v20));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v2,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v18,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v3));
	__m256 v6 = data[6*stride];
	data[6*stride] = acc1;

	// odd samples :: k=7
	__m256 v21 = data[21*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v17,v21));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v2,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v18,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v3,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v19));
	__m256 v7 = data[7*stride];
	data[7*stride] = acc1;

	// even samples :: k=8
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v18,v21));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v3,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v19,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v4));
	__m256 v8 = data[8*stride];
	data[8*stride] = acc1;

	// odd samples :: k=9
	__m256 v22 = data[22*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v18,v22));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v3,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v19,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v4,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v20));
	__m256 v9 = data[9*stride];
	data[9*stride] = acc1;

	// even samples :: k=10
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v19,v22));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v4,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v20,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v5));
	__m256 v10 = data[10*stride];
	data[10*stride] = acc1;

	// odd samples :: k=11
	__m256 v23 = data[23*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v19,v23));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v4,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v20,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v5,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v21));
	__m256 v11 = data[11*stride];
	data[11*stride] = acc1;

	// even samples :: k=12
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v20,v23));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v5,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v21,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v6));
	__m256 v12 = data[12*stride];
	data[12*stride] = acc1;

	// odd samples :: k=13
	__m256 v24 = data[24*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v20,v24));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v5,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v21,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v6,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v22));
	__m256 v13 = data[13*stride];
	data[13*stride] = acc1;

	// even samples :: k=14
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v21,v24));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v6,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v22,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v7));
	__m256 v14 = data[14*stride];
	data[14*stride] = acc1;

	// odd samples :: k=15
	__m256 v25 = data[25*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v21,v25));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v6,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v22,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v7,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v23));
	__m256 v15 = data[15*stride];
	data[15*stride] = acc1;

	// even samples :: k=16
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v22,v25));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v7,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v23,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v8));
	data[16*stride] = acc1;

	// odd samples :: k=17
	__m256 v26 = data[26*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v22,v26));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v7,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v23,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v8,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v24));
	data[17*stride] = acc1;

	// even samples :: k=18
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v23,v26));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v8,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v24,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v9));
	data[18*stride] = acc1;

	// odd samples :: k=19
	__m256 v27 = data[27*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v23,v27));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v8,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v24,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v9,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v25));
	data[19*stride] = acc1;

	// even samples :: k=20
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v24,v27));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v9,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v25,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v10));
	data[20*stride] = acc1;

	// odd samples :: k=21
	__m256 v28 = data[28*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v24,v28));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v9,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v25,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v10,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v26));
	data[21*stride] = acc1;

	// even samples :: k=22
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v25,v28));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v10,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v26,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v11));
	data[22*stride] = acc1;

	// odd samples :: k=23
	__m256 v29 = data[29*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v25,v29));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v10,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v26,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v11,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v27));
	data[23*stride] = acc1;

	// even samples :: k=24
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v26,v29));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v11,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v27,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v12));
	data[24*stride] = acc1;

	// odd samples :: k=25
	__m256 v30 = data[30*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v26,v30));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v11,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v27,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v12,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v28));
	data[25*stride] = acc1;

	// even samples :: k=26
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v27,v30));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v12,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v28,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v13));
	data[26*stride] = acc1;

	// odd samples :: k=27
	__m256 v31 = data[31*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v27,v31));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v12,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v28,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v13,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v29));
	data[27*stride] = acc1;

	// even samples :: k=28
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v28,v31));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v13,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v29,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v14));
	data[28*stride] = acc1;

	// odd samples :: k=29
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v28,v30));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v13,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v29,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v14,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v30));
	data[29*stride] = acc1;

	// even samples :: k=30
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v29,v30));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v14,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v30,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v15));
	data[30*stride] = acc1;

	// odd samples :: k=31
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v29,v29));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v14,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v30,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v15,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v31));
	data[31*stride] = acc1;
}

static void _Us79_AVX_64(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v33 = data[33*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v33,v33));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v1,v1)));
	__m256 v32 = data[32*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v32,v32)));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v0));
	data[0*stride] = acc1;

	// odd samples :: k=1
	__m256 v34 = data[34*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v33,v34));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v1,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v32,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v0,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v32));
	data[1*stride] = acc1;

	// even samples :: k=2
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v32,v34));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v0,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v32,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v1));
	data[2*stride] = acc1;

	// odd samples :: k=3
	__m256 v35 = data[35*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v32,v35));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v0,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v32,v34)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v1,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v33));
	data[3*stride] = acc1;

	// even samples :: k=4
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v32,v35));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v33,v34)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v2));
	__m256 v4 = data[4*stride];
	data[4*stride] = acc1;

	// odd samples :: k=5
	__m256 v36 = data[36*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v32,v36));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v1,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v33,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v2,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v34));
	__m256 v5 = data[5*stride];
	data[5*stride] = acc1;

	// even samples :: k=6
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v33,v36));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v2,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v34,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v3));
	__m256 v6 = data[6*stride];
	data[6*stride] = acc1;

	// odd samples :: k=7
	__m256 v37 = data[37*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v33,v37));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v2,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v34,v36)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v3,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v35));
	__m256 v7 = data[7*stride];
	data[7*stride] = acc1;

	// even samples :: k=8
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v34,v37));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v3,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v35,v36)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v4));
	__m256 v8 = data[8*stride];
	data[8*stride] = acc1;

	// odd samples :: k=9
	__m256 v38 = data[38*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v34,v38));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v3,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v35,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v4,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v36));
	__m256 v9 = data[9*stride];
	data[9*stride] = acc1;

	// even samples :: k=10
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v35,v38));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v4,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v36,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v5));
	__m256 v10 = data[10*stride];
	data[10*stride] = acc1;

	// odd samples :: k=11
	__m256 v39 = data[39*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v35,v39));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v4,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v36,v38)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v5,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v37));
	__m256 v11 = data[11*stride];
	data[11*stride] = acc1;

	// even samples :: k=12
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v36,v39));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v5,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v37,v38)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v6));
	__m256 v12 = data[12*stride];
	data[12*stride] = acc1;

	// odd samples :: k=13
	__m256 v40 = data[40*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v36,v40));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v5,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v37,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v6,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v38));
	__m256 v13 = data[13*stride];
	data[13*stride] = acc1;

	// even samples :: k=14
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v37,v40));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v6,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v38,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v7));
	__m256 v14 = data[14*stride];
	data[14*stride] = acc1;

	// odd samples :: k=15
	__m256 v41 = data[41*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v37,v41));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v6,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v38,v40)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v7,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v39));
	__m256 v15 = data[15*stride];
	data[15*stride] = acc1;

	// even samples :: k=16
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v38,v41));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v7,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v39,v40)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v8));
	__m256 v16 = data[16*stride];
	data[16*stride] = acc1;

	// odd samples :: k=17
	__m256 v42 = data[42*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v38,v42));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v7,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v39,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v8,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v40));
	__m256 v17 = data[17*stride];
	data[17*stride] = acc1;

	// even samples :: k=18
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v39,v42));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v8,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v40,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v9));
	__m256 v18 = data[18*stride];
	data[18*stride] = acc1;

	// odd samples :: k=19
	__m256 v43 = data[43*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v39,v43));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v8,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v40,v42)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v9,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v41));
	__m256 v19 = data[19*stride];
	data[19*stride] = acc1;

	// even samples :: k=20
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v40,v43));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v9,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v41,v42)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v10));
	__m256 v20 = data[20*stride];
	data[20*stride] = acc1;

	// odd samples :: k=21
	__m256 v44 = data[44*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v40,v44));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v9,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v41,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v10,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v42));
	__m256 v21 = data[21*stride];
	data[21*stride] = acc1;

	// even samples :: k=22
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v41,v44));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v10,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v42,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v11));
	__m256 v22 = data[22*stride];
	data[22*stride] = acc1;

	// odd samples :: k=23
	__m256 v45 = data[45*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v41,v45));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v10,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v42,v44)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v11,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v43));
	__m256 v23 = data[23*stride];
	data[23*stride] = acc1;

	// even samples :: k=24
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v42,v45));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v11,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v43,v44)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v12));
	__m256 v24 = data[24*stride];
	data[24*stride] = acc1;

	// odd samples :: k=25
	__m256 v46 = data[46*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v42,v46));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v11,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v43,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v12,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v44));
	__m256 v25 = data[25*stride];
	data[25*stride] = acc1;

	// even samples :: k=26
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v43,v46));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v12,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v44,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v13));
	__m256 v26 = data[26*stride];
	data[26*stride] = acc1;

	// odd samples :: k=27
	__m256 v47 = data[47*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v43,v47));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v12,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v44,v46)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v13,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v45));
	__m256 v27 = data[27*stride];
	data[27*stride] = acc1;

	// even samples :: k=28
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v44,v47));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v13,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v45,v46)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v14));
	__m256 v28 = data[28*stride];
	data[28*stride] = acc1;

	// odd samples :: k=29
	__m256 v48 = data[48*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v44,v48));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v13,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v45,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v14,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v46));
	__m256 v29 = data[29*stride];
	data[29*stride] = acc1;

	// even samples :: k=30
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v45,v48));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v14,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v46,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v15));
	__m256 v30 = data[30*stride];
	data[30*stride] = acc1;

	// odd samples :: k=31
	__m256 v49 = data[49*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v45,v49));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v14,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v46,v48)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v15,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v47));
	__m256 v31 = data[31*stride];
	data[31*stride] = acc1;

	// even samples :: k=32
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v46,v49));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v15,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v47,v48)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v16));
	data[32*stride] = acc1;

	// odd samples :: k=33
	__m256 v50 = data[50*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v46,v50));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v15,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v47,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v16,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v48));
	data[33*stride] = acc1;

	// even samples :: k=34
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v47,v50));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v16,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v48,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v17));
	data[34*stride] = acc1;

	// odd samples :: k=35
	__m256 v51 = data[51*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v47,v51));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v16,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v48,v50)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v17,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v49));
	data[35*stride] = acc1;

	// even samples :: k=36
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v48,v51));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v17,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v49,v50)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v18));
	data[36*stride] = acc1;

	// odd samples :: k=37
	__m256 v52 = data[52*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v48,v52));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v17,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v49,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v18,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v50));
	data[37*stride] = acc1;

	// even samples :: k=38
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v49,v52));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v18,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v50,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v19));
	data[38*stride] = acc1;

	// odd samples :: k=39
	__m256 v53 = data[53*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v49,v53));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v18,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v50,v52)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v19,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v51));
	data[39*stride] = acc1;

	// even samples :: k=40
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v50,v53));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v19,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v51,v52)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v20));
	data[40*stride] = acc1;

	// odd samples :: k=41
	__m256 v54 = data[54*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v50,v54));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v19,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v51,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v20,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v52));
	data[41*stride] = acc1;

	// even samples :: k=42
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v51,v54));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v20,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v52,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v21));
	data[42*stride] = acc1;

	// odd samples :: k=43
	__m256 v55 = data[55*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v51,v55));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v20,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v52,v54)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v21,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v53));
	data[43*stride] = acc1;

	// even samples :: k=44
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v52,v55));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v21,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v53,v54)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v22));
	data[44*stride] = acc1;

	// odd samples :: k=45
	__m256 v56 = data[56*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v52,v56));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v21,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v53,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v22,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v54));
	data[45*stride] = acc1;

	// even samples :: k=46
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v53,v56));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v22,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v54,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v23));
	data[46*stride] = acc1;

	// odd samples :: k=47
	__m256 v57 = data[57*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v53,v57));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v22,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v54,v56)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v23,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v55));
	data[47*stride] = acc1;

	// even samples :: k=48
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v54,v57));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v23,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v55,v56)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v24));
	data[48*stride] = acc1;

	// odd samples :: k=49
	__m256 v58 = data[58*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v54,v58));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v23,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v55,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v24,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v56));
	data[49*stride] = acc1;

	// even samples :: k=50
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v55,v58));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v24,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v56,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v25));
	data[50*stride] = acc1;

	// odd samples :: k=51
	__m256 v59 = data[59*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v55,v59));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v24,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v56,v58)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v25,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v57));
	data[51*stride] = acc1;

	// even samples :: k=52
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v56,v59));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v25,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v57,v58)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v26));
	data[52*stride] = acc1;

	// odd samples :: k=53
	__m256 v60 = data[60*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v56,v60));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v25,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v57,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v26,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v58));
	data[53*stride] = acc1;

	// even samples :: k=54
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v57,v60));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v26,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v58,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v27));
	data[54*stride] = acc1;

	// odd samples :: k=55
	__m256 v61 = data[61*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v57,v61));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v26,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v58,v60)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v27,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v59));
	data[55*stride] = acc1;

	// even samples :: k=56
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v58,v61));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v27,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v59,v60)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v28));
	data[56*stride] = acc1;

	// odd samples :: k=57
	__m256 v62 = data[62*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v58,v62));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v27,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v59,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v28,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v60));
	data[57*stride] = acc1;

	// even samples :: k=58
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v59,v62));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v28,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v60,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v29));
	data[58*stride] = acc1;

	// odd samples :: k=59
	__m256 v63 = data[63*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v59,v63));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v28,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v60,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v29,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v61));
	data[59*stride] = acc1;

	// even samples :: k=60
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v60,v63));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v29,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v61,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v30));
	data[60*stride] = acc1;

	// odd samples :: k=61
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v60,v62));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v29,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v61,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v30,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v62));
	data[61*stride] = acc1;

	// even samples :: k=62
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v61,v62));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v30,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v62,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v31));
	data[62*stride] = acc1;

	// odd samples :: k=63
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v61,v61));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v30,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v62,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v31,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v63));
	data[63*stride] = acc1;
}

static void _Us79_AVX_128(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v65 = data[65*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v65,v65));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v1,v1)));
	__m256 v64 = data[64*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v64,v64)));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v0));
	data[0*stride] = acc1;

	// odd samples :: k=1
	__m256 v66 = data[66*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v65,v66));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v1,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v64,v65)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v0,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v64));
	data[1*stride] = acc1;

	// even samples :: k=2
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v64,v66));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v0,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v64,v65)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v1));
	data[2*stride] = acc1;

	// odd samples :: k=3
	__m256 v67 = data[67*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v64,v67));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v0,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v64,v66)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v1,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v65));
	data[3*stride] = acc1;

	// even samples :: k=4
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v64,v67));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v65,v66)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v2));
	__m256 v4 = data[4*stride];
	data[4*stride] = acc1;

	// odd samples :: k=5
	__m256 v68 = data[68*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v64,v68));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v1,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v65,v67)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v2,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v66));
	__m256 v5 = data[5*stride];
	data[5*stride] = acc1;

	// even samples :: k=6
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v65,v68));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v2,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v66,v67)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v3));
	__m256 v6 = data[6*stride];
	data[6*stride] = acc1;

	// odd samples :: k=7
	__m256 v69 = data[69*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v65,v69));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v2,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v66,v68)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v3,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v67));
	__m256 v7 = data[7*stride];
	data[7*stride] = acc1;

	// even samples :: k=8
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v66,v69));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v3,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v67,v68)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v4));
	__m256 v8 = data[8*stride];
	data[8*stride] = acc1;

	// odd samples :: k=9
	__m256 v70 = data[70*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v66,v70));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v3,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v67,v69)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v4,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v68));
	__m256 v9 = data[9*stride];
	data[9*stride] = acc1;

	// even samples :: k=10
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v67,v70));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v4,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v68,v69)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v5));
	__m256 v10 = data[10*stride];
	data[10*stride] = acc1;

	// odd samples :: k=11
	__m256 v71 = data[71*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v67,v71));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v4,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v68,v70)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v5,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v69));
	__m256 v11 = data[11*stride];
	data[11*stride] = acc1;

	// even samples :: k=12
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v68,v71));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v5,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v69,v70)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v6));
	__m256 v12 = data[12*stride];
	data[12*stride] = acc1;

	// odd samples :: k=13
	__m256 v72 = data[72*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v68,v72));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v5,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v69,v71)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v6,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v70));
	__m256 v13 = data[13*stride];
	data[13*stride] = acc1;

	// even samples :: k=14
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v69,v72));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v6,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v70,v71)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v7));
	__m256 v14 = data[14*stride];
	data[14*stride] = acc1;

	// odd samples :: k=15
	__m256 v73 = data[73*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v69,v73));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v6,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v70,v72)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v7,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v71));
	__m256 v15 = data[15*stride];
	data[15*stride] = acc1;

	// even samples :: k=16
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v70,v73));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v7,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v71,v72)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v8));
	__m256 v16 = data[16*stride];
	data[16*stride] = acc1;

	// odd samples :: k=17
	__m256 v74 = data[74*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v70,v74));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v7,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v71,v73)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v8,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v72));
	__m256 v17 = data[17*stride];
	data[17*stride] = acc1;

	// even samples :: k=18
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v71,v74));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v8,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v72,v73)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v9));
	__m256 v18 = data[18*stride];
	data[18*stride] = acc1;

	// odd samples :: k=19
	__m256 v75 = data[75*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v71,v75));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v8,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v72,v74)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v9,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v73));
	__m256 v19 = data[19*stride];
	data[19*stride] = acc1;

	// even samples :: k=20
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v72,v75));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v9,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v73,v74)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v10));
	__m256 v20 = data[20*stride];
	data[20*stride] = acc1;

	// odd samples :: k=21
	__m256 v76 = data[76*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v72,v76));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v9,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v73,v75)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v10,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v74));
	__m256 v21 = data[21*stride];
	data[21*stride] = acc1;

	// even samples :: k=22
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v73,v76));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v10,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v74,v75)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v11));
	__m256 v22 = data[22*stride];
	data[22*stride] = acc1;

	// odd samples :: k=23
	__m256 v77 = data[77*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v73,v77));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v10,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v74,v76)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v11,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v75));
	__m256 v23 = data[23*stride];
	data[23*stride] = acc1;

	// even samples :: k=24
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v74,v77));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v11,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v75,v76)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v12));
	__m256 v24 = data[24*stride];
	data[24*stride] = acc1;

	// odd samples :: k=25
	__m256 v78 = data[78*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v74,v78));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v11,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v75,v77)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v12,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v76));
	__m256 v25 = data[25*stride];
	data[25*stride] = acc1;

	// even samples :: k=26
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v75,v78));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v12,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v76,v77)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v13));
	__m256 v26 = data[26*stride];
	data[26*stride] = acc1;

	// odd samples :: k=27
	__m256 v79 = data[79*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v75,v79));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v12,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v76,v78)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v13,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v77));
	__m256 v27 = data[27*stride];
	data[27*stride] = acc1;

	// even samples :: k=28
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v76,v79));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v13,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v77,v78)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v14));
	__m256 v28 = data[28*stride];
	data[28*stride] = acc1;

	// odd samples :: k=29
	__m256 v80 = data[80*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v76,v80));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v13,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v77,v79)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v14,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v78));
	__m256 v29 = data[29*stride];
	data[29*stride] = acc1;

	// even samples :: k=30
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v77,v80));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v14,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v78,v79)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v15));
	__m256 v30 = data[30*stride];
	data[30*stride] = acc1;

	// odd samples :: k=31
	__m256 v81 = data[81*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v77,v81));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v14,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v78,v80)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v15,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v79));
	__m256 v31 = data[31*stride];
	data[31*stride] = acc1;

	// even samples :: k=32
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v78,v81));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v15,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v79,v80)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v16));
	__m256 v32 = data[32*stride];
	data[32*stride] = acc1;

	// odd samples :: k=33
	__m256 v82 = data[82*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v78,v82));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v15,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v79,v81)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v16,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v80));
	__m256 v33 = data[33*stride];
	data[33*stride] = acc1;

	// even samples :: k=34
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v79,v82));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v16,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v80,v81)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v17));
	__m256 v34 = data[34*stride];
	data[34*stride] = acc1;

	// odd samples :: k=35
	__m256 v83 = data[83*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v79,v83));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v16,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v80,v82)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v17,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v81));
	__m256 v35 = data[35*stride];
	data[35*stride] = acc1;

	// even samples :: k=36
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v80,v83));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v17,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v81,v82)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v18));
	__m256 v36 = data[36*stride];
	data[36*stride] = acc1;

	// odd samples :: k=37
	__m256 v84 = data[84*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v80,v84));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v17,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v81,v83)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v18,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v82));
	__m256 v37 = data[37*stride];
	data[37*stride] = acc1;

	// even samples :: k=38
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v81,v84));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v18,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v82,v83)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v19));
	__m256 v38 = data[38*stride];
	data[38*stride] = acc1;

	// odd samples :: k=39
	__m256 v85 = data[85*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v81,v85));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v18,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v82,v84)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v19,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v83));
	__m256 v39 = data[39*stride];
	data[39*stride] = acc1;

	// even samples :: k=40
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v82,v85));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v19,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v83,v84)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v20));
	__m256 v40 = data[40*stride];
	data[40*stride] = acc1;

	// odd samples :: k=41
	__m256 v86 = data[86*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v82,v86));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v19,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v83,v85)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v20,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v84));
	__m256 v41 = data[41*stride];
	data[41*stride] = acc1;

	// even samples :: k=42
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v83,v86));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v20,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v84,v85)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v21));
	__m256 v42 = data[42*stride];
	data[42*stride] = acc1;

	// odd samples :: k=43
	__m256 v87 = data[87*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v83,v87));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v20,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v84,v86)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v21,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v85));
	__m256 v43 = data[43*stride];
	data[43*stride] = acc1;

	// even samples :: k=44
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v84,v87));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v21,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v85,v86)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v22));
	__m256 v44 = data[44*stride];
	data[44*stride] = acc1;

	// odd samples :: k=45
	__m256 v88 = data[88*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v84,v88));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v21,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v85,v87)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v22,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v86));
	__m256 v45 = data[45*stride];
	data[45*stride] = acc1;

	// even samples :: k=46
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v85,v88));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v22,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v86,v87)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v23));
	__m256 v46 = data[46*stride];
	data[46*stride] = acc1;

	// odd samples :: k=47
	__m256 v89 = data[89*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v85,v89));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v22,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v86,v88)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v23,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v87));
	__m256 v47 = data[47*stride];
	data[47*stride] = acc1;

	// even samples :: k=48
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v86,v89));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v23,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v87,v88)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v24));
	__m256 v48 = data[48*stride];
	data[48*stride] = acc1;

	// odd samples :: k=49
	__m256 v90 = data[90*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v86,v90));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v23,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v87,v89)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v24,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v88));
	__m256 v49 = data[49*stride];
	data[49*stride] = acc1;

	// even samples :: k=50
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v87,v90));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v24,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v88,v89)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v25));
	__m256 v50 = data[50*stride];
	data[50*stride] = acc1;

	// odd samples :: k=51
	__m256 v91 = data[91*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v87,v91));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v24,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v88,v90)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v25,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v89));
	__m256 v51 = data[51*stride];
	data[51*stride] = acc1;

	// even samples :: k=52
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v88,v91));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v25,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v89,v90)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v26));
	__m256 v52 = data[52*stride];
	data[52*stride] = acc1;

	// odd samples :: k=53
	__m256 v92 = data[92*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v88,v92));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v25,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v89,v91)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v26,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v90));
	__m256 v53 = data[53*stride];
	data[53*stride] = acc1;

	// even samples :: k=54
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v89,v92));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v26,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v90,v91)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v27));
	__m256 v54 = data[54*stride];
	data[54*stride] = acc1;

	// odd samples :: k=55
	__m256 v93 = data[93*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v89,v93));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v26,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v90,v92)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v27,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v91));
	__m256 v55 = data[55*stride];
	data[55*stride] = acc1;

	// even samples :: k=56
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v90,v93));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v27,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v91,v92)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v28));
	__m256 v56 = data[56*stride];
	data[56*stride] = acc1;

	// odd samples :: k=57
	__m256 v94 = data[94*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v90,v94));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v27,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v91,v93)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v28,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v92));
	__m256 v57 = data[57*stride];
	data[57*stride] = acc1;

	// even samples :: k=58
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v91,v94));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v28,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v92,v93)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v29));
	__m256 v58 = data[58*stride];
	data[58*stride] = acc1;

	// odd samples :: k=59
	__m256 v95 = data[95*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v91,v95));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v28,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v92,v94)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v29,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v93));
	__m256 v59 = data[59*stride];
	data[59*stride] = acc1;

	// even samples :: k=60
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v92,v95));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v29,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v93,v94)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v30));
	__m256 v60 = data[60*stride];
	data[60*stride] = acc1;

	// odd samples :: k=61
	__m256 v96 = data[96*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v92,v96));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v29,v32)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v93,v95)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v30,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v94));
	__m256 v61 = data[61*stride];
	data[61*stride] = acc1;

	// even samples :: k=62
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v93,v96));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v30,v32)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v94,v95)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v31));
	__m256 v62 = data[62*stride];
	data[62*stride] = acc1;

	// odd samples :: k=63
	__m256 v97 = data[97*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v93,v97));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v30,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v94,v96)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v31,v32)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v95));
	__m256 v63 = data[63*stride];
	data[63*stride] = acc1;

	// even samples :: k=64
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v94,v97));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v31,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v95,v96)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v32));
	data[64*stride] = acc1;

	// odd samples :: k=65
	__m256 v98 = data[98*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v94,v98));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v31,v34)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v95,v97)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v32,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v96));
	data[65*stride] = acc1;

	// even samples :: k=66
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v95,v98));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v32,v34)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v96,v97)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v33));
	data[66*stride] = acc1;

	// odd samples :: k=67
	__m256 v99 = data[99*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v95,v99));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v32,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v96,v98)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v33,v34)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v97));
	data[67*stride] = acc1;

	// even samples :: k=68
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v96,v99));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v33,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v97,v98)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v34));
	data[68*stride] = acc1;

	// odd samples :: k=69
	__m256 v100 = data[100*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v96,v100));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v33,v36)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v97,v99)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v34,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v98));
	data[69*stride] = acc1;

	// even samples :: k=70
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v97,v100));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v34,v36)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v98,v99)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v35));
	data[70*stride] = acc1;

	// odd samples :: k=71
	__m256 v101 = data[101*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v97,v101));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v34,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v98,v100)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v35,v36)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v99));
	data[71*stride] = acc1;

	// even samples :: k=72
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v98,v101));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v35,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v99,v100)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v36));
	data[72*stride] = acc1;

	// odd samples :: k=73
	__m256 v102 = data[102*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v98,v102));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v35,v38)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v99,v101)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v36,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v100));
	data[73*stride] = acc1;

	// even samples :: k=74
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v99,v102));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v36,v38)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v100,v101)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v37));
	data[74*stride] = acc1;

	// odd samples :: k=75
	__m256 v103 = data[103*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v99,v103));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v36,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v100,v102)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v37,v38)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v101));
	data[75*stride] = acc1;

	// even samples :: k=76
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v100,v103));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v37,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v101,v102)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v38));
	data[76*stride] = acc1;

	// odd samples :: k=77
	__m256 v104 = data[104*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v100,v104));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v37,v40)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v101,v103)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v38,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v102));
	data[77*stride] = acc1;

	// even samples :: k=78
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v101,v104));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v38,v40)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v102,v103)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v39));
	data[78*stride] = acc1;

	// odd samples :: k=79
	__m256 v105 = data[105*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v101,v105));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v38,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v102,v104)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v39,v40)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v103));
	data[79*stride] = acc1;

	// even samples :: k=80
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v102,v105));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v39,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v103,v104)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v40));
	data[80*stride] = acc1;

	// odd samples :: k=81
	__m256 v106 = data[106*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v102,v106));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v39,v42)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v103,v105)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v40,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v104));
	data[81*stride] = acc1;

	// even samples :: k=82
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v103,v106));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v40,v42)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v104,v105)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v41));
	data[82*stride] = acc1;

	// odd samples :: k=83
	__m256 v107 = data[107*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v103,v107));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v40,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v104,v106)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v41,v42)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v105));
	data[83*stride] = acc1;

	// even samples :: k=84
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v104,v107));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v41,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v105,v106)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v42));
	data[84*stride] = acc1;

	// odd samples :: k=85
	__m256 v108 = data[108*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v104,v108));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v41,v44)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v105,v107)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v42,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v106));
	data[85*stride] = acc1;

	// even samples :: k=86
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v105,v108));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v42,v44)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v106,v107)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v43));
	data[86*stride] = acc1;

	// odd samples :: k=87
	__m256 v109 = data[109*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v105,v109));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v42,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v106,v108)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v43,v44)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v107));
	data[87*stride] = acc1;

	// even samples :: k=88
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v106,v109));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v43,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v107,v108)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v44));
	data[88*stride] = acc1;

	// odd samples :: k=89
	__m256 v110 = data[110*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v106,v110));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v43,v46)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v107,v109)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v44,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v108));
	data[89*stride] = acc1;

	// even samples :: k=90
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v107,v110));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v44,v46)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v108,v109)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v45));
	data[90*stride] = acc1;

	// odd samples :: k=91
	__m256 v111 = data[111*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v107,v111));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v44,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v108,v110)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v45,v46)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v109));
	data[91*stride] = acc1;

	// even samples :: k=92
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v108,v111));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v45,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v109,v110)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v46));
	data[92*stride] = acc1;

	// odd samples :: k=93
	__m256 v112 = data[112*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v108,v112));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v45,v48)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v109,v111)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v46,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v110));
	data[93*stride] = acc1;

	// even samples :: k=94
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v109,v112));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v46,v48)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v110,v111)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v47));
	data[94*stride] = acc1;

	// odd samples :: k=95
	__m256 v113 = data[113*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v109,v113));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v46,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v110,v112)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v47,v48)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v111));
	data[95*stride] = acc1;

	// even samples :: k=96
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v110,v113));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v47,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v111,v112)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v48));
	data[96*stride] = acc1;

	// odd samples :: k=97
	__m256 v114 = data[114*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v110,v114));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v47,v50)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v111,v113)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v48,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v112));
	data[97*stride] = acc1;

	// even samples :: k=98
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v111,v114));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v48,v50)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v112,v113)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v49));
	data[98*stride] = acc1;

	// odd samples :: k=99
	__m256 v115 = data[115*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v111,v115));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v48,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v112,v114)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v49,v50)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v113));
	data[99*stride] = acc1;

	// even samples :: k=100
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v112,v115));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v49,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v113,v114)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v50));
	data[100*stride] = acc1;

	// odd samples :: k=101
	__m256 v116 = data[116*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v112,v116));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v49,v52)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v113,v115)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v50,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v114));
	data[101*stride] = acc1;

	// even samples :: k=102
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v113,v116));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v50,v52)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v114,v115)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v51));
	data[102*stride] = acc1;

	// odd samples :: k=103
	__m256 v117 = data[117*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v113,v117));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v50,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v114,v116)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v51,v52)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v115));
	data[103*stride] = acc1;

	// even samples :: k=104
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v114,v117));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v51,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v115,v116)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v52));
	data[104*stride] = acc1;

	// odd samples :: k=105
	__m256 v118 = data[118*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v114,v118));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v51,v54)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v115,v117)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v52,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v116));
	data[105*stride] = acc1;

	// even samples :: k=106
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v115,v118));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v52,v54)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v116,v117)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v53));
	data[106*stride] = acc1;

	// odd samples :: k=107
	__m256 v119 = data[119*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v115,v119));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v52,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v116,v118)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v53,v54)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v117));
	data[107*stride] = acc1;

	// even samples :: k=108
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v116,v119));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v53,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v117,v118)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v54));
	data[108*stride] = acc1;

	// odd samples :: k=109
	__m256 v120 = data[120*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v116,v120));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v53,v56)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v117,v119)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v54,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v118));
	data[109*stride] = acc1;

	// even samples :: k=110
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v117,v120));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v54,v56)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v118,v119)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v55));
	data[110*stride] = acc1;

	// odd samples :: k=111
	__m256 v121 = data[121*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v117,v121));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v54,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v118,v120)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v55,v56)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v119));
	data[111*stride] = acc1;

	// even samples :: k=112
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v118,v121));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v55,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v119,v120)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v56));
	data[112*stride] = acc1;

	// odd samples :: k=113
	__m256 v122 = data[122*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v118,v122));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v55,v58)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v119,v121)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v56,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v120));
	data[113*stride] = acc1;

	// even samples :: k=114
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v119,v122));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v56,v58)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v120,v121)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v57));
	data[114*stride] = acc1;

	// odd samples :: k=115
	__m256 v123 = data[123*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v119,v123));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v56,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v120,v122)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v57,v58)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v121));
	data[115*stride] = acc1;

	// even samples :: k=116
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v120,v123));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v57,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v121,v122)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v58));
	data[116*stride] = acc1;

	// odd samples :: k=117
	__m256 v124 = data[124*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v120,v124));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v57,v60)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v121,v123)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v58,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v122));
	data[117*stride] = acc1;

	// even samples :: k=118
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v121,v124));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v58,v60)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v122,v123)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v59));
	data[118*stride] = acc1;

	// odd samples :: k=119
	__m256 v125 = data[125*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v121,v125));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v58,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v122,v124)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v59,v60)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v123));
	data[119*stride] = acc1;

	// even samples :: k=120
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v122,v125));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v59,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v123,v124)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v60));
	data[120*stride] = acc1;

	// odd samples :: k=121
	__m256 v126 = data[126*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v122,v126));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v59,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v123,v125)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v60,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v124));
	data[121*stride] = acc1;

	// even samples :: k=122
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v123,v126));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v60,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v124,v125)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v61));
	data[122*stride] = acc1;

	// odd samples :: k=123
	__m256 v127 = data[127*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v123,v127));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v60,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v124,v126)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v61,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v125));
	data[123*stride] = acc1;

	// even samples :: k=124
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v124,v127));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v61,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v125,v126)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v62));
	data[124*stride] = acc1;

	// odd samples :: k=125
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v124,v126));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v61,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v125,v127)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v62,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v126));
	data[125*stride] = acc1;

	// even samples :: k=126
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v125,v126));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v62,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v126,v127)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v63));
	data[126*stride] = acc1;

	// odd samples :: k=127
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v125,v125));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v62,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v126,v126)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v63,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v127));
	data[127*stride] = acc1;
}

static void _Us79_AVX_256(__m256* data, int stride)
{
	__m256 acc1;

	// even samples :: k=0
	__m256 v129 = data[129*stride];
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v129,v129));
	__m256 v1 = data[1*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v1,v1)));
	__m256 v128 = data[128*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v128,v128)));
	__m256 v0 = data[0*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v0));
	data[0*stride] = acc1;

	// odd samples :: k=1
	__m256 v130 = data[130*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v129,v130));
	__m256 v2 = data[2*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v1,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v128,v129)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v0,v1)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v128));
	data[1*stride] = acc1;

	// even samples :: k=2
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v128,v130));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v0,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v128,v129)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v1));
	data[2*stride] = acc1;

	// odd samples :: k=3
	__m256 v131 = data[131*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v128,v131));
	__m256 v3 = data[3*stride];
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v0,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v128,v130)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v1,v2)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v129));
	data[3*stride] = acc1;

	// even samples :: k=4
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v128,v131));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v1,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v129,v130)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v2));
	__m256 v4 = data[4*stride];
	data[4*stride] = acc1;

	// odd samples :: k=5
	__m256 v132 = data[132*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v128,v132));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v1,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v129,v131)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v2,v3)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v130));
	__m256 v5 = data[5*stride];
	data[5*stride] = acc1;

	// even samples :: k=6
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v129,v132));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v2,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v130,v131)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v3));
	__m256 v6 = data[6*stride];
	data[6*stride] = acc1;

	// odd samples :: k=7
	__m256 v133 = data[133*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v129,v133));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v2,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v130,v132)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v3,v4)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v131));
	__m256 v7 = data[7*stride];
	data[7*stride] = acc1;

	// even samples :: k=8
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v130,v133));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v3,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v131,v132)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v4));
	__m256 v8 = data[8*stride];
	data[8*stride] = acc1;

	// odd samples :: k=9
	__m256 v134 = data[134*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v130,v134));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v3,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v131,v133)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v4,v5)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v132));
	__m256 v9 = data[9*stride];
	data[9*stride] = acc1;

	// even samples :: k=10
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v131,v134));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v4,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v132,v133)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v5));
	__m256 v10 = data[10*stride];
	data[10*stride] = acc1;

	// odd samples :: k=11
	__m256 v135 = data[135*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v131,v135));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v4,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v132,v134)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v5,v6)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v133));
	__m256 v11 = data[11*stride];
	data[11*stride] = acc1;

	// even samples :: k=12
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v132,v135));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v5,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v133,v134)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v6));
	__m256 v12 = data[12*stride];
	data[12*stride] = acc1;

	// odd samples :: k=13
	__m256 v136 = data[136*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v132,v136));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v5,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v133,v135)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v6,v7)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v134));
	__m256 v13 = data[13*stride];
	data[13*stride] = acc1;

	// even samples :: k=14
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v133,v136));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v6,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v134,v135)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v7));
	__m256 v14 = data[14*stride];
	data[14*stride] = acc1;

	// odd samples :: k=15
	__m256 v137 = data[137*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v133,v137));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v6,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v134,v136)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v7,v8)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v135));
	__m256 v15 = data[15*stride];
	data[15*stride] = acc1;

	// even samples :: k=16
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v134,v137));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v7,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v135,v136)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v8));
	__m256 v16 = data[16*stride];
	data[16*stride] = acc1;

	// odd samples :: k=17
	__m256 v138 = data[138*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v134,v138));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v7,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v135,v137)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v8,v9)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v136));
	__m256 v17 = data[17*stride];
	data[17*stride] = acc1;

	// even samples :: k=18
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v135,v138));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v8,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v136,v137)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v9));
	__m256 v18 = data[18*stride];
	data[18*stride] = acc1;

	// odd samples :: k=19
	__m256 v139 = data[139*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v135,v139));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v8,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v136,v138)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v9,v10)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v137));
	__m256 v19 = data[19*stride];
	data[19*stride] = acc1;

	// even samples :: k=20
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v136,v139));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v9,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v137,v138)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v10));
	__m256 v20 = data[20*stride];
	data[20*stride] = acc1;

	// odd samples :: k=21
	__m256 v140 = data[140*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v136,v140));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v9,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v137,v139)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v10,v11)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v138));
	__m256 v21 = data[21*stride];
	data[21*stride] = acc1;

	// even samples :: k=22
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v137,v140));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v10,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v138,v139)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v11));
	__m256 v22 = data[22*stride];
	data[22*stride] = acc1;

	// odd samples :: k=23
	__m256 v141 = data[141*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v137,v141));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v10,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v138,v140)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v11,v12)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v139));
	__m256 v23 = data[23*stride];
	data[23*stride] = acc1;

	// even samples :: k=24
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v138,v141));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v11,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v139,v140)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v12));
	__m256 v24 = data[24*stride];
	data[24*stride] = acc1;

	// odd samples :: k=25
	__m256 v142 = data[142*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v138,v142));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v11,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v139,v141)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v12,v13)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v140));
	__m256 v25 = data[25*stride];
	data[25*stride] = acc1;

	// even samples :: k=26
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v139,v142));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v12,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v140,v141)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v13));
	__m256 v26 = data[26*stride];
	data[26*stride] = acc1;

	// odd samples :: k=27
	__m256 v143 = data[143*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v139,v143));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v12,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v140,v142)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v13,v14)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v141));
	__m256 v27 = data[27*stride];
	data[27*stride] = acc1;

	// even samples :: k=28
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v140,v143));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v13,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v141,v142)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v14));
	__m256 v28 = data[28*stride];
	data[28*stride] = acc1;

	// odd samples :: k=29
	__m256 v144 = data[144*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v140,v144));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v13,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v141,v143)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v14,v15)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v142));
	__m256 v29 = data[29*stride];
	data[29*stride] = acc1;

	// even samples :: k=30
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v141,v144));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v14,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v142,v143)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v15));
	__m256 v30 = data[30*stride];
	data[30*stride] = acc1;

	// odd samples :: k=31
	__m256 v145 = data[145*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v141,v145));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v14,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v142,v144)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v15,v16)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v143));
	__m256 v31 = data[31*stride];
	data[31*stride] = acc1;

	// even samples :: k=32
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v142,v145));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v15,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v143,v144)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v16));
	__m256 v32 = data[32*stride];
	data[32*stride] = acc1;

	// odd samples :: k=33
	__m256 v146 = data[146*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v142,v146));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v15,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v143,v145)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v16,v17)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v144));
	__m256 v33 = data[33*stride];
	data[33*stride] = acc1;

	// even samples :: k=34
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v143,v146));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v16,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v144,v145)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v17));
	__m256 v34 = data[34*stride];
	data[34*stride] = acc1;

	// odd samples :: k=35
	__m256 v147 = data[147*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v143,v147));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v16,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v144,v146)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v17,v18)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v145));
	__m256 v35 = data[35*stride];
	data[35*stride] = acc1;

	// even samples :: k=36
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v144,v147));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v17,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v145,v146)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v18));
	__m256 v36 = data[36*stride];
	data[36*stride] = acc1;

	// odd samples :: k=37
	__m256 v148 = data[148*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v144,v148));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v17,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v145,v147)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v18,v19)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v146));
	__m256 v37 = data[37*stride];
	data[37*stride] = acc1;

	// even samples :: k=38
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v145,v148));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v18,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v146,v147)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v19));
	__m256 v38 = data[38*stride];
	data[38*stride] = acc1;

	// odd samples :: k=39
	__m256 v149 = data[149*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v145,v149));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v18,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v146,v148)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v19,v20)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v147));
	__m256 v39 = data[39*stride];
	data[39*stride] = acc1;

	// even samples :: k=40
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v146,v149));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v19,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v147,v148)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v20));
	__m256 v40 = data[40*stride];
	data[40*stride] = acc1;

	// odd samples :: k=41
	__m256 v150 = data[150*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v146,v150));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v19,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v147,v149)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v20,v21)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v148));
	__m256 v41 = data[41*stride];
	data[41*stride] = acc1;

	// even samples :: k=42
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v147,v150));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v20,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v148,v149)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v21));
	__m256 v42 = data[42*stride];
	data[42*stride] = acc1;

	// odd samples :: k=43
	__m256 v151 = data[151*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v147,v151));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v20,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v148,v150)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v21,v22)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v149));
	__m256 v43 = data[43*stride];
	data[43*stride] = acc1;

	// even samples :: k=44
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v148,v151));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v21,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v149,v150)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v22));
	__m256 v44 = data[44*stride];
	data[44*stride] = acc1;

	// odd samples :: k=45
	__m256 v152 = data[152*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v148,v152));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v21,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v149,v151)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v22,v23)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v150));
	__m256 v45 = data[45*stride];
	data[45*stride] = acc1;

	// even samples :: k=46
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v149,v152));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v22,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v150,v151)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v23));
	__m256 v46 = data[46*stride];
	data[46*stride] = acc1;

	// odd samples :: k=47
	__m256 v153 = data[153*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v149,v153));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v22,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v150,v152)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v23,v24)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v151));
	__m256 v47 = data[47*stride];
	data[47*stride] = acc1;

	// even samples :: k=48
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v150,v153));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v23,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v151,v152)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v24));
	__m256 v48 = data[48*stride];
	data[48*stride] = acc1;

	// odd samples :: k=49
	__m256 v154 = data[154*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v150,v154));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v23,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v151,v153)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v24,v25)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v152));
	__m256 v49 = data[49*stride];
	data[49*stride] = acc1;

	// even samples :: k=50
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v151,v154));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v24,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v152,v153)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v25));
	__m256 v50 = data[50*stride];
	data[50*stride] = acc1;

	// odd samples :: k=51
	__m256 v155 = data[155*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v151,v155));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v24,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v152,v154)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v25,v26)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v153));
	__m256 v51 = data[51*stride];
	data[51*stride] = acc1;

	// even samples :: k=52
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v152,v155));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v25,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v153,v154)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v26));
	__m256 v52 = data[52*stride];
	data[52*stride] = acc1;

	// odd samples :: k=53
	__m256 v156 = data[156*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v152,v156));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v25,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v153,v155)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v26,v27)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v154));
	__m256 v53 = data[53*stride];
	data[53*stride] = acc1;

	// even samples :: k=54
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v153,v156));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v26,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v154,v155)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v27));
	__m256 v54 = data[54*stride];
	data[54*stride] = acc1;

	// odd samples :: k=55
	__m256 v157 = data[157*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v153,v157));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v26,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v154,v156)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v27,v28)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v155));
	__m256 v55 = data[55*stride];
	data[55*stride] = acc1;

	// even samples :: k=56
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v154,v157));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v27,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v155,v156)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v28));
	__m256 v56 = data[56*stride];
	data[56*stride] = acc1;

	// odd samples :: k=57
	__m256 v158 = data[158*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v154,v158));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v27,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v155,v157)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v28,v29)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v156));
	__m256 v57 = data[57*stride];
	data[57*stride] = acc1;

	// even samples :: k=58
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v155,v158));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v28,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v156,v157)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v29));
	__m256 v58 = data[58*stride];
	data[58*stride] = acc1;

	// odd samples :: k=59
	__m256 v159 = data[159*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v155,v159));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v28,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v156,v158)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v29,v30)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v157));
	__m256 v59 = data[59*stride];
	data[59*stride] = acc1;

	// even samples :: k=60
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v156,v159));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v29,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v157,v158)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v30));
	__m256 v60 = data[60*stride];
	data[60*stride] = acc1;

	// odd samples :: k=61
	__m256 v160 = data[160*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v156,v160));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v29,v32)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v157,v159)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v30,v31)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v158));
	__m256 v61 = data[61*stride];
	data[61*stride] = acc1;

	// even samples :: k=62
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v157,v160));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v30,v32)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v158,v159)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v31));
	__m256 v62 = data[62*stride];
	data[62*stride] = acc1;

	// odd samples :: k=63
	__m256 v161 = data[161*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v157,v161));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v30,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v158,v160)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v31,v32)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v159));
	__m256 v63 = data[63*stride];
	data[63*stride] = acc1;

	// even samples :: k=64
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v158,v161));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v31,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v159,v160)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v32));
	__m256 v64 = data[64*stride];
	data[64*stride] = acc1;

	// odd samples :: k=65
	__m256 v162 = data[162*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v158,v162));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v31,v34)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v159,v161)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v32,v33)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v160));
	__m256 v65 = data[65*stride];
	data[65*stride] = acc1;

	// even samples :: k=66
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v159,v162));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v32,v34)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v160,v161)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v33));
	__m256 v66 = data[66*stride];
	data[66*stride] = acc1;

	// odd samples :: k=67
	__m256 v163 = data[163*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v159,v163));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v32,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v160,v162)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v33,v34)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v161));
	__m256 v67 = data[67*stride];
	data[67*stride] = acc1;

	// even samples :: k=68
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v160,v163));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v33,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v161,v162)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v34));
	__m256 v68 = data[68*stride];
	data[68*stride] = acc1;

	// odd samples :: k=69
	__m256 v164 = data[164*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v160,v164));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v33,v36)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v161,v163)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v34,v35)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v162));
	__m256 v69 = data[69*stride];
	data[69*stride] = acc1;

	// even samples :: k=70
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v161,v164));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v34,v36)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v162,v163)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v35));
	__m256 v70 = data[70*stride];
	data[70*stride] = acc1;

	// odd samples :: k=71
	__m256 v165 = data[165*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v161,v165));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v34,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v162,v164)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v35,v36)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v163));
	__m256 v71 = data[71*stride];
	data[71*stride] = acc1;

	// even samples :: k=72
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v162,v165));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v35,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v163,v164)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v36));
	__m256 v72 = data[72*stride];
	data[72*stride] = acc1;

	// odd samples :: k=73
	__m256 v166 = data[166*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v162,v166));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v35,v38)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v163,v165)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v36,v37)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v164));
	__m256 v73 = data[73*stride];
	data[73*stride] = acc1;

	// even samples :: k=74
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v163,v166));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v36,v38)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v164,v165)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v37));
	__m256 v74 = data[74*stride];
	data[74*stride] = acc1;

	// odd samples :: k=75
	__m256 v167 = data[167*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v163,v167));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v36,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v164,v166)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v37,v38)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v165));
	__m256 v75 = data[75*stride];
	data[75*stride] = acc1;

	// even samples :: k=76
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v164,v167));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v37,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v165,v166)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v38));
	__m256 v76 = data[76*stride];
	data[76*stride] = acc1;

	// odd samples :: k=77
	__m256 v168 = data[168*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v164,v168));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v37,v40)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v165,v167)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v38,v39)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v166));
	__m256 v77 = data[77*stride];
	data[77*stride] = acc1;

	// even samples :: k=78
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v165,v168));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v38,v40)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v166,v167)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v39));
	__m256 v78 = data[78*stride];
	data[78*stride] = acc1;

	// odd samples :: k=79
	__m256 v169 = data[169*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v165,v169));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v38,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v166,v168)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v39,v40)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v167));
	__m256 v79 = data[79*stride];
	data[79*stride] = acc1;

	// even samples :: k=80
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v166,v169));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v39,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v167,v168)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v40));
	__m256 v80 = data[80*stride];
	data[80*stride] = acc1;

	// odd samples :: k=81
	__m256 v170 = data[170*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v166,v170));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v39,v42)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v167,v169)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v40,v41)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v168));
	__m256 v81 = data[81*stride];
	data[81*stride] = acc1;

	// even samples :: k=82
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v167,v170));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v40,v42)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v168,v169)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v41));
	__m256 v82 = data[82*stride];
	data[82*stride] = acc1;

	// odd samples :: k=83
	__m256 v171 = data[171*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v167,v171));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v40,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v168,v170)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v41,v42)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v169));
	__m256 v83 = data[83*stride];
	data[83*stride] = acc1;

	// even samples :: k=84
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v168,v171));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v41,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v169,v170)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v42));
	__m256 v84 = data[84*stride];
	data[84*stride] = acc1;

	// odd samples :: k=85
	__m256 v172 = data[172*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v168,v172));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v41,v44)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v169,v171)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v42,v43)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v170));
	__m256 v85 = data[85*stride];
	data[85*stride] = acc1;

	// even samples :: k=86
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v169,v172));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v42,v44)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v170,v171)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v43));
	__m256 v86 = data[86*stride];
	data[86*stride] = acc1;

	// odd samples :: k=87
	__m256 v173 = data[173*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v169,v173));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v42,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v170,v172)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v43,v44)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v171));
	__m256 v87 = data[87*stride];
	data[87*stride] = acc1;

	// even samples :: k=88
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v170,v173));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v43,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v171,v172)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v44));
	__m256 v88 = data[88*stride];
	data[88*stride] = acc1;

	// odd samples :: k=89
	__m256 v174 = data[174*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v170,v174));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v43,v46)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v171,v173)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v44,v45)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v172));
	__m256 v89 = data[89*stride];
	data[89*stride] = acc1;

	// even samples :: k=90
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v171,v174));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v44,v46)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v172,v173)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v45));
	__m256 v90 = data[90*stride];
	data[90*stride] = acc1;

	// odd samples :: k=91
	__m256 v175 = data[175*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v171,v175));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v44,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v172,v174)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v45,v46)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v173));
	__m256 v91 = data[91*stride];
	data[91*stride] = acc1;

	// even samples :: k=92
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v172,v175));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v45,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v173,v174)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v46));
	__m256 v92 = data[92*stride];
	data[92*stride] = acc1;

	// odd samples :: k=93
	__m256 v176 = data[176*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v172,v176));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v45,v48)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v173,v175)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v46,v47)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v174));
	__m256 v93 = data[93*stride];
	data[93*stride] = acc1;

	// even samples :: k=94
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v173,v176));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v46,v48)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v174,v175)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v47));
	__m256 v94 = data[94*stride];
	data[94*stride] = acc1;

	// odd samples :: k=95
	__m256 v177 = data[177*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v173,v177));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v46,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v174,v176)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v47,v48)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v175));
	__m256 v95 = data[95*stride];
	data[95*stride] = acc1;

	// even samples :: k=96
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v174,v177));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v47,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v175,v176)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v48));
	__m256 v96 = data[96*stride];
	data[96*stride] = acc1;

	// odd samples :: k=97
	__m256 v178 = data[178*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v174,v178));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v47,v50)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v175,v177)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v48,v49)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v176));
	__m256 v97 = data[97*stride];
	data[97*stride] = acc1;

	// even samples :: k=98
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v175,v178));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v48,v50)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v176,v177)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v49));
	__m256 v98 = data[98*stride];
	data[98*stride] = acc1;

	// odd samples :: k=99
	__m256 v179 = data[179*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v175,v179));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v48,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v176,v178)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v49,v50)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v177));
	__m256 v99 = data[99*stride];
	data[99*stride] = acc1;

	// even samples :: k=100
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v176,v179));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v49,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v177,v178)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v50));
	__m256 v100 = data[100*stride];
	data[100*stride] = acc1;

	// odd samples :: k=101
	__m256 v180 = data[180*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v176,v180));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v49,v52)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v177,v179)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v50,v51)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v178));
	__m256 v101 = data[101*stride];
	data[101*stride] = acc1;

	// even samples :: k=102
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v177,v180));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v50,v52)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v178,v179)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v51));
	__m256 v102 = data[102*stride];
	data[102*stride] = acc1;

	// odd samples :: k=103
	__m256 v181 = data[181*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v177,v181));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v50,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v178,v180)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v51,v52)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v179));
	__m256 v103 = data[103*stride];
	data[103*stride] = acc1;

	// even samples :: k=104
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v178,v181));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v51,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v179,v180)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v52));
	__m256 v104 = data[104*stride];
	data[104*stride] = acc1;

	// odd samples :: k=105
	__m256 v182 = data[182*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v178,v182));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v51,v54)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v179,v181)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v52,v53)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v180));
	__m256 v105 = data[105*stride];
	data[105*stride] = acc1;

	// even samples :: k=106
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v179,v182));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v52,v54)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v180,v181)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v53));
	__m256 v106 = data[106*stride];
	data[106*stride] = acc1;

	// odd samples :: k=107
	__m256 v183 = data[183*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v179,v183));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v52,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v180,v182)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v53,v54)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v181));
	__m256 v107 = data[107*stride];
	data[107*stride] = acc1;

	// even samples :: k=108
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v180,v183));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v53,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v181,v182)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v54));
	__m256 v108 = data[108*stride];
	data[108*stride] = acc1;

	// odd samples :: k=109
	__m256 v184 = data[184*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v180,v184));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v53,v56)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v181,v183)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v54,v55)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v182));
	__m256 v109 = data[109*stride];
	data[109*stride] = acc1;

	// even samples :: k=110
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v181,v184));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v54,v56)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v182,v183)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v55));
	__m256 v110 = data[110*stride];
	data[110*stride] = acc1;

	// odd samples :: k=111
	__m256 v185 = data[185*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v181,v185));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v54,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v182,v184)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v55,v56)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v183));
	__m256 v111 = data[111*stride];
	data[111*stride] = acc1;

	// even samples :: k=112
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v182,v185));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v55,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v183,v184)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v56));
	__m256 v112 = data[112*stride];
	data[112*stride] = acc1;

	// odd samples :: k=113
	__m256 v186 = data[186*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v182,v186));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v55,v58)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v183,v185)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v56,v57)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v184));
	__m256 v113 = data[113*stride];
	data[113*stride] = acc1;

	// even samples :: k=114
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v183,v186));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v56,v58)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v184,v185)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v57));
	__m256 v114 = data[114*stride];
	data[114*stride] = acc1;

	// odd samples :: k=115
	__m256 v187 = data[187*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v183,v187));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v56,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v184,v186)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v57,v58)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v185));
	__m256 v115 = data[115*stride];
	data[115*stride] = acc1;

	// even samples :: k=116
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v184,v187));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v57,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v185,v186)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v58));
	__m256 v116 = data[116*stride];
	data[116*stride] = acc1;

	// odd samples :: k=117
	__m256 v188 = data[188*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v184,v188));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v57,v60)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v185,v187)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v58,v59)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v186));
	__m256 v117 = data[117*stride];
	data[117*stride] = acc1;

	// even samples :: k=118
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v185,v188));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v58,v60)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v186,v187)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v59));
	__m256 v118 = data[118*stride];
	data[118*stride] = acc1;

	// odd samples :: k=119
	__m256 v189 = data[189*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v185,v189));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v58,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v186,v188)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v59,v60)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v187));
	__m256 v119 = data[119*stride];
	data[119*stride] = acc1;

	// even samples :: k=120
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v186,v189));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v59,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v187,v188)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v60));
	__m256 v120 = data[120*stride];
	data[120*stride] = acc1;

	// odd samples :: k=121
	__m256 v190 = data[190*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v186,v190));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v59,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v187,v189)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v60,v61)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v188));
	__m256 v121 = data[121*stride];
	data[121*stride] = acc1;

	// even samples :: k=122
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v187,v190));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v60,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v188,v189)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v61));
	__m256 v122 = data[122*stride];
	data[122*stride] = acc1;

	// odd samples :: k=123
	__m256 v191 = data[191*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v187,v191));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v60,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v188,v190)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v61,v62)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v189));
	__m256 v123 = data[123*stride];
	data[123*stride] = acc1;

	// even samples :: k=124
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v188,v191));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v61,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v189,v190)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v62));
	__m256 v124 = data[124*stride];
	data[124*stride] = acc1;

	// odd samples :: k=125
	__m256 v192 = data[192*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v188,v192));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v61,v64)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v189,v191)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v62,v63)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v190));
	__m256 v125 = data[125*stride];
	data[125*stride] = acc1;

	// even samples :: k=126
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v189,v192));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v62,v64)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v190,v191)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v63));
	__m256 v126 = data[126*stride];
	data[126*stride] = acc1;

	// odd samples :: k=127
	__m256 v193 = data[193*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v189,v193));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v62,v65)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v190,v192)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v63,v64)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v191));
	__m256 v127 = data[127*stride];
	data[127*stride] = acc1;

	// even samples :: k=128
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v190,v193));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v63,v65)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v191,v192)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v64));
	data[128*stride] = acc1;

	// odd samples :: k=129
	__m256 v194 = data[194*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v190,v194));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v63,v66)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v191,v193)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v64,v65)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v192));
	data[129*stride] = acc1;

	// even samples :: k=130
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v191,v194));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v64,v66)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v192,v193)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v65));
	data[130*stride] = acc1;

	// odd samples :: k=131
	__m256 v195 = data[195*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v191,v195));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v64,v67)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v192,v194)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v65,v66)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v193));
	data[131*stride] = acc1;

	// even samples :: k=132
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v192,v195));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v65,v67)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v193,v194)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v66));
	data[132*stride] = acc1;

	// odd samples :: k=133
	__m256 v196 = data[196*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v192,v196));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v65,v68)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v193,v195)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v66,v67)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v194));
	data[133*stride] = acc1;

	// even samples :: k=134
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v193,v196));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v66,v68)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v194,v195)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v67));
	data[134*stride] = acc1;

	// odd samples :: k=135
	__m256 v197 = data[197*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v193,v197));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v66,v69)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v194,v196)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v67,v68)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v195));
	data[135*stride] = acc1;

	// even samples :: k=136
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v194,v197));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v67,v69)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v195,v196)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v68));
	data[136*stride] = acc1;

	// odd samples :: k=137
	__m256 v198 = data[198*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v194,v198));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v67,v70)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v195,v197)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v68,v69)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v196));
	data[137*stride] = acc1;

	// even samples :: k=138
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v195,v198));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v68,v70)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v196,v197)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v69));
	data[138*stride] = acc1;

	// odd samples :: k=139
	__m256 v199 = data[199*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v195,v199));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v68,v71)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v196,v198)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v69,v70)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v197));
	data[139*stride] = acc1;

	// even samples :: k=140
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v196,v199));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v69,v71)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v197,v198)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v70));
	data[140*stride] = acc1;

	// odd samples :: k=141
	__m256 v200 = data[200*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v196,v200));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v69,v72)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v197,v199)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v70,v71)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v198));
	data[141*stride] = acc1;

	// even samples :: k=142
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v197,v200));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v70,v72)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v198,v199)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v71));
	data[142*stride] = acc1;

	// odd samples :: k=143
	__m256 v201 = data[201*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v197,v201));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v70,v73)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v198,v200)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v71,v72)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v199));
	data[143*stride] = acc1;

	// even samples :: k=144
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v198,v201));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v71,v73)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v199,v200)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v72));
	data[144*stride] = acc1;

	// odd samples :: k=145
	__m256 v202 = data[202*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v198,v202));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v71,v74)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v199,v201)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v72,v73)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v200));
	data[145*stride] = acc1;

	// even samples :: k=146
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v199,v202));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v72,v74)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v200,v201)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v73));
	data[146*stride] = acc1;

	// odd samples :: k=147
	__m256 v203 = data[203*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v199,v203));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v72,v75)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v200,v202)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v73,v74)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v201));
	data[147*stride] = acc1;

	// even samples :: k=148
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v200,v203));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v73,v75)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v201,v202)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v74));
	data[148*stride] = acc1;

	// odd samples :: k=149
	__m256 v204 = data[204*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v200,v204));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v73,v76)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v201,v203)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v74,v75)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v202));
	data[149*stride] = acc1;

	// even samples :: k=150
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v201,v204));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v74,v76)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v202,v203)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v75));
	data[150*stride] = acc1;

	// odd samples :: k=151
	__m256 v205 = data[205*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v201,v205));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v74,v77)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v202,v204)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v75,v76)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v203));
	data[151*stride] = acc1;

	// even samples :: k=152
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v202,v205));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v75,v77)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v203,v204)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v76));
	data[152*stride] = acc1;

	// odd samples :: k=153
	__m256 v206 = data[206*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v202,v206));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v75,v78)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v203,v205)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v76,v77)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v204));
	data[153*stride] = acc1;

	// even samples :: k=154
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v203,v206));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v76,v78)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v204,v205)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v77));
	data[154*stride] = acc1;

	// odd samples :: k=155
	__m256 v207 = data[207*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v203,v207));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v76,v79)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v204,v206)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v77,v78)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v205));
	data[155*stride] = acc1;

	// even samples :: k=156
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v204,v207));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v77,v79)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v205,v206)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v78));
	data[156*stride] = acc1;

	// odd samples :: k=157
	__m256 v208 = data[208*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v204,v208));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v77,v80)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v205,v207)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v78,v79)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v206));
	data[157*stride] = acc1;

	// even samples :: k=158
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v205,v208));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v78,v80)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v206,v207)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v79));
	data[158*stride] = acc1;

	// odd samples :: k=159
	__m256 v209 = data[209*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v205,v209));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v78,v81)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v206,v208)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v79,v80)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v207));
	data[159*stride] = acc1;

	// even samples :: k=160
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v206,v209));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v79,v81)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v207,v208)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v80));
	data[160*stride] = acc1;

	// odd samples :: k=161
	__m256 v210 = data[210*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v206,v210));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v79,v82)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v207,v209)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v80,v81)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v208));
	data[161*stride] = acc1;

	// even samples :: k=162
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v207,v210));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v80,v82)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v208,v209)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v81));
	data[162*stride] = acc1;

	// odd samples :: k=163
	__m256 v211 = data[211*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v207,v211));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v80,v83)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v208,v210)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v81,v82)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v209));
	data[163*stride] = acc1;

	// even samples :: k=164
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v208,v211));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v81,v83)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v209,v210)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v82));
	data[164*stride] = acc1;

	// odd samples :: k=165
	__m256 v212 = data[212*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v208,v212));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v81,v84)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v209,v211)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v82,v83)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v210));
	data[165*stride] = acc1;

	// even samples :: k=166
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v209,v212));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v82,v84)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v210,v211)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v83));
	data[166*stride] = acc1;

	// odd samples :: k=167
	__m256 v213 = data[213*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v209,v213));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v82,v85)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v210,v212)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v83,v84)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v211));
	data[167*stride] = acc1;

	// even samples :: k=168
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v210,v213));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v83,v85)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v211,v212)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v84));
	data[168*stride] = acc1;

	// odd samples :: k=169
	__m256 v214 = data[214*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v210,v214));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v83,v86)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v211,v213)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v84,v85)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v212));
	data[169*stride] = acc1;

	// even samples :: k=170
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v211,v214));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v84,v86)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v212,v213)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v85));
	data[170*stride] = acc1;

	// odd samples :: k=171
	__m256 v215 = data[215*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v211,v215));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v84,v87)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v212,v214)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v85,v86)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v213));
	data[171*stride] = acc1;

	// even samples :: k=172
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v212,v215));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v85,v87)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v213,v214)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v86));
	data[172*stride] = acc1;

	// odd samples :: k=173
	__m256 v216 = data[216*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v212,v216));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v85,v88)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v213,v215)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v86,v87)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v214));
	data[173*stride] = acc1;

	// even samples :: k=174
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v213,v216));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v86,v88)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v214,v215)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v87));
	data[174*stride] = acc1;

	// odd samples :: k=175
	__m256 v217 = data[217*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v213,v217));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v86,v89)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v214,v216)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v87,v88)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v215));
	data[175*stride] = acc1;

	// even samples :: k=176
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v214,v217));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v87,v89)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v215,v216)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v88));
	data[176*stride] = acc1;

	// odd samples :: k=177
	__m256 v218 = data[218*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v214,v218));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v87,v90)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v215,v217)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v88,v89)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v216));
	data[177*stride] = acc1;

	// even samples :: k=178
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v215,v218));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v88,v90)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v216,v217)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v89));
	data[178*stride] = acc1;

	// odd samples :: k=179
	__m256 v219 = data[219*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v215,v219));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v88,v91)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v216,v218)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v89,v90)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v217));
	data[179*stride] = acc1;

	// even samples :: k=180
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v216,v219));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v89,v91)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v217,v218)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v90));
	data[180*stride] = acc1;

	// odd samples :: k=181
	__m256 v220 = data[220*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v216,v220));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v89,v92)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v217,v219)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v90,v91)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v218));
	data[181*stride] = acc1;

	// even samples :: k=182
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v217,v220));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v90,v92)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v218,v219)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v91));
	data[182*stride] = acc1;

	// odd samples :: k=183
	__m256 v221 = data[221*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v217,v221));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v90,v93)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v218,v220)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v91,v92)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v219));
	data[183*stride] = acc1;

	// even samples :: k=184
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v218,v221));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v91,v93)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v219,v220)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v92));
	data[184*stride] = acc1;

	// odd samples :: k=185
	__m256 v222 = data[222*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v218,v222));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v91,v94)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v219,v221)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v92,v93)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v220));
	data[185*stride] = acc1;

	// even samples :: k=186
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v219,v222));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v92,v94)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v220,v221)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v93));
	data[186*stride] = acc1;

	// odd samples :: k=187
	__m256 v223 = data[223*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v219,v223));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v92,v95)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v220,v222)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v93,v94)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v221));
	data[187*stride] = acc1;

	// even samples :: k=188
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v220,v223));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v93,v95)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v221,v222)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v94));
	data[188*stride] = acc1;

	// odd samples :: k=189
	__m256 v224 = data[224*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v220,v224));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v93,v96)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v221,v223)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v94,v95)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v222));
	data[189*stride] = acc1;

	// even samples :: k=190
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v221,v224));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v94,v96)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v222,v223)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v95));
	data[190*stride] = acc1;

	// odd samples :: k=191
	__m256 v225 = data[225*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v221,v225));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v94,v97)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v222,v224)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v95,v96)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v223));
	data[191*stride] = acc1;

	// even samples :: k=192
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v222,v225));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v95,v97)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v223,v224)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v96));
	data[192*stride] = acc1;

	// odd samples :: k=193
	__m256 v226 = data[226*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v222,v226));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v95,v98)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v223,v225)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v96,v97)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v224));
	data[193*stride] = acc1;

	// even samples :: k=194
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v223,v226));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v96,v98)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v224,v225)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v97));
	data[194*stride] = acc1;

	// odd samples :: k=195
	__m256 v227 = data[227*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v223,v227));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v96,v99)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v224,v226)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v97,v98)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v225));
	data[195*stride] = acc1;

	// even samples :: k=196
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v224,v227));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v97,v99)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v225,v226)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v98));
	data[196*stride] = acc1;

	// odd samples :: k=197
	__m256 v228 = data[228*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v224,v228));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v97,v100)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v225,v227)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v98,v99)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v226));
	data[197*stride] = acc1;

	// even samples :: k=198
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v225,v228));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v98,v100)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v226,v227)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v99));
	data[198*stride] = acc1;

	// odd samples :: k=199
	__m256 v229 = data[229*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v225,v229));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v98,v101)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v226,v228)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v99,v100)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v227));
	data[199*stride] = acc1;

	// even samples :: k=200
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v226,v229));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v99,v101)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v227,v228)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v100));
	data[200*stride] = acc1;

	// odd samples :: k=201
	__m256 v230 = data[230*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v226,v230));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v99,v102)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v227,v229)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v100,v101)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v228));
	data[201*stride] = acc1;

	// even samples :: k=202
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v227,v230));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v100,v102)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v228,v229)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v101));
	data[202*stride] = acc1;

	// odd samples :: k=203
	__m256 v231 = data[231*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v227,v231));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v100,v103)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v228,v230)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v101,v102)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v229));
	data[203*stride] = acc1;

	// even samples :: k=204
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v228,v231));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v101,v103)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v229,v230)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v102));
	data[204*stride] = acc1;

	// odd samples :: k=205
	__m256 v232 = data[232*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v228,v232));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v101,v104)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v229,v231)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v102,v103)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v230));
	data[205*stride] = acc1;

	// even samples :: k=206
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v229,v232));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v102,v104)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v230,v231)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v103));
	data[206*stride] = acc1;

	// odd samples :: k=207
	__m256 v233 = data[233*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v229,v233));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v102,v105)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v230,v232)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v103,v104)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v231));
	data[207*stride] = acc1;

	// even samples :: k=208
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v230,v233));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v103,v105)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v231,v232)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v104));
	data[208*stride] = acc1;

	// odd samples :: k=209
	__m256 v234 = data[234*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v230,v234));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v103,v106)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v231,v233)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v104,v105)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v232));
	data[209*stride] = acc1;

	// even samples :: k=210
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v231,v234));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v104,v106)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v232,v233)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v105));
	data[210*stride] = acc1;

	// odd samples :: k=211
	__m256 v235 = data[235*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v231,v235));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v104,v107)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v232,v234)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v105,v106)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v233));
	data[211*stride] = acc1;

	// even samples :: k=212
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v232,v235));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v105,v107)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v233,v234)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v106));
	data[212*stride] = acc1;

	// odd samples :: k=213
	__m256 v236 = data[236*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v232,v236));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v105,v108)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v233,v235)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v106,v107)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v234));
	data[213*stride] = acc1;

	// even samples :: k=214
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v233,v236));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v106,v108)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v234,v235)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v107));
	data[214*stride] = acc1;

	// odd samples :: k=215
	__m256 v237 = data[237*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v233,v237));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v106,v109)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v234,v236)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v107,v108)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v235));
	data[215*stride] = acc1;

	// even samples :: k=216
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v234,v237));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v107,v109)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v235,v236)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v108));
	data[216*stride] = acc1;

	// odd samples :: k=217
	__m256 v238 = data[238*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v234,v238));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v107,v110)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v235,v237)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v108,v109)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v236));
	data[217*stride] = acc1;

	// even samples :: k=218
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v235,v238));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v108,v110)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v236,v237)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v109));
	data[218*stride] = acc1;

	// odd samples :: k=219
	__m256 v239 = data[239*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v235,v239));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v108,v111)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v236,v238)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v109,v110)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v237));
	data[219*stride] = acc1;

	// even samples :: k=220
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v236,v239));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v109,v111)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v237,v238)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v110));
	data[220*stride] = acc1;

	// odd samples :: k=221
	__m256 v240 = data[240*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v236,v240));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v109,v112)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v237,v239)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v110,v111)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v238));
	data[221*stride] = acc1;

	// even samples :: k=222
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v237,v240));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v110,v112)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v238,v239)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v111));
	data[222*stride] = acc1;

	// odd samples :: k=223
	__m256 v241 = data[241*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v237,v241));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v110,v113)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v238,v240)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v111,v112)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v239));
	data[223*stride] = acc1;

	// even samples :: k=224
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v238,v241));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v111,v113)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v239,v240)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v112));
	data[224*stride] = acc1;

	// odd samples :: k=225
	__m256 v242 = data[242*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v238,v242));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v111,v114)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v239,v241)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v112,v113)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v240));
	data[225*stride] = acc1;

	// even samples :: k=226
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v239,v242));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v112,v114)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v240,v241)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v113));
	data[226*stride] = acc1;

	// odd samples :: k=227
	__m256 v243 = data[243*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v239,v243));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v112,v115)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v240,v242)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v113,v114)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v241));
	data[227*stride] = acc1;

	// even samples :: k=228
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v240,v243));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v113,v115)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v241,v242)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v114));
	data[228*stride] = acc1;

	// odd samples :: k=229
	__m256 v244 = data[244*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v240,v244));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v113,v116)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v241,v243)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v114,v115)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v242));
	data[229*stride] = acc1;

	// even samples :: k=230
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v241,v244));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v114,v116)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v242,v243)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v115));
	data[230*stride] = acc1;

	// odd samples :: k=231
	__m256 v245 = data[245*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v241,v245));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v114,v117)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v242,v244)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v115,v116)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v243));
	data[231*stride] = acc1;

	// even samples :: k=232
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v242,v245));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v115,v117)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v243,v244)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v116));
	data[232*stride] = acc1;

	// odd samples :: k=233
	__m256 v246 = data[246*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v242,v246));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v115,v118)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v243,v245)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v116,v117)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v244));
	data[233*stride] = acc1;

	// even samples :: k=234
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v243,v246));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v116,v118)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v244,v245)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v117));
	data[234*stride] = acc1;

	// odd samples :: k=235
	__m256 v247 = data[247*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v243,v247));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v116,v119)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v244,v246)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v117,v118)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v245));
	data[235*stride] = acc1;

	// even samples :: k=236
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v244,v247));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v117,v119)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v245,v246)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v118));
	data[236*stride] = acc1;

	// odd samples :: k=237
	__m256 v248 = data[248*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v244,v248));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v117,v120)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v245,v247)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v118,v119)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v246));
	data[237*stride] = acc1;

	// even samples :: k=238
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v245,v248));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v118,v120)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v246,v247)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v119));
	data[238*stride] = acc1;

	// odd samples :: k=239
	__m256 v249 = data[249*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v245,v249));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v118,v121)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v246,v248)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v119,v120)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v247));
	data[239*stride] = acc1;

	// even samples :: k=240
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v246,v249));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v119,v121)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v247,v248)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v120));
	data[240*stride] = acc1;

	// odd samples :: k=241
	__m256 v250 = data[250*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v246,v250));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v119,v122)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v247,v249)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v120,v121)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v248));
	data[241*stride] = acc1;

	// even samples :: k=242
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v247,v250));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v120,v122)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v248,v249)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v121));
	data[242*stride] = acc1;

	// odd samples :: k=243
	__m256 v251 = data[251*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v247,v251));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v120,v123)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v248,v250)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v121,v122)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v249));
	data[243*stride] = acc1;

	// even samples :: k=244
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v248,v251));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v121,v123)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v249,v250)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v122));
	data[244*stride] = acc1;

	// odd samples :: k=245
	__m256 v252 = data[252*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v248,v252));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v121,v124)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v249,v251)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v122,v123)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v250));
	data[245*stride] = acc1;

	// even samples :: k=246
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v249,v252));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v122,v124)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v250,v251)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v123));
	data[246*stride] = acc1;

	// odd samples :: k=247
	__m256 v253 = data[253*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v249,v253));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v122,v125)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v250,v252)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v123,v124)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v251));
	data[247*stride] = acc1;

	// even samples :: k=248
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v250,v253));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v123,v125)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v251,v252)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v124));
	data[248*stride] = acc1;

	// odd samples :: k=249
	__m256 v254 = data[254*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v250,v254));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v123,v126)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v251,v253)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v124,v125)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v252));
	data[249*stride] = acc1;

	// even samples :: k=250
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v251,v254));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v124,v126)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v252,v253)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v125));
	data[250*stride] = acc1;

	// odd samples :: k=251
	__m256 v255 = data[255*stride];
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v251,v255));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v124,v127)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v252,v254)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v125,v126)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v253));
	data[251*stride] = acc1;

	// even samples :: k=252
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v252,v255));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v125,v127)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v253,v254)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v126));
	data[252*stride] = acc1;

	// odd samples :: k=253
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v252,v254));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v125,v127)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v253,v255)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v126,v127)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v254));
	data[253*stride] = acc1;

	// even samples :: k=254
	acc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v253,v254));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v126,v127)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v254,v255)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v127));
	data[254*stride] = acc1;

	// odd samples :: k=255
	acc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v253,v253));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v126,v126)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v254,v254)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v127,v127)));
	acc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v255));
	data[255*stride] = acc1;
}

#endif
