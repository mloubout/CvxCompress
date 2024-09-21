#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

/*
 *   This is a rewrite of the ChvCompress code by Ergas et.al.
 *   The purpose of this rewrite is to understand the algorithm.
 *   No attempt has been made at making this code fast.
 */

/*
 *   The following defines establish the wavelet filters used here. The
 *   coefficients are those of the 7-9 tap filters used in the FBI specification,
 *   originally defined by Antonini, and verified by Villasenor to be optimal for
 *   image compression.
 *                                 
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

/*
 * Forward wavelet transform mirrors inputs around 0 and n-1 in the following way:
 *
 * n+1 -> n-3
 * n   -> n-2 
 * n-1
 * n-2
 * n-3
 * .
 * .
 * 2
 * 1
 * 0
 * -1 -> 1
 * -2 -> 2
 *
 * At the lower end, val = abs(val)
 * At the upper end, val = 2*n - 2 - val
 * Note that after lower end mirroring, val may still exceed upper bound.
 * Vice versa, val may be negative after upper end mirroring.
 * To protect against this, the mirroring must be interleaved and chained like this:
 *
 * val = abs(val)
 * val = 2*n - 2 - val
 * val = abs(val)
 * val = 2*n - 2 - val
 */
inline int MIRR(int inp_val, int dim)
{
	int val = inp_val < 0 ? -inp_val : inp_val;
	val = (val >= dim) ? (2*dim-2-val) : val;
	val = val < 0 ? -val : val;
	val = (val >= dim) ? (2*dim-2-val) : val;
	//printf("  -> -> MIRR(%d,%d) = %d\n",inp_val,dim,val);
	return val;
}

static bool Verbose = false;

inline void
Ds79(
	float* p_in,
	float* p_tmp,
	int stride,
	int dim
	)
{
	if (Verbose) printf("Ds79(*,*,stride=%d,dim=%d)\n",stride,dim);
	for (int n = dim;  n >= 2;  n = n-n/2)
	{
		if (Verbose)
		{
			printf("Ds d_inp = %e",p_in[0]);
			for (int i = 1;  i < n;  ++i) printf(", %e",p_in[i*stride]);
			printf("\n");
		}

		// copy inputs to tmp buffer, p_in will be overwritten
		for (int i = 0;  i < n;  ++i) p_tmp[i] = p_in[i*stride];

		int nh = n / 2;
		int nl = n - nh;
		//printf("  -> n=%d, nh=%d, nl=%d\n",n,nh,nl);
		for (int ix = 0;  ix < nl;  ++ix)
		{
			int i0 = 2 * ix;
			int im1 = MIRR(i0-1,n);  int ip1 = MIRR(i0+1,n);
			int im2 = MIRR(i0-2,n);  int ip2 = MIRR(i0+2,n);
			int im3 = MIRR(i0-3,n);  int ip3 = MIRR(i0+3,n);
			int im4 = MIRR(i0-4,n);  int ip4 = MIRR(i0+4,n);
			//printf("d[%d] = al0 * t[%d] + al1 * (t[%d] + t[%d]) + al2 * (t[%d] + t[%d]) + al3 * (t[%d] + t[%d]) + al4 * (t[%d] + t[%d])\n",ix,i0,im1,ip1,im2,ip2,im3,ip3,im4,ip4);
			// sum smallest to largest (most accurate way of summing floats)
			float acc1 = al4 * (p_tmp[im4] + p_tmp[ip4]);
			acc1 += al1 * (p_tmp[im1] + p_tmp[ip1]);
			acc1 += al0 * p_tmp[i0];
			float acc2 = al3 * (p_tmp[im3] + p_tmp[ip3]);
			acc2 += al2 * (p_tmp[im2] + p_tmp[ip2]);
			p_in[ix*stride] = acc1 + acc2;
		}
		for (int ix = 0;  ix < nh;  ++ix)
		{
			int i0 = 2 * ix + 1;
			int im1 = MIRR(i0-1,n);  int ip1 = MIRR(i0+1,n);
			int im2 = MIRR(i0-2,n);  int ip2 = MIRR(i0+2,n);
			int im3 = MIRR(i0-3,n);  int ip3 = MIRR(i0+3,n);
			//printf("d[%d] = ah0 * t[%d] + ah1 * (t[%d] + t[%d]) + ah2 * (t[%d] + t[%d]) + ah3 * (t[%d] + t[%d])\n",nl+ix,i0,im1,ip1,im2,ip2,im3,ip3);
			// sum smallest to largest (most accurate way of summing floats)
			float acc1 = ah3 * (p_tmp[im3] + p_tmp[ip3]);
			acc1 += ah0 * p_tmp[i0];
			float acc2 = ah2 * (p_tmp[im2] + p_tmp[ip2]);
			acc2 += ah1 * (p_tmp[im1] + p_tmp[ip1]);
			p_in[(nl+ix)*stride] = acc1 + acc2;
		}

		if (Verbose)
		{
			printf("Ds d_out = %e",p_in[0]);
			for (int i = 1;  i < n;  ++i) printf(", %e",p_in[i*stride]);
			printf("\n");
		}
	}
	if (Verbose) printf("\n");
}

#define sl0  7.884856164056601e-001f
#define sl1  4.180922732222101e-001f
#define sl2 -4.068941760955800e-002f
#define sl3 -6.453888262893799e-002f

#define sh0  8.526986790094000e-001f
#define sh1 -3.774028556126500e-001f
#define sh2 -1.106244044184200e-001f
#define sh3  2.384946501938001e-002f
#define sh4  3.782845550699501e-002f

/* 
 * Inverse wavelet transform mirroring works like this:
 *
 * Mirroring is done separate for SL and SH coefficients.
 *
 * For SL coefficients, mirroring works like this:
 * nl+1 -> nl-2
 * nl   -> nl-1
 * nl-1
 * nl-2
 * .
 * .
 * 2
 * 1
 * 0
 * -1 -> 1
 * -2 -> 2
 *
 * For SH coefficients, mirroring works like this:
 * n+1 -> n-3
 * n   -> n-2
 * n-1
 * n-2
 * .
 * .
 * nl+2
 * nl+1
 * nl
 * nl-1 -> nl
 * nl-2 -> nl+1
 */
inline int MIRR_SL(int inp_val, int nl)
{
	int val = inp_val;
	val = val < 0 ? -val : val;
	val = (val >= nl) ? (2*nl-1-val) : val;
	val = val < 0 ? -val : val;
	val = (val >= nl) ? (2*nl-1-val) : val;
	val = val < 0 ? -val : val;
	val = (val >= nl) ? (2*nl-1-val) : val;
	return val;
}
inline int MIRR_SH(int inp_val, int nl, int nh)
{
	int val = inp_val - nl;
	val = val < 0 ? -val-1 : val;
	val = (val >= nh) ? (2*nh-2-val) : val;
	val = val < 0 ? -val-1 : val;
	val = (val >= nh) ? (2*nh-2-val) : val;
	val = val < 0 ? -val-1 : val;
	val = (val >= nh) ? (2*nh-2-val) : val;
	return nl + val;
}

inline void
Us79(
	float* p_in,
	float* t,
	int stride,
	int dim
	)
{
	if (Verbose) printf("Us79(*,*,stride=%d,dim=%d)\n",stride,dim);
	int* l = new int[dim];
	int nx = 0;
	for (int n = dim;  n >= 2;  n = n-n/2) {l[nx++] = n;}
	for (int li = nx-1;  li >= 0;  --li)
	{
		int n = l[li];

		if (Verbose)
		{
			printf("Us d_inp = %e",p_in[0]);
			for (int i = 1;  i < n;  ++i) printf(", %e",p_in[i*stride]);
			printf("\n");
		}

		// copy inputs to tmp buffer, p_in will be overwritten
		for (int i = 0;  i < n;  ++i) t[i] = p_in[i*stride];

		int nh = n / 2;
		int nl = n - nh;
		//printf("  -> n=%d, nh=%d, nl=%d\n",n,nh,nl);
		for (int k = 0;  k < nl;  ++k)
		{
			if (Verbose) printf("d[%d] = sl0*t[%d] + sl2*(t[%d]+t[%d]) + sh1*(t[%d]+t[%d]) + sh3*(t[%d]+t[%d])\n",2*k,k,MIRR_SL(k-1,nl),MIRR_SL(k+1,nl),MIRR_SH(nl+k-1,nl,nh),MIRR_SH(nl+k,nl,nh),MIRR_SH(nl+k-2,nl,nh),MIRR_SH(nl+k+1,nl,nh));
			p_in[2*k*stride] = 
				sl0 * t[k] + 
				sl2 * ( t[MIRR_SL(k-1,nl)] + t[MIRR_SL(k+1,nl)] ) +
				sh1 * ( t[MIRR_SH(nl+k-1,nl,nh)] + t[MIRR_SH(nl+k,nl,nh)] ) +
				sh3 * ( t[MIRR_SH(nl+k-2,nl,nh)] + t[MIRR_SH(nl+k+1,nl,nh)] );
		}
		for (int k = 0;  k < nh;  ++k)
		{
			if (Verbose) printf("d[%d] = sl1*(t[%d]+t[%d]) + sl3*(t[%d]+t[%d]) + sh0*t[%d] + sh2*(t[%d]+t[%d]) + sh4*(t[%d]+t[%d])\n",(2*k+1),MIRR_SL(k,nl),MIRR_SL(k+1,nl),MIRR_SL(k-1,nl),MIRR_SL(k+2,nl),nl+k,MIRR_SH(nl+k-1,nl,nh),MIRR_SH(nl+k+1,nl,nh),MIRR_SH(nl+k-2,nl,nh),MIRR_SH(nl+k+2,nl,nh));
			p_in[(2*k+1)*stride] = 
				sl1 * ( t[MIRR_SL(k,nl)] + t[MIRR_SL(k+1,nl)] ) +
				sl3 * ( t[MIRR_SL(k-1,nl)] + t[MIRR_SL(k+2,nl)] ) +
				sh0 * t[nl+k] +
				sh2 * ( t[MIRR_SH(nl+k-1,nl,nh)] + t[MIRR_SH(nl+k+1,nl,nh)] ) +
				sh4 * ( t[MIRR_SH(nl+k-2,nl,nh)] + t[MIRR_SH(nl+k+2,nl,nh)] );
		}

		if (Verbose)
		{
			printf("Us d_out = %e",p_in[0]);
			for (int i = 1;  i < n;  ++i) printf(", %e",p_in[i*stride]);
			printf("\n");
		}
	}
	delete [] l;
	if (Verbose) printf("\n");
}

void Wavelet_Transform_Slow_Forward(
	float* data,
	float* work,
	int bx,
	int by,
	int bz,
	int x0,
	int y0,
	int z0,
	int nx,
	int ny,
	int nz
	)
{
	for (int iz = 0;  iz < bz;  ++iz) {
		if (bx > 1) for (int iy = 0;  iy < by;  ++iy) Ds79(data+((iz+z0)*ny+(iy+y0))*nx+(x0), work,  1, bx);
		if (by > 1) for (int ix = 0;  ix < bx;  ++ix) Ds79(data+((iz+z0)*ny+(y0))*nx+(ix+x0), work, bx, by);
	}
	if (bz > 1) for (int iy = 0;  iy < by;  ++iy) for (int ix = 0;  ix < bx;  ++ix) Ds79(data+((z0)*ny+(iy+y0))*nx+(ix+x0), work, bx*by, bz);
}

void Wavelet_Transform_Slow_Inverse(
	float* data,
	float* work,
	int bx,
	int by,
	int bz,
	int x0,
	int y0,
	int z0,
	int nx,
	int ny,
	int nz
	)
{
	for (int iz = 0;  iz < bz;  ++iz) {
		if (bx > 1) for (int iy = 0;  iy < by;  ++iy) Us79(data+((iz+z0)*ny+(iy+y0))*nx+(x0), work,  1, bx);
		if (by > 1) for (int ix = 0;  ix < bx;  ++ix) Us79(data+((iz+z0)*ny+(y0))*nx+(ix+x0), work, bx, by);
	}
	if (bz > 1) for (int iy = 0;  iy < by;  ++iy) for (int ix = 0;  ix < bx;  ++ix) Us79(data+((z0)*ny+(iy+y0))*nx+(ix+x0), work, bx*by, bz);
}

//
// Code generator. Generates the base AVX and AVX2 implementations of the wavelet forward and inverse transforms.
// Only powers of two are supported for the time being. Would be easy to extend to arbitrary lengths, but this
// wasn't necessary for the compression library.
//

int Find_Index(int* var_prev_idx, int prev_n, int idx)
{
	if (idx >= 0)
		for (int i = 0;  i < prev_n;  ++i)
			if (var_prev_idx[i] == idx)
				return i;
	return -9999;
}

void Print_Load_Line(
		FILE* fp,
		int idx,
		bool tmp_array,
		int* var_prev_idx,
		int num_vars,
		int* var_curr_idx,
		int num_curr
		)
{
	if (tmp_array)
	{
		// count number of variables that need to be loaded
		if (Find_Index(var_prev_idx, num_vars, idx) < 0)
		{
			// find variables that are no longer needed
			for (int i = 0;  i < num_vars;  ++i)
			{
				if (Find_Index(var_curr_idx,num_curr,var_prev_idx[i]) < 0)
				{
					var_prev_idx[i] = idx;
					fprintf(fp,"\tv%d = tmp[%d];\n",i,idx);
					break;
				}
			}
		}
	}
	else
	{
		if (var_prev_idx[idx] < 0)
		{
			var_prev_idx[idx] = idx;
			fprintf(fp,"\t__m256 v%d = data[%d*stride];\n",idx,idx);
		}
	}
}

void Gen_Ds79_Core(FILE* fp, int n, int num_vars, bool avx2)
{
        if (num_vars < 9)
        {
                fprintf(fp,"Warning! num_vars must be 9 or larger.\n");
                fprintf(fp,"         Changing num_vars to 9.\n");
                num_vars = 9;
        }
        bool tmp_array = num_vars < n ? true : false;

	if (n <= 32) fprintf(fp,"static inline "); else fprintf(fp,"static ");
        if (tmp_array)
                fprintf(fp,"void _Ds79_AVX_%d(__m256* data, __m256* tmp, int stride)\n",n);
        else
                fprintf(fp,"void _Ds79_AVX_%d(__m256* data, int stride)\n",n);
        fprintf(fp,"{\n");

	if (tmp_array)
	{
		// copy inputs to tmp buffer, p_in will be overwritten
		fprintf(fp,"\t// copy inputs to tmp buffer\n");
		for (int i = 0;  i < n;  ++i) fprintf(fp,"\ttmp[%d] = data[%d*stride];\n",i,i);
		fprintf(fp,"\n");

		for (int j = 0;  j < (num_vars+7)/8;  ++j)
		{
			fprintf(fp,"\t__m256 ");
			for (int i = j*8;  i < (j+1)*8 && i < num_vars;  ++i)
			{
				if (i == j*8) fprintf(fp,"v%d",i);
				else fprintf(fp,",v%d",i);
			}
			fprintf(fp,";\n");
		}
	}

	int nh = n / 2;
	int nl = n - nh;
	int* var_prev_idx = new int[num_vars];
	for (int i = 0;  i < num_vars;  ++i) var_prev_idx[i] = -1;
	int var_curr_idx[9];
	fprintf(fp,"\t__m256 acc1;\n");
	for (int ix = 0;  ix < nl;  ++ix)
	{
		{
			fprintf(fp,"\n\t// lower band :: ix=%d\n",ix);
			int i0 = 2 * ix;
			int im1 = MIRR(i0-1,n);  int ip1 = MIRR(i0+1,n);
			int im2 = MIRR(i0-2,n);  int ip2 = MIRR(i0+2,n);
			int im3 = MIRR(i0-3,n);  int ip3 = MIRR(i0+3,n);
			int im4 = MIRR(i0-4,n);  int ip4 = MIRR(i0+4,n);
			var_curr_idx[0] = im4;
			var_curr_idx[1] = im3;
			var_curr_idx[2] = im2;
			var_curr_idx[3] = im1;
			var_curr_idx[4] = i0;
			var_curr_idx[5] = ip1;
			var_curr_idx[6] = ip2;
			var_curr_idx[7] = ip3;
			var_curr_idx[8] = ip4;

			Print_Load_Line(fp,im4,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
			Print_Load_Line(fp,ip4,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
			fprintf(fp,"\tacc1 = _mm256_mul_ps(_mm_al4,_mm256_add_ps(v%d,v%d));\n",Find_Index(var_prev_idx,num_vars,im4),Find_Index(var_prev_idx,num_vars,ip4));

                        Print_Load_Line(fp,im3,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
                        Print_Load_Line(fp,ip3,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
			if (avx2)
				fprintf(fp,"\tacc1 = _mm256_fmadd_ps(_mm_al3,_mm256_add_ps(v%d,v%d),acc1);\n",Find_Index(var_prev_idx,num_vars,im3),Find_Index(var_prev_idx,num_vars,ip3));
			else
				fprintf(fp,"\tacc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al3,_mm256_add_ps(v%d,v%d)));\n",Find_Index(var_prev_idx,num_vars,im3),Find_Index(var_prev_idx,num_vars,ip3));

                        Print_Load_Line(fp,im2,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
                        Print_Load_Line(fp,ip2,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
                        if (avx2)
                                fprintf(fp,"\tacc1 = _mm256_fmadd_ps(_mm_al2,_mm256_add_ps(v%d,v%d),acc1);\n",Find_Index(var_prev_idx,num_vars,im2),Find_Index(var_prev_idx,num_vars,ip2));
                        else
                                fprintf(fp,"\tacc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al2,_mm256_add_ps(v%d,v%d)));\n",Find_Index(var_prev_idx,num_vars,im2),Find_Index(var_prev_idx,num_vars,ip2));

			Print_Load_Line(fp,im1,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
			Print_Load_Line(fp,ip1,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
			if (avx2)
				fprintf(fp,"\tacc1 = _mm256_fmadd_ps(_mm_al1,_mm256_add_ps(v%d,v%d),acc1);\n",Find_Index(var_prev_idx,num_vars,im1),Find_Index(var_prev_idx,num_vars,ip1));
                        else
                                fprintf(fp,"\tacc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al1,_mm256_add_ps(v%d,v%d)));\n",Find_Index(var_prev_idx,num_vars,im1),Find_Index(var_prev_idx,num_vars,ip1));

                        Print_Load_Line(fp,i0,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
                        if (avx2)
                                fprintf(fp,"\tacc1 = _mm256_fmadd_ps(_mm_al0, v%d,acc1);\n",Find_Index(var_prev_idx,num_vars,i0));
                        else
                                fprintf(fp,"\tacc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_al0,v%d));\n",Find_Index(var_prev_idx,num_vars,i0));

                        if (!tmp_array) Print_Load_Line(fp,ix,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
                        fprintf(fp,"\tdata[%d*stride] = acc1;\n",ix);
                }

                {
                        fprintf(fp,"\n\t// upper band :: ix=%d\n",ix);
                        int i0 = 2 * ix + 1;
                        int im1 = MIRR(i0-1,n);  int ip1 = MIRR(i0+1,n);
                        int im2 = MIRR(i0-2,n);  int ip2 = MIRR(i0+2,n);
                        int im3 = MIRR(i0-3,n);  int ip3 = MIRR(i0+3,n);
                        var_curr_idx[0] = im3;
                        var_curr_idx[1] = im2;
                        var_curr_idx[2] = im1;
                        var_curr_idx[3] = i0;
                        var_curr_idx[4] = ip1;
                        var_curr_idx[5] = ip2;
                        var_curr_idx[6] = ip3;

                        Print_Load_Line(fp,im3,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
                        Print_Load_Line(fp,ip3,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
                        fprintf(fp,"\tacc1 = _mm256_mul_ps(_mm_ah3,_mm256_add_ps(v%d,v%d));\n",Find_Index(var_prev_idx,num_vars,im3),Find_Index(var_prev_idx,num_vars,ip3));

                        Print_Load_Line(fp,im2,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
                        Print_Load_Line(fp,ip2,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
			if (avx2)
				fprintf(fp,"\tacc1 = _mm256_fmadd_ps(_mm_ah2,_mm256_add_ps(v%d,v%d),acc1);\n",Find_Index(var_prev_idx,num_vars,im2),Find_Index(var_prev_idx,num_vars,ip2));
			else
				fprintf(fp,"\tacc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah2,_mm256_add_ps(v%d,v%d)));\n",Find_Index(var_prev_idx,num_vars,im2),Find_Index(var_prev_idx,num_vars,ip2));

                        Print_Load_Line(fp,im1,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
                        Print_Load_Line(fp,ip1,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
                        if (avx2)
                                fprintf(fp,"\tacc1 = _mm256_fmadd_ps(_mm_ah1,_mm256_add_ps(v%d,v%d),acc1);\n",Find_Index(var_prev_idx,num_vars,im1),Find_Index(var_prev_idx,num_vars,ip1));
                        else
                                fprintf(fp,"\tacc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah1,_mm256_add_ps(v%d,v%d)));\n",Find_Index(var_prev_idx,num_vars,im1),Find_Index(var_prev_idx,num_vars,ip1));

                        Print_Load_Line(fp,i0,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
                        if (avx2)
                                fprintf(fp,"\tacc1 = _mm256_fmadd_ps(_mm_ah0,v%d,acc1);\n",Find_Index(var_prev_idx,num_vars,i0));
                        else
                                fprintf(fp,"\tacc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_ah0,v%d));\n",Find_Index(var_prev_idx,num_vars,i0));

                        if (!tmp_array) Print_Load_Line(fp,nl+ix,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
                        fprintf(fp,"\tdata[%d*stride] = acc1;\n",nl+ix);
                }
        }
        fprintf(fp,"}\n\n");
}

void Gen_Ds79(const char* path, int min_n, int max_n, int num_vars)
{
	FILE* fp = fopen(path, "w");

	fprintf(fp,"/*!\n");
	fprintf(fp," * Don't edit this code, it was automatically generated.\n");
	fprintf(fp," * Base functions for wavelet transforms of length %d to %d.\n",1<<min_n,1<<max_n);
	fprintf(fp," */\n\n");
	fprintf(fp,"#define  SIMDE_ENABLE_NATIVE_ALIASES \n");
	fprintf(fp,"#include <simde/x86/avx.h>  // AVX intrinsics\n\n");

	fprintf(fp,"/*\n");
	fprintf(fp," * Define coefficients for Antonini 7-9 tap filter.\n");
	fprintf(fp," */\n");
	fprintf(fp,"#define al0  8.526986790094000e-001f\n");
	fprintf(fp,"#define al1  3.774028556126500e-001f\n");
	fprintf(fp,"#define al2 -1.106244044184200e-001f\n");
	fprintf(fp,"#define al3 -2.384946501938001e-002f\n");
	fprintf(fp,"#define al4  3.782845550699501e-002f\n\n");

	fprintf(fp,"#define ah0  7.884856164056601e-001f\n");
	fprintf(fp,"#define ah1 -4.180922732222101e-001f\n");
	fprintf(fp,"#define ah2 -4.068941760955800e-002f\n");
	fprintf(fp,"#define ah3  6.453888262893799e-002f\n\n");

	fprintf(fp,"#define _mm_al0 _mm256_set1_ps(al0)\n");
	fprintf(fp,"#define _mm_al1 _mm256_set1_ps(al1)\n");
	fprintf(fp,"#define _mm_al2 _mm256_set1_ps(al2)\n");
	fprintf(fp,"#define _mm_al3 _mm256_set1_ps(al3)\n");
	fprintf(fp,"#define _mm_al4 _mm256_set1_ps(al4)\n\n");

	fprintf(fp,"#define _mm_ah0 _mm256_set1_ps(ah0)\n");
	fprintf(fp,"#define _mm_ah1 _mm256_set1_ps(ah1)\n");
	fprintf(fp,"#define _mm_ah2 _mm256_set1_ps(ah2)\n");
	fprintf(fp,"#define _mm_ah3 _mm256_set1_ps(ah3)\n\n");

	fprintf(fp,"#ifdef __AVX2__\n\n");

	for (int i = min_n;  i <= max_n;  ++i) Gen_Ds79_Core(fp,1<<i,num_vars,true);

	fprintf(fp,"#else\n\n");

	for (int i = min_n;  i <= max_n;  ++i) Gen_Ds79_Core(fp,1<<i,num_vars,false);

	fprintf(fp,"#endif\n");

	fclose(fp);
	
	printf("Wrote Ds79 base code to file %s.\n",path);
}

void Gen_Us79_Core(FILE* fp, int n, int num_vars, bool avx2)
{
        if (num_vars < 9)
        {
                fprintf(fp,"Warning! num_vars must be 9 or larger.\n");
                fprintf(fp,"         Changing num_vars to 9.\n");
                num_vars = 9;
        }
        bool tmp_array = num_vars < n ? true : false;

	if (n <= 32) fprintf(fp,"static inline "); else fprintf(fp,"static ");
        if (tmp_array)
                fprintf(fp,"void _Us79_AVX_%d(__m256* data, __m256* tmp, int stride)\n",n);
        else
                fprintf(fp,"void _Us79_AVX_%d(__m256* data, int stride)\n",n);
        fprintf(fp,"{\n");

	if (tmp_array)
	{
		// copy inputs to tmp buffer, p_in will be overwritten
		fprintf(fp,"\t// copy inputs to tmp buffer\n");
		for (int i = 0;  i < n;  ++i) fprintf(fp,"\ttmp[%d] = data[%d*stride];\n",i,i);
		fprintf(fp,"\n");

		for (int j = 0;  j < (num_vars+7)/8;  ++j)
		{
			fprintf(fp,"\t__m256 ");
			for (int i = j*8;  i < (j+1)*8 && i < num_vars;  ++i)
			{
				if (i == j*8) fprintf(fp,"v%d",i);
				else fprintf(fp,",v%d",i);
			}
			fprintf(fp,";\n");
		}
	}

	int nh = n / 2;
	int nl = n - nh;
	int* var_prev_idx = new int[num_vars];
	for (int i = 0;  i < num_vars;  ++i) var_prev_idx[i] = -1;
	int var_curr_idx[9];
	fprintf(fp,"\t__m256 acc1;\n");
	for (int k = 0;  k < nl;  ++k)
	{
		{
			fprintf(fp,"\n\t// even samples :: k=%d\n",2*k);
			int i0 = k;
			int im1 = MIRR_SH(nl+k-1,nl,nh);  int ip1 = MIRR_SH(nl+k,nl,nh);
			int im2 = MIRR_SL(k-1,nl);        int ip2 = MIRR_SL(k+1,nl);
			int im3 = MIRR_SH(nl+k-2,nl,nh);  int ip3 = MIRR_SH(nl+k+1,nl,nh);
			var_curr_idx[0] = im3;
			var_curr_idx[1] = im2;
			var_curr_idx[2] = im1;
			var_curr_idx[3] = i0;
			var_curr_idx[4] = ip1;
			var_curr_idx[5] = ip2;
			var_curr_idx[6] = ip3;

			/*
			p_in[2*k*stride] = 
				sl0 * t[k] + 
				sl2 * ( t[MIRR_SL(k-1,nl)] + t[MIRR_SL(k+1,nl)] ) +
				sh1 * ( t[MIRR_SH(nl+k-1,nl,nh)] + t[MIRR_SH(nl+k,nl,nh)] ) +
				sh3 * ( t[MIRR_SH(nl+k-2,nl,nh)] + t[MIRR_SH(nl+k+1,nl,nh)] );
				*/

			Print_Load_Line(fp,im3,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
			Print_Load_Line(fp,ip3,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
			fprintf(fp,"\tacc1 = _mm256_mul_ps(_mm_sh3,_mm256_add_ps(v%d,v%d));\n",Find_Index(var_prev_idx,num_vars,im3),Find_Index(var_prev_idx,num_vars,ip3));
			
			Print_Load_Line(fp,im2,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
                        Print_Load_Line(fp,ip2,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
			if (avx2)
				fprintf(fp,"\tacc1 = _mm256_fmadd_ps(_mm_sl2,_mm256_add_ps(v%d,v%d),acc1);\n",Find_Index(var_prev_idx,num_vars,im2),Find_Index(var_prev_idx,num_vars,ip2));
			else
				fprintf(fp,"\tacc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl2,_mm256_add_ps(v%d,v%d)));\n",Find_Index(var_prev_idx,num_vars,im2),Find_Index(var_prev_idx,num_vars,ip2));

			Print_Load_Line(fp,im1,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
                        Print_Load_Line(fp,ip1,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
			if (avx2)
				fprintf(fp,"\tacc1 = _mm256_fmadd_ps(_mm_sh1,_mm256_add_ps(v%d,v%d),acc1);\n",Find_Index(var_prev_idx,num_vars,im1),Find_Index(var_prev_idx,num_vars,ip1));
			else
				fprintf(fp,"\tacc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh1,_mm256_add_ps(v%d,v%d)));\n",Find_Index(var_prev_idx,num_vars,im1),Find_Index(var_prev_idx,num_vars,ip1));

			Print_Load_Line(fp,i0,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
			if (avx2)
				fprintf(fp,"\tacc1 = _mm256_fmadd_ps(_mm_sl0,v%d,acc1);\n",Find_Index(var_prev_idx,num_vars,i0));
			else
				fprintf(fp,"\tacc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl0,v%d));\n",Find_Index(var_prev_idx,num_vars,i0));

			if (!tmp_array) Print_Load_Line(fp,2*k,tmp_array,var_prev_idx,num_vars,var_curr_idx,7);
                        fprintf(fp,"\tdata[%d*stride] = acc1;\n",2*k);
		}
		{
                        fprintf(fp,"\n\t// odd samples :: k=%d\n",2*k+1);
			int i0 = nl+k;
			int im1 = MIRR_SL(k,nl);          int ip1 = MIRR_SL(k+1,nl);
			int im2 = MIRR_SH(nl+k-1,nl,nh);  int ip2 = MIRR_SH(nl+k+1,nl,nh);
			int im3 = MIRR_SL(k-1,nl);        int ip3 = MIRR_SL(k+2,nl);
			int im4 = MIRR_SH(nl+k-2,nl,nh);  int ip4 = MIRR_SH(nl+k+2,nl,nh);
			var_curr_idx[0] = im4;
			var_curr_idx[1] = im3;
			var_curr_idx[2] = im2;
			var_curr_idx[3] = im1;
			var_curr_idx[4] = i0;
			var_curr_idx[5] = ip1;
			var_curr_idx[6] = ip2;
			var_curr_idx[7] = ip3;
			var_curr_idx[8] = ip4;
			/*
			p_in[(2*k+1)*stride] = 
				sl1 * ( t[MIRR_SL(k,nl)] + t[MIRR_SL(k+1,nl)] ) +
				sl3 * ( t[MIRR_SL(k-1,nl)] + t[MIRR_SL(k+2,nl)] ) +
				sh0 * t[nl+k] +
				sh2 * ( t[MIRR_SH(nl+k-1,nl,nh)] + t[MIRR_SH(nl+k+1,nl,nh)] ) +
				sh4 * ( t[MIRR_SH(nl+k-2,nl,nh)] + t[MIRR_SH(nl+k+2,nl,nh)] );
			*/
			Print_Load_Line(fp,im4,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
			Print_Load_Line(fp,ip4,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
			fprintf(fp,"\tacc1 = _mm256_mul_ps(_mm_sh4,_mm256_add_ps(v%d,v%d));\n",Find_Index(var_prev_idx,num_vars,im4),Find_Index(var_prev_idx,num_vars,ip4));

			Print_Load_Line(fp,im3,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
                        Print_Load_Line(fp,ip3,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
			if (avx2)
				fprintf(fp,"\tacc1 = _mm256_fmadd_ps(_mm_sl3,_mm256_add_ps(v%d,v%d),acc1);\n",Find_Index(var_prev_idx,num_vars,im3),Find_Index(var_prev_idx,num_vars,ip3));
			else
				fprintf(fp,"\tacc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl3,_mm256_add_ps(v%d,v%d)));\n",Find_Index(var_prev_idx,num_vars,im3),Find_Index(var_prev_idx,num_vars,ip3));

			Print_Load_Line(fp,im2,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
                        Print_Load_Line(fp,ip2,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
			if (avx2)
				fprintf(fp,"\tacc1 = _mm256_fmadd_ps(_mm_sh2,_mm256_add_ps(v%d,v%d),acc1);\n",Find_Index(var_prev_idx,num_vars,im2),Find_Index(var_prev_idx,num_vars,ip2));
			else
				fprintf(fp,"\tacc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh2,_mm256_add_ps(v%d,v%d)));\n",Find_Index(var_prev_idx,num_vars,im2),Find_Index(var_prev_idx,num_vars,ip2));

			Print_Load_Line(fp,im1,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
                        Print_Load_Line(fp,ip1,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
			if (avx2)
				fprintf(fp,"\tacc1 = _mm256_fmadd_ps(_mm_sl1,_mm256_add_ps(v%d,v%d),acc1);\n",Find_Index(var_prev_idx,num_vars,im1),Find_Index(var_prev_idx,num_vars,ip1));
			else
				fprintf(fp,"\tacc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sl1,_mm256_add_ps(v%d,v%d)));\n",Find_Index(var_prev_idx,num_vars,im1),Find_Index(var_prev_idx,num_vars,ip1));

			Print_Load_Line(fp,i0,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
			if (avx2)
				fprintf(fp,"\tacc1 = _mm256_fmadd_ps(_mm_sh0,v%d,acc1);\n",Find_Index(var_prev_idx,num_vars,i0));
			else
				fprintf(fp,"\tacc1 = _mm256_add_ps(acc1,_mm256_mul_ps(_mm_sh0,v%d));\n",Find_Index(var_prev_idx,num_vars,i0));

			if (!tmp_array) Print_Load_Line(fp,2*k+1,tmp_array,var_prev_idx,num_vars,var_curr_idx,9);
                        fprintf(fp,"\tdata[%d*stride] = acc1;\n",2*k+1);
		}
	}
        fprintf(fp,"}\n\n");
}

void Gen_Us79(const char* path, int min_n, int max_n, int num_vars)
{
	FILE* fp = fopen(path, "w");

	fprintf(fp,"/*!\n");
	fprintf(fp," * Don't edit this code, it was automatically generated.\n");
	fprintf(fp," * Base functions for wavelet transforms of length %d to %d.\n",1<<min_n,1<<max_n);
	fprintf(fp," */\n");

	fprintf(fp,"/*\n");
	fprintf(fp," * Define coefficients for Antonini 7-9 tap filter.\n");
	fprintf(fp," */\n");
	fprintf(fp,"#define sl0  7.884856164056601e-001f\n");
	fprintf(fp,"#define sl1  4.180922732222101e-001f\n");
	fprintf(fp,"#define sl2 -4.068941760955800e-002f\n");
	fprintf(fp,"#define sl3 -6.453888262893799e-002f\n");

	fprintf(fp,"#define sh0  8.526986790094000e-001f\n");
	fprintf(fp,"#define sh1 -3.774028556126500e-001f\n");
	fprintf(fp,"#define sh2 -1.106244044184200e-001f\n");
	fprintf(fp,"#define sh3  2.384946501938001e-002f\n");
	fprintf(fp,"#define sh4  3.782845550699501e-002f\n\n");

	fprintf(fp,"#define _mm_sl0 _mm256_set1_ps(sl0)\n");
	fprintf(fp,"#define _mm_sl1 _mm256_set1_ps(sl1)\n");
	fprintf(fp,"#define _mm_sl2 _mm256_set1_ps(sl2)\n");
	fprintf(fp,"#define _mm_sl3 _mm256_set1_ps(sl3)\n");

	fprintf(fp,"#define _mm_sh0 _mm256_set1_ps(sh0)\n\n");
	fprintf(fp,"#define _mm_sh1 _mm256_set1_ps(sh1)\n");
	fprintf(fp,"#define _mm_sh2 _mm256_set1_ps(sh2)\n");
	fprintf(fp,"#define _mm_sh3 _mm256_set1_ps(sh3)\n");
	fprintf(fp,"#define _mm_sh4 _mm256_set1_ps(sh4)\n\n");

	fprintf(fp,"#ifdef __AVX2__\n\n");

	for (int i = min_n;  i <= max_n;  ++i) Gen_Us79_Core(fp,1<<i,num_vars,true);

	fprintf(fp,"#else\n\n");

	for (int i = min_n;  i <= max_n;  ++i) Gen_Us79_Core(fp,1<<i,num_vars,false);

	fprintf(fp,"#endif\n");

	fclose(fp);
	
	printf("Wrote Us79 base code to file %s.\n",path);
}
