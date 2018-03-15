//
// This code implements a fast 3D wavelet transform based on Antonini's 7-9 tap filter.
// This code performs the same functions as the subroutine ChvDs79 in Ergas et.al. ChvCompress code.
// In order to get maximum performance out of the code, only block sizes that are powers of two
// between 8 and 256 are supported.
// This means the smallest supported block is 8x8x8 and the largest is 256x256x256.
//

// Base functions for wavelet transform. These were auto-generated.
#include "Ds79_Base.cpp"
#include "Us79_Base.cpp"

void Wavelet_Transform_Fast_Forward(
	__m256* work,
	__m256* tmp,
	int bx,
	int by,
	int bz
	)
{
	int _mm_bx = bx >> 2;
	int _mm256_bx = bx >> 3;
	int _mm256_stride_y = _mm256_bx;
	for (int iz = 0;  iz < bz;  ++iz)
        {
		// x
		for (int iy = 0;  iy < by;  iy+=8)
        	{
			__m128* data = ((__m128*)work) + (iz*by + iy) * _mm_bx;
			for (int ix = 0;  ix < _mm_bx;  ++ix)
			{
				__m128 v0 = data[ix];
				__m128 v1 = data[ix+_mm_bx];
				__m128 v2 = data[ix+2*_mm_bx];
				__m128 v3 = data[ix+3*_mm_bx];
				_MM_TRANSPOSE4_PS(v0,v1,v2,v3);
				__m128 v4 = data[ix+4*_mm_bx];
				__m128 v5 = data[ix+5*_mm_bx];
				__m128 v6 = data[ix+6*_mm_bx];
				__m128 v7 = data[ix+7*_mm_bx];
				_MM_TRANSPOSE4_PS(v4,v5,v6,v7);
				tmp[ix*4] = _mm256_insertf128_ps(_mm256_castps128_ps256(v0),v4,1);
				tmp[ix*4+1] = _mm256_insertf128_ps(_mm256_castps128_ps256(v1),v5,1);
				tmp[ix*4+2] = _mm256_insertf128_ps(_mm256_castps128_ps256(v2),v6,1);
				tmp[ix*4+3] = _mm256_insertf128_ps(_mm256_castps128_ps256(v3),v7,1);
			}
		
			if (bx >= 256) _Ds79_AVX_256(tmp, 1);
			if (bx >= 128) _Ds79_AVX_128(tmp, 1);
			if (bx >= 64) _Ds79_AVX_64(tmp, 1);
			if (bx >= 32) _Ds79_AVX_32(tmp, 1);
			if (bx >= 16) _Ds79_AVX_16(tmp, 1);
			_Ds79_AVX_8(tmp,1);
			_Ds79_AVX_4(tmp,1);
			_Ds79_AVX_2(tmp,1);

			for (int ix = 0;  ix < _mm_bx;  ++ix)
			{
				__m256 lv0 = tmp[ix*4];
				__m256 lv1 = tmp[ix*4+1];
				__m256 lv2 = tmp[ix*4+2];
				__m256 lv3 = tmp[ix*4+3];
				__m128 v0 = _mm256_extractf128_ps(lv0,0);
				__m128 v1 = _mm256_extractf128_ps(lv1,0);
				__m128 v2 = _mm256_extractf128_ps(lv2,0);
				__m128 v3 = _mm256_extractf128_ps(lv3,0);
				_MM_TRANSPOSE4_PS(v0,v1,v2,v3);
				__m128 v4 = _mm256_extractf128_ps(lv0,1);
				__m128 v5 = _mm256_extractf128_ps(lv1,1);
				__m128 v6 = _mm256_extractf128_ps(lv2,1);
				__m128 v7 = _mm256_extractf128_ps(lv3,1);
				_MM_TRANSPOSE4_PS(v4,v5,v6,v7);
				data[ix] = v0;
				data[ix+_mm_bx] = v1;
				data[ix+2*_mm_bx] = v2;
				data[ix+3*_mm_bx] = v3;
				data[ix+4*_mm_bx] = v4;
				data[ix+5*_mm_bx] = v5;
				data[ix+6*_mm_bx] = v6;
				data[ix+7*_mm_bx] = v7;
			}
		}
	
		// y
		__m256* data = work + iz*by*_mm256_bx;
		for (int ix = 0;  ix < _mm256_bx;  ++ix)
		{
			if (by >= 256) _Ds79_AVX_256(data+ix, _mm256_stride_y);
			if (by >= 128) _Ds79_AVX_128(data+ix, _mm256_stride_y);
			if (by >= 64) _Ds79_AVX_64(data+ix, _mm256_stride_y);
			if (by >= 32) _Ds79_AVX_32(data+ix, _mm256_stride_y);
			if (by >= 16) _Ds79_AVX_16(data+ix, _mm256_stride_y);
			_Ds79_AVX_8(data+ix, _mm256_stride_y);
			_Ds79_AVX_4(data+ix, _mm256_stride_y);
			_Ds79_AVX_2(data+ix, _mm256_stride_y);
		}
	}
	
	// z
	if (bz > 1)
	{
		int _mm256_stride_z = by * _mm256_bx;
		for (int iy = 0;  iy < by;  ++iy)
		{
			__m256* data = work + iy*_mm256_bx;
			if ((bx*by) >= 1024 && ((bx*by)&1023) == 0)
			{
				// z stride is a multiple of 4096 bytes.
				// this is bad because all values along the Z axis maps to the same cache page.
				// this reduces effective L1 cache size from 32KB to 0.5KB and has negative effects
				// on L2 cache as well. We prevent this by copying the inputs to a temporary buffer.
				for (int ix = 0;  ix < _mm256_bx;  ++ix)
				{
					for (int iz = 0;  iz < bz;  ++iz) tmp[iz] = data[iz*_mm256_stride_z+ix];

					if (bz >= 256) _Ds79_AVX_256(tmp, 1);
					if (bz >= 128) _Ds79_AVX_128(tmp, 1);
					if (bz >= 64) _Ds79_AVX_64(tmp, 1);
					if (bz >= 32) _Ds79_AVX_32(tmp, 1);
					if (bz >= 16) _Ds79_AVX_16(tmp, 1);
					// smallest block size is 8, so no need to test for anything less than that.
					_Ds79_AVX_8(tmp, 1);
					_Ds79_AVX_4(tmp, 1);
					_Ds79_AVX_2(tmp, 1);

					for (int iz = 0;  iz < bz;  ++iz) data[iz*_mm256_stride_z+ix] = tmp[iz];
				}
			}
			else
			{
				for (int ix = 0;  ix < _mm256_bx;  ++ix)
				{
					if (bz >= 256) _Ds79_AVX_256(data+ix, _mm256_stride_z);
					if (bz >= 128) _Ds79_AVX_128(data+ix, _mm256_stride_z);
					if (bz >= 64) _Ds79_AVX_64(data+ix, _mm256_stride_z);
					if (bz >= 32) _Ds79_AVX_32(data+ix, _mm256_stride_z);
					if (bz >= 16) _Ds79_AVX_16(data+ix, _mm256_stride_z);
					// smallest block size is 8, so no need to test for anything less than that.
					_Ds79_AVX_8(data+ix, _mm256_stride_z);
					_Ds79_AVX_4(data+ix, _mm256_stride_z);
					_Ds79_AVX_2(data+ix, _mm256_stride_z);
				}
			}
		}
	}
}

void Wavelet_Transform_Fast_Inverse(
	__m256* work,
	__m256* tmp,
	int bx,
	int by,
	int bz
	)
{
	int _mm_bx = bx >> 2;
	int _mm256_bx = bx >> 3;
	int _mm256_stride_y = _mm256_bx;
	for (int iz = 0;  iz < bz;  ++iz)
        {
		// x
		for (int iy = 0;  iy < by;  iy+=8)
        	{
			__m128* data = ((__m128*)work) + (iz*by + iy) * _mm_bx;
			for (int ix = 0;  ix < _mm_bx;  ++ix)
			{
				__m128 v0 = data[ix];
				__m128 v1 = data[ix+_mm_bx];
				__m128 v2 = data[ix+2*_mm_bx];
				__m128 v3 = data[ix+3*_mm_bx];
				_MM_TRANSPOSE4_PS(v0,v1,v2,v3);
				__m128 v4 = data[ix+4*_mm_bx];
				__m128 v5 = data[ix+5*_mm_bx];
				__m128 v6 = data[ix+6*_mm_bx];
				__m128 v7 = data[ix+7*_mm_bx];
				_MM_TRANSPOSE4_PS(v4,v5,v6,v7);
				tmp[ix*4] = _mm256_insertf128_ps(_mm256_castps128_ps256(v0),v4,1);
				tmp[ix*4+1] = _mm256_insertf128_ps(_mm256_castps128_ps256(v1),v5,1);
				tmp[ix*4+2] = _mm256_insertf128_ps(_mm256_castps128_ps256(v2),v6,1);
				tmp[ix*4+3] = _mm256_insertf128_ps(_mm256_castps128_ps256(v3),v7,1);
			}

			_Us79_AVX_2(tmp,1);
			_Us79_AVX_4(tmp,1);
			_Us79_AVX_8(tmp,1);
			if (bx >= 16) _Us79_AVX_16(tmp,1);
			if (bx >= 32) _Us79_AVX_32(tmp,1);
			if (bx >= 64) _Us79_AVX_64(tmp,1);
			if (bx >= 128) _Us79_AVX_128(tmp,1);
			if (bx >= 256) _Us79_AVX_256(tmp,1);

			for (int ix = 0;  ix < _mm_bx;  ++ix)
			{
				__m256 lv0 = tmp[ix*4];
				__m256 lv1 = tmp[ix*4+1];
				__m256 lv2 = tmp[ix*4+2];
				__m256 lv3 = tmp[ix*4+3];
				__m128 v0 = _mm256_extractf128_ps(lv0,0);
				__m128 v1 = _mm256_extractf128_ps(lv1,0);
				__m128 v2 = _mm256_extractf128_ps(lv2,0);
				__m128 v3 = _mm256_extractf128_ps(lv3,0);
				_MM_TRANSPOSE4_PS(v0,v1,v2,v3);
				__m128 v4 = _mm256_extractf128_ps(lv0,1);
				__m128 v5 = _mm256_extractf128_ps(lv1,1);
				__m128 v6 = _mm256_extractf128_ps(lv2,1);
				__m128 v7 = _mm256_extractf128_ps(lv3,1);
				_MM_TRANSPOSE4_PS(v4,v5,v6,v7);
				data[ix] = v0;
				data[ix+_mm_bx] = v1;
				data[ix+2*_mm_bx] = v2;
				data[ix+3*_mm_bx] = v3;
				data[ix+4*_mm_bx] = v4;
				data[ix+5*_mm_bx] = v5;
				data[ix+6*_mm_bx] = v6;
				data[ix+7*_mm_bx] = v7;
			}
		}

		// y
		__m256* data = work + iz*by*_mm256_bx;
		for (int ix = 0;  ix < _mm256_bx;  ++ix)
		{
			_Us79_AVX_2(data+ix, _mm256_stride_y);
			_Us79_AVX_4(data+ix, _mm256_stride_y);
			_Us79_AVX_8(data+ix, _mm256_stride_y);
			if (by >= 16) _Us79_AVX_16(data+ix, _mm256_stride_y);
			if (by >= 32) _Us79_AVX_32(data+ix, _mm256_stride_y);
			if (by >= 64) _Us79_AVX_64(data+ix, _mm256_stride_y);
			if (by >= 128) _Us79_AVX_128(data+ix, _mm256_stride_y);
			if (by >= 256) _Us79_AVX_256(data+ix, _mm256_stride_y);
		}
	}
	
	// z
	if (bz > 1)
	{
		int _mm256_stride_z = by * _mm256_bx;
		for (int iy = 0;  iy < by;  ++iy)
		{
			__m256* data = work + iy*_mm256_bx;
			if ((bx*by) >= 1024 && ((bx*by)&1023) == 0)
			{
				// z stride is a multiple of 4096 bytes.
				// this is bad because all values along the Z axis maps to the same cache page.
				// this reduces effective L1 cache size from 32KB to 0.5KB and has negative effects
				// on L2 cache as well. We prevent this by copying the inputs to a temporary buffer.
				for (int ix = 0;  ix < _mm256_bx;  ++ix)
				{
					for (int iz = 0;  iz < bz;  ++iz) tmp[iz] = data[iz*_mm256_stride_z+ix];

					_Us79_AVX_2(tmp, 1);
					_Us79_AVX_4(tmp, 1);
					_Us79_AVX_8(tmp, 1);
					if (bz >= 16) _Us79_AVX_16(tmp, 1);
					if (bz >= 32) _Us79_AVX_32(tmp, 1);
					if (bz >= 64) _Us79_AVX_64(tmp, 1);
					if (bz >= 128) _Us79_AVX_128(tmp, 1);
					if (bz >= 256) _Us79_AVX_256(tmp, 1);

					for (int iz = 0;  iz < bz;  ++iz) data[iz*_mm256_stride_z+ix] = tmp[iz];
				}
			}
			else
			{
				for (int ix = 0;  ix < _mm256_bx;  ++ix)
				{
					_Us79_AVX_2(data+ix, _mm256_stride_z);
					_Us79_AVX_4(data+ix, _mm256_stride_z);
					_Us79_AVX_8(data+ix, _mm256_stride_z);
					if (bz >= 16) _Us79_AVX_16(data+ix, _mm256_stride_z);
					if (bz >= 32) _Us79_AVX_32(data+ix, _mm256_stride_z);
					if (bz >= 64) _Us79_AVX_64(data+ix, _mm256_stride_z);
					if (bz >= 128) _Us79_AVX_128(data+ix, _mm256_stride_z);
					if (bz >= 256) _Us79_AVX_256(data+ix, _mm256_stride_z);
				}
			}
		}
	}
}

