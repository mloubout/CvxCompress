#include "Block_Copy.hxx"

/*!
 * Efficiently copy a block of dimension bx*by*bz from volume with dimension nx*ny*nz starting at location x0,y0,z0.
 *
 * Arguments:
 *
 * data  Pointer to source volume
 * x0    .
 * y0    Start copy at this location
 * z0    .
 * nx    .
 * ny    Source volume dimensions
 * nz    .
 * work  Pointer to destination block
 * bx    . Must be multiple of 4.
 * by    Destination block dimensions.
 * bz    .
 *
 */
void Copy_To_Block(
	float* data,
	int x0,
	int y0,
	int z0,
	int nx,
	int ny,
	int nz,
	__m128* work,
	int bx,
	int by,
	int bz
	)
{
	int _mm_bx = bx >> 2;
	int nclipx = (x0+bx) - nx;
	int clipped_mm_bx = nclipx > 0 ? (bx - nclipx) >> 2 : _mm_bx;
	int clipped_bx = nclipx > 0 ? (bx - nclipx) : bx;
	int y_stop = y0+by < ny ? y0+by : ny;
	int iz;
	for (iz = z0;  iz < z0+bz && iz < nz;  ++iz)
	{
		int iy = y0;
		__m128* src = (__m128*)(data + (long)(iz*ny + iy)*nx + x0);
		__m128* dst = work + (long)((iz-z0)*by + (iy-y0))*_mm_bx;
		if (nclipx <= 0)
		{
			if (_mm_bx == 1)
			{
				for (;  iy < y_stop;  ++iy)
				{
					//printf("....copying iz=%d,iy=%d\n",iz,iy);
					*(dst++) = _mm_loadu_ps((float*)(src));
					src = (__m128*)(((float*)src) + nx);
				}
			}
			else if (_mm_bx == 2)
			{
				for (;  iy < y_stop;  ++iy)
				{
					//printf("....copying iz=%d,iy=%d\n",iz,iy);
					*(dst++) = _mm_loadu_ps((float*)(src));
					*(dst++) = _mm_loadu_ps((float*)(src+1));
					src = (__m128*)(((float*)src) + nx);
				}
			}
			else if (_mm_bx == 4)
			{
				for (;  iy < y_stop;  ++iy)
				{
					//printf("....copying iz=%d,iy=%d\n",iz,iy);
					*(dst++) = _mm_loadu_ps((float*)(src));
					*(dst++) = _mm_loadu_ps((float*)(src+1));
					*(dst++) = _mm_loadu_ps((float*)(src+2));
					*(dst++) = _mm_loadu_ps((float*)(src+3));
					src = (__m128*)(((float*)src) + nx);
				}
			}
			else
			{
				for (;  iy < y_stop;  ++iy)
				{
					//printf("....copying iz=%d,iy=%d\n",iz,iy);
					int ix;
					for (ix = 0;  ix < clipped_mm_bx;  ++ix) *(dst++) = _mm_loadu_ps((float*)(src+ix));
					src = (__m128*)(((float*)src) + nx);
				}
			}
		}
		else
		{
			for (;  iy < y_stop;  ++iy)
			{
				//printf("....copying iz=%d,iy=%d\n",iz,iy);
				int ix;
				for (ix = 0;  ix < clipped_mm_bx;  ++ix) dst[ix] = _mm_loadu_ps((float*)(src+ix));
				for (ix=ix*4;  ix < clipped_bx;  ++ix) ((float*)dst)[ix] = ((float*)src)[ix];
				for (;  ix < bx;  ++ix) ((float*)dst)[ix] = 0.0f;
				src = (__m128*)(((float*)src) + nx);
				dst += _mm_bx;
			}
		}
		for (;  iy < y0+by;  ++iy)
		{
			//printf("....zeroing iz=%d,iy=%d\n",iz,iy);
			for (int ix = 0;  ix < _mm_bx;  ++ix) dst[ix] = _mm_set1_ps(0.0f);
			dst += _mm_bx;
		}
	}
	for (;  iz < z0+bz;  ++iz)
	{
		//printf("...zeroing iz=%d\n",iz);
		__m128* dst = work + by*_mm_bx*(long)(iz-z0);
		for (int ixy = 0;  ixy < by*_mm_bx;  ++ixy) dst[ixy] = _mm_set1_ps(0.0f);
	}
}

/*!
 * Efficiently copy a block of dimensions bx*by*bz to destination volume starting at location x0,y0,z0.
 *
 * Arguments:
 *
 * work  Pointer to destination block
 * bx    . Must be multiple of 4.
 * by    Destination block dimensions.
 * bz    .
 * data  Pointer to source volume
 * x0    .
 * y0    Start copy at this location
 * z0    .
 * nx    .
 * ny    Source volume dimensions
 * nz    .
 *
 */
void Copy_From_Block(
	__m128* work,
	int bx,
	int by,
	int bz,
	float* data,
	int x0,
	int y0,
	int z0,
	int nx,
	int ny,
	int nz
	)
{
	int _mm_bx = bx >> 2;
	int nclipx = (x0+bx) - nx;
	int clipped_mm_bx = nclipx > 0 ? (bx - nclipx) >> 2 : bx >> 2;
	int clipped_bx = nclipx > 0 ? (bx - nclipx) : bx;
	int y_stop = y0+by < ny ? y0+by : ny;
	for (int iz = z0;  iz < z0+bz && iz < nz;  ++iz)
	{
		int iy = y0;
		__m128i* dst = (__m128i*)(data + (long)(iz*ny + iy)*nx + x0);
		__m128i* src = (__m128i*)work + (long)((iz-z0)*by + (iy-y0))*_mm_bx;
		if (nclipx <= 0)
		{
			if (_mm_bx == 1)
			{
				for (;  iy < y_stop;  ++iy)
				{
					_mm_storeu_si128(dst, _mm_load_si128(src++));
					dst = (__m128i*)(((float*)dst) + nx);
				}
			}
			else if (_mm_bx == 2)
			{
				for (;  iy < y_stop;  ++iy)
				{
					_mm_storeu_si128(dst, _mm_load_si128(src++));
					_mm_storeu_si128(dst+1, _mm_load_si128(src++));
					dst = (__m128i*)(((float*)dst) + nx);
				}
			}
			else if (_mm_bx == 4)
			{
				for (;  iy < y_stop;  ++iy)
				{
					_mm_storeu_si128(dst, _mm_load_si128(src++));
					_mm_storeu_si128(dst+1, _mm_load_si128(src++));
					_mm_storeu_si128(dst+2, _mm_load_si128(src++));
					_mm_storeu_si128(dst+3, _mm_load_si128(src++));
					dst = (__m128i*)(((float*)dst) + nx);
				}
			}
			else
			{
				for (;  iy < y_stop;  ++iy)
				{
					int ix;
					for (ix = 0;  ix < clipped_mm_bx;  ++ix) _mm_storeu_si128(dst+ix, _mm_load_si128(src++));
					dst = (__m128i*)(((float*)dst) + nx);
				}
			}
		}
		else
		{
			for (;  iy < y_stop;  ++iy)
			{
				int ix;
				for (ix = 0;  ix < clipped_mm_bx;  ++ix) _mm_storeu_si128(dst+ix, _mm_load_si128(src+ix));
				for (ix=ix*4;  ix < clipped_bx;  ++ix) ((float*)dst)[ix] = ((float*)src)[ix];
				dst = (__m128i*)(((float*)dst) + nx);
				src += _mm_bx;
			}
		}
	}
}
