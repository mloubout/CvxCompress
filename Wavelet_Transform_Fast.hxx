#ifndef CVX_WAVELET_TRANSFORM_FAST_HXX
#define CVX_WAVELET_TRANSFORM_FAST_HXX

#include <immintrin.h>

/*!
 * Perform forward wavelet transform.
 * Arguments:
 * work - pointer to the block you want to transform. must be aligned on 32 byte boundary.
 * tmp  - temporary buffer used internally. must be at least 8*MAX(bx,by,bz) floats large and be aligned on 32 byte boundary.
 * bx   - x block size (number of floats)
 * by   - y block size (number of floats)
 * bz   - z block size (number of floats)
 *
 */
void Wavelet_Transform_Fast_Forward(
	__m256* work,
	__m256* tmp,
	int bx,
	int by,
	int bz
	);

/*!
 * Perform inverse wavelet transform.
 * Arguments:
 * work - pointer to the block you want to transform. must be aligned on 32 byte boundary.
 * tmp  - temporary buffer used internally. must be at least 8*MAX(bx,by,bz) floats large and be aligned on 32 byte boundary.
 * bx   - x block size (number of floats)
 * by   - y block size (number of floats)
 * bz   - z block size (number of floats)
 *
 */
void Wavelet_Transform_Fast_Inverse(
	__m256* work,
	__m256* tmp,
	int bx,
	int by,
	int bz
	);

#endif
