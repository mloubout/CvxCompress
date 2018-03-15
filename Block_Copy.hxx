#ifndef CVX_CVXCOMPRESS_BLOCK_COPY_HXX
#define CVX_CVXCOMPRESS_BLOCK_COPY_HXX

#include <xmmintrin.h>

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
	);
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
	);

#endif
