#ifndef CVX_CVXCOMPRESS_BLOCK_COPY_HXX
#define CVX_CVXCOMPRESS_BLOCK_COPY_HXX

#ifndef MY_AVX_DEFINED
	#define SIMDE_ENABLE_NATIVE_ALIASES
	#include <simde/x86/avx512.h>  // SSE intrinsics
#endif

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
