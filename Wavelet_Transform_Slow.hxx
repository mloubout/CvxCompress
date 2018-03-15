#ifndef CVX_TJHC_WAVELET_TRANSFORM_SLOW
#define CVX_TJHC_WAVELET_TRANSFORM_SLOW

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
	);

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
	);

void Gen_Ds79(const char* path, int min_n, int max_n, int num_vars);
void Gen_Us79(const char* path, int min_n, int max_n, int num_vars);

#endif
