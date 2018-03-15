#ifndef CVXCOMPRESS_READ_RAW_VOLUME_HXX
#define CVXCOMPRESS_READ_RAW_VOLUME_HXX

void 
Read_Raw_Volume(
		const char* filename, 
		int& nx, 
		int& ny, 
		int& nz, 
		float*& vol
	       );

#endif

