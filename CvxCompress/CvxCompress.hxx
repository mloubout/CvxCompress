#ifndef CVX_CVXCOMPRESS_HXX
#define CVX_CVXCOMPRESS_HXX

/*!
 * This class implements lossy compression for seismic wavefields.
 * It uses Antonini's 7-9 tap filter wavelet transform.
 * It was derived from the ChvCompress code by Ergas et.al.
 * Recoded with AVX & AVX2 intrinisics for maximum throughput.
 * Words of praise or heaps of scorn can be directed at thor.johnsen@chevron.com.
 */
class CvxCompress
{
public:
	CvxCompress();
	virtual ~CvxCompress();

	/*!< Compress a 3D wavefield using a given block size */
	float Compress(
			float scale,
			float* vol,
			int nx,
			int ny,
			int nz,
			int bx,
			int by,
			int bz,
			unsigned int* compressed,
			long& compressed_length 
		      );
	/*!< Decompress a 3D wavefield that was compressed with Compress(...) method */
	void Decompress(
			float*& vol,
			int& nx,
			int& ny,
			int& nz,
			unsigned int* compressed,
			long compressed_length 
			);

	static int Min_BX() {return  8;}  /*!< Get minimum X block size. Will always be a power of two.*/
	static int Max_BX() {return 256;}  /*!< Get maximum X block size. Will always be a power of two.*/
	static int Min_BY() {return  8;}  /*!< Get minimum Y block size. Will always be a power of two.*/
	static int Max_BY() {return 256;}  /*!< Get maximum Y block size. Will always be a power of two.*/
	static int Min_BZ() {return  8;}  /*!< Get minimum Z block size. Will always be a power of two.*/
	static int Max_BZ() {return 256;}  /*!< Get maximum Z block size. Will always be a power of two.*/

	bool Run_Module_Tests(bool verbose, bool exhaustive_throughput_tests);  /*!< Execute module tests.*/

private:
	bool _Valid_Block_Size(int bx, int by, int bz);
};

#endif

