#ifndef CVX_SEISMOD_ASCII2EBCDIC
#define CVX_SEISMOD_ASCII2EBCDIC
#ifdef __cplusplus
extern "C"
{
#endif

void Convert_ASCII_To_EBCDIC(
	char* ASCII_str,
	char* EBCDIC_str,
	int num_chars
	);

#ifdef __cplusplus
}
#endif
#endif
