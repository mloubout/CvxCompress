#ifndef CVX_MI_UTILS_SWAPBYTES_H
#define CVX_MI_UTILS_SWAPBYTES_H

#ifdef __cplusplus
extern "C"
{
#endif

void swap2bytes(short *i2, long n);
void swap4bytes(int *i4, long n);
void swap8bytes(long* i8, long n);

#ifdef __cplusplus
}
#endif

#endif

