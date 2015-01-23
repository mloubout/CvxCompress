#ifndef CVX_MI_UTILS_SWAPBYTES_H
#define CVX_MI_UTILS_SWAPBYTES_H

#ifdef __cplusplus
extern "C"
{
#endif

void swap2bytes(short *i2, int n);
void swap4bytes(int *i4, int n);
void swap8bytes(long* i8, int n);

#ifdef __cplusplus
}
#endif

#endif

