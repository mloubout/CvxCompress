/** Performs byte swapping using Joe Stefani's swap4bytes.c **/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void swapbytes_ (float *temp, int *n)
{
 void swap4(int*, int);
 swap4((int*)temp, *n);
}


/***** swap4bytes:  dcba --> abcd *****/
void swap4(int *i4, int n)
{
  int k, i, a, b, c, d, bmask = 16711680, cmask = 65280, dmask = 255;
  for(k=0; k<n; k++)
  { i = i4[k];
    a =  i << 24;          
    b = (i << 8)  & bmask;
    c = (i >> 8) & cmask;  
    d = (i >> 24) & dmask;
    i4[k] = a | b | c | d ;
  }
}
