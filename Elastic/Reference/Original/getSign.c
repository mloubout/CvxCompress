/** Gets sign of a float **/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

float getSign (float x)
{
  if (x > 0) 
     return 1;
  else if (x < 0) 
      return -1;
  else
    return 0;
}

