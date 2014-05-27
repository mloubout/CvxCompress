#include <stdio.h>

void swap_endian(float* v)
{
	int* iv = (int*)v;
	*iv = 
		(((*iv) >> 24) & 255) | 
		(((*iv) >>  8) & 65280) | 
		(((*iv) & 65280) <<  8) |
		(((*iv) & 255) << 24);
}

int main(int argc, char** argv)
{
	FILE* fp = fopen(argv[1], "r");
	if (fp != 0L)
	{
		long idx = 0;
		float min = 1e38f;
		float max = -1e38f;
		float f[1024];
		for (
				size_t nread = fread(f, sizeof(float), 1024, fp);  
				nread > 0;  
				nread = fread(f, sizeof(float), 1024, fp))
		{
		 	for (int i = 0;  i < nread;  ++i)
			{
				swap_endian(f+i);
				float v = f[i];
				if (v < min) min = v;
				if (v > max) max = v;
				//if (v == 0.0f)
				//{
				//	printf("Perfect ZERO at idx %ld\n",idx+i);
				//}
			}
			idx += nread;
		}
		fclose(fp);
		printf("min = %f, max = %f\n",min,max);
	}
	return 0;
}
