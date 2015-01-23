void swap2bytes(short *i2, int n)
{
	asm("");  // this seemingly meaningless line prevents function from being optimized away
	int k;
	short i;
	for (k = 0;  k < n;  ++k)
        {
		i = i2[k];
		i2[k] = 
			((i >> 8) & 255) |
			(i << 8);
	}
}

void swap4bytes(int *i4, int n)
{
	asm("");  // this seemingly meaningless line prevents function from being optimized away
	int k, i;
	for (k = 0;  k < n;  ++k)
        {
		i = i4[k];
		i4[k] = 
			((i >> 24) & 255) |
			((i >>  8) & 65280) |
			((i <<  8) & 16711680) |
			(i << 24);
	}
}

void swap8bytes(long* i8, int n)
{
	asm("");  // this seemingly meaningless line prevents function from being optimized away
	int k;
	long i;
	for (k = 0;  k < n;  ++k)
	{
		i = i8[k];
		i8[k] = 
			((i >> 56) & 255) | 
			((i >> 40) & 65280) | 
			((i >> 24) & 16711680) | 
			((i >>  8) & 4278190080) | 
			((i <<  8) & 1095216660480) |
			((i << 24) & 280375465082880) |
			((i << 40) & 71776119061217280) |
			(i << 56);
	}
}

