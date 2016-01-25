#include <assert.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include "Run_Length_Escape_Codes.hxx"

// un-comment if you want debug printouts during encoding
//#define DEBUG_ENCODE
//#define DEBUG_DECODE

// un-comment if you want individual byte mem read and writes
//#define BYTEIO

inline void EncodeRLE_Slow(int& rle, char* dst, int& bytepos)
{
	if (rle > 0)
	{
		if (rle == 1)
		{
#ifdef DEBUG_ENCODE
			printf("BYTE encode %d at index %d\n",(char)0,bytepos);
#endif
			dst[bytepos++] = (char)0;
		}
		/*
		else if (rle == 2)
		{
			dst[bytepos++] = (char)0;
			dst[bytepos++] = (char)0;
		}
		*/
		else if (rle < 256)
		{
#ifdef DEBUG_ENCODE
			printf("RLESC1 rle = %d at index %d\n",rle,bytepos);
#endif
			dst[bytepos++] = (char)RLESC1;
			dst[bytepos++] = (unsigned char)rle;
		}
		else
		{
#ifdef DEBUG_ENCODE
			printf("RLESC3 rle = %d at index %d\n",rle,bytepos);
#endif
			dst[bytepos++] = (char)RLESC3;
			dst[bytepos++] = rle & 0xFF;
			dst[bytepos++] = (rle >> 8) & 0xFF;
			dst[bytepos++] = (rle >> 16) & 0xFF;
		}
		rle = 0;
	}
}

/*
 * Run length encode a quantized block. 
 * The encoded block is usually smaller, but can be up to 25% larger after encoding.
 * Please ensure that work block (compressed ptr) is at least 5*sizeof(float)*num/4 bytes.
 */
void Run_Length_Encode_Slow(float scale, float* vals, int num, unsigned long* compressed, int& bytepos)
{
	int rle = 0;
	char* dst = (char*)compressed;
	for (int i = 0;  i < num;  ++i)
	{
		float fval = scale*vals[i];
		int ival = (int)fval;
		if (ival == 0)
		{
			++rle;
		}
		else
		{
			EncodeRLE_Slow(rle,dst,bytepos);
			if (ival > VLESC2 && ival < RLESC3)
			{
#ifdef DEBUG_ENCODE
				printf("BYTE encode %d at index %d\n",(char)ival,bytepos);
#endif
				dst[bytepos++] = (char)ival;
			}
			else if (ival >= -32768 && ival <= 32767)
			{
#ifdef DEBUG_ENCODE
				printf("VLESC2 encode %d at index %d\n",ival,bytepos);
#endif
				dst[bytepos++] = (char)VLESC2;
#ifndef BYTEIO
				*((short*)(dst+bytepos)) = ival & 0xFFFF;
				bytepos += 2;
#else
				dst[bytepos++] = ival & 0xFF;
				dst[bytepos++] = (ival >> 8) & 0xFF;
#endif
			}
			else if (ival >= -8388608 && ival <= 8388607)
			{
#ifdef DEBUG_ENCODE
				printf("VLESC3 encode %d at index %d\n",ival,bytepos);
#endif
				dst[bytepos++] = (char)VLESC3;
#ifndef BYTEIO
				*((int*)(dst+bytepos)) = ival & 0xFFFFFF;
				bytepos += 3;
#else
				dst[bytepos++] = ival & 0xFF;
				dst[bytepos++] = (ival >> 8) & 0xFF;
				dst[bytepos++] = (ival >> 16) & 0xFF;
#endif
			}
			else
			{
#ifdef DEBUG_ENCODE
                                printf("VLESC4 encode %e at index %d\n",fval,bytepos);
#endif
                                dst[bytepos++] = (char)VLESC4;
#ifndef BYTEIO
				*((float*)(dst+bytepos)) = fval;
				bytepos += 4;
#else
				unsigned int ifval = *((int*)&fval);
                                dst[bytepos++] = ifval & 0xFF;
                                dst[bytepos++] = (ifval >> 8) & 0xFF;
                                dst[bytepos++] = (ifval >> 16) & 0xFF;
                                dst[bytepos++] = (ifval >> 24) & 0xFF;
#endif
			}
		}
	}
	EncodeRLE_Slow(rle,dst,bytepos);
}

int Run_Length_Decode_Slow(float scale, float* vals, int num_expected_vals, unsigned long* compressed)
{
#ifdef DEBUG_DECODE
	printf("*** WTF! ***\n");
#endif
	int num = 0;
	char* p = (char*)compressed;
	float scalefac = 1.0f / scale;
	for (;  num < num_expected_vals;  ++p)
	{
		int ival = *p;
		if (ival == RLESC1)
		{
			int rle = ((unsigned char*)p)[1];
#ifdef DEBUG_DECODE
			printf("  RLESC1 rle=%d, num=%d, num_expected_vals=%d\n",rle,num,num_expected_vals);
			assert(num+rle <= num_expected_vals);
#endif
			for (int j = 0;  j < rle;  ++j) vals[num+j] = 0.0f;
			num += rle;
			p += 1;
		}
		else if (ival == RLESC3)
		{
			int rle = (((unsigned char*)p)[3] << 16) | (((unsigned char*)p)[2] << 8) | ((unsigned char*)p)[1];
#ifdef DEBUG_DECODE
			printf("  RLESC3 rle=%d, num=%d, num_expected_vals=%d\n",rle,num,num_expected_vals);
			assert(num+rle <= num_expected_vals);
#endif
			for (int j = 0;  j < rle;  ++j) vals[num+j] = 0.0f;
			num += rle;
			p += 3;
		}
		else if (ival == VLESC2)
		{
			int quant = (p[2] << 8) | ((unsigned char*)p)[1];
#ifdef DEBUG_DECODE
			printf("  VLESC2 quant=%d, num=%d\n",quant,num);
			assert(num < num_expected_vals);
#endif
			vals[num++] = (float)quant * scalefac;
			p += 2;
		}
		else if (ival == VLESC3)
		{
			int quant = (p[3] << 16) | (((unsigned char*)p)[2] << 8) | ((unsigned char*)p)[1];
#ifdef DEBUG_DECODE
			printf("  VLESC3 quant=%d, num=%d\n",quant,num);
			assert(num < num_expected_vals);
#endif
			vals[num++] = (float)quant * scalefac;
			p += 3;
		}
		else if (ival == VLESC4)
		{
			unsigned int fval = (p[4] << 24) | (((unsigned char*)p)[3] << 16) | (((unsigned char*)p)[2] << 8) | ((unsigned char*)p)[1];
#ifdef DEBUG_DECODE
			printf("  VLESC4 fval=%e, num=%d\n",*((float*)&fval),num);
			assert(num < num_expected_vals);
#endif
			vals[num++] = *((float*)&fval) * scalefac;
			p += 4;
		}
		else
		{
#ifdef DEBUG_DECODE
			printf("  BYTE=%d, num=%d\n",ival,num);
			assert(num < num_expected_vals);
#endif
			vals[num++] = (float)ival * scalefac;
		}
	}
	return num;
}

bool Run_Length_Encode_Compare(unsigned long* compressed, int bytepos, unsigned long* compressed2, int bytepos2)
{
	bool retval = false;
	printf("bytepos = %d, bytepos2 = %d\n",bytepos,bytepos2);
	if (bytepos == bytepos2)
	{
		retval = true;
		char* p = (char*)compressed;
		char* p2 = (char*)compressed2;
		for (int i = 0;  i < bytepos;  ++i)
		{
#ifdef DEBUG_ENCODE
			printf("p[%4d] = %4d, p2[%4d] = %4d\n",i,p[i],i,p2[i]);
#endif
			if (p[i] != p2[i])
			{
#ifndef DEBUG_ENCODE
				printf("p[%4d] = %4d, p2[%4d] = %4d\n",i,p[i],i,p2[i]);
#endif
				printf("Arrays differ at byte %d\n",i);
				retval = false;
			}
		}
	}
	if (retval) printf("Arrays are identical\n"); else printf("Arrays differ\n");
	return retval;
}

