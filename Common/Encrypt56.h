#ifndef CVX_COMMON_ENCRYPT56_H
#define CVX_COMMON_ENCRYPT56_H

#include <unistd.h>
#include <cstdlib>
#include <crypt.h>

static long __key56 = 0x78487AE9123981B6;

static struct crypt_data Setkey56(long key56)
{
	struct crypt_data data;
        data.initialized = 0;
        char key[64];
        for (int i = 0;  i < 64;  ++i)
        {
		key[i] = ((key56 & 1) == 1) ? 1 : 0;
		key56 = key56 >> 1;
        }
	setkey_r(key,&data);
	return data;
}

static long Encrypt56(float v1, float v2)
{
	struct crypt_data data = Setkey56(__key56);
	char block[64];
	long lval = 0;
	char* l = (char*)&lval;
	char* p = (char*)&v1;
	l[0] = p[0];
	l[1] = p[1];
	l[2] = p[2];
	l[3] = p[3];
	p = (char*)&v2;
	l[4] = p[0];
	l[5] = p[1];
	l[6] = p[2];
	l[7] = p[3];
	for (int i = 0;  i < 64;  ++i)
	{
		block[i] = ((lval & 1) == 1) ? 1 : 0;
		lval = lval >> 1;
	}
	encrypt_r(block,0,&data);
	lval = 0;
	for (int i = 0;  i < 64;  ++i)
	{
		lval |= block[63-i];
		if (i < 63) lval = lval << 1;
	}
	return lval;
}

static long Encrypt56(int v1, int v2)
{
	return Encrypt56(*((float*)&v1),*((float*)&v2));
}

static void Decrypt56(long encrypted, float& v1, float& v2)
{
	long lval = encrypted;
	struct crypt_data data = Setkey56(__key56);
	v1 = 0.0f;
	v2 = 0.0f;
	char block[64];
	for (int i = 0;  i < 64;  ++i)
        {
		block[i] = ((lval & 1) == 1) ? 1 : 0;
		lval = lval >> 1;
	}
	encrypt_r(block,1,&data);
	lval = 0;
	for (int i = 0;  i < 64;  ++i)
        {
                lval |= block[63-i];
                if (i < 63) lval = lval << 1;
        }
	char* l = (char*)&lval;
        char* p = (char*)&v1;
	p[0] = l[0];
	p[1] = l[1];
	p[2] = l[2];
	p[3] = l[3];
	p = (char*)&v2;
	p[0] = l[4];
        p[1] = l[5];
        p[2] = l[6];
        p[3] = l[7];
}

static void Decrypt56(long encrypted, int& v1, int& v2)
{
	float f1,f2;
	Decrypt56(encrypted,f1,f2);
	v1 = *((int*)&f1);
	v2 = *((int*)&f2);
}

#endif

