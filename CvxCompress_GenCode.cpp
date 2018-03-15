#include <stdio.h>
#include "CvxCompress.hxx"
#include "Wavelet_Transform_Slow.hxx"

static int Find_Pow2(int val)
{
	int cnt = -1;
	while (val > 0)
	{
		val = val >> 1;
		++cnt;
	}
	return cnt;
}

int main(int argc, char* argv[])
{
	const char* path_Ds79 = "Ds79_Base.cpp";
	const char* path_Us79 = "Us79_Base.cpp";
#define MIN(a,b) (a<b?a:b)
#define MAX(a,b) (a>b?a:b)
	int min_bs = 2;
	int max_bs = MAX(CvxCompress::Max_BX(),MAX(CvxCompress::Max_BY(),CvxCompress::Max_BZ()));
	Gen_Ds79(path_Ds79,Find_Pow2(min_bs),Find_Pow2(max_bs),max_bs);
	Gen_Us79(path_Us79,Find_Pow2(min_bs),Find_Pow2(max_bs),max_bs);
#undef MAX
#undef MIN
	return 0;
}
