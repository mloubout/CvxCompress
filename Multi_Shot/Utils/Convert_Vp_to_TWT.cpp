#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <libgen.h>
#include <swapbytes.h>

void Print_Usage(const char* cmd)
{
	printf("\n");
	printf("Usage : %s <svp-bin-file> <dz> <zwb>\n",cmd);
	printf("<svp-bin-file> File containing the svp binary data.\n");
	printf("<dz> Vertical cell size in meters.\n");
	printf("<zwb> Water depth in meters.\n");
	printf("\n");
}

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		Print_Usage(argv[0]);
		return -1;
	}

	std::string filepath = (std::string)(argv[1]);
	std::string filename = (std::string)basename((char*)filepath.c_str());

	int ffid = 0;
	if (sscanf(filename.c_str(),"ffid_%d.src.Vp.bin",&ffid) < 1)
	{
		printf("Error! Unable to determine FFID from filename %s.\n",argv[1]);
		Print_Usage(argv[0]);
                return -8;
	}

	double dz = atof(argv[2]);
	if (dz <= 0.0)
	{
		printf("Error! <dz> Must be positive number.\n");
		Print_Usage(argv[0]);
		return -2;
	}
		
	double zwb = atof(argv[3]);
	int iwb = (int)trunc(zwb/dz);
	if (iwb <= 0)
	{
		printf("Error! <zwb> Must be positive number >= <dz>.\n");
		Print_Usage(argv[0]);
		return -3;
	}
	
	struct stat st;
	if (stat(argv[1],&st) != 0)
	{
		printf("Error! Unable to stat file %s.\n",argv[1]);
		Print_Usage(argv[0]);
		return -4;
	}
	int nsamp = st.st_size / sizeof(float);
	if (nsamp <= 0)
	{
		printf("Error! File %s is empty.\n",argv[1]);
		Print_Usage(argv[0]);
                return -5;
	}
	if (iwb >= nsamp)
	{
		printf("Error! <zwb> is beyond depth of <svp-bin-file> %s.\n",argv[1]);
		Print_Usage(argv[0]);
                return -6;
	}

	float* samples = new float[nsamp];
	FILE* fp = fopen(filepath.c_str(),"rb");
	int nread = fread((void*)samples,sizeof(float),nsamp,fp);
	fclose(fp);
	if (nread < nsamp)
	{
		printf("Error! Unable to read all %d values from <svp-bin-file> %s.\n",nsamp,argv[1]);
		Print_Usage(argv[0]);
                return -7;
	}
	swap4bytes((int*)samples,nsamp);

	double owt = 0.0;
	for (int i = 0;  i < iwb;  ++i)
	{
		owt = owt + (dz / samples[i]);
	}
	double avg_Vp = ((double)iwb * dz) / owt;
	double twt = 2.0 * zwb / avg_Vp;

	printf("%d %.2f %.2lf\n",ffid,zwb,twt*1e3);

	delete [] samples;
	return 0;
}
