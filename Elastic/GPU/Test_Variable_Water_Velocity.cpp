#include <cassert>
#include <cstdio>
#include <string>
#include <Variable_Water_Velocity.hxx>
#include <Voxet.hxx>
#include <Voxet_Property.hxx>
#include <Voxet_Memory_Mapper.hxx>
#include <Global_Coordinate_System.hxx>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int main(int argc, char* argv[])
{
	if (argc != 5)
	{
		printf("Usage : %s <voxet> <Vs-property> <uvw-mapping> <wbIdx-file>\n",argv[0]);
		return -1;
	}

	const int log_level = 4;
	setenv("TZ", "UTC", 1);

	Voxet* voxet = new Voxet(log_level,argv[1]);
	Global_Coordinate_System* gcs = voxet->Get_Global_Coordinate_System();
	gcs->Set_Transpose(argv[3]);
	Voxet_Property* Vs_Prop = voxet->Get_Property_By_Moniker(argv[2]);
	if (Vs_Prop == 0L)
	{
		printf("Error! %s is not a property in voxet %s.\n",argv[2],argv[1]);
		printf("This voxet has the following properties:\n");
		for (int i = 0;  i < voxet->Get_Number_Of_Properties();  ++i)
		{
			Voxet_Property* prop = voxet->Get_Property_By_Index(i);
			printf("%s\n",prop->Get_Moniker());
		}
		return -2;
	}
	long nxny = (long)gcs->Get_NX() * (long)gcs->Get_NY();
	
	Voxet_Memory_Mapper* mapper = new Voxet_Memory_Mapper();
	float* Vs = mapper->Get_Memory_Mapped_File(Vs_Prop->Get_Full_Path());
	if (Vs == 0L)
	{
		printf("Unable to memory map file %s\n",Vs_Prop->Get_Full_Path());
		return -3;
	}

	printf("Creating interpolation volumes.\n");
	std::map<time_t,std::string> VwInterp;
	const long VwNz = 161;
	float* buf = new float[VwNz];
	for (int i = 0;  i < 4;  ++i)
	{
		for (int j = 0;  j < VwNz;  ++j) buf[j] = (float)(i+1);
		char str[256];
		sprintf(str,"VwInterp%d",i);
		FILE* fp = fopen(str,"wb");
		for (int y = 0;  y < gcs->Get_NY();  ++y)
		{
			for (int x = 0;  x < gcs->Get_NX();  ++x)
			{
				for (int z = 0;  z < VwNz;  ++z)
				{
					buf[z] = (float)(i+1);// + (float)x + (float)y + (float)z;
				}
				fwrite(buf,sizeof(float),VwNz,fp);
			}
		}
		fclose(fp);
		struct tm st;
		st.tm_year = 116;
		st.tm_mon = 2+i;
		st.tm_mday = 1;
		st.tm_hour = 0;
		st.tm_min = 1;
		st.tm_sec = 0;
		time_t start_time = mktime(&st);
		VwInterp[start_time] = (std::string)str;
	}
	delete [] buf;

	for (std::map<time_t,std::string>::iterator it = VwInterp.begin();  it != VwInterp.end();  ++it)
	{
		printf("VwInterp %s",asctime(localtime(&(it->first))));
	}
	
	printf("Reading in interpolation volumes.\n");
	Variable_Water_Velocity* vwv = new Variable_Water_Velocity(gcs,Vs,VwInterp);
	int* wbIdx = vwv->Get_Water_Bottom();
	FILE* fp = fopen(argv[4],"w");
	if (fp == 0L)
	{
		printf("Unable to open %s for writing.\n",argv[3]);
	}
	for (int y = 0;  y < gcs->Get_NY();  ++y)
	{
		for (int x = 0;  x < gcs->Get_NX();  ++x)
		{
			fprintf(fp,"%d %d %d\n",x,y,wbIdx[y*gcs->Get_NX()+x]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	printf("Done!\n");

	for (int i = -1;  i < 5;  ++i)
	{
		struct tm st;
                st.tm_year = 116;
                st.tm_mon = 2+i;
                st.tm_mday = 15;
		st.tm_hour = 0;
		st.tm_min = 1;
		st.tm_sec = 0;
                time_t start_time = mktime(&st);
		vwv->Set_Shot_Time(start_time);
		int x = gcs->Get_NX()/2;
		int y = gcs->Get_NY()/2;
		int z = 5;
		float Vw = vwv->Compute_Velocity_Increment(x,y,z);// - (float)x - (float)y - (float)z;
		printf("Velocity increment is %e at time %s",Vw,asctime(&st));
	}

	return 0;
}
