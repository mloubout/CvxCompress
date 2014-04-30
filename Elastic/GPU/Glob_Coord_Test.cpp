#include <stdio.h>

#include "Elastic_Modeling_Job.hxx"
#include "Voxet.hxx"
#include "Global_Coordinate_System.hxx"

int main(int argc, char** argv)
{
	int log_level = 4;
	
	Elastic_Modeling_Job* eJob = new Elastic_Modeling_Job(log_level,argv[1]);
	if (!eJob->Is_Valid())
	{
		printf("Invalid parameter file.\n");
		return -1;
	}

	Voxet* voxet = eJob->Get_Voxet();
	voxet->Dump();
	printf("\n");

	Global_Coordinate_System* gcs = voxet->Get_Global_Coordinate_System();
	printf("\n");

	double g0 = 266840.0;
	double g1 = 7716425.0;
	double g2 = -12.5;
	printf("Global coordinate [%lf,%lf,%lf]\n",g0,g1,g2);

	double nlu,nlv,nlw;
	gcs->Convert_Global_To_Normalized_Local(g0,g1,g2,nlu,nlv,nlw);
	printf("Normalized local [%lf,%lf,%lf]\n",nlu,nlv,nlw);

	double lu,lv,lw;
	gcs->Convert_Normalized_Local_To_Local(nlu,nlv,nlw,lu,lv,lw);
	printf("Local [%lf,%lf,%lf]\n",lu,lv,lw);

	double filu,filv,filw;
	gcs->Convert_Normalized_Local_To_Fractional_Local_Index(nlu,nlv,nlw,filu,filv,filw);
	printf("Fractional index [%lf,%lf,%lf]\n",filu,filv,filw);

	int ilu,ilv,ilw;
	gcs->Convert_Fractional_Local_Index_To_Local_Index(filu,filv,filw,ilu,ilv,ilw);
	printf("Local index [%d,%d,%d]\n",ilu,ilv,ilw);

	int ix,iy,iz;
	gcs->Convert_Local_Index_To_Transposed_Index(ilu,ilv,ilw,ix,iy,iz);
	printf("Transposed index [%d,%d,%d]\n",ix,iy,iz);

	double h0,h1,h2;
	gcs->Convert_Transposed_Index_To_Global(ix,iy,iz,h0,h1,h2);
	printf("Approximate global from transposed index [%lf,%lf,%lf]\n",h0,h1,h2);
	
	double x,y,z;
	gcs->Convert_Global_To_Transposed_Fractional_Index(g0,g1,g2,x,y,z);
	printf("Transposed fractional index [%lf,%lf,%lf]\n",x,y,z);

	double j0,j1,j2;
	gcs->Convert_Transposed_Fractional_Index_To_Global(x,y,z,j0,j1,j2);
	printf("Exact global from transposed fractional index [%lf,%lf,%lf]\n",j0,j1,j2);

	return 0;
}

