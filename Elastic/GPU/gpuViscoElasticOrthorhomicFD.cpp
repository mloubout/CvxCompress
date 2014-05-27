#include <stdio.h>
#include <unistd.h>

#include <cuda_runtime_api.h>

#include "Elastic_Propagator.hxx"
#include "Elastic_Modeling_Job.hxx"
#include "Elastic_Shot.hxx"
#include "Elastic_Buffer.hxx"

int main(int argc, char* argv[])
{
	printf("\n%s v0.8\nVariable density visco-elastic orthorhombic finite difference forward modeling.\n\n",argv[0]);

#ifdef GPU_DEBUG
	printf("Warning! This executable was compiled with -DGPU_DEBUG.\n\n");
	printf("GPU kernels will be serialized and performance will be slow.\n\n");
#endif

	int log_level = 4;

	Elastic_Modeling_Job* job = new Elastic_Modeling_Job(log_level, argv[1]);
	if (job->Is_Valid())
	{
		Elastic_Propagator* prop = new Elastic_Propagator(job);
		prop->Print_Graphical();

		printf("Allocating host memory...\n");
		prop->Allocate_Host_Memory(false, true);
		printf("Done!\n");

		prop->Read_Earth_Model();
		job->Write_Earth_Model_XZ_Slice("slices/xz_slice",job->Get_Propagation_NY()/2);
		job->Write_Earth_Model_XY_Slice("slices/xy_slice",100);

		printf("Allocating device memory...\n");
		prop->Allocate_Device_Memory();
		printf("Done!\n");

		//prop->Set_WF_Value(6,500,400,400,1.0f);
		//prop->Set_WF_Value(7,500,400,400,1.0f);
		//prop->Set_WF_Value(8,500,400,400,1.0f);
		//prop->Set_WF_Value(9,403,400,407,1.0f);
		//prop->Set_WF_Value(10,403,400,407,1.0f);
		//prop->Set_WF_Value(11,403,400,407,1.0f);
	
		/*	
		prop->Set_WF_Value(0,401,400,400,1.0f);
		prop->Set_WF_Value(1,401,400,400,1.0f);
		prop->Set_WF_Value(2,401,400,400,1.0f);
		*/

		//prop->Set_WF_Value(0,200,200,200,1.0f);
		//prop->Set_WF_Value(1,200,200,200,2.0f);
		//prop->Set_WF_Value(2,200,200,200,3.0f);
		for (int iShot = 0;  iShot < job->Get_Number_Of_Shots();  ++iShot)
		{
			/*
			Elastic_Shot* shot = job->Get_Shot_By_Index(iShot);
			while (!prop->Propagate_One_Block(prop->Get_Total_Number_Of_Timesteps(), shot));
			break;

			prop->Prepare_For_Propagation(shot);
			int nn = prop->Get_Total_Number_Of_Timesteps();
			for (int i = 0;  i < 200;  ++i)
			{
				while (!prop->Propagate_One_Block((i+1)*nn,shot));
				char path[4096];
				sprintf(path,"Gorgon_P_XZ_Slice_%04d",(i+1)*nn);
				job->Write_XZ_Slice(path, 3, (int)(shot->Get_Propagation_Source_Y()/2.0));
			}
			//prop->Propagate_Shot(shot);
			*/
			
			Elastic_Shot* shot = job->Get_Shot_By_Index(iShot);
			prop->Propagate_Shot(shot);
		}

		printf("Freeing host memory...\n");
		prop->Free_Host_Memory();
		printf("Done!\n");

		printf("Freeing device memory...\n");
		prop->Free_Device_Memory();
		printf("Done!\n");

		delete prop;
	}
	delete job;
	return 0;
}
