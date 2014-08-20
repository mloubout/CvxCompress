#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include <cuda_runtime_api.h>

#include "Elastic_Propagator.hxx"
#include "Elastic_Modeling_Job.hxx"
#include "Elastic_Shot.hxx"
#include "Elastic_Buffer.hxx"

int main(int argc, char* argv[])
{
	printf("\n%s v0.9\nVariable density visco-elastic orthorhombic finite difference forward modeling.\n\n",argv[0]);
	if (argc == 1)
	{
		printf("Usage : %s <parmfile> [options]\n",argv[0]);
		printf("Valid options are:\n");
		printf("xz_model	Output X-Z crossplot of earth model intersecting source location.\n");
		printf("xz_slices	Output X-Z crossplot of P wavefield intersecting source location during propagation.\n");
		printf("source_wavelet	Output filtered and time integrated source wavelet.\n");
		return -1;
	}

#ifdef GPU_DEBUG
	printf("Warning! This executable was compiled with -DGPU_DEBUG.\n\n");
	printf("GPU kernels will be serialized and performance will be slow.\n\n");
#endif

	int log_level = 3;

	Elastic_Modeling_Job* job = new Elastic_Modeling_Job(log_level, argv[1]);
	if (job->Is_Valid())
	{
		bool xz_model = false, xz_slices = false, source_wavelet = false;
		for (int argi = 2;  argi < argc;  ++argi)
		{
			if (strcmp(argv[argi], "xz_model") == 0) xz_model = true;
			if (strcmp(argv[argi], "xz_slices") == 0) xz_slices = true;
			if (strcmp(argv[argi], "source_wavelet") == 0) source_wavelet = true;
		}

		Elastic_Propagator* prop = new Elastic_Propagator(job);
		if (prop != 0L)
		{
			prop->Configure();

			prop->Read_Earth_Model();
			if (job->Lower_Q_Seafloor_Enabled()) job->Lower_Q_Seafloor();
			if (xz_model && job->Get_Number_Of_Shots() > 0)
			{
				Elastic_Shot* shot = job->Get_Shot_By_Index(0);
				job->Write_Earth_Model_XZ_Slice("slices/xz_slice",(int)round(shot->Get_Propagation_Source_Y()));
				//job->Write_Earth_Model_XY_Slice("slices/xy_slice",100);
			}

			for (int iShot = 0;  iShot < job->Get_Number_Of_Shots();  ++iShot)
			{
				Elastic_Shot* shot = job->Get_Shot_By_Index(iShot);
				prop->Propagate_Shot(shot,source_wavelet,xz_slices);
			}

			printf("Freeing host memory...\n");
			prop->Free_Host_Memory();
			printf("Done!\n");

			printf("Freeing device memory...\n");
			prop->Free_Device_Memory();
			printf("Done!\n");

			delete prop;
		}
	}
	delete job;
	return 0;
}

