#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Vanilla_Visco_Elastic_CPU_Propagator.hxx"

void _generate_ricker_wavelet(double dt, double fmax, int* tsrc, double* stf)
{
        double fpeak = fmax / 2.0;
        double tshift = 5.0 * sqrt(1.5) / (fpeak * 3.1415926535897932384626433832795);
        //double tshift = 0.107;  // HACK
        *tsrc = (int)round((2.0 * tshift) / dt);
        for (int i = 0;  i < *tsrc;  ++i)
        {
                double t = (double)i * dt - tshift;
                double arg = -9.8696044010893586188344909998762 * fpeak * fpeak * t * t;
                stf[i] = (1.0 + 2.0 * arg) * exp(arg);
        }
        printf("SOURCE TYPE 2 : imax=%d\n",*tsrc);
}

void Compute_Time_Integrated_Source_Wavelet(
        double* stf,
        double* stf_int,
        int len,
        double dt
        )
{
        double dc = stf[0];
        for (int i = 1;  i < len;  ++i) dc += stf[i];
        dc /= (double)len;
        stf_int[0] = (stf[0] - dc) * dt;
        for (int i = 1;  i < len;  ++i) stf_int[i] = stf_int[i-1] + (stf[i] - dc) * dt;
}

int main(int argc, char* argv[])
{
	int bx = 64;
	int by = 32;
	if (argc == 3)
	{
		bx = atoi(argv[1]);
		by = atoi(argv[2]);
	}
	printf("Using block size %d by %d\n",bx,by);

	//Vanilla_Visco_Elastic_CPU_Propagator::Run_Unit_Tests();

	const int nx = 1760;
	const int ny = 1760;
	const int nz = 1040;
	const float dx = 6.25f;	
	const float dy = 6.25f;	
	const float dz = 6.25f;
	const float fq = 25.0f;	

	const float dti = 4.5e-4f;
	const float fmax = 60.0f;

	const int num_steps = 17550;

	const int stf_len = 32768;
	double* stf = new double[stf_len];
	memset((void*)stf, 0, stf_len*sizeof(double));
	int imax = 0;
	_generate_ricker_wavelet(dti,fmax,&imax,stf);
	double* stf_int = new double[stf_len];
	Compute_Time_Integrated_Source_Wavelet(stf,stf_int,stf_len,dti);

	Vanilla_Visco_Elastic_CPU_Propagator* prop = new Vanilla_Visco_Elastic_CPU_Propagator(nx,ny,nz,dx,dy,dz,fq);
	prop->Set_Tile_Size(bx,by);

	// populate earth model
	printf("Populating earth model...\n");
	for (int z = 0;  z < prop->Get_NZ();  ++z)
	{
		for (int y = 0;  y < prop->Get_NY();  ++y)
		{
			for (int x = 0;  x < prop->Get_NX();  ++x)
			{
				float depth = (float)z * dz;
				if (depth < 1000.0f)
				{
					prop->Set_Vp(x,y,z,1500.0f);
					prop->Set_Vs(x,y,z,0.0f);
					prop->Set_Rho(x,y,z,1000.0f);
					prop->Set_Q(x,y,z,1e6f);
				}
				else if (depth < 2000.0f)
				{
					prop->Set_Vp(x,y,z,2500.0f);
					prop->Set_Vs(x,y,z,1000.0f);
					prop->Set_Rho(x,y,z,2000.0f);
					prop->Set_Q(x,y,z,1e6f);
				}
				else
				{
					prop->Set_Vp(x,y,z,3500.0f);
					prop->Set_Vs(x,y,z,2000.0f);
					prop->Set_Rho(x,y,z,2500.0f);
					prop->Set_Q(x,y,z,1e6f);
				}
			}
		}
	}

	// insert pressure spike in the middle of volume
	int srcx = nx/2;
	int srcy = ny/2;
	int srcz = 20;

	// compute internal constants needed for propagation
	printf("Prepare for propagation.\n");
	prop->Prepare_For_Propagation();
	
	// do 100 timesteps
	printf("Propagation.\n");
	for (int i = 0;  i < num_steps;  ++i)
	{
		if ((i%10) == 0 && i > 0)
		{
			// write X-Z slice
			char str[4096];
			sprintf(str, "slices/xz_slice_y=%d_t=%d",ny/2,i);
			FILE* fp = fopen(str, "w");
			if (fp != 0L)
			{
				printf("Writing timestep to %s\n",str);
				for (int z = 0;  z < prop->Get_NZ();  ++z)
				{
					for (int x = 0;  x < prop->Get_NX();  ++x)
					{
						fprintf(fp, "%d %d %e\n",x,z,(prop->Txx(x,ny/2,z) + prop->Tyy(x,ny/2,z) + prop->Tzz(x,ny/2,z)) / 3.0f);
					}
					fprintf(fp, "\n");
				}
				fclose(fp);
			}
		}

		// insert pressure source
		if (i < imax)
		{
			prop->Set_Txx(srcx,srcy,srcz,prop->Txx(srcx,srcy,srcz)+0.33333333f*stf_int[i]);
			prop->Set_Tyy(srcx,srcy,srcz,prop->Tyy(srcx,srcy,srcz)+0.33333333f*stf_int[i]);
			prop->Set_Tzz(srcx,srcy,srcz,prop->Tzz(srcx,srcy,srcz)+0.33333333f*stf_int[i]);
		}

		double mpt_per_second = prop->Timestep(dti);
		printf("Timestep %4d/%4d - %.0f Mpt/s\n",i+1,num_steps,mpt_per_second);
	}

	return 0;
}

