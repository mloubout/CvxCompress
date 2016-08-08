#include <assert.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <ElaOrthModel.hxx>
#include <Voxet.hxx>
#include <Global_Coordinate_System.hxx>

using namespace std;

int main(int argc, char* argv[])
{
	char buf[1000000];
	ifstream fs(argv[1]);
	fs.read(buf,1000000);
	ElaOrthModel* m = new ElaOrthModel(buf);

	int nsamp;
	double sample_rate;
	m->Get_Trace_Arguments(nsamp,sample_rate);

	Voxet* voxet = m->Get_Voxet();
	Global_Coordinate_System* gcs = voxet->Get_Global_Coordinate_System();

	int nx = gcs->Get_NX();
	int ny = gcs->Get_NY();
	int nz = gcs->Get_NZ();

	if (argc == 2)
	{
		int nTraces = 4;
		int nShots = 10;

		float* traces = new float[nsamp*nTraces];
		float* recx = new float[nTraces]; 
		float* recy = new float[nTraces]; 
		float* recz = new float[nTraces];

		srand48(time(0L));
		for (int iShot = 0;  iShot < nShots;  ++iShot)
		{
			double slx = drand48() * (double)nx * gcs->Get_DX();
			double sly = drand48() * (double)ny * gcs->Get_DY();
			double slz = 2.0 * gcs->Get_DZ();
			double soux,souy,souz;
			gcs->Convert_Transposed_Fractional_Index_To_Global(slx/gcs->Get_DX(),sly/gcs->Get_DY(),slz/gcs->Get_DZ(),soux,souy,souz);
			double rcx,rcy,rcz;
			gcs->Convert_Transposed_Fractional_Index_To_Global((slx-1000)/gcs->Get_DX(),(sly-1000)/gcs->Get_DY(),slz/gcs->Get_DZ(),rcx,rcy,rcz);
			recx[0]=rcx;  recy[0]=rcy;  recz[0]=rcz;
			gcs->Convert_Transposed_Fractional_Index_To_Global((slx+1000)/gcs->Get_DX(),(sly-1000)/gcs->Get_DY(),slz/gcs->Get_DZ(),rcx,rcy,rcz);
			recx[1]=rcx;  recy[1]=rcy;  recz[1]=rcz;
			gcs->Convert_Transposed_Fractional_Index_To_Global((slx+1000)/gcs->Get_DX(),(sly+1000)/gcs->Get_DY(),slz/gcs->Get_DZ(),rcx,rcy,rcz);
			recx[2]=rcx;  recy[2]=rcy;  recz[2]=rcz;
			gcs->Convert_Transposed_Fractional_Index_To_Global((slx-1000)/gcs->Get_DX(),(sly+1000)/gcs->Get_DY(),slz/gcs->Get_DZ(),rcx,rcy,rcz);
			recx[3]=rcx;  recy[3]=rcy;  recz[3]=rcz;
			for (int iTrc = 0;  iTrc < nTraces;  ++iTrc)
			{
				printf("rec[%d] = %.2f, %.2f, %.2f\n",iTrc,recx[iTrc],recy[iTrc],recz[iTrc]);
			}
			m->runShot(nTraces,soux,souy,souz,recx,recy,recz,traces);
		}
	}
	else if (argc == 3)
	{
		int nShots = 10;
		for (int iShot = 0;  iShot < nShots;  ++iShot)
		{
			char str[1024];
			ifstream fs2(argv[2]);
			int line = 0;
			int nTraces;
			double x, y, z;
			float soux, souy, souz;
			float *traces, *recx, *recy, *recz;
			for (fs2.getline(str,1024);  !fs2.eof();  fs2.getline(str,1024))
			{
				if (line == 0)
				{
					// source coordinates
					assert(sscanf(str, "sou x,y,e = %lf %lf %lf", &x, &y, &z) == 3);
					soux = x;
					souy = y;
					souz = z;
				}
				else if (line == 1)
				{
					// number of receivers
					assert(sscanf(str, "num_rec = %d", &nTraces) == 1);
					assert(nTraces > 0);
					traces = new float[nsamp*nTraces];
					recx = new float[nTraces]; 
					recy = new float[nTraces]; 
					recz = new float[nTraces];
				}
				else
				{
					// one receiver
					int iTrc;
					assert(sscanf(str, "rec i, x,y,e = %d %lf %lf %lf", &iTrc, &x, &y, &z) == 4);
					recx[iTrc] = x;
					recy[iTrc] = y;
					recz[iTrc] = z;
				}
				++line;
			}
			m->runShot(nTraces,soux,souy,souz,recx,recy,recz,traces);
			delete [] traces;
			delete [] recx;
			delete [] recy;
			delete [] recz;
		}
	}
	return 0;
}
