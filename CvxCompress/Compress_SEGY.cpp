#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <SEGY_File.h>
#include <CvxCompress.hxx>
#include <map>

double Compute_Distance(SEGY_Trace_Header* hdr0, SEGY_Trace_Header* hdr1, bool common_shot)
{
	if (common_shot)
	{
		double delta_x = hdr0->Get_Rec_X() - hdr1->Get_Rec_X();
		double delta_y = hdr0->Get_Rec_Y() - hdr1->Get_Rec_Y();
		return sqrt(delta_x*delta_x+delta_y*delta_y);
	}
	else
	{
		double delta_x = hdr0->Get_Src_X() - hdr1->Get_Src_X();
		double delta_y = hdr0->Get_Src_Y() - hdr1->Get_Src_Y();
		return sqrt(delta_x*delta_x+delta_y*delta_y);
	}
}

int main(int argc, char* argv[])
{
	const double min_SNR = 10.0;
	const double max_SNR = 120.0;

	if (argc != 7)
	{
		printf("\nUsage: %s <segy-file> <segy-output-file> <SNR> <loc-RMS-flag> <bsamp> <btrace>\n",argv[0]);
		printf("<segy-file> is a valid SEGY file with at least read permission.\n");
		printf("<segy-output-file> is a valid file with write permission.\n");
		printf("<SNR> is the desired signal-to-noise ratio after compression, expressed in dB.\n");
		printf("<SNR> must be >= %.0fdB and <= %.0fdB.\n",min_SNR,max_SNR);
		printf("<loc-RMS-flag> is either 1 or 0. 1->Use local RMS.\n");
		printf("<bsamp> number of samples in block. Must be power of two between 8 and 256.\n");
		printf("<btrace> number of traces in block. Must be power of two between 8 and 256.\n");
		printf("\n");
		return -1;
	}

	CvxCompress* compressor = new CvxCompress();

	int bsamp = atoi(argv[5]);
	int btrace = atoi(argv[6]);
	if (!compressor->Is_Valid_Block_Size(bsamp,btrace,1))
	{
		printf("\nError: <bsamp>=%d, <btrace>=%d is not a valid block size.\n",bsamp,btrace);
		return -7;
	}

	bool use_local_RMS = false;
	if (strcmp(argv[4],"1") == 0) use_local_RMS = true;
	printf("Using %s RMS.\n",use_local_RMS?"local":"global");

	double SNR = atof(argv[3]);
	if (SNR < min_SNR)
	{
		printf("\nError: <SNR> is too small (%.2fdB). Must be >= %.2fdB.\n\n",SNR,min_SNR);
		return -3;
	}
	if (SNR > max_SNR)
	{
		printf("\nError: <SNR> is too large (%.2fdB). Must be <= %.2fdB.\n\n",SNR,max_SNR);
		return -4;
	}

	SEGY_File* segy = new SEGY_File(argv[1],false);
	if (!segy->Is_Valid())
	{
		printf("\nError: %s is not a valid SEGY file.\n\n",argv[1]);
		return -2;
	}

	int ntrace = segy->Get_Number_Of_Traces();
	if (ntrace < 1)
	{
		printf("\nError: SEGY file %s has no traces.\n\n",argv[1]);
		return -5;
	}
	printf("File has %d traces.\n",ntrace);
	SEGY_Trace_Header* hdr0 = segy->Get_Trace(0)->Get_Trace_Header();	
	SEGY_Trace_Header* hdr1 = segy->Get_Trace(1)->Get_Trace_Header();
	bool common_shot = (hdr0->Get_Rec_X() == hdr1->Get_Rec_X() && hdr0->Get_Rec_Y() == hdr1->Get_Rec_Y()) ? false : true;
	printf("This is a %s gather.\nTrying to determine frame size.\n",common_shot?"shot":"receiver");

	int prev_Trc = 0;
	double distance_estimate = Compute_Distance(hdr0,hdr1,common_shot);
	std::map<int,int> votes;
	for (int iTrc = 2;  iTrc < ntrace;  ++iTrc)
	{
		hdr0 = segy->Get_Trace(iTrc-1)->Get_Trace_Header();
		hdr1 = segy->Get_Trace(iTrc)->Get_Trace_Header();
		double distance = Compute_Distance(hdr0,hdr1,common_shot);
		if (distance > distance_estimate * 2.0)
		{
			int nTrc = iTrc - prev_Trc;
			prev_Trc = iTrc;
			if (votes.find(nTrc) != votes.end())
				votes[nTrc]++;
			else
				votes[nTrc] = 1;
		}
	}
	int frame_size = 0;
	int frame_size_votes = 0;
	for (auto const &ent : votes)
	{
		int nTrc = ent.first;
		int nvotes = ent.second;
		if (nvotes > frame_size_votes)
		{
			frame_size = nTrc;
			frame_size_votes = nvotes;
		}
	}
	if (frame_size <= 0)
	{
		printf("Unable to determine frame size.\n");
	}
	else
	{
		printf("Most likely frame size is %d traces.\n",frame_size);
		if ((ntrace % frame_size) == 0)
		{
			printf("File had exactly %d frames.\n",ntrace/frame_size);
		}
		else
		{
			printf("WARNING! Number of traces is not a multiple of the likely frame size.\n");
		}
	}

	int nsamp = segy->Get_Trace(0)->Get_Trace_Header()->Get_NSAMP();

	// allocate single buffer for all traces
	float* samples;
	posix_memalign((void**)&samples, 64, (long)sizeof(float)*(long)ntrace*(long)nsamp);
	if (samples == 0L)
	{
		printf("\nError: Out of memory, could not allocate trace buffer.\n\n");
		return -6;
	}
	// allocate compress buffer
	unsigned int* compressed;
	posix_memalign((void**)&compressed, 64, (long)sizeof(float)*(long)ntrace*(long)nsamp);
	if (compressed == 0L)
	{
		printf("\nError: Out of memory, could not allocate compress buffer.\n\n");
		return -6;
	}
	// allocate decompress buffer
	float* decompressed;
	posix_memalign((void**)&decompressed, 64, (long)sizeof(float)*(long)ntrace*(long)nsamp);
	if (decompressed == 0L)
	{
		printf("\nError: Out of memory, could not allocate decompress buffer.\n\n");
		return -6;
	}
	
	printf("Reading SEGY file.\n");
	for (int iTrc = 0;  iTrc < ntrace;  ++iTrc)
	{
		SEGY_Trace* trace = segy->Get_Trace(iTrc);
		long idx = (long)iTrc * (long)nsamp;
		memcpy(samples+idx,trace->Get_Samples(),sizeof(float)*nsamp);
	}
	
	printf("Compressing. This is an iterative process to find parameters that yield %.2fdB SNR.\n",SNR);
	double desired_error = pow(10.0,-(SNR/20.0));
	double scale = desired_error;
	double SNR_diff = 0.0;
	long length = 0;
	do
	{
		float ratio = compressor->Compress(scale,samples,nsamp,ntrace,1,bsamp,btrace,1,use_local_RMS,compressed,length);
		compressor->Decompress(decompressed,nsamp,ntrace,1,compressed,length);
		double acc0=0.0, acc1=0.0;
#pragma omp parallel for reduction(+:acc0,acc1)
		for (int iTrc = 0;  iTrc < ntrace;  ++iTrc)
		{
			long idx = (long)iTrc * (long)nsamp;
			for (long isamp = 0;  isamp < nsamp;  ++isamp)
			{
				double val0 = samples[idx+isamp];
				double val1 = decompressed[idx+isamp] - val0;
				acc0 += val0 * val0;
				acc1 += val1 * val1;
			}
		}
		double actual_error = sqrt(acc1)/sqrt(acc0);
		double actual_SNR = -20.0 * log10(actual_error);
		SNR_diff = fabs(SNR-actual_SNR);
		printf(" -> scale=%.4e, ratio %.2f:1, achieved SNR %.2fdB\n",scale,ratio,actual_SNR);
		scale *= (desired_error / actual_error);
	} while (SNR_diff > 4e-3);

	printf("Writing decompressed data to %s.\n",argv[2]);
	for (int iTrc = 0;  iTrc < ntrace;  ++iTrc)
        {
                SEGY_Trace* trace = segy->Get_Trace(iTrc);
		long idx = (long)iTrc * (long)nsamp;
		memcpy(trace->Get_Samples(),decompressed+idx,sizeof(float)*nsamp);
	}
	segy->Write(argv[2]);

	delete compressor;
	free(decompressed);
	free(compressed);
	free(samples);
	return 0;
}
