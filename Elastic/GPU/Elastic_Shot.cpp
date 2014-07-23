#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "Elastic_Shot.hxx"
#include "Elastic_Propagator.hxx"
#include "Elastic_Modeling_Job.hxx"
#include "Elastic_Buffer.hxx"
#include "Elastic_SEGY_File.hxx"
#include "Elastic_Pipeline.hxx"
#include "Extract_Receiver_Values.hxx"
#include "Voxet.hxx"
#include "Global_Coordinate_System.hxx"
#include "DFT.hxx"
#include "gpuAssert.h"

//#define RESAMPLE_DEBUG

Elastic_Shot::Elastic_Shot(int log_level, Elastic_Modeling_Job* job, int souidx, double x, double y, double z)
{
	_log_level = log_level;
	_job = job;
	_ordertime = 2;  // hardcode for now.
	_souintrp = Point;
	_souidx = souidx;
	_x = x;
	_y = y;
	_z = z;
	_soutype = 0;
	_ampl1 = 0.0;
	_ampl2 = 0.0;
	_ampl3 = 0.0;
	_wavetype = 0;
	_max_freq = 10.0;
	_wavelet_path = 0L;
	_filter_order = 7;
	_tsrc = 0;
	_stf = new double[65536];
	_stf_int = new double[65536];
	_segy_files = 0L;
	_num_segy_files = 0;
        _totSteps = 0;
        _nBlks = 0;
        _num_pipes = 0;
        _h_rcv_loc = 0L;
        _h_pinned_rcv_loc = 0L;
        _h_rcv_binned = 0L;
        _h_rcv_trcidx = 0L;
        _h_rcv_trcflag = 0L;
        _h_rcv_loc_size_f = 0L;
	_num_traces = 0;
	_h_trace_rcv_x = 0L;
	_h_trace_rcv_y = 0L;
	_h_trace_rcv_z = 0L;
	_h_trace_iline = 0L;
	_h_trace_xline = 0L;
	_h_trace_trcens = 0L;
	_h_trace_tshift = 0L;
	_h_trace_sample_rate = 0L;
	_h_trace_in = 0L;
	_h_trace_out = 0L;
	_h_trace_flag = 0L;
	_h_trace_idx_in = 0L;
	_h_trace_idx_in_nn = 0L;
	_h_trace_idx_out = 0L;
	_h_trace_iFile = 0L;
	_h_trace_touched = 0L;
	_h_trace_nsamp_in = 0L;
	_h_trace_nsamp_out = 0L;
	_h_trace_out_idxM = 0L;
	_h_trace_out_sinc_coeffs = 0L;
}

Elastic_Shot::~Elastic_Shot()
{
	Free_Trace_Resample_Buffers();
	if (_stf != 0L) delete [] _stf;
	if (_stf_int != 0L) delete [] _stf_int;
	if (_num_segy_files > 0)
	{
		for (int i = 0;  i < _num_segy_files;  ++i)
		{
			delete _segy_files[i];
			for (int idx = 0;  idx < _h_trace_nsamp_out[i];  ++idx)
			{
				delete [] _h_trace_out_sinc_coeffs[i][idx];
			}
			delete [] _h_trace_out_sinc_coeffs[i];
			delete [] _h_trace_out_idxM[i];
		}
		delete [] _h_trace_out_sinc_coeffs;
		delete [] _h_trace_out_idxM;
		delete [] _h_trace_nsamp_in;
		delete [] _h_trace_nsamp_out;
		delete [] _segy_files;
	}
	if (_wavelet_path != 0L) free(_wavelet_path);
}

double Elastic_Shot::Get_Propagation_Source_X() 
{
	return _x - _job->Get_Propagation_X0();
}

double Elastic_Shot::Get_Propagation_Source_Y() 
{
	return _y - _job->Get_Propagation_Y0();
}

double Elastic_Shot::Get_Propagation_Source_Z() 
{
	return _z - _job->Get_Propagation_Z0();
}

double Elastic_Shot::Get_Propagation_Time()
{
	double propagation_time = 0.0;
	for (int iFile = 0;  iFile < _num_segy_files;  ++iFile)
	{
		double file_propagation_time = _segy_files[iFile]->Get_Timeshift() + _segy_files[iFile]->Get_Record_Length();	
		if (file_propagation_time > propagation_time) propagation_time = file_propagation_time;
	}
	return propagation_time;
}

void Elastic_Shot::Add_Receiver_Array(int nrec,	double* rec_x, double* rec_y, double* rec_z, int* iline, int* xline, int* trcens ) {
	
	for (int iFile = 0;  iFile < _num_segy_files;  ++iFile){
		_segy_files[iFile]->Add_Receiver_Array(nrec,rec_x,rec_y,rec_z, iline, xline, trcens);
	}
}

const char* Elastic_Shot::Get_Source_Type_String()
{
	if (_soutype == Source_Type_Force)
	{
		return "Force";
	}
	else if (_soutype == Source_Type_Velocity)
	{
		return "Velocity";
	}
	else if (_soutype == Source_Type_Pressure)
	{
		return "Pressure";
	}
	else
	{
		return "Unknown";
	}
}

void Elastic_Shot::Prepare_Source_Wavelet(double dt, bool debug_output_source_wavelet)
{
	if (_wavetype == 1)
	{
		_src(dt, _max_freq, 1, 0L, &_tsrc, _stf);
		//for (int i = 0;  i < _tsrc;  ++i) _stf[i] = 0.0f;
		//_stf[3] = 1.0f;
	}
	else if (_wavetype == 2)
	{
		_generate_ricker_wavelet(dt, _max_freq, &_tsrc, _stf);
	}
	else if (_wavetype == 3)
	{
		_src(dt, _max_freq, 2, _wavelet_path, &_tsrc, _stf);
	}
	// find peak
	int ipeak = 0;
	float peak_val = _stf[0];
	for (int i = 1;  i < _tsrc;  ++i) if (_stf[i] > peak_val) {ipeak = i;  peak_val = _stf[i];}
	int ichoplo = ipeak;
	while (_stf[ichoplo-1] > 0.0f) --ichoplo;
	int ichophi = ipeak;
	int nzc = 0;  // number of zero crossings
	while (nzc < 2)
	{
		if (_stf[ichophi] > 0.0f && _stf[ichophi+1] <= 0.0f)
		{
			++nzc;
		}
		++ichophi;
	};
	printf("ipeak = %d, ichoplo = %d, ichophi = %d\n",ipeak,ichoplo,ichophi);
	//for (int i = ichoplo-1;  i >= 0;  --i) _stf[i] = 0.0f;
	//for (int i = ichophi+1;  i < _tsrc;  ++i) _stf[i] = 0.0f;

	Compute_Time_Integrated_Source_Wavelet(_log_level,_stf,_stf_int,_tsrc,dt);
	if (debug_output_source_wavelet)
	{
		FILE* fp = fopen("filtered.txt", "w");
		if (fp != 0L)
		{
			printf("Writing filtered source wavelet to file filtered.txt\n");
			for (int i = 0;  i < _tsrc;  ++i) fprintf(fp, "%e %e\n",(double)i*dt,_stf[i]);
			fclose(fp);
		}
		fp = fopen("filtered_int.txt", "w");
		if (fp != 0L)
		{
			printf("Writing filtered and time integrated source wavelet to file filtered_int.txt\n");
			for (int i = 0;  i < _tsrc;  ++i) fprintf(fp,"%e %e\n",(double)i*dt,_stf_int[i]);
			fclose(fp);
		}
	}
}

bool Elastic_Shot::Read_Source_Wavelet_From_File(const char* wavelet_path, double max_freq, int filter_order)
{
	FILE* fp = fopen(wavelet_path, "rb");
	if (fp != 0L)
	{
		fclose(fp);
		_wavelet_path = strdup(wavelet_path);
		_max_freq = max_freq;
		_filter_order = filter_order;
		_wavetype = 3;
		return false;
	}
	else
	{
		return true;
	}
}

bool Elastic_Shot::Use_Builtin_Source_Wavelet(const char* wavetype, double max_freq, const char* parmfile_path, int line_num)
{
	if (strcmp(wavetype, "gaussian") == 0)
	{
		_wavetype = 1;
		_max_freq = max_freq;
		return false;
	}
	else if (strcmp(wavetype, "ricker") == 0)
	{
		_wavetype = 2;
		_max_freq = max_freq;
		return false;
	}
	else
	{
		printf("%s (line %d): Error - SOURCE_WAVELET unknown wavelet type %s.\n", parmfile_path, line_num, wavetype);
		return true;
	}
}

void Elastic_Shot::_generate_ricker_wavelet(double dt, double fmax, int* tsrc, double* stf)
{
	if (_log_level >= 4)
	{
		printf("_generate_ricker_wavelet(dt=%e, fmax=%e, tsrc=%s, stf=%s)\n",dt,fmax,tsrc!=0L?"ok":"nil",stf!=0L?"ok":"nil");
		fflush(stdout);
	}
	double fpeak = fmax / 2.0;
	double tshift = 5.0 * sqrt(1.5) / (fpeak * 3.1415926535897932384626433832795);
//	double tshift = 0.107;  // HACK
	*tsrc = (int)round((2.0 * tshift) / dt);
	for (int i = 0;  i < *tsrc;  ++i)
	{
		double t = (double)i * dt - tshift;
		double arg = -9.8696044010893586188344909998762 * fpeak * fpeak * t * t;
		stf[i] = (1.0 + 2.0 * arg) * exp(arg);
	}
	if (_log_level >= 4) printf("SOURCE TYPE 2 : imax=%d\n",*tsrc);
}

void Elastic_Shot::_src(double dt, double fmax, int type, char* stfname, int* tsrc, double* stf)
{
        if (_log_level >= 4)
        {
                printf("src(_log_level=%d, dt=%e, fmax=%e, type=%d, stfname=%s, tsrc=%s, stf=%s)\n",_log_level,dt,fmax,type,stfname!=0L?stfname:"nil",tsrc!=0L?"ok":"nil",stf!=0L?"ok":"nil");
                fflush(stdout);
        }
        if (type == 1)
        {
                /* type==1: first derivative of a Gaussian, with linear extension tapers */
                int i, imax, imaxhalf, ntap;
                double w0, wmax, ts,t,rt, wt, wt2, diff;

                wmax = 2.0*M_PI*fmax;
                /* Note:   tsrc = ts/dt = 2/(gam*khmax)*(rt*rt)*(Vmax/Vmin) */

                /* if(type==1) */  /* only one type for now */
                {
                        rt = 3.571625; /* guarantees SourceAmplitude(tmax)   = 0.01*MaxAmplitude
                                           and guarantees SpectralAmplitude(fmax) = 0.01*MaxSpectrum
                                           and: wmax/w0 = rt = w0*ts/2  */
                        w0 = wmax/rt;  /* w0i = 1./w0; */
                        ts = 2.0*rt/w0;  /* total source time */
                        imax = (int)(ts/dt) + 1;

                        for(i=0;i<imax;i++)
                        { t=i*dt-0.5*ts;
                                stf[i]  = -sqrt(M_E)*w0*t*exp(-0.5*w0*w0*t*t);
				//printf("stf[%d] = %e\n",i,stf[i]);
                        }

                        /* taper (linearly extend) front and back ends */
                        /* front end */
                        diff = stf[1]-stf[0];
                        ntap = (int)(fabs(stf[0]/diff));
                        for(i=imax-1; i>=0; i--) stf[i+ntap] = stf[i];  // shift
                        for(i=ntap-1; i>=0; i--) stf[i] = stf[i+1] - diff; // taper
                        imax += ntap;

                        /* back end: */
                        diff = stf[imax-1]-stf[imax-2];
                        ntap = (int)(fabs(stf[imax-1]/diff));
                        for(i=0; i<ntap; i++)  stf[imax+i] = stf[imax+i-1] + diff; // taper
                        imax += ntap;
                }

                *tsrc = imax;

                if (_log_level >= 4) printf("SOURCE TYPE 1 : imax=%d\n",imax);

                // for(i=0; i<imax; i++) stf[i]=0.0f; stf[0] = 1.0f;

                // for(i=0; i<imax; i++) printf("%d  %f\n", i,stf[i]);
        }
        else if (type == 2)
        {
                /* type==2: user provided source function from file */
                int nfine;
                double dtfine;
                FILE* stffile = fopen(stfname,"r");
                if (stffile == 0L)
                {
                        printf("ERROR! src(...) - Source wavelet file '%s' cannot be read.\n",stfname);
                        fflush(stdout);
                        exit(-1);
                }
                fscanf(stffile,"%d %lf", &nfine, &dtfine);
                double* stffine = (double*)malloc((nfine+1)*sizeof(double));
                for(int i=0; i<nfine; i++) fscanf(stffile,"%lf", &stffine[i]);
                stffine[nfine] = 0.;
		
		double f_cut = Butterworth_Low_Pass_Filter_Find_Fcut_From_Fmax(_log_level,fmax,_filter_order);
		if (_log_level >= 4) printf("Using f_cut = %.2lfHz for f_max = %.2lf\n",f_cut,fmax);
		Apply_Butterworth_Low_Pass_Filter(_log_level,stffine,stffine,nfine,dtfine,f_cut,_filter_order);

                int imax = (int)((nfine-1)*dtfine/dt) + 1;
                double absmax = -1e37;
                for(int i=0; i<imax; i++)
                {
                        double t = i*dt;
                        int tfine = (int)(t/dtfine);
                        double frac = t/dtfine - tfine;
                        double val = (1.-frac)*stffine[tfine] + frac*stffine[tfine+1];
                        stf[i] = val;
                        double absval = val < 0.0 ? -val : val;
                        if (absval > absmax) absmax = absval;
                }
                *tsrc = imax;

                if (_log_level >= 4) printf("SOURCE TYPE 2 : nfine=%d, dtfine=%e, imax=%d, absmax=%e\n",nfine,dtfine,imax,absmax);

                for(int i=0; i<imax; i++)
                {
                        stf[i] /= absmax;
                        if (_log_level >= 5) printf("stf[%d] = %e\n",i,stf[i]);
                }
        }
}

bool Elastic_Shot::_Range_Intersects(int x0, int x1, int x_lo, int x_hi)
{
        if ((x_lo < x0 && x_hi >= x0) || (x_lo >= x0 && x_lo <= x1))
                return true;
        else
                return false;
}

bool Elastic_Shot::_Receiver_Intersects(Elastic_Interpolation_t interpolation_method, int x0, int x1, int y0, int y1, float recx, float recy)
{
	if (interpolation_method == Point)
	{
		int ix = (int)trunc(recx);
		int iy = (int)trunc(recy);
		return _Range_Intersects(x0,x1,ix,ix) && _Range_Intersects(y0,y1,iy,iy);
	}
	else if (interpolation_method == Trilinear)
	{
		int ix = (int)lrintf(recx);
		int iy = (int)lrintf(recy);
		return _Range_Intersects(x0,x1,ix,ix+1) && _Range_Intersects(y0,y1,iy,iy+1);
	}
	else if (interpolation_method == Sinc)
	{
		int ix = (int)lrintf(recx) + 1;
		int iy = (int)lrintf(recy) + 1;
		return _Range_Intersects(x0,x1,ix-3,ix+4) && _Range_Intersects(y0,y1,iy-3,iy+4);
	}
	else
	{
		return false;
	}
}

void Elastic_Shot::Start_Extract_Receiver_Values_From_Device(
	Elastic_Propagator* prop, 
	Elastic_Pipeline* pipe,
	int device_id, 
	int* block_offsets, 
	int* timesteps, 
	int* num_rx,
	int* flags,
	int num_blocks, 
	float** d_rxloc_block, 
	float* d_rxres,
	float* h_rxres
	)
{
	// copy rxloc to device
	int max_block_offset = block_offsets[num_blocks-1];
	if (max_block_offset >= 0)
	{
		float* src = _h_pinned_rcv_loc[pipe->Get_ID()] + ((float*)(_h_rcv_binned[pipe->Get_ID()][max_block_offset]) - (float*)(_h_rcv_binned[pipe->Get_ID()][0]));
		int length = _Comp_RxLoc_Length(_h_rcv_binned[pipe->Get_ID()][max_block_offset]);
		if (length > 0) gpuErrchk( cudaMemcpyAsync(d_rxloc_block[num_blocks-1], src, length, cudaMemcpyHostToDevice, prop->Get_Receiver_Stream(device_id)) );
	}
}

void Elastic_Shot::Extract_Receiver_Values_From_Device(
	Elastic_Propagator* prop, 
	Elastic_Pipeline* pipe,
	int device_id, 
	int* block_offsets, 
	int* timesteps, 
	int* num_rx,
	int* flags,
	int num_blocks, 
	float** d_rxloc_block, 
	float* d_rxres,
	float* h_rxres
	)
{
	// copy rxloc to device
	int max_block_offset = block_offsets[num_blocks-1];
	if (max_block_offset >= 0)
	{
		/*
		float* src = _h_pinned_rcv_loc[pipe->Get_ID()] + ((float*)(_h_rcv_binned[pipe->Get_ID()][max_block_offset]) - (float*)(_h_rcv_binned[pipe->Get_ID()][0]));
		int length = _Comp_RxLoc_Length(_h_rcv_binned[pipe->Get_ID()][max_block_offset]);
		gpuErrchk( cudaMemcpyAsync(d_rxloc_block[num_blocks-1], src, length, cudaMemcpyHostToDevice, prop->Get_Receiver_Stream(device_id)) );
		*/

		// extract
		cudaSetDevice(device_id);

		int _arg_is_pv[16];
		float* _arg_cmp[16];
		int _arg_x0[16];
		int _arg_y0[16];
		int _arg_nx[16];
		int _arg_ny[16];
		int _arg_nz[16];
		float* _arg_rxloc[16];
		int _arg_nn[16];
		float* _arg_rxres[16];
		int _nk = 0;

		int rxres_offset = 0;
		for (int iBlk = 0;  iBlk < num_blocks;  ++iBlk)
		{
			num_rx[iBlk] = 0;
			flags[iBlk] = 0;
			int curr_timestep = timesteps[iBlk];
			if (curr_timestep > 0)
			{
				bool stop = false;
				int curr_blk_offset = block_offsets[iBlk];
				for (int ibuf = 0;  !stop && ibuf < pipe->Get_Number_Of_Buffers();  ++ibuf)
				{
					Elastic_Buffer* buffer = pipe->Get_Buffer(ibuf);
					if (buffer->Get_Device_ID() == device_id && buffer->Is_Compute() && !buffer->Is_Partial_Compute())
					{
						int rx_block_offset = buffer->Get_Block_Offset(1,0);
						if (rx_block_offset == curr_blk_offset)
						{
							float* rxloc = d_rxloc_block[iBlk];
							float* rxres = d_rxres + rxres_offset;

							int allowed_flags = buffer->Is_Particle_Velocity() ? 14 : 1;

							// count number of receivers that can be extracted from this buffer
							// this means matching buffer type (particle velocity or stresses)
							// with file selection flags
							int nn = _Comp_RxRes_Length(_h_rcv_binned[pipe->Get_ID()][curr_blk_offset], allowed_flags) / sizeof(float);
							if (nn > 0)
							{
								num_rx[iBlk] = _Comp_RxRes_Length(_h_rcv_binned[pipe->Get_ID()][curr_blk_offset], 15) / sizeof(float);
								flags[iBlk] = allowed_flags;

								// clip search region to output region of pipe
								// otherwise, we will get two contributions for cells on the fringe

								int x0 = rx_block_offset * prop->Get_Block_Size_X();
								int y0 = pipe->Get_Y0();
								Elastic_Modeling_Job* job = prop->Get_Job();
								int nx = prop->Get_Block_Size_X();
								int ny = pipe->Get_Width();
								int nz = buffer->Get_Z1() - buffer->Get_Z0() + 1;

								int one_wf_size_f = 4 * nz;
								int one_y_size_f = one_wf_size_f * 6;
								float* cmp_block = (float*)(buffer->Get_Block_By_Offset(rx_block_offset,0)) + (pipe->Get_Y0() - buffer->Get_Y0()) * one_y_size_f;

								_arg_is_pv[_nk] = buffer->Is_Particle_Velocity() ? 1 : 0;
								_arg_cmp[_nk] = (float*)cmp_block;
								_arg_x0[_nk] = x0;
								_arg_y0[_nk] = y0;
								_arg_nx[_nk] = nx;
								_arg_ny[_nk] = ny;
								_arg_nz[_nk] = nz;
								_arg_rxloc[_nk] = rxloc;
								_arg_nn[_nk] = nn;
								_arg_rxres[_nk] = rxres;
								++_nk;

								if (_nk == 16)
								{
									Host_Extract_Receiver_Values(
										prop->Get_Receiver_Stream(device_id),
										_arg_is_pv,_arg_cmp,_arg_x0,_arg_y0,_arg_nx,_arg_ny,_arg_nz,
										-_job->Get_Propagation_Z0(),
										_arg_rxloc,_arg_nn,_arg_rxres,
										_nk
										);
									_nk = 0;
								}

								// prepare for next iteration
								rxres_offset += nn;
							}

							stop = true;
						}
					}
				}
			}
		}
		if (_nk > 0)
		{
			Host_Extract_Receiver_Values(
					prop->Get_Receiver_Stream(device_id),
					_arg_is_pv,_arg_cmp,_arg_x0,_arg_y0,_arg_nx,_arg_ny,_arg_nz,
					-_job->Get_Propagation_Z0(),
					_arg_rxloc,_arg_nn,_arg_rxres,
					_nk
					);
			_nk = 0;
		}

		// copy rxres from device
		if (rxres_offset > 0)
		{
			//printf("h_rxres = %p, d_rxres = %p, rxres_offset = %d\n",h_rxres,d_rxres,rxres_offset);
			gpuErrchk( cudaMemcpyAsync(h_rxres, d_rxres, rxres_offset*sizeof(float), cudaMemcpyDeviceToHost, prop->Get_Receiver_Stream(device_id)) );
					//printf("RxRes :: Copied %d bytes from device %d\n",rxres_offset*sizeof(float),device_id);
		}
	}
}

void Elastic_Shot::DEMUX_Receiver_Values(
	Elastic_Pipeline* pipe,
	int* block_offsets, 
	int* timesteps, 
	int* num_rx,
	int* flags,
	int num_blocks, 
	float* h_rxres
	)
{
	//printf("Elastic_Shot::DEMUX_Receiver_Values - start\n");
	float** rxres = new float*[num_blocks];
	for (int iBlk = 0, rxres_offset = 0;  iBlk < num_blocks;  ++iBlk) 
	{
		int curr_timestep = timesteps[iBlk];
		if (curr_timestep > 0)
		{
			rxres[iBlk] = h_rxres + rxres_offset;

			// advance rxres pointer for next block
			int allowed_flags = flags[iBlk];
			int curr_blk_offset = block_offsets[iBlk];
			int* trcflag = _h_rcv_trcflag[pipe->Get_ID()][curr_blk_offset];

			int nn = 0;
			for (int ii = 0;  ii < num_rx[iBlk];  ++ii)
				if (trcflag[ii] & allowed_flags)
					++nn;
			rxres_offset += nn;
		}
		else
		{
			rxres[iBlk] = 0L;
		}
	}

	/*
	int* jj_max = new int[num_blocks];
	int** jj2iTrc = new int*[num_blocks];
#pragma omp parallel for
	for (int iBlk = 0;  iBlk < num_blocks;  ++iBlk)
	{
		jj_max[iBlk] = 0;
		jj2iTrc[iBlk] = 0L;
		if (num_rx[iBlk] > 0)
		{
			jj2iTrc[iBlk] = new int[num_rx[iBlk]];
			int curr_timestep = timesteps[iBlk];
			if (curr_timestep > 0)
			{
				int curr_blk_offset = block_offsets[iBlk];
				int allowed_flags = flags[iBlk];
				int* trcidx = _h_rcv_trcidx[pipe->Get_ID()][curr_blk_offset];
				int* trcflag = _h_rcv_trcflag[pipe->Get_ID()][curr_blk_offset];
				int jj = 0;
				for (int ii = 0;  ii < num_rx[iBlk];  ++ii)
				{
					if (trcflag[ii] & allowed_flags)
					{
						jj2iTrc[iBlk][jj] = trcidx[ii];
						++jj;
					}
				}
				jj_max[iBlk] = jj;
			}
		}
	}
	for (int iBlk = 0;  iBlk < num_blocks;  ++iBlk)
        {
		if (jj_max[iBlk] > 0)
		{
#pragma omp parallel for
			for (int jj = 0;  jj < jj_max[iBlk];  ++jj)
			{
				int iTrc = jj2iTrc[iBlk][jj];
				int iFile = _h_trace_iFile[iTrc];
				int nsamp_out = _h_trace_nsamp_out[iFile];
				if (_h_trace_idx_out[iTrc] < nsamp_out)
				{
					int nsamp_in = _h_trace_nsamp_in[iFile];
					int idx_in = timesteps[iBlk] - _h_trace_idx_in[iTrc];
					if (idx_in >= 0 && idx_in < nsamp_in)
					{
						_h_trace_touched[iTrc] = true;
						_h_trace_in[iTrc][idx_in] += rxres[iBlk][jj];
						if (idx_in >= _h_trace_idx_in_nn[iTrc]) _h_trace_idx_in_nn[iTrc] = idx_in + 1;
					}
				}
			}
			delete [] jj2iTrc[iBlk];
		}
	}
	delete [] jj2iTrc;
	delete [] jj_max;
	*/

#pragma omp parallel for
	for (int iBlk = 0;  iBlk < num_blocks;  ++iBlk)
	{
		int curr_timestep = timesteps[iBlk];
		if (curr_timestep > 0)
		{
			int curr_blk_offset = block_offsets[iBlk];
			int allowed_flags = flags[iBlk];
			if (num_rx[iBlk] > 0)
			{
				int* trcidx = _h_rcv_trcidx[pipe->Get_ID()][curr_blk_offset];
				int* trcflag = _h_rcv_trcflag[pipe->Get_ID()][curr_blk_offset];
				for (int ii = 0, jj = 0;  ii < num_rx[iBlk];  ++ii)
				{
					if (trcflag[ii] & allowed_flags)
					{
						int iTrc = trcidx[ii];
						int iFile = _h_trace_iFile[iTrc];
						int nsamp_out = _h_trace_nsamp_out[iFile];
						if (_h_trace_idx_out[iTrc] < nsamp_out)
						{
							int nsamp_in = _h_trace_nsamp_in[iFile];
							int idx_in = curr_timestep - _h_trace_idx_in[iTrc];
							if (idx_in >= 0 && idx_in < nsamp_in)
							{
								_h_trace_touched[iTrc] = true;
								_h_trace_in[iTrc][idx_in] += rxres[iBlk][jj];
								if (idx_in >= _h_trace_idx_in_nn[iTrc]) _h_trace_idx_in_nn[iTrc] = idx_in + 1;
#ifdef RESAMPLE_DEBUG
								if (iTrc == 827)
								{
									printf("iTrc=%d, iPipe=%d, iBlk=%d, flags=%d, curr_timestep=%d, curr_blk_offset=%d, _h_trace_idx_out=%d, _h_trace_idx_in=%d, idx_in=%d, _h_trace_idx_in_nn=%d, _h_trace_in=%f\n",iTrc,pipe->Get_ID(),iBlk,allowed_flags,curr_timestep,curr_blk_offset,_h_trace_idx_out[iTrc],_h_trace_idx_in[iTrc],idx_in,_h_trace_idx_in_nn[iTrc],_h_trace_in[iTrc][idx_in]);
								}
#endif
							}
						}
						++jj;
					}
				}
			}
		}
	}

	delete [] rxres;
	//printf("Elastic_Shot::DEMUX_Receiver_Values - end\n");
}

void Elastic_Shot::Resample_Receiver_Traces()
{
	int num_threads = 0;
#pragma omp parallel
	{
		num_threads = omp_get_num_threads();
	}
	int trc_per_thread = (_num_traces + num_threads - 1) / num_threads;

#pragma omp parallel for
	for (int iThr = 0;  iThr < num_threads;  ++iThr)
	{
		int min_iTrc = iThr * trc_per_thread;
		int max_iTrc = min_iTrc + trc_per_thread;
		if (max_iTrc > _num_traces) max_iTrc = _num_traces;
		for (int iTrc = min_iTrc;  iTrc < max_iTrc;  ++iTrc)
		{
			if (_h_trace_touched[iTrc])
			{
				_h_trace_touched[iTrc] = false;
				int iFile = _h_trace_iFile[iTrc];
				int nsamp_out = _h_trace_nsamp_out[iFile];
				bool ding_dong = false, jabba_dabba_doo = false;
				do
				{
					jabba_dabba_doo = false;
					int idx_out = _h_trace_idx_out[iTrc];
					if (idx_out < nsamp_out)
					{
						float* coeffs = _h_trace_out_sinc_coeffs[iFile][idx_out];
						int idxM = _h_trace_out_idxM[iFile][idx_out] - _h_trace_idx_in[iTrc];
						if (idxM >= 0)
						{
							int idx0 = idxM - 3;
							int idx1 = idx0 + 7;
							if (idx1 < _h_trace_idx_in_nn[iTrc] - 12)
							{
								jabba_dabba_doo = true;
								ding_dong = true;
								float val_out = 0.0f;
								if (idx0 < 0)
								{
#ifdef RESAMPLE_DEBUG
									if (iTrc == 827) printf("trace no %d - idxM = %d, idx0 = %d - linear\n",iTrc,idxM,idx0);
#endif
									// linear interpolation
									val_out = _h_trace_in[iTrc][idxM+1] * coeffs[4] + _h_trace_in[iTrc][idxM] * coeffs[3];
								}
								else
								{
#ifdef RESAMPLE_DEBUG
									if (iTrc == 827)
									{
										printf("trace no %d - idxM = %d, idx0 = %d - sinc(%f",iTrc,idxM,idx0,coeffs[0]);
										for (int i = -2;  i <= 3;  ++i) printf(",%f",coeffs[i+3]);
										printf(",%f)\n",coeffs[7]);
									}
#endif
									for (int i = -3;  i <= 4;  ++i)
									{
										val_out += _h_trace_in[iTrc][idxM+i] * coeffs[i+3];
									}
								}	
								_h_trace_out[iTrc][idx_out] = val_out;
								++_h_trace_idx_out[iTrc];
							}
						}
					}
				} while (jabba_dabba_doo);
				if (ding_dong && _h_trace_idx_out[iTrc] < nsamp_out)
				{
					// processed at least one output sample, so advance input buffer
					int idx0 = _h_trace_out_idxM[iFile][_h_trace_idx_out[iTrc]] - 3 - _h_trace_idx_in[iTrc];
					if (idx0 > 0)
					{
#ifdef RESAMPLE_DEBUG
						if (iTrc == 827) printf("...advancing _h_trace_idx_in[%d] from %d to ",iTrc,_h_trace_idx_in[iTrc]);
#endif
						for (int i = 0;  i < _h_trace_idx_in_nn[iTrc] - idx0;  ++i)
						{
							_h_trace_in[iTrc][i] = _h_trace_in[iTrc][i+idx0];
						}
						for (int i = _h_trace_idx_in_nn[iTrc] - idx0;  i < _h_trace_idx_in_nn[iTrc];  ++i)
						{
							_h_trace_in[iTrc][i] = 0.0f;
						}
						_h_trace_idx_in_nn[iTrc] -= idx0;
						_h_trace_idx_in[iTrc] += idx0;
#ifdef RESAMPLE_DEBUG
						if (iTrc == 827) printf("%d\n",_h_trace_idx_in[iTrc]);
#endif
					}
				}
			}
		}
	}
}

void Elastic_Shot::Write_SEGY_Files()
{
	printf("Elastic_Shot::Write_SEGY_Files\n");
	Global_Coordinate_System* gcs = _job->Get_Voxet()->Get_Global_Coordinate_System();
	for (int iFile = 0;  iFile < _num_segy_files;  ++iFile)
	{
		for (int iSel = 0;  iSel < 4;  ++iSel)
		{
			int flag = 1 << iSel;
			
			// count number of traces belonging to this file with this flag
			int count = 0;
			for (int iTrc = 0;  iTrc < _num_traces;  ++iTrc)
				if (_h_trace_iFile[iTrc] == iFile && _h_trace_flag[iTrc] == flag)
					++count;
			char filename[4096];
			printf("%s has %d traces.\n",_segy_files[iFile]->Get_Full_Path(filename,flag),count);

			// gather output traces and write to file
			if (count > 0)
			{
				double srcx, srcy, srcz;
				gcs->Convert_Transposed_Fractional_Index_To_Global(_x,_y,_z,srcx,srcy,srcz);
				float** traces = new float*[count];
				double* recx = new double[count];
				double* recy = new double[count];
				double* recz = new double[count];
				int *iline = new int[count];
				int *xline = new int[count];
				int *trcens = new int[count];
				for (int iTrc = 0, ii = 0;  iTrc < _num_traces;  ++iTrc)
				{
	                                if (_h_trace_iFile[iTrc] == iFile && _h_trace_flag[iTrc] == flag)
					{
						traces[ii] = _h_trace_out[iTrc];
						gcs->Convert_Transposed_Fractional_Index_To_Global(_h_trace_rcv_x[iTrc],_h_trace_rcv_y[iTrc],_h_trace_rcv_z[iTrc],recx[ii],recy[ii],recz[ii]);
						iline[ii] = _h_trace_iline[iTrc];
						xline[ii] = _h_trace_xline[iTrc];
						trcens[ii] = _h_trace_trcens[iTrc];
						++ii;
					}
				}
				int nsamp = _h_trace_nsamp_out[iFile];
				_segy_files[iFile]->Write_SEGY_File(traces,srcx,srcy,srcz,recx,recy,recz,iline,xline,trcens,count,nsamp,flag);
				delete [] traces;
				delete [] recx;
				delete [] recy;
				delete [] recz;
				delete [] iline;
				delete [] xline;
				delete [] trcens;
			}
		}
	}
}

// compute the length of this rxloc buffer in bytes
unsigned long Elastic_Shot::_Comp_RxLoc_Length(float* rxloc)
{
	if (rxloc != 0L)
	{
		int* iloc = (int*)rxloc;
		int num_files = iloc[0];
		unsigned long length = (unsigned long)(1+2*num_files);
		for (int iFile = 0;  iFile < num_files;  ++iFile)
		{
			int nn = iloc[1+2*iFile];
			length += 3 * nn;
		}
		return length * sizeof(float);
	}
	else
	{
		return 0L;
	}
}

// compute the maximum length needed for this rxres buffer in bytes.
// receivers are filtered on flags
unsigned long Elastic_Shot::_Comp_RxRes_Length(float* rxloc, int flags)
{
	if (rxloc != 0L)
	{
		int* iloc = (int*)rxloc;
		int num_files = iloc[0];
		unsigned long length = 0;
		for (int iFile = 0;  iFile < num_files;  ++iFile)
		{
			int curr_nn = iloc[1+2*iFile];
			int curr_flags = iloc[2+2*iFile] & flags;
			int num_wf = 0;
			if (curr_flags & 1) ++num_wf;
			if (curr_flags & 2) ++num_wf;
			if (curr_flags & 4) ++num_wf;
			if (curr_flags & 8) ++num_wf;
			length += num_wf * curr_nn;
		}
		return length * sizeof(float);
	}
	else
	{
		return 0L;
	}
}

// Calculate the maximum size of the buffer on the GPU side used to hold rx locations
// for all full compute buffers
void Elastic_Shot::Calculate_RX_Locations_And_Results_Size(Elastic_Propagator* prop, int device_id, int& rxloc_size, int& rxres_size, int& num_full_compute)
{
	if (_num_pipes <= 0) _Create_Receiver_Transfer_Buffers(prop);
	
	// which pipeline does device belong to?
	Elastic_Pipeline* pipe = 0L;
	for (int iPipe = 0;  pipe == 0L && iPipe < prop->Get_Number_Of_Pipelines();  ++iPipe)
	{
		Elastic_Pipeline* curr_pipe = prop->Get_Pipeline(iPipe);
		for (int iDev = 0;  pipe == 0L && iDev < curr_pipe->Get_Device_Count();  ++iDev)
		{
			if (curr_pipe->Get_All_Device_IDs()[iDev] == device_id)
			{
				pipe = curr_pipe;
			}	
		}
	}

	// determine how many full compute buffers this device has
	num_full_compute = 0;
	for (int ibuf = 0;  ibuf < pipe->Get_Number_Of_Buffers();  ++ibuf)
	{
		Elastic_Buffer* buffer = pipe->Get_Buffer(ibuf);
		if (buffer->Get_Device_ID() == device_id && buffer->Is_Compute() && !buffer->Is_Partial_Compute()) ++ num_full_compute;
	}
	
	// loop over all the block offsets and calculate maximum values for rxloc and rxres buffers
	rxloc_size = 0;
	rxres_size = 0;
	for (int iBlk = 0;  iBlk < _nBlks;  ++iBlk)
	{
		// rxloc buffer is allocated individually for each block offset
		float* rxloc = _h_rcv_binned[pipe->Get_ID()][iBlk];
		int curr_rxloc_size = _Comp_RxLoc_Length(rxloc) * num_full_compute;
		if (curr_rxloc_size > rxloc_size) rxloc_size = curr_rxloc_size;

		// rxres buffer is allocated in one chunk for all block offsets
		int curr_rxres_size = 0;
		for (int iDelta = 0;  iDelta < num_full_compute;  ++iDelta)
		{
			int curr_blk = (iBlk + iDelta) % _nBlks;
			float* rxloc2 = _h_rcv_binned[pipe->Get_ID()][curr_blk];
			int curr_rxres_pressure_size = _Comp_RxRes_Length(rxloc2, 1);
			int curr_rxres_stress_size = _Comp_RxRes_Length(rxloc2, 14);
			curr_rxres_size += curr_rxres_pressure_size > curr_rxres_stress_size ? curr_rxres_pressure_size : curr_rxres_stress_size;
		}
		if (curr_rxres_size > rxres_size) rxres_size = curr_rxres_size;
	}

	//printf("Elastic_Shot::Calculate_RX_Locations_And_Results_Size - num_full_compute=%d, rxloc_size=%d, rxres_size=%d\n",num_full_compute,rxloc_size,rxres_size);
}

double Elastic_Shot::Find_Deepest_Receiver()
{
	double deepest = 0.0;
	for (int iFile = 0;  iFile < _num_segy_files;  ++iFile)
        {
                double *rcv_x, *rcv_y, *rcv_z;
                int num_rx = _segy_files[iFile]->Compute_Receiver_Locations(rcv_x, rcv_y, rcv_z);
		for (int iRx = 0;  iRx < num_rx;  ++iRx)
		{
			if (rcv_z[iRx] > deepest) deepest = rcv_z[iRx];
		}
		delete [] rcv_x;
		delete [] rcv_y;
		delete [] rcv_z;
	}
	return deepest;
}

void Elastic_Shot::_Create_Receiver_Transfer_Buffers(Elastic_Propagator* prop)
{
        //
        // create binned receiver locations for all files
        //

	double dti = prop->Get_Internal_Sample_Rate();
        _totSteps = prop->Get_Total_Number_Of_Timesteps();
        _nBlks = prop->Get_Number_Of_Blocks();
        _num_pipes = prop->Get_Number_Of_Pipelines();

        _h_rcv_loc = new float*[_num_pipes];
        _h_rcv_binned = new float**[_num_pipes];
        _h_rcv_trcidx = new int**[_num_pipes];
        _h_rcv_trcflag = new int**[_num_pipes];
        _h_rcv_loc_size_f = new int[_num_pipes];
        for (int iPipe = 0;  iPipe < _num_pipes;  ++iPipe)
        {
                _h_rcv_loc[iPipe] = 0L;
                _h_rcv_binned[iPipe] = new float*[_nBlks];
                _h_rcv_trcidx[iPipe] = new int*[_nBlks];
                _h_rcv_trcflag[iPipe] = new int*[_nBlks];
                _h_rcv_loc_size_f[iPipe] = 0;
                for (int iBlk = 0;  iBlk < _nBlks;  ++iBlk)
                {
                        _h_rcv_binned[iPipe][iBlk] = 0L;
                        _h_rcv_trcidx[iPipe][iBlk] = 0L;
                        _h_rcv_trcflag[iPipe][iBlk] = 0L;
                }
        }

	// determine total number of traces in all files
	_num_traces = 0;
	for (int iFile = 0;  iFile < _num_segy_files;  ++iFile)
	{
		double *rcv_x, *rcv_y, *rcv_z;
		int num_rx = _segy_files[iFile]->Compute_Receiver_Locations(rcv_x, rcv_y, rcv_z);
		if (num_rx > 0)
		{
			int num_wf = 0;
			int flags = _segy_files[iFile]->Get_Selection_Flags();
			for (int iSel = 0;  iSel < 4;  ++iSel) if (flags & (1 << iSel)) ++num_wf;
			_num_traces += num_wf * num_rx;
		}
		delete [] rcv_x;
		delete [] rcv_y;
		delete [] rcv_z;
	}

	_h_trace_rcv_x = new double[_num_traces];
	_h_trace_rcv_y = new double[_num_traces];
	_h_trace_rcv_z = new double[_num_traces];
	_h_trace_iline = new int[_num_traces];
	_h_trace_xline = new int[_num_traces];
	_h_trace_trcens = new int[_num_traces];
	for (int iFile = 0, trace_no_idx = 0;  iFile < _num_segy_files;  ++iFile)
        {
                double *rcv_x, *rcv_y, *rcv_z;
		int *il, *xl, *trce;
                int num_rx = _segy_files[iFile]->Compute_Receiver_Locations(rcv_x, rcv_y, rcv_z, il, xl, trce);
		int flags = _segy_files[iFile]->Get_Selection_Flags();
		for (int iRx = 0;  iRx < num_rx;  ++iRx)
                {
			for (int iSel = 0;  iSel < 4;  ++iSel)
			{
				if (flags & (1 << iSel))
				{
					_h_trace_rcv_x[trace_no_idx] = rcv_x[iRx];
					_h_trace_rcv_y[trace_no_idx] = rcv_y[iRx];
					_h_trace_rcv_z[trace_no_idx] = rcv_z[iRx];
					_h_trace_iline[trace_no_idx] = il[iRx];
					_h_trace_xline[trace_no_idx] = xl[iRx];
					_h_trace_trcens[trace_no_idx] = trce[iRx];
					++trace_no_idx;
				}
			}
		}
		delete [] rcv_x;
		delete [] rcv_y;
		delete [] rcv_z;
		delete [] il;
		delete [] xl;
		delete [] trce;
	}

        for (int iPipe = 0;  iPipe < _num_pipes;  ++iPipe)
        {
                Elastic_Pipeline* pipe = prop->Get_Pipeline(iPipe);
                int y0 = pipe->Get_Y0() + _job->Get_Propagation_Y0();
                int y1 = pipe->Get_Y1() + _job->Get_Propagation_Y0();

                // ..count total number of receiver locations in this pipe
                int tot_rx = 0, tot_traces = 0;
                for (int iFile = 0;  iFile < _num_segy_files;  ++iFile)
                {
                        double *rcv_x, *rcv_y, *rcv_z;
                        int num_rx = _segy_files[iFile]->Compute_Receiver_Locations(rcv_x, rcv_y, rcv_z);
			int flags = _segy_files[iFile]->Get_Selection_Flags();
			Elastic_Interpolation_t interpolation_method = _segy_files[iFile]->Get_Interpolation_Method();
			int num_wf = 0;
			if (flags & 1) ++num_wf;
			if (flags & 2) ++num_wf;
			if (flags & 4) ++num_wf;
			if (flags & 8) ++num_wf;
                        if (num_rx > 0)
                        {
                                for (int iBlk = 0;  iBlk < _nBlks;  ++iBlk)
                                {
                                        int x0 = iBlk * prop->Get_Block_Size_X() + _job->Get_Propagation_X0();
                                        int x1 = x0 + prop->Get_Block_Size_X() - 1;
                                        for (int iRx = 0;  iRx < num_rx;  ++iRx)
                                        {
                                                if (_Receiver_Intersects(interpolation_method,x0,x1,y0,y1,rcv_x[iRx],rcv_y[iRx]))
                                                {
                                                        ++tot_rx;
							tot_traces+=num_wf;
                                                }
                                        }
                                }
                                delete [] rcv_x;
                                delete [] rcv_y;
                                delete [] rcv_z;
                        }
                }

                // ..bin receiver locations
                _h_rcv_loc_size_f[iPipe] = _nBlks * (1 + 2 * _num_segy_files) + 3 * tot_rx;
		_h_rcv_loc[iPipe] = new float[_h_rcv_loc_size_f[iPipe]];
                float* floc = _h_rcv_loc[iPipe];
                int* iloc = (int*)floc;
                int* trcidx = new int[tot_traces];
		int* trcflag = new int[tot_traces];
                for (int iBlk = 0;  iBlk < _nBlks;  ++iBlk)
                {
                        int x0 = iBlk * prop->Get_Block_Size_X() + _job->Get_Propagation_X0();
                        int x1 = x0 + prop->Get_Block_Size_X() - 1;

                        // ..count receivers in this block
                        _h_rcv_binned[iPipe][iBlk] = floc;
                        iloc[0] = _num_segy_files;
                        _h_rcv_trcidx[iPipe][iBlk] = trcidx;
			_h_rcv_trcflag[iPipe][iBlk] = trcflag;
                        int blk_offset = 0;
			int trace_blk_offset = 0;
			int rxres_length_f = 0;
			int trace_no_idx = 0;
                        for (int iFile = 0, trace_no = 0;  iFile < _num_segy_files;  ++iFile)
                        {
				int flags = _segy_files[iFile]->Get_Selection_Flags();
				Elastic_Interpolation_t interpolation_method = _segy_files[iFile]->Get_Interpolation_Method();
				bool receiver_ghost_enabled = !_job->Freesurface_Enabled() & _job->Receiver_Ghost_Enabled();
                                double *rcv_x, *rcv_y, *rcv_z;
                                int num_rx = _segy_files[iFile]->Compute_Receiver_Locations(rcv_x, rcv_y, rcv_z);
				int pr_rx = (flags & 1) ? num_rx : 0;
                                iloc[1+2*iFile] = 0;
                                // 0001 -> P
                                // 0010 -> Vx
                                // 0100 -> Vy
                                // 1000 -> Vz
				iloc[2+2*iFile] = flags | ((int)interpolation_method << 16) | ((receiver_ghost_enabled ? 1 : 0) << 31);
                                if (num_rx > 0)
                                {
                                        int file_offset = 0;
                                        for (int iRx = 0;  iRx < num_rx;  ++iRx)
                                        {
                                                if (_Receiver_Intersects(interpolation_method,x0,x1,y0,y1,rcv_x[iRx],rcv_y[iRx]))
                                                {
                                                        int off = 1 + 2 * _num_segy_files + 3 * (blk_offset + file_offset);
                                                        floc[off  ] = rcv_x[iRx] - _job->Get_Propagation_X0();
                                                        floc[off+1] = rcv_y[iRx] - _job->Get_Propagation_Y0();
                                                        floc[off+2] = rcv_z[iRx] - _job->Get_Propagation_Z0();
							++file_offset;
						
							for (int iSel = 0, trace_no_off = 0;  iSel < 4;  ++iSel)
							{
								if (flags & (1 << iSel))
								{
									trcidx[trace_no_idx] = trace_no + trace_no_off;
									++trace_no_off;
									trcflag[trace_no_idx] = 1 << iSel;
									++trace_no_idx;
								}
							}
                                                }
						for (int iSel = 0;  iSel < 4;  ++iSel) if (flags & (1 << iSel)) ++trace_no;
                                        }
                                        iloc[1+2*iFile] = file_offset;
                                        blk_offset += file_offset;
                                }
                                delete [] rcv_z;
                                delete [] rcv_y;
                                delete [] rcv_x;
                        }
			//printf("iBlk = %d, trace_no_idx = %d\n",iBlk,trace_no_idx);
                        int glob_offset_inc = 1 + 2 * _num_segy_files + 3 * blk_offset;
                        floc = floc + glob_offset_inc;
                        iloc = iloc + glob_offset_inc;
                        trcidx = trcidx + trace_no_idx;
			trcflag = trcflag + trace_no_idx;
                }
        }

#ifdef DEBUG_TMJ
	FILE* fp = fopen("/panfs07/esdrd/tjhc/ELA_on_GPU/rx_locs.txt", "w");
	if (fp != 0L)
	{
		fprintf(fp,"\n");
		for (int iPipe = 0;  iPipe < _num_pipes;  ++iPipe)
		{
			Elastic_Pipeline* pipe = prop->Get_Pipeline(iPipe);
			int y0 = pipe->Get_Y0();
			int y1 = pipe->Get_Y1();
			fprintf(fp,"P I P E   %d\n\n",iPipe);

			for (int iBlk = 0;  iBlk < _nBlks;  ++iBlk)
			{
				int x0 = iBlk * prop->Get_Block_Size_X();
				int x1 = x0 + prop->Get_Block_Size_X() - 1;
				fprintf(fp,"   B L O C K   %d\n\n",iBlk);

				float* floc = _h_rcv_binned[iPipe][iBlk];
				int* iloc = (int*)floc;
				int* trcidx = _h_rcv_trcidx[iPipe][iBlk];
				int* trcflag = _h_rcv_trcflag[iPipe][iBlk];

				int num_files = iloc[0];
				int bulk_offset = 1 + 2 * num_files;

				for (int iFile = 0, trace_no_idx = 0;  iFile < num_files;  ++iFile)
				{
					int num_rx = iloc[1+2*iFile];
					int flags = iloc[1+2*iFile+1];
					fprintf(fp,"      F I L E   %d - %d receivers, %s%s%s%s selected\n\n",iFile,num_rx,(flags&1)?"P ":"",(flags&2)?"Vx ":"",(flags&4)?"Vy ":"",(flags&8)?"Vz ":"");
					for (int iRx = 0;  iRx < num_rx;  ++iRx)
					{
						float x = floc[bulk_offset+3*iRx  ];
						float y = floc[bulk_offset+3*iRx+1];
						float z = floc[bulk_offset+3*iRx+2];
						fprintf(fp,"         %.2f %.2f %.2f :: trace",x,y,z);
						for (int iSel = 0;  iSel < 4;  ++iSel)
						{
							if (flags & (1 << iSel))
							{
								char* FlagStr = "";
								if (trcflag[trace_no_idx] & 1) FlagStr = "P";
								else if (trcflag[trace_no_idx] & 2) FlagStr = "Vx";
								else if (trcflag[trace_no_idx] & 4) FlagStr = "Vy";
								else if (trcflag[trace_no_idx] & 8) FlagStr = "Vz";
								fprintf(fp," %d(%s)",trcidx[trace_no_idx],FlagStr);
								++trace_no_idx;
							}
						}
						fprintf(fp,"\n");
					}
					bulk_offset += 3 * num_rx;
					fprintf(fp,"\n");
				}
			}
		}
		fclose(fp);
	}
	else
	{
		printf("Warning! Unable to write to %s.\n","/panfs07/esdrd/tjhc/ELA_on_GPU/rx_locs.txt");
	}
#endif
}

void Elastic_Shot::Free_Trace_Resample_Buffers()
{
	if (_num_traces > 0)
	{
		for (int iTrc = 0;  iTrc < _num_traces;  ++iTrc)
		{
			delete [] _h_trace_in[iTrc];
			delete [] _h_trace_out[iTrc];
		}
		delete [] _h_trace_in;
		_h_trace_in = 0L;
		delete [] _h_trace_out;
		_h_trace_out = 0L;
		delete [] _h_trace_flag;
		_h_trace_flag = 0L;
		delete [] _h_trace_rcv_x;
		_h_trace_rcv_x = 0L;
		delete [] _h_trace_rcv_y;
		_h_trace_rcv_y = 0L;
		delete [] _h_trace_rcv_z;
		_h_trace_rcv_z = 0L;
		delete [] _h_trace_iline;
		_h_trace_iline = 0L;
		delete [] _h_trace_xline;
		_h_trace_xline = 0L;
		delete [] _h_trace_trcens;
		_h_trace_trcens = 0L;
		delete [] _h_trace_tshift;
		_h_trace_tshift = 0L;
		delete [] _h_trace_sample_rate;
		_h_trace_sample_rate = 0L;
		delete [] _h_trace_idx_in;
		_h_trace_idx_in = 0L;
		delete [] _h_trace_idx_in_nn;
		_h_trace_idx_in_nn = 0L;
		delete [] _h_trace_idx_out;
		_h_trace_idx_out = 0L;
		delete [] _h_trace_iFile;
		_h_trace_iFile = 0L;
		delete [] _h_trace_touched;
		_h_trace_touched = 0L;
		_num_traces = 0;
	}
	if (_num_pipes > 0)
	{
		for (int iPipe = 0;  iPipe < _num_pipes;  ++iPipe)
		{
			if (_h_rcv_loc[iPipe] != 0L) delete [] _h_rcv_loc[iPipe];
			if (_h_rcv_binned[iPipe] != 0L) delete [] _h_rcv_binned[iPipe];
			delete [] _h_rcv_trcidx[iPipe][0];
			delete [] _h_rcv_trcidx[iPipe];
			delete [] _h_rcv_trcflag[iPipe][0];
			delete [] _h_rcv_trcflag[iPipe];
		}
		delete [] _h_rcv_loc;
		delete [] _h_rcv_binned;
		delete [] _h_rcv_trcidx;
		delete [] _h_rcv_trcflag;
		_h_rcv_loc = 0L;
		_h_rcv_binned = 0L;
		_h_rcv_trcidx = 0L;
		_h_rcv_trcflag = 0L;
		_num_pipes = 0;
	}
}

void Elastic_Shot::Create_Trace_Resample_Buffers(Elastic_Propagator* prop)
{
	double dti = prop->Get_Internal_Sample_Rate();

	// create structures needed for temporal resampling
	_h_trace_tshift = new double[_num_traces];
	_h_trace_sample_rate = new double[_num_traces];
	_h_trace_in = new float*[_num_traces];
	_h_trace_out = new float*[_num_traces];
	_h_trace_flag = new int[_num_traces];
	_h_trace_idx_in = new int[_num_traces];
	_h_trace_idx_in_nn = new int[_num_traces];
	_h_trace_idx_out = new int[_num_traces];
	_h_trace_iFile = new int[_num_traces];
	_h_trace_touched = new bool[_num_traces];
	_h_trace_nsamp_in = new int[_num_segy_files];
	_h_trace_nsamp_out = new int[_num_segy_files];
	_h_trace_out_idxM = new int*[_num_segy_files];
	_h_trace_out_sinc_coeffs = new float**[_num_segy_files];
	int nsamp_in = _totSteps + 12;
	for (int iFile = 0, iTrc = 0;  iFile < _num_segy_files;  ++iFile)
        {
		double tshift = _segy_files[iFile]->Get_Timeshift();
		double sample_rate = _segy_files[iFile]->Get_Sample_Rate();
		int nsamp_out = (int)round(_segy_files[iFile]->Get_Record_Length() / sample_rate);
		int flags = _segy_files[iFile]->Get_Selection_Flags();

		// create trace info
                double *rcv_x, *rcv_y, *rcv_z;
                int num_rx = _segy_files[iFile]->Compute_Receiver_Locations(rcv_x, rcv_y, rcv_z);
                if (num_rx > 0)
                {
			for (int iRx = 0;  iRx < num_rx;  ++iRx)
			{
				for (int iSel = 0;  iSel < 4;  ++iSel)
				{
					if (flags & (1 << iSel))
					{
						_h_trace_tshift[iTrc] = tshift;
						_h_trace_sample_rate[iTrc] = sample_rate;
						_h_trace_flag[iTrc] = 1 << iSel;
						double t0 = tshift + sample_rate;
						int cnt = 0;
						for (double tcurr = (double)cnt * dti;  tcurr < t0;  tcurr = (double)cnt * dti) ++cnt;
						cnt -= 4;
						if (cnt < 0) cnt = 0;
#ifdef TRACE_DEBUG
						printf("iFile=%d, iRx=%d, tshift=%f, dti=%f, cnt=%d\n",iFile,iRx,tshift,dti,cnt);
#endif
						_h_trace_idx_in[iTrc] = cnt;
						_h_trace_idx_in_nn[iTrc] = 0;
						_h_trace_idx_out[iTrc] = 1;
						_h_trace_in[iTrc] = new float[nsamp_in];
						for (int i = 0;  i < nsamp_in;  ++i) _h_trace_in[iTrc][i] = 0.0f;
						_h_trace_out[iTrc] = new float[nsamp_out];
						for (int i = 0;  i < nsamp_out;  ++i) _h_trace_out[iTrc][i] = 0.0f;
						_h_trace_iFile[iTrc] = iFile;
						++iTrc;
					}
				}
			}
		}
		delete [] rcv_x;
		delete [] rcv_y;
		delete [] rcv_z;

		// create sinc coefficients for each file
		_h_trace_nsamp_in[iFile] = nsamp_in;
		_h_trace_nsamp_out[iFile] = nsamp_out;
		_h_trace_out_idxM[iFile] = new int[nsamp_out];
		_h_trace_out_sinc_coeffs[iFile] = new float*[nsamp_out];
		_h_trace_out_sinc_coeffs[iFile][0] = 0L;
		_h_trace_out_idxM[iFile][0] = 0;
		for (int idx = 1;  idx < nsamp_out;  ++idx)
		{
			double t_out = (double)idx * sample_rate + tshift;
			double i_out = trunc(t_out/dti);
			double remainder = (t_out - i_out * dti) / sample_rate;
		
			_h_trace_out_idxM[iFile][idx] = (int)i_out;
			_h_trace_out_sinc_coeffs[iFile][idx] = new float[8];
			int idx0 = _h_trace_out_idxM[iFile][idx] - 3;
			if (idx0 < 0 || remainder == 0.0)
			{
				// linear interpolation
				_h_trace_out_sinc_coeffs[iFile][idx][0] = 0.0f;
				_h_trace_out_sinc_coeffs[iFile][idx][1] = 0.0f;
				_h_trace_out_sinc_coeffs[iFile][idx][2] = 0.0f;
				_h_trace_out_sinc_coeffs[iFile][idx][3] = 1.0 - remainder;
				_h_trace_out_sinc_coeffs[iFile][idx][4] = remainder;
				_h_trace_out_sinc_coeffs[iFile][idx][5] = 0.0f;
				_h_trace_out_sinc_coeffs[iFile][idx][6] = 0.0f;
				_h_trace_out_sinc_coeffs[iFile][idx][7] = 0.0f;
			}
			else
			{
				for (int i = -3;  i <= 4;  ++i)
				{
					double x = 3.1415926535897932384626433832795 * ((double)i - remainder);
					_h_trace_out_sinc_coeffs[iFile][idx][i+3] = x == 0.0 ? 1.0 : (float)(sin(x)/x);
				}
			}
		}
	}
}

void Elastic_Shot::Allocate_Pinned_Host_Memory(Elastic_Propagator* prop)
{
	if (_num_pipes <= 0) _Create_Receiver_Transfer_Buffers(prop);
	
	_h_pinned_rcv_loc = new float*[_num_pipes];
	for (int iPipe = 0;  iPipe < prop->Get_Number_Of_Pipelines();  ++iPipe)
	{
		// allocate pinned buffer for rx locations 
		gpuErrchk( cudaHostAlloc((void**)&(_h_pinned_rcv_loc[iPipe]), (size_t)_h_rcv_loc_size_f[iPipe] * sizeof(float), cudaHostAllocDefault) );
		// copy from non-pinned buffer.
		memcpy(_h_pinned_rcv_loc[iPipe], _h_rcv_loc[iPipe], _h_rcv_loc_size_f[iPipe] * sizeof(float));
	}
}

void Elastic_Shot::Free_Pinned_Host_Memory(Elastic_Propagator* prop)
{
	if (_h_pinned_rcv_loc != 0l)
	{
		for (int iPipe = 0;  iPipe < prop->Get_Number_Of_Pipelines();  ++iPipe)
		{
			if (_h_pinned_rcv_loc[iPipe] != 0L)
			{
				gpuErrchk( cudaFreeHost(_h_pinned_rcv_loc[iPipe]) );
				_h_pinned_rcv_loc[iPipe] = 0L;
			}
		}
		delete [] _h_pinned_rcv_loc;
		_h_pinned_rcv_loc = 0L;
	}
}

void Elastic_Shot::Add_SEGY_File(Elastic_SEGY_File* segy_file)
{
	if (_segy_files == 0L)
	{
		_segy_files = new Elastic_SEGY_File*[1];
		_segy_files[0] = segy_file;
		_num_segy_files = 1;
	}
	else
	{
		Elastic_SEGY_File** new_files = new Elastic_SEGY_File*[_num_segy_files+1];
		for (int i = 0;  i < _num_segy_files;  ++i) new_files[i] = _segy_files[i];
		new_files[_num_segy_files] = segy_file;
		delete [] _segy_files;
		_segy_files = new_files;
		++_num_segy_files;
	}
}

Elastic_SEGY_File* Elastic_Shot::Get_SEGY_File(int segy_file_idx)
{
	for (int i = 0;  i < _num_segy_files;  ++i)
	{
		if (_segy_files[i]->Get_File_Index() == segy_file_idx)
		{
			return _segy_files[i];
		}
	}
	return 0L;
}

