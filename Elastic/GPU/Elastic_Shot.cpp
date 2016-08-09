#include <cassert>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <xmmintrin.h>
#include <omp.h>
#include <Elastic_Shot.hxx>
#include <Elastic_Propagator.hxx>
#include <Elastic_Modeling_Job.hxx>
#include <Elastic_Buffer.hxx>
#include <Elastic_SEGY_File.hxx>
#include <Elastic_Pipeline.hxx>
#include <Extract_Receiver_Values.hxx>
#include <Voxet.hxx>
#include <Global_Coordinate_System.hxx>
#include <DFT.hxx>
#include <gpuAssert.h>

//#define RESAMPLE_DEBUG 0
//#define DEBUG_TMJ

Elastic_Shot::Elastic_Shot(int log_level, Elastic_Modeling_Job* job, int souidx, double x, double y, double z)
{
	_log_level = log_level;
	_job = job;
	_ordertime = 2;  // hardcode for now.
	_souintrp = Trilinear; // 10/09/14 - Changed default from point to trilinear
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
	_dti = 0.0;
        _h_rcv_loc = 0L;
        _h_pinned_rcv_loc = 0L;
        _h_rcv_binned = 0L;
        _h_rcv_trcidx = 0L;
        _h_rcv_trcflag = 0L;
        _h_rcv_loc_size_f = 0L;
	_num_traces = 0;
	_h_traces_hdr = 0L;
	_h_traces = 0L;
	/*
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
	*/
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
		}
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

void Elastic_Shot::Add_Receiver_Array(int nrec,	double* rec_x, double* rec_y, double* rec_z, int* iline, int* xline, int* trcens, int* rec_ffid, time_t* acqtime, int* usec ) {
	
	for (int iFile = 0;  iFile < _num_segy_files;  ++iFile){
		_segy_files[iFile]->Add_Receiver_Array(nrec,rec_x,rec_y,rec_z, iline, xline, trcens, rec_ffid, acqtime, usec);
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

	Compute_Time_Integrated_Source_Wavelet(_log_level,_stf,_stf_int,_tsrc,dt);
	if (debug_output_source_wavelet)
	{
		for (int iFile = 0;  iFile < _num_segy_files;  ++iFile)
		{
			double srcx, srcy, srcz;
			Global_Coordinate_System* gcs = _job->Get_Voxet()->Get_Global_Coordinate_System();
			gcs->Convert_Transposed_Fractional_Index_To_Global(_x,_y,_z,srcx,srcy,srcz);
			_segy_files[iFile]->Write_Source_Wavelet_To_SEGY_File(_job->Is_Vwxyzt(),_stf,_stf_int,dt,_tsrc,srcx,srcy,srcz,_il,_xl);
		}
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

bool Elastic_Shot::Read_Source_Wavelet_From_File(const char* wavelet_path)
{
	FILE* fp = fopen(wavelet_path, "rb");
	if (fp != 0L)
	{
		fclose(fp);
		_wavelet_path = strdup(wavelet_path);
		_max_freq = -1.0;
		_filter_order = -1;
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
	double fpeak = fmax / 2.77;
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
                int i, imax, ntap;
                double w0, wmax, ts,t,rt, diff;

                wmax = 2.0*M_PI*fmax;
                /* Note:   tsrc = ts/dt = 2/(gam*khmax)*(rt*rt)*(Vmax/Vmin) */

                /* if(type==1) */  /* only one type for now */
                {
                        rt = 3.571625; /* guarantees SourceAmplitude(tmax)   = 0.01*MaxAmplitude
                                           and guarantees SpectralAmplitude(fmax) = 0.01*MaxSpectrum
                                           and: wmax/w0 = rt = w0*ts/2  */
                        w0 = wmax/rt;  /* w0i = 1./w0; */
                        ts = 2.0*rt/w0;  /* total source time */
                        imax = (int)(2.0*ts/dt) + 1;

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
		for(int i=0; i<nfine; i++) stffine[i] = 0.0;
                for(int i=0; i<nfine; i++) fscanf(stffile,"%lf", &stffine[i]);
	
		if (_filter_order > 0)
		{	
			double f_cut = Butterworth_Low_Pass_Filter_Find_Fcut_From_Fmax(_log_level,fmax,_filter_order);
			if (_log_level >= 4) printf("Using f_cut = %.2lfHz for f_max = %.2lf\n",f_cut,fmax);
			Apply_Butterworth_Low_Pass_Filter(_log_level,stffine,stffine,nfine,dtfine,f_cut,_filter_order);
		}
		else
		{
			if (_log_level >= 3) printf("Warning! No filter will be applied to this source wavelet, this may cause dispersion.\n");
		}

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
		int ix = (int)lrintf(recx);
		int iy = (int)lrintf(recy);
		return _Range_Intersects(x0,x1,ix,ix) && _Range_Intersects(y0,y1,iy,iy);
	}
	else if (interpolation_method == Trilinear)
	{
		int ix = (int)truncf(recx);
		int iy = (int)truncf(recy);
		// TMJ, -1,+2 range is because Vx and Vy receivers are staggered horizontally. DON'T CHANGE TO 0,+1
		return _Range_Intersects(x0,x1,ix-1,ix+2) && _Range_Intersects(y0,y1,iy-1,iy+2);
	}
	else if (interpolation_method == Sinc)
	{
		int ix = (int)lrintf(recx);
		int iy = (int)lrintf(recy);
		// TMJ, all receivers can use the same range, we adjust sinc coefficients to account for staggered grids
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
		if (length > 0) 
		{
			gpuErrchk( cudaMemcpyAsync(d_rxloc_block[num_blocks-1], src, length, cudaMemcpyHostToDevice, prop->Get_Receiver_Stream(device_id)) );
			prop->Add_H2D(length);
		}
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
								int nx = prop->Get_Block_Size_X();
								int ny = pipe->Get_Width();
								int nz = buffer->Get_Z1() - buffer->Get_Z0() + 1;

								int one_wf_size_f = nx * nz;
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
			prop->Add_D2H(rxres_offset*sizeof(float));
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

	// TMJ 12/18/2014
	// Blocks cannot be processed in parallel.
	// Two adjacent blocks may contribute to the same trace.
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
					if ((ii&15) == 0)
					{
						_mm_prefetch((char*)(trcflag+ii+64),_MM_HINT_T0);
						_mm_prefetch((char*)(trcidx+ii+64),_MM_HINT_T0);
					}
					if (trcflag[ii] & allowed_flags)
					{
						int iTrc = trcidx[ii];
						Elastic_Resampled_Trace* trace = _h_traces[iTrc];
						trace->Add_To_Trace(curr_timestep, rxres[iBlk][jj]);
						//printf("Add_To_Trace(iTrc=%d, timestep=%d, value=%e\n",iTrc,curr_timestep,rxres[iBlk][jj]);
						++jj;
					}
				}
			}
		}
	}

	delete [] rxres;
	//printf("Elastic_Shot::DEMUX_Receiver_Values - end\n");
}

float Elastic_Shot::Compute_Reciprocal_Scale_Factor(int flag, float srcx, float srcy, float srcz, float recx, float recy, float recz)
{
	int isx = (int)round(srcx) - _job->Get_Propagation_X0();
	int isy = (int)round(srcy) - _job->Get_Propagation_Y0();
	int isz = (int)round(srcz) - _job->Get_Propagation_Z0();
	int irx = (int)round(recx) - _job->Get_Propagation_X0();
	int iry = (int)round(recy) - _job->Get_Propagation_Y0();
	int irz = (int)round(recz) - _job->Get_Propagation_Z0();

	//printf("Elastic_Shot::Compute_Reciprocal_Scale_Factor - isx=%d, isy=%d, isz=%d, irx=%d, iry=%d, irz=%d\n",isx,isy,isz,irx,iry,irz);

	bool silent=true, error=false, combined_error=false;
	float rhoA = _job->Get_Earth_Model_Attribute(_job->Attr_Idx_Density, irx, iry, irz, silent, error) * 1000.0f;
	if (error) combined_error = true;
	float rhoB = _job->Get_Earth_Model_Attribute(_job->Attr_Idx_Density, isx, isy, isz, silent, error) * 1000.0f;
	if (error) combined_error = true;
	float VpA = _job->Get_Earth_Model_Attribute(_job->Attr_Idx_Vp, irx, iry, irz, silent, error);
	if (error) combined_error = true;
	float VpB = _job->Get_Earth_Model_Attribute(_job->Attr_Idx_Vp, isx, isy, isz, silent, error);
	if (error) combined_error = true;
	float VsA = _job->Get_Earth_Model_Attribute(_job->Attr_Idx_Vs, irx, iry, irz, silent, error);
	if (error) combined_error = true;
	float VsB = _job->Get_Earth_Model_Attribute(_job->Attr_Idx_Vs, isx, isy, isz, silent, error);
	if (error) combined_error = true;

	// return zero if an error was returned by earth model. zero scaling factor kills the trace.
	if (combined_error) return 0.0f;

	float kA = (VpA * VpA - 1.333333333f * VsA * VsA) * rhoA;
	float kB = (VpB * VpB - 1.333333333f * VsB * VsB) * rhoB;

	float scalefac = 1.0f;
	if (Get_Source_Type() == Source_Type_Pressure)
	{
		if (flag == 1)
		{
			//  P(A) <- P(B)
			scalefac = (kB/kA);
		}
		else if (flag == 2)
		{
			// Vx(A) <- P(B)
			scalefac = -(rhoA*kB);
		}
		else if (flag == 4)
		{
			// Vy(A) <- P(B)
			scalefac = -(rhoA*kB);
		}
		else if (flag == 8)
		{
			// Vz(A) <- P(B)
			scalefac = -(rhoA*kB);
		}
		else
		{
			printf("Elastic_Shot::Compute_Reciprocal_Scale_Factor - Cannot use reciprocity for this receiver type.\n");
			return 0.0f;
		}
	}
	else
	{
		if (Get_Amplitude1() != 0.0f && Get_Amplitude2() == 0.0f && Get_Amplitude3() == 0.0f)
		{
			// Vx source
			if (flag == 1)
			{
				//  P(A) <- Vx(B)
				scalefac = -1.0f / (rhoB*kA);
			}
			else if (flag == 2)
			{
				// Vx(A) <- Vx(B)
				scalefac = rhoA / rhoB;
			}
			else if (flag == 4)
			{
				// Vy(A) <- Vx(B)
				scalefac = rhoA / rhoB;
			}
			else if (flag == 8)
			{
				// Vz(A) <- Vx(B)
				scalefac = rhoA / rhoB;
			}
			else
			{
				printf("Elastic_Shot::Compute_Reciprocal_Scale_Factor - Cannot use reciprocity for this receiver type.\n");
				return 0.0f;
			}
		}
		else if (Get_Amplitude1() == 0.0f && Get_Amplitude2() != 0.0f && Get_Amplitude3() == 0.0f)
		{
			// Vy source
			if (flag == 1)
			{
				//  P(A) <- Vy(B)
				scalefac = -1.0f / (rhoB*kA);
			}
			else if (flag == 2)
			{
				// Vx(A) <- Vy(B)
				scalefac = rhoA / rhoB;
			}
			else if (flag == 4)
			{
				// Vy(A) <- Vy(B)
				scalefac = rhoA / rhoB;
			}
			else if (flag == 8)
			{
				// Vz(A) <- Vy(B)
				scalefac = rhoA / rhoB;
			}
			else
			{
				printf("Elastic_Shot::Compute_Reciprocal_Scale_Factor - Cannot use reciprocity for this receiver type.\n");
				return 0.0f;
			}
                }
                else if (Get_Amplitude1() == 0.0f && Get_Amplitude2() == 0.0f && Get_Amplitude3() != 0.0f)
                {
                        // Vz source
			if (flag == 1)
			{
				//  P(A) <- Vz(B)
				scalefac = -1.0f / (rhoB*kA);
			}
			else if (flag == 2)
			{
				// Vx(A) <- Vz(B)
				scalefac = rhoA / rhoB;
			}
			else if (flag == 4)
			{
				// Vy(A) <- Vz(B)
				scalefac = rhoA / rhoB;
			}
			else if (flag == 8)
			{
				// Vz(A) <- Vz(B)
				scalefac = rhoA / rhoB;
			}
			else
			{
				printf("Elastic_Shot::Compute_Reciprocal_Scale_Factor - Cannot use reciprocity for this receiver type.\n");
				return 0.0f;
			}
                }
		else
		{
			printf("Elastic_Shot::Compute_Reciprocal_Scale_Factor - Cannot use reciprocity for this velocity source.\n");
			return 0.0f;
		}
	}
	return scalefac;
}

void Elastic_Shot::Copy_Traces_To_External_Buffer(
		int nTraces,
		int nsamp,
		float* samples
		)
{
	//printf("_num_traces = %d, nTraces = %d\n",_num_traces,nTraces);
	assert(nTraces == _num_traces);
	for (int iTrc = 0;  iTrc < _num_traces;  ++iTrc)
        {
		float* dst = samples + nsamp * iTrc;
		float* src = _h_traces[iTrc]->Get_Samples();
		int nsamp_src = _h_traces[iTrc]->Get_NSAMP();
		if (nsamp_src >= nsamp)
		{
			// source is at least same size as destination, copy using destination length
			memcpy(dst,src,nsamp*sizeof(float));
		}
		else
		{
			// source is shorter than destination, copy and zero pad
			memcpy(dst,src,nsamp_src*sizeof(float));
			memset(dst+nsamp_src,0,(nsamp-nsamp_src)*sizeof(float));
		}
	}
}

void Elastic_Shot::Write_SEGY_Files()
{
	/*
        // TMJ 01/07/15 Debug outputs
	printf("Elastic_Shot::Write_SEGY_Files\n");
	float min_before =  1e37f;
	float max_before = -1e37f;
	float min_after  =  1e37f;
	float max_after  = -1e37f;
#pragma omp parallel for
	for (int iTrace = 0;  iTrace < _num_traces;  ++iTrace) 
	{
		float local_min_before, local_max_before, local_min_after, local_max_after;
		_h_traces[iTrace]->Correct_Amplitudes(local_min_before,local_max_before,local_min_after,local_max_after);
#pragma omp critical
		{
			if (local_min_before < min_before) min_before = local_min_before;
			if (local_max_before > max_before) max_before = local_max_before;
			if (local_min_after < min_after) min_after = local_min_after;
			if (local_max_after > max_after) max_after = local_max_after;
		}
	}
	printf("MIN MAX :: before=[%e,%e], after=[%e,%e]\n",min_before,max_before,min_after,max_after);
	*/
	char EBCDIC_Header[3200];
	memset(EBCDIC_Header, 32, 3200);  // fill with space
	if (_job->Get_EBCDIC_Header_Filename() != 0L)
	{
		FILE* fp = fopen((char*)(_job->Get_EBCDIC_Header_Filename()), "r");
		if (fp != 0L)
		{
			char buf[4096];
			for (int i = 0;  i < 40 && !feof(fp);  ++i)
			{
				fgets(buf,4096,fp);
				int len = strlen(buf);
				len = len > 80 ? 80 : len;  // truncate string to max 80 characters
				for (int j = 0;  j < len;  ++j)
				{
					// replace non-printable characters with space
					EBCDIC_Header[i*80+j] = isprint(buf[j]) ? buf[j] : ' ';
				}
			}
			fclose(fp);
		}
		printf("\nFormatted EBCDIC Header will be included in SEGY files.\n");
		printf("--------------------------------------------------------------------------------\n");
		printf("12345678901234567890123456789012345678901234567890123456789012345678901234567890\n");
		for (int i = 0;  i < 40;  ++i)
		{
			char buf[81];
			for (int j = 0;  j < 80;  ++j) buf[j] = EBCDIC_Header[i*80+j];
			buf[80] = 0;
			printf("%s\n",buf);
		}
		printf("--------------------------------------------------------------------------------\n\n");
	}
	Global_Coordinate_System* gcs = _job->Get_Voxet()->Get_Global_Coordinate_System();
	for (int iFile = 0;  iFile < _num_segy_files;  ++iFile)
	{
		for (int iSel = 0;  iSel < 4;  ++iSel)
		{
			int flag = 1 << iSel;
			
			// count number of traces belonging to this file with this flag
			int count = 0;
			for (int iTrc = 0;  iTrc < _num_traces;  ++iTrc)
				if (_h_traces_hdr[iTrc]->Get_File_Number() == iFile && _h_traces_hdr[iTrc]->Get_Flags() == flag)
					++count;
			char filename[4096];
			printf("%s has %d traces.\n",_segy_files[iFile]->Get_Full_Path(filename,flag),count);

			// gather output traces and write to file
			if (count > 0)
			{
				double srcx, srcy, srcz;
				gcs->Convert_Transposed_Fractional_Index_To_Global(_x,_y,_z,srcx,srcy,srcz);
				float src_model_water_depth, src_model_water_Vp;
				_job->Compute_Model_Water_Depth_And_Avg_Vp((int)round(_x-_job->Get_Propagation_X0()),(int)round(_y-_job->Get_Propagation_Y0()),src_model_water_depth,src_model_water_Vp);
				double src_bath_z;
				// HACK: src_bath_z contains the true bathymetry water depth at the source location.
				// For Agbami job, this is alway 1.5 cells shallower than the Z source placement.
				src_bath_z = fabs(srcz) - 1.5f * gcs->Get_DZ();
				float** traces = new float*[count];
				double* recx = new double[count];
				double* recy = new double[count];
				double* recz = new double[count];
				int *iline = new int[count];
				int *xline = new int[count];
				int *trcens = new int[count];
				short* compon = new short[count];
				int* irec = new int[count];
				time_t *acqtime = new time_t[count];
				int* usec = new int[count];
				float* rec_model_water_depth = new float[count];
				float* rec_model_water_Vp = new float[count];
				float* rec_bath_z = new float[count];
				int nsamp = 0;
				for (int iTrc = 0, ii = 0;  iTrc < _num_traces;  ++iTrc)
				{
					if (_h_traces_hdr[iTrc]->Get_File_Number() == iFile && _h_traces_hdr[iTrc]->Get_Flags() == flag)
					{
						nsamp = _h_traces[iTrc]->Get_NSAMP();
						traces[ii] = _h_traces[iTrc]->Get_Samples();
						gcs->Convert_Transposed_Fractional_Index_To_Global(
							_h_traces_hdr[iTrc]->Get_Location_X(),_h_traces_hdr[iTrc]->Get_Location_Y(),_h_traces_hdr[iTrc]->Get_Location_Z(),
							recx[ii],recy[ii],recz[ii]
							);
						_job->Compute_Model_Water_Depth_And_Avg_Vp(
							(int)round(_h_traces_hdr[iTrc]->Get_Location_X()-_job->Get_Propagation_X0()),(int)round(_h_traces_hdr[iTrc]->Get_Location_Y()-_job->Get_Propagation_Y0()),
							rec_model_water_depth[ii],rec_model_water_Vp[ii]
							);
						// HACK: rec_bath_z contains the true bathymetry water depth at the receiver location.
						// For Agbami job, this is always 1.5 cells shallower than the Z receiver placement.
						rec_bath_z[ii] = fabs(recz[ii]) - 1.5f * gcs->Get_DZ();
						iline[ii] = _h_traces_hdr[iTrc]->Get_Inline();
						xline[ii] = _h_traces_hdr[iTrc]->Get_Crossline();
						trcens[ii] = _h_traces_hdr[iTrc]->Get_Trace_Ensemble();
						if (_h_traces_hdr[iTrc]->Is_Pressure())
						{
							compon[ii] = 1;
						}
						else if (_h_traces_hdr[iTrc]->Is_Vz())
						{
							compon[ii] = 2;
						}
						else if (_h_traces_hdr[iTrc]->Is_Vx())
                                                {
                                                        compon[ii] = 3;
                                                }
						else if (_h_traces_hdr[iTrc]->Is_Vy())
                                                {
                                                        compon[ii] = 4;
                                                }
						irec[ii] = _h_traces_hdr[iTrc]->Get_Receiver_FFID();
						acqtime[ii] = _h_traces_hdr[iTrc]->Get_Shot_Time();
						usec[ii] = _h_traces_hdr[iTrc]->Get_Shot_Time_usec();
						++ii;
					}
				}
				_segy_files[iFile]->Write_SEGY_File(
					traces,EBCDIC_Header,_job->Is_Vwxyzt(),
					srcx,srcy,srcz,_il,_xl,recx,recy,recz,iline,xline,trcens,compon,irec,acqtime,usec,
					src_model_water_depth,src_model_water_Vp,src_bath_z,rec_model_water_depth,rec_model_water_Vp,rec_bath_z,count,nsamp,flag
					);
				delete [] rec_bath_z;
				delete [] rec_model_water_Vp;
				delete [] rec_model_water_depth;
				delete [] traces;
				delete [] recx;
				delete [] recy;
				delete [] recz;
				delete [] iline;
				delete [] xline;
				delete [] trcens;
				delete [] irec;
				delete [] compon;
				delete [] acqtime;
				delete [] usec;
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
                int num_rx = _segy_files[iFile]->Compute_Receiver_Locations_NO_COPY(rcv_x, rcv_y, rcv_z);
		for (int iRx = 0;  iRx < num_rx;  ++iRx)
		{
			if (rcv_z[iRx] > deepest) deepest = rcv_z[iRx];
		}
	}
	return deepest;
}

void Elastic_Shot::_Create_Receiver_Transfer_Buffers(Elastic_Propagator* prop)
{
        //
        // create binned receiver locations for all files
        //

        _totSteps = prop->Get_Total_Number_Of_Timesteps();
        _nBlks = prop->Get_Number_Of_Blocks();
        _num_pipes = prop->Get_Number_Of_Pipelines();

        _h_rcv_loc = new float*[_num_pipes];
        _h_rcv_binned = new float**[_num_pipes];
        _h_rcv_trcidx = new int**[_num_pipes];
        _h_rcv_trcflag = new int**[_num_pipes];
        _h_rcv_loc_size_f = new int[_num_pipes];
#pragma omp parallel for
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

	struct timespec ts1;
	clock_gettime(CLOCK_REALTIME, &ts1);

	// determine total number of traces in all files
	_num_traces = 0;
	for (int iFile = 0;  iFile < _num_segy_files;  ++iFile)
	{
		double *rcv_x, *rcv_y, *rcv_z;
		int num_rx = _segy_files[iFile]->Compute_Receiver_Locations_NO_COPY(rcv_x, rcv_y, rcv_z);
		if (num_rx > 0)
		{
			int num_wf = 0;
			int flags = _segy_files[iFile]->Get_Selection_Flags();
			for (int iSel = 0;  iSel < 4;  ++iSel) if (flags & (1 << iSel)) ++num_wf;
			_num_traces += num_wf * num_rx;
		}
	}

	struct timespec ts2;
	clock_gettime(CLOCK_REALTIME, &ts2);

	_h_traces_hdr = new Elastic_Trace_Header*[_num_traces];
	for (int iFile = 0, trace_no_idx = 0;  iFile < _num_segy_files;  ++iFile)
        {
                double *rcv_x, *rcv_y, *rcv_z;
		int *il, *xl, *rec_ffid, *trce;
		time_t *acqt;
		int* usec;
		int num_rx = _segy_files[iFile]->Compute_Receiver_Locations_NO_COPY(rcv_x, rcv_y, rcv_z, il, xl, trce, rec_ffid, acqt, usec);
		double tshift = _segy_files[iFile]->Get_Timeshift();
		double sample_rate = _segy_files[iFile]->Get_Sample_Rate();
		int nsamp_out = (int)round(_segy_files[iFile]->Get_Record_Length() / sample_rate) + 1;
		int flags = _segy_files[iFile]->Get_Selection_Flags();

		int num_wf = 0;
		for (int iSel = 0;  iSel < 4;  ++iSel) if (flags & (1 << iSel)) ++num_wf;
		if (num_wf > 0)
		{
			for (int iRx = 0;  iRx < num_rx;  ++iRx)
			{
				int curr_trace_no_idx = trace_no_idx + iRx * num_wf;
				for (int iSel = 0;  iSel < 4;  ++iSel)
				{
					if (flags & (1 << iSel))
					{
						// To-Do: Adjust & verify input start-time for velocity receiver. Make sure _dti is valid.
						_h_traces_hdr[curr_trace_no_idx] = new Elastic_Trace_Header(
							tshift,				// trace start time
							sample_rate,			// trace sample rate
							nsamp_out,			// # samples in trace, length of sinc operator
							iFile,1<<iSel,rcv_x[iRx],rcv_y[iRx],rcv_z[iRx],il[iRx],xl[iRx],trce[iRx],rec_ffid[iRx],acqt[iRx],usec[iRx]
							);
						++curr_trace_no_idx;
					}
				}
			}
			trace_no_idx += num_rx * num_wf;
		}
	}

	struct timespec ts3;
	clock_gettime(CLOCK_REALTIME, &ts3);

        for (int iPipe = 0;  iPipe < _num_pipes;  ++iPipe)
        {
                Elastic_Pipeline* pipe = prop->Get_Pipeline(iPipe);
		int y0 = pipe->Get_Y0();
		int y1 = pipe->Get_Y1();

                int tot_rx = 0, tot_traces = 0;

		// bin receivers first
		int **binRx_num = new int*[_num_segy_files];
		double ***binRx_x = new double**[_num_segy_files];
		double ***binRx_y = new double**[_num_segy_files];
		double ***binRx_z = new double**[_num_segy_files];
		int ***binRx_trcno = new int**[_num_segy_files];
		for (int iFile = 0, traceno = 0;  iFile < _num_segy_files;  ++iFile)
		{
			int flags = _segy_files[iFile]->Get_Selection_Flags();
			int num_wf = 0;
			if (flags & 1) ++num_wf;
			if (flags & 2) ++num_wf;
			if (flags & 4) ++num_wf;
			if (flags & 8) ++num_wf;
			Elastic_Interpolation_t interpolation_method = _segy_files[iFile]->Get_Interpolation_Method();
			double *rcv_x, *rcv_y, *rcv_z;
			int num_rx = _segy_files[iFile]->Compute_Receiver_Locations_NO_COPY(rcv_x, rcv_y, rcv_z);
			binRx_num[iFile] = new int[_nBlks];
			for (int iBlk = 0;  iBlk < _nBlks;  ++iBlk) binRx_num[iFile][iBlk] = 0;
			for (int iRx = 0;  iRx < num_rx;  ++iRx)
			{
				int ix = (int)truncf((float)(rcv_x[iRx] - _job->Get_Propagation_X0()));
				int ixBlk = ix / prop->Get_Block_Size_X();
				int iBlk_min = ixBlk-2;
				if (iBlk_min < 0) iBlk_min = 0;
				int iBlk_max = ixBlk+2;
				if (iBlk_max >= _nBlks) iBlk_max = _nBlks-1;
				//printf("rcv[%d][%d] = %f,%f,%f :: iBlk=[%d,%d]\n",iFile,iRx,rcv_x[iRx],rcv_y[iRx],rcv_z[iRx],iBlk_min,iBlk_max);
				for (int iBlk = iBlk_min;  iBlk <= iBlk_max;  ++iBlk)
				{
					int x0 = iBlk * prop->Get_Block_Size_X();
					int x1 = x0 + prop->Get_Block_Size_X() - 1;
					//printf(" => x0=%d, x1=%d\n",x0,x1);
					if (_Receiver_Intersects(interpolation_method,x0,x1,y0,y1,rcv_x[iRx]-_job->Get_Propagation_X0(),rcv_y[iRx]-_job->Get_Propagation_Y0()))
					{
						++binRx_num[iFile][iBlk];
						++tot_rx;
						tot_traces += num_wf;
					}
				}
			}
			binRx_x[iFile] = new double*[_nBlks];
			binRx_y[iFile] = new double*[_nBlks];
			binRx_z[iFile] = new double*[_nBlks];
			binRx_trcno[iFile] = new int*[_nBlks];
			for (int iBlk = 0;  iBlk < _nBlks;  ++iBlk)
			{
				//printf("binRx_num[%d][%d] = %d\n",iFile,iBlk,binRx_num[iFile][iBlk]);
				if (binRx_num[iFile][iBlk] > 0)
				{
					binRx_x[iFile][iBlk] = new double[binRx_num[iFile][iBlk]];
					binRx_y[iFile][iBlk] = new double[binRx_num[iFile][iBlk]];
					binRx_z[iFile][iBlk] = new double[binRx_num[iFile][iBlk]];
					binRx_trcno[iFile][iBlk] = new int[binRx_num[iFile][iBlk]];
				}
				else
				{
					binRx_x[iFile][iBlk] = 0L;
					binRx_y[iFile][iBlk] = 0L;
					binRx_z[iFile][iBlk] = 0L;
					binRx_trcno[iFile][iBlk] = 0L;
				}
				binRx_num[iFile][iBlk] = 0;
			}
			for (int iRx = 0;  iRx < num_rx;  ++iRx, traceno+=num_wf)
			{
				int ix = (int)truncf((float)(rcv_x[iRx] - _job->Get_Propagation_X0()));
				int ixBlk = ix / prop->Get_Block_Size_X();
				int iBlk_min = ixBlk-2;
				if (iBlk_min < 0) iBlk_min = 0;
				int iBlk_max = ixBlk+2;
				if (iBlk_max >= _nBlks) iBlk_max = _nBlks-1;
				for (int iBlk = iBlk_min;  iBlk <= iBlk_max;  ++iBlk)
				{
					int x0 = iBlk * prop->Get_Block_Size_X();
					int x1 = x0 + prop->Get_Block_Size_X() - 1;
					if (_Receiver_Intersects(interpolation_method,x0,x1,y0,y1,rcv_x[iRx]-_job->Get_Propagation_X0(),rcv_y[iRx]-_job->Get_Propagation_Y0()))
					{
						binRx_x[iFile][iBlk][binRx_num[iFile][iBlk]] = rcv_x[iRx];
						binRx_y[iFile][iBlk][binRx_num[iFile][iBlk]] = rcv_y[iRx];
						binRx_z[iFile][iBlk][binRx_num[iFile][iBlk]] = rcv_z[iRx];
						binRx_trcno[iFile][iBlk][binRx_num[iFile][iBlk]] = traceno;
						++binRx_num[iFile][iBlk];
					}
				}
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
                        _h_rcv_binned[iPipe][iBlk] = floc;
                        iloc[0] = _num_segy_files;
                        _h_rcv_trcidx[iPipe][iBlk] = trcidx;
                        _h_rcv_trcflag[iPipe][iBlk] = trcflag;
                        int blk_offset = 0;
                        int trace_no_idx = 0;
                        for (int iFile = 0;  iFile < _num_segy_files;  ++iFile)
                        {
                                int flags = _segy_files[iFile]->Get_Selection_Flags();
                                Elastic_Interpolation_t interpolation_method = _segy_files[iFile]->Get_Interpolation_Method();
                                bool receiver_ghost_enabled = !_job->Freesurface_Enabled() && _job->Receiver_Ghost_Enabled();
				double *rcv_x = binRx_x[iFile][iBlk];
				double *rcv_y = binRx_y[iFile][iBlk];
				double *rcv_z = binRx_z[iFile][iBlk];
				int* traceno = binRx_trcno[iFile][iBlk];
				int num_rx = binRx_num[iFile][iBlk];
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
						int off = 1 + 2 * _num_segy_files + 3 * (blk_offset + file_offset);
						floc[off  ] = rcv_x[iRx] - _job->Get_Propagation_X0();
						floc[off+1] = rcv_y[iRx] - _job->Get_Propagation_Y0();
						floc[off+2] = rcv_z[iRx] - _job->Get_Propagation_Z0();
						++file_offset;

						for (int iSel = 0, trace_no_off = traceno[iRx];  iSel < 4;  ++iSel)
						{
							if (flags & (1 << iSel))
							{
								trcidx[trace_no_idx] = trace_no_off++;
								trcflag[trace_no_idx] = 1 << iSel;
								++trace_no_idx;
							}
						}
					}
                                        iloc[1+2*iFile] = file_offset;
                                        blk_offset += file_offset;
                                }
                        }
                        //printf("iBlk = %d, trace_no_idx = %d\n",iBlk,trace_no_idx);
                        int glob_offset_inc = 1 + 2 * _num_segy_files + 3 * blk_offset;
                        floc = floc + glob_offset_inc;
                        iloc = iloc + glob_offset_inc;
                        trcidx = trcidx + trace_no_idx;
                        trcflag = trcflag + trace_no_idx;
                }

		for (int iFile = 0;  iFile < _num_segy_files;  ++iFile)
		{
			for (int iBlk = 0;  iBlk < _nBlks;  ++iBlk)
			{
				int num_rx = binRx_num[iFile][iBlk];
				if (num_rx > 0)
				{
					delete [] binRx_x[iFile][iBlk];
					delete [] binRx_y[iFile][iBlk];
					delete [] binRx_z[iFile][iBlk];
					delete [] binRx_trcno[iFile][iBlk];
				}
			}
			delete [] binRx_x[iFile];
			delete [] binRx_y[iFile];
			delete [] binRx_z[iFile];
			delete [] binRx_trcno[iFile];
			delete [] binRx_num[iFile];
		}
		delete [] binRx_x;
		delete [] binRx_y;
		delete [] binRx_z;
		delete [] binRx_trcno;
		delete [] binRx_num;
        }

#ifdef DEBUG_TMJ
	FILE* fp = fopen("/cpfs/lfs01/ESDRD/tjhc/rx_locs.txt", "w");
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
					Elastic_Interpolation_t interpolation_method = (Elastic_Interpolation_t)((flags >> 16) & 3);
					char* str_intrp = (char*)ToString_Elastic_Interpolation_t(interpolation_method);
					bool receiver_ghost_enabled = ((flags & (1 << 31)) != 0) ? true : false;
					fprintf(fp,"      F I L E   %d - %d receivers, RX Ghost %s, %s interpolation, selected %s%s%s%s\n\n",iFile,num_rx,receiver_ghost_enabled?"enabled":"disabled",str_intrp,(flags&1)?"P ":"",(flags&2)?"Vx ":"",(flags&4)?"Vy ":"",(flags&8)?"Vz ":"");
					for (int iRx = 0;  iRx < num_rx;  ++iRx)
					{
						float x = floc[bulk_offset+3*iRx  ];
						float y = floc[bulk_offset+3*iRx+1];
						float z = floc[bulk_offset+3*iRx+2];
						fprintf(fp,"         %e %e %e :: trace",x,y,z);
						for (int iSel = 0;  iSel < 4;  ++iSel)
						{
							if (flags & (1 << iSel))
							{
								char* FlagStr = (char*)"";
								if (trcflag[trace_no_idx] & 1) FlagStr = (char*)"P";
								else if (trcflag[trace_no_idx] & 2) FlagStr = (char*)"Vx";
								else if (trcflag[trace_no_idx] & 4) FlagStr = (char*)"Vy";
								else if (trcflag[trace_no_idx] & 8) FlagStr = (char*)"Vz";
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

void Elastic_Shot::Create_Trace_Resample_Buffers(Elastic_Propagator* prop)
{
	double dti = prop->Get_Internal_Sample_Rate();
	_dti = dti;

	_h_traces = new Elastic_Resampled_Trace*[_num_traces];
	for (int iTrc = 0;  iTrc < _num_traces;  ++iTrc)
	{
		_h_traces[iTrc] = new Elastic_Resampled_Trace(
			_h_traces_hdr[iTrc]->Is_Velocity() ? (_dti/2.0) : 0.0,	// wavefield start time
			_dti,							// wavefield timestep
			_h_traces_hdr[iTrc]->Get_Start_Time(),			// trace start time
			_h_traces_hdr[iTrc]->Get_Sample_Rate(),			// trace sample rate
			_h_traces_hdr[iTrc]->Get_NSAMP(),			// # samples in trace
			2							// length of sinc operator, ==2 yields linear interpolation
			);
	}
}

void Elastic_Shot::Free_Trace_Resample_Buffers()
{
	if (_num_traces > 0)
	{
		if (_h_traces != 0L)
		{
			for (int iTrc = 0;  iTrc < _num_traces;  ++iTrc)
                        {
				delete _h_traces[iTrc];
			}
			delete [] _h_traces;
			_h_traces = 0L;
		}
		if (_h_traces_hdr != 0L)
		{
			for (int iTrc = 0;  iTrc < _num_traces;  ++iTrc)
                        {
                                delete _h_traces_hdr[iTrc];
                        }
                        delete [] _h_traces_hdr;
                        _h_traces_hdr = 0L;
		}
		_num_traces = 0;
	}
	if (_num_pipes > 0)
	{
		for (int iPipe = 0;  iPipe < _num_pipes;  ++iPipe)
		{
			if (_h_rcv_loc != 0L) delete [] _h_rcv_loc[iPipe];
			if (_h_rcv_binned != 0L) delete [] _h_rcv_binned[iPipe];
			if (_h_rcv_trcidx != 0L)
			{
				if (_h_rcv_trcidx[iPipe] != 0L) delete [] _h_rcv_trcidx[iPipe][0];
				delete [] _h_rcv_trcidx[iPipe];
			}
			if (_h_rcv_trcflag != 0L)
			{
				if (_h_rcv_trcflag[iPipe] != 0L) delete [] _h_rcv_trcflag[iPipe][0];
				delete [] _h_rcv_trcflag[iPipe];
			}
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

