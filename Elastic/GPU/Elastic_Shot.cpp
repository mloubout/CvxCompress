#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "Elastic_Shot.hxx"
#include "Elastic_Propagator.hxx"
#include "Elastic_Modeling_Job.hxx"
#include "Elastic_Buffer.hxx"
#include "Elastic_SEGY_File.hxx"

Elastic_Shot::Elastic_Shot(int log_level, Elastic_Modeling_Job* job, int souidx, double x, double y, double z)
{
	_log_level = log_level;
	_job = job;
	_ordertime = 2;  // hardcode for now.
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
	_tsrc = 0;
	_stf = new double[4096];
	_propagation_time = 0.0;
	_segy_files = 0L;
	_num_segy_files = 0;
}

Elastic_Shot::~Elastic_Shot()
{
	if (_stf != 0L) delete [] _stf;
	if (_segy_files != 0L)
	{
		for (int i = 0;  i < _num_segy_files;  ++i)
		{
			delete _segy_files[i];
		}
		delete [] _segy_files;
	}
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

bool Elastic_Shot::Set_Propagation_Time(double PropTime, const char* parmfile_path, int line_num)
{
	if (PropTime <= 0.0)
	{
		printf("%s (line %d): Error - PROPAGATION_TIME is in seconds and must be positive number.\n", parmfile_path, line_num);
		return true;
	}
	else
	{
		_propagation_time = PropTime;
		return false;
	}
}

void Elastic_Shot::Prepare_Source_Wavelet(double dt)
{
	if (_wavetype == 1)
	{
		_src(dt, _max_freq, 1, 0L, &_tsrc, _stf);
		//for (int i = 0;  i < _tsrc;  ++i) _stf[i] = 0.0f;
		//_stf[3] = 1.0f;
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
	else
	{
		printf("%s (line %d): Error - SOURCE_WAVELET unknown wavelet type %s.\n", parmfile_path, line_num, wavetype);
		return true;
	}
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
				printf("stf[%d] = %e\n",i,stf[i]);
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
                        if (_log_level >= 4) printf("stf[%d] = %e\n",i,stf[i]);
                }
        }
}

void Elastic_Shot::Shift_Receiver_Transfer_Buffers()
{
	for (int iFile = 0;  iFile < _num_segy_files;  ++iFile)
	{
		_segy_files[iFile]->Shift_Receiver_Transfer_Buffers();
	}
}

void Elastic_Shot::Create_Receiver_Transfer_Buffers(Elastic_Propagator* prop)
{
	for (int iFile = 0;  iFile < _num_segy_files;  ++iFile)
	{
		_segy_files[iFile]->Create_Receiver_Transfer_Buffers(prop);
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

