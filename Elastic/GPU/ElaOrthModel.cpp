#include <cstring>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <map>
#include <ElaOrthModel.hxx>
#include <Elastic_Modeling_Job.hxx>
#include <Elastic_Propagator.hxx>
#include <Elastic_Shot.hxx>
#include <Voxet.hxx>
#include <Voxet_Property.hxx>
#include <Voxet_Memory_Mapper.hxx>
#include <Global_Coordinate_System.hxx>
#include <Variable_Water_Velocity.hxx>

using namespace std;

static void Format_Voxet_Property(
	ostringstream& oss,
	map<string,string>& key_value_pairs,
	map<string,double>& min,
	map<string,double>& max,
	string key
	)
{
	string val = key_value_pairs[key];
	char* endptr;
	double dval = strtod(val.c_str(),&endptr);
	if (dval == 0.0 && *endptr != 0)
	{
		// value is a string
		if (min.find(val) != min.end())
		{
			oss << "PROPERTY " << key << " = " << val << " " << min[val] << " " << max[val] << endl;
		}
		else
		{
			oss << "PROPERTY " << key << " = " << val << " 0 0" << endl;
		}
	}
	else
	{
		// value is a valid constant
		oss << "PROPERTY " << key << " = " << val << endl;
	}
}

// translate "yes", "1", "true" to enabled.
// translate "no", "0", "false" to disabled.
// anything else translates to bad_value.
static string Translate_Boolean_Choice(string val)
{
	if (val == "yes" || val == "1" || val == "true")
	{
		return (string)"enabled";
	}
	else if (val == "no" || val == "0" || val == "true")
	{
		return (string)"disabled";
	}
	else
	{
		return (string)"bad_value";
	}
}

static string Translate_Anchor(string val)
{
	if (val == "source_and_receivers")
	{
		return (string)"Source_And_Receivers";
	}
	else if (val == "source")
	{
		return (string)"Source";
	}
	else if (val == "receivers")
	{
		return (string)"Receivers";
	}
	else if (val == "model_origin")
	{
		return (string)"Volume";
	}
	else
	{
		return (string)"bad_value";
	}
}

ElaOrthModel::ElaOrthModel(
		char* parmString
		)
{
	_job = 0L;
	_parmString = strdup(parmString);
	_nsamp = 0;
}

void ElaOrthModel::Get_Trace_Arguments(
	int& samples_per_trace,
	double& sample_rate_in_seconds
	)
{
	map<string,string> key_value_pairs;
	istringstream iss(_parmString);
	char buf[1024];
        for (iss.getline(buf,1024);  !iss.eof();  iss.getline(buf,1024))
        {
                char key[1024];
                char val[1024];
                if (sscanf(buf,"%s = %s\n",key,val) == 2)
                {
                        key_value_pairs[(string)key] = (string)val;
                }
	}
	samples_per_trace = atoi(key_value_pairs["N_SAMPLES_PER_TRACE"].c_str());
	sample_rate_in_seconds = atof(key_value_pairs["SAMPLE_INTERVAL"].c_str()) / 1000.0;
}

Voxet* ElaOrthModel::Get_Voxet()
{
	if (_job == 0L) _job = _Create_Job_Object(0.0f,0.0f,0.0f);
	return _job->Get_Voxet();
}

Elastic_Modeling_Job* ElaOrthModel::_Create_Job_Object(float soux, float souy, float souz)
{
	map<string,string> key_value_pairs;

	// insert defaults
	key_value_pairs["Transpose"] = "zxy";
	key_value_pairs["Anchor"] = "source_and_receivers";
	key_value_pairs["Model_Type"] = "Isotropic";
	key_value_pairs["Vp"] = "p_velocity";
	key_value_pairs["Density"] = "density";
	key_value_pairs["Q"] = "1e6";
	key_value_pairs["FQ"] = "10";
	key_value_pairs["Scholte_Attenuation"] = "no";
	key_value_pairs["Extend_Model"] = "yes";
	key_value_pairs["Freesurface"] = "yes";
	key_value_pairs["Source_Ghost"] = "yes";
	key_value_pairs["Receiver_Ghost"] = "yes";
	key_value_pairs["Timestep_Multiplier"] = "1.0";
	key_value_pairs["Spatial_Order"] = "8";
	key_value_pairs["Spatial_Interpolator"] = "Sinc";
	key_value_pairs["Source_Type"] = "P";
	key_value_pairs["Receiver_Type"] = "P";
	key_value_pairs["Gather_Type"] = "Common_Shot";
	key_value_pairs["Source_Time_Shift"] = "0.0";
	key_value_pairs["Verbosity"] = "2";
	
	// overwrite defaults with values from parmString
	istringstream iss(_parmString);
	char buf[1024];
        for (iss.getline(buf,1024);  !iss.eof();  iss.getline(buf,1024))
        {
		char* splt_pt = strstr(buf,"=");
		if (splt_pt == 0L)
		{
			printf("Error! Invalid key-value pair (%s)\n",buf);
			exit(-1);
		}
		else
		{
			char* key = buf;
			char* val = splt_pt + 1;
			while (splt_pt[-1] == ' ' && splt_pt > key) --splt_pt;
			*splt_pt = '\0';
			while (*val == ' ') ++val;
			printf("KEY = '%s' :: VALUE = '%s'\n",key,val);
			key_value_pairs[key] = val;
		}
	}

	// translate voxet attribute keys to fit model type
	if (key_value_pairs["Model_Type"] == (string)"Isotropic")
	{
		key_value_pairs["Vs"] = "0";
		key_value_pairs["Dip"] = "0";
		key_value_pairs["Azimuth"] = "0";
		key_value_pairs["Rake"] = "0";
		key_value_pairs["Delta1"] = "0";
		key_value_pairs["Delta2"] = "0";
		key_value_pairs["Delta3"] = "0";
		key_value_pairs["Epsilon1"] = "0";
		key_value_pairs["Epsilon2"] = "0";
		key_value_pairs["Gamma1"] = "0";
		key_value_pairs["Gamma2"] = "0";
	}
	else if (key_value_pairs["Model_Type"] == (string)"Acoustic VTI")
	{
		key_value_pairs["Vs"] = "0";
		key_value_pairs["Dip"] = "0";
		key_value_pairs["Azimuth"] = "0";
		key_value_pairs["Rake"] = "0";
		key_value_pairs["Delta2"] = key_value_pairs["Delta1"];
		key_value_pairs["Delta3"] = "0";
		key_value_pairs["Epsilon1"] = "0";
		key_value_pairs["Epsilon2"] = key_value_pairs["Epsilon1"];
		key_value_pairs["Gamma1"] = "0";
		key_value_pairs["Gamma2"] = "0";
	}
	else if (key_value_pairs["Model_Type"] == (string)"Acoustic TTI")
	{
		key_value_pairs["Vs"] = "0";
		key_value_pairs["Rake"] = "0";
		key_value_pairs["Delta2"] = key_value_pairs["Delta1"];
		key_value_pairs["Delta3"] = "0";
		key_value_pairs["Epsilon1"] = "0";
		key_value_pairs["Epsilon2"] = key_value_pairs["Epsilon1"];
		key_value_pairs["Gamma1"] = "0";
		key_value_pairs["Gamma2"] = "0";
	}
	else if (key_value_pairs["Model_Type"] == (string)"Elastic Isotropic")
	{
		key_value_pairs["Dip"] = "0";
		key_value_pairs["Azimuth"] = "0";
		key_value_pairs["Rake"] = "0";
		key_value_pairs["Delta1"] = "0";
		key_value_pairs["Delta2"] = "0";
		key_value_pairs["Delta3"] = "0";
		key_value_pairs["Epsilon1"] = "0";
		key_value_pairs["Epsilon2"] = "0";
		key_value_pairs["Gamma1"] = "0";
		key_value_pairs["Gamma2"] = "0";
	}
	else if (key_value_pairs["Model_Type"] == (string)"Elastic VTI")
	{
		key_value_pairs["Dip"] = "0";
		key_value_pairs["Azimuth"] = "0";
		key_value_pairs["Rake"] = "0";
		key_value_pairs["Delta2"] = key_value_pairs["Delta1"];
		key_value_pairs["Delta3"] = "0";
		key_value_pairs["Epsilon1"] = "0";
		key_value_pairs["Epsilon2"] = key_value_pairs["Epsilon1"];
		key_value_pairs["Gamma1"] = "0";
		key_value_pairs["Gamma2"] = "0";
	}
	else if (key_value_pairs["Model_Type"] == (string)"Elastic TTI")
	{
		key_value_pairs["Rake"] = "0";
		key_value_pairs["Delta2"] = key_value_pairs["Delta1"];
		key_value_pairs["Delta3"] = "0";
		key_value_pairs["Epsilon1"] = "0";
		key_value_pairs["Epsilon2"] = key_value_pairs["Epsilon1"];
		key_value_pairs["Gamma1"] = "0";
		key_value_pairs["Gamma2"] = "0";
	}
	else if (key_value_pairs["Model_Type"] == (string)"Elastic VOr")
	{
		key_value_pairs["Dip"] = "0";
		key_value_pairs["Azimuth"] = "0";
		key_value_pairs["Rake"] = "0";
	}
	else if (key_value_pairs["Model_Type"] == (string)"Elastic TOr")
	{
		// all fields are valid for this mode, so do nothing
	}
	
	// generate parmfile from key value pairs
	ostringstream oss;
	oss << "USE VOXET " << key_value_pairs["Voxet"] << endl;
	oss << "TRANSPOSE UVW = " << key_value_pairs["Transpose"] << endl;
	oss << "MODEL_TYPE = " << key_value_pairs["Model_Type"] << endl;
	Format_Voxet_Property(oss,key_value_pairs,_min_prop_val,_max_prop_val,"Vp");
	Format_Voxet_Property(oss,key_value_pairs,_min_prop_val,_max_prop_val,"Vs");
	Format_Voxet_Property(oss,key_value_pairs,_min_prop_val,_max_prop_val,"Density");
	Format_Voxet_Property(oss,key_value_pairs,_min_prop_val,_max_prop_val,"Q");
	oss << "SET FQ = " << key_value_pairs["FQ"] << endl;
	oss << "LOWER_Q_ALONG_SEAFLOOR " << Translate_Boolean_Choice(key_value_pairs["Scholte_Attenuation"]) << endl;
	Format_Voxet_Property(oss,key_value_pairs,_min_prop_val,_max_prop_val,"Dip");
	Format_Voxet_Property(oss,key_value_pairs,_min_prop_val,_max_prop_val,"Azimuth");
	Format_Voxet_Property(oss,key_value_pairs,_min_prop_val,_max_prop_val,"Rake");
	Format_Voxet_Property(oss,key_value_pairs,_min_prop_val,_max_prop_val,"Delta1");
	Format_Voxet_Property(oss,key_value_pairs,_min_prop_val,_max_prop_val,"Delta2");
	Format_Voxet_Property(oss,key_value_pairs,_min_prop_val,_max_prop_val,"Delta3");
	Format_Voxet_Property(oss,key_value_pairs,_min_prop_val,_max_prop_val,"Epsilon1");
	Format_Voxet_Property(oss,key_value_pairs,_min_prop_val,_max_prop_val,"Epsilon2");
	Format_Voxet_Property(oss,key_value_pairs,_min_prop_val,_max_prop_val,"Gamma1");
	Format_Voxet_Property(oss,key_value_pairs,_min_prop_val,_max_prop_val,"Gamma2");
	oss << "PROPAGATE_ORIGIN = " << Translate_Anchor(key_value_pairs["Anchor"]) << endl;
	oss << "PROPAGATE_X = -" << key_value_pairs["Aperture_X"] << " " << key_value_pairs["Aperture_X"] << " local" << endl;
	oss << "PROPAGATE_Y = -" << key_value_pairs["Aperture_Y"] << " " << key_value_pairs["Aperture_Y"] << " local" << endl;
	oss << "PROPAGATE_Z = 0 100 %" << endl;
	oss << "PROPAGATE_EXTEND_MODEL " << Translate_Boolean_Choice(key_value_pairs["Extend_Model"]) << endl;
	oss << "NABC_SDX = " << key_value_pairs["ABC_X"] << " local extend" << endl;
	oss << "NABC_SDY = " << key_value_pairs["ABC_Y"] << " local extend" << endl;
	oss << "NABC_TOP = " << key_value_pairs["ABC_TOP"] << " local extend" << endl;
	oss << "NABC_BTM = " << key_value_pairs["ABC_BTM"] << " local extend" << endl;
	oss << "FREESURFACE = " << Translate_Boolean_Choice(key_value_pairs["Freesurface"]) << endl;
	oss << "SOURCE_GHOST = " << Translate_Boolean_Choice(key_value_pairs["Source_Ghost"]) << endl;
	oss << "RECEIVER_GHOST = " << Translate_Boolean_Choice(key_value_pairs["Receiver_Ghost"]) << endl;
	oss << "COURANT_FACTOR = " << key_value_pairs["Timestep_Multiplier"] << endl;
	oss << "SPATIAL_ORDER = " << key_value_pairs["Spatial_Order"] << endl;
	oss << "SHOT 1 SOURCE_LOCATION " << soux << " " << souy << " " << souz << " global" << endl;
	if (key_value_pairs["Source_Type"] == "P")
	{
		oss << "SHOT 1 SOURCE_TYPE PRESSURE" << endl;
		oss << "SHOT 1 SOURCE_AMPLITUDE 1.0" << endl;
	}
	else if (key_value_pairs["Source_Type"] == "Vx")
	{
		oss << "SHOT 1 SOURCE_TYPE VELOCITY" << endl;
		oss << "SHOT 1 SOURCE_AMPLITUDE 1.0 0.0 0.0" << endl;
	}
	else if (key_value_pairs["Source_Type"] == "Vy")
	{
		oss << "SHOT 1 SOURCE_TYPE VELOCITY" << endl;
		oss << "SHOT 1 SOURCE_AMPLITUDE 0.0 1.0 0.0" << endl;
	}
	else if (key_value_pairs["Source_Type"] == "Vz")
	{
		oss << "SHOT 1 SOURCE_TYPE VELOCITY" << endl;
		oss << "SHOT 1 SOURCE_AMPLITUDE 0.0 0.0 1.0" << endl;
	}
	oss << "SHOT 1 SOURCE_INTERPOLATION " << key_value_pairs["Spatial_Interpolator"] << endl;
	double fmax;
	if (sscanf(key_value_pairs["Source_Wavelet"].c_str(), "ricker %lf", &fmax) == 1)
	{
		oss << "SHOT 1 SOURCE_WAVELET Ricker " << fmax << endl;
		double fpeak = fmax / 2.77;
		double tshift = 5.0 * sqrt(1.5) / (fpeak * 3.1415926535897932384626433832795);
		char str[256];
		sprintf(str,"%.14e",tshift*1e3);
		key_value_pairs["Source_Time_Shift"] = str;
	}
	else
	{
		oss << "SHOT 1 SOURCE_WAVELET FILE " << key_value_pairs["Source_Wavelet"] << " DONT_FILTER" << endl;
	}
	_nsamp = atoi(key_value_pairs["N_SAMPLES_PER_TRACE"].c_str());
	double sample_rate = atof(key_value_pairs["SAMPLE_INTERVAL"].c_str()) / 1000.0;
	double start_time = atof(key_value_pairs["Source_Time_Shift"].c_str()) / 1000.0;
	double end_time = start_time + (double)(_nsamp-1) * sample_rate;
	char sample_rate_str[1024];
	sprintf(sample_rate_str,"%.6f",sample_rate);
	char start_str[1024];
	sprintf(start_str,"%.6f",start_time);
	char end_str[1024];
	sprintf(end_str,"%.6f",end_time);
	oss << "SHOT 1 SEGY_FILE 9999 FILE dummy " << sample_rate_str << " " <<  start_str << " " << end_str << " " << key_value_pairs["Receiver_Type"] << endl;
	//oss << "SHOT 1 SEGY_FILE 9999 RECEIVER_LOCATIONS 9999 RANGE_X 0 0 1 local" << endl;
	//oss << "SHOT 1 SEGY_FILE 9999 RECEIVER_LOCATIONS 9999 RANGE_Y 0 0 1 local" << endl;
	//oss << "SHOT 1 SEGY_FILE 9999 RECEIVER_LOCATIONS 9999 RANGE_Z 0 0 1 local" << endl;
	
	cout << oss.str() << endl;
	istringstream iss2(oss.str());
	int verbosity = atoi(key_value_pairs["Verbosity"].c_str());
	_mapper = new Voxet_Memory_Mapper();
	_Vwxyzt = new Variable_Water_Velocity();
	Elastic_Modeling_Job* job = new Elastic_Modeling_Job(verbosity,"ElaOrthoModel",_mapper,_Vwxyzt,iss2);
	for (int iProp = 0;  iProp < job->Get_Number_of_Voxet_Properties();  ++iProp)
	{
		Voxet_Property* prop = job->Get_Voxet_Property(iProp);
		if (prop != 0L)
		{
			string moniker = (string)(prop->Get_Moniker());
			double min_val = prop->Get_Min();
			double max_val = prop->Get_Max();
			_min_prop_val[moniker] = min_val;
			_max_prop_val[moniker] = max_val;
			printf(" => Stored Voxet_Property(%s,min=%e,max=%e)\n",moniker.c_str(),min_val,max_val);
		}
	}
	return job;
}

ElaOrthModel::~ElaOrthModel()
{
	delete _mapper;
	delete _Vwxyzt;
	delete _job;
	free(_parmString);
}

void ElaOrthModel::runWorker()
{
	// entry point for MPI worker processes
	// do nothing for now
}

int ElaOrthModel::runShot(
		int	nTraces,
		float	soux,
		float	souy,
		float	souz,
		float*	recx,
		float*	recy,
		float*	recz,
		float*	samples
		)
{
	Elastic_Modeling_Job* job = _Create_Job_Object(soux,souy,souz);
	//return 0;  // testing

	Global_Coordinate_System* gcs = job->Get_Voxet()->Get_Global_Coordinate_System();

	Elastic_Shot* shot = job->Get_Shot_By_Index(0);
	double isx = shot->Get_Source_X();
	double isy = shot->Get_Source_Y();
	double isz = shot->Get_Source_Z();
	if (isx < 0.0 || isx > gcs->Get_NX() || isy < 0.0 || isy > gcs->Get_NY() || isz < 0.0 || isz > gcs->Get_NZ())
	{
		printf("ElaOrthModel::runShot - Error! Source lies outside model bounds.\nAborting, all output traces will be zeroed out.\n");
		memset((void*)samples,0,(long)sizeof(float)*(long)nTraces*(long)_nsamp);
		return -1;
	}
	else
	{
		double *drecx = new double[nTraces];
		double *drecy = new double[nTraces];
		double *drecz = new double[nTraces];
		int* iline = new int[nTraces];
		int* xline = new int[nTraces];
		int* trcens = new int[nTraces];
		int* rec_ffid = new int[nTraces];
		time_t* acqtime = new time_t[nTraces];
		int* usec = new int[nTraces];
		time_t ima = time(0L);  // TO-DO: Copy shot time from SeisSpace trace header.
		for (int iRec = 0;  iRec < nTraces;  ++iRec)
		{
			gcs->Convert_Global_To_Transposed_Fractional_Index(recx[iRec],recy[iRec],recz[iRec],drecx[iRec],drecy[iRec],drecz[iRec]);
			iline[iRec] = iRec + 1;
			xline[iRec] = 1;
			trcens[iRec] = 1;
			rec_ffid[iRec] = 0;
			acqtime[iRec] = ima;
			usec[iRec] = 0;
		}
		Elastic_Shot* shot = job->Get_Shot_By_Index(0);
		shot->Add_Receiver_Array(nTraces,drecx,drecy,drecz,iline,xline,trcens,rec_ffid,acqtime,usec);
		delete [] drecx;
		delete [] drecy;
		delete [] drecz;
		delete [] iline;
		delete [] xline;
		delete [] trcens;
		delete [] rec_ffid;
		delete [] acqtime;
		job->Compute_Subvolume();

		/*
		char job_dim_str[256];
		sprintf(job_dim_str,"%dx%dx%d",job->Get_Propagation_NX(),job->Get_Propagation_NY(),job->Get_Propagation_NZ());
		if (_gpu_device_params.find((string)job_dim_str) != _gpu_device_params.end())
		{
			GPU_Runtime_Parameters* params = _gpu_device_params[(string)job_dim_str];
			job->Set_Number_Of_GPU_Pipes(params->Get_Number_Of_GPU_Pipes());
			job->Set_Steps_Per_GPU(params->Get_Steps_Per_GPU());
			job->Set_GPU_Devices(params->Get_GPU_Devices(),params->Get_Number_Of_GPU_Devices());
			printf(" ** Reusing GPU device parameters from previous run **\n");
			printf(" ** %s\n",job_dim_str);
			printf(" ** Number_Of_GPU_Pipes = %d\n",params->Get_Number_Of_GPU_Pipes());
			printf(" ** Steps_Per_GPU = %d\n",params->Get_Steps_Per_GPU());
			printf(" ** Device_IDs = ");
			for (int iDev = 0;  iDev < params->Get_Number_Of_GPU_Devices();  ++iDev)
			{
				if (iDev == 0)
					printf("%d",params->Get_GPU_Devices()[iDev]);
				else
					printf(",%d",params->Get_GPU_Devices()[iDev]);
			}
			printf("\n");
		}
		*/

		Elastic_Propagator* prop = new Elastic_Propagator(job);
		if (prop != 0L)
		{
			prop->Configure();
			prop->Read_Earth_Model();
			if (job->Lower_Q_Seafloor_Enabled()) job->Lower_Q_Seafloor(job->Scholte_Only());

			prop->Propagate_Shot_Dont_Write_SEGY(shot,false,false);
			shot->Copy_Traces_To_External_Buffer(nTraces,_nsamp,samples);

			/*
			if (_gpu_device_params.find((string)job_dim_str) == _gpu_device_params.end())
			{
				GPU_Runtime_Parameters* params = new GPU_Runtime_Parameters(prop->Get_Number_Of_Pipelines(),prop->Get_Steps_Per_GPU(),prop->Get_Number_Of_GPU_Devices(),prop->Get_GPU_Devices());
				_gpu_device_params[(string)job_dim_str] = params;
				printf(" ** Storing GPU device parameters for reuse in succeeding runs **\n");
				printf(" ** %s\n",job_dim_str);
				printf(" ** Number_Of_GPU_Pipes = %d\n",params->Get_Number_Of_GPU_Pipes());
				printf(" ** Steps_Per_GPU = %d\n",params->Get_Steps_Per_GPU());
				printf(" ** Device_IDs = ");
				for (int iDev = 0;  iDev < params->Get_Number_Of_GPU_Devices();  ++iDev)
				{
					if (iDev == 0)
						printf("%d",params->Get_GPU_Devices()[iDev]);
					else
						printf(",%d",params->Get_GPU_Devices()[iDev]);
				}
				printf("\n");
			}
			*/
			prop->Release_Resources_After_Propagation(shot);

			prop->Free_Host_Memory();
			prop->Free_Device_Memory();
			delete prop;
		}
		delete job;
		return 0;
	}
}

