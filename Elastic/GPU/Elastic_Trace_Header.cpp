#include "Elastic_Trace_Header.hxx"

Elastic_Trace_Header::Elastic_Trace_Header(
		double start_time,
		double sample_rate,
		int nsamp,
		int iFile,
                int flags,
                double loc_x,
                double loc_y,
                double loc_z,
                int iline,
                int xline,
                int trcens
                )
{
	_start_time = start_time;
	_sample_rate = sample_rate;
	_nsamp = nsamp;
	_iFile = iFile;
	_flags = flags;
	_loc_x = loc_x;
	_loc_y = loc_y;
	_loc_z = loc_z;
	_iline = iline;
	_xline = xline;
	_trcens = trcens;
}

Elastic_Trace_Header::~Elastic_Trace_Header()
{
}

