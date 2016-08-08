#ifndef CVX_ELASTIC_TRACE_HEADER_HXX
#define CVX_ELASTIC_TRACE_HEADER_HXX

#include <time.h>

//
// This class is a container for trace meta-data normally found in the SEGY trace header.
//

class Elastic_Trace_Header
{
public:
	Elastic_Trace_Header(
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
		int trcens,
		int rec_ffid,
		time_t acqtime,
		int usec
		);
	virtual ~Elastic_Trace_Header();

	inline double Get_Start_Time() {return _start_time;}
	inline double Get_Sample_Rate() {return _sample_rate;}
	inline int Get_NSAMP() {return _nsamp;}

	inline int Get_File_Number() {return _iFile;}

	inline int Get_Flags() {return _flags;}
	inline bool Is_Pressure() {return (_flags & 1) != 0 ? true : false;}
	inline bool Is_Vx() {return (_flags & 2) != 0 ? true : false;}
	inline bool Is_Vy() {return (_flags & 4) != 0 ? true : false;}
	inline bool Is_Vz() {return (_flags & 8) != 0 ? true : false;}
	inline bool Is_Velocity() {return Is_Vx() | Is_Vy() | Is_Vz();}

	inline double Get_Location_X() {return _loc_x;}
	inline double Get_Location_Y() {return _loc_y;}
	inline double Get_Location_Z() {return _loc_z;}
	
	inline int Get_Inline() {return _iline;}
	inline int Get_Crossline() {return _xline;}
	inline int Get_Trace_Ensemble() {return _trcens;}

	inline int Get_Receiver_FFID() {return _rec_ffid;}
	inline time_t Get_Shot_Time() {return _acqtime;}
	inline int Get_Shot_Time_usec() {return _usec;}

protected:
	double _start_time;
	double _sample_rate;
	int _nsamp;
	int _iFile;
	int _flags;
	double _loc_x;
	double _loc_y;
	double _loc_z;
	int _iline;
	int _xline;
	int _trcens;
	int _rec_ffid;
	time_t _acqtime;
	int _usec;
};

#endif
