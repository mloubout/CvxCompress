#ifndef CVX_ESDRD_MI_TMJ_ELASTIC_SEGY_FILE_HXX
#define CVX_ESDRD_MI_TMJ_ELASTIC_SEGY_FILE_HXX

#include "Elastic_Interpolation.hxx"
#include "Elastic_Gather_Type.hxx"

class Elastic_SEGY_File_Receiver_Range;
class Elastic_Propagator;

class Elastic_SEGY_File
{
public:
	Elastic_SEGY_File(
		int fileidx,
		const char* base_filename,
		double sample_rate,
		double tshift,
		double reclen,
		bool do_P,
		bool do_Vx,
		bool do_Vy,
		bool do_Vz
		);
	~Elastic_SEGY_File();

	bool Is_Valid() {return _Is_Valid;}

	int Get_File_Index() {return _fileidx;}
	void Set_File_Index(int fileidx) {_fileidx=fileidx;}

	int Get_SeqNo() {return _seqno;}
	void Set_SeqNo(int seqno) {_seqno=seqno;}

	int Get_GunSeq() {return _gunseq;}
	void Set_GunSeq(int gunseq) {_gunseq=gunseq;}

	Elastic_Gather_Type_t Get_Gather_Type() {return _gather_type;}
	void Set_Gather_Type(Elastic_Gather_Type_t gather_type) {_gather_type = gather_type;}

	Elastic_Interpolation_t Get_Interpolation_Method() {return _interpolation_method;}
	void Set_Interpolation_Method(Elastic_Interpolation_t interpolation_method) {_interpolation_method = interpolation_method;}
	
	double Get_Sample_Rate() {return _sample_rate;}
	double Get_Timeshift() {return _tshift;}
	double Get_Record_Length() {return _reclen;}

        int Get_Selection_Flags();
	bool Do_P() {return _do_P;}
	bool Do_Vx() {return _do_Vx;}
	bool Do_Vy() {return _do_Vy;}
	bool Do_Vz() {return _do_Vz;}

	const char* Get_Base_Filename() {return _base_filename;}
	const char* Get_Full_Path(char* buf, int flag);

	void Add_Receiver_Range_X(
		int range_idx,
		double start,
		double end,
		double interval
		);
	void Add_Receiver_Range_Y(
		int range_idx,
		double start,
		double end,
		double interval
		);
	void Add_Receiver_Range_Z(
		int range_idx,
		double start,
		double end,
		double interval
		);

	void Add_Receiver_Array(
		int nrec, 
		double* rec_x,
		double* rec_y,
		double* rec_z,
		int* iline,
		int* xline,
		int* trcens,
		int* rec_ffid,
		time_t* acqtime,
		int* usec
		);

	int Compute_Receiver_Locations(
		double*& rcv_x,
		double*& rcv_y,
		double*& rcv_z
		);
	int Compute_Receiver_Locations(
		double*& rcv_x,
		double*& rcv_y,
		double*& rcv_z,
		int*& iline,
		int*& xline,
		int*& trcens,
		int*& rec_ffid,
		time_t*& acqtime,
		int*& usec
		);
	int Compute_Receiver_Locations_NO_COPY(
		double*& rcv_x,
		double*& rcv_y,
		double*& rcv_z
		);
	int Compute_Receiver_Locations_NO_COPY(
		double*& rcv_x,
		double*& rcv_y,
		double*& rcv_z,
		int*& iline,
		int*& xline,
		int*& trcens,
		int*& rec_ffid,
		time_t*& acqtime,
		int*& usec
		);
	
	void Write_Source_Wavelet_To_SEGY_File(
			bool Is_Vwxyzt,
			double* filtered,
			double* filtered_int,
			double sample_rate,
			int nsamp,
			double srcx,
			double srcy,
			double srcz,
			int src_il,
			int src_xl
			);
	void Write_SEGY_File(
			float** traces,
			char* EBCDIC_Header,
			bool Is_Vwxyzt,
			double srcx,
			double srcy,
			double srcz,
			int src_il,
			int src_xl,
			double* recx,
			double* recy,
			double* recz,
			int* rec_il,
			int* rec_xl,
			int* trcens,
			short* compon,
			int* rec_ffid,
			time_t* acqtime,
			int* usec,
			float src_model_water_depth,
			float src_model_water_Vp,
			float src_bath_z,
			float* rec_model_water_depth,
			float* rec_model_water_Vp,
			float* rec_bath_z,
			int num_traces,
			int nsamp,
			int flag
			);
	void Write_SEGY_File(
			const char* filename,
			double sample_rate,
			Elastic_Gather_Type_t gather_type,
			bool Is_Vwxyzt,
			int file_idx,
			int seqno,
			int gunseq,
			double start_time,
			float** traces,
			char* EBCDIC_Header,
			double srcx,
			double srcy,
			double srcz,
			int src_il,
			int src_xl,
			double* recx,
			double* recy,
			double* recz,
			int* rec_il,
			int* rec_xl,
			int* trcens,
			short* compon,
			int* rec_ffid,
			time_t* acqtime,
			int* usec,
			float src_model_water_depth,
			float src_model_water_Vp,
			float src_bath_z,
			float* rec_model_water_depth,
			float* rec_model_water_Vp,
			float* rec_bath_z,
			int num_traces,
			int nsamp
			);

	void printRec();
private:
	bool _Is_Valid;
	int _fileidx;
	int _seqno;
	int _gunseq;
	char* _base_filename;
	double _sample_rate;
	double _tshift;
	double _reclen;
	bool _do_P;
	bool _do_Vx;
	bool _do_Vy;
	bool _do_Vz;
	Elastic_Interpolation_t _interpolation_method;

	Elastic_Gather_Type_t _gather_type;

	Elastic_SEGY_File_Receiver_Range** _rcv_ranges;
	int _num_rcv_ranges;
	double* _h_user_rcv_x;
	double* _h_user_rcv_y;
	double* _h_user_rcv_z;
	int* _h_user_iline;
	int* _h_user_xline;
	int* _h_user_trcens;
	int* _h_user_rec_ffid;
	time_t* _h_user_acqtime;
	int* _h_user_usec;
	int _num_user_rcv;

	Elastic_SEGY_File_Receiver_Range* _Get_Receiver_Range(int range_idx);
};

#endif

