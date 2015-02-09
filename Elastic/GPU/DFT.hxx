#ifndef CVX_SEISMOD_ELASTIC_DFT_HXX
#define CVX_SEISMOD_ELASTIC_DFT_HXX

void cmplx_DFT(
	double* real_in,
	double* imag_in,
	double* real_out,
	double* imag_out,
	int len
	);

void Apply_Butterworth_Low_Pass_Filter(
	int log_level,
	double* signal_in,
	double* signal_out,
	int len,
	double dt,  // in seconds
	double f_cut,
	int filter_order
	);

double Butterworth_Low_Pass_Filter_Find_Fcut_From_Fmax(
	int log_level,
	double f_max,
	int filter_order
	);

void Compute_Time_Integrated_Source_Wavelet(
	int log_level,
	double* stf,
	double* stf_int,
	int len,
	double dt
	);

#endif

