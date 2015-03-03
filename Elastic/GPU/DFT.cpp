#include <math.h>
#include <stdio.h>

//
// Real-to-complex discrete fourier transform implementation.
// Yes, I know there are many, many alternative implementations.
// This one is neither fast nor brilliant, but speed is not an issue
// for this part of the code and there are benefits to not having
// yet another 3rd party library dependency.
//

// in and out cannot overlap.
void cmplx_DFT(
	double* real_in,
	double* imag_in,
	double* real_out,
	double* imag_out,
	int len
	)
{
	for (int k = 0;  k < len;  ++k)
	{
		double real = 0.0;
		double imag = 0.0;
		for (int n = 0;  n < len;  ++n)
		{
			double arg = -6.283185307179586476925286766559 * (double)(k * n) / (double)len;
			double euler_real, euler_imag;
			sincos(arg, &euler_imag, &euler_real);
			//double euler_real = cos(arg);
			//double euler_imag = sin(arg);
			real += (real_in[n] * euler_real - imag_in[n] * euler_imag);
			imag += (real_in[n] * euler_imag + imag_in[n] * euler_real);
		}
		real_out[k] = real;
		imag_out[k] = imag;
	}
}

double Compute_Butterworth_Gain_At_Frequency(
	double f,
	double f_cut,
	int filter_order
	)
{
	return sqrt(1.0 / (1.0 + pow(f/f_cut,2.0*(double)filter_order)));
}

void Compute_Time_Integrated_Source_Wavelet(
	int log_level,
	double* stf,
	double* stf_int,
	int len,
	double dt
	)
{
	double dc = stf[0];
	for (int i = 1;  i < len;  ++i) dc += stf[i];
	dc /= (double)len;
	stf_int[0] = (stf[0] - dc) * dt;
	for (int i = 1;  i < len;  ++i) stf_int[i] = stf_int[i-1] + (stf[i] - dc) * dt;
}

void Apply_Butterworth_Low_Pass_Filter(
	int log_level,
	double* signal_in,
	double* signal_out,
	int len,
	double dt,  // in seconds
	double f_cut,
	int filter_order
	)
{
	int filt_len = (len + 1) & (~1);  // make it lowest even number larger than or equal to len
	double* real1 = new double[filt_len];
	double* imag1 = new double[filt_len];
	double* real2 = new double[filt_len];
	double* imag2 = new double[filt_len];
	for (int i = 0;  i < filt_len;  ++i)
	{
		real1[i] = i < len ? signal_in[i] : 0.0;
		imag1[i] = 0.0;
	}
	cmplx_DFT(real1,imag1,real2,imag2,filt_len);
	int nk = (filt_len/2)+1;
	double f_ny = 0.5 / dt;
	for (int k = 0;  k < nk;  ++k)
	{
		double freq = f_ny * (double)k / (double)(nk-1);
		double bw_gain = Compute_Butterworth_Gain_At_Frequency(freq,f_cut,filter_order) / (double)filt_len;
		if (log_level >= 5) printf("bw_gain is %.6lf at f=%.2f\n",bw_gain*(double)filt_len,freq);
		real2[k] *= bw_gain;
		imag2[k] *= -bw_gain;
		if (k > 0)
		{
			real2[filt_len-k] *= bw_gain;
			imag2[filt_len-k] *= -bw_gain;
		}
	}
	cmplx_DFT(real2,imag2,real1,imag1,filt_len);
	for (int i = 0;  i < len;  ++i) signal_out[i] = real1[i];
	delete [] imag2;
	delete [] real2;
	delete [] imag1;
	delete [] real1;
}

double Butterworth_Low_Pass_Filter_Find_Fcut_From_Fmax(
	int log_level,
	double f_max,
	int filter_order
	)
{
	// TMJ 06/19/14
	// binary search for f_max.
	// there is an analytical solution for sure, but this numerical approach is the path of least resistance for me.
	double f0 = 0.0;
	double f1 = f_max;
	double f_cut;
	bool done = false;
	do
	{
		f_cut = (f0 + f1) / 2.0;
		double bw_att = -20.0 * log10(Compute_Butterworth_Gain_At_Frequency(f_max,f_cut,filter_order));
		if (log_level >= 5) printf("f0=%lf, f1=%lf, f_cut=%lf, bw_att=%lf\n",f0,f1,f_cut,bw_att);
		if (bw_att < 39.99)
		{
			// signal is too strong
			f1 = f_cut;
		}
		else if (bw_att > 40.01)
		{
			// signal is too weak
			f0 = f_cut;
		}
		else
		{
			done = true;
		}
	} while (!done);
	return f_cut;
}

/*
int main(int argc, char* arv[])
{
	int filter_order = 7;
	double f_cut = Butterworth_Low_Pass_Filter_Find_Fcut_From_Fmax(80.0f,filter_order);
	printf("f_max=80Hz requires f_cut=%.2fHz for filter order %d\n",f_cut,filter_order);

	double real1[8];
	double imag1[8];
	double real2[8];
	double imag2[8];
	
	for (int i = 0;  i < 8;  ++i) {real1[i] = imag1[i] = real2[i] = imag2[i] = 0.0f;}
	for (int i = 0;  i < 8;  ++i)
	{
		double arg = 6.283185307179586476925286766559f * (double)i / 8.0f;
		real1[i] = sinf(arg);
	}
	cmplx_DFT(real1,imag1,real2,imag2,8);
	for (int i = 0;  i < 8;  ++i)
	{
		printf("DFT1[%d] = %f + i * %f\n",i,real2[i],imag2[i]);
	}
	for (int i = 0;  i < 8;  ++i)
	{
		real2[i] = 
		imag2[i] = -imag2[i];
	}
	cmplx_DFT(real2,imag2,real1,imag1,8);
	for (int i = 0;  i < 8;  ++i)
        {
                printf("DFT2[%d] = %f + i * %f\n",i,real1[i],imag1[i]);
        }
}
*/

