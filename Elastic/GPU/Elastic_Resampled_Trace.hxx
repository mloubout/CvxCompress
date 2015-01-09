#ifndef CVX_ELASTIC_RESAMPLED_TRACE_HXX
#define CVX_ELASTIC_RESAMPLED_TRACE_HXX

#include <math.h>
#include "Elastic_Trace.hxx"

class Elastic_Resampled_Trace : public Elastic_Trace
{
	public:
		Elastic_Resampled_Trace(
				double input_start_time,
				double input_sample_rate,
				double output_start_time,
				double output_sample_rate,
				int nsamp,
				int sinc_length
				);
		virtual ~Elastic_Resampled_Trace();

		inline void Add_To_Trace(int input_index, float delta)
		{
			double lowest_input_time = _input_start_time + (double)(input_index - _sinc_pos_idx - 1) * _input_sample_rate;
			double highest_input_time = _input_start_time + (double)(input_index - _sinc_neg_idx + 1) * _input_sample_rate;
			int lowest_out_index = (int)floor((lowest_input_time - _start_time) * _rcp_sample_rate);
			lowest_out_index = lowest_out_index > 0 ? lowest_out_index : 0;
			int highest_out_index = (int)ceil((highest_input_time - _start_time) * _rcp_sample_rate);
			highest_out_index = highest_out_index < _nsamp ? highest_out_index : _nsamp - 1;
			for (int out_index = lowest_out_index;  out_index <= highest_out_index;  ++out_index)
			{
				// convert sample output index to time
				double out_time = _start_time + (double)out_index * _sample_rate;
				// center sample input index (i.e. index of input sample that is the "center" in the sinc).
				int center_input_idx = (int)trunc((out_time - _input_start_time) * _rcp_input_sample_rate);
				// loop over all input samples covered by sinc.
				int idx_0 = center_input_idx + _sinc_neg_idx;
				if (_sinc_length <= 2 || idx_0 < 0)
				{
					// linear interpolation when sinc length is two samples or we don't have enough samples yet to do a sinc
					if (input_index >= center_input_idx && input_index <= center_input_idx+1)
					{
						double inp_sample_distance = (_input_start_time + (double)input_index * _input_sample_rate - out_time) * _rcp_input_sample_rate;
						double coeff = 1.0 - fabs(inp_sample_distance);
						_samples[out_index] += (float)coeff * delta;
					}
				}
				else
				{
					// sinc interpolation
					int idx_1 = center_input_idx + _sinc_pos_idx;
					if (input_index >= idx_0 && input_index <= idx_1)
					{
						// voila!
						_samples[out_index] += Compute_Sinc_Coefficient(input_index,out_time) * delta;
					}
				}
			}
		}

		void Correct_Amplitudes(float& min_before, float& max_before, float& min_after, float& max_after);

	protected:
		int _sinc_length;
		int _sinc_neg_idx;
		int _sinc_pos_idx;
		double _input_start_time;
		double _input_sample_rate;

		double _rcp_sample_rate;
		double _rcp_input_sample_rate;

		inline float Compute_Sinc_Coefficient(
				int idx,
				double out_time
				)
		{
			// compute distance from output time
			double inp_sample_distance = (_input_start_time + (double)idx * _input_sample_rate - out_time) / _input_sample_rate;
			double x = inp_sample_distance * 3.1415926535897932384626433832795;
			double coeff = inp_sample_distance == 0.0 ? 1.0 : sin(x) / x;
			return (float)coeff;
		}
};

#endif
