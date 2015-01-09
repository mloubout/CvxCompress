#include <math.h>
#include <stdio.h>
#include "Elastic_Trace.hxx"
#include "Elastic_Resampled_Trace.hxx"

Elastic_Resampled_Trace::Elastic_Resampled_Trace(
		double input_start_time,
		double input_sample_rate,
		double output_start_time,
		double output_sample_rate,
		int nsamp,
		int sinc_length
		)
: Elastic_Trace(output_start_time,output_sample_rate,nsamp)
{
	_sinc_length = sinc_length;
	_sinc_neg_idx = -((_sinc_length - 1) / 2);
	_sinc_pos_idx = _sinc_neg_idx + _sinc_length - 1;
	_input_start_time = input_start_time;
	_input_sample_rate = input_sample_rate;
	_rcp_sample_rate = 1.0 / output_sample_rate;
	_rcp_input_sample_rate = 1.0 / _input_sample_rate;
}

Elastic_Resampled_Trace::~Elastic_Resampled_Trace()
{
}

void Elastic_Resampled_Trace::Correct_Amplitudes(float& min_before, float& max_before, float& min_after, float& max_after)
{
	if (_sinc_length > 2)
	{
		min_before =  1e37f;
		max_before = -1e37f;
		min_after  =  1e37f;
		max_after  = -1e37f;
		for (int out_index = 0;  out_index < _nsamp;  ++out_index)
		{
			if (_samples[out_index] < min_before) min_before = _samples[out_index];
			if (_samples[out_index] > max_before) max_before = _samples[out_index];

			// convert sample output index to time
			double out_time = _start_time + (double)out_index * _sample_rate;
			// center sample input index (i.e. index of input sample that is the "center" in the sinc).
			int center_input_idx = (int)trunc((out_time - _input_start_time) * _rcp_input_sample_rate);
			int idx_0 = center_input_idx + _sinc_neg_idx;
			if (idx_0 >= 0)
			{
				// loop over all input samples covered by sinc.
				int idx_1 = center_input_idx + _sinc_pos_idx;
				float weight = 0.0f;
				for (int idx = idx_0;  idx <= idx_1;  ++idx)
				{
					float coeff = Compute_Sinc_Coefficient(idx,out_time);
					weight += coeff;
				}
				_samples[out_index] /= weight;
			}

			if (_samples[out_index] < min_after) min_after = _samples[out_index];
			if (_samples[out_index] > max_after) max_after = _samples[out_index];
		}
	}
	else
	{
		min_before = 0.0f;
		max_before = 0.0f;
		min_after  = 0.0f;
		max_after  = 0.0f;
	}
}

