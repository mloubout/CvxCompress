#include <string.h>
#include "Elastic_Trace.hxx"

Elastic_Trace::Elastic_Trace(
		double start_time,
		double sample_rate,
		int nsamp
		)
{
	_start_time = start_time;
	_sample_rate = sample_rate;
	_nsamp = nsamp;
	_samples = new float[_nsamp];
	memset((void*)_samples, 0, _nsamp*sizeof(float));
}

Elastic_Trace::~Elastic_Trace()
{
	delete [] _samples;
}

float Elastic_Trace::Get_Sample(int index)
{
	if (index >= 0 && index < _nsamp)
	{
		return _samples[index];
	}
	return 0.0f;
}

void Elastic_Trace::Set_Sample(int index, float new_value)
{
	if (index >= 0 && index < _nsamp)
        {
		_samples[index] = new_value;
	}
}

void Elastic_Trace::Add_To_Sample(int index, float delta)
{
	if (index >= 0 && index < _nsamp)
        {
		_samples[index] += delta;
	}
}

