#ifndef CVX_ELASTIC_TRACE_HXX
#define CVX_ELASTIC_TRACE_HXX

class Elastic_Trace
{
public:
	Elastic_Trace(
		double start_time,
		double sample_rate,
		int nsamp
		);
	virtual ~Elastic_Trace();

	inline int Get_NSAMP() {return _nsamp;}
	inline float* Get_Samples() {return _samples;}

	float Get_Sample(int index);
	void Set_Sample(int index, float new_value);
	void Add_To_Sample(int index, float delta);

protected:
	double _start_time;
	double _sample_rate;

	int _nsamp;
	float* _samples;
};

#endif

