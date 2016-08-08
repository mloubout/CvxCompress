#ifndef CVX_ELAORTHMODEL_HXX
#define CVX_ELAORTHMODEL_HXX

#include <map>
#include <string>

class GPU_Runtime_Parameters
{
	public:
		GPU_Runtime_Parameters(int num_pipes, int steps_per_GPU, int num_devices, const int* device_Ids)
		{
			_num_devices = num_devices;
			_device_Ids = new int[_num_devices];
			for (int i = 0;  i < _num_devices;  ++i) _device_Ids[i] = device_Ids[i];
			_num_pipes = num_pipes;
			_steps_per_GPU = steps_per_GPU;
		}
		virtual ~GPU_Runtime_Parameters()
		{
			delete [] _device_Ids;
		}

		int Get_Number_Of_GPU_Pipes() {return _num_pipes;}
		int Get_Steps_Per_GPU() {return _steps_per_GPU;}
		const int* Get_GPU_Devices() {return _device_Ids;}
		int Get_Number_Of_GPU_Devices() {return _num_devices;}

	private:
		int _num_pipes;
		int _steps_per_GPU;
		int _num_devices;
		int* _device_Ids;
};
class Voxet;
class Voxet_Memory_Mapper;
class Elastic_Modeling_Job;
class Variable_Water_Velocity;
class ElaOrthModel
{
	public:
		ElaOrthModel(
			char* parmString
		);

		virtual ~ElaOrthModel();

		void Get_Trace_Arguments(
				int& samples_per_trace,
				double& sample_rate_in_seconds
				);
		Voxet* Get_Voxet();

		void runWorker();

		int runShot(
			int	nTraces,
			float	soux,
			float	souy,
			float	souz,
			float*	recx,
			float*	recy,
			float*	recz,
			float*	samples
			);

	private:
		int _nsamp;
		char* _parmString;
		Voxet_Memory_Mapper* _mapper;
		Variable_Water_Velocity* _Vwxyzt;
		Elastic_Modeling_Job* _job;
		Elastic_Modeling_Job* _Create_Job_Object(float soux, float souy, float souz);
		std::map<std::string,double> _min_prop_val;
		std::map<std::string,double> _max_prop_val;
		std::map<std::string,GPU_Runtime_Parameters*> _gpu_device_params;
};
#endif
