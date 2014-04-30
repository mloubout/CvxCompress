#include <stdio.h>
#include "Elastic_Buffer.hxx"
#include "Elastic_Pipeline.hxx"
#include "Elastic_Propagator.hxx"
#include "Elastic_Shot.hxx"

Elastic_Pipeline::Elastic_Pipeline(
		Elastic_Propagator* prop,
		int pipe_id,
		int pipe_y0,
		int pipe_y1,
		int pipe_z0,
		int pipe_z1
		)
{
	_prop = prop;
	_pipe_id = pipe_id;
	_pipe_y0 = pipe_y0;
	_pipe_y1 = pipe_y1;
	_pipe_z0 = pipe_z0;
	_pipe_z1 = pipe_z1;
	_num_buffers = 0;
	_buffers = new Elastic_Buffer*[1000];
	_num_devices = 0;
	_device_IDs = 0L;
	_d_Mem = 0L;
}

Elastic_Pipeline::~Elastic_Pipeline()
{
	Free_Device_Memory();
	if (_buffers != 0L)
	{
		for (int i = 0;  i < _num_buffers;  ++i)
		{
			if (_buffers[i] != 0L)
			{
				delete _buffers[i];
				_buffers[i] = 0L;
			}
		}
		delete [] _buffers;
	}
}

void Elastic_Pipeline::Append_Buffer(Elastic_Buffer* new_buffer)
{
	_buffers[_num_buffers] = new_buffer;
	++_num_buffers;
}

void Elastic_Pipeline::Add_EM_Buffer(Elastic_Buffer* new_buffer)
{
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		if (_buffers[i]->Is_Compute() && _buffers[i]->Get_Device_ID() == new_buffer->Get_Device_ID())
		{
			_buffers[i]->Add_EM_Buffer(new_buffer);
		}
	}
	Append_Buffer(new_buffer);
}

void Elastic_Pipeline::Shift_Buffers()
{
	for (int i = 0;  i < _num_buffers;  ++i)
        {
                _buffers[i]->Shift_Buffer();
	}
}

void Elastic_Pipeline::Launch_Data_Transfers()
{
	for (int i = 0;  i < _num_buffers;  ++i)
        {
                _buffers[i]->Launch_Data_Transfers();
	}
}

void Elastic_Pipeline::Launch_Simple_Copy_Kernel()
{
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		_buffers[i]->Launch_Simple_Copy_Kernel();
	}
}

void Elastic_Pipeline::Launch_Compute_Kernel(float dti, Elastic_Shot* shot)
{
	for (int i = 0;  i < _num_buffers;  ++i)
        {
                _buffers[i]->Launch_Compute_Kernel(false, dti, shot);
        }
}

int Elastic_Pipeline::Get_Number_Of_Buffers()
{
	return _num_buffers;
}

Elastic_Buffer* Elastic_Pipeline::Get_Buffer(int index)
{
	if (index >= 0 && index < _num_buffers)
	{
		return _buffers[index];
	}
	else if (index < 0 && (_num_buffers+index) >= 0)
	{
		return _buffers[_num_buffers+index];
	}
	else
	{
		return 0L;
	}
}

int Elastic_Pipeline::Get_Device_Count()
{
	int count = 1, curr_device_id = _buffers[0]->Get_Device_ID();
	for (int i = 1;  i < _num_buffers;  ++i)
	{
		if (_buffers[i]->Get_Device_ID() != curr_device_id)
		{
			++count;
			curr_device_id = _buffers[i]->Get_Device_ID();
		}
	}
	return count;
}

int Elastic_Pipeline::Get_Total_Number_Of_Timesteps()
{
	int max_timesteps = -1000000000;
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		if (_buffers[i]->Is_Compute() && _buffers[i]->Get_Relative_Timestep() > max_timesteps)
		{
			max_timesteps = _buffers[i]->Get_Relative_Timestep();
		}
	}
	return max_timesteps;
}

int Elastic_Pipeline::Get_Total_Block_Offset()
{
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		if (_buffers[i]->Is_Device2Host())
		{
			return _buffers[i]->Get_Max_Relative_Block_Offset();
		}
	}
	return 0;
}

int Elastic_Pipeline::Get_Input_Block_Offset(int iteration)
{
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		if (_buffers[i]->Is_Host2Device())
		{
			if (_buffers[i]->Get_Block_Timestep(0,iteration) >= 0)
			{
				return _buffers[i]->Get_Block_Offset(0,iteration);
			}
		}
	}
	return -1;
}

int Elastic_Pipeline::Get_Output_Block_Offset(int iteration)
{
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		if (_buffers[i]->Is_Device2Host())
		{
			if (_buffers[i]->Get_Block_Timestep(-1,iteration) > 0)
			{
				return _buffers[i]->Get_Block_Offset(-1,iteration);
			}
		}
	}
	return -1;
}

int Elastic_Pipeline::Get_Input_Block_Timestep(int iteration)
{
	for (int i = 0;  i < _num_buffers;  ++i)
        {
                if (_buffers[i]->Is_Host2Device())
                {
			return _buffers[i]->Get_Block_Timestep(0,iteration);
		}
	}
	return -1;
}

int Elastic_Pipeline::Get_Output_Block_Timestep(int iteration)
{
	for (int i = 0;  i < _num_buffers;  ++i)
        {
                if (_buffers[i]->Is_Device2Host())
                {
			return _buffers[i]->Get_Block_Timestep(-1,iteration);
		}
	}
	return -1;
}

double Elastic_Pipeline::Get_Workload()
{
	double acc = 0.0;
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		acc += _buffers[i]->Get_Workload();
	}
	return acc;
}

double Elastic_Pipeline::Get_Minimum_Workload()
{
	return _prop->Get_Minimum_Workload();
}

double Elastic_Pipeline::Get_Workload(int device_id)
{
	double acc = 0.0;
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		if (_buffers[i]->Get_Device_ID() == device_id)
		{
			acc += _buffers[i]->Get_Workload();
		}
	}
	return acc;
}

double Elastic_Pipeline::Get_Computational_Overhead(int device_id)
{
	return (double)Get_Device_Count() * Get_Workload(device_id) / Get_Minimum_Workload();
}

unsigned long Elastic_Pipeline::Compute_Device_Memory_Requirement(int device_id)
{
	unsigned long acc = 0;
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		if (_buffers[i]->Get_Device_ID() == device_id)
		{
			acc += _buffers[i]->Compute_Device_Memory_Requirement();
		}
	}
	return acc;
}

void Elastic_Pipeline::Print_Graphical(int device_id)
{
	printf("device_id = %d\n",device_id);
	// find block offset range
	int min_block_offset = 1000000000, max_block_offset = -1000000000;
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		if (device_id == _buffers[i]->Get_Device_ID())
		{
			if (_buffers[i]->Get_Min_Relative_Block_Offset() < min_block_offset) min_block_offset = _buffers[i]->Get_Min_Relative_Block_Offset();
			if (_buffers[i]->Get_Max_Relative_Block_Offset() > max_block_offset) max_block_offset = _buffers[i]->Get_Max_Relative_Block_Offset();
		}
	}

	// print bX
	printf("  bX  ");
	for (int i = min_block_offset;  i <= max_block_offset;  ++i) printf("%3d ",i);
	double size_MB = (double)Compute_Device_Memory_Requirement(device_id) / (double)(1024.0 * 1024.0);
	printf(" | %.0fMB, Workload=%.2f%%\n",size_MB,100.0*Get_Computational_Overhead(device_id));

	// print buffers
	char buf[256];
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		if (device_id == _buffers[i]->Get_Device_ID())
		{
			printf("%4s ",_buffers[i]->Get_Name_String(buf));
			for (int j = min_block_offset;  j <= max_block_offset;  ++j)
			{
				if (_buffers[i]->Block_Is_Included_By_Relative_Offset(j))
				{
					if (_buffers[i]->Block_Is_Host2Device_By_Relative_Offset(j))
					{
						printf(" H2D");
					}
					else if (_buffers[i]->Block_Is_Input_By_Relative_Offset(j))
					{
						printf("   I");
					}
					else if (_buffers[i]->Block_Is_Partial_Compute_By_Relative_Offset(j))
					{
						printf("  PC");
					}
					else if (_buffers[i]->Block_Is_Compute_By_Relative_Offset(j))
					{
						printf("   C");
					}
					else if (_buffers[i]->Block_Is_Device2Host_By_Relative_Offset(j))
					{
						printf(" D2H");
					}
					else if (Block_Is_Output_By_Relative_Offset(_buffers[i], j))
					{
						printf("   O");
					}
					else
					{
						printf("   X");
					}
				}
				else
				{
					printf("    ");
				}
			}
			printf("  | y=[%4d,%4d]",_buffers[i]->Get_Y0(),_buffers[i]->Get_Y1());
			if (_buffers[i]->Is_Input()) 
				printf("  inp_y=[%4d,%4d]",_buffers[i]->Get_Inp_Y0(),_buffers[i]->Get_Inp_Y1());
			else
				printf("                   ");
			if (_buffers[i]->Is_Compute()) 
				printf("  cmp_y=[%4d,%4d]",_buffers[i]->Get_Cmp_Y0(),_buffers[i]->Get_Cmp_Y1());
			else
				printf("                   ");
			printf("\n");
		}
	}
}

void Elastic_Pipeline::Print_Graphical()
{
	// this routine and others assume buffers with same device_id are stored consecutively
	printf("\nPipe %d :: y=[%4d,%4d]\n",_pipe_id+1,_pipe_y0,_pipe_y1);
	int curr_device_id = -1;
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		int device_id = _buffers[i]->Get_Device_ID();
		if (device_id != curr_device_id)
		{
			curr_device_id = device_id;
			Print_Graphical(device_id);
			printf("\n");
		}
	}
}

bool Elastic_Pipeline::Block_Is_Output_By_Relative_Offset(Elastic_Buffer* buffer, int relative_block_offset)
{
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		if (_buffers[i]->Is_Input() && _buffers[i]->Get_Source_Buffer() == buffer && _buffers[i]->Get_Source_Block_Relative_Offset() == relative_block_offset)
		{
			return true;
		}
	}
	return false;
}

void Elastic_Pipeline::Free_Device_Memory()
{
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		_buffers[i]->Free_Device_Blocks();
	}
	if (_d_Mem != 0L)
        {
                for (int i = 0;  i < _num_devices;  ++i)
                {
                        cudaSetDevice(_device_IDs[i]);
                        if (_d_Mem[i] != 0L) cudaFree(_d_Mem[i]);
                }
                delete [] _d_Mem;
                _d_Mem = 0L;
        }
}

int* Elastic_Pipeline::Get_All_Device_IDs()
{
	if (_device_IDs == 0L)
	{
		_Compile_Device_IDs();
	}
	return _device_IDs;
}

void Elastic_Pipeline::_Compile_Device_IDs()
{
	if (_device_IDs != 0L)
	{
		delete [] _device_IDs;
		_device_IDs = 0L;
	}
	_num_devices = 0;
	if (_num_buffers > 0)
	{
		_device_IDs = new int[_num_buffers];
		for (int i = 0;  i < _num_buffers;  ++i)
		{
			int device_id = _buffers[i]->Get_Device_ID();
			bool found = false;
			for (int j = 0;  j < _num_devices && !found;  ++j)
			{
				if (_device_IDs[j] == device_id)
				{
					found = true;
				}
			}
			if (!found)
			{
				_device_IDs[_num_devices] = device_id;
				++_num_devices;
			}
		}
	}
}

void Elastic_Pipeline::Allocate_Device_Memory()
{
	Free_Device_Memory();
	_Compile_Device_IDs();
	if (_num_devices > 0)
	{
		_d_Mem = new void*[_num_devices];
		for (int i = 0;  i < _num_devices;  ++i)
		{
			int device_id = _device_IDs[i];
			unsigned long reqd_mem = Compute_Device_Memory_Requirement(device_id);
			double reqd_mem_MB = (double)reqd_mem / 1048576.0;
			cudaSetDevice(device_id);
			cudaError_t err = cudaMalloc(&(_d_Mem[i]),reqd_mem);
			if (err == cudaSuccess)
			{
				cudaMemset(_d_Mem[i], 0, reqd_mem);  // zero block, always do this
				printf("cudaMalloc (device %d) :: ALLOCATED %.2f MB device memory\n",device_id,reqd_mem_MB);
				unsigned long offset = 0;
				for (int j = 0;  j < _num_buffers;  ++j)
				{
					if (_buffers[j]->Get_Device_ID() == device_id)
					{
						offset = _buffers[j]->Allocate_Device_Blocks(_d_Mem[i],offset);
					}
				}
			}
			else
			{
				_d_Mem[i] = 0L;
				printf("cudaMalloc (device %d) :: FAILED TO ALLOCATE %.2f MB device memory\n",device_id,reqd_mem_MB);
			}
		}

		for (int i = 0;  i < _num_buffers;  ++i)
		{
			_buffers[i]->Enable_Peer_Access();
		}
	}	
}

