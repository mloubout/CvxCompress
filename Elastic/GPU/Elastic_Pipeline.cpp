#include <stdio.h>
#include <stdlib.h>
#include "Elastic_Buffer.hxx"
#include "Elastic_Pipeline.hxx"
#include "Elastic_Propagator.hxx"
#include "Elastic_Shot.hxx"
#include "Elastic_Modeling_Job.hxx"
#include "gpuAssert.h"

Elastic_Pipeline::Elastic_Pipeline(
		int log_level,
		Elastic_Propagator* prop,
		int pipe_id,
		int pipe_y0,
		int pipe_y1,
		int pipe_z0,
		int pipe_z1
		)
{
	_log_level = log_level;
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

	_optseq_launch = 0L;
	_optseq_data_in = 0L;
	_optseq_data_out = 0L;
	_num_optseq_launch = 0;
	_num_optseq_data_in = 0;
	_num_optseq_data_out = 0;

	_d_Mem = 0L;

	_d_RxLoc = 0L;
	_d_RxLoc_block = 0L;
	_h_RxLoc_block_Offset = 0L;
	_h_RxLoc_num_blocks = 0L;

	_d_RxRes = 0L;

	_h_RxRes_curr = 0L;
	_h_RxRes_curr_block_offset = 0L;
	_h_RxRes_curr_timestep = 0L;
	_h_RxRes_curr_num_rx = 0L;
	_h_RxRes_curr_flags = 0L;
	_h_RxRes_curr_num_blocks = 0L;
	
	_h_RxRes_prev = 0L;
	_h_RxRes_prev_block_offset = 0L;
	_h_RxRes_prev_timestep = 0L;
	_h_RxRes_prev_num_rx = 0L;
	_h_RxRes_prev_flags = 0L;
	_h_RxRes_prev_num_blocks = 0L;
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

void Elastic_Pipeline::Reset()
{
	for (int i = 0;  i < _num_buffers;  ++i)
	{
		_buffers[i]->Reset();
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

void Elastic_Pipeline::_Shift_Rx_Blocks()
{
	//printf("Elastic_Pipeline::_Shift_Rx_Blocks - start\n");
	for (int iDev = 0;  iDev < _num_devices;  ++iDev)
	{
		int nn = _h_RxLoc_num_blocks[iDev];
		void* first_d_RxLoc_block = _d_RxLoc_block[iDev][0];
		for (int ib = 1;  ib < nn;  ++ib)
		{
			_d_RxLoc_block[iDev][ib-1] = _d_RxLoc_block[iDev][ib];
			_h_RxLoc_block_Offset[iDev][ib-1] = _h_RxLoc_block_Offset[iDev][ib];
		}
		_d_RxLoc_block[iDev][nn-1] = first_d_RxLoc_block;
		_h_RxLoc_block_Offset[iDev][nn-1] = -1;
	}

	// flip pinned receiver results transfer blocks
	void** tmp = _h_RxRes_curr;
	_h_RxRes_curr = _h_RxRes_prev;
	_h_RxRes_prev = tmp;

	// shift meta data
	if (_h_RxRes_prev_num_blocks != 0L)
	{
		for (int iDev = 0;  iDev < _num_devices;  ++iDev)
	        {
			delete [] _h_RxRes_prev_block_offset[iDev];
			delete [] _h_RxRes_prev_timestep[iDev];
			delete [] _h_RxRes_prev_num_rx[iDev];
			delete [] _h_RxRes_prev_flags[iDev];
		}
		delete [] _h_RxRes_prev_block_offset;
		delete [] _h_RxRes_prev_timestep;
		delete [] _h_RxRes_prev_num_rx;
		delete [] _h_RxRes_prev_flags;
		delete [] _h_RxRes_prev_num_blocks;
	}

	_h_RxRes_prev_block_offset = _h_RxRes_curr_block_offset;
	_h_RxRes_prev_timestep = _h_RxRes_curr_timestep;
	_h_RxRes_prev_num_rx = _h_RxRes_curr_num_rx;
	_h_RxRes_prev_flags = _h_RxRes_curr_flags;
	_h_RxRes_prev_num_blocks = _h_RxRes_curr_num_blocks;

	_h_RxRes_curr_block_offset = 0L;
	_h_RxRes_curr_timestep = 0L;
	_h_RxRes_curr_num_rx = 0L;
	_h_RxRes_curr_flags = 0L;
	_h_RxRes_curr_num_blocks = 0L;
	//printf("Elastic_Pipeline::_Shift_Rx_Blocks - end\n");
}

void Elastic_Pipeline::Launch_Receiver_Extraction_Kernels(Elastic_Shot* shot)
{
#pragma omp parallel for
	for (int iDev = 0;  iDev < _num_devices;  ++iDev)
	{
		int device_id = _device_IDs[iDev];
		int num_full_compute = _h_RxRes_curr_num_blocks[iDev];
		shot->Extract_Receiver_Values_From_Device(
				_prop, this, device_id,
				_h_RxRes_curr_block_offset[iDev],_h_RxRes_curr_timestep[iDev],_h_RxRes_curr_num_rx[iDev],_h_RxRes_curr_flags[iDev],num_full_compute,
				(float**)(_d_RxLoc_block[iDev]),
				(float*)(_d_RxRes[iDev]),
				(float*)(_h_RxRes_curr[iDev]));
	}
}

void Elastic_Pipeline::Launch_Receiver_Data_Transfers(Elastic_Shot* shot)
{
	_Shift_Rx_Blocks();

	_h_RxRes_curr_block_offset = new int*[_num_devices];
	_h_RxRes_curr_timestep = new int*[_num_devices];
	_h_RxRes_curr_num_rx = new int*[_num_devices];
	_h_RxRes_curr_flags = new int*[_num_devices];
	_h_RxRes_curr_num_blocks = new int[_num_devices];

#pragma omp parallel for
	for (int iDev = 0;  iDev < _num_devices;  ++iDev)
	{
		int device_id = _device_IDs[iDev];

		int num_full_compute = _h_RxLoc_num_blocks[iDev];
		if (num_full_compute > 0)
		{
			// find first full compute buffer in pipeline
			Elastic_Buffer* first_full_compute_buffer = 0L;
			for (int iBuf = 0;  iBuf < _num_buffers;  ++iBuf)
			{
				Elastic_Buffer* buffer = _buffers[iBuf];
				if (buffer->Get_Device_ID() == device_id && buffer->Is_Compute() && !buffer->Is_Partial_Compute())
				{
					if (!buffer->Get_M1_Buffer()->Is_Compute() || (buffer->Get_M1_Buffer()->Is_Compute() && buffer->Get_M1_Buffer()->Is_Partial_Compute()))
					{
						first_full_compute_buffer = buffer;
						break;
					}
				}
			}
			//char bufname[4096];
			//printf("First full compute buffer is %s\n",first_full_compute_buffer->Get_Name_String(bufname));
			_h_RxLoc_block_Offset[iDev][num_full_compute-1] = first_full_compute_buffer->Get_Block_Offset(1,0);

			_h_RxRes_curr_block_offset[iDev] = new int[num_full_compute];
			_h_RxRes_curr_timestep[iDev] = new int[num_full_compute];
			_h_RxRes_curr_num_rx[iDev] = new int[num_full_compute];
			_h_RxRes_curr_flags[iDev] = new int[num_full_compute];
			_h_RxRes_curr_num_blocks[iDev] = num_full_compute;
	
			for (int ib = 0;  ib < num_full_compute;  ++ib)
			{
				// defaults
				_h_RxRes_curr_block_offset[iDev][ib] = -1;
				_h_RxRes_curr_timestep[iDev][ib] = -1;
				_h_RxRes_curr_num_rx[iDev][ib] = 0;
				_h_RxRes_curr_flags[iDev][ib] = 0;

				// match compute buffers with RxLoc buffers
				int rxloc_block_offset = _h_RxLoc_block_Offset[iDev][ib];
				if (rxloc_block_offset >= 0)
				{
					for (int iBuf = 0;  iBuf < _num_buffers;  ++iBuf)
					{
						Elastic_Buffer* buffer = _buffers[iBuf];
						int cmp_buffer_block_offset = buffer->Get_Block_Offset(1,0);
						if (buffer->Get_Device_ID() == device_id && buffer->Is_Compute() && !buffer->Is_Partial_Compute() && rxloc_block_offset == cmp_buffer_block_offset)
						{
							// found it
							_h_RxRes_curr_block_offset[iDev][ib] = rxloc_block_offset;
							_h_RxRes_curr_timestep[iDev][ib] = buffer->Get_Block_Timestep(1,0);
							break;
						}
					}
				}
			}

			//printf("Launch_Extraction_Kernels :: num_full_compute=%d\n",num_full_compute);
			//for (int ib = 0;  ib < num_full_compute;  ++ib)
			//{
			//	printf("  Block %d :: Timestep=%d, #RX=%d, Flags=%d\n",_h_RxRes_curr_block_offset[iDev][ib],_h_RxRes_curr_timestep[iDev][ib],_h_RxRes_curr_num_rx[iDev][ib],_h_RxRes_curr_flags[iDev][ib]);
			//}

			shot->Start_Extract_Receiver_Values_From_Device(
				_prop, this, device_id,
				_h_RxRes_curr_block_offset[iDev],_h_RxRes_curr_timestep[iDev],_h_RxRes_curr_num_rx[iDev],_h_RxRes_curr_flags[iDev],num_full_compute,
				(float**)(_d_RxLoc_block[iDev]),
				(float*)(_d_RxRes[iDev]),
				(float*)(_h_RxRes_curr[iDev]));
		}
		else
		{
			_h_RxRes_curr_block_offset[iDev] = 0L;
			_h_RxRes_curr_timestep[iDev] = 0L;
			_h_RxRes_curr_num_rx[iDev] = 0L;
			_h_RxRes_curr_flags[iDev] = 0L;
			_h_RxRes_curr_num_blocks[iDev] = 0;
		}
	}
}

void Elastic_Pipeline::DEMUX_Receiver_Values_For_One_Device(Elastic_Shot* shot, int device_index)
{
	if (device_index >= 0 && device_index < _num_devices)
	{
		if (_h_RxRes_prev_num_blocks != 0L && _h_RxRes_prev_num_blocks[device_index] > 0)
		{
			shot->DEMUX_Receiver_Values(
					this,
					_h_RxRes_prev_block_offset[device_index],
					_h_RxRes_prev_timestep[device_index],
					_h_RxRes_prev_num_rx[device_index],
					_h_RxRes_prev_flags[device_index],
					_h_RxRes_prev_num_blocks[device_index],
					(float*)(_h_RxRes_prev[device_index]));
		}
	}
	else
	{
		printf("Error! Elastic_Pipeline::DEMUX_Receiver_Values_For_One_Device - Illegal device index, was %d, should be >= 0 and < %d\n",device_index,_num_devices);
		exit(-1);
	}
}

void Elastic_Pipeline::DEMUX_Receiver_Values(Elastic_Shot* shot)
{
#pragma omp parallel for
	for (int iDev = 0;  iDev < _num_devices;  ++iDev)
	{
		DEMUX_Receiver_Values_For_One_Device(shot, iDev);
	}
}

void Elastic_Pipeline::Launch_Data_Transfers()
{
	/*
	for (int i = 0;  i < _num_optseq_data_in;  ++i)
	{
		_buffers[_optseq_data_in[i]]->Launch_Data_Transfers();
	}
	for (int i = 0;  i < _num_optseq_data_out;  ++i)
	{
		_buffers[_optseq_data_out[i]]->Launch_Data_Transfers();
	}
	*/
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

void Elastic_Pipeline::Launch_Compute_Kernel(float dti, Elastic_Shot* shot, int num_z)
{
	for (int i = 0;  i < _num_optseq_launch;  ++i)
	{
		_buffers[_optseq_launch[i]]->Launch_Compute_Kernel(false, dti, shot, num_z);
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

void Elastic_Pipeline::Compute_RX_Device_Memory_Requirement(int device_id, int& rxloc_size, int& rxres_size)
{
	rxloc_size = 0;
	rxres_size = 0;
	Elastic_Modeling_Job* job = _prop->Get_Job();
	for (int idx = 0;  idx < job->Get_Number_Of_Shots();  ++idx)
	{
		Elastic_Shot* shot = job->Get_Shot_By_Index(idx);
		int curr_rxloc_size, curr_rxres_size, num_full_compute;
		shot->Calculate_RX_Locations_And_Results_Size(_prop, device_id, curr_rxloc_size, curr_rxres_size, num_full_compute);
		if (curr_rxloc_size > rxloc_size) rxloc_size = curr_rxloc_size;
		if (curr_rxres_size > rxres_size) rxres_size = curr_rxres_size;
	}
}

bool Elastic_Pipeline::Verify_All_Devices_Have_Enough_Memory()
{
	int* device_ids = Get_All_Device_IDs();
	//printf("_num_devices = %d\n",_num_devices);
 	for (int i = 0;   i < _num_devices;  ++i)
	{
		int device_id = device_ids[i];
		unsigned long reqd_memory = Compute_Device_Memory_Requirement(device_id);
		int rxloc_size, rxres_size;
		Compute_RX_Device_Memory_Requirement(device_id,rxloc_size,rxres_size);
		reqd_memory += (unsigned long)rxloc_size;
		reqd_memory += (unsigned long)rxres_size;
		cudaSetDevice(device_id);
		size_t free, total;
		cudaError_t err = cudaMemGetInfo(&free,&total);
		//printf("device %d :: free %.2f MB, reqd %.2f MB\n",device_id,(double)free/1048576.0,(double)reqd_memory/1048576.0);
                if (err == cudaSuccess)
                {
			if (free < reqd_memory)
			{
				printf("Device %d is short %ld bytes.\n",device_id,reqd_memory-free);
				return false;
			}
		}
		else
		{
			printf("Unable to check free memory for device %d.\n",device_id);
			return false;
		}
	}
	return true;
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
	unsigned long reqd_mem = Compute_Device_Memory_Requirement(device_id);
	int rxloc_size, rxres_size;
	Compute_RX_Device_Memory_Requirement(device_id, rxloc_size, rxres_size);
	unsigned long tot_reqd_mem = reqd_mem + (unsigned long)(rxloc_size + rxres_size);
	double size_MB = (double)tot_reqd_mem / (double)(1024.0 * 1024.0);
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
	delete [] _optseq_launch;
	delete [] _optseq_data_in;
	delete [] _optseq_data_out;
	_optseq_launch = 0L;
	_optseq_data_in = 0L;
	_optseq_data_out = 0L;
	_num_optseq_launch = 0;
	_num_optseq_data_in = 0;
	_num_optseq_data_out = 0;
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
	if (_h_RxRes_curr != 0L)
	{
		for (int i = 0;  i < _num_devices;  ++i)
                {
			if (_h_RxRes_curr[i] != 0L) cudaFreeHost(_h_RxRes_curr[i]);
		}
		delete [] _h_RxRes_curr;
		_h_RxRes_curr = 0L;
	}
	if (_h_RxRes_prev != 0L)
	{
		for (int i = 0;  i < _num_devices;  ++i)
                {
			if (_h_RxRes_prev[i] != 0L) cudaFreeHost(_h_RxRes_prev[i]);
		}
		delete [] _h_RxRes_prev;
		_h_RxRes_prev = 0L;
	}
	if (_d_RxLoc != 0L)
	{
		delete [] _d_RxLoc;
		_d_RxLoc = 0L;
	}
	if (_d_RxRes != 0L)
	{
		delete [] _d_RxRes;
		_d_RxRes = 0L;
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

void Elastic_Pipeline::Allocate_RxLoc_Buffer(Elastic_Shot* shot)
{
	//printf("Elastic_Pipeline::Allocate_RxLoc_Buffer - start\n");
	if (_num_devices > 0)
	{
		_d_RxLoc_block = new void**[_num_devices];
                _h_RxLoc_block_Offset = new int*[_num_devices];
                _h_RxLoc_num_blocks = new int[_num_devices];
		for (int i = 0;  i < _num_devices;  ++i)
		{
			int device_id = _device_IDs[i];
			//printf("_d_RxLoc[%d] = %p\n",device_id,_d_RxLoc[i]);
			int rxloc_size, rxres_size, num_full_compute;
			shot->Calculate_RX_Locations_And_Results_Size(_prop, device_id, rxloc_size, rxres_size, num_full_compute);
			_h_RxLoc_num_blocks[i] = num_full_compute;
			if (num_full_compute > 0)
			{
				_d_RxLoc_block[i] = new void*[num_full_compute];
				_d_RxLoc_block[i][0] = _d_RxLoc[i];
				//printf("_d_RxLoc_block[%d][%d] = %p\n",device_id,0,_d_RxLoc_block[i][0]);
				_h_RxLoc_block_Offset[i] = new int[num_full_compute];
				_h_RxLoc_block_Offset[i][0] = -1;
				for (int ib = 1;  ib < num_full_compute;  ++ib)
				{
					_d_RxLoc_block[i][ib] = ((char*)_d_RxLoc_block[i][ib-1]) + (rxloc_size / num_full_compute);
					_h_RxLoc_block_Offset[i][ib] = -1;
					//printf("_d_RxLoc_block[%d][%d] = %p\n",device_id,ib,_d_RxLoc_block[i][ib]);
				}
			}
			else
			{
				_d_RxLoc_block[i] = 0L;
				_h_RxLoc_block_Offset[i] = 0L;
			}
		}
	}
	//printf("Elastic_Pipeline::Allocate_RxLoc_Buffer - end\n");
}

void Elastic_Pipeline::Free_RxLoc_Buffer(Elastic_Shot* shot)
{
	if (_d_RxLoc_block != 0L)
	{
		for (int i = 0;  i < _num_devices;  ++i)
                {
			delete [] _d_RxLoc_block[i];
		}
		delete [] _d_RxLoc_block;
		delete [] _h_RxLoc_block_Offset;
		delete [] _h_RxLoc_num_blocks;
		_d_RxLoc_block = 0L;
		_h_RxLoc_block_Offset = 0L;
		_h_RxLoc_num_blocks = 0L;
	}
}

bool Elastic_Pipeline::Allocate_Device_Memory()
{
	bool success = false;
	Free_Device_Memory();
	_Compile_Device_IDs();
	if (_num_devices > 0)
	{
		success = true;
		_d_Mem = new void*[_num_devices];
		_d_RxLoc = new void*[_num_devices];
		_d_RxRes = new void*[_num_devices];
		_h_RxRes_curr = new void*[_num_devices];
		_h_RxRes_prev = new void*[_num_devices];
		for (int i = 0;  i < _num_devices;  ++i)
		{
			int device_id = _device_IDs[i];
			unsigned long reqd_mem = Compute_Device_Memory_Requirement(device_id);
			int rxloc_size, rxres_size;
			Compute_RX_Device_Memory_Requirement(device_id, rxloc_size, rxres_size);
			unsigned long tot_reqd_mem = reqd_mem + (unsigned long)(rxloc_size + rxres_size);
			double tot_reqd_mem_MB = (double)tot_reqd_mem / 1048576.0;
			cudaSetDevice(device_id);
			cudaError_t err = cudaMalloc(&(_d_Mem[i]),tot_reqd_mem);
			if (err == cudaSuccess)
			{
				cudaMemset(_d_Mem[i], 0, tot_reqd_mem);  // zero block, always do this
				if (_log_level >= 4) printf("cudaMalloc (device %d) :: ALLOCATED %.2f MB device memory\n",device_id,tot_reqd_mem_MB);
				unsigned long offset = 0;
				for (int j = 0;  j < _num_buffers;  ++j)
				{
					if (_buffers[j]->Get_Device_ID() == device_id)
					{
						offset = _buffers[j]->Allocate_Device_Blocks(_d_Mem[i],offset);
					}
				}
				_d_RxLoc[i] = ((char*)_d_Mem[i]) + reqd_mem;
				_d_RxRes[i] = ((char*)_d_RxLoc[i]) + rxloc_size;
				gpuErrchk( cudaHostAlloc(&(_h_RxRes_curr[i]), rxres_size, cudaHostAllocDefault) );
				gpuErrchk( cudaHostAlloc(&(_h_RxRes_prev[i]), rxres_size, cudaHostAllocDefault) );
			}
			else
			{
				success = false;
				_d_Mem[i] = 0L;
				_d_RxLoc[i] = 0L;
				_d_RxRes[i] = 0L;
				_h_RxRes_curr[i] = 0L;
				_h_RxRes_prev[i] = 0L;
				if (_log_level >= 4) printf("cudaMalloc (device %d) :: FAILED TO ALLOCATE %.2f MB device memory\n",device_id,tot_reqd_mem_MB);
			}
		}

		if (success)
		{
			for (int i = 0;  i < _num_buffers;  ++i)
			{
				_buffers[i]->Enable_Peer_Access();
			}

			// determine optimal iteration sequences for the CUDA kernel launches and data transfers
			// ..go wide before deep
			int* idx_launch = new int[_num_devices];
			int* idx_data_in = new int[_num_devices];
			int* idx_data_out = new int[_num_devices];
			for (int i = 0;  i < _num_devices;  ++i)
			{
				idx_launch[i] = 0;
				idx_data_in[i] = 0;
				idx_data_out[i] = 0;
			}
			_optseq_launch = new int[_num_buffers];
			_optseq_data_in = new int[_num_buffers];
			_optseq_data_out = new int[_num_buffers];
			_num_optseq_launch = 0;
			_num_optseq_data_in = 0;
			_num_optseq_data_out = 0;
			bool not_done;
			do
			{
				not_done = false;
				for (int iDev = 0;  iDev < _num_devices;  ++iDev)
				{
					int device_id = _device_IDs[iDev];
					// kernel launches
					for (int i = idx_launch[iDev];  i < _num_buffers;  ++i)
					{
						if (_buffers[i]->Get_Device_ID() == device_id && _buffers[i]->Is_Compute())
						{
							_optseq_launch[_num_optseq_launch++] = i;
							idx_launch[iDev] = i+1;
							//char namestr[4096];
							//printf("OPTIMAL LAUNCH :: %s (device %d)\n",_buffers[i]->Get_Name_String(namestr),device_id);
							not_done = true;
							break;
						}
					}
					// data in
					for (int i = idx_data_in[iDev];  i < _num_buffers;  ++i)
					{
						if (_buffers[i]->Get_Device_ID() == device_id && _buffers[i]->Is_Input())
						{
							_optseq_data_in[_num_optseq_data_in++] = i;
							idx_data_in[iDev] = i+1;
							not_done = true;
							break;
						}
					}
					// data out
					for (int i = idx_data_out[iDev];  i < _num_buffers;  ++i)
					{
						if (_buffers[i]->Get_Device_ID() == device_id && _buffers[i]->Is_Device2Host())
						{
							_optseq_data_out[_num_optseq_data_out++] = i;
							idx_data_out[iDev] = i+1;
							not_done = true;
							break;
						}
					}
				}
			} while (not_done);
			delete [] idx_data_out;
			delete [] idx_data_in;
			delete [] idx_launch;
		}
		else
		{
			_optseq_launch = 0L;
			_optseq_data_in = 0L;
			_optseq_data_out = 0L;
			_num_optseq_launch = 0;
			_num_optseq_data_in = 0;
			_num_optseq_data_out = 0;
		}
	}
	return success;
}

