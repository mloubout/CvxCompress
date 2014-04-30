#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "Elastic_SEGY_File.hxx"
#include "Elastic_SEGY_File_Receiver_Range.hxx"
#include "Elastic_Buffer.hxx"
#include "Elastic_Pipeline.hxx"
#include "Elastic_Propagator.hxx"
#include <cuda_runtime_api.h>
#include "gpuAssert.h"

Elastic_SEGY_File::Elastic_SEGY_File(
		int fileidx,
		const char* base_filename,
		double sample_rate,
		double tshift,
		double reclen,
		bool do_P,
		bool do_Vx,
		bool do_Vy,
		bool do_Vz
		)
{
	_Is_Valid = false;
	_fileidx = fileidx;
	_base_filename = strdup(base_filename);
	_sample_rate = sample_rate;
	_tshift = tshift;
	_reclen = reclen;
	_do_P = do_P;
	_do_Vx = do_Vx;
	_do_Vy = do_Vy;
	_do_Vz = do_Vz;
	_rcv_ranges = 0L;
	_num_rcv_ranges = 0;
	_nWF = 0;
	if (_do_P) ++_nWF;
	if (_do_Vx) ++_nWF;
	if (_do_Vy) ++_nWF;
	if (_do_Vz) ++_nWF;
	_totSteps = 0;
	_h_transfer = 0L;
	_h_rx = 0L;
	_nBlks = 0;
	_num_pipes = 0;
	_rcv_x = 0L;
	_rcv_y = 0L;
	_rcv_z = 0L;
	_rcv_n = 0L;
	_rcv_i = 0L;
	_rcv_nn = 0L;
	_max_rx = 0L;
	_tot_rx = 0L;
	_h_prev_ts = 0L;
	_h_curr_ts = 0L;
	_device2pipe = 0L;
	_d_rcv_loc = 0L;
	_d_rcv_x = 0L;
	_d_rcv_y = 0L;
	_d_rcv_z = 0L;
	_d_transfer = 0L;
	_max_device_id = -1;
}

Elastic_SEGY_File::~Elastic_SEGY_File()
{
	if (_base_filename != 0L) free(_base_filename);
	_Destroy_Buffers();
}

void Elastic_SEGY_File::Add_Receiver_Range_X(
		int range_idx,
		double start,
		double end,
		double interval
		)
{
	_Get_Receiver_Range(range_idx)->Add_X(start,end,interval);
}

void Elastic_SEGY_File::Add_Receiver_Range_Y(
		int range_idx,
		double start,
		double end,
		double interval
		)
{
	_Get_Receiver_Range(range_idx)->Add_Y(start,end,interval);
}

void Elastic_SEGY_File::Add_Receiver_Range_Z(
		int range_idx,
		double start,
		double end,
		double interval
		)
{
	_Get_Receiver_Range(range_idx)->Add_Z(start,end,interval);
}

void Elastic_SEGY_File::_Destroy_Buffers()
{
	if (_nBlks > 0 && _num_pipes > 0)
	{
		for (int iBlk = 0;  iBlk < _nBlks;  ++iBlk)
		{
			delete [] _h_rx[iBlk];
			for (int iPipe = 0;  iPipe < _num_pipes;  ++iPipe)
			{
				if (_rcv_n[iBlk][iPipe] > 0)
				{
					delete [] _rcv_x[iBlk][iPipe];
					delete [] _rcv_y[iBlk][iPipe];
					delete [] _rcv_z[iBlk][iPipe];
				}
			}
			delete [] _rcv_x[iBlk];
			delete [] _rcv_y[iBlk];
			delete [] _rcv_z[iBlk];
			delete [] _rcv_n[iBlk];
			delete [] _rcv_i[iBlk];
		}
		delete [] _device2pipe;
		delete [] _h_prev_ts;
		delete [] _h_curr_ts;
		delete [] _h_rx;
		delete [] _rcv_x;
		delete [] _rcv_y;
		delete [] _rcv_z;
		delete [] _rcv_n;
		delete [] _rcv_i;
		delete [] _rcv_nn;
		delete [] _max_rx;
		delete [] _tot_rx;
	}
	_totSteps = 0;
	_device2pipe = 0L;
	_h_prev_ts = 0L;
	_h_curr_ts = 0L;
	if (_h_transfer != 0L) gpuErrchk( cudaFreeHost(_h_transfer) );
	_h_transfer = 0L;
	_h_rx = 0L;
	_nBlks = 0;
	_num_pipes = 0;
	_rcv_x = 0L;
	_rcv_y = 0L;
	_rcv_z = 0L;
	_rcv_n = 0L;
	_rcv_i = 0L;
	_rcv_nn = 0L;
	_max_rx = 0L;
	_tot_rx = 0L;
	if (_max_device_id >= 0)
	{
		for (int device_id = 0;  device_id <= _max_device_id;  ++device_id)
		{
			if (_d_rcv_loc[device_id] != 0L)
			{
				cudaSetDevice(device_id);
				gpuErrchk( cudaFree(_d_rcv_loc[device_id]) );
				_d_rcv_loc[device_id] = 0L;
			}
			if (_d_rcv_x[device_id] != 0L) {delete [] _d_rcv_x[device_id];  _d_rcv_x[device_id]=0L;}
			if (_d_rcv_y[device_id] != 0L) {delete [] _d_rcv_y[device_id];  _d_rcv_y[device_id]=0L;}
			if (_d_rcv_z[device_id] != 0L) {delete [] _d_rcv_z[device_id];  _d_rcv_z[device_id]=0L;}
			if (_d_transfer[device_id] != 0L)
			{
				for (int iStep = 0;  iStep < _totSteps;  ++iStep)
				{
					for (int iWF = 0;  iWF < 4;  ++iWF)
					{
						delete [] _d_transfer[device_id][iStep][iWF];
					}
					delete [] _d_transfer[device_id][iStep];
				}
				delete [] _d_transfer[device_id];
				_d_transfer[device_id] = 0L;
			}
		}
		delete [] _d_rcv_loc;
		delete [] _d_rcv_x;
		delete [] _d_rcv_y;
		delete [] _d_rcv_z;
		delete [] _d_transfer;
	}
	_max_device_id = -1;
}

void Elastic_SEGY_File::Shift_Receiver_Transfer_Buffers()
{
	for (int iBlk = 0;  iBlk < _nBlks;  ++iBlk)
	{
		_h_prev_ts[iBlk] = _h_curr_ts[iBlk];
		// leave _h_curr_ts unchanged
	}
}

bool Elastic_SEGY_File::Create_New_Device_To_Host_Transfer(
	int device_id,
	int block_number,
	int timestep,
	float*& dst_buf,
	int& dst_size
	)
{
	if (block_number >= 0 && block_number < _nBlks)
        {
		int pipe_id = _device2pipe[device_id];
		dst_buf = _h_rx[block_number][timestep%_totSteps] + _rcv_i[block_number][pipe_id];
		dst_size = _rcv_n[block_number][pipe_id];
		_h_curr_ts[block_number] = timestep;
		return true;
	}
	return false;
}

void Elastic_SEGY_File::Create_Receiver_Transfer_Buffers(
	Elastic_Propagator* prop
	)
{
	// delete old buffer structures
	_Destroy_Buffers();

	// create new buffer structures
	float *rcv_x, *rcv_y, *rcv_z;
	int num_rx = Compute_Receiver_Locations(rcv_x, rcv_y, rcv_z);
	if (num_rx > 0)
	{
		// sort into blocks and pipes
		_totSteps = prop->Get_Total_Number_Of_Timesteps();
		_nBlks = prop->Get_Number_Of_Blocks();
		_num_pipes = prop->Get_Number_Of_Pipelines();
		_rcv_x = new float**[_nBlks];
		_rcv_y = new float**[_nBlks];
		_rcv_z = new float**[_nBlks];
		_rcv_n = new int*[_nBlks];
		_rcv_i = new int*[_nBlks];
		_rcv_nn = new int[_nBlks];
		_trc_idx = new int**[_nBlks];
		_max_rx = new int[_num_pipes];
		_tot_rx = new int[_num_pipes];
		_max_device_id = 0;
		for (int iPipe = 0;  iPipe < _num_pipes;  ++iPipe)
		{
			_max_rx[iPipe] = 0;
			_tot_rx[iPipe] = 0;
			Elastic_Pipeline* pipe = prop->Get_Pipeline(iPipe);
			for (int iDev = 0;  iDev < pipe->Get_Device_Count();  ++iDev) if (pipe->Get_All_Device_IDs()[iDev] > _max_device_id) _max_device_id = pipe->Get_All_Device_IDs()[iDev];
		}
		_device2pipe = new int[_max_device_id+1];
		for (int i = 0;  i <= _max_device_id;  ++i) _device2pipe[i] = -1;
		for (int iPipe = 0;  iPipe < _num_pipes;  ++iPipe)
                {
			Elastic_Pipeline* pipe = prop->Get_Pipeline(iPipe);
			for (int iDev = 0;  iDev < pipe->Get_Device_Count();  ++iDev)
			{
				int device_id = pipe->Get_All_Device_IDs()[iDev];
				_device2pipe[device_id] = iPipe;
			}
		}
		for (int iBlk = 0;  iBlk < _nBlks;  ++iBlk)
		{
			_rcv_nn[iBlk] = 0;
			_rcv_x[iBlk] = new float*[_num_pipes];
			_rcv_y[iBlk] = new float*[_num_pipes];
			_rcv_z[iBlk] = new float*[_num_pipes];
			_rcv_n[iBlk] = new int[_num_pipes];
			_rcv_i[iBlk] = new int[_num_pipes];
			_trc_idx[iBlk] = new int*[_num_pipes];
			for (int iPipe = 0, rxIdx = 0;  iPipe < _num_pipes;  ++iPipe)
			{
				Elastic_Pipeline* pipe = prop->Get_Pipeline(iPipe);

				// count number of receivers that fall into this block
				int x0 = iBlk * prop->Get_Block_Size_X();
				int x1 = x0 + prop->Get_Block_Size_X() - 1;
				int y0 = pipe->Get_Y0();
				int y1 = pipe->Get_Y1();
				
				int rxCnt = 0;
				for (int iRx = 0;  iRx < num_rx;  ++iRx)
				{
					int ix = (int)round(rcv_x[iRx]) + 1;
					int iy = (int)round(rcv_y[iRx]) + 1;
					int iz = (int)round(rcv_z[iRx]) + 1;
					if (ix >= x0 && ix <= x1 && iy >= y0 && iy <= y1)
					{
						++rxCnt;
					}
				}
				if (rxCnt > 0)
				{
					_rcv_nn[iBlk] += rxCnt;
					_rcv_n[iBlk][iPipe] = rxCnt;
					_rcv_i[iBlk][iPipe] = rxIdx;
					rxIdx += rxCnt;
					_rcv_x[iBlk][iPipe] = new float[rxCnt];
					_rcv_y[iBlk][iPipe] = new float[rxCnt];
					_rcv_z[iBlk][iPipe] = new float[rxCnt];
					_trc_idx[iBlk][iPipe] = new int[rxCnt];
					if (rxCnt > _max_rx[iPipe]) _max_rx[iPipe] = rxCnt;
					_tot_rx[iPipe] += rxCnt;
				}
				else
				{
					_rcv_n[iBlk][iPipe] = 0;
					_rcv_i[iBlk][iPipe] = -1;
					_rcv_x[iBlk][iPipe] = 0L;
					_rcv_y[iBlk][iPipe] = 0L;
					_rcv_z[iBlk][iPipe] = 0L;
					_trc_idx[iBlk][iPipe] = 0L;
				}
				rxCnt = 0;
				for (int iRx = 0;  iRx < num_rx;  ++iRx)
                                {
                                        int ix = (int)round(rcv_x[iRx]) + 1;
                                        int iy = (int)round(rcv_y[iRx]) + 1;
                                        int iz = (int)round(rcv_z[iRx]) + 1;
                                        if (ix >= x0 && ix <= x1 && iy >= y0 && iy <= y1)
                                        {
						_rcv_x[iBlk][iPipe][rxCnt] = rcv_x[iRx];
						_rcv_y[iBlk][iPipe][rxCnt] = rcv_y[iRx];
						_rcv_z[iBlk][iPipe][rxCnt] = rcv_z[iRx];
						_trc_idx[iBlk][iPipe][rxCnt] = iRx;
						++rxCnt;
					}
				}
			}
		}

		delete [] rcv_x;
		delete [] rcv_y;
		delete [] rcv_z;

		// create pinned input transfer buffer
		int HostTransferBufSize = _nWF * num_rx * prop->Get_Total_Number_Of_Timesteps() * sizeof(float);
		gpuErrchk( cudaHostAlloc((void**)&_h_transfer, HostTransferBufSize, cudaHostAllocDefault) );
		_h_rx = new float**[_nBlks];
		_h_prev_ts = new int[_nBlks];
		_h_curr_ts = new int[_nBlks];
		for (int iBlk = 0, hidx = 0;  iBlk < _nBlks;  ++iBlk)
		{
			_h_rx[iBlk] = new float*[prop->Get_Total_Number_Of_Timesteps()];
			for (int iStep = 0;  iStep < prop->Get_Total_Number_Of_Timesteps();  ++iStep)
			{
				_h_rx[iBlk][iStep] = _h_transfer + hidx;
				hidx += _rcv_nn[iBlk];
			}
		}

		// create work buffers for each device
		_d_rcv_loc = new float*[_max_device_id+1];
		_d_rcv_x = new float**[_max_device_id+1];
		_d_rcv_y = new float**[_max_device_id+1];
		_d_rcv_z = new float**[_max_device_id+1];
		_d_transfer = new float****[_max_device_id+1];  // _d_transfers[device_id][timestep][0->P,1->Vx,2->Vy,3->Vz][0->comp,1->transfer_out]
		for (int i = 0;  i <= _max_device_id;  ++i)
		{
			_d_rcv_x[i] = 0L;
			_d_rcv_y[i] = 0L;
			_d_rcv_z[i] = 0L;
			_d_rcv_loc[i] = 0L;
			_d_transfer[i] = 0L;
		}
		for (int iPipe = 0;  iPipe < _num_pipes;  ++iPipe)
                {
			// create tempoary host buffer and pack receiver locations into it.
			// this buffer is copied to each device.
			int RxLocPipeBufSize = _tot_rx[iPipe] * 3 * sizeof(float);
			float* h_rcv_locs = 0L;
			gpuErrchk( cudaHostAlloc((void**)&h_rcv_locs, RxLocPipeBufSize, cudaHostAllocDefault) );
			for (int iBlk = 0, idx = 0;  iBlk < _nBlks;  ++iBlk)
			{
				for (int iRx = 0;  iRx < _rcv_n[iBlk][iPipe];  ++iRx) h_rcv_locs[idx+iRx] = _rcv_x[iBlk][iPipe][iRx];		idx+=_rcv_n[iBlk][iPipe];
				for (int iRx = 0;  iRx < _rcv_n[iBlk][iPipe];  ++iRx) h_rcv_locs[idx+iRx] = _rcv_y[iBlk][iPipe][iRx];           idx+=_rcv_n[iBlk][iPipe];
				for (int iRx = 0;  iRx < _rcv_n[iBlk][iPipe];  ++iRx) h_rcv_locs[idx+iRx] = _rcv_z[iBlk][iPipe][iRx];           idx+=_rcv_n[iBlk][iPipe];
			}

			// allocate buffer on device side, copy over flattened rx locations, create ptr structure on host side
			Elastic_Pipeline* pipe = prop->Get_Pipeline(iPipe);
			for (int iDev = 0;  iDev < pipe->Get_Device_Count();  ++iDev)
			{
				int device_id = pipe->Get_All_Device_IDs()[iDev];
				if (_d_transfer[device_id] == 0L)
				{
					_d_transfer[device_id] = new float***[_totSteps+1];
					for (int iStep = 0;  iStep <= _totSteps;  ++iStep)
					{
						_d_transfer[device_id][iStep] = new float**[4];
						for (int iWF = 0;  iWF < 4;  ++iWF)
						{
							_d_transfer[device_id][iStep][iWF] = new float*[2];
							for (int ii = 0;  ii < 2;  ++ii)
							{
								_d_transfer[device_id][iStep][iWF][ii] = 0L;
							}
						}
					}
				}
				
				// figure out how much device memory is needed for rx transfer buffer(s)
				int RxTransferBufSize = 0;
				for (int iBuf = 0;  iBuf < pipe->Get_Number_Of_Buffers();  ++iBuf)
				{
					Elastic_Buffer* buffer = pipe->Get_Buffer(iBuf);
					if (buffer->Get_Device_ID() == device_id && buffer->Is_Compute() && !buffer->Is_Partial_Compute())
					{
						if (buffer->Is_Particle_Velocity())
						{
							if (_do_P) ++RxTransferBufSize;
						}
						else
						{
							if (_do_Vx) ++RxTransferBufSize;
							if (_do_Vy) ++RxTransferBufSize;
							if (_do_Vz) ++RxTransferBufSize;
						}
					}
				}
				RxTransferBufSize = RxTransferBufSize * 2 * _max_rx[iPipe] * sizeof(float);

				cudaSetDevice(device_id);
				int RxTotBufSize = RxLocPipeBufSize + RxTransferBufSize;
				gpuErrchk( cudaMalloc((void**)&(_d_rcv_loc[device_id]), RxTotBufSize) );
				gpuErrchk( cudaMemcpy((void*)_d_rcv_loc[device_id], (const void*)h_rcv_locs, RxLocPipeBufSize, cudaMemcpyHostToDevice) );
				gpuErrchk( cudaMemset((void*)((char*)(_d_rcv_loc[device_id])+RxLocPipeBufSize), 0, RxTransferBufSize) );

				_d_rcv_x[device_id] = new float*[_nBlks];
				_d_rcv_y[device_id] = new float*[_nBlks];
				_d_rcv_z[device_id] = new float*[_nBlks];
				for (int iBlk = 0, idx = 0;  iBlk < _nBlks;  ++iBlk)
				{
					_d_rcv_x[device_id][iBlk] = _d_rcv_loc[device_id] + idx;		idx+=_rcv_n[iBlk][iPipe];
					_d_rcv_y[device_id][iBlk] = _d_rcv_loc[device_id] + idx;		idx+=_rcv_n[iBlk][iPipe];
					_d_rcv_z[device_id][iBlk] = _d_rcv_loc[device_id] + idx;		idx+=_rcv_n[iBlk][iPipe];
				}

				float* d_rx_out = (float*)((char*)(_d_rcv_loc[device_id]) + RxLocPipeBufSize);
				for (int iBuf = 0, idx = 0;  iBuf < pipe->Get_Number_Of_Buffers();  ++iBuf)
                                {
                                        Elastic_Buffer* buffer = pipe->Get_Buffer(iBuf);
                                        if (buffer->Get_Device_ID() == device_id && buffer->Is_Compute() && !buffer->Is_Partial_Compute())
                                        {
                                                if (buffer->Is_Particle_Velocity())
                                                {
                                                        if (_do_P)
							{
								_d_transfer[device_id][buffer->Get_Relative_Timestep()][0][0] = d_rx_out + idx;	idx+=_max_rx[iPipe];
								_d_transfer[device_id][buffer->Get_Relative_Timestep()][0][1] = d_rx_out + idx;	idx+=_max_rx[iPipe];
							}
                                                }
                                                else
                                                {
                                                        if (_do_Vx)
							{
								_d_transfer[device_id][buffer->Get_Relative_Timestep()][1][0] = d_rx_out + idx;	idx+=_max_rx[iPipe];
								_d_transfer[device_id][buffer->Get_Relative_Timestep()][1][1] = d_rx_out + idx;	idx+=_max_rx[iPipe];
							}
                                                        if (_do_Vy)
							{
								_d_transfer[device_id][buffer->Get_Relative_Timestep()][2][0] = d_rx_out + idx;	idx+=_max_rx[iPipe];
								_d_transfer[device_id][buffer->Get_Relative_Timestep()][2][1] = d_rx_out + idx;	idx+=_max_rx[iPipe];
							}
                                                        if (_do_Vz)
							{
								_d_transfer[device_id][buffer->Get_Relative_Timestep()][3][0] = d_rx_out + idx;	idx+=_max_rx[iPipe];
								_d_transfer[device_id][buffer->Get_Relative_Timestep()][3][1] = d_rx_out + idx;	idx+=_max_rx[iPipe];
							}
                                                }
                                        }
                                }
			}

			// release pinned memory
			gpuErrchk( cudaFreeHost((void*)h_rcv_locs) );
		}
	}
}

int Elastic_SEGY_File::Compute_Receiver_Locations(
		float*& rcv_x,
		float*& rcv_y,
		float*& rcv_z
		)
{
	int num = 0;
	rcv_x = 0L;
	rcv_y = 0L;
	rcv_z = 0L;
	for (int i = 0;  i < _num_rcv_ranges;  ++i)
	{
		float *x,*y,*z;
		int nn = _rcv_ranges[i]->Compute_Receiver_Locations(x,y,z);
		if (nn > 0)
		{
			if (num == 0)
			{
				rcv_x = x;
				rcv_y = y;
				rcv_z = z;
				num = nn;
			}
			else
			{
				float* tmp = new float[num+nn];
				for (int j = 0;  j < num;  ++j) tmp[j] = rcv_x[j];
				for (int j = 0;  j < nn;  ++j) tmp[j+num] = x[j];
				delete [] rcv_x;
				delete [] x;
				rcv_x = tmp;

				tmp = new float[num+nn];
				for (int j = 0;  j < num;  ++j) tmp[j] = rcv_y[j];
				for (int j = 0;  j < nn;  ++j) tmp[j+num] = y[j];
				delete [] rcv_y;
				delete [] y;
				rcv_y = tmp;

				tmp = new float[num+nn];
				for (int j = 0;  j < num;  ++j) tmp[j] = rcv_z[j];
				for (int j = 0;  j < nn;  ++j) tmp[j+num] = z[j];
				delete [] rcv_z;
				delete [] z;
				rcv_z = tmp;

				num = num + nn;
			}
		}
	}
	return num;
}

Elastic_SEGY_File_Receiver_Range* Elastic_SEGY_File::_Get_Receiver_Range(int range_idx)
{
	if (_num_rcv_ranges == 0)
	{
		// add 1st range
		_rcv_ranges = new Elastic_SEGY_File_Receiver_Range*[1];
		_rcv_ranges[0] = new Elastic_SEGY_File_Receiver_Range(range_idx);
		_num_rcv_ranges = 1;
	}
	else
	{
		for (int i = 0;  i < _num_rcv_ranges;  ++i)
		{
			if (_rcv_ranges[i]->Get_Range_Idx() == range_idx)
			{
				// return existing range
				return _rcv_ranges[i];
			}
		}
		// add new range
		Elastic_SEGY_File_Receiver_Range** new_ranges = new Elastic_SEGY_File_Receiver_Range*[_num_rcv_ranges+1];
		for (int i = 0;  i < _num_rcv_ranges;  ++i) new_ranges[i] = _rcv_ranges[i];
		new_ranges[_num_rcv_ranges] = new Elastic_SEGY_File_Receiver_Range(range_idx);
		delete [] _rcv_ranges;
		_rcv_ranges = new_ranges;
		++_num_rcv_ranges;
	}
	return _rcv_ranges[_num_rcv_ranges-1];
}

