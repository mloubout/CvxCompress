#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>

#include "Elastic_Buffer.hxx"
#include "Elastic_Propagator.hxx"
#include "Elastic_Pipeline.hxx"
#include "Elastic_Modeling_Job.hxx"
#include "Elastic_Shot.hxx"
#include "Voxet.hxx"
#include "Global_Coordinate_System.hxx"

#include "Propagate_Stresses.hxx"
#include "Propagate_Particle_Velocities.hxx"

// forward declaration of global call that starts cuda copy kernel
void Host_Simple_Copy_Kernel(
	cudaStream_t cmp_stream,
	void* d_dst,
	void* d_src,
	int nx,
	int ny,
	int nz,
	size_t one_y_size
	);

// Create earth model buffer
Elastic_Buffer::Elastic_Buffer(
		Elastic_Propagator* prop,
		Elastic_Pipeline* pipe,
		int device_id,
		int y0,
		int y1,
		int z0,
		int z1,
		int num_blocks,
		int block_offset,
		Elastic_Buffer* src,
		int dst_block_id
		)
{
	_prop = prop;
	_pipe = pipe;
	_blocks = 0L;
	_current_block_offset = 0L;
	_current_block_timestep = 0L;
	_device_id = device_id;
	_timestep = -1;
	_Is_PV = false;
	_y0 = y0;
	_y1 = y1;
	_z0 = z0;
	_z1 = z1;
	_inp_y0 = y0;
	_inp_y1 = y1;
	_cmp_y0 = 0;
	_cmp_y1 = 0;
	_Is_Model = true;
	_Is_Input = true;
	_Is_Device2Host = false;
	_Is_Compute = false;
	_Is_Partial = false;
	_allocated_bytes = 0;
	_failed_to_allocate_bytes = 0;
	_num_blocks = num_blocks;
	_block_offset = block_offset;
	_cmp_block_id = 0;
	_src = src;
	_dst_block_id = dst_block_id;
	_inp_m2 = 0L;
	_inp_m1 = 0L;
}

// Create input buffer
Elastic_Buffer::Elastic_Buffer(
		Elastic_Propagator* prop,
		Elastic_Pipeline* pipe,
		int device_id,
		int timestep,
		bool Is_PV,
		int y0,
		int y1,
		int z0,
		int z1,
		int num_blocks,
		int block_offset,
		Elastic_Buffer* src,
		int dst_block_id
		)
{
	_prop = prop;
	_pipe = pipe;
	_blocks = 0L;
	_current_block_offset = 0L;
	_current_block_timestep = 0L;
	_device_id = device_id;
	_timestep = timestep;
	_Is_PV = Is_PV;
	_y0 = y0;
	_y1 = y1;
	_z0 = z0;
	_z1 = z1;
	_inp_y0 = y0;
	_inp_y1 = y1;
	_cmp_y0 = 0;
	_cmp_y1 = 0;
	_Is_Model = false;
	_Is_Input = true;
	_Is_Device2Host = false;
	_Is_Compute = false;
	_Is_Partial = false;
	_num_blocks = num_blocks;
	_block_offset = block_offset;
	_cmp_block_id = 0;			
	_src = src;
	_dst_block_id = dst_block_id;
	_inp_m2 = 0L;
	_inp_m1 = 0L;
}

// Create compute buffer
Elastic_Buffer::Elastic_Buffer(
		Elastic_Propagator* prop,
		Elastic_Pipeline* pipe,
		int device_id,
		int timestep,
		bool Is_PV,
		int y0,
		int y1,
		int cmp_y0,
		int cmp_y1,
		int z0,
		int z1,
		int num_blocks,
		int block_offset,
		int cmp_block_id,
		Elastic_Buffer* inp_m2,
		Elastic_Buffer* inp_m1,
		Elastic_Buffer* src,
		int dst_block_id
		)
{
	_prop = prop;
	_pipe = pipe;
	_blocks = 0L;
	_current_block_offset = 0L;
	_current_block_timestep = 0L;
	_device_id = device_id;
	_timestep = timestep;
	_Is_PV = Is_PV;
	_y0 = y0;
	_y1 = y1;
	_z0 = z0;
	_z1 = z1;
	_cmp_y0 = cmp_y0;
	_cmp_y1 = cmp_y1;
	_Is_Model = false;
	_Is_Compute = true;
	if (cmp_y0 > y0)
	{
		_Is_Input = true;
		_Is_Partial = true;
		_inp_y0 = y0;
		_inp_y1 = cmp_y0 - 1;
	}
	else
	{
		_Is_Input = false;
		_Is_Partial = false;
		_inp_y0 = 0;
		_inp_y1 = 0;
	}
	_Is_Device2Host = false;
	_num_blocks = num_blocks;
	_block_offset = block_offset;
	_cmp_block_id = cmp_block_id;
	_inp_m2 = inp_m2;
	_inp_m1 = inp_m1;
	_src = src;
	_dst_block_id = dst_block_id;
}

Elastic_Buffer::~Elastic_Buffer()
{
	Free_Device_Blocks();
}

char* Elastic_Buffer::Get_Name_String(char* buf)
{
	if (_Is_Model)
	{
		sprintf(buf, "E_M");
	}
	else
	{
		sprintf(buf, "%2d-%s", _timestep, _Is_PV?"V":"S");
	}
	return buf;
}

int Elastic_Buffer::Get_Device_ID()
{
	return _device_id;
}

void Elastic_Buffer::Add_EM_Buffer(Elastic_Buffer* new_buffer)
{
	_em = new_buffer;
}

bool Elastic_Buffer::Is_Input()
{
	return _Is_Input;
}

Elastic_Buffer* Elastic_Buffer::Get_Source_Buffer()
{
	return _src;
}

// get input buffer once removed
Elastic_Buffer* Elastic_Buffer::Get_M1_Buffer()
{
	return _inp_m1;
}

// get input buffer twice removed
Elastic_Buffer* Elastic_Buffer::Get_M2_Buffer()
{
	return _inp_m2;
}

int Elastic_Buffer::Get_Source_Block_Relative_Offset()
{
	return Get_Relative_Block_Offset(_dst_block_id);
}

int Elastic_Buffer::Get_Min_Relative_Block_Offset()
{
	return _block_offset - _num_blocks + 1;
}

int Elastic_Buffer::Get_Max_Relative_Block_Offset()
{
	return _block_offset;
}

bool Elastic_Buffer::Block_Is_Included_By_Relative_Offset(int relative_block_offset)
{
	if (relative_block_offset >= Get_Min_Relative_Block_Offset() && relative_block_offset <= Get_Max_Relative_Block_Offset())
	{
		return true;
	}
	else
	{
		return false;
	}
}

int Elastic_Buffer::Get_Relative_Block_Offset(int block_id)
{
	return _block_offset - ((block_id < 0) ? (_num_blocks + block_id) : block_id);
}

int Elastic_Buffer::Get_Block_Offset(int block_id, int iteration)
{
	if (block_id >= 0 && block_id < _num_blocks)
	{
		int boff = _current_block_offset[block_id];
		int ts = _current_block_timestep[block_id];
		Advance_Block_Offset(iteration,boff,ts);
		return boff;
	}
	else if (block_id < 0)
	{
		block_id = _num_blocks + block_id;
		if (block_id >= 0 && block_id < _num_blocks)
		{
			int boff = _current_block_offset[block_id];
			int ts = _current_block_timestep[block_id];
			Advance_Block_Offset(iteration,boff,ts);
			return boff;
		}
		else
		{
			return -1;
		}
	}
}

int Elastic_Buffer::Get_Block_Timestep(int block_id, int iteration)
{
	if (block_id >= 0 && block_id < _num_blocks)
        {
		int boff = _current_block_offset[block_id];
		int ts = _current_block_timestep[block_id];
		Advance_Block_Offset(iteration,boff,ts);
                return ts;
        }
        else
        {
                block_id = _num_blocks + block_id;
                if (block_id >= 0 && block_id < _num_blocks)
		{
			int boff = _current_block_offset[block_id];
			int ts = _current_block_timestep[block_id];
			Advance_Block_Offset(iteration,boff,ts);
			return ts;
		}
                else
                {
                        return -1;
                }
        }
}

bool Elastic_Buffer::Block_Is_Valid(int block_id, int iteration)
{
	return Get_Block_Offset(block_id,iteration) >= 0 ? true : false;
}

bool Elastic_Buffer::Block_Is_Input_By_Relative_Offset(int relative_block_offset)
{
	if (_Is_Input && relative_block_offset == _block_offset)
	{
		return true;
	}
	else
	{
		return false;
	}	
}

bool Elastic_Buffer::Block_Is_Compute_By_Relative_Offset(int relative_block_offset)
{
	if (_Is_Compute && relative_block_offset == Get_Relative_Block_Offset(_cmp_block_id))
	{
		return true;
	}
	else
	{
		return false;
	}
}

void* Elastic_Buffer::Get_Block_By_Offset(int block_offset, int iteration)
{
	//char buf[4096];
	for (int i = 0;  i < _num_blocks;  ++i)
	{
		//printf("%s :: block=%d, block_offset=%d - looking for %d\n",Get_Name_String(buf),i,Get_Block_Offset(i,iteration),block_offset);
		if (Get_Block_Offset(i,iteration) == block_offset)
		{
			//printf("%s :: FOUND BLOCK\n",Get_Name_String(buf));
			return _blocks[i];
		}
	}
	//printf("%s :: NO BLOCK FOUND!\n",Get_Name_String(buf));
	//exit(0);
	return 0L;
}

bool Elastic_Buffer::Block_Is_Device2Host_By_Relative_Offset(int relative_block_offset)
{
	return _Is_Device2Host && Get_Relative_Block_Offset(-1) == relative_block_offset;
}

bool Elastic_Buffer::Block_Is_Host2Device_By_Relative_Offset(int relative_block_offset)
{
	return Block_Is_Input_By_Relative_Offset(relative_block_offset) && (_src == 0L);
}

bool Elastic_Buffer::Set_Is_Device2Host(int flag)
{
	_Is_Device2Host = flag;
}

bool Elastic_Buffer::Is_Device2Host()
{
	return _Is_Device2Host;
}

bool Elastic_Buffer::Is_Host2Device()
{
	return _Is_Input && (_src == 0L);
}

bool Elastic_Buffer::Is_Compute()
{
	return _Is_Compute;
}

bool Elastic_Buffer::Block_Is_Partial_Compute_By_Relative_Offset(int relative_block_offset)
{
	return _Is_Partial && Block_Is_Compute_By_Relative_Offset(relative_block_offset);
}

bool Elastic_Buffer::Is_Partial_Compute()
{
	return _Is_Partial;
}

int Elastic_Buffer::Get_Relative_Timestep()
{
	return _timestep;
}

int Elastic_Buffer::Get_Y0()
{
	return _y0;
}

int Elastic_Buffer::Get_Y1()
{
	return _y1;
}

int Elastic_Buffer::Get_Inp_Y0()
{
	return _inp_y0;
}

int Elastic_Buffer::Get_Inp_Y1()
{
	return _inp_y1;
}

int Elastic_Buffer::Get_Cmp_Y0()
{
	return _cmp_y0;
}

int Elastic_Buffer::Get_Cmp_Y1()
{
	return _cmp_y1;
}

int Elastic_Buffer::Get_Z0()
{
	return _z0;
}

int Elastic_Buffer::Get_Z1()
{
	return _z1;
}

double Elastic_Buffer::Get_Workload()
{
	if (_Is_Compute)
	{
		return (double)(_cmp_y1-_cmp_y0+1) * Get_Relative_Cost();
	}
	else
	{
		return 0;
	}
}

unsigned long Elastic_Buffer::Get_Bytes_Per_Cell()
{
	if (_Is_Model)
		return 16;
	else
		return 24;
}

unsigned long Elastic_Buffer::Compute_Device_Memory_Block_Size()
{
	return (unsigned long)_prop->Get_Block_Size_X() * Get_Bytes_Per_Cell() * (unsigned long)(_y1-_y0+1) * (unsigned long)(_z1-_z0+1);
}

unsigned long Elastic_Buffer::Compute_Device_Memory_Requirement()
{
	return Compute_Device_Memory_Block_Size() * (unsigned long)_num_blocks;
}

void Elastic_Buffer::Add_To_YRange(int& min_y, int& max_y)
{
	if (_Is_Compute)
	{
		if (_cmp_y0 < min_y) min_y = _cmp_y0;
		if (_cmp_y1 > max_y) max_y = _cmp_y1;
	}
}

double Elastic_Buffer::Get_Relative_Cost()
{
	return _prop->Get_Relative_Cost(_Is_PV);
}

cudaStream_t Elastic_Buffer::Get_Compute_Stream()
{
	return _prop->Get_Compute_Stream(_device_id);
}

cudaStream_t Elastic_Buffer::Get_Input_Stream()
{
	return _prop->Get_Input_Stream(_device_id);
}

cudaStream_t Elastic_Buffer::Get_Output_Stream()
{
	return _prop->Get_Output_Stream(_device_id);
}

void Elastic_Buffer::_Find_Non_Zeros(char* dst, size_t len)
{
	size_t len_f = len / 4;
	for (size_t i = 0;  i < len_f;  ++i)
	{
		if (((float*)dst)[i] != 0.0f)
		{
			char name[256];
			printf("%s :: val at idx=%ld is %f\n",Get_Name_String(name),i,((float*)dst)[i]);
		}
	}
}

void Elastic_Buffer::Launch_Input_Transfers()
{
	if (_Is_Input)
	{
		size_t one_y_size = (size_t)((_z1 - _z0 + 1) * _prop->Get_Block_Size_X() * Get_Bytes_Per_Cell());
		if (_src == 0L)
		{
			// host-2-device
			size_t dst_off_y = (size_t)(_inp_y0 - _y0);
			size_t dst_len_y = (size_t)(_inp_y1 - _inp_y0 + 1);
			size_t src_off_y = (size_t)_inp_y0;
			_prop->Add_H2D(dst_len_y * one_y_size);
			int block_offset = Get_Block_Offset(_dst_block_id,0);
			//if (!_Is_Model)	_Find_Non_Zeros((char*)_prop->Get_Host_Block(block_offset, _Is_Model, _Is_PV, true) + src_off_y * one_y_size, dst_len_y * one_y_size);
			//printf("H2D H=%p D=%p\n",(char*)_prop->Get_Host_Block(block_offset, _Is_Model, _Is_PV, true) + src_off_y * one_y_size,(char*)Get_Block_By_Offset(block_offset,0) + dst_off_y * one_y_size);
			cudaMemcpyAsync(
				(char*)Get_Block_By_Offset(block_offset,0) + dst_off_y * one_y_size,
				(char*)_prop->Get_Host_Block(block_offset, _Is_Model, _Is_PV, true) + src_off_y * one_y_size,
				dst_len_y * one_y_size,
				cudaMemcpyHostToDevice,
				Get_Input_Stream()
				);
			//printf("Host to Device %2d :: %ld bytes\n",_device_id,dst_len_y * one_y_size);
		}
		else
		{
			// device-2-device
			size_t dst_off_y = (size_t)(_inp_y0 - _y0);
			size_t dst_len_y = (size_t)(_inp_y1 - _inp_y0 + 1);
			size_t src_off_y = (size_t)(_inp_y0 - _src->Get_Y0());
			int block_offset = Get_Block_Offset(_dst_block_id,0);
			//cudaSetDevice(_src->Get_Device_ID());
			cudaMemcpyPeerAsync(
				(char*)Get_Block_By_Offset(block_offset,0) + dst_off_y * one_y_size, _device_id, 
				(char*)_src->Get_Block_By_Offset(block_offset,0) + src_off_y * one_y_size, _src->Get_Device_ID(),
				dst_len_y * one_y_size,
				_src->Get_Output_Stream()
				);
			//printf("Device %2d to Device %2d :: %ld bytes\n",_src->Get_Device_ID(),_device_id,dst_len_y * one_y_size);
		}
	}
}

void Elastic_Buffer::Launch_Output_Transfers()
{
	if (_Is_Device2Host && (_Is_Model || Get_Block_Timestep(-1,0) > 0))
	{
		size_t one_y_size = (size_t)((_z1 - _z0 + 1) * _prop->Get_Block_Size_X() * Get_Bytes_Per_Cell());
		size_t dst_off_y = (size_t)_pipe->Get_Y0();
		size_t dst_len_y = (size_t)(_pipe->Get_Y1() - _pipe->Get_Y0() + 1);
		size_t src_off_y = (size_t)(_pipe->Get_Y0() - _y0);
		_prop->Add_D2H(dst_len_y * one_y_size);
		int block_offset = Get_Block_Offset(-1,0);
		cudaMemcpyAsync(
			(char*)_prop->Get_Host_Block(block_offset, _Is_Model, _Is_PV, false) + dst_off_y * one_y_size,
			(char*)Get_Block_By_Offset(block_offset,0) + src_off_y * one_y_size,
			dst_len_y * one_y_size,
			cudaMemcpyDeviceToHost,
			Get_Output_Stream()
			);
		//printf("Device %2d to Host :: %ld bytes\n",_device_id,dst_len_y * one_y_size);
	}
}

// launch a simple GPU kernel that copies data from src to dst buffer
void Elastic_Buffer::Launch_Compute_Kernel(bool Simple_Copy, float dti, Elastic_Shot* shot)
{
	if (_Is_Compute)
	{
		// return offset is in the range [0,_prop->Get_Number_Of_Blocks()-1]
		// compute block is always rightmost block
		int cmp_block_offset = Get_Block_Offset(0,0);
		int cmp_block_timestep = Get_Block_Timestep(0,0);
		void* cmp_block = Get_Block_By_Offset(cmp_block_offset,0);

		void* em_block = _em->Get_Block_By_Offset(cmp_block_offset,0);
		// if cmp_block_offset is 0, this call will return nil, which is intentional
		void* m1L_block = _inp_m1->Get_Block_By_Offset(cmp_block_offset-1,0);
		void* m1C_block = _inp_m1->Get_Block_By_Offset(cmp_block_offset,0);
		// if cmp block_offset is _prop->Get_Number_Of_Blocks()-1, this call will return nil, which is intentional
		void* m1R_block = _inp_m1->Get_Block_By_Offset(cmp_block_offset+1,0);
		
		void* m2C_block = _inp_m2->Get_Block_By_Offset(cmp_block_offset,0);

		int one_y_size = (int)((_z1 - _z0 + 1) * _prop->Get_Block_Size_X() * Get_Bytes_Per_Cell());
		int iem_offset = (int)(_cmp_y0 - _em->Get_Y0());
		int im2_offset = (int)(_cmp_y0 - _inp_m2->Get_Y0());
		int im1_offset = (int)(_cmp_y0 - _inp_m1->Get_Y0());
		int cmp_offset = (int)(_cmp_y0 - _y0);
		int ny = _cmp_y1 - _cmp_y0 + 1;
		int nz = _z1 - _z0 + 1;	

		int one_wf_size = one_y_size / 6;
		int em_one_word_size = one_wf_size;
		int em_one_y_size = em_one_word_size * 4;

		bool has_low_YHalo  = (_cmp_y0 - _inp_m1->Get_Y0()) >= 4 ? true : false;
		bool has_high_YHalo = (_inp_m1->Get_Y1() - _cmp_y1) >= 4 ? true : false;

		Elastic_Modeling_Job* job = _prop->Get_Job();
		bool isosphere = job->Use_Isotropic_Sphere_During_Source_Injection() && shot->Inject_Source(cmp_block_timestep);

#ifdef GPU_DEBUG
		char buf1[256];
		char buf2[256];
		char buf3[256];
		if (false)//if (strcmp(Get_Name_String(buf1), " 1-S") == 0)
		{
			printf("_em = %s, _inp_m1 = %s, _inp_m2 = %s\n",_em->Get_Name_String(buf1),_inp_m1->Get_Name_String(buf2),_inp_m2->Get_Name_String(buf3));
			printf("%4s Launch_Compute_Kernel (device=%d) :: cmp_block(offset=%d, timestep=%d), halos:lo:hi=%s:%s, em=%s, cmp=%s, m1=%s%s%s, m2=%s, iem_offset=%d, im2_offset=%d, im1_offset=%d, cmp_offset=%d\n",
					Get_Name_String(buf1), _device_id, cmp_block_offset, cmp_block_timestep,
					has_low_YHalo ? "Y" : "N", has_high_YHalo ? "Y" : "N",
					em_block == 0L ? "N" : "Y",
					cmp_block == 0L ? "N" : "Y",
					m1L_block == 0L ? "N" : "Y", m1C_block == 0L ? "N" : "Y", m1R_block == 0L ? "N" : "Y",
					m2C_block == 0L ? "N" : "Y",
					one_y_size,iem_offset,im2_offset,im1_offset,cmp_offset
			      );
			printf("ny=%d, nz=%d, em_one_y_size=%d, one_y_size=%d\n",ny,nz,em_one_y_size,one_y_size);
			printf("m1L_block = %p\n",m1L_block);
			printf("m1C_block = %p\n",m1C_block);
			printf("m1R_block = %p\n",m1R_block);
		}
#endif

		cudaSetDevice(_device_id);
		if (Simple_Copy)
		{
			Host_Simple_Copy_Kernel(
					Get_Compute_Stream(),
					(char*)cmp_block + cmp_offset * one_y_size,
					(char*)m2C_block + im2_offset * one_y_size,
					_prop->Get_Block_Size_X(),ny,nz,one_y_size
					);
		}
		else if (_Is_PV)
		{
			// launch P-V kernel
			Host_Propagate_Particle_Velocities_Kernel(
				cmp_block_timestep,
				Get_Compute_Stream(),
				4,
				cmp_block_offset * _prop->Get_Block_Size_X(),
				_cmp_y0,
				_cmp_y1,
				_inp_m1->Get_Y0(),
                                _inp_m1->Get_Y1(),
				job->Get_Propagation_NX(),
				job->Get_Propagation_NY(),
				job->Get_Propagation_NZ(),
				dti,
				(unsigned int*)((char*)em_block + iem_offset * em_one_y_size),
				(float*)((char*)cmp_block + cmp_offset * one_y_size),
				m1L_block != 0L ? (float*)((char*)m1L_block + im1_offset * one_y_size) : 0L,
				(float*)((char*)m1C_block + im1_offset * one_y_size),
				m1R_block != 0L ? (float*)((char*)m1R_block + im1_offset * one_y_size) : 0L,
				(float*)((char*)m2C_block + im2_offset * one_y_size),
				_C0,
				_C1,
				_C2,
				_C3,
				1.0 / job->Get_DX(),
				1.0 / job->Get_DY(),
				1.0 / job->Get_DZ(),
				job->Get_Vpvert_Avg_Top(),
				job->Get_Vpvert_Avg_Bot(),
				job->Get_NABC_SDX(),
				job->Get_NABC_SDY(),
				job->Get_NABC_TOP(),
				job->Get_NABC_BOT(),
				job->Get_IsoOrEarth_Model_Attribute_Min(job->Attr_Idx_Q,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Q,isosphere),
				job->Get_FQ(),
				job->Get_IsoOrEarth_Model_Attribute_Min(job->Attr_Idx_Density,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Density,isosphere),
				one_y_size,
				shot->Inject_Source(cmp_block_timestep),
				shot->Get_Source_Interpolation_Method(),
				shot->Get_Source_Is_Force(),
				shot->Get_Source_Is_Velocity(),
				shot->Get_Amplitude1(),
				shot->Get_Amplitude2(),
				shot->Get_Amplitude3(),
				shot->Get_Source_Wavelet_Sample(cmp_block_timestep),
				shot->Get_Propagation_Source_X(),
				shot->Get_Propagation_Source_Y(),
				shot->Get_Propagation_Source_Z()
				);
		}
		else
		{
			// launch S-T kernel
			float swav = 
				cmp_block_timestep > 1 ? 
				(shot->Get_Source_Wavelet_Sample(cmp_block_timestep-1) + shot->Get_Source_Wavelet_Sample(cmp_block_timestep))/2.0f : 
				shot->Get_Source_Wavelet_Sample(cmp_block_timestep)/2.0f;
			Host_Propagate_Stresses_Orthorhombic_Kernel(
				cmp_block_timestep,
				Get_Compute_Stream(),
				4,
				cmp_block_offset * _prop->Get_Block_Size_X(),
				_cmp_y0,
				_cmp_y1,
				_inp_m1->Get_Y0(),
                                _inp_m1->Get_Y1(),
				job->Get_Propagation_NX(),
				job->Get_Propagation_NY(),
				job->Get_Propagation_NZ(),
				dti,
				(unsigned int*)((char*)em_block + iem_offset * em_one_y_size),
				(float*)((char*)cmp_block + cmp_offset * one_y_size),
				m1L_block != 0L ? (float*)((char*)m1L_block + im1_offset * one_y_size) : 0L,
				(float*)((char*)m1C_block + im1_offset * one_y_size),
				m1R_block != 0L ? (float*)((char*)m1R_block + im1_offset * one_y_size) : 0L,
				(float*)((char*)m2C_block + im2_offset * one_y_size),
				_C0,
				_C1,
				_C2,
				_C3,
				1.0 / job->Get_DX(),
				1.0 / job->Get_DY(),
				1.0 / job->Get_DZ(),
				job->Get_Vpvert_Avg_Top(),
				job->Get_Vpvert_Avg_Bot(),
				job->Get_NABC_SDX(),
				job->Get_NABC_SDY(),
				job->Get_NABC_TOP(),
				job->Get_NABC_BOT(),
				job->Get_IsoOrEarth_Model_Attribute_Min  (job->Attr_Idx_Vp,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Vp,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Min  (job->Attr_Idx_Vs,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Vs,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Min  (job->Attr_Idx_Density,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Density,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Min  (job->Attr_Idx_Dip,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Dip,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Min  (job->Attr_Idx_Azimuth,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Azimuth,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Min  (job->Attr_Idx_Rake,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Rake,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Min  (job->Attr_Idx_Delta1,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Delta1,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Min  (job->Attr_Idx_Delta2,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Delta2,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Min  (job->Attr_Idx_Delta3,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Delta3,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Min  (job->Attr_Idx_Epsilon1,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Epsilon1,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Min  (job->Attr_Idx_Epsilon2,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Epsilon2,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Min  (job->Attr_Idx_Gamma1,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Gamma1,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Min  (job->Attr_Idx_Gamma2,isosphere),
				job->Get_IsoOrEarth_Model_Attribute_Range(job->Attr_Idx_Gamma2,isosphere),
				one_y_size,
				shot->Inject_Source(cmp_block_timestep),
				shot->Get_Source_Interpolation_Method(),
				shot->Get_Source_Is_Pressure(),
				shot->Get_Amplitude1(),
				swav,
				shot->Get_Propagation_Source_X(),
				shot->Get_Propagation_Source_Y(),
				shot->Get_Propagation_Source_Z()
				);
		}
	}
}

void Elastic_Buffer::Launch_Simple_Copy_Kernel()
{
	Launch_Compute_Kernel(true, 0.0f, 0L);
}

// compute block offset for future or past iteration
// iteration == 1 computes block offset for the next iteration,
// iteration == -1 computes block offset for the previous iteration etc.
void Elastic_Buffer::Advance_Block_Offset(int iteration, int& block_offset, int& timestep)
{
	block_offset += iteration;
	while (block_offset >= _prop->Get_Number_Of_Blocks())
	{
		block_offset -= _prop->Get_Number_Of_Blocks();
		timestep += _prop->Get_Total_Number_Of_Timesteps();
	}
	while (block_offset < 0)
	{
		block_offset += _prop->Get_Number_Of_Blocks();
		timestep -= _prop->Get_Total_Number_Of_Timesteps();
	}
}

void Elastic_Buffer::Shift_Buffer()
{
	// shift buffer one position to the left 
	void* tmp = _blocks[_num_blocks-1];
	for (int i = _num_blocks-1;  i > 0;  --i) _blocks[i] = _blocks[i-1];
	_blocks[0] = tmp;
	for (int i = 0;  i < _num_blocks;  ++i)
	{
		Advance_Block_Offset(1, _current_block_offset[i], _current_block_timestep[i]);
	}
}

void Elastic_Buffer::Launch_Data_Transfers()
{
	cudaSetDevice(_device_id);
	Launch_Input_Transfers();
	Launch_Output_Transfers();
}

void Elastic_Buffer::Free_Device_Blocks()
{
	if (_blocks != 0L)
	{
		for (int i = 0;  i < _num_blocks;  ++i)
		{
			_blocks[i] = 0L;
		}
		delete [] _blocks;
		_blocks = 0L;
	}
	if (_current_block_offset != 0L)
	{
		delete [] _current_block_offset;
		_current_block_offset = 0L;
	}
	if (_current_block_timestep != 0L)
	{
		delete [] _current_block_timestep;
		_current_block_timestep = 0L;
	}
}

// Reset block off offset and timestep counters, thus preparing for next shot.
void Elastic_Buffer::Reset()
{
	for (int i = 0;  i < _num_blocks;  ++i)
        {
		_current_block_offset[i] = Get_Relative_Block_Offset(i);
		_current_block_timestep[i] = _timestep;
	}
}

unsigned long Elastic_Buffer::Allocate_Device_Blocks(void* d_Mem, unsigned long offset)
{
	Free_Device_Blocks();
	cudaSetDevice(_device_id);
	size_t blkSize = (size_t)Compute_Device_Memory_Block_Size();
	_blocks = new void*[_num_blocks];
	_current_block_offset = new int[_num_blocks];
	_current_block_timestep = new int[_num_blocks];
	size_t my_offset = offset;
	for (int i = 0;  i < _num_blocks;  ++i)
	{
		_blocks[i] = (char*)d_Mem + my_offset;
		my_offset += blkSize;
		_current_block_offset[i] = Get_Relative_Block_Offset(i);
		_current_block_timestep[i] = _timestep;
		if (_current_block_offset[i] < 0)
		{
			_current_block_offset[i] += _prop->Get_Number_Of_Blocks();
			_current_block_timestep[i] -= _prop->Get_Total_Number_Of_Timesteps();
		}
	}
	return (unsigned long)my_offset;
}

void Elastic_Buffer::Enable_Peer_Access()
{
	if (_src != 0L) _prop->Enable_Peer_Access(_device_id, _src->Get_Device_ID());
}

