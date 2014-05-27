#ifndef CVX_SEISMOD_ELASTIC_BUFFER
#define CVX_SEISMOD_ELASTIC_BUFFER

#include <cuda_runtime_api.h>

class Elastic_Propagator;
class Elastic_Pipeline;
class Elastic_Shot;

class Elastic_Buffer
{
public:
	// Create earth model buffer
	Elastic_Buffer(
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
		);

	// Create input buffer
	Elastic_Buffer(
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
		);

	// Create compute buffer
	Elastic_Buffer(
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
		);
	
	~Elastic_Buffer();
	
	char* Get_Name_String(char* buf);

	int Get_Device_ID();

	void Reset();

	void Add_EM_Buffer(Elastic_Buffer* new_buffer);

	//
	// Buffer level methods
	//

	// set Is_Device2Host flag to supplied value
	bool Set_Is_Device2Host(int flag);
	bool Is_Device2Host();
	bool Is_Host2Device();

	bool Is_Input();
	bool Is_Compute();
	bool Is_Partial_Compute();

	bool Is_Particle_Velocity() {return _Is_PV;}

	//
	// Relative offset block level methods.
	//

	bool Block_Is_Device2Host_By_Relative_Offset(int relative_block_offset);
	bool Block_Is_Host2Device_By_Relative_Offset(int relative_block_offset);

	bool Block_Is_Input_By_Relative_Offset(int relative_block_offset);
	bool Block_Is_Compute_By_Relative_Offset(int relative_block_offset);
	bool Block_Is_Partial_Compute_By_Relative_Offset(int relative_block_offset);

	Elastic_Buffer* Get_Source_Buffer();
	int Get_Source_Block_Relative_Offset();

	Elastic_Buffer* Get_M1_Buffer();
	Elastic_Buffer* Get_M2_Buffer();

	int Get_Relative_Timestep();

	int Get_Min_Relative_Block_Offset();
	int Get_Max_Relative_Block_Offset();

	// returns true if a block with the relative block offset is included in this buffer.
	bool Block_Is_Included_By_Relative_Offset(int relative_block_offset);

	// get block offset relative to start of pipeline.
	// this number is always <= zero.
	int Get_Relative_Block_Offset(int block_id);

	//
	// Current block offset methods.
	// 

	// return true if block is valid.
	// blocks will remain invalid until the pipeline has filled up.
	// blocks that appear earlier in the pipeline will become valid first.
	bool Block_Is_Valid(int block_id, int iteration);

	// get block offset for current, past or future iteration.
	// this value will change with every call to Propagate_One_Block.
	// negative number until Block_Is_Valid(...) returns true.
	// after that, a value in the range [0,_prop->Get_Number_Of_Blocks-1]
	int Get_Block_Offset(int block_id, int iteration);

	// get timestep for current, past or future iteration.
	// this value will change with every call to Propagate_One_Block.
	// negative number until Block_Is_Valid(...) returns true.
	int Get_Block_Timestep(int block_id, int iteration);

	// compute block offset and timestep for past or future iteration.
	// calling this function with iteration == 0 does nothing.
	void Advance_Block_Offset(int iteration, int& block_offset, int& block_timestep);

	// get block by it's block offset. this is sometimes more convenient when you are looking for a match to a block
	// in another buffer.
	void* Get_Block_By_Offset(int block_offset, int iteration);

	int Get_Y0();
	int Get_Y1();

	int Get_Inp_Y0();
	int Get_Inp_Y1();

	int Get_Cmp_Y0();
	int Get_Cmp_Y1();

	int Get_Z0();
	int Get_Z1();

	double Get_Workload();

	unsigned long Get_Bytes_Per_Cell();

	// compute size of one block in bytes
	unsigned long Compute_Device_Memory_Block_Size();

	// compute aggregate size of all blocks in bytes
	unsigned long Compute_Device_Memory_Requirement();

	void Add_To_YRange(int& min_y, int& max_y);

	void Free_Device_Blocks();
	unsigned long Allocate_Device_Blocks(void* d_Mem, unsigned long offset);
	void Enable_Peer_Access();

	cudaStream_t Get_Compute_Stream();
	cudaStream_t Get_Input_Stream();
	cudaStream_t Get_Output_Stream();

	void Launch_Data_Transfers();
	void Launch_Simple_Copy_Kernel();

	/* Original coefficients 
	static const float _C0 = 1225.0f/1024.0f;
	static const float _C1 = -245.0f/3072.0f;
	static const float _C2 = 49.0f/5120.0f;
	static const float _C3 = -5.0f/7168.0f; */
	
	// Optimized (for lower dispersion) coefficients 
	static const float _C0 = 1.1850912100109303f;
	static const float _C1 = -0.0835270299926924f;
	static const float _C2 = 0.016837760894350576f;
	static const float _C3 = -0.0027386181103430177f;

private:
	friend class Elastic_Pipeline;

	Elastic_Propagator* _prop;	// Root object
	Elastic_Pipeline* _pipe;	// Parent object
	int _device_id;			// GPU device id
	int _timestep;
	bool _Is_PV;
	bool _Is_Model;			// true if this is an earth model buffer
	bool _Is_Input;			// true if this is an input buffer
	bool _Is_Device2Host;		// true if content of this buffer gets sent back to host
	bool _Is_Compute;		// true if this is a compute buffer
	bool _Is_Partial;		// true if this is a partial compute buffer
	int _y0;			// y coordinates of buffer
	int _y1;
	int _z0;			// z coordinates of buffer
	int _z1;
	int _inp_y0;			// input y coordinates, only valid if _Is_Input is true
	int _inp_y1;
	int _cmp_y0;			// compute y coordinate, only valid if _Is_Compute is true
	int _cmp_y1;
	int _num_blocks;		// number of blocks in this buffer
	int _block_offset;		// block offset of the rightmost block (block_idx == 0)
	int _cmp_block_id;		// computed results go to this block
	Elastic_Buffer* _em;		// Earth model buffer
	Elastic_Buffer* _inp_m2;	// input buffer, two steps removed - must be same device id as this
	Elastic_Buffer* _inp_m1;	// input buffer, one step removed - must be same device id as this
	Elastic_Buffer* _src;		// source buffer for copy operations - must be different device id
	int _dst_block_id;		// copy to this buffer into the block with id equal to _dst_block_id

	unsigned long _allocated_bytes;	// total number of allocated bytes on device
	unsigned long _failed_to_allocate_bytes;	// total number of bytes that failed to allocate on device

	int* _current_block_offset;	// current block offset for each block
	int* _current_block_timestep;	// current timestep for each block
	void** _blocks;			// pointer to the blocks in device memory

	void _Find_Non_Zeros(char* dst, size_t len);
	double Get_Relative_Cost();

	void Shift_Buffer();

	void Launch_Input_Transfers();
	void Launch_Output_Transfers();
	void Launch_Compute_Kernel(bool Simple_Copy, float dti, Elastic_Shot* shot);
};

#endif

