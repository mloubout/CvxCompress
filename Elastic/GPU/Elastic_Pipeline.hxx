#ifndef CVX_SEISMOD_ELASTIC_PIPELINE
#define CVX_SEISMOD_ELASTIC_PIPELINE

class Elastic_Buffer;
class Elastic_Propagator;
class Elastic_Shot;

class Elastic_Pipeline
{
public:
	Elastic_Pipeline(
		int log_level,
		Elastic_Propagator* prop,
		int pipe_id,
		int pipe_y0,
		int pipe_y1,
		int pipe_z0,
		int pipe_z1
		);

	~Elastic_Pipeline();

	int Get_ID() {return _pipe_id;}

	void Reset();

	void Append_Buffer(Elastic_Buffer* new_buffer);
	void Add_EM_Buffer(Elastic_Buffer* new_buffer);
	
	int Get_Y0() {return _pipe_y0;}
	int Get_Y1() {return _pipe_y1;}
	int Get_Width() {return _pipe_y1 - _pipe_y0 + 1;}

	int Get_Number_Of_Buffers();	
	Elastic_Buffer* Get_Buffer(int index);

	// get the number of devices contributing to this pipeline
	int Get_Device_Count();

	int* Get_All_Device_IDs();

	// get the total number of timesteps a block is propagated by this pipeline
	int Get_Total_Number_Of_Timesteps();

	// get the total spatial shift for this pipeline.
	// note that this is a negative number.
	int Get_Total_Block_Offset();

	// Compute total memory requirement for one device
	unsigned long Compute_Device_Memory_Requirement(int device_id);
	void Compute_RX_Device_Memory_Requirement(int device_id, int& rxloc_size, int& rxres_size);

	bool Verify_All_Devices_Have_Enough_Memory();

	// get the input and output block offsets for the pinned buffers for current, past or future iteration.
	// returns -1 if block offsets are not valid yet.
	int Get_Input_Block_Offset(int iteration);
	int Get_Output_Block_Offset(int iteration);

	int Get_Input_Block_Timestep(int iteration);
	int Get_Output_Block_Timestep(int iteration);

	// Get total workload of this pipeline
	double Get_Workload();
	
	// Get total workload of this pipeline without the halo overhead
	double Get_Minimum_Workload();

	// Get workload for one device in this pipeline
	double Get_Workload(int device_id);
	
	// Get computational overhead for one device in this pipeline.
	// A value of 1.0 means no overhead, 1.5 means 50% overhead.
	double Get_Computational_Overhead(int device_id);

	void Print_Graphical(int device_id);
	void Print_Graphical();

	bool Allocate_Device_Memory();
	void Free_Device_Memory();

	void Allocate_RxLoc_Buffer(Elastic_Shot* shot);
	void Free_RxLoc_Buffer(Elastic_Shot* shot);

	void Launch_Receiver_Data_Transfers(Elastic_Shot* shot);
	void Launch_Receiver_Extraction_Kernels(Elastic_Shot* shot);
	void Launch_Data_Transfers();	
	void Launch_Simple_Copy_Kernel();
	void Launch_Compute_Kernel(float dti, Elastic_Shot* shot, int num_z);

	void DEMUX_Receiver_Values_For_One_Device(Elastic_Shot* shot, int device_index);
	void DEMUX_Receiver_Values(Elastic_Shot* shot);

private:
	friend class Elastic_Propagator;
	Elastic_Propagator* _prop;

	int _log_level;

	int _pipe_id;
	int _pipe_y0;
	int _pipe_y1;
	int _pipe_z0;
	int _pipe_z1;

	int _num_buffers;
	Elastic_Buffer** _buffers;
	void Shift_Buffers();

	int* _optseq_launch;
	int* _optseq_data_in;
	int* _optseq_data_out;
	int _num_optseq_launch;
	int _num_optseq_data_in;
	int _num_optseq_data_out;

	int _num_devices;
	int* _device_IDs;

	void** _d_Mem;				// pointer to ALL the device memory allocated

	void** _d_RxLoc;			// points to device memory reserved for receiver locations
	void*** _d_RxLoc_block;			// receiver locations are organized into blocks
	int** _h_RxLoc_block_Offset;		// block_offset for each receiver location block
	int* _h_RxLoc_num_blocks;

	void** _d_RxRes;			// points to device memory reserved for receiver results

	void** _h_RxRes_curr;			// points to host memory reserved for receiver results
	int** _h_RxRes_curr_block_offset;	// block offsets for receiver results
	int** _h_RxRes_curr_timestep;		// timesteps for receiver results
	int** _h_RxRes_curr_num_rx;		// number of receivers for receiver results
	int** _h_RxRes_curr_flags;		// flag for receiver results
	int* _h_RxRes_curr_num_blocks;		// number of blocks
	
	void** _h_RxRes_prev;			// points to host memory reserved for receiver results
	int** _h_RxRes_prev_block_offset;	// block offsets for receiver results
	int** _h_RxRes_prev_timestep;		// timesteps for receiver results
	int** _h_RxRes_prev_num_rx;		// number of receivers for receiver results
	int** _h_RxRes_prev_flags;		// flag for receiver results
	int* _h_RxRes_prev_num_blocks;		// number of blocks

	void _Shift_Rx_Blocks();

	// get a list containing the device IDs of every buffer.
	// each device ID appears only once in this list.
	void _Compile_Device_IDs();

	bool Block_Is_Output_By_Relative_Offset(Elastic_Buffer* buffer, int relative_block_offset);
};

#endif

