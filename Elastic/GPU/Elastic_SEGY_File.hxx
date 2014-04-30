#ifndef CVX_ESDRD_MI_TMJ_ELASTIC_SEGY_FILE_HXX
#define CVX_ESDRD_MI_TMJ_ELASTIC_SEGY_FILE_HXX

class Elastic_SEGY_File_Receiver_Range;
class Elastic_Propagator;

class Elastic_SEGY_File
{
public:
	Elastic_SEGY_File(
		int fileidx,
		const char* base_filename,
		double sample_rate,
		double tshift,
		double reclen,
		bool do_P,
		bool do_Vx,
		bool do_Vy,
		bool do_Vz
		);
	~Elastic_SEGY_File();

	bool Is_Valid() {return _Is_Valid;}

	int Get_File_Index() {return _fileidx;}

	void Add_Receiver_Range_X(
		int range_idx,
		double start,
		double end,
		double interval
		);
	void Add_Receiver_Range_Y(
		int range_idx,
		double start,
		double end,
		double interval
		);
	void Add_Receiver_Range_Z(
		int range_idx,
		double start,
		double end,
		double interval
		);

	int Compute_Receiver_Locations(
		float*& rcv_x,
		float*& rcv_y,
		float*& rcv_z
		);	

	void Shift_Receiver_Transfer_Buffers();

	void Create_Receiver_Transfer_Buffers(
			Elastic_Propagator* prop
			);

	bool Create_New_Device_To_Host_Transfer(
			int device_id,
			int block_number,
			int timestep,
			float*& dst_buf,
			int& dst_size
			);

private:
	bool _Is_Valid;
	int _fileidx;
	char* _base_filename;
	double _sample_rate;
	double _tshift;
	double _reclen;
	bool _do_P;
	bool _do_Vx;
	bool _do_Vy;
	bool _do_Vz;

	Elastic_SEGY_File_Receiver_Range** _rcv_ranges;
	int _num_rcv_ranges;

	Elastic_SEGY_File_Receiver_Range* _Get_Receiver_Range(int range_idx);

	void _Destroy_Buffers();

	int _totSteps;			// total number of timesteps performed by entire pipeline
	int _nWF;			// number of output wavefields
	int _nBlks;			// number of blocks, i.e. Elastic_Propagator::Get_Number_Of_Blocks()
	int _num_pipes;			// number of pipelines, i.e. Elastic_Propagator::Get_Number_Of_Pipelines()
	float*** _rcv_x;		// receiver locations _rcv_x[Block][Pipeline][idx]
	float*** _rcv_y;		// receiver locations
	float*** _rcv_z;		// receiver locations
	int** _rcv_n;			// number of receivers _rcv_n[Block][Pipeline]
	int** _rcv_i;			// receiver location transfer buffer index _rcv_i[Block][Pipeline]. used to calculate destination ptr for device-to-host transfer.
	int* _rcv_nn;			// number of receivers per block. _rcv_nn[Block]
	int*** _trc_idx;		// trace index (order traces should appear in SEGY file). _trc_idx[Block][Pipeline][idx]
	int* _max_rx;			// maximum receiver count for all blocks in a pipeline. _max_rx[Pipeline]
	int* _tot_rx;			// total receiver count for each pipeline. _tot_rx[Pipeline]

	float* _h_transfer;		// pinned host transfer buffer
	float*** _h_rx;			// ptrs to _h_transfer. _h_rx[Block][Timestep] 
	int* _h_prev_ts;		// timestep for previous buffer transfer.
	int* _h_curr_ts;		// timestep for current (still ongoing) buffer transfer.
	int* _device2pipe;		// mapping from device id to pipe id.

	int _max_device_id;		// highest device id among the utilized GPUs
	float** _d_rcv_loc;		// device memory buffer for receiver locations and device-to-host transfer. _d_rcv_loc[device_id]
	float*** _d_rcv_x;		// device memory ptr for rcv_x. _d_rcv_x[device_id][block]
	float*** _d_rcv_y;		// device memory ptr for rcv_y. _d_rcv_x[device_id][block]
	float*** _d_rcv_z;		// device memory ptr for rcv_z. _d_rcv_x[device_id][block]
	float***** _d_transfer;		// device memory ptr for device-to-host transfer. _d_transfer[device_id][timestep][wf(0->P,1->Vx,2->Vy,3->Vz)][usage(0->comp,1->transfer_out)]
};

#endif

