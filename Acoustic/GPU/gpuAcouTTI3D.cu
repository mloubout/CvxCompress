#include <stdio.h>
#include <exception>
#include <omp.h>
#include <unistd.h>
#include <cuda_runtime.h>

#include "Parmfile_Reader.hxx"
#include "Acoustic_Earth_Model.hxx"

#define MB 1048576

//
// Trial implementation of "scanner" GPU kernel implementation.
// Scanner direction is Y.
// Axis ordering is X-Z-Y.
//

//
// simple kernel that just copies from src to dst.
//
__global__ void SimpleCopyKernel(int device, int layer_id, float* d_dst, float* d_src, int halo_ny, int y_offset, int dimz)
{
	//const unsigned long stride_x = 1;
	const unsigned long stride_y = (unsigned long)dimz * 16;
	const unsigned long stride_z = 32;

	int thr_X = threadIdx.x;
	int thr_Y = threadIdx.y + blockIdx.y * blockDim.y;
	int thr_Z = threadIdx.z + blockIdx.z * blockDim.z;

	float2* dst = (float2*)d_dst + (thr_X + (thr_Y + y_offset) * stride_y + thr_Z * stride_z);
	float2* src = (float2*)d_src + (thr_X + (thr_Y + halo_ny) * stride_y + thr_Z * stride_z);

	for (int iZ = 0;  iZ < dimz/2;  iZ+=gridDim.z*blockDim.z)
	{
		float2 val = src[iZ*stride_z];
		dst[iZ*stride_z] = val;
	}
}

#include "Acoustic_Kernels_Sponges.cpp"

//#define NAN_DESU_KA
#include "TTIDenQ_T2_GPU.cu"

class TimeStepLayer
{
public:

	TimeStepLayer(int arg_device, char* arg_name, int arg_is_compute_layer, int arg_is_p_layer, int arg_timestep, int arg_bytes_per_cell, int arg_total_num_timesteps, int arg_total_num_stripes, int arg_y0)
	{
		device = arg_device;
		d_buf = 0L;

		name = 0L;
		setName(arg_name);

		num_layers_sharing_device_resources = 0;
		Layers_Sharing_Device_Resources = 0L;
	
		is_compute_layer = arg_is_compute_layer;
		is_p_layer = arg_is_p_layer;
		is_output_layer = 0;
		timestep = arg_timestep;
		bytes_per_cell = arg_bytes_per_cell;
		total_num_timesteps = arg_total_num_timesteps;
		total_num_stripes = arg_total_num_stripes;
		y0 = arg_y0;
		stripe_offset = 0;
		input_off_x = 0;
	}

	~TimeStepLayer()
	{
		//printf("~TimeSteplayer\n");
		fflush(stdout);
		if (name != 0L) delete [] name;
		if (d_buf != 0L)
		{
			cudaFree(d_buf);
			d_buf = 0L;
		}
		if (Layers_Sharing_Device_Resources != 0L && this == Layers_Sharing_Device_Resources[0])
		{
			delete [] Layers_Sharing_Device_Resources;
		}
		Layers_Sharing_Device_Resources = 0L;
	}

	void setName(char* arg_name)
	{
		if (name != 0L) delete [] name;
		name = 0L;

		if (arg_name != 0L)
		{
			int len = strlen(arg_name);
			if (len > 0)
			{
				name = new char[len+1];
				memcpy((void*)name, (void*)arg_name, len);
				name[len] = 0;
			}
		}
	}

	unsigned long Compute_Storage_Requirement_Single_Stripe_Num_Cells()
	{
		if (num_layers_sharing_device_resources > 0)
		{
			if (this == Layers_Sharing_Device_Resources[0])
			{
				unsigned long maxReqMem = 0;
				for (int i = 0;  i < num_layers_sharing_device_resources;  ++i)
				{
					unsigned long ReqMem = Layers_Sharing_Device_Resources[i]->Private_Compute_Storage_Requirement_Single_Stripe_Num_Cells();
					if (ReqMem > maxReqMem) maxReqMem = ReqMem;
				}
				return maxReqMem;
			}
			else
			{
				// this layer shares device resources with one or more other layers, so it doesn't have any requirements of it's own.
				return 0L;
			}
		}
		else
		{
			return Private_Compute_Storage_Requirement_Single_Stripe_Num_Cells();
		}
	}

	unsigned long Compute_Storage_Requirement_Single_Stripe_Bytes()
	{
		return Compute_Storage_Requirement_Single_Stripe_Num_Cells() * (unsigned long)bytes_per_cell;
	}

	unsigned long Compute_Storage_Requirement()
	{
		return Compute_Storage_Requirement_Single_Stripe_Bytes() * (unsigned long)num_stripes;
	}

	void Print()
	{
		printf("%s\t%2d stripes, storage %d by %d by %d",name,num_stripes,16,output_dim_y,output_dim_z);
		if (is_compute_layer)
		{
			printf(", timestep %d, compute offset %d,%d,%d, compute dimension %d by %d by %d, compute out %d,%d,%d\n",timestep,compute_off_x,compute_off_y,compute_off_z,16,compute_dim_y,compute_dim_z,0,compute_out_y,compute_out_z);
		}
		else
		{
			printf(", timestep %d, input offset %d,%d,%d\n",timestep,input_off_x,input_off_y,input_off_z);
		}
	}

	//
	// Shift the stripes in the layer. Each layer is a sliding window, implemented as a ring buffer.
	// Shifting the ring buffer has the same effect as moving the sliding window one position.
	//
	void Move_Sliding_Window(int logLevel)
	{
		if (num_layers_sharing_device_resources > 0)
		{
			Layers_Sharing_Device_Resources[0]->Private_Move_Sliding_Window(logLevel);
		}
		else
		{
			Private_Move_Sliding_Window(logLevel);
		}
	}

	//
	// Get an indexed stripe from the layer. Each layer is a sliding window with num_stripes stripes in it.
	// Each stripe holds 16 consecutive X locations.
	// stripe index can be in range [-num_stripes,num_stripes>.
	// Negative index count from end of sliding window, so -1 is the last stripe in the window.
	// Last stripe is always populated most recent.
	//
	void* Get_Stripe(int stripe_number)
	{
		if (num_layers_sharing_device_resources > 0)
                {
			return Layers_Sharing_Device_Resources[0]->Private_Get_Stripe(stripe_number);
		}
		else
		{
			return Private_Get_Stripe(stripe_number);
		}
	}

	// Returns non-zero (true) if allocation succeeds.
	int Allocate_Device_Resources()
	{
		if (num_layers_sharing_device_resources > 0)
                {
			if (this == Layers_Sharing_Device_Resources[0])
			{
                        	return Layers_Sharing_Device_Resources[0]->Private_Allocate_Device_Resources();
			}
			else
			{
				// no device resources to allocate, so will never fail.
				return 1;
			}
                }
                else
                {
                        return Private_Allocate_Device_Resources();
                }
	}

	//
	// Compute Y extent this layer.
	// Y coordinates are not clipped.
	//
	void Compute_Y_Extent(
		int pipeline_net_y0,
		int pipeline_net_y1,
		int& layer_y0,
		int& layer_y1
		)
	{
		if (is_compute_layer)
		{
			layer_y0 = pipeline_net_y0 + compute_off_y - compute_out_y;
			layer_y1 = layer_y0 + output_dim_y - 1;
		}
		else
		{
			layer_y0 = pipeline_net_y0 + input_off_y;
			layer_y1 = layer_y0 + output_dim_y - 1;
		}
	}

	//
	// Compute parameters when this layer is used as src input in a stripe copy.
	// src_y0, src_y1 - y range of stripe to be copied.
	// src_off_y - copy should start at this offset in stripe.
	// d_src - address of stripe to be copied in device memory.
	//
	void Compute_Source_Parameters(
		int pipeline_net_y0,
                int pipeline_net_y1,
		int dimy,
		int stripe_id,
		int& src_y0,
		int& src_y1,
		int& src_off_y,
		void*& d_src
		)
	{
		if (is_compute_layer)
		{
			src_off_y = compute_out_y;
			src_y0 = pipeline_net_y0 + compute_off_y - compute_out_y;
			src_y1 = src_y0 + compute_dim_y - 1;
			d_src = Get_Stripe(stripe_id);
		}
		else
		{
			src_off_y = 0;
			src_y0 = pipeline_net_y0 + input_off_y;
			src_y1 = src_y0 + output_dim_y - 1;
			if (src_y0 < 0)
			{
				src_off_y = -src_y0;	
				src_y0 = 0;
			}
			else
			{
				src_off_y = 0;
			}
			if (src_y1 >= dimy) src_y1 = dimy - 1;
			d_src = Get_Stripe(stripe_id);
		}
	}

	char* name;

	int num_stripes;	// number of stripes held by this layer
	int total_num_stripes;	// total number of stripes in propagation volume.
	int bytes_per_cell;	// bytes of storage required per cell.

	int is_compute_layer;	// 0->input layer, non-zero means compute layer.
	int is_p_layer;		// non-zero if this is a p wavefield layer.
	int is_output_layer;	// last layer, so need to generate some outputs.

	int timestep;		// timestep of this layer. first timestep is 1, second is 2 and so forth.
	int total_num_timesteps;	// total number of timesteps for entire pipeline

	int input_off_x;	// x offset of first cell in input buffer.
	int input_off_y;	// y offset of first cell in input buffer.
	int input_off_z;	// z offset of first cell in input buffer.

	int compute_dim_y;	// number of y cells to compute by kernel
	int compute_dim_z;	// number of z cells to compute by kernel

	int compute_out_y;	// compute kernel starts writing results at this y offset.
	int compute_out_z;	// compute kernel starts writing results at this z offset.

	int y0;			// y0 of pipeline. y0 + compute_off_y = y coordinate of first y in input block.

	int compute_off_x;	// x offset of this layer relative to input from host block.
	int compute_off_y;	// y offset of this layer relative to input from host block.
	int compute_off_z;	// z offset of this layer relative to input from host block.

	int output_dim_y;	// size of output buffer in y.
	int output_dim_z;	// size of output buffer in z.

	int curr_input_stripe;	// current absolute input stripe. 0 -> 0th stripe of 0th timestep. total_num_stripes -> 0th stripe of 1st timestep, 2*total_num_stripes -> 0th stripe of 2nd timestep etc.

	int curr_input_stripe_number_1;	// current stripe number in range [0,total_num_stripes> [first stripe]
	int curr_input_timestep_1;	// current input timestep. After total_num_stripes have been processed, this value is incremented by one and curr_input_stripe_number is reset to 0.
	int curr_input_stripe_number_2;	// current stripe number in range [0,total_num_stripes> [second stripe]
	int curr_input_timestep_2;	// current input timestep. After total_num_stripes have been processed, this value is incremented by one and curr_input_stripe_number is reset to 0.

	int curr_output_stripe_number_1;
	int curr_output_timestep_1;
	int curr_output_stripe_number_2;
	int curr_output_timestep_2;

private:
	int device;

	// device memory

	float* d_buf;
	int stripe_offset;

	unsigned long Private_Compute_Storage_Requirement_Single_Stripe_Num_Cells()
	{
		return (unsigned long)16 * (unsigned long)output_dim_y * (unsigned long)output_dim_z;
	}

	void Private_Move_Sliding_Window(int logLevel)
	{
		stripe_offset = (stripe_offset + 2) % num_stripes;

		// update slice indexes so we know which part of the model we are working on.

		++curr_input_stripe;

		curr_input_stripe_number_1 = (curr_input_stripe + total_num_stripes) % total_num_stripes;
		curr_input_timestep_1 = ((curr_input_stripe + total_num_stripes) / total_num_stripes) - 1;
		curr_input_timestep_1 = curr_input_timestep_1 * total_num_timesteps + timestep;

		curr_output_stripe_number_1 = (curr_input_stripe - 2 + total_num_stripes) % total_num_stripes;
		curr_output_timestep_1 = ((curr_input_stripe - 2 + total_num_stripes) / total_num_stripes) - 1;
		curr_output_timestep_1 = curr_output_timestep_1 * total_num_timesteps + timestep;

		++curr_input_stripe;

		curr_input_stripe_number_2 = (curr_input_stripe + total_num_stripes) % total_num_stripes;
		curr_input_timestep_2 = ((curr_input_stripe + total_num_stripes) / total_num_stripes) - 1;
		curr_input_timestep_2 = curr_input_timestep_2 * total_num_timesteps + timestep;

		curr_output_stripe_number_2 = (curr_input_stripe - 2 + total_num_stripes) % total_num_stripes;
		curr_output_timestep_2 = ((curr_input_stripe - 2 + total_num_stripes) / total_num_stripes) - 1;
		curr_output_timestep_2 = curr_output_timestep_2 * total_num_timesteps + timestep;

		if (logLevel >= 5)
		{
			printf("device=%d %s :: inp1 is %d,%d - out1 is %d,%d :: inp2 is %d,%d - out2 is %d,%d\n",device,name,curr_input_timestep_1,curr_input_stripe_number_1,curr_output_timestep_1,curr_output_stripe_number_1,curr_input_timestep_2,curr_input_stripe_number_2,curr_output_timestep_2,curr_output_stripe_number_2);
		}
	}

	void* Private_Get_Stripe(int stripe_number)
	{
		void* retval = 0L;
		if (d_buf != 0L)
		{
			if (stripe_number >= -num_stripes && stripe_number < num_stripes)
			{
				int idx = (stripe_number + stripe_offset + num_stripes) % num_stripes;
				//printf("Get_Stripe :: dev=%d:%s, stripe_number=%d, stripe_offset=%d, num_stripes=%d, idx=%d\n",device,name,stripe_number,stripe_offset,num_stripes,idx);
				retval = d_buf + Compute_Storage_Requirement_Single_Stripe_Num_Cells() * (unsigned long)2 * (unsigned long)idx;
				//unsigned long off_buf = (unsigned long)((char*)retval - (char*)d_buf);
				//printf("off_buf = %ld\n",off_buf);
			}
		}
		if (retval == 0L)
		{
			printf("Warning! Get_Stripe return value is NIL\n");
		}
		return retval;
	}

	int Private_Allocate_Device_Resources()
	{
		unsigned long ReqMem = Compute_Storage_Requirement();
		if (ReqMem > 0)
		{
			// NB! Assumes cudaSetDevice is called by (parent) SingleDevicePropagator
			cudaError_t err = cudaMalloc((void**)&d_buf, ReqMem);
			if (err != cudaSuccess)
			{
				printf("Error! Failed to allocate device memory for layer %s on device %d\n",name,device);
				return 0;
			}
			else
			{
				cudaMemset(d_buf, 0, ReqMem);
				printf("Allocated %ld bytes of device memory for layer %s on device %d.\n",ReqMem,name,device);
				/*
				printf("Press ENTER to continue...");
				fflush(stdout);
				char buf[1024];
				gets(buf);
				*/
				return 1;
			}
		}
		else
		{
			return 0;
		}
	}

	friend class DevicePropagator;
	int num_layers_sharing_device_resources;
	TimeStepLayer** Layers_Sharing_Device_Resources;
};

class SingleDevicePropagator
{
public:
	//
	// Create device propagator to be placed in front of pipeline element
	//
	SingleDevicePropagator(
			int arg_Kernel_Type,
			SingleDevicePropagator* arg_Succeeding_Propagator,
			char* arg_device_name,
			int arg_device,
			int arg_num_timesteps,
			int arg_total_num_timesteps,
			int arg_total_num_stripes,
			int arg_y0				// y0 of pipeline
			)
	{
		Succeeding_Propagator = arg_Succeeding_Propagator;
		halo_ny = Succeeding_Propagator->halo_ny;
		exclude_halo_ylo = Succeeding_Propagator->exclude_halo_ylo;
		exclude_halo_yhi = Succeeding_Propagator->exclude_halo_yhi;
		bY = Succeeding_Propagator->layers[3]->compute_dim_y + (!exclude_halo_ylo ? halo_ny : 0) + (!exclude_halo_yhi ? halo_ny : 0);
		Create_Configuration(
			arg_Kernel_Type,
			arg_device_name,
			arg_device,
			arg_num_timesteps,
			arg_total_num_timesteps,
			arg_total_num_stripes,
			Succeeding_Propagator->OTflag,
			bY,
			Succeeding_Propagator->layers[3]->compute_dim_z,
			halo_ny,
			exclude_halo_ylo,
			exclude_halo_yhi,
			arg_y0
			);
	}

	SingleDevicePropagator(
			int arg_Kernel_Type,
			char* arg_device_name,
			int arg_device, 		// not needed?
			int arg_num_timesteps,		// number of timesteps to be performed on this GPU
			int arg_total_num_timesteps,	// total number of timesteps for entire pipeline
			int arg_total_num_stripes,	// total number of stripes in propagation volume
			int arg_OTflag,			// temporal order, either 2 or 4
			int arg_bY,			// output block size Y
			int arg_bZ,			// output block size Z
			int arg_halo_ny,		// number of halo cells Y, one side
			int arg_exclude_halo_ylo,	// flag: true means halo should be excluded on low Y side for each output
			int arg_exclude_halo_yhi,
			int arg_y0			// y0 of pipeline
			)
	{
		Succeeding_Propagator = 0L;
		Create_Configuration(arg_Kernel_Type,arg_device_name,arg_device,arg_num_timesteps,arg_total_num_timesteps,arg_total_num_stripes,arg_OTflag,arg_bY,arg_bZ,arg_halo_ny,arg_exclude_halo_ylo,arg_exclude_halo_yhi,arg_y0);	
	}

	~SingleDevicePropagator()
	{
		//printf("~SingleDevicePropagator\n");
		fflush(stdout);
		cudaSetDevice(device);
		if (layers != 0L)
		{
			for (int i = 0;  i < num_layers;  ++i) delete layers[i];
			delete [] layers;
			layers = 0L;
		}
		num_layers = 0;
		if (device_name != 0L) delete [] device_name;
		if (input_stream != 0)
		{
			cudaStreamDestroy(input_stream);
			input_stream = 0;
		}
		if (compute_stream != 0)
		{
			cudaStreamDestroy(compute_stream);
			compute_stream = 0;
		}
		if (output_stream != 0)
		{
			cudaStreamDestroy(output_stream);
			output_stream = 0;
		}
		if (d_temp1 != 0L) cudaFree(d_temp1);
		if (d_temp2 != 0L) cudaFree(d_temp2);
		if (d_temp3 != 0L) cudaFree(d_temp3);
		if (d_spg_x != 0L) cudaFree(d_spg_x);
		if (d_spg_y != 0L) cudaFree(d_spg_y);
		if (d_spg_z != 0L) cudaFree(d_spg_z);
		Free_Receiver_Locations();
	}

	void Put_Source_Wavelet(
		float* src,
		int len,
		int arg_xsrc,
		int arg_ysrc,
		int arg_zsrc,
		int arg_zsrcghost
		)
	{
		h_src_wavelet = src;
		src_len = len;
		xsrc = arg_xsrc;
		ysrc = arg_ysrc;
		zsrc = arg_zsrc;
		zsrcghost = arg_zsrcghost;
	}

	int Put_Receiver_Locations(
		int* arg_rcx,
		int* arg_rcy,
		int* arg_rcz,
		int arg_recnum,
		int dimx,
		int y0,
		int y1,
		int recghost_Flag,
		int z_freesurface
		)
	{
		int retval = -1;

		Free_Receiver_Locations();

		int nx = dimx / 16;
		h_totrecnx = arg_recnum;

		_recghost_Flag = recghost_Flag;
		_z_freesurface = z_freesurface;

		// determine how many receivers are covered by this pipeline
		int totcnt = 0;
		int* bincnt = new int[nx];
		for (int i = 0;  i < nx;  ++i) bincnt[i] = 0;
		for (int i = 0;  i < arg_recnum;  ++i)
		{
			int x = arg_rcx[i];
			int y = arg_rcy[i];
			int bin = x / 16;	
			if (y >= y0 && y <= y1)
			{
				++totcnt;
				++bincnt[bin];
			}
		}
		int maxbincnt = bincnt[0];
		for (int i = 1;  i < nx;  ++i)
		{
			if (bincnt[i] > maxbincnt) maxbincnt = bincnt[i];
		}
		h_rcoutmaxbinsize = (unsigned long)maxbincnt;
		printf("Number of bins = %d, totcnt = %d, maxbincnt = %d\n",nx,totcnt,maxbincnt);

		if (totcnt > 0)
		{
			unsigned long rcbufSize = (unsigned long)12 * (unsigned long)totcnt;
			cudaError_t err1 = cudaMalloc((void**)&d_rcbuf, rcbufSize);
			if (err1 != cudaSuccess)
			{
				printf("Error! Failed to allocate device memory for receiver buffer on device %d\n",device);
				d_rcbuf = 0L;
				retval = 0;
			}
			else
			{
				printf("Allocated %ld bytes of device memory for receiver buffer on device %d.\n",rcbufSize,device);
			}

			if (retval != 0)
			{
				unsigned long rcoutSize = (unsigned long)(16 * num_timesteps) * (unsigned long)maxbincnt;  // four buffers per timestep, elements are floats
				cudaError_t err2 = cudaMalloc((void**)&d_rcout, rcoutSize);
				if (err2 != cudaSuccess)
				{
					printf("Error! Failed to allocate device memory for receiver output buffer on device %d\n",device);
					d_rcout = 0L;
					retval = 0;
				}
				else
				{
					printf("Allocated %ld bytes of device memory for receiver buffer on device %d.\n",rcoutSize,device);
				}
				d_rcoutflip = d_rcout;
				d_rcoutflop = d_rcout + (2 * num_timesteps * maxbincnt);

				if (retval != 0)
				{
					h_rcbuf = new int[3*totcnt];
					unsigned long rcoutHostSize = rcoutSize;
					cudaError_t err3 = cudaHostAlloc((void**)&h_rcout, rcoutHostSize, cudaHostAllocDefault);
					if (err3 != cudaSuccess)
					{
						printf("Error! Failed to allocate host memory for receiver output buffer on device %d\n",device);
						h_rcout = 0L;
						retval = 0;
					}	
					else
					{
						printf("Allocated %ld bytes of host memory for receiver output buffer on device %d\n",rcoutHostSize,device);
					}
					h_rcoutflop = h_rcout;
					h_rcoutflap = h_rcout + (2 * num_timesteps * maxbincnt);

					if (retval != 0)
					{
						d_rcx = new int*[nx];
						d_rcy = new int*[nx];
						d_rcz = new int*[nx];

						h_rcx = new int*[nx];
						h_rcy = new int*[nx];
						h_rcz = new int*[nx];
						h_idx = new int*[nx];

						h_recnum = new int[nx];
						h_recnx = nx;

						h_rcoutstatus = new int[12*num_timesteps];
						h_rcoutstatusflip = h_rcoutstatus;
						h_rcoutstatusflop = h_rcoutstatus + 4*num_timesteps;
						h_rcoutstatusflap = h_rcoutstatus + 8*num_timesteps;
						for (int i = 0;  i < 12*num_timesteps;  ++i)
						{
							h_rcoutstatus[i] = -1;
						}

						unsigned long boffset1 = 0L;
						unsigned long boffset2 = 0L;
						for (int iX = 0;  iX < nx;  ++iX)
						{
							int x0 = iX * 16;
							int x1 = x0 + 15;

							int cnt = 0;
							for (int i = 0;  i < arg_recnum;  ++i)
							{
								int x = arg_rcx[i];
								int y = arg_rcy[i];
								if (x >= x0 && x <= x1 && y >= y0 && y <= y1)
								{
									++cnt;
								}
							}
							if (cnt > 0)
							{
								printf("bin %d :: cnt = %d\n",iX,cnt);

								d_rcx[iX] = d_rcbuf + boffset1;			boffset1 += (unsigned long)cnt;
								d_rcy[iX] = d_rcbuf + boffset1;			boffset1 += (unsigned long)cnt;
								d_rcz[iX] = d_rcbuf + boffset1;			boffset1 += (unsigned long)cnt;

								h_rcx[iX] = h_rcbuf + boffset2;			boffset2 += (unsigned long)cnt;
								h_rcy[iX] = h_rcbuf + boffset2;			boffset2 += (unsigned long)cnt;
								h_rcz[iX] = h_rcbuf + boffset2;			boffset2 += (unsigned long)cnt;
								h_idx[iX] = new int[cnt];

								h_recnum[iX] = cnt;

								for (int i = 0, idx = 0;  i < arg_recnum;  ++i)
								{
									int x = arg_rcx[i];
									int y = arg_rcy[i];
									int z = arg_rcz[i];
									if (x >= x0 && x <= x1 && y >= y0 && y <= y1)
									{
										h_rcx[iX][idx] = x;
										h_rcy[iX][idx] = y;
										h_rcz[iX][idx] = z;
										h_idx[iX][idx] = i;
										++idx;
									}
								}
							}
							else
							{
								d_rcx[iX] = 0L;
								d_rcy[iX] = 0L;
								d_rcz[iX] = 0L;

								h_rcx[iX] = 0L;
								h_rcy[iX] = 0L;
								h_rcz[iX] = 0L;

								h_recnum[iX] = 0;
							}
						}

						cudaMemcpy(d_rcbuf, h_rcbuf, (unsigned long)totcnt * 12, cudaMemcpyHostToDevice);
					}
				}
			}
		}

		delete [] bincnt;
		if (retval == 0)
		{
			printf("Calling Free_Receiver_Locations();\n");
			Free_Receiver_Locations();
		}

		printf("Exiting Put_Receiver_Locations()\n");
		return retval;
	}

	int Put_Constant_Sponges(
		int dimx,
		int dimy,
		int dimz
		)
	{
		int max = dimx;
		if (dimy > max) max = dimy;
		if (dimz > max) max = dimz;
		float* buf = new float[max];
		for (int i = 0;  i < max;  ++i) buf[i] = 1.0f;
		int retval = Put_Sponges(buf,dimx,buf,dimy,buf,dimz);
		delete [] buf;
		return retval;
	}

	int Put_Sponges(
		float* h_spg_x,
		int dimx,
		float* h_spg_y,
		int dimy,
		float* h_spg_z,
		int dimz
		)
	{
		int retval = 1;

		cudaSetDevice(device);

		if (h_spg_x != 0L)
		{
			// free X sponge in device memory
			if (d_spg_x != 0L) cudaFree(d_spg_x);

			// allocate X sponge in device memory and copy over
			unsigned long spgxSize = (unsigned long)dimx * (unsigned long)sizeof(float);
			unsigned long spgxSizePad = (unsigned long)16 * (unsigned long)sizeof(float);
			cudaError_t err1 = cudaMalloc((void**)&d_spg_x, spgxSize + spgxSizePad);
			if (err1 != cudaSuccess)
			{
				printf("Error! Failed to allocate device memory for sponge_x buffer on device %d\n",device);
				d_spg_x = 0L;
				retval = 0;
			}
			else
			{
				printf("Allocated %ld bytes of device memory for sponge_x buffer on device %d.\n",spgxSize+spgxSizePad,device);
				cudaMemcpy(d_spg_x, h_spg_x, spgxSize, cudaMemcpyHostToDevice);
				cudaMemcpy(d_spg_x+dimx, h_spg_x, spgxSizePad, cudaMemcpyHostToDevice);
				// note to self: Copying the padding with a separate cudaMemcpy is really inefficient, but only done once.
			}
		}
	
		if (h_spg_y != 0L && retval != 0)
		{
			// free Y sponge in device memory
			if (d_spg_y != 0L) cudaFree(d_spg_y);

			// allocate Y sponge in device memory and copy over
			unsigned long spgySize = (unsigned long)dimy * (unsigned long)sizeof(float);
			cudaError_t err1 = cudaMalloc((void**)&d_spg_y, spgySize);
			if (err1 != cudaSuccess)
			{
				printf("Error! Failed to allocate device memory for sponge_y buffer on device %d\n",device);
				d_spg_y = 0L;
				retval = 0;
			}
			else
			{
				printf("Allocated %ld bytes of device memory for sponge_y buffer on device %d.\n",spgySize,device);
				cudaMemcpy(d_spg_y, h_spg_y, spgySize, cudaMemcpyHostToDevice);
			}
		}

		if (h_spg_z != 0L && retval != 0)
		{
			// free Z sponge in device memory
			if (d_spg_z != 0L) cudaFree(d_spg_z);

			// allocate Z sponge in device memory and copy over
			unsigned long spgzSize = (unsigned long)dimz * (unsigned long)sizeof(float);
			cudaError_t err1 = cudaMalloc((void**)&d_spg_z, spgzSize);
			if (err1 != cudaSuccess)
			{
				printf("Error! Failed to allocate device memory for sponge_z buffer on device %d\n",device);
				d_spg_z = 0L;
				retval = 0;
			}
			else
			{
				printf("Allocated %ld bytes of device memory for sponge_z buffer on device %d.\n",spgzSize,device);
				cudaMemcpy(d_spg_z, h_spg_z, spgzSize, cudaMemcpyHostToDevice);
			}
		}

		return retval;
	}

	int Allocate_Device_Resources(int dimx, int dimy, int dimz)
	{
		int retval = 1;
		cudaSetDevice(device);
		cudaDeviceProp devProps;
		cudaGetDeviceProperties(&devProps, device);
		cc_major = devProps.major;
		cc_minor = devProps.minor;	
		canAccessPeer = 0;
		if (Succeeding_Propagator != 0L)
		{
			cudaError_t err = cudaDeviceCanAccessPeer(&canAccessPeer, device, Succeeding_Propagator->Get_Device_ID());
			if (canAccessPeer)
			{
				cudaDeviceEnablePeerAccess(Succeeding_Propagator->Get_Device_ID(), 0);
				printf("Peer-2-peer access enabled between device %d and %d\n",device,Succeeding_Propagator->Get_Device_ID());
			}
			else
			{
				printf("No peer-2-peer access possible between device %d and %d\n",device,Succeeding_Propagator->Get_Device_ID());
			}
		}
		for (int i = 0;  i < num_layers;  ++i)
		{
			if (layers[i]->Allocate_Device_Resources() == 0)
			{
				retval = 0;
				break;
			}
		}
		if (retval != 0)
		{
			cudaError_t err1 = cudaStreamCreate(&compute_stream);
			if (err1 != cudaSuccess)
			{
				compute_stream = 0;
				retval = 0;
			}
			else
			{
				cudaError_t err2 = cudaStreamCreate(&input_stream);
				if (err2 != cudaSuccess)
				{
					input_stream = 0;
					retval = 0;
				}
				else
				{
					cudaError_t err3 = cudaStreamCreate(&output_stream);
					if (err3 != cudaSuccess)
					{
						output_stream = 0;
						retval = 0;
					}
				}
			}
		}
		if (retval != 0)
		{
			// allocate temporary buffers.
			if (layers != 0L)
			{
				unsigned long tempSize1 = Compute_Storage_Requirements_of_Temporary_Buffer(128);
				cudaError_t err1 = cudaMalloc((void**)&d_temp1, tempSize1);
				if (err1 != cudaSuccess)
				{
					printf("Error! Failed to allocate device memory for temporary buffer #1 on device %d\n",device);
					d_temp1 = 0L;
					retval = 0;
				}
				else
				{
					printf("Allocated %ld bytes of device memory for temporary buffer #1 on device %d.\n",tempSize1,device);
					cudaMemset(d_temp1, 0, tempSize1);

					unsigned long tempSize2 = tempSize1;
					cudaError_t err2 = cudaMalloc((void**)&d_temp2, tempSize2);
					if (err2 != cudaSuccess)
					{
						printf("Error! Failed to allocate device memory for temporary buffer #2 on device %d\n",device);
						d_temp2 = 0L;
						retval = 0;
					}
					else
					{
						printf("Allocated %ld bytes of device memory for temporary buffer #2 on device %d.\n",tempSize2,device);
						cudaMemset(d_temp2, 0, tempSize2);

						unsigned long tempSize3 = Compute_Storage_Requirements_of_Temporary_Buffer(96);
						cudaError_t err3 = cudaMalloc((void**)&d_temp3, tempSize3);
						if (err3 != cudaSuccess)
						{
							printf("Error! Failed to allocate device memory for temporary buffer #3 on device %d\n",device);
							d_temp3 = 0L;
							retval = 0;
						}
						else
						{
							printf("Allocated %ld bytes of device memory for temporary buffer #3 on device %d.\n",tempSize3,device);
							cudaMemset(d_temp3, 0, tempSize3);
						}
					}
				}
			}
		}
		if (retval != 0)
		{
			retval = Put_Constant_Sponges(dimx,dimy,dimz);
		}
		return retval;
	}

	// Get total number of bytes required for all temporary buffers.
	unsigned long Compute_Storage_Requirements_of_Temporary_Buffer()
	{
		return Compute_Storage_Requirements_of_Temporary_Buffer(128) + 
			Compute_Storage_Requirements_of_Temporary_Buffer(128) +
			Compute_Storage_Requirements_of_Temporary_Buffer(96);
	}

	unsigned long Compute_Storage_Requirements_of_Temporary_Buffer(int floats_per_cell)
	{
		int max_compute_dim_y, max_compute_dim_z;
		Find_Largest_Compute_Dimensions(max_compute_dim_y, max_compute_dim_z);
		printf("Largest compute dimensions : 16 by %d by %d\n",max_compute_dim_y,max_compute_dim_z);
		return floats_per_cell * sizeof(float) * (max_compute_dim_y+18) * max_compute_dim_z;
	}

	void Find_Largest_Compute_Dimensions(int& max_compute_dim_y, int& max_compute_dim_z)
	{
		max_compute_dim_y = 0;
		max_compute_dim_z = 0;
		for (int i = 0;  i < num_layers;  ++i)
		{
			if (layers[i]->is_compute_layer)
			{
				if (layers[i]->compute_dim_y > max_compute_dim_y) max_compute_dim_y = layers[i]->compute_dim_y;
				if (layers[i]->compute_dim_z > max_compute_dim_z) max_compute_dim_z = layers[i]->compute_dim_z;
			}
		}
	}

	//
	// Computational overhead is a factor yielding how much extra work the halos add to the compute kernel launches.
	// A factor of 1.0 is optimal and says NO overhead was added by halos.
	// Returns 0.0 if layers have not been instantiated yet.
	//
	double Compute_Computational_Overhead()
	{
		if (num_layers > 0)
		{
			int net_compute_dim_x, net_compute_dim_y, net_compute_dim_z;
			Get_Net_Compute_Dimensions(net_compute_dim_x, net_compute_dim_y, net_compute_dim_z);
			unsigned long netCompute = (unsigned long)net_compute_dim_y * (unsigned long)net_compute_dim_z;
			unsigned long cellAcc = 0, num_compute_layers = 0;
			for (int i = 0;  i < num_layers;  ++i)
			{
				if (layers[i]->is_compute_layer)
				{
					++num_compute_layers;
					cellAcc += (unsigned long)layers[i]->compute_dim_y * (unsigned long)layers[i]->compute_dim_z;
				}
			}
			netCompute *= num_compute_layers;
			return (double)cellAcc / (double)netCompute;
		}
		else
		{
			return 0.0;
		}
	}

	double Estimate_Computational_Cost()
	{
		double cost = 0.0;
		for (int i = 0;  i < num_layers;  ++i)
		{
			if (layers[i]->is_compute_layer)
			{
				cost += (double)layers[i]->compute_dim_y * (double)layers[i]->compute_dim_z;
			}
		}
		return cost;
	}

	void Print()
	{
		if (layers != 0L)
		{
			printf("device %d :: %s\n",device,device_name);
			for (int i = 0;  i < num_layers;  ++i) layers[i]->Print();
			printf("Output X Offset is %d\n",Get_Output_X_Offset());
			printf("%.2f GB device memory required, computational overhead is %.2f, bytes per cell per timestep is %.2f\n",
				Compute_Memory_Required_GB(),
				Compute_Computational_Overhead(),
				Compute_Bytes_Per_Cell_Per_Timestep()
				);
		}
		else
		{
			printf("Error! SingleDevicePropagator object not initialized properly!\n");
		}
	}

	void Compute_Y_Extent(
		int pipeline_net_y0,
		int pipeline_net_y1,
		int layer_id,
		int& layer_y0,
		int& layer_y1
		)
	{
		TimeStepLayer* layer = Get_Layer(layer_id);
		if (layer != 0L)
		{
			layer->Compute_Y_Extent(pipeline_net_y0,pipeline_net_y1,layer_y0,layer_y1);
		}
		else
		{
			layer_y0 = 0;
			layer_y1 = 0;
		}
	}

	void Compute_Source_Parameters(
		int pipeline_net_y0,
                int pipeline_net_y1,
		int dimy,
		int layer_id,
		int stripe_id,
		int& src_y0,
		int& src_y1,
		int& src_off_y,
		void*& d_src
		)
	{
		TimeStepLayer* layer = Get_Layer(layer_id);
                if (layer != 0L)
                {
			layer->Compute_Source_Parameters(pipeline_net_y0,pipeline_net_y1,dimy,stripe_id,src_y0,src_y1,src_off_y,d_src);
		}
		else
		{
			src_y0 = 0;
			src_y1 = 0;
			src_off_y = 0;
			d_src = 0L;
		}
	}

	void* Get_Stripe(
		int layer_id,
		int stripe_id
		)
	{
		TimeStepLayer* layer = Get_Layer(layer_id);
		if (layer != 0L)
		{
			return layer->Get_Stripe(stripe_id);
		}
		else
		{
			return 0L;
		}
	}

	TimeStepLayer* Get_Layer(int layer_id)
	{
		if (layer_id >= -num_layers && layer_id < num_layers)
		{
			int id = (layer_id + num_layers) % num_layers;
			return layers[id];
		}
		return 0L;
	}

	int Get_Number_of_Layers()
	{
		return num_layers;
	}

	int Get_Device_ID()
	{
		return device;
	}

	int Get_Number_of_Timesteps()
	{
		return num_timesteps;
	}

	int Get_X_Delay()
	{
		return Get_Output_X_Offset() / -16;
	}

	int Get_Output_X_Offset()
	{
		// Same offset for OTflag 2 and 4.
		return layers[num_layers-1]->compute_off_x-32;
	}

	void Get_Offsets(int& compute_off_x, int& input_off_y)
	{
		compute_off_x = layers[num_layers-1]->compute_off_x;
		input_off_y = layers[2]->input_off_y;
	}

	void Get_Net_Compute_Dimensions(int& net_compute_dim_x, int& net_compute_dim_y, int& net_compute_dim_z)
	{
		if (Succeeding_Propagator != 0L)
		{
			Succeeding_Propagator->Get_Net_Compute_Dimensions(net_compute_dim_x,net_compute_dim_y,net_compute_dim_z);
		}
		else
		{
			net_compute_dim_x = 16;
			net_compute_dim_y = bY;
			net_compute_dim_z = bZ;
		}
	}

	unsigned long Compute_Memory_Required()
	{
		unsigned long memReq = 0;
		for (int i = 0;  i < num_layers;  ++i)
		{
			memReq += layers[i]->Compute_Storage_Requirement();
		}
		return memReq + Compute_Storage_Requirements_of_Temporary_Buffer();
	}
	
	double Compute_Memory_Required_GB()
	{
		return (double)Compute_Memory_Required() / 1073741824.0;
	}

	double Compute_Bytes_Per_Cell_Per_Timestep()
	{
		int net_compute_dim_x, net_compute_dim_y, net_compute_dim_z;
		Get_Net_Compute_Dimensions(net_compute_dim_x,net_compute_dim_y,net_compute_dim_z);

		int model_y0 = layers[0]->input_off_y;
		int model_y1 = layers[0]->output_dim_y + layers[0]->input_off_y;
		int m1_y0 = layers[1]->input_off_y;
		int m1_y1 = layers[1]->output_dim_y + layers[1]->input_off_y;
		int p0_y0 = layers[2]->input_off_y;
		int p0_y1 = layers[2]->output_dim_y + layers[2]->input_off_y;

		//printf("model y0=%d, y1=%d\n",model_y0,model_y1);
		//printf("m1 y0=%d, y1=%d\n",m1_y0,m1_y1);
		//printf("p0 y0=%d, y1=%d\n",p0_y0,p0_y1);
		
		if (exclude_halo_ylo)
		{
			if (model_y0 < 0) model_y0 = 0;
			if (m1_y0 < 0) m1_y0 = 0;
			if (p0_y0 < 0) p0_y0 = 0;
		}
		
		//printf("model y0=%d, y1=%d\n",model_y0,model_y1);
		//printf("m1 y0=%d, y1=%d\n",m1_y0,m1_y1);
		//printf("p0 y0=%d, y1=%d\n",p0_y0,p0_y1);

		if (exclude_halo_yhi)
		{
			if (model_y1 > net_compute_dim_y) model_y1 = net_compute_dim_y;
			if (m1_y1 > net_compute_dim_y) m1_y1 = net_compute_dim_y;
			if (p0_y1 > net_compute_dim_y) p0_y1 = net_compute_dim_y;
		}

		//printf("model y0=%d, y1=%d\n",model_y0,model_y1);
		//printf("m1 y0=%d, y1=%d\n",m1_y0,m1_y1);
		//printf("p0 y0=%d, y1=%d\n",p0_y0,p0_y1);

		int net_input_model_y = model_y1 - model_y0;
		int net_input_m1_y = m1_y1 - m1_y0;
		int net_input_p0_y = p0_y1 - p0_y0;

		double bytes_per_cell_per_timestep = 
			((double)net_input_model_y * (double)layers[0]->bytes_per_cell + 
			 (double)net_input_m1_y * (double)layers[1]->bytes_per_cell +
			 (double)net_input_p0_y * (double)layers[2]->bytes_per_cell) / 
			((double)net_compute_dim_y * (double)num_timesteps);
		return bytes_per_cell_per_timestep;
	}

	void CopyFromHostToDevice(
		int dst_layer_id,
		int dst_stripe_id,
		void* src_stripe,
		int pipeline_net_y0,
		int pipeline_net_y1,
		int dimy
		)
	{
		int dst_y0, dst_y1;
		Compute_Y_Extent(pipeline_net_y0,pipeline_net_y1,dst_layer_id,dst_y0,dst_y1);

		// clip against host volume dimensions

		int dst_off_y = 0;
		if (dst_y0 < 0)
		{
			dst_off_y = -dst_y0;
			dst_y0 = 0;
		}
		if (dst_y1 >= dimy)
		{
			dst_y1 = dimy - 1;
		}
		int dst_dim_y = dst_y1 - dst_y0 + 1;

		int src_off_y = dst_y0;

		cudaSetDevice(device);
                float* d_dst = (float*)Get_Stripe(dst_layer_id,dst_stripe_id) + (unsigned long)32 * (unsigned long)dst_off_y * (unsigned long)bZ;
                float* h_src = (float*)src_stripe + (unsigned long)32 * (unsigned long)src_off_y * (unsigned long)bZ;
                unsigned long count = (unsigned long)128 * (unsigned long)dst_dim_y * (unsigned long)bZ;
                cudaMemcpyAsync(d_dst,h_src,count,cudaMemcpyHostToDevice,input_stream);	
		h2d_byte_count += count;

		//printf("CopyFromHostToDevice :: dst_device=%d, dst_layer=%d, dst_stripe_id=%d, y0=%d, y1=%d, src_off_y=%d, dst_off_dy=%d, dst_dim_y=%d, count=%ld\n",device,dst_layer_id,dst_stripe_id,pipeline_net_y0,pipeline_net_y1,src_off_y,dst_off_y,dst_dim_y,count);	
	}
		
	void CopyFromDeviceToHost(
		int src_layer_id,
		int src_stripe_id,
		void* dst_stripe,
		int pipeline_net_y0,
		int pipeline_net_y1
		)
	{
		int src_y0, src_y1;
		Compute_Y_Extent(pipeline_net_y0,pipeline_net_y1,src_layer_id,src_y0,src_y1);

		int src_off_y = 0;
		if (src_y0 < pipeline_net_y0)
		{
			src_off_y = (pipeline_net_y0 - src_y0);
			src_y0 = pipeline_net_y0;
		}
		if (src_y1 > pipeline_net_y1)
		{
			src_y1 = pipeline_net_y1;
		}
		int src_dim_y = src_y1 - src_y0 + 1;

                cudaSetDevice(device);
                float* h_dst = (float*)dst_stripe + (unsigned long)32 * (unsigned long)src_y0 * (unsigned long)bZ;
                float* d_src = (float*)Get_Stripe(src_layer_id,src_stripe_id) + (unsigned long)32 * (unsigned long)src_off_y * (unsigned long)bZ;
                unsigned long count = (unsigned long)128 * (unsigned long)src_dim_y * (unsigned long)bZ;
                cudaMemcpyAsync(h_dst,d_src,count,cudaMemcpyDeviceToHost,output_stream);
		d2h_byte_count += count;

		//printf("CopyFromDeviceToHost :: src_layer_id=%d, src_stripe_id=%d, y0=%d, y1=%d, src_off_y=%d, src_dim_y=%d, dst_off_dy=%d, count=%ld\n",src_layer_id,src_stripe_id,pipeline_net_y0,pipeline_net_y1,src_off_y,src_dim_y,src_y0,count);
	}

	void CopyFromDeviceToDevice(
		SingleDevicePropagator* src_prop,
		int src_layer_id,
		int src_stripe_id,
		int dst_layer_id,
		int dst_stripe_id,
		int pipeline_net_y0,
                int pipeline_net_y1,
		int dimy
		)
	{
		int dst_y0, dst_y1;
		Compute_Y_Extent(pipeline_net_y0,pipeline_net_y1,dst_layer_id,dst_y0,dst_y1);

		// clip destination against global volume 

		int dst_off_y = 0;
		if (dst_y0 < 0)
		{
			dst_off_y = -dst_y0;
			dst_y0 = 0;
		}
		if (dst_y1 >= dimy)
		{
			dst_y1 = dimy - 1;
		}

		// get source parameters

		void* src_stripe = src_prop->Get_Stripe(src_layer_id,src_stripe_id);

		int src_y0, src_y1;
		src_prop->Compute_Y_Extent(pipeline_net_y0,pipeline_net_y1,src_layer_id,src_y0,src_y1);

		// clip source against destination

		int src_off_y = 0;
                if (src_y0 < dst_y0)
                {
			src_off_y = (dst_y0 - src_y0);
			src_y0 = dst_y0;
                }
                if (src_y1 > dst_y1)
                {
                        src_y1 = dst_y1;
                }
		int src_dim_y = src_y1 - src_y0 + 1;
		
		//printf("CopyFromDeviceToDevice :: src_device=%d, src_layer_id=%d, src_stripe_id=%d, dst_device=%d, dst_layer_id=%d, dst_stripe_id=%d, y0=%d, y1=%d, src_off_y=%d, src_dim_y=%d, dst_off_dy=%d, src_y0=%d, dst_y0=%d\n",src_prop->Get_Device_ID(),src_layer_id,src_stripe_id,device,dst_layer_id,dst_stripe_id,pipeline_net_y0,pipeline_net_y1,src_off_y,src_dim_y,dst_off_y,src_y0,dst_y0);

		cudaSetDevice(src_prop->Get_Device_ID());
		float* d_dst = (float*)Get_Stripe(dst_layer_id,dst_stripe_id) + (unsigned long)32 * (unsigned long)dst_off_y * (unsigned long)bZ;
		float* d_src = (float*)src_stripe + (unsigned long)32 * (unsigned long)src_off_y * (unsigned long)bZ;
		unsigned long count = (unsigned long)128 * (unsigned long)src_dim_y * (unsigned long)bZ;
		cudaMemcpyPeerAsync(d_dst,device,d_src,src_prop->Get_Device_ID(),count,src_prop->output_stream);
	}

	void Start_Receiver_Data_Flow()
	{
		if (h_recnx > 0)
		{
			cudaSetDevice(device);
			for (int i = 0;  i < num_timesteps;  ++i)
			{
				int stripe_1 = h_rcoutstatusflop[4*i];
				int timestep_1 = h_rcoutstatusflop[4*i+1];
				if (stripe_1 >= 0 && timestep_1 >= 0)
				{
					//printf("Asynchronous transfer of receiver values for stripe %d timestep %d\n",stripe_1,timestep_1);
					float* src = d_rcoutflop + (unsigned long)(2*i) * h_rcoutmaxbinsize;
					float* dst = h_rcoutflop + (unsigned long)(2*i) * h_rcoutmaxbinsize;
					unsigned long count = h_recnum[stripe_1] * 4;
					cudaMemcpyAsync(dst,src,count,cudaMemcpyDeviceToHost,output_stream);
				}

				int stripe_2 = h_rcoutstatusflop[4*i+2];
				int timestep_2 = h_rcoutstatusflop[4*i+3];
				if (stripe_2 >= 0 && timestep_2 >= 0)
				{
					//printf("Asynchronous transfer of receiver values for stripe %d timestep %d\n",stripe_2,timestep_2);
					float* src = d_rcoutflop + (unsigned long)(2*i+1) * h_rcoutmaxbinsize;
					float* dst = h_rcoutflop + (unsigned long)(2*i+1) * h_rcoutmaxbinsize;
					unsigned long count = h_recnum[stripe_2] * 4;
					cudaMemcpyAsync(dst,src,count,cudaMemcpyDeviceToHost,output_stream);
				}
			}
		}
	}

	void Start_Data_Flow(
		void* h_em_inp1, 
		void* h_em_inp2, 
		void* h_pq_inp1, 
		void* h_pq_inp2, 
		void* h_prev_pq_inp1, 
		void* h_prev_pq_inp2, 
		void* h_pq_out1, 
		void* h_pq_out2, 
		void* h_prev_pq_out1,
		void* h_prev_pq_out2,
		int y0,
		int y1,
		int dimy,
		int is_first,
		int is_last,
		SingleDevicePropagator* preceeding_propagator
		)
	{
		if (is_first)
		{
			//printf("FIRST!\n");
			//Print();

			// first device, so input from host
			CopyFromHostToDevice(0,-2,h_em_inp1,y0,y1,dimy);
			CopyFromHostToDevice(1,-2,h_prev_pq_inp1,y0,y1,dimy);
			CopyFromHostToDevice(2,-2,h_pq_inp1,y0,y1,dimy);

			CopyFromHostToDevice(0,-1,h_em_inp2,y0,y1,dimy);
			CopyFromHostToDevice(1,-1,h_prev_pq_inp2,y0,y1,dimy);
			CopyFromHostToDevice(2,-1,h_pq_inp2,y0,y1,dimy);
		}
		else
		{
			//printf("P2P!\n");
			//Print();

			// peer-2-peer transfer
			int src_num_layers = preceeding_propagator->Get_Number_of_Layers();
			int src_m1_layer_id = OTflag == 2 ? src_num_layers-2 : src_num_layers - 3;
			int src_p0_layer_id = src_num_layers - 1;

			CopyFromDeviceToDevice(preceeding_propagator,0              ,0,0,-2,y0,y1,dimy);   // earth model
			CopyFromDeviceToDevice(preceeding_propagator,src_m1_layer_id,0,1,-2,y0,y1,dimy);
			CopyFromDeviceToDevice(preceeding_propagator,src_p0_layer_id,0,2,-2,y0,y1,dimy);

			CopyFromDeviceToDevice(preceeding_propagator,0              ,1,0,-1,y0,y1,dimy);   // earth model
			CopyFromDeviceToDevice(preceeding_propagator,src_m1_layer_id,1,1,-1,y0,y1,dimy);
			CopyFromDeviceToDevice(preceeding_propagator,src_p0_layer_id,1,2,-1,y0,y1,dimy);
		}
		if (is_last)
		{
			//printf("LAST!\n");
			//Print();

			// last device, so also output final result to host,
			// but only after pipeline fills up.

			int src_num_layers = Get_Number_of_Layers();
			int src_p0_layer_id = OTflag == 2 ? src_num_layers - 2 : src_num_layers - 3;
			int src_p1_layer_id = src_num_layers - 1;

			if (h_pq_out1 != 0L && h_prev_pq_out1 != 0L)
			{
				CopyFromDeviceToHost(src_p0_layer_id,0,h_prev_pq_out1,y0,y1);
				CopyFromDeviceToHost(src_p1_layer_id,0,h_pq_out1,y0,y1);
			}

			if (h_pq_out2 != 0L && h_prev_pq_out2 != 0L)
			{
				CopyFromDeviceToHost(src_p0_layer_id,1,h_prev_pq_out2,y0,y1);
				CopyFromDeviceToHost(src_p1_layer_id,1,h_pq_out2,y0,y1);
			}
		}
	}

	//
	// Copy inputs into outputs. Used for debugging to verify data flow.
	//
	void LaunchSimpleCopyKernel()
	{
		cudaSetDevice(device);
		for (int i = 0;  i < num_timesteps;  ++i)
		{
			if (OTflag == 2)
			{
				int src_layer = 2 + i;
				int dst_layer = src_layer + 1;

				LaunchSimpleCopyKernelSingleLayer(i == 0 ? -5 : -3, src_layer, -2, dst_layer);
				LaunchSimpleCopyKernelSingleLayer(i == 0 ? -4 : -2, src_layer, -1, dst_layer);
			}
			else if (OTflag == 4)
			{
				// Apq[+1]
				int Apq_src_layer = 2 + i * 2;
				int Apq_dst_layer = Apq_src_layer + 1;

				LaunchSimpleCopyKernelSingleLayer(i == 0 ? -5 : -3, Apq_src_layer, -2, Apq_dst_layer);
				LaunchSimpleCopyKernelSingleLayer(i == 0 ? -4 : -2, Apq_src_layer, -1, Apq_dst_layer);

				// pq[+1]
				int pq_dst_layer = Apq_dst_layer + 1;

				LaunchSimpleCopyKernelSingleLayer(-3, Apq_dst_layer, -2, pq_dst_layer);
				LaunchSimpleCopyKernelSingleLayer(-2, Apq_dst_layer, -1, pq_dst_layer);
			}
		}
	}

	void CollectReceiverValues(float* muxed_traces)
	{
		if (h_recnx > 0)
		{
			for (int i = 0;  i < num_timesteps;  ++i)
			{
#pragma omp parallel sections num_threads(2)
				{
#pragma omp section
				{
				int stripe_1 = h_rcoutstatusflap[4*i];
				if (stripe_1 >= 0)
				{
					int imgstep_1 = h_rcoutstatusflap[4*i+1];
					unsigned long offset_1 = (unsigned long)h_totrecnx * (unsigned long)imgstep_1;
					float* mux_1 = muxed_traces + offset_1;
					int* h_idx_1 = h_idx[stripe_1];
					float* h_out_1 = h_rcoutflap + (2*i)*h_rcoutmaxbinsize;
					int recnum_1 = h_recnum[stripe_1];

					if (recnum_1 > 0)
					{
						//printf("Collecting receiver values for stripe %d timestep %d - recnum_1=%d, totrecnx=%d, offset_1=%ld\n",stripe_1,imgstep_1,recnum_1,h_totrecnx,offset_1);
						for (int j = 0;  j < recnum_1;  ++j)
						{
							int dst_idx = h_idx_1[j];
							mux_1[dst_idx] = h_out_1[j];
							//printf("...%d <- %d -- %e\n",dst_idx,j,h_out_1[j]);
						}
					}
				}
				}
#pragma omp section
				{
				int stripe_2 = h_rcoutstatusflap[4*i+2];
				if (stripe_2 >= 0)
				{
					int imgstep_2 = h_rcoutstatusflap[4*i+3];
					unsigned long offset_2 = (unsigned long)h_totrecnx * (unsigned long)imgstep_2;
					float* mux_2 = muxed_traces + offset_2;
					int* h_idx_2 = h_idx[stripe_2];
					float* h_out_2 = h_rcoutflap + (2*i+1)*h_rcoutmaxbinsize;
					int recnum_2 = h_recnum[stripe_2];

					if (recnum_2 > 0)
					{
						//printf("Collecting receiver values for stripe %d timestep %d - recnum_2=%d, totrecnx=%d, offset_2=%ld\n",stripe_2,imgstep_2,recnum_2,h_totrecnx,offset_2);
						for (int j = 0;  j < recnum_2;  ++j)
						{
							int dst_idx = h_idx_2[j];
							mux_2[dst_idx] = h_out_2[j];
							//printf("...%d <- %d -- %e\n",dst_idx,j,h_out_2[j]);
						}
					}
				}
				}
				}
			}
		}
	}

	void FlipFlopFlap()
	{
		float* tmp = h_rcoutflap;
		h_rcoutflap = h_rcoutflop;
		h_rcoutflop = tmp;

		tmp = d_rcoutflop;
		d_rcoutflop = d_rcoutflip;
		d_rcoutflip = tmp;

		int* itmp = h_rcoutstatusflap;
		h_rcoutstatusflap = h_rcoutstatusflop;
		h_rcoutstatusflop = h_rcoutstatusflip;
		h_rcoutstatusflip = itmp;
	}

	void LaunchPropagationKernel(int logLevel, int* mux_imgstep, float* mux_time_delta, int receiver_num_timesteps)
	{
		cudaSetDevice(device);
		for (int i = 0;  i < num_timesteps;  ++i)
		{
			if (OTflag == 2)
			{
				int src_layer = 2 + i;

				// check if inputs are valid (i.e. pipeline has reached this layer)
				if (layers[src_layer+1]->curr_input_timestep_2 > 0)
				{
					/*
					int left_edge = layers[src_layer+1]->curr_input_stripe_number_1 == 0 ? 1 : 0;
					int right_edge = layers[src_layer+1]->curr_input_stripe_number_2 == (total_num_stripes - 1) ? 1 : 0;
					if (left_edge) printf("LEFT EDGE\n");
					if (right_edge) printf("RIGHT EDGE\n");
					*/

					const float* d_prev_pq1 = (const float*)(layers[src_layer-1]->Get_Stripe(0));
					const float* d_prev_pq2 = (const float*)(layers[src_layer-1]->Get_Stripe(1));

					const float* d_curr_pq0 = (const float*)(layers[src_layer]->Get_Stripe(i==0?-6:-4));
					const float* d_curr_pq1 = (const float*)(layers[src_layer]->Get_Stripe(i==0?-5:-3));
					const float* d_curr_pq2 = (const float*)(layers[src_layer]->Get_Stripe(i==0?-4:-2));
					const float* d_curr_pq3 = (const float*)(layers[src_layer]->Get_Stripe(i==0?-3:-1));

					float* d_next_pq1 = (float*)(layers[src_layer+1]->Get_Stripe(-2));
					float* d_next_pq2 = (float*)(layers[src_layer+1]->Get_Stripe(-1));

					const int* d_em0 = (const int*)(layers[0]->Get_Stripe(-(6+i)));
					const int* d_em1 = (const int*)(layers[0]->Get_Stripe(-(5+i)));
					const int* d_em2 = (const int*)(layers[0]->Get_Stripe(-(4+i)));
					const int* d_em3 = (const int*)(layers[0]->Get_Stripe(-(3+i)));

					int dimy = layers[src_layer+1]->compute_dim_y;
					int dimz = layers[src_layer+1]->compute_dim_z;
				
					int   em_y0 = layers[0]->input_off_y;
					int prev_y0 = layers[src_layer-1]->is_compute_layer ? layers[src_layer-1]->compute_off_y - layers[src_layer-1]->compute_out_y : layers[src_layer-1]->input_off_y;
					int curr_y0 = layers[src_layer  ]->is_compute_layer ? layers[src_layer  ]->compute_off_y - layers[src_layer  ]->compute_out_y : layers[src_layer  ]->input_off_y;
					int next_y0 = layers[src_layer+1]->compute_off_y - 9;

					int   em_y_offset = next_y0 -   em_y0;
					int prev_y_offset = next_y0 - prev_y0;
					int curr_y_offset = next_y0 - curr_y0;
					int next_y_offset = layers[src_layer+1]->compute_out_y;

					// x-y corner of current stripe
					int x0_1 = layers[src_layer+1]->curr_input_stripe_number_1 * 16;
					int x0_2 = layers[src_layer+1]->curr_input_stripe_number_2 * 16;
					int y0 = layers[src_layer+1]->y0 + layers[src_layer+1]->compute_off_y;
					if (logLevel >= 5)
					{
						printf("device=%d, layer=%d, dimy = %d, dimz = %d, em_y_offset = %d, prev_y_offset = %d, curr_y_offset = %d, next_y_offset = %d, x0_1 = %d, x0_2 = %d, y0 = %d\n",device,i,dimy,dimz,em_y_offset,prev_y_offset,curr_y_offset,next_y_offset,x0_1,x0_2,y0);
						fflush(stdout);
					}

					// figure out if source should be injected
					int stripe_xsrc_1 = xsrc - x0_1;
					int stripe_ysrc_1 = ysrc - y0;
					int stripe_zsrc_1, stripe_zsrcghost_1;
					float src_val_1;
					if (src_len > 0 && layers[src_layer+1]->curr_input_timestep_1 < src_len && stripe_xsrc_1 >= 0 && stripe_xsrc_1 < 16 && stripe_ysrc_1 >= 0 && stripe_ysrc_1 < dimy)
					{
						// inject source
						stripe_zsrc_1 = zsrc;
						stripe_zsrcghost_1 = zsrcghost;
						src_val_1 = h_src_wavelet[layers[src_layer+1]->curr_input_timestep_1];
					}
					else
					{
						// don't inject source
						stripe_xsrc_1 = -1;
						stripe_ysrc_1 = -1;
						stripe_zsrc_1 = -1;
						stripe_zsrcghost_1 = -1;
						src_val_1 = 0.0f;
					}

					int stripe_xsrc_2 = xsrc - x0_2;
					int stripe_ysrc_2 = ysrc - y0;
					int stripe_zsrc_2, stripe_zsrcghost_2;
					float src_val_2;
					if (src_len > 0 && layers[src_layer+1]->curr_input_timestep_2 < src_len && stripe_xsrc_2 >= 0 && stripe_xsrc_2 < 16 && stripe_ysrc_2 >= 0 && stripe_ysrc_2 < dimy)
					{
						// inject source
						stripe_zsrc_2 = zsrc;
						stripe_zsrcghost_2 = zsrcghost;
						src_val_2 = h_src_wavelet[layers[src_layer+1]->curr_input_timestep_2];
					}
					else
					{
						// don't inject source
						stripe_xsrc_2 = -1;
						stripe_ysrc_2 = -1;
						stripe_zsrc_2 = -1;
						stripe_zsrcghost_2 = -1;
						src_val_2 = 0.0f;
					}

					// determine if we want receiver outputs
					int *d_rcx_1 = 0L, *d_rcy_1 = 0L, *d_rcz_1 = 0L;
					float *d_out_1 = 0L;
					int recnum_1 = 0;
					float tfrac1 = 0.0f;
					if (h_recnx > 0)
					{
						h_rcoutstatusflip[4*i  ] = -1;	// stripe #
						h_rcoutstatusflip[4*i+1] = -1;	// timestep
						int timestep_1 = layers[src_layer+1]->curr_input_timestep_1;
						if (timestep_1 >= 0 && timestep_1 < receiver_num_timesteps)
						{
							int imgstep_1 = mux_imgstep[timestep_1];
							if (imgstep_1 >= 0)
							{
								// this is an output step
								int stripe_1 = layers[src_layer+1]->curr_input_stripe_number_1;
								recnum_1 = h_recnum[stripe_1];
								if (recnum_1 > 0)
								{
									tfrac1 = mux_time_delta[timestep_1];
									h_rcoutstatusflip[4*i  ] = stripe_1;
									h_rcoutstatusflip[4*i+1] = imgstep_1;
									d_rcx_1 = d_rcx[stripe_1];
									d_rcy_1 = d_rcy[stripe_1];
									d_rcz_1 = d_rcz[stripe_1];
									d_out_1 = d_rcoutflip + (unsigned long)(2*i)*h_rcoutmaxbinsize;
								}
							}
						}
					}
					
					int *d_rcx_2 = 0L, *d_rcy_2 = 0L, *d_rcz_2 = 0L;
					float *d_out_2 = 0L;
					int recnum_2 = 0;
					float tfrac2 = 0.0f;
					if (h_recnx > 0)
					{
						h_rcoutstatusflip[4*i+2] = -1;	// stripe #
						h_rcoutstatusflip[4*i+3] = -1;	// timestep
						int timestep_2 = layers[src_layer+1]->curr_input_timestep_2;
						if (timestep_2 >= 0 && timestep_2 < receiver_num_timesteps)
						{
							int imgstep_2 = mux_imgstep[timestep_2];
							if (imgstep_2 >= 0)
							{
								// this is an output step
								int stripe_2 = layers[src_layer+1]->curr_input_stripe_number_2;
								recnum_2 = h_recnum[stripe_2];
								if (recnum_2 > 0)
								{
									tfrac2 = mux_time_delta[timestep_2];
									h_rcoutstatusflip[4*i+2] = stripe_2;
									h_rcoutstatusflip[4*i+3] = imgstep_2;
									d_rcx_2 = d_rcx[stripe_2];
									d_rcy_2 = d_rcy[stripe_2];
									d_rcz_2 = d_rcz[stripe_2];
									d_out_2 = d_rcoutflip + (unsigned long)(2*i+1)*h_rcoutmaxbinsize;
								}
							}
						}
					}

					TTIDenQ_T2_GPU_Timestep(
							Kernel_Type,
							compute_stream,
							d_prev_pq1,d_prev_pq2,
							d_curr_pq0,d_curr_pq1,d_curr_pq2,d_curr_pq3,
							d_next_pq1,d_next_pq2,
							d_em0,d_em1,d_em2,d_em3,
							d_spg_x,
							d_spg_y,
							d_spg_z,
							prev_y_offset,curr_y_offset,next_y_offset,em_y_offset,
							d_temp1,d_temp2,d_temp3,
							x0_1,x0_2,y0,dimy+18,dimz,
							(float)dt,
							dx_A1,dx_A2,dx_A3,dx_A4,dx_A5,
							dz_A1,dz_A2,dz_A3,dz_A4,dz_A5,
							dz_E1,
							Q_min,Q_scaler,
							density_min,density_scaler,
							dip_min,dip_scaler,
							azm_min,azm_scaler,
							C44C33_min,C44C33_scaler,
							Vel_min,Vel_scaler,
							Del_min,Del_scaler,
							Eps_min,Eps_scaler,
							stripe_xsrc_1,stripe_ysrc_1,stripe_zsrc_1,stripe_zsrcghost_1,src_val_1,
							stripe_xsrc_2,stripe_ysrc_2,stripe_zsrc_2,stripe_zsrcghost_2,src_val_2,
							tfrac1,d_rcx_1,d_rcy_1,d_rcz_1,d_out_1,recnum_1,
							tfrac2,d_rcx_2,d_rcy_2,d_rcz_2,d_out_2,recnum_2,
							_recghost_Flag,_z_freesurface
							);
				}
			}
		}
	}

	void WaitForInputStream()
	{
		cudaSetDevice(device);
		cudaStreamSynchronize(input_stream);
	}

	void WaitForComputeStream()
	{
		cudaSetDevice(device);
		cudaStreamSynchronize(compute_stream);
		cudaError_t err = cudaGetLastError();
		if (err != cudaSuccess)
		{
			printf("CUDA error : %s\n",cudaGetErrorString(err));
		}
	}

	void WaitForOutputStream()
	{
		cudaSetDevice(device);
		cudaStreamSynchronize(output_stream);
	}

	void Move_Sliding_Window(int logLevel)
	{
		for (int i = 0;  i < num_layers;  ++i)
		{
			layers[i]->Move_Sliding_Window(logLevel);
		}
	}

	unsigned long Get_H2D_Byte_Count()
	{
		return h2d_byte_count;
	}

	unsigned long Get_D2H_Byte_Count()
	{
		return d2h_byte_count;
	}

private:
	void LaunchSimpleCopyKernelSingleLayer(int src_stripe, int src_layer, int dst_stripe, int dst_layer)
	{
		float* d_src = (float*)(layers[src_layer]->Get_Stripe(src_stripe));
		float* d_dst = (float*)(layers[dst_layer]->Get_Stripe(dst_stripe));

		dim3 blockShape(32,1,4);
		dim3 gridShape(1,layers[dst_layer]->compute_dim_y,1);

		int y_offset = layers[dst_layer]->compute_out_y;
		int dimz = layers[dst_layer]->compute_dim_z;
		//printf("device=%d dst_layer=%d,'%s' compute_out_y=%d blockShape=%d,%d,%d gridShape=%d,%d,%d compute_dim_y=%d, compute_dim_z=%d\n",device,dst_layer,layers[dst_layer]->name,y_offset,blockShape.x,blockShape.y,blockShape.z,gridShape.x,gridShape.y,gridShape.z,layers[dst_layer]->compute_dim_y,dimz);

		SimpleCopyKernel<<<gridShape,blockShape,0,compute_stream>>>(device,dst_layer,d_dst,d_src,halo_ny,y_offset,dimz);
		//printf("y_offset=%d\n",y_offset);

		/*
		cudaStreamSynchronize(compute_stream);
		cudaError_t err = cudaGetLastError();
		if (err != cudaSuccess)
		{
			printf("CUDA ERROR (device %d, layer %d) : %s\n",device,dst_layer,cudaGetErrorString(err));
		}
		*/
	}

	char* strdup(char* src)
	{
		if (src == 0L) return strdup("nil");
		int len = strlen(src);
		if (len <= 0) return strdup("nil");
		char* dst = new char[len+1];
		for (int i = 0;  i < len;  ++i) dst[i] = src[i];
		dst[len] = 0;
		return dst;
	}

	void Free_Receiver_Locations()
	{
		_recghost_Flag = 0;
		_z_freesurface = 0;
		cudaSetDevice(device);
		if (d_rcbuf != 0L)
		{
			cudaFree(d_rcbuf);
			d_rcbuf = 0L;
		}
		if (d_rcout != 0L)
		{
			cudaFree(d_rcout);
			d_rcout = 0L;
		}
		d_rcoutflip = 0L;
		d_rcoutflop = 0L;
		if (d_rcx != 0L)
		{
			delete [] d_rcx;
			d_rcx = 0L;
		}
		if (d_rcy != 0L)
		{
			delete [] d_rcy;
			d_rcy = 0L;
		}
		if (d_rcz != 0L)
		{
			delete [] d_rcx;
			d_rcz = 0L;
		}
		if (h_rcx != 0L)
		{
			delete [] h_rcx;
			h_rcx = 0L;
		}
		if (h_rcy != 0L)
		{
			delete [] h_rcy;
			h_rcy = 0L;
		}
		if (h_rcz != 0L)
		{
			delete [] h_rcx;
			h_rcz = 0L;
		}
		if (h_idx != 0L)
		{
			for (int i = 0;  i < h_recnx;  ++i)
			{
				if (h_idx[i] != 0L) delete [] h_idx[i];
			}
			delete [] h_idx;
			h_idx = 0L;
		}
		if (h_recnum != 0L)
		{
			delete [] h_recnum;
			h_recnum = 0L;
		}
		if (h_rcout != 0L)
		{
			cudaFreeHost(h_rcout);
			h_rcout = 0L;
		}
		h_rcoutflop = 0L;
		h_rcoutflap = 0L;
		h_rcoutmaxbinsize = 0;
		h_recnx = 0;
		if (h_rcoutstatus != 0L)
		{
			delete [] h_rcoutstatus;
			h_rcoutstatus = 0L;
		}
		h_rcoutstatusflip = 0L;
		h_rcoutstatusflop = 0L;
		h_rcoutstatusflap = 0L;
	}

	void Create_Configuration(
			int arg_Kernel_Type,
			char* arg_device_name,
			int arg_device, 		// not needed?
			int arg_num_timesteps,		// number of timesteps to be performed on this GPU
			int arg_total_num_timesteps,	// total number of timesteps for entire pipeline
			int arg_total_num_stripes,	// total number of stripes in propagation volume
			int arg_OTflag,			// temporal order, either 2 or 4
			int arg_bY,			// output block size Y
			int arg_bZ,			// output block size Z
			int arg_halo_ny,		// number of halo cells Y, one side
			int arg_exclude_halo_ylo,	// flag: true means halo should be excluded on low Y side for each output
			int arg_exclude_halo_yhi,
			int arg_y0			// y0 of pipeline
			)
	{
		h2d_byte_count = 0;
		d2h_byte_count = 0;

		Kernel_Type = arg_Kernel_Type;

		d_temp1 = 0L;
		d_temp2 = 0L;
		d_temp3 = 0L;

		d_spg_x = 0L;
		d_spg_y = 0L;
		d_spg_z = 0L;
	
		h_src_wavelet = 0L;
		src_len = 0;
		xsrc = -1;
		ysrc = -1;
		zsrc = -1;
		zsrcghost = -1;

		_recghost_Flag = 0;
		_z_freesurface = 0;

		d_rcbuf = 0L;
		d_rcx = 0L;
		d_rcy = 0L;
		d_rcz = 0L;
		d_rcout = 0L;

		h_rcbuf = 0L;
		h_rcx = 0L;
		h_rcy = 0L;
		h_rcz = 0L;
		h_idx = 0L;
		h_rcout = 0L;
		h_rcoutmaxbinsize = 0;

		h_recnum = 0L;
		h_recnx = 0;
		h_totrecnx = 0;

		h_rcoutstatus = 0L;
		h_rcoutstatusflip = 0L;
		h_rcoutstatusflop = 0L;

		// zero out propagation parameters
		// use defaults that will not generate NaN values.

		dt = 1e-3;

		const float deg2rad = 3.1415926536f / 180.0f;

		dip_min = 45.0f * deg2rad;
		dip_scaler = 0.0f;

		azm_min = 135.0f * deg2rad;
		azm_scaler = 0.0f;

		density_min = 1.0f;
		density_scaler = 0.0f;

		Q_min = 1.0f;
		Q_scaler = 0.0f;

		Vel_min = 2500.0f;
		Vel_scaler = 0.0f;

		C44C33_min = 0.01f;
		C44C33_scaler = 0.0f;

		Del_min = 0.021043f;
		Del_scaler = 0.0f;

		Eps_min = 0.094139f;
		Eps_scaler = 0.0f;

		// copy all arguments.
		device_name = strdup(arg_device_name);
		device = arg_device;
		canAccessPeer = 0;
		cc_major = 0;
		cc_minor = 0;
		num_timesteps = arg_num_timesteps;
		total_num_timesteps = arg_total_num_timesteps;
		total_num_stripes = arg_total_num_stripes;
		OTflag = arg_OTflag;
		bY = arg_bY;
		bZ = arg_bZ;
		halo_ny = arg_halo_ny;
		exclude_halo_ylo = arg_exclude_halo_ylo;
		exclude_halo_yhi = arg_exclude_halo_yhi;

		// determine configuration
		input_stream = 0;
		compute_stream = 0;
		output_stream = 0;
		num_layers = 0;
		layers = 0L;  // default, in case of failure
		if (OTflag == 2)
		{
			num_layers = 3 + num_timesteps;
			layers = new TimeStepLayer*[num_layers];
			
			layers[0] = new TimeStepLayer(device,"Model",0,0,0,8,arg_total_num_timesteps,arg_total_num_stripes,arg_y0);
                        layers[1] = new TimeStepLayer(device," pq[-1]",0,1,-1,8,arg_total_num_timesteps,arg_total_num_stripes,arg_y0);
                        layers[2] = new TimeStepLayer(device," pq[ 0]",0,1,0,8,arg_total_num_timesteps,arg_total_num_stripes,arg_y0);

			char buf[1024];
			int compute_dim_y = bY, input_off_y = 0;
			for (int iStep = num_timesteps;  iStep >= 1;  --iStep)
			{
				int layer_num = 2 + iStep;

				sprintf(buf, " pq[%2d]", iStep);
				layers[layer_num] = new TimeStepLayer(device,buf,1,1,iStep,8,arg_total_num_timesteps,arg_total_num_stripes,arg_y0);

				layers[layer_num]->input_off_x = 0;  // for now, adjusted later.
				layers[layer_num]->input_off_y = input_off_y - halo_ny;
				layers[layer_num]->input_off_z = 0;

				layers[layer_num]->compute_off_x = (iStep+2) * -16;
				layers[layer_num]->compute_off_y = layers[layer_num]->input_off_y + halo_ny;
				layers[layer_num]->compute_off_z = 0;

				layers[layer_num]->compute_dim_y = compute_dim_y;
				layers[layer_num]->compute_dim_z = bZ;

				layers[layer_num]->compute_out_y = exclude_halo_ylo ? (iStep == num_timesteps ? 0 : halo_ny) : 0;
				layers[layer_num]->compute_out_z = 0;

				layers[layer_num]->curr_input_stripe = -4 - iStep;

				if (iStep == num_timesteps)
				{
					layers[layer_num]->num_stripes = 4;
					layers[layer_num]->output_dim_y = compute_dim_y;
					layers[layer_num]->output_dim_z = bZ;
					layers[layer_num]->is_output_layer = 1;
				}
				else if (iStep == num_timesteps-1)
				{
					layers[layer_num]->num_stripes = 5;
					layers[layer_num]->output_dim_y = compute_dim_y + (exclude_halo_ylo ? halo_ny : 0) + (exclude_halo_yhi ? halo_ny : 0);
					layers[layer_num]->output_dim_z = bZ;
				}
				else
				{
					layers[layer_num]->num_stripes = 4;
					layers[layer_num]->output_dim_y = compute_dim_y + (exclude_halo_ylo ? halo_ny : 0) + (exclude_halo_yhi ? halo_ny : 0);
					layers[layer_num]->output_dim_z = bZ;
				}

				if (!exclude_halo_ylo)
				{
					compute_dim_y += halo_ny;
					input_off_y -= halo_ny;
				}
				if (!exclude_halo_yhi)
				{
					compute_dim_y += halo_ny;
				}
			}

			layers[2]->num_stripes = 6;
			layers[2]->curr_input_stripe = -2;
			layers[2]->input_off_y = layers[3]->input_off_y;
			layers[2]->input_off_z = layers[3]->input_off_z;
			layers[2]->output_dim_y = layers[3]->compute_dim_y + 2 * halo_ny;
			layers[2]->output_dim_z = layers[3]->compute_dim_z;

			layers[1]->num_stripes = 5;
			layers[1]->curr_input_stripe = -2;
			layers[1]->input_off_y = layers[3]->compute_off_y;
			layers[1]->input_off_z = layers[3]->compute_off_z;
			layers[1]->output_dim_y = layers[3]->compute_dim_y;
			layers[1]->output_dim_z = layers[3]->compute_dim_z;

			layers[0]->num_stripes = 6 + num_timesteps;
			layers[0]->curr_input_stripe = -2;
			layers[0]->input_off_y = layers[3]->compute_off_y - 5;
			layers[0]->input_off_z = layers[3]->compute_off_z;
			layers[0]->output_dim_y = layers[3]->compute_dim_y + halo_ny;
			layers[0]->output_dim_z = layers[3]->compute_dim_z;
		}
		else if (OTflag == 4)
		{
			num_layers = 3 + 2 * num_timesteps;
			layers = new TimeStepLayer*[num_layers];

			layers[0] = new TimeStepLayer(device,"Model",0,0,0,8,arg_total_num_timesteps,arg_total_num_stripes,arg_y0);
			layers[1] = new TimeStepLayer(device," pq[-1]",0,-1,1,8,arg_total_num_timesteps,arg_total_num_stripes,arg_y0);
			layers[2] = new TimeStepLayer(device," pq[ 0]",0,0,1,8,arg_total_num_timesteps,arg_total_num_stripes,arg_y0);

			char buf[1024];
			int compute_dim_y = bY, input_off_y = 0;
			for (int iStep = num_timesteps;  iStep >= 1;  --iStep)
			{
				int layer_num = 2 + iStep * 2;

				sprintf(buf, "Apq[%2d]", iStep);
				layers[layer_num-1] = new TimeStepLayer(device,buf,1,0,iStep,8,arg_total_num_timesteps,arg_total_num_stripes,arg_y0);

				sprintf(buf, " pq[%2d]", iStep);
				layers[layer_num  ] = new TimeStepLayer(device,buf,1,1,iStep,8,arg_total_num_timesteps,arg_total_num_stripes,arg_y0);

				// process pq and Apq in reverse order

				// pq
				layers[layer_num]->input_off_y = input_off_y - halo_ny;
				layers[layer_num]->input_off_z = 0;

				layers[layer_num]->compute_off_x = (iStep+1) * -32;
				layers[layer_num]->compute_off_y = layers[layer_num]->input_off_y + halo_ny;
				layers[layer_num]->compute_off_z = 0;

				layers[layer_num]->compute_dim_y = compute_dim_y;
				layers[layer_num]->compute_dim_z = bZ;

				layers[layer_num]->compute_out_y = exclude_halo_ylo ? (iStep == num_timesteps ? 0 : halo_ny) : 0;
				layers[layer_num]->compute_out_z = 0;

				layers[layer_num]->curr_input_stripe = -3 - 2 * iStep;

				if (iStep == num_timesteps)
				{
					layers[layer_num]->num_stripes = 4;
					layers[layer_num]->output_dim_y = compute_dim_y;
					layers[layer_num]->output_dim_z = bZ;
				}
				else if (iStep == num_timesteps-1)
				{
					layers[layer_num]->num_stripes = 6;
					layers[layer_num]->output_dim_y = compute_dim_y + (exclude_halo_ylo ? halo_ny : 0) + (exclude_halo_yhi ? halo_ny : 0);
					layers[layer_num]->output_dim_z = bZ;
				}
				else
				{
					layers[layer_num]->num_stripes = 6;
					layers[layer_num]->output_dim_y = compute_dim_y + (exclude_halo_ylo ? halo_ny : 0) + (exclude_halo_yhi ? halo_ny : 0);
					layers[layer_num]->output_dim_z = bZ;
				}

				if (!exclude_halo_ylo)
				{
					compute_dim_y += halo_ny;
					input_off_y -= halo_ny;
				}
				if (!exclude_halo_yhi)
				{
					compute_dim_y += halo_ny;
				}

				// Apq
				layers[layer_num-1]->input_off_y = input_off_y - halo_ny;
				layers[layer_num-1]->input_off_z = 0;

				layers[layer_num-1]->compute_off_x = (iStep+1) * -32 + 16;
				layers[layer_num-1]->compute_off_y = layers[layer_num-1]->input_off_y + halo_ny;
				layers[layer_num-1]->compute_off_z = 0;

				layers[layer_num-1]->compute_dim_y = compute_dim_y;
				layers[layer_num-1]->compute_dim_z = bZ;

				layers[layer_num-1]->compute_out_y = exclude_halo_ylo ? halo_ny : 0;
				layers[layer_num-1]->compute_out_z = 0;

				layers[layer_num-1]->num_stripes = 4;
				layers[layer_num-1]->output_dim_y = compute_dim_y + (exclude_halo_ylo ? halo_ny : 0) + (exclude_halo_yhi ? halo_ny : 0);
				layers[layer_num-1]->output_dim_z = bZ;

				layers[layer_num-1]->curr_input_stripe = -4 - 2 * iStep;

				if (!exclude_halo_ylo)
				{
					compute_dim_y += halo_ny;
					input_off_y -= halo_ny;
				}
				if (!exclude_halo_yhi)
				{
					compute_dim_y += halo_ny;
				}
			}

			layers[2]->num_stripes = 6;
			layers[2]->curr_input_stripe = -2;
			layers[2]->input_off_y = layers[3]->input_off_y;
			layers[2]->input_off_z = layers[3]->input_off_z;
			layers[2]->output_dim_y = layers[3]->compute_dim_y + 2 * halo_ny;
			layers[2]->output_dim_z = layers[3]->compute_dim_z;

			layers[1]->num_stripes = 8;
			layers[1]->curr_input_stripe = -2;
			layers[1]->input_off_y = layers[4]->compute_off_y;
			layers[1]->input_off_z = layers[4]->compute_off_z;
			layers[1]->output_dim_y = layers[4]->compute_dim_y;
			layers[1]->output_dim_z = layers[4]->compute_dim_z;

			layers[0]->num_stripes = 6 + 2 * num_timesteps;
			layers[0]->curr_input_stripe = -2;
			layers[0]->input_off_y = layers[3]->compute_off_y - 5;
			layers[0]->input_off_z = layers[3]->compute_off_z;
			layers[0]->output_dim_y = layers[3]->compute_dim_y + halo_ny;
			layers[0]->output_dim_z = layers[3]->compute_dim_z;
		}

		if (Succeeding_Propagator != 0L)
		{
			if (!exclude_halo_ylo)
			{
				int succ_input_off_y, succ_compute_off_x;
				Succeeding_Propagator->Get_Offsets(succ_compute_off_x,succ_input_off_y); 
				for (int i =0;  i < num_layers;  ++i)
				{
					if (layers[i]->is_compute_layer)
					{
						layers[i]->compute_off_y += succ_input_off_y; 
					}
					else
					{
						layers[i]->input_off_y += succ_input_off_y;
					}
				}
			}

			Succeeding_Propagator->Adjust_X_Offsets(this->Get_Output_X_Offset(),layers[num_layers-1]->timestep);
		}
	}

	void Adjust_X_Offsets(int adjust_off_x, int adjust_timesteps)
	{
		int adjust_off_stripe = adjust_off_x / 16;
		for (int i =0;  i < num_layers;  ++i)
		{
			if (layers[i]->is_compute_layer) 
			{
				layers[i]->compute_off_x += adjust_off_x;
			}
			else /* input layer */
			{
				layers[i]->input_off_x += adjust_off_x;
			}
			layers[i]->curr_input_stripe += adjust_off_stripe;
			layers[i]->timestep += adjust_timesteps;
		}
		if (Succeeding_Propagator != 0L) Succeeding_Propagator->Adjust_X_Offsets(adjust_off_x,adjust_timesteps);
	}

	void Initialize_Earth_Model_Decompression(
		double arg_dt,
		float arg_dip_scaler,
		float arg_dip_min,
		float arg_azm_scaler,
		float arg_azm_min,
		float arg_density_scaler,
		float arg_density_min,
		float arg_Q_scaler,
		float arg_Q_min,
		float arg_Vel_scaler,
		float arg_Vel_min,
		float arg_C44C33_scaler,
		float arg_C44C33_min,
		float arg_Del_scaler,
		float arg_Del_min,
		float arg_Eps_scaler,
		float arg_Eps_min
		)
	{
		dt = arg_dt;
		dip_scaler = arg_dip_scaler;
		dip_min = arg_dip_min;
		azm_scaler = arg_azm_scaler;
		azm_min = arg_azm_min;
		density_scaler = arg_density_scaler;
		density_min = arg_density_min;
		Q_scaler = arg_Q_scaler;
		Q_min = arg_Q_min;
		Vel_scaler = arg_Vel_scaler;
		Vel_min = arg_Vel_min;
		C44C33_scaler = arg_C44C33_scaler;
		C44C33_min = arg_C44C33_min;
		Del_scaler = arg_Del_scaler;
		Del_min = arg_Del_min;
		Eps_scaler = arg_Eps_scaler;
		Eps_min = arg_Eps_min;
	}

	friend class DevicePropagator;
	friend class SinglePipelinePropagator;

	// arguments
	SingleDevicePropagator* Succeeding_Propagator;

	int Kernel_Type;
	char* device_name;
	int device;
	int canAccessPeer;
	int cc_major;
	int cc_minor;
	int num_timesteps;
	int total_num_timesteps;
	int total_num_stripes;
	int OTflag;
	int bY;
	int bZ;
	int halo_ny;
	int exclude_halo_ylo;
	int exclude_halo_yhi;

	// configuration
	int num_layers;
	TimeStepLayer** layers;

	// device resources
	cudaStream_t input_stream;
	cudaStream_t compute_stream;
	cudaStream_t output_stream;

	// monitoring
	unsigned long h2d_byte_count;
	unsigned long d2h_byte_count;

	//
	// GPU specific code. Eventually move to subclass.
	//

	double dt;

	float dip_scaler;
	float dip_min;
	float azm_scaler;
	float azm_min;
	float density_scaler;
	float density_min;
	float Q_scaler;
	float Q_min;
	float Vel_scaler;
	float Vel_min;
	float C44C33_scaler;
	float C44C33_min;
	float Del_scaler;
	float Del_min;
	float Eps_scaler;
	float Eps_min;

	float* d_temp1;
	float* d_temp2;
	float* d_temp3;

	//
	// d_spg_x
	//
	// In some situations, a block can contain two stripes that are on opposite ends of the volume.
	// One is on the extreme right (x=[dimx-16,dimx-1]), the other on the extreme left (x=[0,15]).
	// To handle this situation, d_spg_x has dimx+16 values, the last 16 are for x=[0,15] in the above case.
	//	
	float* d_spg_x;
	float* d_spg_y;
	float* d_spg_z;

	int src_len;
	float* h_src_wavelet;
	int xsrc;
	int ysrc;
	int zsrc;
	int zsrcghost;

	int* d_rcbuf;
	int** d_rcx;
	int** d_rcy;
	int** d_rcz;
	float* d_rcout;
	float* d_rcoutflip;
	float* d_rcoutflop;

	int _recghost_Flag;
	int _z_freesurface;
	int* h_rcbuf;
	int** h_rcx;
	int** h_rcy;
	int** h_rcz;
	int** h_idx;
	float* h_rcout;
	float* h_rcoutflop;
	float* h_rcoutflap;
	unsigned long h_rcoutmaxbinsize;

	int* h_recnum;
	int h_recnx;
	int h_totrecnx;

	int* h_rcoutstatus;
	int* h_rcoutstatusflip;
	int* h_rcoutstatusflop;
	int* h_rcoutstatusflap;

	void Compute_Stencils(int OTflag, float dx, float dy, float dz)
	{
		float a1, a2, a3, a4, a5; 
		//float b1, b2, b3, b4;
		//float c1, c2, c3;
		//float d1, d2;
		float e1;

		if(OTflag==4)
		{
			a1 =  1.250394163714f;  a2 = -0.119656543874f;  a3 = 0.031206223579f;
			a4 = -0.009128136972f;  a5 =  0.001882183398f;
		}
		else
		{
			a1 =  1.248489029341f;  a2 = -0.120133754290f;  a3 = 0.031688119039f;
			a4 = -0.008048796917f;  a5 =  0.001090357653f;
		}
		//b1 =  1.231650129521f;  b2 = -0.103861125624f;  b3 = 0.020166542235f;
		//b4 = -0.002985637689f;
		//c1 =  1.199634495725f;  c2 = -0.080370339530f;  c3 = 0.008295304573f;
		//d1 =  1.134389630713f;  d2 = -0.044796543571f;
		e1 =  0.5f;

		dx_A1 = a1/dx;  dx_A2 = a2/dx;  dx_A3 = a3/dx;  dx_A4 = a4/dx;  dx_A5 = a5/dx;

		dy_A1 = a1/dy;  dy_A2 = a2/dy;  dy_A3 = a3/dy;  dy_A4 = a4/dy;  dy_A5 = a5/dy;

		dz_A1 = a1/dz;  dz_A2 = a2/dz;  dz_A3 = a3/dz;  dz_A4 = a4/dz;  dz_A5 = a5/dz;
		dz_E1 = e1/dz;
	}

	float dx_A1, dx_A2, dx_A3, dx_A4, dx_A5;
	float dy_A1, dy_A2, dy_A3, dy_A4, dy_A5;
	float dz_A1, dz_A2, dz_A3, dz_A4, dz_A5;
	float dz_E1;
};

class SinglePipelinePropagator
{
public:
	SinglePipelinePropagator(
		int arg_Kernel_Type,
		int arg_dimx, 
		int arg_y0,
		int arg_y1, 
		int arg_dimy,
		int arg_dimz, 
		int arg_first_device,
		int arg_device_count,
		int arg_exclude_halo_ylo,
		int arg_exclude_halo_yhi,
		int arg_num_timesteps_per_device,
		int arg_halo_ny,
		int arg_OTflag,
		float dx,
		float dy,
		float dz
		)
	{
		// copy arguments
		dimx = arg_dimx;
		y0 = arg_y0;
		y1 = arg_y1;
		dimy = arg_dimy;
		dimz = arg_dimz;
		first_device = arg_first_device;
		device_count = arg_device_count;
		exclude_halo_ylo = arg_exclude_halo_ylo;
		exclude_halo_yhi = arg_exclude_halo_yhi;
		num_timesteps_per_device = arg_num_timesteps_per_device;
		halo_ny = arg_halo_ny;
		OTflag = arg_OTflag;

		int total_num_stripes = dimx / 16;
		int total_num_timesteps = num_timesteps_per_device * device_count;

		// create device propagators
		int bY = y1 - y0 + 1;
		props = new SingleDevicePropagator*[device_count];
		for (int iDev = device_count - 1;  iDev >= 0;  --iDev)
		{
			int device_no = first_device + iDev;
			char device_name[1024];
			sprintf(device_name, "device %d", device_no);
			if (iDev == device_count - 1)
			{
				props[iDev] = new SingleDevicePropagator(
						arg_Kernel_Type,
						device_name,device_no,
						num_timesteps_per_device,total_num_timesteps,total_num_stripes,OTflag,
						bY,dimz,halo_ny,exclude_halo_ylo,exclude_halo_yhi,y0
						);
			}
			else
			{
				props[iDev] = new SingleDevicePropagator(arg_Kernel_Type,props[iDev+1],device_name,device_no,num_timesteps_per_device,total_num_timesteps,total_num_stripes,y0);
			}
			props[iDev]->Compute_Stencils(OTflag,dx,dy,dz);
		}
	}

	~SinglePipelinePropagator()
	{
		//printf("~SinglePipelinePropagator\n");
		fflush(stdout);
		if (props != 0L)
		{	
			for (int i = 0;  i < device_count;  ++i)
			{
				delete props[i];
				props[i] = 0L;
			}
			delete [] props;
			props = 0L;
		}
	}

	void Put_Source_Wavelet(
		float* src,
		int len,
		int xsrc,
		int ysrc,
		int zsrc,
		int zsrcghost
		)
	{
		if (props != 0L)
                {
                        for (int i = 0;  i < device_count;  ++i)
                        {
				props[i]->Put_Source_Wavelet(src,len,xsrc,ysrc,zsrc,zsrcghost);
			}
		}
	}

	int Put_Receiver_Locations(
		int* arg_rcx,
		int* arg_rcy,
		int* arg_rcz,
		int arg_recnum,
		int dimx,
		int recghost_Flag,
		int z_freesurface
		)
	{
		int retval = -1;
		if (props != 0L)
                {
                        for (int i = 0;  i < device_count;  ++i)
                        {
				retval = props[i]->Put_Receiver_Locations(arg_rcx,arg_rcy,arg_rcz,arg_recnum,dimx,y0,y1,recghost_Flag,z_freesurface) ? retval : 0;
                        }
                }
		return retval;
	}

	void Put_Sponges(
		float* h_spg_x,
		int dimx,
		float* h_spg_y,
		int dimy,
		float* h_spg_z,
		int dimz
		)
	{
		if (props != 0L)
		{
			for (int i = 0;  i < device_count;  ++i)
			{
				props[i]->Put_Sponges(h_spg_x,dimx,h_spg_y,dimy,h_spg_z,dimz);
			}
		}
	}

	void Initialize_Earth_Model_Decompression(
		double dt,
		float dip_scaler,
		float dip_min,
		float azm_scaler,
		float azm_min,
		float density_scaler,
		float density_min,
		float Q_scaler,
		float Q_min,
		float Vel_scaler,
		float Vel_min,
		float C44C33_scaler,
		float C44C33_min,
		float Del_scaler,
		float Del_min,
		float Eps_scaler,
		float Eps_min
		)
	{
                if (props != 0L)
                {
                        for (int i = 0;  i < device_count;  ++i)
                        {
				props[i]->Initialize_Earth_Model_Decompression(dt,dip_scaler,dip_min,azm_scaler,azm_min,density_scaler,density_min,Q_scaler,Q_min,Vel_scaler,Vel_min,C44C33_scaler,C44C33_min,Del_scaler,Del_min,Eps_scaler,Eps_min);
			}
		}		
	}
	
	void Move_Sliding_Window(int logLevel)
	{
		for (int iDev = 0;  iDev < device_count;  ++iDev)
                {
                        props[iDev]->Move_Sliding_Window(logLevel);
		}
	}

	//
	// All inputs point to the current input / output stripe.
	//
	/*
	void Start_Data_Flow(
		void* h_em_inp, 
		void* h_pq_inp, 
		void* h_prev_pq_inp, 
		void* h_pq_out, 
		void* h_prev_pq_out
		)
	{
		for (int iDev = 0;  iDev < device_count;  ++iDev)
		{
			if (iDev == 0)
			{
				// first device, so input from host
				props[iDev]->CopyFromHostToDevice(0,-1,h_em_inp,y0,y1,dimy);
				props[iDev]->CopyFromHostToDevice(1,-1,h_pq_inp,y0,y1,dimy);
				props[iDev]->CopyFromHostToDevice(2,-1,h_prev_pq_inp,y0,y1,dimy);
			}
			else
			{
				// peer-2-peer transfer
				int src_num_layers = props[iDev-1]->Get_Number_of_Layers();
				int src_m1_layer_id = OTflag == 2 ? src_num_layers-2 : src_num_layers - 3;
				int src_p0_layer_id = src_num_layers - 1;
				props[iDev]->CopyFromDeviceToDevice(props[iDev-1],0,0,0,-1,y0,y1,dimy);   // earth model
				props[iDev]->CopyFromDeviceToDevice(props[iDev-1],src_m1_layer_id,0,1,-1,y0,y1,dimy);
				props[iDev]->CopyFromDeviceToDevice(props[iDev-1],src_p0_layer_id,0,2,-1,y0,y1,dimy);
			}
			if (iDev == device_count - 1)
			{
				// last device, so also output final result to host,
				// but only after pipeline fills up.
				if (h_pq_out != 0L && h_prev_pq_out != 0L)
				{
					int src_num_layers = props[iDev]->Get_Number_of_Layers();
					int src_p0_layer_id = OTflag == 2 ? src_num_layers - 2 : src_num_layers - 3;
					int src_p1_layer_id = src_num_layers - 1;
					props[iDev]->CopyFromDeviceToHost(src_p0_layer_id,0,h_prev_pq_out,y0,y1);
					props[iDev]->CopyFromDeviceToHost(src_p1_layer_id,0,h_pq_out,y0,y1);
				}
			}
		}
	}

	void Start_Input_Stream()
	{
		for (int iDev = 0;  iDev < device_count;  ++iDev)
                {
			int is_first = iDev == 0 ? 1 : 0;
			int is_last = iDev == device_count-1 ? 1 : 0;
			props[iDev]->Start_Input_Stream(is_first, is_last);
		}
	}

	void LaunchSimpleCopyKernel()
	{
		for (int iDev = 0;  iDev < device_count;  ++iDev)
		{
			props[iDev]->LaunchSimpleCopyKernel();
		}
	}

	void Start_Output_Stream()
	{
		for (int iDev = 0;  iDev < device_count;  ++iDev)
                {
			int is_first = iDev == 0 ? 1 : 0;
			int is_last = iDev == device_count-1 ? 1 : 0;
			props[iDev]->Start_Output_Stream(is_first, is_last);
		}
	}

	void WaitForInputStream()
	{
		for (int iDev = 0;  iDev < device_count;  ++iDev)
		{
			props[iDev]->WaitForInputStream();
		}
	}

	void WaitForComputeStream()
	{
		for (int iDev = 0;  iDev < device_count;  ++iDev)
		{
			props[iDev]->WaitForComputeStream();
		}
	}

	void WaitForOutputStream()
	{
		for (int iDev = 0;  iDev < device_count;  ++iDev)
		{
			props[iDev]->WaitForOutputStream();
		}
	}
	*/

	void Compute_Y_Extent(
		int device_id,
		int layer_id,
		int& layer_y0,
		int& layer_y1
		)
	{
		SingleDevicePropagator* prop = Get_Device(device_id);
		if (prop != 0L)
		{
			prop->Compute_Y_Extent(y0,y1,layer_id,layer_y0,layer_y1);
		}
		else
		{
			layer_y0 = 0;
			layer_y1 = 0;
		}
	}

	SingleDevicePropagator* Get_Device(int device_id)
	{
		if (device_id >= -device_count && device_id < device_count)
		{
			int device = (device_id + device_count) % device_count;
			return props[device];
		}
		return 0L;
	}

	//
	// Returns the X delay, in number of stripes, for the entire pipeline.
	// This is a function of the temporal order and the total number of timesteps.
	//
	int Get_X_Delay()
	{
		return props[device_count-1]->Get_X_Delay();
	}

	int Get_Total_Number_of_Timesteps()
	{
		int acc_ts = 0;
		for (int iDev = 0;  iDev < device_count;  ++iDev)
                {
			acc_ts += props[iDev]->Get_Number_of_Timesteps();
		}
		return acc_ts;
	}

	int Get_Y0()
	{
		return y0;
	}
	
	int Get_Y1()
	{
		return y1;
	}

	int Allocate_Device_Resources(int dimx, int dimy, int dimz)
	{
		int retval = 1;
		for (int iDev = 0;  retval != 0 && iDev < device_count;  ++iDev)
		{
			if (props[iDev]->Allocate_Device_Resources(dimx,dimy,dimz) == 0)
			{	
				retval = 0;
			}
		}
		return retval;
	}

	double Compute_Max_Memory_Required_GB()
	{
		double max_mem_GB = 0.0;
		for (int iDev = 0;  iDev < device_count;  ++iDev)
		{
			double mem_GB = props[iDev]->Compute_Memory_Required_GB();
			if (mem_GB > max_mem_GB) max_mem_GB = mem_GB;
		}
		return max_mem_GB;
	}

	double Compute_Max_Bytes_Per_Cell_Per_Timestep()
	{
		double retval = 0.0;
		for (int iDev = 0;  iDev < device_count;  ++iDev)
		{
			double bytes_per_cell_per_timestep = props[iDev]->Compute_Bytes_Per_Cell_Per_Timestep();
			if (bytes_per_cell_per_timestep > retval) retval = bytes_per_cell_per_timestep;
		}
		return retval;
	}

	double Estimate_Computational_Cost()
	{
		double cost = 0.0;
		for (int iDev = 0;  iDev < 1;/*device_count;*/  ++iDev)
                {
			cost += props[iDev]->Estimate_Computational_Cost();
		}
		return cost;
	}

	unsigned long Get_H2D_Byte_Count()
	{
		unsigned long acc_count = 0;
		for (int iDev = 0;  iDev < device_count;  ++iDev)
		{
			acc_count += props[iDev]->Get_H2D_Byte_Count();
		}
		return acc_count;
	}

	unsigned long Get_D2H_Byte_Count()
	{
		unsigned long acc_count = 0;
		for (int iDev = 0;  iDev < device_count;  ++iDev)
		{
			acc_count += props[iDev]->Get_D2H_Byte_Count();
		}
		return acc_count;
	}

	void Print()
	{
		for (int iDev = 0;  iDev < device_count;  ++iDev)
		{
			props[iDev]->Print();
			printf("\n");
		}
	}

	//
	// Get device propagator.
	// iDev is in range [0,Get_Device_Count()-1]
	//
	SingleDevicePropagator* Get_Device_Propagator(int iDev, int& is_first, int& is_last, SingleDevicePropagator*& preceeding_propagator, int& pipeline_y0, int& pipeline_y1)
	{
		pipeline_y0 = y0;
		pipeline_y1 = y1;
		if (iDev >= 0 && iDev < device_count)
		{
			if (iDev == 0)
			{
				is_first = 1;
				preceeding_propagator = 0L;
			}
			else
			{
				is_first = 0;
				preceeding_propagator = props[iDev-1];
			}
			is_last = (iDev == device_count - 1) ? 1 : 0;
			return props[iDev];
		}
		is_first = 0;
		is_last = 0;
		preceeding_propagator = 0L;
		return 0L;
	}

	int Get_Device_Count()
	{
		return device_count;
	}

private:
	int dimx;
	int y0;
	int y1; 
	int dimy;
	int dimz;
	int first_device;
	int device_count;
	int exclude_halo_ylo;
	int exclude_halo_yhi;
	int num_timesteps_per_device;
	int halo_ny;
	int OTflag;

	friend class DevicePropagator;

	SingleDevicePropagator** props;
};

class DevicePropagator
{
public:
	DevicePropagator(
		int arg_Kernel_Type,
		int arg_dimx, 
		int arg_dimy, 
		int arg_dimz, 
		int arg_first_device,
		int arg_device_count,
		int arg_num_y_pipelines,
		int arg_num_timesteps_per_device,
		int arg_halo_ny,
		int arg_OTflag,
		float arg_dx,
		float arg_dy,
		float arg_dz,
		float arg_dt
		)
	{
		is_valid = 0;
		dimx = arg_dimx;
		dimy = arg_dimy;
		dimz = arg_dimz;
		first_device = arg_first_device;
		device_count = arg_device_count;
		num_y_pipelines = arg_num_y_pipelines;
		num_timesteps_per_device = arg_num_timesteps_per_device;
		halo_ny = arg_halo_ny;
		OTflag = arg_OTflag;
		dx = arg_dx;
		dy = arg_dy;
		dz = arg_dz;
		dt = arg_dt;

		single_device_multi_track = (device_count == 1 && num_y_pipelines > 1) ? 1 : 0;

		if (!single_device_multi_track && (device_count < 1 || device_count < num_y_pipelines || (device_count % num_y_pipelines) != 0))
		{
			printf("device_count must be a multiple of num_y_pipelines and cannot be zero!\n");
			pipelines = 0L;
			return;
		}

		h_pinned = 0;
		h_pq = 0L;
		h_prev_pq = 0L;
		h_em_deallocate = 0;
		h_em = 0L;
		h_src_wavelet = 0L;

		Kernel_Type = arg_Kernel_Type;

		rcx = 0L;
		rcy = 0L;
		rcz = 0L;
		nsamp = 0;
		recnum = 0;
		muxed_traces = 0L;
		rec_out_steps = 0L;
		time_delta = 0L;
		receiver_num_timesteps = 0;
		mux_imgstep = 0L;
		mux_time_delta = 0L;

		// create GPU configuration
		is_valid = 1;

		int num_y_per_pipeline = (dimy + num_y_pipelines - 1) / num_y_pipelines;
		pipe_dimy = new int[num_y_pipelines];
		for (int iY = 0;  iY < num_y_pipelines;  ++iY)
		{
			pipe_dimy[iY] = num_y_per_pipeline;
		}
		pipelines = 0L;
		shared = 0L;

		Create_Configuration(arg_Kernel_Type, pipe_dimy, num_y_pipelines);
		Optimize();

		size_of_one_stripe_in_floats = (unsigned long)32 * (unsigned long)dimy * (unsigned long)dimz;
		total_num_stripes = dimx / 16;

		curr_input_slice = -2;

		curr_input_stripe_number = total_num_stripes - 2;
		curr_input_timestep = -1;
	}

	~DevicePropagator()
	{
		//printf("~DevicePropagator\n");
		fflush(stdout);
		Free_Receiver_Locations();
		if (shared != 0L) delete [] shared;
		if (pipelines != 0L)
		{
			for (int iY = 0;  iY < num_y_pipelines;  ++iY)
			{
				delete pipelines[iY];
				pipelines[iY] = 0L;
			}
			delete [] pipelines;
			pipelines = 0L;
		}
		if (h_pq != 0L) cudaFreeHost((void*)h_pq);
		//printf("Freed h_pq\n");
		fflush(stdout);

		if (h_prev_pq != 0L) cudaFreeHost((void*)h_prev_pq);
		//printf("Freed h_prev_pq\n");
		fflush(stdout);

		if (h_em_deallocate != 0 && h_em != 0L) cudaFreeHost((void*)h_em);
		//printf("Freed h_em\n");
		fflush(stdout);

		if (h_src_wavelet != 0L) free((void*)h_src_wavelet);
		//printf("Freed h_src_wavelet\n");
		fflush(stdout);
	}

	void Optimize()
	{
		double* cost = new double[num_y_pipelines];
		for (int iter = 0;  iter < 1000;  ++iter)
		{
			Estimate_Computational_Cost(cost);

			// add 1 cell to cheapest pipeline
			int cheapest_idx = 0;
			double cheapest_cost = cost[0];
			for (int iY = 1;  iY < num_y_pipelines;  ++iY)
			{
				if (cost[iY] < cheapest_cost)
				{
					cheapest_idx = iY;
					cheapest_cost = cost[iY];
				}
			}
			pipe_dimy[cheapest_idx] += 1;
			
			// subtract 1 cell from most expensive pipeline
			int expensive_idx = 0;
			double expensive_cost = cost[0];
			for (int iY = 1;  iY < num_y_pipelines;  ++iY)
                        {
                                if (cost[iY] > expensive_cost)
                                {
                                        expensive_idx = iY;
                                        expensive_cost = cost[iY];
                                }
                        }
			pipe_dimy[expensive_idx] -= 1;

			Create_Configuration(Kernel_Type, pipe_dimy, num_y_pipelines);
		}
	}

	// estimate computational cost of each pipeline.
	// values are normalized so that "cheapest" pipeline has cost 1.0
	void Estimate_Computational_Cost(double* cost)
	{
		for (int iY = 0;  iY < num_y_pipelines;  ++iY)
		{
			cost[iY] = pipelines[iY]->Estimate_Computational_Cost();
		}
		
		double min = cost[0];
		for (int iY = 1;  iY < num_y_pipelines;  ++iY)
		{
			if (cost[iY] < min) min = cost[iY];
		}

		for (int iY = 0;  iY < num_y_pipelines;  ++iY)
                {
			cost[iY] /= min;
		}
	}

	int Allocate_Host_Memory(int pinned, int functionality_test, int* em)
	{
		int nx = dimx / 16;

		unsigned long pqNCells = (unsigned long)nx * (unsigned long)16 * (unsigned long)dimy * (unsigned long)dimz;

		cudaError_t err1 = cudaHostAlloc((void**)&h_pq, (size_t)(pqNCells * (unsigned long)8), cudaHostAllocDefault);
		if (err1 == cudaSuccess)
		{
			cudaError_t err2 = cudaHostAlloc((void**)&h_prev_pq, (size_t)(pqNCells * (unsigned long)8), cudaHostAllocDefault);
			if (err2 == cudaSuccess)
			{
				cudaError_t err3;
				if (em != 0L)
				{
					h_em_deallocate = 0;
					h_em = em;
					err3 = cudaSuccess;
				}
				else
				{
					h_em_deallocate = 1;
					err3 = cudaHostAlloc((void**)&h_em, (size_t)(pqNCells * (unsigned long)8), cudaHostAllocDefault);
				}
				if (err3 == cudaSuccess)
				{
					if (functionality_test == 1)
					{
						// fill with pattern that can be verified
						for (unsigned long i = 0;  i < pqNCells*2;  ++i)
						{
							((float*)h_pq)[i] = (float)i;
							((float*)h_prev_pq)[i] = (float)i;
							if (h_em_deallocate) ((float*)h_em)[i] = (float)i;
						}
					}
					else if (functionality_test == 2)
					{
						// anisotropic half space
					}
					else
					{
						// fill with zeros
						memset((void*)h_pq, 0, pqNCells * (unsigned long)8);
						memset((void*)h_prev_pq, 0, pqNCells * (unsigned long)8);
						if (h_em_deallocate) memset((void*)h_em, 0, pqNCells * (unsigned long)8);
					}
					return 1;
				}
			}
		}
		return 0;
	}

	//
	// Inputs:
	// iY		Y coordinate of cross section
	// p_or_q_flag	0->p, 1->q
	//
	// Return values:
	//  0		No error
	// -1		Bad arguments.
	//
	// data array must have at least dimx * dimz items.
	//
	int XZ_Cross_Section(
		int iY,
		int p_or_q_flag,
		float* data
		)
	{
		if (iY >= 0 && iY < dimy && p_or_q_flag >= 0 && p_or_q_flag < 2)
		{
			unsigned long stride_z = 32;
			unsigned long stride_y = stride_z * (unsigned long)dimz;
			int nx = dimx / 16;
			for (int iiX = 0;  iiX < nx;  ++iiX)
			{
				int iX0 = iiX * 16;
				float* wf_src_x = h_pq + (unsigned long)iY * stride_y + (unsigned long)iiX * stride_y * (unsigned long)dimy;
				for (int iZ = 0;  iZ < dimz;  ++iZ)
				{
					float* wf_src_z = wf_src_x + iZ * stride_z;
					float* wf_dst_z = data + iZ * dimx;
					for (int iX = 0;  iX < 16;  ++iX)
					{
						wf_dst_z[iX+iX0] = wf_src_z[iX*2+p_or_q_flag];
					}
				}
			}
			return 0;
		}
		return -1;
	}

	int XY_Cross_Section(
		int iZ,
		int p_or_q_flag,
		float* data
		)
	{
		if (iZ >= 0 && iZ < dimz && p_or_q_flag >= 0 && p_or_q_flag < 2)
		{
			unsigned long stride_z = 32;
                        unsigned long stride_y = stride_z * (unsigned long)dimz;
                        int nx = dimx / 16;
                        for (int iiX = 0;  iiX < nx;  ++iiX)
                        {
                                int iX0 = iiX * 16;
                                float* wf_src_x = h_pq + (unsigned long)iZ * stride_z + (unsigned long)iiX * stride_y * (unsigned long)dimy;
                                for (int iY = 0;  iY < dimy;  ++iY)
                                {
                                        float* wf_src_y = wf_src_x + (unsigned long)iY * stride_y;
                                        float* wf_dst_y = data + (unsigned long)iY * dimx;
                                        for (int iX = 0;  iX < 16;  ++iX)
                                        {
                                                wf_dst_y[iX+iX0] = wf_src_y[iX*2+p_or_q_flag];
                                        }
                                }
                        }
                        return 0;
                }
                return -1;
	}

	void AddToPQ(
		int iX,
		int iY,
		int iZ,
		float delta_p,
		float delta_q
		)
	{
		const unsigned long stride_z = 32;
		const unsigned long stride_y = stride_z * (unsigned long)dimz;
		const unsigned long stride_x = stride_y * (unsigned long)dimy;

		unsigned long iXvol = iX / 16;
		unsigned long iiX = iX & 15;
		unsigned long idx = iXvol * stride_x + iY * stride_y + iZ * stride_z + iiX * 2;
		h_pq[idx  ] += delta_p;
		h_pq[idx+1] += delta_q;
	}

	// get multiplexed traces.
	float* Get_Muxed_Traces()
	{
		return muxed_traces;
	}

	// get number of samples per traces.
	int Get_NSAMP()
	{
		return nsamp;
	}

	//
	// Put source wavelet on GPU cards. 
	// Source wavelet has been resampled to fit internal stepping rate.
	// Set zsrcghost = -1 if you want no ghost source inserted.
	//
	void Put_Source_Wavelet(
		float* src,
		int len,
		int xsrc,
		int ysrc,
		int zsrc,
		int zsrcghost
		)
	{
		if (h_src_wavelet != 0L) free((void*)h_src_wavelet);
		h_src_wavelet = (float*)malloc(len * sizeof(float));
		memcpy((void*)h_src_wavelet, (void*)src, len * sizeof(float));		

		if (pipelines != 0L)
		{
			for (int iY = 0;  iY < num_y_pipelines;  ++iY)
			{
				pipelines[iY]->Put_Source_Wavelet(h_src_wavelet,len,xsrc,ysrc,zsrc,zsrcghost);
			}
		}
	}

	int Put_Receiver_Locations(
		int* arg_rcx,
		int* arg_rcy,
		int* arg_rcz,
		int arg_recnum,
		double dtout,			// sample rate of output traces in seconds
		double reclen,			// length of one trace in seconds
		int recghost_Flag,
		int z_freesurface,
		int& out_receiver_num_timesteps	// total number of timesteps wavefields must be propagated to fill up traces.
		)
	{
		int retval = -1;
		Free_Receiver_Locations();
		if (pipelines != 0L)
                {
                        for (int iY = 0;  iY < num_y_pipelines;  ++iY)
                        {
                                retval = pipelines[iY]->Put_Receiver_Locations(arg_rcx,arg_rcy,arg_rcz,arg_recnum,dimx,recghost_Flag,z_freesurface) ? retval : 0;
			}
		}
		if (retval != 0)
		{
			recnum = arg_recnum;
			rcx = new int[arg_recnum];
			rcy = new int[arg_recnum];
			rcz = new int[arg_recnum];
			memcpy((void*)rcx, (void*)arg_rcx, arg_recnum * sizeof(int));
			memcpy((void*)rcy, (void*)arg_rcy, arg_recnum * sizeof(int));
			memcpy((void*)rcz, (void*)arg_rcz, arg_recnum * sizeof(int));
			nsamp = (int)((reclen / dtout) + 0.5) + 1;
			muxed_traces = new float[(unsigned long)nsamp*(unsigned long)arg_recnum];
			rec_out_steps = new int[nsamp];
			time_delta = new float[nsamp];

			const int tstartsim = 0;
			const int tstartrec = 0;
			int t = tstartsim;
			for(int ntout=0, itout=0; ntout < nsamp; ++t)
			{
				float tout = tstartsim*dt + itout * dtout;
				float tfrac = t + 1.0f - tout/dt;
				if (t == 0) tfrac = 1.0f;
				//printf("t=%d, itout=%d, tout=%f, tfrac=%f\n",t,itout,tout,tfrac);
				if (tfrac > 0.0f)// && tfrac <= 1.0f)
				{
					if (t >= tstartrec)
					{
						//printf("NEW :: t = %d, tfrac = %f, tout = %f\n",t,tfrac,tout);
						rec_out_steps[ntout] = t;
						time_delta[ntout] = tfrac;
						++ntout;
					}
					++itout;
				}
			}
			/*
			int j = 0;
			for (int i = 0;  i < nsamp;  ++j)
			{
				double t = (double)j * dt;  // dt = internal timestep
				double tout = (double)i * dtout;
				if (t >= tout)
				{
					rec_out_steps[i] = j;
					time_delta[i] = (t - tout) / dt;
					++i;
				}
			}
			*/

			receiver_num_timesteps = t;
			out_receiver_num_timesteps = receiver_num_timesteps;
			mux_imgstep = new int[receiver_num_timesteps];
			mux_time_delta = new float[receiver_num_timesteps];
			for (int i = 0;  i < receiver_num_timesteps;  ++i)
			{
				mux_imgstep[i] = -1;
				mux_time_delta[i] = 0.0f;
			}
			for (int i = 0;  i < nsamp;  ++i)
                        {
				int k = rec_out_steps[i];
				mux_imgstep[k] = i;
				mux_time_delta[k] = time_delta[i];
				//printf("imaging step %d - timestep %d - delta_t = %e\n",i,k,time_delta[i]);
			}
		}
		return retval;
	}

	void Put_Sponges(
		float* h_spg_x,
		float* h_spg_y,
		float* h_spg_z
		)
	{
		if (pipelines != 0L)
		{
			for (int iY = 0;  iY < num_y_pipelines;  ++iY)
			{
				pipelines[iY]->Put_Sponges(h_spg_x,dimx,h_spg_y,dimy,h_spg_z,dimz);
			}
		}
	}

	//
	// This method re-arranges padded earth model formed for CPU kernels into acceptable format for GPU kernels.
	// Re-arranged earth model is written to h_em.
	//
	void Put_Earth_Model(
		unsigned int* PadDenAng,
		unsigned int* PadVelAnis,
		int xh,
		int yh,
		int zh,
		double dt,
		float Vel_min,
		float Vel_binsize,
		float Eps_min,
		float Eps_binsize,
		float Del_min,
		float Del_binsize,
		float Den_min,
		float Den_binsize,
		float Buoy_min,
		float Buoy_binsize,
		float Dip_min,
		float Dip_binsize,
		float Azm_min,
		float Azm_binsize,
		float Q_min,
		float Q_binsize,
		float c44c33_min,
		float c44c33_binsize
		)
	{
		unsigned long nx = dimx / 16;
		
		// DenAng and VelAnis are organized unsigned longo nx volumes of 16 by dimz by dimy cells.
		unsigned long stride_z = 16;
		unsigned long stride_y = dimz * stride_z;
		unsigned long stride_vol = stride_y * dimy;

		// PadDenAng and PadVelAnis are organized unsigned longo one volume of (dimx+2*xh) by (dimy+2*yh) by (dimz+2*zh) cells.
		unsigned long dimxh = dimx + 2 * xh;
		unsigned long dimyh = dimy + 2 * yh;

		Initialize_Earth_Model_Decompression(
			dt,
			Dip_binsize,
			Dip_min,
			Azm_binsize,
			Azm_min,
			Den_binsize,
			Den_min,
			Q_binsize,
			Q_min,
			Vel_binsize,
			Vel_min,
			c44c33_binsize,
			c44c33_min,
			Del_binsize,
			Del_min,
			Eps_binsize,
			Eps_min
			);

		for (unsigned long ii = 0;  ii < nx;  ++ii)
		{
			unsigned long iX0 = ii * 16;
			for (unsigned long iY = 0;  iY < dimy;  ++iY)
			{
				for (unsigned long iZ = 0;  iZ < dimz;  ++iZ)
				{
					for (unsigned long iiX = 0;  iiX < 16;  ++iiX)
					{
						unsigned long iX = iX0 + iiX;
						unsigned long pad_idx = (iX + xh) + (iY + yh) * dimxh + (iZ + zh) * dimxh * dimyh;
						unsigned long idx = 2 * stride_vol * ii + iX  + iY * 2 * stride_y + iZ * stride_z;
						h_em[idx] = PadDenAng[pad_idx];
						h_em[idx+stride_y] = PadVelAnis[pad_idx];
					}
				}
			}
		}
	}

	void Check_Host_Memory()
	{
                int nx = dimx / 16;
                unsigned long pqNCells = (unsigned long)nx * (unsigned long)16 * (unsigned long)dimy * (unsigned long)dimz;
		for (unsigned long i = 0;  i < pqNCells*2-10;  ++i)
		{
			//printf("pq[%ld] = %f, prev_pq[%ld] = %f\n",i,((float*)h_pq)[i],i,((float*)h_prev_pq)[i]);
			/*
			int match = 1;
			for (unsigned long ii = 0;  ii < 10;  ++ii)
			{
				float val = (float)ii;
				if (val != ((float*)h_pq)[i+ii] || val != ((float*)h_prev_pq)[i+ii] || val != ((float*)h_em)[i+ii]) match = 0;
			}

			if (match)
			{
				printf("Found match at index %ld\n",i);
				exit(0);
			}
			*/

			float val = (float)i;
			float pq = ((float*)h_pq)[i];
			float prev_pq = ((float*)h_prev_pq)[i];
			float em = ((float*)h_em)[i];
			if (pq != val)
			{
				printf("Error! h_pq[%ld] is %f, expected %f\n",i,pq,val);
				exit(0);
			}
			if (prev_pq != val)
			{
				printf("Error! h_prev_pq[%ld] is %f, expected %f\n",i,prev_pq,val);
				exit(0);
			}
			if (em != val)
			{
				printf("Error! h_em[%ld] is %f, expected %f\n",i,em,val);
				exit(0);
			}
		}
		printf("No errors found in host memory!\n");
	}

	int Allocate_Device_Resources()
	{
		int retval = 0;
		if (is_valid)
		{
			retval = 1;
			for (int iY = 0;  retval != 0 && iY < num_y_pipelines;  ++iY)
                        {
				pipelines[iY]->Allocate_Device_Resources(dimx,dimy,dimz);
			}
		}
		return retval;
	}

	unsigned long Get_H2D_Byte_Count()
	{
		unsigned long acc_count = 0;
		for (int iY = 0;  iY < num_y_pipelines;  ++iY)
		{
			acc_count += pipelines[iY]->Get_H2D_Byte_Count();
		}
		return acc_count;
	}

	unsigned long Get_D2H_Byte_Count()
	{
		unsigned long acc_count = 0;
		for (int iY = 0;  iY < num_y_pipelines;  ++iY)
		{
			acc_count += pipelines[iY]->Get_D2H_Byte_Count();
		}
		return acc_count;
	}

	//
	// Get device propagator and properties.
	// iDev is in range [0,this.Get_Device_Count()-1], latter returns total device count for all pipelines.
	//
	SingleDevicePropagator* Get_Device_Propagator(int iDev, int& is_first, int& is_last, SingleDevicePropagator*& preceeding_propagator, int& pipeline_y0, int& pipeline_y1)
	{
		if (iDev >= 0 && iDev < Get_Device_Count())
		{
			for (int iY = 0;  iY < num_y_pipelines;  ++iY)
			{
				if (iDev < pipelines[iY]->Get_Device_Count())
				{
					SingleDevicePropagator* retval = pipelines[iY]->Get_Device_Propagator(iDev,is_first,is_last,preceeding_propagator,pipeline_y0,pipeline_y1);
					return retval;
				}
				else
				{
					iDev -= pipelines[iY]->Get_Device_Count();
				}
			}
		}
		is_first = 0;
		is_last = 0;
		preceeding_propagator = 0L;
		return 0L;
	}

	int Get_Device_Count()
	{
		int total_device_count = 0;
		for (int iY = 0;  iY < num_y_pipelines;  ++iY)
		{
			total_device_count += pipelines[iY]->Get_Device_Count();
		}
		return total_device_count;
	}

	void Initialize_Earth_Model_Decompression(
		double dt,
		float dip_scaler,
		float dip_min,
		float azm_scaler,
		float azm_min,
		float density_scaler,
		float density_min,
		float Q_scaler,
		float Q_min,
		float Vel_scaler,
		float Vel_min,
		float C44C33_scaler,
		float C44C33_min,
		float Del_scaler,
		float Del_min,
		float Eps_scaler,
		float Eps_min
		)
	{
		printf("dt = %e\n",dt);
		printf("dip_scaler = %e\n",dip_scaler);
		printf("dip_min = %e\n",dip_min);
		printf("azm_scaler = %e\n",azm_scaler);
		printf("azm_min = %e\n",azm_min);
		printf("density_scaler = %e\n",density_scaler);
		printf("density_min = %e\n",density_min);
		printf("Q_scaler = %e\n",Q_scaler);
		printf("Q_min = %e\n",Q_min);
		printf("Vel_scaler = %e\n",Vel_scaler);
		printf("Vel_min = %e\n",Vel_min);
		printf("C44C33_scaler = %e\n",C44C33_scaler);
		printf("C44C33_min = %e\n",C44C33_min);
		printf("Del_scaler = %e\n",Del_scaler);
		printf("Del_min = %e\n",Del_min);
		printf("Eps_scaler = %e\n",Eps_scaler);
		printf("Eps_min = %e\n",Eps_min);
		for (int iY = 0;  iY < num_y_pipelines;  ++iY)
		{
			pipelines[iY]->Initialize_Earth_Model_Decompression(dt,dip_scaler,dip_min,azm_scaler,azm_min,density_scaler,density_min,Q_scaler,Q_min,Vel_scaler,Vel_min,C44C33_scaler,C44C33_min,Del_scaler,Del_min,Eps_scaler,Eps_min);
		}
	}

	// 
	// Send in the next two stripes (next 16+16 consecutive X locations).
	// Propagate through all the timesteps.
	// Take out resulting stripes.
	//
	int Process_Two_Stripes(int logLevel, int num_timesteps, int& output_timestep)
	{
		if (pipelines != 0L)
		{
			// move sliding window for all layers in all pipelines.
			for (int iY = 0;  iY < num_y_pipelines;  ++iY)
			{
				pipelines[iY]->Move_Sliding_Window(logLevel);
			}

			// move sliding window for host.
                        Host_Move_Sliding_Window(logLevel);

			unsigned long curr_input_host_offset1 = size_of_one_stripe_in_floats * (unsigned long)curr_input_stripe_number;
			float* h_curr_inp_pq1 = h_pq + curr_input_host_offset1;
			float* h_curr_inp_prev_pq1 = h_prev_pq + curr_input_host_offset1;
			int* h_curr_inp_em1 = h_em + curr_input_host_offset1;
				
			int output_1 = curr_output_timestep > 0 && (curr_output_timestep - Get_Total_Number_of_Timesteps()) < num_timesteps ? 1 : 0;
			unsigned long curr_output_host_offset1 = size_of_one_stripe_in_floats * (unsigned long)curr_output_stripe_number;
			float* h_curr_out_pq1 = output_1 ? h_pq + curr_output_host_offset1 : 0L;
			float* h_curr_out_prev_pq1 = output_1 ? h_prev_pq + curr_output_host_offset1 : 0L;

			//printf("input :: %d,%d - output :: %d,%d (%s)\n",curr_input_timestep,curr_input_stripe_number,curr_output_timestep,curr_output_stripe_number,output_1?"yes":"no");

                        Host_Move_Sliding_Window(logLevel);

			unsigned long curr_input_host_offset2 = size_of_one_stripe_in_floats * (unsigned long)curr_input_stripe_number;
			float* h_curr_inp_pq2 = h_pq + curr_input_host_offset2;
			float* h_curr_inp_prev_pq2 = h_prev_pq + curr_input_host_offset2;
			int* h_curr_inp_em2 = h_em + curr_input_host_offset2;
				
			int output_2 = curr_output_timestep > 0 && (curr_output_timestep - Get_Total_Number_of_Timesteps()) < num_timesteps ? 1 : 0;
			unsigned long curr_output_host_offset2 = size_of_one_stripe_in_floats * (unsigned long)curr_output_stripe_number;
			float* h_curr_out_pq2 = output_2 ? h_pq + curr_output_host_offset2 : 0L;
			float* h_curr_out_prev_pq2 = output_2 ? h_prev_pq + curr_output_host_offset2 : 0L;
			
			//printf("input :: %d,%d - output :: %d,%d (%s)\n",curr_input_timestep,curr_input_stripe_number,curr_output_timestep,curr_output_stripe_number,output_2?"yes":"no");
			for (int iDev = 0;  iDev < Get_Device_Count();  ++iDev)
			{
				int is_first, is_last, pipeline_y0, pipeline_y1;
				SingleDevicePropagator* preceeding_propagator;
				SingleDevicePropagator* propagator = Get_Device_Propagator(iDev,is_first,is_last,preceeding_propagator,pipeline_y0,pipeline_y1);
				
				if (propagator != 0L)
				{
					//propagator->LaunchSimpleCopyKernel();
					propagator->LaunchPropagationKernel(logLevel,mux_imgstep,mux_time_delta,receiver_num_timesteps);
					propagator->Start_Data_Flow(
						h_curr_inp_em1,
						h_curr_inp_em2,
						h_curr_inp_pq1,
						h_curr_inp_pq2,
						h_curr_inp_prev_pq1,
						h_curr_inp_prev_pq2,
						h_curr_out_pq1,
						h_curr_out_pq2,
						h_curr_out_prev_pq1,
						h_curr_out_prev_pq2,
						pipeline_y0,pipeline_y1,dimy,is_first,is_last,preceeding_propagator);
					propagator->Start_Receiver_Data_Flow();
					propagator->CollectReceiverValues(muxed_traces);
				}
			}
			//printf("\n\n\n");

			for (int iDev = 0;  iDev < Get_Device_Count();  ++iDev)
                        {
                                int is_first, is_last, pipeline_y0, pipeline_y1;
                                SingleDevicePropagator* preceeding_propagator;
                                SingleDevicePropagator* propagator = Get_Device_Propagator(iDev,is_first,is_last,preceeding_propagator,pipeline_y0,pipeline_y1);

				if (propagator != 0L)
                                {
					propagator->WaitForInputStream();
					propagator->WaitForOutputStream();
					propagator->WaitForComputeStream();
					propagator->FlipFlopFlap();
				}
			}

			output_timestep = curr_output_timestep;
			return curr_output_timestep - Get_Total_Number_of_Timesteps() >= num_timesteps ? 1 : 0;
		}

		output_timestep = -1;
		return 1;
	}

	void CompleteReceiverValues()
	{
		// FLOP :: GPU->Host 
		// FLAP :: CollectReceiverValues
		for (int iDev = 0;  iDev < Get_Device_Count();  ++iDev)
		{
			int is_first, is_last, pipeline_y0, pipeline_y1;
			SingleDevicePropagator* preceeding_propagator;
			SingleDevicePropagator* propagator = Get_Device_Propagator(iDev,is_first,is_last,preceeding_propagator,pipeline_y0,pipeline_y1);

			if (propagator != 0L)
			{
				propagator->Start_Receiver_Data_Flow();
				propagator->CollectReceiverValues(muxed_traces);
			}
		}
	
		// FLOP :: Wait for GPU->Host transfer to finish
		for (int iDev = 0;  iDev < Get_Device_Count();  ++iDev)
		{
			int is_first, is_last, pipeline_y0, pipeline_y1;
			SingleDevicePropagator* preceeding_propagator;
			SingleDevicePropagator* propagator = Get_Device_Propagator(iDev,is_first,is_last,preceeding_propagator,pipeline_y0,pipeline_y1);

			if (propagator != 0L)
			{
				propagator->WaitForOutputStream();
				propagator->FlipFlopFlap();
			}
		}

		// FLOP :: CollectReceiverValues
		for (int iDev = 0;  iDev < Get_Device_Count();  ++iDev)
                {
                        int is_first, is_last, pipeline_y0, pipeline_y1;
                        SingleDevicePropagator* preceeding_propagator;
                        SingleDevicePropagator* propagator = Get_Device_Propagator(iDev,is_first,is_last,preceeding_propagator,pipeline_y0,pipeline_y1);

                        if (propagator != 0L)
                        {
                                propagator->CollectReceiverValues(muxed_traces);
                        }
                }
	}

	double Compute_Max_Memory_Required_GB()
	{
		double max_mem_GB = 0.0;
		if (is_valid)
		{
			for (int iY = 0;  iY < num_y_pipelines;  ++iY)
			{
				double mem_GB = pipelines[iY]->Compute_Max_Memory_Required_GB();
				if (mem_GB > max_mem_GB) max_mem_GB = mem_GB;
			}
		}
		return max_mem_GB;
	}

	double Compute_Max_Bytes_Per_Cell_Per_Timestep()
	{
		double retval = 0.0;
		if (is_valid)
		{
			for (int iY = 0;  iY < num_y_pipelines;  ++iY)
                        {
				double bytes_per_cell_per_timestep = pipelines[iY]->Compute_Max_Bytes_Per_Cell_Per_Timestep();
				if (bytes_per_cell_per_timestep > retval) retval = bytes_per_cell_per_timestep;
			}
		}
		return retval;
	}

	//
	// Get the total number of timesteps. This is the depth of the pipeline(s).
	//
	int Get_Total_Number_of_Timesteps()
	{
		int max_ts = 0;
		for (int iY = 0;  iY < num_y_pipelines;  ++iY)
		{
			int curr_ts = pipelines[iY]->Get_Total_Number_of_Timesteps();
			if (curr_ts > max_ts) max_ts = curr_ts;
		}
		return max_ts;	
	}

	// wavefield contains this many stripes.
	int Total_Number_of_Stripes()
	{
		return total_num_stripes;
	}

	void Print()
	{
		if (is_valid)
		{
			for (int iPipeline = 0;  iPipeline < num_y_pipelines;  ++iPipeline)
			{
				printf("\nT R A C K   %d   -   y=[%d,%d]\n\n",iPipeline+1,pipelines[iPipeline]->Get_Y0(),pipelines[iPipeline]->Get_Y1());
				pipelines[iPipeline]->Print();
			}
		}
		else
		{
			printf("Error! Invalid device configuration!\n");
		}
	}

private:
	void Host_Move_Sliding_Window(int logLevel)
	{
		++curr_input_slice;

		curr_input_stripe_number = (curr_input_slice + total_num_stripes) % total_num_stripes;
		curr_input_timestep = (((curr_input_slice + total_num_stripes) / total_num_stripes) - 1) * Get_Total_Number_of_Timesteps();

		curr_output_stripe_number = (curr_input_slice - pipelines[0]->Get_X_Delay() + 2 * total_num_stripes) % total_num_stripes;
		curr_output_timestep = (((curr_input_slice - pipelines[0]->Get_X_Delay() + 2 * total_num_stripes) / total_num_stripes) - 1) * Get_Total_Number_of_Timesteps();

		if (logLevel >= 5)
		{
			printf(">>> input = ts:%d,%d -- output = ts:%d,%d\n",curr_input_timestep,curr_input_stripe_number,curr_output_timestep,curr_output_stripe_number);
		}
	}

	void Free_Receiver_Locations()
	{
		if (rcx != 0L) delete [] rcx;
		if (rcy != 0L) delete [] rcy;
		if (rcz != 0L) delete [] rcz;
		if (muxed_traces != 0L) delete [] muxed_traces;
		if (rec_out_steps != 0L) delete [] rec_out_steps;
		if (time_delta != 0L) delete [] time_delta;
		if (mux_imgstep != 0L) delete [] mux_imgstep;
		if (mux_time_delta != 0L) delete [] mux_time_delta;
		nsamp = 0;
		recnum = 0;
	}

	void Create_Configuration(int Kernel_Type, int* pipe_dimy, int num_y_pipelines)
	{
		if (pipelines != 0L)
		{
			for (int iY = 0;  iY < num_y_pipelines;  ++iY)
			{
				if (pipelines[iY] != 0L) delete pipelines[iY];
			}
			delete [] pipelines;
		}
		if (shared != 0L) delete [] shared;

		int devices_per_y_pipeline = single_device_multi_track ? 1 : device_count / num_y_pipelines;
		//int num_y_per_pipeline = (dimy + num_y_pipelines - 1) / num_y_pipelines;
		
		int y0 = 0;
		pipelines = new SinglePipelinePropagator*[num_y_pipelines];
		for (int iY = 0;  iY < num_y_pipelines;  ++iY)
		{
			int y1 = y0 + pipe_dimy[iY] - 1;
			int exclude_halo_ylo = (y0 == 0 ? 1 : 0);
			int exclude_halo_yhi = 0;
			if (y1 >= (dimy-1))
			{
				y1 = dimy - 1;
				exclude_halo_yhi = 1;
			}

			int device_id = single_device_multi_track ? first_device : first_device+devices_per_y_pipeline*iY;

			pipelines[iY] = new SinglePipelinePropagator(
				Kernel_Type,
				dimx,y0,y1,dimy,dimz,
				device_id,devices_per_y_pipeline,
				exclude_halo_ylo,exclude_halo_yhi,
				num_timesteps_per_device,halo_ny,OTflag,dx,dy,dz);

			y0 = y1 + 1;
		}
		int nstep = pipelines[0]->Get_Total_Number_of_Timesteps();
		if ((nstep & 1) == 0)
		{
			// even number of steps
			curr_input_slice = -1;
		}
		else
		{
			// odd number of steps
			curr_input_slice = -2;
		}

		if (single_device_multi_track)
		{
			for (int iLayer = 0;  iLayer < pipelines[0]->props[0]->Get_Number_of_Layers();  ++iLayer)
			{
				shared = new TimeStepLayer*[num_y_pipelines];
				for (int iY = 0;  iY < num_y_pipelines;  ++iY)
	                        {
					shared[iY] = pipelines[iY]->props[0]->layers[iLayer];
					pipelines[iY]->props[0]->layers[iLayer]->num_layers_sharing_device_resources = num_y_pipelines;
					pipelines[iY]->props[0]->layers[iLayer]->Layers_Sharing_Device_Resources = shared;
				}
			}
		}
	}

	int is_valid;

	int Kernel_Type;
	int dimx;
	int dimy;
	int dimz;
	int num_y_pipelines;
	int first_device;
	int device_count;
	int num_timesteps_per_device;
	int halo_ny;
	int OTflag;

	float dx;
	float dy;
	float dz;
	float dt;

	int* pipe_dimy;
	TimeStepLayer** shared;

	int single_device_multi_track;
	SinglePipelinePropagator** pipelines;

	int h_pinned;		// 0->regular, swappable memory, 1->pinned memory
	float* h_pq;		// p and q wavefields interleaved
	float* h_prev_pq;	// p and q wavefields from previous timestep

	int h_em_deallocate;	// de-allocate h_em if this field is true (non-zero).
	int* h_em;		// compressed earth model

	float* h_src_wavelet;		// source wavelet

	int* rcx;			// receiver coordinates (in cell indexes)
	int* rcy;
	int* rcz;
	int recnum;			// number of receiver locations
	int nsamp;			// number of samples in output traces
	float* muxed_traces;		// multiplexed traces, nsamp * recnum values.
	int* rec_out_steps;		// list of all the receiver output timesteps.
	float* time_delta;		// time delta in fractions of dt_out (output sampling rate).
	int receiver_num_timesteps;	// minimum number of timesteps wavefields must be propagated to fill up receiver traces.
	int* mux_imgstep;		// dimension is total_num_timesteps. Contains offset into muxed_traces array for each timestep. Steps that are not output steps have this set to -1.
	float* mux_time_delta;		// dimension is total_num_timesteps. Contains time delta used for linear interpolation for each timestep. Steps that are not output steps have this set to 0.

	int total_num_stripes;		// total number of stripes. == dimx / 16;
	unsigned long size_of_one_stripe_in_floats;

	int curr_input_slice;	// current absolute input slice. 0 -> 0th slice of 0th timestep. total_num_stripes -> 0th slice of 1st timestep, 2*total_num_stripes -> 0th slice of 2nd timestep etc.

	int curr_input_stripe_number;	// current stripe number in range [0,total_num_stripes>
	int curr_input_timestep;	// current input timestep. After total_num_stripes have been processed, this value is incremented by one and curr_input_stripe_number is reset to 0.

	int curr_output_stripe_number;
	int curr_output_timestep;
};

/***** swap2bytes ba --> ab *****/
void swap2bytes(short *i2, int n)
{
        int i;
        short a,b;
        for (i=0; i<n; i++)
        {
                a = i2[i] << 8;
                b = (i2[i] >> 8) & 255;
                i2[i] = a | b;
        }
}


/***** swap4bytes:  dcba --> abcd *****/
void swap4bytes(int *i4, int n)
{
  int k, i, a, b, c, d, bmask = 16711680, cmask = 65280, dmask = 255;
  for(k=0; k<n; k++)
  { i = i4[k];
    a =  i << 24;          b = (i << 8)  & bmask;
    c = (i >> 8) & cmask;  d = (i >> 24) & dmask;
    i4[k] = a | b | c | d ;
  }
}

/***** CONVERT ENTIRE SHOT PROFILE TO SEGY & WRITE TO DISK
       Transpose 3D shot profile output from p[t][y][x] to p[y][x][t]
       (i.e., from planes of constant time to trace columns);
       Optionally swap data & header bytes from linux little endian
       (for ProMAX readability if data and/or headers need byte swapping); *****/
void writeSEGY(float* seisdata, char *segyfilename, char *hdrstring, int swapflag,
                int ffid, float srcx, float srcy, float srcz,
                float dt, float dtout, float timestartrec, int nsamp, int nrec,
                float recxstart, int nx_inline, float dx_inline,
                float recystart, int ny_crossline, float dy_crossline,
                float reczstart, int nzrec, float dzrec)
{
        /* SEGY DATA TYPES ***/
        char reel_id_hdr1[3200];
        memset((void*)reel_id_hdr1, 0, 3200);

        struct
        {
                int jobid;
                int lineid;
                int reelid;
                short ntrc_per_record;
                short nauxtrc;
                short dtreel;
                short dtfield;
                short nsampreel;
                short nsampfield;
                short datafmt;
                short cmpfold;
                short sortcode;
                char skip[370];
        } reel_id_hdr2;
        memset((void*)&reel_id_hdr2, 0, sizeof(reel_id_hdr2));

        struct
        {
                int trcseqno;
                int skip0;
                int isrc;
                int ichan;
                int skip1;
                int cmpbin;
                int trcensemb;
                short code;
                char skip3[6];
                int offset;
                int recelev;
                int elevatsrc;
                int srcdepth;
                char skip4[16];
                short scalar1;
                short scalar2;
                int srcx;
                int srcy;
                int recx;
                int recy;
                short lenunit;
                char skip5[18];
                short tstartrec;
                char skip6[4];
                short nsamp;
                short dtmicro;
                char skip7[82];
                float cmp_x;
                float cmp_y;
                int iline_no;
                int xline_no;
                float xoff;
                float yoff;
                float azim;
                char skip8[12];
        } trc_id_hdr;
        memset((void*)&trc_id_hdr, 0, sizeof(trc_id_hdr));
        /* cmp_x starts at byte position 201 */

        short one2=1, five2=5, nrec2, dtmicro2, nsamp2, fold2, trc_sortcode2, tstartrec2;
        short neg100 = -100;
        int one=1, one4=1, elevatsrc, recelev, srcdepth, xsrc, ysrc, xrec, yrec, izrec;
        int ichn, ichan, trcseq, trcens;
        int xline, iline, xinline, ycrossline, offset;
        int t, trcseqno, trcensemb, dtmicro;
        long zxyline;
        long ynx, currtrace;
        long ptr, nsegyvol;
        float recx, recy, recz, cmpx, cmpy, xoff, yoff, azim;
        float *segyvol, *trace;
        void swap2bytes(short*, int), swap4bytes(int*, int);
        FILE *segyfile;
        reel_id_hdr2.jobid = 0;
        trc_id_hdr.trcseqno = 0;

        /*** OPEN OUTPUT FILE ***/
        if( !(segyfile = fopen(segyfilename,"w")) )
        {
                fprintf(stderr,"Cannot open %s\n",segyfilename);
                exit(0);
        }

        /*** ALLOCATE MEM ***/
        trace = (float *)malloc(nsamp*sizeof(float));
        nsegyvol = nrec*(nsamp*sizeof(float) + 240) + 3600;  /* 240=trc header, 3600=main+reel */
        segyvol = (float *)malloc(nsegyvol);  /* Rec array */

        /*** FILL REEL ID HEADER 1 ***/
        /* Write Time of Day */
        /*status = gettimeofday(&time, &tz);
          secs = time.tv_sec;
          timenow = ctime(&secs);
          strcat(reel_id_hdr1, timenow);
         */
        strcpy(reel_id_hdr1, hdrstring);
        memcpy(segyvol, &reel_id_hdr1, 3200);
        ptr = 800;  /* 800 = number of 4-byte words in 3200 byte master header */

        /*** FILL REEL ID HEADER 2 ***/
        dtmicro = (int)(1000000.*dtout + 0.5);
        trc_sortcode2 = 1;  /* as recorded, no sorting */
        fold2 = 1;
        one2 = one;       one4 = one;
        nrec2 = nrec;     dtmicro2 = dtmicro;        nsamp2 = nsamp;

        if(swapflag)
        {
                swap2bytes(&one2, 1);        swap2bytes(&five2, 1);  swap4bytes(&one4, 1);
                swap2bytes(&nrec2, 1);       swap2bytes(&dtmicro2, 1);
                swap2bytes(&nsamp2, 1);      swap2bytes(&trc_sortcode2, 1);
                swap2bytes(&fold2, 1);
        }
        reel_id_hdr2.jobid = reel_id_hdr2.lineid = reel_id_hdr2.reelid = one4;
        reel_id_hdr2.ntrc_per_record = nrec2;
        reel_id_hdr2.dtreel = reel_id_hdr2.dtfield = dtmicro2;
        reel_id_hdr2.nsampreel = reel_id_hdr2.nsampfield = nsamp2;
        reel_id_hdr2.datafmt = five2;
        reel_id_hdr2.cmpfold = fold2;
        reel_id_hdr2.sortcode = trc_sortcode2;
        memcpy(&segyvol[ptr], &reel_id_hdr2, 400);
        ptr += 100;  /* 100 = number of 4-byte words in 400 byte master header */


        /*** FILL SOURCE-RELATED PART OF TRACE HEADER ***/
        elevatsrc = 0;
        srcdepth = (int)(100.*srcz);
        xsrc = (int)(100.*srcx);
        ysrc = (int)(100.*srcy);
        tstartrec2 = (int)(timestartrec*1000. + 0.5);
        if(swapflag)
        {
                swap4bytes(&ffid, 1);
                swap4bytes(&elevatsrc, 1);     swap4bytes(&srcdepth, 1);
                swap4bytes(&xsrc, 1);          swap4bytes(&ysrc, 1);
                swap2bytes(&tstartrec2, 1);
                swap2bytes(&neg100, 1);
        }
        trc_id_hdr.isrc = ffid;
        trc_id_hdr.elevatsrc = elevatsrc; trc_id_hdr.srcdepth = srcdepth;
        trc_id_hdr.srcx = xsrc;           trc_id_hdr.srcy = ysrc;
        trc_id_hdr.nsamp = nsamp2;
        trc_id_hdr.tstartrec = tstartrec2;
        trc_id_hdr.dtmicro = dtmicro2;
        trc_id_hdr.scalar1 = neg100;
        trc_id_hdr.scalar2 = neg100;

	const float r2d = 57.29577951308232311f;


        /*** READ IN SEISMIC DATA (slow axis is Time, med axis is Y, fast axis is X)
          AND WRITE OUT TRACE HEADER + DATA TRACE ***/
        for(izrec = 0, trcseqno=ichn=0; izrec<nzrec; izrec++)
        {
                zxyline = (long)izrec *  (long)nx_inline * (long)ny_crossline;
                recz = reczstart + izrec*dzrec;
                recelev = -(int)(100.*recz);  if(swapflag) swap4bytes(&recelev, 1);
                trc_id_hdr.recelev = recelev;

                for(ycrossline=0; ycrossline<ny_crossline; ycrossline++)
                {
                        ynx = zxyline + (long)ycrossline * (long)nx_inline;
                        recy = recystart + ycrossline*dy_crossline;
                        yrec = (int)(100.*recy);  xline = ycrossline+1;
                        if(swapflag) { swap4bytes(&yrec, 1); swap4bytes(&xline, 1); }
                        trc_id_hdr.recy = yrec;
                        trc_id_hdr.iline_no = xline; /* yes, this is correct */

                        for(xinline=0, trcensemb=1;    xinline<nx_inline;   xinline++, trcensemb++)
                        {
                                currtrace = (long)ynx + (long)xinline;

                                trcseq = trcseqno++;       ichan = ichn++;      trcens = trcensemb;
                                recx = recxstart + xinline*dx_inline;
                                xrec = (int)(100.*recx);
                                xoff = recx - srcx;        yoff = recy - srcy;
                                cmpx = 0.5*(srcx + recx);  cmpy = 0.5*(srcy + recy);
                                iline = xinline+1;
                                offset = (int)(sqrtf(yoff*yoff + xoff*xoff) + 0.5);
                                azim = r2d*atan2f(yoff, xoff);

                                if(swapflag)
                                {
                                        swap4bytes(&trcseq, 1);  swap4bytes(&ichan, 1);  swap4bytes(&trcens, 1);
                                        swap4bytes(&xrec, 1);
                                        swap4bytes((int*)(&cmpx), 1); swap4bytes((int*)(&cmpy), 1);
                                        swap4bytes(&iline, 1);
                                        swap4bytes((int*)(&xoff), 1); swap4bytes((int*)(&yoff), 1);
                                        swap4bytes(&offset, 1);       swap4bytes((int*)(&azim), 1);
                                }

                                /* Assign & Write Trace Header */
                                trc_id_hdr.trcseqno = trcseq;
                                trc_id_hdr.ichan = ichan;
                                trc_id_hdr.trcensemb = trcens;
                                trc_id_hdr.offset = offset;
                                trc_id_hdr.recx = xrec;
                                trc_id_hdr.cmp_x = cmpx;
                                trc_id_hdr.cmp_y = cmpy;
                                trc_id_hdr.xline_no = iline; /* yes, this is correct */
                                trc_id_hdr.xoff = xoff;
                                trc_id_hdr.yoff = yoff;
                                trc_id_hdr.azim = azim;
                                memcpy(&segyvol[ptr], &trc_id_hdr, 240);
                                ptr += 60;  /* 60 = number of 4-byte words in 240 byte trace header */

                                /* Read one trace into trace[] array and swapbytes */
                                for(t=0; t<nsamp; t++)  trace[t] = seisdata[currtrace + (long)t * (long)nrec];
                                if(swapflag)  swap4bytes((int*)trace,nsamp);

                                /* Write Trace Vector to memory */
                                memcpy(&segyvol[ptr], trace, nsamp*sizeof(float));
                                ptr += nsamp;
                        }
                }
        }

        /* WRITE EVERYTHING TO DISK IN ONE GO */
        fwrite(segyvol, nsegyvol, 1, segyfile);

        fclose(segyfile);

        free((void*)trace);
        free((void*)segyvol);
}

void src(int logLevel, float dt, float fmax, int type, char* stfname, int* tsrc, float* stf)
{
        if (logLevel >= 4)
        {
                printf("src(logLevel=%d, dt=%e, fmax=%e, type=%d, stfname=%s, tsrc=%s, stf=%s)\n",logLevel,dt,fmax,type,stfname!=0L?stfname:"nil",tsrc!=0L?"ok":"nil",stf!=0L?"ok":"nil");
                fflush(stdout);
        }
        if (type == 1)
        {
                /* type==1: first derivative of a Gaussian, with linear extension tapers */
                int i, imax, ntap;
                float w0, wmax, ts,t,rt, diff;

                wmax = 2.0f*M_PI*fmax;
                /* Note:   tsrc = ts/dt = 2/(gam*khmax)*(rt*rt)*(Vmax/Vmin) */

                /* if(type==1) */  /* only one type for now */
                {
                        rt = 3.571625f; /* guarantees SourceAmplitude(tmax)   = 0.01*MaxAmplitude
                                           and guarantees SpectralAmplitude(fmax) = 0.01*MaxSpectrum
                                           and: wmax/w0 = rt = w0*ts/2  */
                        w0 = wmax/rt;  /* w0i = 1./w0; */
                        ts = 2.0f*rt/w0;  /* total source time */
                        imax = (int)(ts/dt) + 1;

                        for(i=0;i<imax;i++)
                        { t=i*dt-0.5f*ts;
                                stf[i] = -sqrtf(M_E)*w0*t*exp(-0.5f*w0*w0*t*t);
                        }

                        /* taper (linearly extend) front and back ends */
                        /* front end */
                        diff = stf[1]-stf[0];
                        ntap = (int)(fabs(stf[0]/diff));
                        for(i=imax-1; i>=0; i--) stf[i+ntap] = stf[i];  /* shift */
                        for(i=ntap-1; i>=0; i--) stf[i] = stf[i+1] - diff; /* taper */
                        imax += ntap;

                        /* back end: */
                        diff = stf[imax-1]-stf[imax-2];
                        ntap = (int)(fabs(stf[imax-1]/diff));
                        for(i=0; i<ntap; i++)  stf[imax+i] = stf[imax+i-1] + diff; /* taper */
                        imax += ntap;
                }

                *tsrc = imax;

                if (logLevel >= 4) printf("SOURCE TYPE 1 : imax=%d\n",imax);

                // for(i=0; i<imax; i++) stf[i]=0.0f; stf[0] = 1.0f;

                // for(i=0; i<imax; i++) printf("%d  %f\n", i,stf[i]);
        }
        else if (type == 2)
        {
                /* type==2: user provided source function from file */
                int nfine;
                float dtfine;
                FILE* stffile = fopen(stfname,"r");
                if (stffile == 0L)
                {
                        printf("ERROR! src(...) - Source wavelet file '%s' cannot be read.\n",stfname);
                        fflush(stdout);
                        exit(-1);
                }
                fscanf(stffile,"%d %f", &nfine, &dtfine);
                float* stffine = (float*)malloc((nfine+1)*sizeof(float));
                for(int i=0; i<nfine; i++) fscanf(stffile,"%f", &stffine[i]);
                stffine[nfine] = 0.;

                int imax = (int)((nfine-1)*dtfine/dt) + 1;
                float absmax = -1e37f;
                for(int i=0; i<imax; i++)
                {
                        float t = i*dt;
                        int tfine = (int)(t/dtfine);
                        float frac = t/dtfine - tfine;
                        float val = (1.-frac)*stffine[tfine] + frac*stffine[tfine+1];
                        stf[i] = val;
                        float absval = val < 0.0f ? -val : val;
                        if (absval > absmax) absmax = absval;
                }
                *tsrc = imax;

                if (logLevel >= 4) printf("SOURCE TYPE 2 : nfine=%d, dtfine=%e, imax=%d, absmax=%e\n",nfine,dtfine,imax,absmax);

                for(int i=0; i<imax; i++)
                {
                        stf[i] /= absmax;
                        if (logLevel >= 4) printf("stf[%d] = %e\n",i,stf[i]);
                }
        }
}

int main(int argc, char* argv[])
{
	//
	// Determine number of devices and create one thread for each device.
	// Use OpenMP for threading.
	//
	int debug = 0;

	printf("\nAcoustic TTI FDTD GPU v0.9\n\n");
	if (argc != 6)
	{
		printf("Usage: %s <parmfile> <first_device> <device_count> <num_pipelines> <timesteps_per_device>\n",argv[0]);
		return 0;
	}

	char hostname[4096];
	gethostname(hostname,4096);
	printf("Running on host %s\n",hostname);

	int num_y_pipelines = 1;
	int max_timesteps_per_device = 4;

#ifdef NO_GPU
	printf("\nNO_GPU OPTION!\n\n");
	int first_device = 0;
	int device_count = 1;
#else
	int first_device = 0;
	int device_count = 0;
	cudaGetDeviceCount(&device_count);
	float min_dev_glob_mem_gb = 0.0f;
#endif
	if (device_count == 0)
	{
		printf("No CUDA capable devices detected.\n");
		return -1;
	}
	else
	{
		int requested_device_count = -1, requested_first_device = -1, requested_pipelines = -1, requested_max_timesteps_per_device = -1;
		if (argc == 6)
		{
			requested_first_device = atoi(argv[2]);
			requested_device_count = atoi(argv[3]);
			requested_pipelines = atoi(argv[4]);
			requested_max_timesteps_per_device = atoi(argv[5]);
			if (requested_max_timesteps_per_device < 0)
			{
				printf("ENABLING DEBUG CROSSPLOT OUTPUTS - LONGER RUNTIME.\n");
				requested_max_timesteps_per_device = -requested_max_timesteps_per_device;
				debug = 1;
			}
			if (
					requested_device_count >  0 && requested_device_count <= device_count && 
					requested_first_device >= 0 && requested_first_device <  device_count &&
					requested_pipelines    >  0 && requested_pipelines    <= device_count && (device_count % requested_pipelines) == 0 &&
					requested_max_timesteps_per_device > 0
			   )
			{
				printf("User Selected :: First Device = %d, Device Count = %d, Pipelines = %d, MAX Timesteps per Device = %d\n",requested_first_device,requested_device_count,requested_pipelines,requested_max_timesteps_per_device);
				first_device = requested_first_device;
				device_count = requested_device_count;
				num_y_pipelines = requested_pipelines;	
				max_timesteps_per_device = requested_max_timesteps_per_device;
			}
			else
			{
				printf("ERROR! Invalid device selection (requested first device = %d, requested device count = %d)\n",requested_first_device,requested_device_count);
				return -1;
			}
		}
		else
		{
			num_y_pipelines = device_count;
			printf("Using all %d devices with %d pipelines.\n",device_count,num_y_pipelines);
		}

#ifndef NO_GPU
		for (int i = 0;  i < device_count;  ++i)
		{
			int device = i + first_device;
			cudaDeviceProp devProps;
			cudaGetDeviceProperties(&devProps, device);
			printf("device %d :: %s, CC=%d.%d, asyncEngineCount=%d, multiProcessorCount=%d, memoryBusWidth=%d, memoryClockRate=%d, sharedMemPerBlock=%d, totalGlobalMem=%ld\n",device,devProps.name,devProps.major,devProps.minor,devProps.asyncEngineCount,devProps.multiProcessorCount,devProps.memoryBusWidth,devProps.memoryClockRate,devProps.sharedMemPerBlock,devProps.totalGlobalMem);
			float glob_mem_gb = (float)devProps.totalGlobalMem * 1e-9f;
			if (min_dev_glob_mem_gb == 0.0f || glob_mem_gb < min_dev_glob_mem_gb) min_dev_glob_mem_gb = glob_mem_gb;
		}
#endif
	}

	Parmfile_Reader* parmfile = new Parmfile_Reader(argv[1]);

	const int sky_z = 0;
	const float sky_factor = 1.0f;

	const int absorb_z = parmfile->absorbz0;

	int Log_Level = 4;
	Acoustic_Earth_Model* earth_model = new Acoustic_Earth_Model(
			Log_Level,
			parmfile->Kernel_Type,
			1,                                                      // Target_Type == GPU
			parmfile->vpname,
			parmfile->dnname,
			parmfile->epsetaname,
			parmfile->deltaname,
			parmfile->dipdxname,
			parmfile->azmdyname,
			parmfile->Qname,
			parmfile->eta_Flag,
			parmfile->dipxdipy_Flag,
			parmfile->degrees_Flag,
			parmfile->swap_Flag,
			parmfile->VsoVp0,
			parmfile->newfmax,
			parmfile->fast_Axis,
			parmfile->med_Axis,
			parmfile->slow_Axis,
			parmfile->nx,
			parmfile->ny,
			parmfile->nz,
			parmfile->sub_xoff,
			parmfile->sub_xoff + parmfile->sub_nx - 1,
			parmfile->sub_yoff,
			parmfile->sub_yoff + parmfile->sub_ny - 1,
			parmfile->sub_zoff,
			parmfile->sub_zoff + parmfile->sub_nz - 1,
			sky_z,
			sky_factor,
			absorb_z,
			parmfile->dh,
			parmfile->dz
			);
	if (!earth_model->Is_Valid())
	{
		printf("INVALID EARTH MODEL - Exiting\n");
		return -1;
	}
	if (parmfile->Kernel_Type == 0)
	{
		printf("Isotropic GPU kernel not supported. Use VTI kernel with anisotropic parameters set to zero.\nExiting!\n");
		return -2;
	}

	int OTflag = 2;
	float fq = 0.0f;
	earth_model->Create_Compressed_Earth_Model(OTflag,parmfile->gamfac,fq,parmfile->dtout);

	float dx = parmfile->dh;
	float dy = parmfile->dh;
	float dz = parmfile->dz;
	double dt = earth_model->Get_DT();

	const int logLevel = 4;

	int num_timesteps_per_device = 1;
	const int dimx = earth_model->Get_NX();
	const int dimy = earth_model->Get_NY();
	const int dimz = earth_model->Get_NZ() + sky_z + absorb_z;
	const int halo_ny = 9;

	printf("NX = %d, Actual_NX = %d\n",dimx,earth_model->Get_Actual_NX());
	printf("NY = %d, Actual_NY = %d\n",dimy,earth_model->Get_Actual_NY());
	printf("NZ = %d, Actual_NZ = %d\n",dimz,earth_model->Get_Actual_NZ());

	unsigned long volSize = (unsigned long)dimx * (unsigned long)dimy * (unsigned long)dimz;

	DevicePropagator* prop = 0L;
	DevicePropagator* new_prop = 0L;
	for (
			new_prop = new DevicePropagator(parmfile->Kernel_Type,dimx,dimy,dimz,first_device,device_count,num_y_pipelines,num_timesteps_per_device,halo_ny,OTflag,dx,dy,dz,dt);
			num_timesteps_per_device <= max_timesteps_per_device && new_prop->Compute_Max_Memory_Required_GB() < (min_dev_glob_mem_gb * 0.9f) && (prop == 0L || new_prop->Compute_Max_Bytes_Per_Cell_Per_Timestep() < prop->Compute_Max_Bytes_Per_Cell_Per_Timestep());
			new_prop = new DevicePropagator(parmfile->Kernel_Type,dimx,dimy,dimz,first_device,device_count,num_y_pipelines,num_timesteps_per_device,halo_ny,OTflag,dx,dy,dz,dt)
	    )
	{
		if (prop != 0L) delete prop;
		prop = new_prop;
		printf("num_timesteps_per_device = %d -> %.2f GB required, bytes per cell per timestep = %.2f\n",num_timesteps_per_device,prop->Compute_Max_Memory_Required_GB(),prop->Compute_Max_Bytes_Per_Cell_Per_Timestep());
		++num_timesteps_per_device;
	}
	if (new_prop != 0L) delete new_prop;
	if (prop != 0L)
	{
		prop->Print();
		//return 0;

		int *PadVelAnis, *PadDenAng;
		int ***VelAnis, ***DenAng;
		earth_model->Get_Compressed_Earth_Model(PadVelAnis,PadDenAng,VelAnis,DenAng);

		prop->Initialize_Earth_Model_Decompression(
			earth_model->Get_DT(),
			earth_model->Get_Dip_Scaler(),
			earth_model->Get_Dip_Min(),
			earth_model->Get_Azm_Scaler(),
			earth_model->Get_Azm_Min(),
			earth_model->Get_Den_Scaler(),
			earth_model->Get_Den_Min(),
			earth_model->Get_Q_Scaler(),
			earth_model->Get_Q_Min(),
			earth_model->Get_Vp_Scaler(),
			earth_model->Get_Vp_Min(),
			parmfile->VsoVp0 * parmfile->VsoVp0,
			0.0f,
			earth_model->Get_Dta_Scaler(),
			earth_model->Get_Dta_Min(),
			earth_model->Get_Eps_Scaler(),
			earth_model->Get_Eps_Min()
			);

		prop->Allocate_Device_Resources();
		printf("Device memory allocated\n");

		float* spgx = new float[dimx];
		float* spgy = new float[dimy];
		float* spgz = new float[dimz];
		Compute_Sponges(
				parmfile->spongecoeff_x,
				parmfile->spongecoeff_y,
				parmfile->spongecoeff_z_lo,
				parmfile->spongecoeff_z_hi,
				parmfile->spongewidth_x,
				parmfile->spongewidth_y,
				parmfile->spongewidth_z_lo,
				parmfile->spongewidth_z_hi,
				parmfile->absorbz0_Flag,
				dimx,
				earth_model->Get_NX() - earth_model->Get_Actual_NX(),
				dimy,
				earth_model->Get_NY() - earth_model->Get_Actual_NY(),
				dimz,
				earth_model->Get_NZ() - earth_model->Get_Actual_NZ(),
				spgx,
				spgy,
				spgz
			       );
		prop->Put_Sponges(spgx,spgy,spgz);
		delete [] spgx;
		delete [] spgy;
		delete [] spgz;

		prop->Allocate_Host_Memory(1,0,PadDenAng);
		printf("Host memory allocated\n");

		int xsrc = (int)round(parmfile->srcx / parmfile->dh);
		int ysrc = (int)round(parmfile->srcy / parmfile->dh);
		int zsrc = (int)round(parmfile->srcz / parmfile->dz) + sky_z + absorb_z;
		int zsrcghost = parmfile->srcghost_Flag ? 2 * (sky_z + absorb_z) - zsrc : -1;

		double dtout = parmfile->dtout;
		double reclen = parmfile->maxtime;

		int num_timesteps = (int)((reclen / dtout) + 0.5);

		int nrecx = 0, nrecy = 0, nrecz = 0;
		for (int irecx = parmfile->xrecstart;  irecx <= parmfile->xrecend;  irecx += parmfile->xrecstride) ++nrecx;
		for (int irecy = parmfile->yrecstart;  irecy <= parmfile->yrecend;  irecy += parmfile->yrecstride) ++nrecy;
		for (int irecz = parmfile->zrecstart;  irecz <= parmfile->zrecend;  irecz += parmfile->zrecstride) ++nrecz;

		int recnum = nrecx * nrecy * nrecz;
		printf("nrecx=%d, nrecy=%d, nrecz=%d, recnum=%d\n",nrecx,nrecy,nrecz,recnum);
		int* rcx = new int[recnum];
		int* rcy = new int[recnum];
		int* rcz = new int[recnum];
		for (int irecz = parmfile->zrecstart, idx = 0;  irecz <= parmfile->zrecend;  irecz += parmfile->zrecstride)
		{
			for (int irecy = parmfile->yrecstart;  irecy <= parmfile->yrecend;  irecy += parmfile->yrecstride)
			{
				for (int irecx = parmfile->xrecstart;  irecx <= parmfile->xrecend;  irecx += parmfile->xrecstride)
				{
					rcx[idx] = irecx;
					rcy[idx] = irecy;
					rcz[idx] = irecz + sky_z + absorb_z;
					++idx;
				}
			}
		}
		prop->Put_Receiver_Locations(rcx,rcy,rcz,recnum,dtout,reclen,parmfile->recghost_Flag,sky_z+absorb_z,num_timesteps);
		delete [] rcx;
		delete [] rcy;
		delete [] rcz;

		if (debug)
		{
			earth_model->Write_XZ_Slice_Gnuplot(
					ysrc,
					"/users/tjhc/slices/em_vp.dat",
					"/users/tjhc/slices/em_eps.dat",
					"/users/tjhc/slices/em_dta.dat",
					"/users/tjhc/slices/em_c44c33.dat",
					"/users/tjhc/slices/em_den.dat",
					"/users/tjhc/slices/em_dip.dat",
					"/users/tjhc/slices/em_azm.dat",
					"/users/tjhc/slices/em_Q.dat"
					);

			earth_model->Write_XZ_Slice_Gnuplot(
					ysrc+1,
					"/users/tjhc/slices/em_p1_vp.dat",
					"/users/tjhc/slices/em_p1_eps.dat",
					"/users/tjhc/slices/em_p1_dta.dat",
					"/users/tjhc/slices/em_p1_c44c33.dat",
					"/users/tjhc/slices/em_p1_den.dat",
					"/users/tjhc/slices/em_p1_dip.dat",
					"/users/tjhc/slices/em_p1_azm.dat",
					"/users/tjhc/slices/em_p1_Q.dat"
					);
		}
		printf("Number of output timesteps = %d\n",num_timesteps);

		float fmax = earth_model->Get_FMAX();

		float* stf = new float[2000];
		int tsrc = 0;
		src(4, dt, fmax, 1, "nil", &tsrc, stf);
		prop->Put_Source_Wavelet(stf,tsrc,xsrc,ysrc,zsrc,zsrcghost);
		//prop->AddToPQ(xsrc,ysrc,zsrc,1.0f,1.0f);

		if (debug)
		{
			char segywaveletfilename[1024];
			sprintf(segywaveletfilename, "%ssource_wavelet_%05d.segy", parmfile->seisname, parmfile->id);
			writeSEGY(
					stf, segywaveletfilename, "GPU variable density acoustic TTI source wavelet", 1, parmfile->id,
					xsrc*dx, ysrc*dy, zsrc*dz, dt, dt, 0.0f, tsrc, 1, 
					xsrc*dx, 1, dx,
					ysrc*dy, 1, dy,
					zsrc*dz, 1, dz
				 );
			printf("Source wavelet written to %s.\n",segywaveletfilename);
		}

		float* xz = (float*)malloc(dimx*dimz*sizeof(float));
		float* xy = (float*)malloc(dimx*dimy*sizeof(float));

		struct timespec before1;
		clock_gettime(CLOCK_REALTIME, &before1);
	
		unsigned long h2d_byte_count = prop->Get_H2D_Byte_Count();
		unsigned long d2h_byte_count = prop->Get_D2H_Byte_Count();
		struct timespec before2;
		clock_gettime(CLOCK_REALTIME, &before2);

		int curr_timestep = 0, done = 0, ocnt = 0;
		while (!done)
		{
			int output_timestep;
			done = prop->Process_Two_Stripes(logLevel, num_timesteps, output_timestep);
			if (output_timestep > curr_timestep)
			{
				if (curr_timestep > 0)
				{
					// finished one iteration, i.e. X timesteps
					struct timespec after2;
					clock_gettime(CLOCK_REALTIME, &after2);
					double elapsed_time2 = (double)after2.tv_sec + (double)after2.tv_nsec * 1e-9 - ((double)before2.tv_sec + (double)before2.tv_nsec * 1e-9);
					unsigned long net_h2d_byte_count = prop->Get_H2D_Byte_Count() - h2d_byte_count;
					unsigned long net_d2h_byte_count = prop->Get_D2H_Byte_Count() - d2h_byte_count;
					double data_rate_h2d_GB_sec = (double)net_h2d_byte_count * 1e-9 / elapsed_time2;
					double data_rate_d2h_GB_sec = (double)net_d2h_byte_count * 1e-9 / elapsed_time2;
					double data_rate_GCell_sec = (double)volSize * (double)prop->Get_Total_Number_of_Timesteps() * 1e-9 / elapsed_time2;
					int j0 = 2 * curr_timestep - output_timestep + 1;
					int j1 = curr_timestep;

					printf("Timesteps %4d to %4d - %.2f s, H2D %.2f GB/s, D2H %.2f GB/s, %.0f MCell/s\n",j0,j1,elapsed_time2,data_rate_h2d_GB_sec,data_rate_d2h_GB_sec,data_rate_GCell_sec*1e3);
		
					// output every 4th iteration
					if (debug)
					{
						++ocnt;
						if (ocnt >= 4)
						{
							ocnt = 0;

							prop->XZ_Cross_Section(ysrc,0,xz);
							char fname[256];
							sprintf(fname,"/users/tjhc/slices/GPU_slice_y=%d_t=%d.dat",ysrc,j1);
							FILE* fp = fopen(fname, "w");
							if (fp != 0L)
							{
								for (int iZ = 0;  iZ < dimz;  ++iZ)
								{
									for (int iX = 0;  iX < dimx;  ++iX)
									{
										fprintf(fp,"%d %d %e\n",iX,iZ,xz[iZ*dimx+iX]);
									}
									fprintf(fp,"\n");
								}
								fclose(fp);
							}

							prop->XY_Cross_Section(zsrc,0,xy);
							sprintf(fname,"/users/tjhc/slices/GPU_slice_z=%d_t=%d.dat",zsrc,j1);
							fp = fopen(fname, "w");
							if (fp != 0L)
							{
								for (int iY = 0;  iY < dimy;  ++iY)
								{
									for (int iX = 0;  iX < dimx;  ++iX)
									{
										fprintf(fp,"%d %d %e\n",iX,iY,xy[iY*dimx+iX]);
									}
									fprintf(fp,"\n");
								}
								fclose(fp);
							}
						}
					}

					h2d_byte_count = prop->Get_H2D_Byte_Count();
					d2h_byte_count = prop->Get_D2H_Byte_Count();
					clock_gettime(CLOCK_REALTIME, &before2);
				}
				curr_timestep = output_timestep;
			}
		}
		prop->CompleteReceiverValues();

		struct timespec after1;
		clock_gettime(CLOCK_REALTIME, &after1);
		double elapsed_time1 = (double)after1.tv_sec + (double)after1.tv_nsec * 1e-9 - ((double)before1.tv_sec + (double)before1.tv_nsec * 1e-9);

		char segyfilename[1024];
		sprintf(segyfilename, "%s%05d.segy", parmfile->seisname, parmfile->id);
		writeSEGY(
			prop->Get_Muxed_Traces(), segyfilename, "GPU variable density acoustic TTI propagation", 1, parmfile->id, 
			xsrc*dx, ysrc*dy, zsrc*dz, dt, dtout, 0.0f, prop->Get_NSAMP(), recnum, 
			parmfile->xrecstart*dx, nrecx, parmfile->xrecstride*dx,
			parmfile->yrecstart*dy, nrecy, parmfile->yrecstride*dy,
			parmfile->zrecstart*dz, nrecz, parmfile->zrecstride*dz
			);

		for (int irecx = parmfile->xrecstart;  irecx <= parmfile->xrecend;  irecx += parmfile->xrecstride) ++nrecx;
		for (int irecy = parmfile->yrecstart;  irecy <= parmfile->yrecend;  irecy += parmfile->yrecstride) ++nrecy;
		for (int irecz = parmfile->zrecstart;  irecz <= parmfile->zrecend;  irecz += parmfile->zrecstride) ++nrecz;

		//prop->Check_Host_Memory();

		printf("%d timesteps took %.2fs - %.0fms on average per timestep.\n",num_timesteps,elapsed_time1,1e3f*elapsed_time1/(double)(num_timesteps));
		fflush(stdout);
		//gets(buf);

		delete prop;
	}
}
