#ifndef CVX_ESDRD_MI_TMJ_EXTRACT_RECEIVER_VALUES_HXX
#define CVX_ESDRD_MI_TMJ_EXTRACT_RECEIVER_VALUES_HXX

void 
Host_Extract_Receiver_Values(
	cudaStream_t stream,
	int* is_vp,
        float** cmp,
        int* x0,
        int* y0,
        int* nx,
        int* ny,
        int* nz,
        float** rcv_binned,
        int* num_rx,
        float** res,
	int num_kernels
        );

#endif

