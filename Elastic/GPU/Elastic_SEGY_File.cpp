#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "swapbytes.h"
#include "Elastic_SEGY_File.hxx"
#include "Elastic_SEGY_File_Receiver_Range.hxx"
#include "Elastic_Buffer.hxx"
#include "Elastic_Pipeline.hxx"
#include "Elastic_Propagator.hxx"
#include "Elastic_Gather_Type.hxx"
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
	_gather_type = Common_Shot_Gather;
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
	_interpolation_method = Trilinear;

	_num_user_rcv = 0;
	_h_user_rcv_x = 0L;
	_h_user_rcv_y = 0L;
	_h_user_rcv_z = 0L;
	_h_user_iline = 0L;
	_h_user_xline = 0L;
	_h_user_trcens = 0L;
}

Elastic_SEGY_File::~Elastic_SEGY_File()
{
	delete [] _h_user_rcv_x;
	_h_user_rcv_x = 0L;
	delete [] _h_user_rcv_y;
	_h_user_rcv_y = 0L;
	delete [] _h_user_rcv_z;
	_h_user_rcv_z = 0L;
	delete [] _h_user_iline;
	_h_user_iline = 0L;
	delete [] _h_user_xline;
	_h_user_xline = 0L;
	delete [] _h_user_trcens;
	_h_user_trcens = 0L;

	if (_base_filename != 0L) free(_base_filename);
}

int Elastic_SEGY_File::Get_Selection_Flags()
{
        int flags = 0;
        if (_do_P ) flags |= 1;
        if (_do_Vx) flags |= 2;
        if (_do_Vy) flags |= 4;
        if (_do_Vz) flags |= 8;
        return flags;
}

const char* Elastic_SEGY_File::Get_Full_Path(char* buf, int flag)
{
	const char* trace_type_str;
	if (flag == 1)
		trace_type_str = "P";
	else if (flag == 2)
		trace_type_str = "Vx";
	else if (flag == 4)
		trace_type_str = "Vy";
	else if (flag == 8)
		trace_type_str = "Vz";
	else
		trace_type_str = "??";

	sprintf(buf,"%s_%05d_%s.segy",_base_filename,_fileidx,trace_type_str);
	return buf;
}

void Elastic_SEGY_File::Write_Source_Wavelet_To_SEGY_File(
	double* filtered,
	double* filtered_int,
	double sample_rate,
	int nsamp,
	double srcx,
	double srcy,
        double srcz
	)
{
	char buf[4096];
	sprintf(buf,"%s_%05d_Source_Wavelet.segy",_base_filename,_fileidx);
	float* traces[2];
	traces[0] = new float[nsamp];
	traces[1] = new float[nsamp];
	for (int i = 0;  i < nsamp;  ++i) {traces[0][i] = (float)filtered[i]; traces[1][i] = (float)filtered_int[i];}
	double recx[2], recy[2], recz[2];
	recx[0] = recx[1] = srcx;
	recy[0] = recy[1] = srcy;
	recz[0] = recz[1] = srcz;
	int iline[2], xline[2], trcens[2];
	iline[0] = iline[1] = 1;
	xline[0] = xline[1] = 1;
	trcens[0] = trcens[1] = 1;
	float rec_model_water_depth[2], rec_model_Vp[2], rec_bath_z[2];
	rec_model_water_depth[0] = rec_model_water_depth[1] = 0.0f;
	rec_model_Vp[0] = rec_model_Vp[1] = 0.0f;
	rec_bath_z[0] = rec_bath_z[1] = 0.0f;
	this->Write_SEGY_File((const char*)buf,sample_rate,Common_Shot_Gather,_fileidx,0.0,&traces[0],srcx,srcy,srcz,&recx[0],&recy[0],&recz[0],&iline[0],&xline[0],&trcens[0],0.0f,0.0f,0.0f,&rec_model_water_depth[0],&rec_model_Vp[0],&rec_bath_z[0],2,nsamp);
	delete [] traces[0];
	delete [] traces[1];
}

void Elastic_SEGY_File::Write_SEGY_File(
	float** traces,
	double srcx,
	double srcy,
	double srcz,
	double* recx,
	double* recy,
	double* recz,
	int* iline,
	int* xline,
	int* trcens,
	float src_model_water_depth,
	float src_model_water_Vp,
	float src_bath_z,
	float* rec_model_water_depth,
	float* rec_model_water_Vp,
	float* rec_bath_z,
	int num_traces,
	int nsamp,
	int flag
	)
{
	char filename[4096];
	Get_Full_Path(filename,flag);
	this->Write_SEGY_File(
		filename,_sample_rate,Get_Gather_Type(),_fileidx,_tshift,traces,
		srcx,srcy,srcz,recx,recy,recz,iline,xline,trcens,
		src_model_water_depth,src_model_water_Vp,src_bath_z,rec_model_water_depth,rec_model_water_Vp,rec_bath_z,num_traces,nsamp
		);
}

void Elastic_SEGY_File::Write_SEGY_File(
	const char* filename,
	double sample_rate,
	Elastic_Gather_Type_t gather_type,
	int file_idx,
	double start_time,
	float** traces,
	double srcx,
	double srcy,
	double srcz,
	double* recx,
	double* recy,
	double* recz,
	int* iline,
	int* xline,
	int* trcens,
	float src_model_water_depth,
	float src_model_water_Vp,
	float src_bath_z,
	float* rec_model_water_depth,
	float* rec_model_water_Vp,
	float* rec_bath_z,
	int num_traces,
	int nsamp
	)
{
	const int swapflag = 1;
	const double r2d = 57.295779513082320876798154814105;

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
		// 0-3
                int trcseqno;

		// 4-7
                int skip0;

		// 8-11
                int isrc;

		// 12-15
                int ichan;

		// 16-19
                int skip1;

		// 20-23
                int cmpbin;

		// 24-27
                int trcensemb;

		// 28-29
                short code;

		// 30-35
                char skip2[6];

		// 36-39
                int offset;

		// 40-43
                int recelev;

		// 44-47
                int elevatsrc;

		// 48-51
                int srcdepth;

		// 52-59
		char skip3[8];

		// 60-63
		int srcwaterdepth;

		// 64-67
		int recwaterdepth;

		// 68-69
                short scalar1;

		// 70-71
                short scalar2;

		// 72-75
                int srcx;

		// 76-79
                int srcy;

		// 80-83
                int recx;
	
		// 84-87
                int recy;
	
		// 88-89
                short lenunit;
		
		// 90-97
		char skip4[8];

		// 98-99
		short srcstatic;

		// 100-101
		short recstatic;

		// 102-107
                char skip5[6];

		// 108-109
                short tstartrec;

		// 110-113
                char skip6[4];

		// 114-115
                short nsamp;

		// 116-117
                short dtmicro;

		// 118-119
		char skip7[2];

		// 120-127
		double sou_xd;

		// 128-135
		double sou_yd;

		// 136-143
		double rec_xd;

		// 144-151
		double rec_yd;

		// 152-179
		char skip8[28];

		// 180-183
                int cmp_x;

		// 184-187
                int cmp_y;
	
		// 188-191
                int iline_no;

		// 192-195
                int xline_no;

		// 196-199
                int shot_point;

		// 200-201
		short scalar3;

		// 202-213
		char skip9[12];

		// 214-215
		short scalar4;

		// 216-227
		float src_water_model_depth;
		float src_water_bathymetry_depth;
		float src_water_Vp;
		
		// 228-239
		float rec_water_model_depth;
		float rec_water_bathymetry_depth;
		float rec_water_Vp;
        } trc_id_hdr;
        memset((void*)&trc_id_hdr, 0, sizeof(trc_id_hdr));
	/* cmp_x starts at byte position 201 */

	FILE* fp = fopen(filename, "wb");
	if (fp != 0L)
	{
		/*** FILL REEL ID HEADER 1 ***/
		char hdrstring[4096];
		sprintf(hdrstring,"Variable density elastic tilted orthorhombic seismic wave propagation.\nGPU code v0.9\n");
		strcpy(reel_id_hdr1, hdrstring);
		fwrite((void*)reel_id_hdr1,1,3200,fp);

		/*** FILL REEL ID HEADER 2 ***/
		int dtmicro = (int)(1000000.*sample_rate + 0.5);
		short trc_sortcode2 = 1;  /* as recorded, no sorting */
		short fold2 = 1;
		short one2 = 1;       
		int one4 = 1;
		short five2 = 5;
		short nrec2 = num_traces;     
		short dtmicro2 = dtmicro;
		short nsamp2 = nsamp;

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
		fwrite((void*)&reel_id_hdr2,1,400,fp);

		float* trcbuf = new float[nsamp];
		for (int iTrc = 0;  iTrc < num_traces;  ++iTrc)
		{
			double sx, sy, sz, rx, ry, rz;
			float src_water_model_depth;
			float src_water_bathymetry_depth;
			float src_water_Vp;
			float rec_water_model_depth;
			float rec_water_bathymetry_depth;
			float rec_water_Vp;
			if (gather_type == Common_Receiver_Gather)
			{
				// common receiver gather generated through reciprocity
				// flip source and receiver locations
				sx = recx[iTrc];
				sy = recy[iTrc];
				sz = recz[iTrc];
				rx = srcx;
				ry = srcy;
				rz = srcz;
				rec_water_model_depth = src_model_water_depth;
				rec_water_bathymetry_depth = src_bath_z;
				rec_water_Vp = src_model_water_Vp;
				src_water_model_depth = rec_model_water_depth[iTrc];
				src_water_bathymetry_depth = rec_bath_z[iTrc];
				src_water_Vp = rec_model_water_Vp[iTrc];
			}
			else
			{
				sx = srcx;
				sy = srcy;
				sz = srcz;
				rx = recx[iTrc];
				ry = recy[iTrc];
				rz = recz[iTrc];
				src_water_model_depth = src_model_water_depth;
				src_water_bathymetry_depth = src_bath_z;
				src_water_Vp = src_model_water_Vp;
				rec_water_model_depth = rec_model_water_depth[iTrc];
				rec_water_bathymetry_depth = rec_bath_z[iTrc];
				rec_water_Vp = rec_model_water_Vp[iTrc];
			}
		
			int recwaterdepthbathymetry = (int)(10.*rec_water_bathymetry_depth);
			int srcwaterdepthbathymetry = (int)(10.*src_water_bathymetry_depth);
			int srcwaterVp = (int)(10.*src_water_Vp);
			int recwaterVp = (int)(10.*rec_water_Vp);
			int srcwaterdepth = (int)(100.*src_water_model_depth);
			int recwaterdepth = (int)(100.*rec_water_model_depth);

			/*** FILL SOURCE-RELATED PART OF TRACE HEADER ***/
			int ffid = file_idx;
			int elevatsrc = 0;
			int srcdepth = (int)(100.*sz);
			int xsrc = (int)(100.*sx);
			int ysrc = (int)(100.*sy);
			short tstartrec2 = (int)(start_time*1000. + 0.5);
			short neg10 = -10;
			short neg100 = -100;
			short srcstatic = recwaterdepthbathymetry;
			short recstatic = recwaterVp;
			if(swapflag)
			{
				swap4bytes(&ffid, 1);
				swap4bytes(&elevatsrc, 1);     swap4bytes(&srcdepth, 1);
				swap4bytes(&srcwaterdepth, 1);
				swap4bytes(&recwaterdepth, 1);
				swap2bytes(&srcstatic, 1);
				swap2bytes(&recstatic, 1);
				swap4bytes(&xsrc, 1);          swap4bytes(&ysrc, 1);
				swap2bytes(&tstartrec2, 1);
				swap2bytes(&neg100, 1);
				swap2bytes(&neg10, 1);

				swap4bytes((int*)&src_water_bathymetry_depth, 1);
				swap4bytes((int*)&src_water_model_depth, 1);
				swap4bytes((int*)&src_water_Vp, 1);
				swap4bytes((int*)&rec_water_bathymetry_depth, 1);
				swap4bytes((int*)&rec_water_model_depth, 1);
				swap4bytes((int*)&rec_water_Vp, 1);
			}
			trc_id_hdr.isrc = ffid;
			trc_id_hdr.elevatsrc = elevatsrc; 
			trc_id_hdr.srcdepth = srcdepth;
			trc_id_hdr.srcwaterdepth = srcwaterdepth;
			trc_id_hdr.recwaterdepth = recwaterdepth;
			//trc_id_hdr.srcstatic = srcstatic;
			//trc_id_hdr.recstatic = recstatic;
			trc_id_hdr.src_water_bathymetry_depth = src_water_bathymetry_depth;
			trc_id_hdr.src_water_model_depth = src_water_model_depth;
			trc_id_hdr.src_water_Vp = src_water_Vp;
			trc_id_hdr.rec_water_bathymetry_depth = rec_water_bathymetry_depth;
			trc_id_hdr.rec_water_model_depth = rec_water_model_depth;
			trc_id_hdr.rec_water_Vp = rec_water_Vp;
			trc_id_hdr.srcx = xsrc;           
			trc_id_hdr.srcy = ysrc;
			trc_id_hdr.nsamp = nsamp2;
			trc_id_hdr.tstartrec = tstartrec2;
			trc_id_hdr.dtmicro = dtmicro2;
			trc_id_hdr.scalar1 = neg100;
			trc_id_hdr.scalar2 = neg100;
			trc_id_hdr.scalar3 = neg100;
			trc_id_hdr.scalar4 = neg10;

			int recelev = -(int)(100.*rz);
			if(swapflag) swap4bytes(&recelev, 1);
			trc_id_hdr.recelev = recelev;

			int yrec = (int)(100.*ry);
			int xl = xline[iTrc];
			if(swapflag) { swap4bytes(&yrec, 1); swap4bytes(&xl, 1); }
			trc_id_hdr.recy = yrec;
                        trc_id_hdr.iline_no = xl; /* yes, this is correct */

			int trcseq = iTrc+1;       
			int ichan = iTrc+1;      
			int trce = trcens[iTrc];
			int xrec = (int)(100.*rx);
			double xoff = rx - sx;        
			double yoff = ry - sy;
			double cmpx = 0.5*(sx + rx);  
			double cmpy = 0.5*(sy + ry);
			int il = iline[iTrc];
			int offset = (int)round(sqrt(yoff*yoff + xoff*xoff));
			double azim = r2d*atan2(yoff, xoff);

			int xcmp = (int)(100.*cmpx);
			int ycmp = (int)(100.*cmpy);
			if(swapflag)
			{
				swap4bytes(&trcseq, 1);  swap4bytes(&ichan, 1);  swap4bytes(&trce, 1);
				swap4bytes(&xrec, 1);
				swap4bytes((int*)(&xcmp), 1); swap4bytes((int*)(&ycmp), 1);
				swap4bytes(&il, 1);
				swap4bytes(&offset, 1);

				swap8bytes((long*)&sx, 1);
				swap8bytes((long*)&sy, 1);
				swap8bytes((long*)&rx, 1);
				swap8bytes((long*)&ry, 1);
			}
			trc_id_hdr.sou_xd = sx;
			trc_id_hdr.sou_yd = sy;
			trc_id_hdr.rec_xd = rx;
			trc_id_hdr.rec_yd = ry;

			/* Assign & Write Trace Header */
			trc_id_hdr.trcseqno = trcseq;
			trc_id_hdr.ichan = ichan;
			trc_id_hdr.trcensemb = trce;
			trc_id_hdr.offset = offset;
			trc_id_hdr.recx = xrec;
			trc_id_hdr.cmp_x = xcmp;
			trc_id_hdr.cmp_y = ycmp;
			trc_id_hdr.xline_no = il; /* yes, this is correct */
			//trc_id_hdr.xoff = xoff;
			//trc_id_hdr.yoff = yoff;
			//trc_id_hdr.azim = azim;
			fwrite((void*)&trc_id_hdr,1,240,fp);

			if (swapflag)
			{
				memcpy((void*)trcbuf, (void*)(traces[iTrc]), nsamp*sizeof(float));
				swap4bytes((int*)trcbuf,nsamp);
				fwrite((void*)trcbuf,4,nsamp,fp);
			}
			else
			{
				fwrite((void*)traces[iTrc],4,nsamp,fp);
			}
		}
		delete [] trcbuf;

		fclose(fp);
	}
}

/*Add an array of receivers that may not have a fixed range. Range is not set here!
*/
void Elastic_SEGY_File::Add_Receiver_Array(
		int nrec,
		double* rcv_x,
		double* rcv_y,
		double* rcv_z,
		int* iline,
		int* xline,
		int* trcens
)
{

	_h_user_rcv_x = new double[nrec]; 
	_h_user_rcv_y = new double[nrec]; 
	_h_user_rcv_z = new double[nrec];
	_h_user_iline = new int[nrec];
	_h_user_xline = new int[nrec];
	_h_user_trcens = new int[nrec];
	_num_user_rcv = nrec;

	for (int i=0;i<_num_user_rcv;i++) {
		_h_user_rcv_x[i]=rcv_x[i];
		_h_user_rcv_y[i]=rcv_y[i];
		_h_user_rcv_z[i]=rcv_z[i];
		_h_user_iline[i]=iline[i];
		_h_user_xline[i]=xline[i];
		_h_user_trcens[i]=trcens[i];
	}

}
void Elastic_SEGY_File::printRec() {
	for (int i=0;i<_num_user_rcv;i++) {
		printf("printRec::%d %12.3f\t%12.3f\t%12.3f\n",i+1,_h_user_rcv_x[i],_h_user_rcv_y[i],_h_user_rcv_z[i]);
	}
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

int Elastic_SEGY_File::Compute_Receiver_Locations(
		double*& rcv_x,
		double*& rcv_y,
		double*& rcv_z
		)
{
	double *rx, *ry, *rz;
	int *il, *xl, *te;
	int num = Compute_Receiver_Locations_NO_COPY(rx,ry,rz,il,xl,te);
	if (num > 0)
	{
		rcv_x = new double[num];
		rcv_y = new double[num];
		rcv_z = new double[num];
		for (int i = 0;  i < num;  ++i)
		{
			rcv_x[i] = rx[i];
			rcv_y[i] = ry[i];
			rcv_z[i] = rz[i];
		}
	}
	else
	{
		rcv_x = 0L;
		rcv_y = 0L;
		rcv_z = 0L;
	}
	return num;
}

int Elastic_SEGY_File::Compute_Receiver_Locations_NO_COPY(
		double*& rcv_x,
		double*& rcv_y,
		double*& rcv_z
		)
{
	int *iline, *xline, *trcens;
	int num_rx = Compute_Receiver_Locations_NO_COPY(rcv_x,rcv_y,rcv_z,iline,xline,trcens);
	return num_rx;
}

int Elastic_SEGY_File::Compute_Receiver_Locations(
		double*& rcv_x,
		double*& rcv_y,
		double*& rcv_z,
		int*& iline,
		int*& xline,
		int*& trcens
		)
{
	double *rx, *ry, *rz;
	int *il, *xl, *te;
	int num = Compute_Receiver_Locations_NO_COPY(rx,ry,rz,il,xl,te);
	if (num > 0)
	{
		rcv_x = new double[num];
		rcv_y = new double[num];
		rcv_z = new double[num];
		iline = new int[num];
		xline = new int[num];
		trcens = new int[num];
		for (int i = 0;  i < num;  ++i)
		{
			rcv_x[i] = rx[i];
			rcv_y[i] = ry[i];
			rcv_z[i] = rz[i];
			iline[i] = il[i];
			xline[i] = xl[i];
			trcens[i] = te[i];
		}
	}
	else
	{
		rcv_x = 0L;
		rcv_y = 0L;
		rcv_z = 0L;
		iline = 0L;
		xline = 0L;
		trcens = 0L;
	}
	return num;
}

int Elastic_SEGY_File::Compute_Receiver_Locations_NO_COPY(
		double*& rcv_x,
		double*& rcv_y,
		double*& rcv_z,
		int*& iline,
		int*& xline,
		int*& trcens
		)
{
	int num = 0;
	rcv_x = 0L;
	rcv_y = 0L;
	rcv_z = 0L;
	iline = 0L;
	xline = 0L;
	trcens = 0L;

	if (_h_user_rcv_x == 0L) { //Array of receivers has not been specified by user, create arrays from range in parmfile

		for (int i = 0;  i < _num_rcv_ranges;  ++i)
		{
			double *x,*y,*z;
			int *il,*xl,*trce;
			int nn = _rcv_ranges[i]->Compute_Receiver_Locations(x,y,z,il,xl,trce);
			if (nn > 0)
			{
				if (num == 0)
				{
					rcv_x = x;
					rcv_y = y;
					rcv_z = z;
					iline = il;
					xline = xl;
					trcens = trce;
					num = nn;
				}
				else
				{
					double* tmp = new double[num+nn];
					for (int j = 0;  j < num;  ++j) tmp[j] = rcv_x[j];
					for (int j = 0;  j < nn;  ++j) tmp[j+num] = x[j];
					delete [] rcv_x;
					delete [] x;
					rcv_x = tmp;

					tmp = new double[num+nn];
					for (int j = 0;  j < num;  ++j) tmp[j] = rcv_y[j];
					for (int j = 0;  j < nn;  ++j) tmp[j+num] = y[j];
					delete [] rcv_y;
					delete [] y;
					rcv_y = tmp;

					tmp = new double[num+nn];
					for (int j = 0;  j < num;  ++j) tmp[j] = rcv_z[j];
					for (int j = 0;  j < nn;  ++j) tmp[j+num] = z[j];
					delete [] rcv_z;
					delete [] z;
					rcv_z = tmp;

					int* itmp = new int[num+nn];
					for (int j = 0;  j < num;  ++j) itmp[j] = iline[j];
					for (int j = 0;  j < nn;  ++j) itmp[j+num] = il[j];
					delete [] iline;
					delete [] il;
					iline = itmp;

					itmp = new int[num+nn];
					for (int j = 0;  j < num;  ++j) itmp[j] = xline[j];
					for (int j = 0;  j < nn;  ++j) itmp[j+num] = xl[j];
					delete [] xline;
					delete [] xl;
					xline = itmp;

					itmp = new int[num+nn];
					for (int j = 0;  j < num;  ++j) itmp[j] = trcens[j];
					for (int j = 0;  j < nn;  ++j) itmp[j+num] = trce[j];
					delete [] trcens;
					delete [] trce;
					trcens = itmp;

					num = num + nn;
				}
			}
		}

		_h_user_rcv_x = rcv_x;
		_h_user_rcv_y = rcv_y;
		_h_user_rcv_z = rcv_z;
		_h_user_iline = iline;
		_h_user_xline = xline;
		_h_user_trcens = trcens;
		_num_user_rcv = num;
	} 
	else { //Array of receivers has been specified by the user
		rcv_x = _h_user_rcv_x;		
		rcv_y = _h_user_rcv_y;		
		rcv_z = _h_user_rcv_z;
		iline = _h_user_iline;
		xline = _h_user_xline;
		trcens = _h_user_trcens;
		num = _num_user_rcv;		
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

