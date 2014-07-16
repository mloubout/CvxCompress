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
	_interpolation_method = Trilinear;

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

/***** swap2bytes ba --> ab *****/
void Elastic_SEGY_File::swap2bytes(short *i2, int n)
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
void Elastic_SEGY_File::swap4bytes(int *i4, int n)
{
	int k, i, a, b, c, d, bmask = 16711680, cmask = 65280, dmask = 255;
	for(k=0; k<n; k++)
	{ i = i4[k];
		a =  i << 24;          b = (i << 8)  & bmask;
		c = (i >> 8) & cmask;  d = (i >> 24) & dmask;
		i4[k] = a | b | c | d ;
	}
}

char* Elastic_SEGY_File::Get_Full_Path(char* buf, int flag)
{
	char* trace_type_str;
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
	int num_traces,
	int nsamp,
	int flag
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

	char filename[4096];
	Get_Full_Path(filename,flag);
	FILE* fp = fopen(filename, "wb");
	if (fp != 0L)
	{
		/*** FILL REEL ID HEADER 1 ***/
		char hdrstring[4096];
		sprintf(hdrstring,"Variable density elastic tilted orthorhombic seismic wave propagation.\nGPU code v0.9\n");
		strcpy(reel_id_hdr1, hdrstring);
		fwrite((void*)reel_id_hdr1,1,3200,fp);

		/*** FILL REEL ID HEADER 2 ***/
		int dtmicro = (int)(1000000.*_sample_rate + 0.5);
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

		/*** FILL SOURCE-RELATED PART OF TRACE HEADER ***/
		int ffid = _fileidx;
		int elevatsrc = 0;
		int srcdepth = (int)(100.*srcz);
		int xsrc = (int)(100.*srcx);
		int ysrc = (int)(100.*srcy);
		short tstartrec2 = (int)(_tshift*1000. + 0.5);
		short neg100 = -100;
		if(swapflag)
		{
			swap4bytes(&ffid, 1);
			swap4bytes(&elevatsrc, 1);     swap4bytes(&srcdepth, 1);
			swap4bytes(&xsrc, 1);          swap4bytes(&ysrc, 1);
			swap2bytes(&tstartrec2, 1);
			swap2bytes(&neg100, 1);
		}
		trc_id_hdr.isrc = ffid;
		trc_id_hdr.elevatsrc = elevatsrc; 
		trc_id_hdr.srcdepth = srcdepth;
		trc_id_hdr.srcx = xsrc;           
		trc_id_hdr.srcy = ysrc;
		trc_id_hdr.nsamp = nsamp2;
		trc_id_hdr.tstartrec = tstartrec2;
		trc_id_hdr.dtmicro = dtmicro2;
		trc_id_hdr.scalar1 = neg100;
		trc_id_hdr.scalar2 = neg100;

		float* trcbuf = new float[nsamp];
		for (int iTrc = 0;  iTrc < num_traces;  ++iTrc)
		{
			int recelev = -(int)(100.*recz[iTrc]);
			if(swapflag) swap4bytes(&recelev, 1);
			trc_id_hdr.recelev = recelev;

			int yrec = (int)(100.*recy[iTrc]);
			int xl = xline[iTrc];
			if(swapflag) { swap4bytes(&yrec, 1); swap4bytes(&xl, 1); }
			trc_id_hdr.recy = yrec;
                        trc_id_hdr.iline_no = xl; /* yes, this is correct */

			int trcseq = iTrc+1;       
			int ichan = iTrc+1;      
			int trce = trcens[iTrc];
			int xrec = (int)(100.*recx[iTrc]);
			double xoff = recx[iTrc] - srcx;        
			double yoff = recy[iTrc] - srcy;
			double cmpx = 0.5*(srcx + recx[iTrc]);  
			double cmpy = 0.5*(srcy + recy[iTrc]);
			int il = iline[iTrc];
			int offset = (int)round(sqrt(yoff*yoff + xoff*xoff));
			double azim = r2d*atan2(yoff, xoff);

			if(swapflag)
			{
				swap4bytes(&trcseq, 1);  swap4bytes(&ichan, 1);  swap4bytes(&trce, 1);
				swap4bytes(&xrec, 1);
				swap4bytes((int*)(&cmpx), 1); swap4bytes((int*)(&cmpy), 1);
				swap4bytes(&il, 1);
				swap4bytes((int*)(&xoff), 1); swap4bytes((int*)(&yoff), 1);
				swap4bytes(&offset, 1);       swap4bytes((int*)(&azim), 1);
			}

			/* Assign & Write Trace Header */
			trc_id_hdr.trcseqno = trcseq;
			trc_id_hdr.ichan = ichan;
			trc_id_hdr.trcensemb = trce;
			trc_id_hdr.offset = offset;
			trc_id_hdr.recx = xrec;
			trc_id_hdr.cmp_x = cmpx;
			trc_id_hdr.cmp_y = cmpy;
			trc_id_hdr.xline_no = il; /* yes, this is correct */
			trc_id_hdr.xoff = xoff;
			trc_id_hdr.yoff = yoff;
			trc_id_hdr.azim = azim;
			fwrite((void*)&trc_id_hdr,1,240,fp);

			if (swapflag)
			{
				for (int i = 0;  i < nsamp;  ++i)
				{
					float v = traces[iTrc][i];
					swap4bytes((int*)&v, 1);
					trcbuf[i] = v;
				}
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

	printf("***num_user_rcv= %d\n",_num_user_rcv);
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
	int *iline, *xline, *trcens;
	int num_rx = Compute_Receiver_Locations(rcv_x,rcv_y,rcv_z,iline,xline,trcens);
	delete [] iline;
	delete [] xline;
	delete [] trcens;
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
	int num = 0;
	rcv_x = 0L;
	rcv_y = 0L;
	rcv_z = 0L;
	iline = 0L;
	xline = 0L;
	trcens = 0L;

	if (_h_user_rcv_x == 0L) { //Array of receivers has not been specified by user, create arrays from range in parmfile
		printf("Receiver arrays created by fd tool from ranges in param file\n");

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
	} else { //Array of receivers has been specified by the user
		printf("Receiver Arrays provided by user\n");
		rcv_x = new double[_num_user_rcv];
		rcv_y = new double[_num_user_rcv];
		rcv_z = new double[_num_user_rcv];
		iline = new int[_num_user_rcv];
		xline = new int[_num_user_rcv];
		trcens = new int[_num_user_rcv];
		for (int i=0;i<_num_user_rcv;i++) {
			rcv_x[i]=_h_user_rcv_x[i];
			rcv_y[i]=_h_user_rcv_y[i];
			rcv_z[i]=_h_user_rcv_z[i];
			iline[i]=_h_user_iline[i];
			xline[i]=_h_user_xline[i];
			trcens[i]=_h_user_trcens[i];
		}
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

