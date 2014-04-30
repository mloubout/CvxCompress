#ifndef PARMFILE_READER_HXX
#define PARMFILE_READER_HXX

class Parmfile_Reader
{
	public:
		Parmfile_Reader(
				const char* filename
			       );
		virtual ~Parmfile_Reader();

		int Is_Valid() {return _Is_Valid;}

		char seisname[256];
		char vpname[256];
		char epsetaname[256];
		char deltaname[256];
		char dnname[256];
		char dipdxname[256];
		char azmdyname[256];
		char Qname[256];
		char stfname[256];
		int id;
		int eta_Flag;
		int dipxdipy_Flag;
		int degrees_Flag;
		int swap_Flag;
		int smoothanis_Flag;
		int isorad;
		int absorbz0_Flag;
		int absorbz0;
		int srcghost_Flag;
		int recghost_Flag;
		int nx, ny, nz;
		int sub_xoff, sub_yoff, sub_zoff;
		int sub_nx, sub_ny, sub_nz;
		float VsoVp0;
		float dh;
		float dz;
		float stretchfacz;
		float srcx, srcy, srcz;
		float maxtime;
		float timestartrec;
		float dtout;
		float newfmax;
		float gamfac;
		int sourcetype;
		int OTflag;
		int xrecstart, xrecend, xrecstride;
		int yrecstart, yrecend, yrecstride;
		int zrecstart, zrecend, zrecstride;
		int fast_Axis, med_Axis, slow_Axis;
		int Kernel_Type; // 0->ISO, 1->VTI, 2->TTI
		int spongewidth_x, spongewidth_y, spongewidth_z_lo, spongewidth_z_hi;
		float spongeendval;
                float spongecoeff_x;
                float spongecoeff_y;
                float spongecoeff_z_lo;
                float spongecoeff_z_hi;

	private:
		int _Is_Valid;
};

#endif
