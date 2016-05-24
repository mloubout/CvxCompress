#ifndef SEGY_HEADER_STRUCTS_H
#define SEGY_HEADER_STRUCTS_H

struct reel_id_hdr2_t
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
};

struct trc_id_hdr_t
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
};

struct custom_trc_id_hdr_t
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
};

#endif

