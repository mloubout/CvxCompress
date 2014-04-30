// Function:  writeSEGY.c

//--------------------------------------------------------------------------
/* writeSEGY.c is a stripped down version of Joe Stefani's code EFDSEGYMULT.

   o Note!:
     Could not get this function to work as C++ with the Intel compiler. 

   o Program converts ascii traces to SEGY:
     - Transposes ascii profile output from W[t][x] to W[x][t]
     - Optionally window data by offset;
     - Optionally swap data bytes and/or header bytes, as required by
       SeisSpace.
     - Output in SEGY format.
   
   o Last Modified:  
     01-25-13  Adapted from /users/knih/NumericalCodes/Ssg2d/Svn/trunk/
               segywrite2d.c

   o Comments:
     To enable 0.01 accuracies in source and receiver locations in 
       SeisSpace, use the X,Y scalars (Z values will be correct), by 
       selecting "Yes" to "Use the coordinate scalar?" in the SEG-Y Input 
       module. If this is not selected, the X,Y locations will be 100x too 
       large.
     To use iline_no and xline_no SEGY header info in SeisSpace, select 
       "Yes" to "Remap SEG-Y header values?" in SEG-Y Input module 
       (/iline_no,,4I,,209/xline_no,,4I,,213/).
*/
//--------------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
//#include "mpi.h"

#define  r2d 57.29578

// Define Structures:
typedef struct
{ float x;  float y; float z; }  LOC;

//..SEGY data types 
char reel_id_hdr1[3200];

struct
{ int jobid;  int lineid;  int reelid;  short ntrc_per_record;  short nauxtrc;
  short dtreel;   short dtfield;   short nsampreel;   short nsampfield;
  short datafmt;  short cmpfold;   short sortcode;    char skip[370];
} reel_id_hdr2;

struct
{ int trcseqno;   int skip0;     int  isrc;       int ichan;
  int skip1;      int cmpbin;    int skip2;       short code;   char skip3[6];
  int offset;     int recelev;
  int elevatsrc;  int srcdepth;    char skip4[16];
  short scalar1;  short scalar2;
  int srcx;       int srcy;        int recx;      int recy;     short lenunit;
  char skip5[24]; short nsamp;     short dtmicro; char skip6[82];
  float cmp_x;    float cmp_y;   int iline_no;    int xline_no;
  float xoff;     float yoff;    float azim;
  char skip8[12];
} trc_id_hdr;    // note: cmp_x starts at byte position 201
  
void local2global(double O1, double O2, double O3,
                 double U1, double U2, double U3,
                 double V1, double V2, double V3,
                 double W1, double W2, double W3,
                 char* orderaxes,
                 int nx, int ny, int nz,
                 float dx, float dy, float dz,
		 float x, float y, float z,
		 double* gx, double* gy, double* gz
		)
{
	double _x = (double)x / ((double)(nx-1) * (double)dx);
	double _y = (double)y / ((double)(ny-1) * (double)dy);
	double _z = (double)z / ((double)(nz-1) * (double)dz);
	if (strncmp(orderaxes,"xyz",3) == 0)
	{
		*gx = _x*U1 + _y*V1 + _z*W1 + O1;
		*gy = _x*U2 + _y*V2 + _z*W2 + O2;
		*gz = _x*U3 + _y*V3 + _z*W3 + O3;
	}
	else if (strncmp(orderaxes,"zxy",3) == 0)
	{
		*gx = _x*V1 + _y*W1 + _z*U1 + O1;
		*gy = _x*V2 + _y*W2 + _z*U2 + O2;
		*gz = _x*V3 + _y*W3 + _z*U3 + O3;
	}
	else if (strncmp(orderaxes,"zyx",3) == 0)
	{
		*gx = _x*W1 + _y*V1 + _z*U1 + O1;
		*gy = _x*W2 + _y*V2 + _z*U2 + O2;
		*gz = _x*W3 + _y*V3 + _z*U3 + O3;
	}

	printf("local2global :: (%.2f,%.2f,%.2f) ==> (%.0f,%.0f,%.2f)\n",x,y,z,*gx,*gy,*gz);
}

double getsign(double v)
{
	if (v < 0.0)
		return -1.0;
	else
		return 1.0;
}

//--------------------------------------------------------------------------

void writesegy_ (float *p, int *iline, int *xline, int *ksho, 
  	         char *filename1, float *srcminx, float *srcminy, 
                 float *srcminz , float *xrec, float *yrec, float *zrec, 
	         int *nrec, int *nsampy, float *dt, 
	         int *irec_offset,
		 double* O1, double* O2, double* O3,
		 double* U1, double* U2, double* U3,
		 double* V1, double* V2, double* V3,
		 double* W1, double* W2, double* W3,
		 char* orderaxes,
		 int* nx, int* ny, int* nz,
		 float* dx, float* dy, float* dz
		) 
{

  // Define Variables & Arrays:
  char profname[100];

  int swapflag, nsrc, isrc, irec;
  int irecstart, irecend, nrecout, srcx, srcy, srcz, recx, recy, recz, 
      recelev;
  double srcx_, srcy_, srcz_, recx_, recy_, recz_;
  int t, trcseqno=1, dtmicro;
  int one=1, one4=1, elevatsrc, srcdepth, offset, isr, ichan, trcseq;

  int nsamp = *nsampy;  // this gets used a lot, so let's assign it
  long itr;              // index to break up p vector into rec. traces

  short one2=1, five2=5, nrec2, dtmicro2, nsamp2, fold2, trc_sortcode2;

  short neg100 = -100;  // preserves 0.01 accuracies in sou & rec locations
  int pos100 = 100;

  //float offminout, offmaxout;    // !! disabled windowing feature
  float dsrc;
  float zsurface;
  float *Ptrace;  
  float **matrix(int,int);

  double cmpx, cmpy, xoff, yoff, zoff, azim;
  int ilineno, xlineno;

  // Define Functions:
  void swap2bytes(short *i2, int n), swap4bytes(int*, int);
  void free_matrix(float**, int);

  LOC *src, *rec;

  // SEGY Output File Stuff:
  FILE *Psegy;
  strcpy(profname,filename1);
  Psegy = fopen(profname,"w");    

  printf("O=(%.0f,%.0f,%.0f) - U=(%.0f,%.0f,%.0f) - V=(%.0f,%.0f,%.0f) - W=(%.0f,%.0f,%.0f)\n",*O1,*O2,*O3,*U1,*U2,*U3,*V1,*V2,*V3,*W1,*W2,*W3);
  printf("nx=%d, ny=%d, nz=%d, dx=%.2f, dy=%.2f, dz=%.2f, orderaxes=%s, strlen(orderaxes)=%d\n",*nx,*ny,*nz,*dx,*dy,*dz,orderaxes,strlen(orderaxes));

  memset((void*)reel_id_hdr1, 0, 3200);
  memset((void*)&reel_id_hdr2, 0,  400);

  /*
  // !! Debugging Print-out Stuff:
  printf("hello from inside segywrite3d.c \n");
  printf("filename1 = %s \n ", filename1);
  printf("srcminx = %e \n ", *srcminx);
  printf("srcminz = %e \n ", *srcminz);
  printf("xrec = %e \n ", *xrec);
  printf("nrec = %d \n ", *nrec);
  printf("zrec = %e \n ", *zrec);
  printf("nsamp = %d \n ", nsamp);
  printf("dt = %e \n ", *dt);
  int it, ir;
  for(ir=0; ir<*nrec; ir++) {
    for(it=0; it<nsamp; it++) {
      itr = it + nsamp*ir;    // index to break up p vector into rec. traces
      printf("p[itr] = %e \n", p[itr]);
    }
  } 
 */

  // Set Some Variables as Constants:
  nsrc = 1;  // hardwired to 1 source
  dsrc = 0.;
  isrc = 0;

  swapflag  = 1;   // swapflag = 1 for Linux (little endian)

  src = (LOC *)malloc((nsrc+4)*sizeof(LOC));
  rec = (LOC *)malloc(*nrec*sizeof(LOC));


  src[isrc].x = *srcminx;
  src[isrc].y = *srcminy;
  src[isrc].z = *srcminz;

  for(irec=0; irec<*nrec; irec++) {
    rec[irec].x = xrec[irec];
    rec[irec].y = yrec[irec];
    rec[irec].z = zrec[irec];
  }

  fold2 = 1;
  nrecout = *nrec;  // !! hardwired 

  // Allocate space for one vector:
  Ptrace = (float *)malloc(nsamp*sizeof(float)); 

  // DEALING WITH ONLY ONE SOURCE: 
  isrc=0;

 /*** WRITE REEL ID HEADER 1 ***/
  strcpy(reel_id_hdr1, "C 1 CLIENT     PRESSURE DATA IN SEGY FORMAT\n");
  fwrite(reel_id_hdr1, 3200, 1, Psegy);

 /*** WRITE REEL ID HEADER 2 ***/
  trc_sortcode2 = 1;  /* as recorded, no sorting */
  dtmicro = (int)(1000000.*(*dt) + 0.5);
  one2 = one;                  
  one4 = one;
  nrec2 = nrecout;     
  dtmicro2 = dtmicro;        
  nsamp2 = nsamp;
  if(swapflag) {
    swap2bytes(&one2, 1);       
    swap2bytes(&five2, 1);       
    swap4bytes(&one4, 1);
    swap2bytes(&nrec2, 1);       
    swap2bytes(&dtmicro2, 1);
    swap2bytes(&nsamp2, 1);     
    swap2bytes(&trc_sortcode2, 1);
    swap2bytes(&fold2, 1);
    swap2bytes(&neg100, 1);
  }
  reel_id_hdr2.jobid = reel_id_hdr2.lineid = reel_id_hdr2.reelid = one4;
  reel_id_hdr2.ntrc_per_record = nrec2;
  reel_id_hdr2.dtreel = reel_id_hdr2.dtfield = dtmicro2;
  reel_id_hdr2.nsampreel = reel_id_hdr2.nsampfield = nsamp2;
  //reel_id_hdr2.datafmt = one2;  // for IBM format
  reel_id_hdr2.datafmt = five2;   // for IEEE format
  reel_id_hdr2.cmpfold = fold2;
  reel_id_hdr2.sortcode = trc_sortcode2;
  fwrite(&reel_id_hdr2, 400, 1, Psegy);

  /* TREAT CURRENT SOURCE */
  irecstart = 0;       // !! hardwired 
  irecend = *nrec-1;   // !! hardwired 

  /* LOOP OVER RECEIVERS */
  // write trace ID headers & traces: 
  for(irec=irecstart; irec<=irecend; irec++) {
    printf("irec = %d / %d\n",irec,irecend);
    fflush(stdout);

    // source stuff:
    zsurface = 0;
    isr = *ksho;  // ffid = isr, i.e., ffid = shot#
    elevatsrc = 0;  // hardwired for marine data   


    // assign & write trace header:
    trcseq = trcseqno++;          
    ichan = irec + 1 + *irec_offset;

    // compute all global coordinates as double's
    printf("SRC local2global\n");
    local2global(
	*O1,*O2,*O3,*U1,*U2,*U3,*V1,*V2,*V3,*W1,*W2,*W3,orderaxes,*nx,*ny,*nz,*dx,*dy,*dz,
	src[isrc].x,src[isrc].y,src[isrc].z,
	&srcx_,&srcy_,&srcz_
	);
    printf("REC local2global\n");
    local2global(
	*O1,*O2,*O3,*U1,*U2,*U3,*V1,*V2,*V3,*W1,*W2,*W3,orderaxes,*nx,*ny,*nz,*dx,*dy,*dz,
        rec[irec].x,rec[irec].y,rec[irec].z,  
        &recx_,&recy_,&recz_
        );

    xoff = recx_ - srcx_;        
    yoff = recy_ - srcy_;
    zoff = -recz_ - srcz_;   // minus sign on recz_ to convert from an elevation to a depth

    cmpx = 0.5*(srcx_ + recx_);  
    cmpy = 0.5*(srcy_ + recy_);

    srcdepth = (int)round(100.0*(srcz_ - zsurface));
    srcx = (int)round(100.0*srcx_);	
    srcy = (int)round(100.0*srcy_);	
    srcz = (int)round(100.0*srcz_);	
    recx = (int)round(100.0*recx_);	
    recy = (int)round(100.0*recy_);	
    recz = (int)round(100.0*recz_);
    recelev = recz;	

    ilineno = *iline;
    xlineno = *xline;

    // signed offset:
    xoff = (int)(xoff > 0.0 ? xoff + 0.5 : xoff - 0.5);
    yoff = (int)(yoff > 0.0 ? yoff + 0.5 : yoff - 0.5);
    offset = sqrt(pow(xoff, 2.) + pow(yoff, 2.));
    if (xoff < 0.0) {offset = -offset;}    // sign on offset comes from xoff


    azim = r2d*atan2f(yoff, xoff);

    
    if(swapflag) {
      swap4bytes(&elevatsrc, 1);     
      swap4bytes(&srcdepth, 1);
      swap4bytes(&srcx, 1);          
      swap4bytes(&srcy, 1);          
      swap4bytes(&isr, 1);
      swap4bytes(&offset, 1);
      swap4bytes(&recelev,1);        
      swap4bytes(&recx, 1);
      swap4bytes(&recy, 1);
      swap4bytes(&ilineno, 1);
      swap4bytes(&xlineno, 1);
      swap4bytes(&trcseq, 1);        
      swap4bytes(&ichan, 1);
      swap4bytes((int*)(&xoff), 1);
      swap4bytes((int*)(&yoff), 1);
      swap4bytes((int*)(&cmpx), 1);
      swap4bytes((int*)(&cmpy), 1);
      swap4bytes((int*)(&azim), 1);
    }

    trc_id_hdr.trcseqno = trcseq;
    trc_id_hdr.isrc = isr;
    trc_id_hdr.ichan = ichan;
    //trc_id_hdr.cmpbin = cmpbin;
    trc_id_hdr.code = one2;
    trc_id_hdr.offset = offset;
    trc_id_hdr.recelev = recelev;
    trc_id_hdr.elevatsrc = elevatsrc;
    trc_id_hdr.srcdepth = srcdepth;
    trc_id_hdr.srcx = srcx;
    trc_id_hdr.srcy = srcy;
    trc_id_hdr.recx = recx;
    trc_id_hdr.recy = recy;
    trc_id_hdr.iline_no = ilineno;
    trc_id_hdr.xline_no = xlineno;
    trc_id_hdr.lenunit = one2;
    trc_id_hdr.nsamp = nsamp2;
    trc_id_hdr.dtmicro = dtmicro2;
    trc_id_hdr.scalar1 = neg100;
    trc_id_hdr.scalar2 = neg100;
    trc_id_hdr.xoff = xoff;
    trc_id_hdr.yoff = yoff;
    trc_id_hdr.cmp_x = cmpx;
    trc_id_hdr.cmp_y = cmpy;
    trc_id_hdr.azim = azim;

    fwrite(&trc_id_hdr, 240, 1, Psegy);

    /* Write Trace Vector */
    for(t=0; t<nsamp; t++) {
      //printf("irec = %d \n", irec);    // !! debugging printf's
      //printf("t = %d \n", t);
      itr = (long)t + (long)nsamp * (long)irec;    // index to break up p vector into rec. traces
      Ptrace[t] = p[itr];  
    } 
    if(swapflag) swap4bytes((int*)(Ptrace), nsamp);  // !!! Recast Ptrace to int !!!
    fwrite(Ptrace, sizeof(float), nsamp, Psegy); 

  }  // end of trace loop 


  fclose(Psegy);
  free((char*)src);    
  free((char*)rec);
  free((char*)Ptrace);

}  // end of writeSEGY.c 

// !! Uncomment-out the line below when using Portland Group 64 bit compiler:
// !!PG!! }  // close of extern "C" bracket 



/***** swap2bytes ba --> ab *****/
void swap2bytes(short *i2, int n)
{
  int i;
  short a,b;
  for (i=0; i<n; i++) {
    a = i2[i] << 8;     b = (i2[i] >> 8) & 255;    i2[i] = a | b; 
  }
}


/***** swap4bytes:  dcba --> abcd *****/
void swap4bytes(int *i4, int n)
{
  int k, i, a, b, c, d, bmask = 16711680, cmask = 65280, dmask = 255;
  for(k=0; k<n; k++) {
    i = i4[k];
    a =  i << 24;          
    b = (i << 8)  & bmask;
    c = (i >> 8) & cmask;  
    d = (i >> 24) & dmask;
    i4[k] = a | b | c | d ;
  }
}


/***** STORAGE MANAGEMENT FUNCTIONS ADAPTED FROM NUMERICAL RECIPES IN C *****/
float **matrix(int nt, int nx)   /* allocate a float matrix ntXnx */
{ 
  int i;   
  float **m;
  m=(float **) malloc((unsigned) nt*sizeof(float*));
  if (!m) {
    fprintf(stderr,"Memory allocation failure in matrix.\n");  exit(0); }

  for(i=0;i<nt;i++) {
    m[i]=(float *) malloc((unsigned) nx*sizeof(float));
    if (!m[i]) {
      fprintf(stderr,"Memory allocation failure in matrix.\n");  exit(0); }
  }
  return m;
}

void free_matrix(float **m, int nt)
/* free float matrix allocated with matrix() */
{ 
  int i;
  for(i=nt-1;i>=0;i--) free((char*) m[i]);
  free((char*) m);
}

