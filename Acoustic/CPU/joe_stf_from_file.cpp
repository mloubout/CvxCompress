void src(float dt, float fmax, int type, int* tsrc, float* stf)
{ /* type==1: first derivative of a Gaussian, with linear extension tapers
     type==2: future specification */
int i, imax, imaxhalf, ntap, nfine, tfine;
float w0, wmax, ts,t,rt, wt, wt2, diff, dtfine;
float frac, *stffine;
FILE *stffile;

wmax = 2.0f*M_PI*fmax;
/* Note:   tsrc = ts/dt = 2/(gam*khmax)*(rt*rt)*(Vmax/Vmin) */

if(type==1)
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

else // type==2  source time function provided
{ stffile = fopen("stf","r");
  fscanf(stffile,"%d %f", &nfine, &dtfine);
  stffine = (float*)malloc((nfine+1)*sizeof(float));
  for(i=0; i<nfine; i++) fscanf(stffile,"%f", &stffine[i]);
  stffine[nfine] = 0.;

  imax = (int)((nfine-1)*dtfine/dt) + 1;
  for(i=0; i<imax; i++)
  { t = i*dt;
    tfine = (int)(t/dtfine);
    frac = t/dtfine - tfine;
    stf[i] = (1.-frac)*stffine[tfine] + frac*stffine[tfine+1];
  }

}

