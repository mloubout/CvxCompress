#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Parmfile_Reader.hxx"

Parmfile_Reader::Parmfile_Reader(
		const char* filename
	       )
{
	_Is_Valid = 0;
	FILE* parmfile = fopen(filename, "r");
        if (parmfile == 0L)
        {
                fprintf(stderr,"Cannot find parmfile\n");
        }
        else
        {
                fscanf(parmfile,"%s %d %d %d %d %d %s %s %d %s %f %s %s %s %d %d %d %s", seisname, &id, &Kernel_Type, &fast_Axis, &med_Axis, &slow_Axis,
                                vpname, epsetaname, &eta_Flag, deltaname, &VsoVp0, dnname, dipdxname, azmdyname, &dipxdipy_Flag, &degrees_Flag, &swap_Flag, Qname);
                fscanf(parmfile,"%d %d %d %d %d", &smoothanis_Flag, &isorad, &absorbz0_Flag, &srcghost_Flag, &recghost_Flag);
                fscanf(parmfile,"%d %f", &spongewidth_x, &spongeendval);
                fscanf(parmfile,"%f %f %f %d %d %d", &dh, &dz, &stretchfacz, &nx, &ny, &nz);
                fscanf(parmfile,"%d %d %d %d %d %d", &sub_xoff, &sub_yoff, &sub_zoff, &sub_nx, &sub_ny, &sub_nz);
                fscanf(parmfile,"%f %f %f %d %s", &srcx, &srcy, &srcz, &sourcetype, stfname);
                fscanf(parmfile,"%d %d %d", &xrecstart, &xrecend, &xrecstride);
                fscanf(parmfile,"%d %d %d", &yrecstart, &yrecend, &yrecstride);
                fscanf(parmfile,"%d %d %d", &zrecstart, &zrecend, &zrecstride);
                fscanf(parmfile,"%f %f %f %f %f %d", &maxtime, &timestartrec, &dtout,
                                &newfmax, &gamfac, &OTflag /* &nthread */ );
                fclose(parmfile);

                spongewidth_y = spongewidth_x;
                spongewidth_z_hi = (int)(((float)spongewidth_x * dh / dz) + 0.5f);
                spongewidth_z_lo = spongewidth_z_hi;
		absorbz0 = absorbz0_Flag ? (int)ceilf(fabsf(srcz/dz)) + 2 + spongewidth_z_lo : 0;

                // compute sponge coefficient

                spongecoeff_x = (1.0f - spongeendval) / (float)((spongewidth_x-1) * (spongewidth_x-1));
                spongecoeff_y = (1.0f - spongeendval) / (float)((spongewidth_y-1) * (spongewidth_y-1));
                spongecoeff_z_lo = (1.0f - spongeendval) / (float)((spongewidth_z_lo-1) * (spongewidth_z_lo-1));
                spongecoeff_z_hi = (1.0f - spongeendval) / (float)((spongewidth_z_hi-1) * (spongewidth_z_hi-1));

                // Specifies how data is ordered in input files.
                if (
                        fast_Axis < 0 || fast_Axis > 2 ||
                        med_Axis  < 0 || med_Axis  > 2 ||
                        slow_Axis < 0 || slow_Axis > 2 ||
                        (fast_Axis + med_Axis + slow_Axis) != 3
                        )
                {
                        printf("Invalid axis ordering : %d %d %d\nAborting!\n",fast_Axis,med_Axis,slow_Axis);
                }
		else
		{
			_Is_Valid = 1;
		}
	}
}

Parmfile_Reader::~Parmfile_Reader()
{
}

