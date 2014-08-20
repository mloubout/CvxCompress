#include "Elastic_SEGY_File_Receiver_Range.hxx"

Elastic_SEGY_File_Receiver_Range::Elastic_SEGY_File_Receiver_Range(int range_idx)
{
	_range_idx = range_idx;
	_has_x = false;
	_has_y = false;
	_has_z = false;
}

Elastic_SEGY_File_Receiver_Range::~Elastic_SEGY_File_Receiver_Range()
{
}

void Elastic_SEGY_File_Receiver_Range::Add_X(double start, double end, double interval)
{
	_has_x = true;
	_x_start = start;
	_x_end = end;
	_x_interval = interval;
}

void Elastic_SEGY_File_Receiver_Range::Add_Y(double start, double end, double interval)
{
	_has_y = true;
	_y_start = start;
	_y_end = end;
	_y_interval = interval;
}

void Elastic_SEGY_File_Receiver_Range::Add_Z(double start, double end, double interval)
{
	_has_z = true;
	_z_start = start;
	_z_end = end;
	_z_interval = interval;
}

int Elastic_SEGY_File_Receiver_Range::Compute_Receiver_Locations(
		double*& rcv_x,
		double*& rcv_y,
		double*& rcv_z,
		int*& iline,
		int*& xline,
		int*& trcens
		)
{
	if (Is_Valid())
	{
		// add 1e-5 of interval to prevent round off error from eliminating last point
		int nx = (int)((_x_end + _x_interval * 1e-5 - _x_start) / _x_interval) + 1;

		// add 1e-5 of interval to prevent round off error from eliminating last point
		int ny = (int)((_y_end + _y_interval * 1e-5 - _y_start) / _y_interval) + 1;

		// add 1e-5 of interval to prevent round off error from eliminating last point
		int nz = (int)((_z_end + _z_interval * 1e-5 - _z_start) / _z_interval) + 1;

		rcv_x = new double[nx*ny*nz];
		rcv_y = new double[nx*ny*nz];
		rcv_z = new double[nx*ny*nz];
		iline = new int[nx*ny*nz];
		xline = new int[nx*ny*nz];
		trcens = new int[nx*ny*nz];
#pragma omp paralle for
		for (int iz = 0;  iz < nz;  ++iz)
		{
			for (int iy = 0;  iy < ny;  ++iy)
			{
				for (int ix = 0;  ix < nx;  ++ix)
				{
					int idx = ix + iy * nx + iz * nx * ny;
					rcv_x[idx] = _x_start + (double)ix * _x_interval;
					rcv_y[idx] = _y_start + (double)iy * _y_interval;
					rcv_z[idx] = _z_start + (double)iz * _z_interval;
					// TO-DO :: Add proper offsets to these inline, xline numbers
					xline[idx] = ix;
					iline[idx] = iy;
					trcens[idx] = ix+1;
				}
			}
		}
		return nx*ny*nz;
	}
	rcv_x = 0L;
	rcv_y = 0L;
	rcv_z = 0L;
	iline = 0L;
	xline = 0L;
	trcens = 0L;
	return 0;
}

