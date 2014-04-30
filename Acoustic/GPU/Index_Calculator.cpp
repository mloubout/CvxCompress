#include <stdio.h>
#include <stdlib.h>
#include "Index_Calculator.hxx"

Index_Calculator::Index_Calculator(int Target_Platform, int num_fields, int nx, int ny, int nz, int sky_z, int absorb_z, int xh, int yh, int zh)
{
	_Target_Platform = Target_Platform;
	_num_fields = num_fields;
	_nx = nx;
	_ny = ny;
	_nz = nz;
	_sky_z = sky_z;
	_absorb_z = absorb_z;
	_xh = xh;
	_yh = yh;
	_zh = zh;
}

Index_Calculator::~Index_Calculator()
{
}

unsigned long Index_Calculator::Get_Field_Offset()
{
	if (_Target_Platform == 1)
	{
		unsigned long stride_z = 16;
		unsigned long stride_field = stride_z * (_nz + _sky_z + _absorb_z);
		return stride_field;
	}
	else
	{
		return 0;
	}
}

