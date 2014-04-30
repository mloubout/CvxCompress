#include <stdio.h>
#include <stdlib.h>
#include "UVW_XYZ_Calculator.hxx"

UVW_XYZ_Calculator::UVW_XYZ_Calculator(int fast_Axis, int med_Axis, int slow_Axis)
{
	_fast_Axis = fast_Axis;
	_med_Axis = med_Axis;
	_slow_Axis = slow_Axis;
}

UVW_XYZ_Calculator::~UVW_XYZ_Calculator()
{
}

int UVW_XYZ_Calculator::Is_Valid()
{
	if (
                                _fast_Axis >= 0 && _fast_Axis <= 2 &&
                                _med_Axis  >= 0 && _med_Axis  <= 2 &&
                                _slow_Axis >= 0 && _slow_Axis <= 2 &&
                                ((1 << _fast_Axis) | (1 << _med_Axis) || (1 << _slow_Axis)) == 7
                   )
                {
                        return 1;
                }
                else
                {
			return 0;
                }
}

void UVW_XYZ_Calculator::Compute_XYZ_From_UVW(int u, int v, int w, int& x, int& y, int& z)
{
	switch (_fast_Axis)
	{
		case 0:
			x = u;
			break;
		case 1:
			y = u;
			break;
		case 2:
			z = u;
			break;
		default:
			fprintf(stderr, "ERROR! Bad value for U!\n");
			exit(-1);
	}
	switch (_med_Axis)
	{
		case 0:
			x = v;
			break;
		case 1:
			y = v;
			break;
		case 2:
			z = v;
			break;
		default:
			fprintf(stderr, "ERROR! Bad value for V!\n");
			exit(-1);
	}
	switch (_slow_Axis)
	{
		case 0:
			x = w;
			break;
		case 1:
			y = w;
			break;
		case 2:
			z = w;
			break;
		default:
			fprintf(stderr, "ERROR! Bad value for W!\n");
			exit(-1);
	}
}

void UVW_XYZ_Calculator::Compute_UVW_From_XYZ(int x, int y, int z, int& u, int& v, int& w)
{
	switch (_fast_Axis)
	{
		case 0:
			u = x;
			break;
		case 1:
			u = y;
			break;
		case 2:
			u = z;
			break;
		default:
			fprintf(stderr, "ERROR! Bad value for U!\n");
			exit(-1);
	}
	switch (_med_Axis)
	{
		case 0:
			v = x;
			break;
		case 1:
			v = y;
			break;
		case 2:
			v = z;
			break;
		default:
			fprintf(stderr, "ERROR! Bad value for V!\n");
			exit(-1);
	}
	switch (_slow_Axis)
	{
		case 0:
			w = x;
			break;
		case 1:
			w = y;
			break;
		case 2:
			w = z;
			break;
		default:
			fprintf(stderr, "ERROR! Bad value for W!\n");
			exit(-1);
	}
}

