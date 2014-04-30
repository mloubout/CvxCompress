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
                                ((1 << _fast_Axis) | (1 << _med_Axis) | (1 << _slow_Axis)) == 7
                   )
                {
                        return 1;
                }
                else
                {
			return 0;
                }
}

