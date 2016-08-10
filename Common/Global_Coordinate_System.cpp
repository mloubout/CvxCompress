#include <math.h>
#include <stdio.h>
#include <string.h>
#include <Global_Coordinate_System.hxx>

Global_Coordinate_System::Global_Coordinate_System(const char* path)
{
	_NU = _NV = _NW = 0;
	_DU = _DV = _DW = 0.0;
	_Read(path);
	_transpose = 5;  // default is xyz
	if (_DU == 0.0 && _NU != 0)
	{
		_DU = sqrt(_U0*_U0+_U1*_U1+_U2*_U2) * (_LUMAX-_LUMIN) / (double)(_NU - 1);
		_DV = sqrt(_V0*_V0+_V1*_V1+_V2*_V2) * (_LVMAX-_LVMIN) / (double)(_NV - 1);
		_DW = sqrt(_W0*_W0+_W1*_W1+_W2*_W2) * (_LWMAX-_LWMIN) / (double)(_NW - 1);
	}
	else if (_NU == 0 && _DU != 0.0)
	{
		_NU = (int)round(sqrt(_U0*_U0+_U1*_U1+_U2*_U2) * (_LUMAX-_LUMIN) / _DU) + 1;
		_NV = (int)round(sqrt(_V0*_V0+_V1*_V1+_V2*_V2) * (_LVMAX-_LVMIN) / _DV) + 1;
		_NW = (int)round(sqrt(_W0*_W0+_W1*_W1+_W2*_W2) * (_LWMAX-_LWMIN) / _DW) + 1;
	}
	printf("_DU=%.7e, _DV=%.7f, _DW=%.7f\n",_DU,_DV,_DW);

	double lenU = sqrt(_U0*_U0 + _U1*_U1 + _U2*_U2);
	double lenV = sqrt(_V0*_V0 + _V1*_V1 + _V2*_V2);
	double lenW = sqrt(_W0*_W0 + _W1*_W1 + _W2*_W2);
	double u0 = _U0 / lenU;
	double u1 = _U1 / lenU;
	double u2 = _U2 / lenU;
	double v0 = _V0 / lenV;
	double v1 = _V1 / lenV;
	double v2 = _V2 / lenV;
	double w0 = _W0 / lenW;
	double w1 = _W1 / lenW;
	double w2 = _W2 / lenW;

	double UonV = fabs(u0*v0 + u1*v1 + u2*v2);
	double UonW = fabs(u0*w0 + u1*w1 + u2*w2);
	double VonW = fabs(v0*w0 + v1*w1 + v1*w2);
	if (UonV > 1e-7 || UonW > 1e-7 || VonW > 1e-7)
	{
		printf("WARNING! This is not an orthogonal coordinate system.\n");
		printf("         Mapping between global and local coordinates will introduce errors.\n");
		double maxU = UonV > UonW ? UonV : UonW;
		double maxV = UonV > VonW ? UonV : VonW;
		double maxW = UonW > VonW ? UonW : VonW;
		if (maxU > 1e-7) printf("         Error per 100km for U axis is approximately %.2lfm.\n",maxU*1e5);
		if (maxV > 1e-7) printf("         Error per 100km for V axis is approximately %.2lfm.\n",maxV*1e5);
		if (maxW > 1e-7) printf("         Error per 100km for W axis is approximately %.2lfm.\n",maxW*1e5);
	}

	// override these so that meaning of local coordinate is what users expect.
	//_LUMIN = 0.0;
	//_LUMAX = _DU * (double)(_NU - 1);
	//_LVMIN = 0.0;
	//_LVMAX = _DV * (double)(_NV - 1);
	//_LWMIN = 0.0;
	//_LWMAX = _DW * (double)(_NW - 1);
}

Global_Coordinate_System::~Global_Coordinate_System()
{
}

int Global_Coordinate_System::Get_NU()
{
	return _NU;
}

int Global_Coordinate_System::Get_NV()
{
	return _NV;
}

int Global_Coordinate_System::Get_NW()
{
	return _NW;
}

double Global_Coordinate_System::Get_DU()
{
	return _DU;
}

double Global_Coordinate_System::Get_DV()
{
	return _DV;
}

double Global_Coordinate_System::Get_DW()
{
	return _DW;
}

int Global_Coordinate_System::Get_NX()
{
	int nx,ny,nz;
	Convert_Local_Index_To_Transposed_Index(_NU,_NV,_NW,nx,ny,nz);
	return nx;
}

int Global_Coordinate_System::Get_NY()
{
	int nx,ny,nz;
	Convert_Local_Index_To_Transposed_Index(_NU,_NV,_NW,nx,ny,nz);
	return ny;
}

int Global_Coordinate_System::Get_NZ()
{
	int nx,ny,nz;
	Convert_Local_Index_To_Transposed_Index(_NU,_NV,_NW,nx,ny,nz);
	return nz;
}

double Global_Coordinate_System::Get_DX()
{
	double dx,dy,dz;
	Convert_Fractional_Local_Index_To_Transposed_Fractional_Index(_DU,_DV,_DW,dx,dy,dz);
	return dx;
}

double Global_Coordinate_System::Get_DY()
{
	double dx,dy,dz;
	Convert_Fractional_Local_Index_To_Transposed_Fractional_Index(_DU,_DV,_DW,dx,dy,dz);
	return dy;
}

double Global_Coordinate_System::Get_DZ()
{
	double dx,dy,dz;
	Convert_Fractional_Local_Index_To_Transposed_Fractional_Index(_DU,_DV,_DW,dx,dy,dz);
	return dz;
}

void Global_Coordinate_System::Convert_Global_To_Normalized_Local(
		double g0,
		double g1,
		double g2,
		double &nlu,
		double &nlv,
		double &nlw
		)
{
        double lg0 = g0 - _O0;
        double lg1 = g1 - _O1;
        double lg2 = g2 - _O2;

        double lenU = _U0*_U0 + _U1*_U1 + _U2*_U2;
        double lenV = _V0*_V0 + _V1*_V1 + _V2*_V2;
        double lenW = _W0*_W0 + _W1*_W1 + _W2*_W2;

        nlu = (lg0*_U0 + lg1*_U1 + lg2*_U2) / lenU;
        nlv = (lg0*_V0 + lg1*_V1 + lg2*_V2) / lenV;
        nlw = (lg0*_W0 + lg1*_W1 + lg2*_W2) / lenW;
}

void Global_Coordinate_System::Convert_Normalized_Local_To_Global(
		double nlu,
		double nlv,
		double nlw,
		double &g0,
		double &g1,
		double &g2
		)
{
	g0 = _O0 + _U0 * nlu + _V0 * nlv + _W0 * nlw;
	g1 = _O1 + _U1 * nlu + _V1 * nlv + _W1 * nlw;
	g2 = _O2 + _U2 * nlu + _V2 * nlv + _W2 * nlw;
}

void Global_Coordinate_System::Convert_Local_To_Normalized_Local(
		double lu,
		double lv,
		double lw,
		double &nlu,
		double &nlv,
		double &nlw
		)
{
	nlu = lu / _DU;
	nlv = lv / _DV;
	nlw = lw / _DW;
}

void Global_Coordinate_System::Convert_Normalized_Local_To_Local(
		double nlu,
		double nlv,
		double nlw,
		double &lu,
		double &lv,
		double &lw
		)
{
	lu = nlu * _DU;
	lv = nlv * _DV;
	lw = nlw * _DW;
}

void Global_Coordinate_System::Convert_Normalized_Local_To_Fractional_Local_Index(
		double nlu,
		double nlv,
		double nlw,
		double &filu,
		double &filv,
		double &filw
		)
{
	filu = (nlu - _LUMIN) * (double)(_NU - 1) / (_LUMAX - _LUMIN);
	filv = (nlv - _LVMIN) * (double)(_NV - 1) / (_LVMAX - _LVMIN);
	filw = (nlw - _LWMIN) * (double)(_NW - 1) / (_LWMAX - _LWMIN);
}

void Global_Coordinate_System::Convert_Fractional_Local_Index_To_Normalized_Local(
		double filu,
		double filv,
		double filw,
		double &nlu,
		double &nlv,
		double &nlw
		)
{
	nlu = filu * (_LUMAX - _LUMIN) / (double)(_NU - 1) + _LUMIN;
	nlv = filv * (_LVMAX - _LVMIN) / (double)(_NV - 1) + _LVMIN;
	nlw = filw * (_LWMAX - _LWMIN) / (double)(_NW - 1) + _LWMIN;
}

void Global_Coordinate_System::Convert_Fractional_Local_Index_To_Local_Index(
		double filu,
		double filv,
		double filw,
		int& ilu,
		int& ilv,
		int& ilw
		)
{
	ilu = (int)round(filu);
        ilv = (int)round(filv);
        ilw = (int)round(filw);
}

void Global_Coordinate_System::Convert_Local_Index_To_Fractional_Local_Index(
		int ilu,
		int ilv,
		int ilw,
		double &filu,
		double &filv,
		double &filw
		)
{
	filu = (double)ilu;
	filv = (double)ilv;
	filw = (double)ilw;
}

void Global_Coordinate_System::Convert_Normalized_Local_To_Local_Index(
		double nlu,
		double nlv,
		double nlw,
		int &ilu,
		int &ilv,
		int &ilw
		)
{
	double filu, filv, filw;
	Convert_Normalized_Local_To_Fractional_Local_Index(nlu,nlv,nlw,filu,filv,filw);
	Convert_Fractional_Local_Index_To_Local_Index(filu,filv,filw,ilu,ilv,ilw);
}

void Global_Coordinate_System::Convert_Local_Index_To_Normalized_Local(
		int ilu,
		int ilv,
		int ilw,
		double &nlu,
		double &nlv,
		double &nlw
		)
{
	double filu,filv,filw;
	Convert_Local_Index_To_Fractional_Local_Index(ilu,ilv,ilw,filu,filv,filw);
	Convert_Fractional_Local_Index_To_Normalized_Local(filu,filv,filw,nlu,nlv,nlw);
}

void Global_Coordinate_System::Convert_Global_To_Local(
		double g0, 
		double g1, 
		double g2,
		double &lu, 
		double &lv,
		double &lw
		)
{
	double nlu, nlv, nlw;
	Convert_Global_To_Normalized_Local(g0,g1,g2,nlu,nlv,nlw);
	Convert_Normalized_Local_To_Local(nlu,nlv,nlw,lu,lv,lw);
}

void Global_Coordinate_System::Convert_Local_To_Global(
		double lu,
		double lv,
		double lw,
		double &g0,
		double &g1,
		double &g2
		)
{
	double nlu, nlv, nlw;
	Convert_Local_To_Normalized_Local(lu,lv,lw,nlu,nlv,nlw);
	Convert_Normalized_Local_To_Global(nlu,nlv,nlw,g0,g1,g2);
}

void Global_Coordinate_System::Convert_Global_To_Local_Index(
		double g0,
		double g1,
		double g2,
		int &ilu,
		int &ilv,
		int &ilw
		)
{
	double nlu, nlv, nlw;
	Convert_Global_To_Normalized_Local(g0,g1,g2,nlu,nlv,nlw);
	Convert_Normalized_Local_To_Local_Index(nlu,nlv,nlw,ilu,ilv,ilw);
}

void Global_Coordinate_System::Convert_Local_Index_To_Global(
		int ilu,
		int ilv,
		int ilw,
		double &g0,
		double &g1,
		double &g2
		)
{
	double nlu, nlv, nlw;
	Convert_Local_Index_To_Normalized_Local(ilu,ilv,ilw,nlu,nlv,nlw);
	Convert_Normalized_Local_To_Global(nlu,nlv,nlw,g0,g1,g2);
}

bool Global_Coordinate_System::Set_Transpose(const char* transpose)
{
	int i = _Transpose_String2Idx(transpose);
	if (i >= 0)
	{
		_transpose = i;
		return true;
	}
	return false;
}

bool Global_Coordinate_System::U_Is_Z()
{
	if (
			_transpose == _Transpose_String2Idx("zxy") || 
			_transpose == _Transpose_String2Idx("zyx"))
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool Global_Coordinate_System::V_Is_Z()
{
	if (
			_transpose == _Transpose_String2Idx("xzy") || 
			_transpose == _Transpose_String2Idx("yzx"))
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool Global_Coordinate_System::W_Is_Z()
{
	if (
			_transpose == _Transpose_String2Idx("xyz") || 
			_transpose == _Transpose_String2Idx("yxz"))
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Global_Coordinate_System::Convert_Local_To_Transposed_Fractional_Index(
		double lu,
                double lv,
                double lw,
                double &x,
                double &y,
                double &z
                )
{
	double nlu, nlv, nlw;
        Convert_Local_To_Normalized_Local(lu,lv,lw,nlu,nlv,nlw);
        double filu, filv, filw;
        Convert_Normalized_Local_To_Fractional_Local_Index(nlu,nlv,nlw,filu,filv,filw);
	Convert_Fractional_Local_Index_To_Transposed_Fractional_Index(filu,filv,filw,x,y,z);
}

void Global_Coordinate_System::Convert_Transposed_Fractional_Index_To_Local(
		double x,
		double y,
		double z,
		double &lu,
		double &lv,
		double &lw
		)
{
	if (_Transpose_Is_Valid())
	{
		double filu,filv,filw;
		Convert_Transposed_Fractional_Index_To_Fractional_Local_Index(x,y,z,filu,filv,filw);
		double nlu, nlv, nlw;
		Convert_Fractional_Local_Index_To_Normalized_Local(filu,filv,filw,nlu,nlv,nlw);
		Convert_Normalized_Local_To_Local(nlu,nlv,nlw,lu,lv,lw);
	}
	else
	{
		lu = 0.0;
		lv = 0.0;
		lw = 0.0;
	}
}

void Global_Coordinate_System::Convert_Fractional_Local_Index_To_Transposed_Fractional_Index(
		double filu,
		double filv,
		double filw,
		double &x,
		double &y,
		double &z
		)
{
	switch (_transpose)
	{
		case 0: // zyx
			z = filu;
			y = filv;
			x = filw;
			break;
		case 1: // zxy
			z = filu;
			x = filv;
			y = filw;
			break;
		case 2: // yzx
			y = filu;
			z = filv;
			x = filw;
			break;
		case 3: // yxz
			y = filu;
			x = filv;
			z = filw;
			break;
		case 4: // xzy
			x = filu;
			z = filv;
			y = filw;
			break;
		case 5: // xyz
			x = filu;
			y = filv;
			z = filw;
			break;
		default:
			x = 0.0f;
			y = 0.0f;
			z = 0.0f;
			break;
	}
}

void Global_Coordinate_System::Convert_Transposed_Fractional_Index_To_Fractional_Local_Index(
		double x,
		double y,
		double z,
		double &filu,
		double &filv,
		double &filw
		)
{
	switch (_transpose)
	{
		case 0: // zyx
			filu = z;
			filv = y;
			filw = x;
			break;
		case 1: // zxy
			filu = z;
			filv = x;
			filw = y;
			break;
		case 2: // yzx
			filu = y;
			filv = z;
			filw = x;
			break;
		case 3: // yxz
			filu = y;
			filv = x;
			filw = z;
			break;
		case 4: // xzy
			filu = x;
			filv = z;
			filw = y;
			break;
		case 5: // xyz
			filu = x;
			filv = y;
			filw = z;
			break;
		default:
			filu = 0.0;
			filv = 0.0;
			filw = 0.0;
			break;
	}
}

void Global_Coordinate_System::Convert_Global_To_Transposed_Fractional_Index(
		double g0,
		double g1,
		double g2,
		double &x,
		double &y,
		double &z
		)
{
	double nlu, nlv, nlw;
	Convert_Global_To_Normalized_Local(g0,g1,g2,nlu,nlv,nlw);
	double filu, filv, filw;
	Convert_Normalized_Local_To_Fractional_Local_Index(nlu,nlv,nlw,filu,filv,filw);
	Convert_Fractional_Local_Index_To_Transposed_Fractional_Index(filu,filv,filw,x,y,z);
}

void Global_Coordinate_System::Convert_Transposed_Fractional_Index_To_Global(
		double x,
		double y,
		double z,
		double &g0,
		double &g1,
		double &g2
		)
{
	
	double filu, filv, filw;
	Convert_Transposed_Fractional_Index_To_Fractional_Local_Index(x,y,z,filu,filv,filw);
	double nlu, nlv, nlw;
	Convert_Fractional_Local_Index_To_Normalized_Local(filu,filv,filw,nlu,nlv,nlw);
	Convert_Normalized_Local_To_Global(nlu,nlv,nlw,g0,g1,g2);
}

void Global_Coordinate_System::Convert_Global_To_Transposed_Index(
		double g0,
		double g1,
		double g2,
		int &ix,
		int &iy,
		int &iz
		)
{
	int ilu, ilv, ilw;
	Convert_Global_To_Local_Index(g0,g1,g2,ilu,ilv,ilw);
	Convert_Local_Index_To_Transposed_Index(ilu,ilv,ilw,ix,iy,iz);
}

void Global_Coordinate_System::Convert_Transposed_Index_To_Global(
		int ix,
		int iy,
		int iz,
		double &g0,
		double &g1,
		double &g2
		)
{
	int ilu, ilv, ilw;
	Convert_Transposed_Index_To_Local_Index(ix,iy,iz,ilu,ilv,ilw);
	Convert_Local_Index_To_Global(ilu,ilv,ilw,g0,g1,g2);
}

void Global_Coordinate_System::_Read(const char* path)
{
	FILE* fp = fopen(path, "r");
	if (fp != 0L)
	{
		char str[1024];
		for (char* s = fgets(str, 1024, fp);  s != 0L;  s = fgets(str, 1024, fp))
		{
			char token[1024];
			double v0, v1, v2;
			if (sscanf(s, "%s %lf %lf %lf", token, &v0, &v1, &v2) == 4)
			{
				if (strcmp(token, "AXIS_O") == 0)
				{
					_O0 = v0;
					_O1 = v1;
					_O2 = v2;
				}
				else if(strcmp(token, "AXIS_U") == 0)
				{
					_U0 = v0;
					_U1 = v1;
					_U2 = v2;
				}
				else if(strcmp(token, "AXIS_V") == 0)
				{
					_V0 = v0;
					_V1 = v1;
					_V2 = v2;
				}
				else if(strcmp(token, "AXIS_W") == 0)
				{
					_W0 = v0;
					_W1 = v1;
					_W2 = v2;
				}
				else if(strcmp(token, "AXIS_MIN") == 0)
				{
					_LUMIN = v0;
					_LVMIN = v1;
					_LWMIN = v2;
				}
				else if(strcmp(token, "AXIS_MAX") == 0)
				{
					_LUMAX = v0;
					_LVMAX = v1;
					_LWMAX = v2;
				}
				else if(strcmp(token, "AXIS_D") == 0)
				{
					_DU = v0;
					_DV = v1;
					_DW = v2;
				}
			}
			int i0, i1, i2;
			if (sscanf(s, "%s %d %d %d", token, &i0, &i1, &i2) == 4)
                        {
                                if (strcmp(token, "AXIS_N") == 0)
                                {
                                        _NU = i0;
                                        _NV = i1;
                                        _NW = i2;
                                }
			}
		}
		fclose(fp);
	}
};

int Global_Coordinate_System::_Transpose_String2Idx(const char* str)
{
	const char* s = 0L;
	for (int i = 0;  (s = _Transpose_Idx2String(i)) != 0L;  ++i)
	{
		if (strcmp(s, str) == 0)
		{
			return i;
		}
	}
	return -1;	
}

const char* Global_Coordinate_System::_Transpose_Idx2String(int index)
{
	switch (index)
	{
		case 0:
			return "zyx";
		case 1:
			return "zxy";
		case 2:
			return "yzx";
		case 3:
			return "yxz";
		case 4:
			return "xzy";
		case 5:
			return "xyz";
		default:
			return 0L;
	}
}

const char* Global_Coordinate_System::Get_Axis_Labels()
{
	switch (_transpose)
        {
                case 0:
                        return "Z Y X";
                case 1:
                        return "Z X Y";
                case 2:
                        return "Y Z X";
                case 3:
                        return "Y X Z";
                case 4:
                        return "X Z Y";
                case 5:
                        return "X Y Z";
                default:
                        return 0L;
        }
}

bool Global_Coordinate_System::_Transpose_Is_Valid()
{
	return _Transpose_Idx2String(_transpose) != 0L ? true : false;
}

void Global_Coordinate_System::Dump()
{
	printf("Global_Coordinate_System instance %p\n",this);
	printf("_NU=%d, _NV=%d, _NW=%d, _DU=%lf, _DV=%lf, _DW=%lf\n",_NU,_NV,_NW,_DU,_DV,_DW);
	printf("_O=[%f,%f,%f]\n",_O0,_O1,_O2);
	printf("_U=[%f,%f,%f]\n",_U0,_U1,_U2);
	printf("_V=[%f,%f,%f]\n",_V0,_V1,_V2);
	printf("_W=[%f,%f,%f]\n",_W0,_W1,_W2);
	printf("_MIN=[%f,%f,%f]\n",_LUMIN,_LVMIN,_LWMIN);
	printf("_MAX=[%f,%f,%f]\n",_LUMAX,_LVMAX,_LWMAX);
	printf("_transpose = %d (uvw -> %s)\n",_transpose,_Transpose_Idx2String(_transpose));
	printf("nx=%d, ny=%d, nz=%d, dx=%lf, dy=%lf, dz=%lf\n",
		Get_NX(),Get_NY(),Get_NZ(),Get_DX(),Get_DY(),Get_DZ()
		);
}

