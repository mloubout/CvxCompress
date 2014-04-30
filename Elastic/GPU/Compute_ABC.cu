//
// This method computes deta(x,y,z) on-the-fly.
//

__device__ 
static float 
Compute_ABC(
	int x,
	int y,
	int z,
	int vol_nx,
	int vol_ny,
	int vol_nz,
	int nabc_top,
	int nabc_bot,
	int nabc_sdx,
	int nabc_sdy,
	float vpvert_avtop,
	float vpvert_avbot,
	float inv_DX,
	float inv_DY,
	float inv_DZ
	)
{
	bool Is_West = x < nabc_sdx;
	bool Is_East = x >= vol_nx - nabc_sdx;	

	bool Is_South = y < nabc_sdy;
	bool Is_North = y >= vol_ny - nabc_sdy;
	
	bool Is_Top = z < nabc_top;
	bool Is_Bot = z >= vol_nz - nabc_bot;

	bool Is_X = Is_West || Is_East;
	bool Is_Y = Is_South || Is_North;
	bool Is_Z = Is_Top || Is_Bot;

	if ( !Is_X && !Is_Y && !Is_Z )
	{
		// no dampening
		return 0.0f;
	}
	else
	{
		if (Is_Z)
		{
			float zr = Is_Top ? (float)(nabc_top - z - 1) / (float)nabc_top : (float)(nabc_bot - vol_nz + z) / (float)nabc_bot;
			float deta_max = Is_Top ? vpvert_avtop * 18.0f * inv_DZ / (float)nabc_top : vpvert_avbot * 18.0f * inv_DZ / (float)nabc_bot;
			if (Is_X && Is_Y)
			{
				// cube
				float xr = Is_West ? (float)(nabc_sdx - x - 1) / (float)nabc_sdx : (float)(nabc_sdx - vol_nx + x) / (float)nabc_sdx;
				float yr = Is_South ? (float)(nabc_sdy - y - 1) / (float)nabc_sdy : (float)(nabc_sdy - vol_ny + y) / (float)nabc_sdy;
				xr = xr * xr * xr;
				xr = xr * xr;
				yr = yr * yr * yr;
				yr = yr * yr;
				zr = zr * zr * zr;
				zr = zr * zr;
				return powf(deta_max*deta_max*(xr+yr+zr),1.0f/3.0f);
			}
			else if (Is_Y)
			{
				// south or north beam
				float yr = Is_South ? (float)(nabc_sdy - y - 1) / (float)nabc_sdy : (float)(nabc_sdy - vol_ny + y) / (float)nabc_sdy;
				yr = yr * yr;
				yr = yr * yr;
				zr = zr * zr;
				zr = zr * zr;
				return sqrtf(deta_max*deta_max*(yr+zr));
			}	
			else if (Is_X)
			{
				// west or east beam
				float xr = Is_West ? (float)(nabc_sdx - x - 1) / (float)nabc_sdx : (float)(nabc_sdx - vol_nx + x) / (float)nabc_sdx;
				xr = xr * xr;
				xr = xr * xr;
				zr = zr * zr;
				zr = zr * zr;
				return sqrtf(deta_max*deta_max*(xr+zr));
			}
			else
			{
				// top or bottom plate
				return deta_max*zr*zr;
			}
		}
		else if (Is_Y) // south or north plate
		{
			float yr = Is_South ? (float)(nabc_sdy - y - 1) / (float)nabc_sdy : (float)(nabc_sdy - vol_ny + y) / (float)nabc_sdy;
			float deta_maxsdy_top = vpvert_avtop * 18.0f * inv_DY / (float)nabc_sdy;
			float deta_maxsdy_bot = vpvert_avbot * 18.0f * inv_DY / (float)nabc_sdy;
			float deta_grady = (deta_maxsdy_bot - deta_maxsdy_top) / (float)(vol_nz-1+nabc_bot);
			float deta_maxsdy = deta_maxsdy_top + deta_grady * (float)(z-nabc_top);
			return deta_maxsdy*yr*yr;
		}
		else // west or east plate
		{
			float xr = Is_West ? (float)(nabc_sdx - x - 1) / (float)nabc_sdx : (float)(nabc_sdx - vol_nx + x) / (float)nabc_sdx;
			float deta_maxsdx_top = vpvert_avtop * 18.0f * inv_DX / (float)nabc_sdx;
			float deta_maxsdx_bot = vpvert_avbot * 18.0f * inv_DX / (float)nabc_sdx;
			float deta_gradx = (deta_maxsdx_bot - deta_maxsdx_top) / (float)(vol_nz-1+nabc_bot);
			float deta_maxsdx = deta_maxsdx_top + deta_gradx * (float)(z-nabc_top);
			return deta_maxsdx*xr*xr;
		}
	}
}

