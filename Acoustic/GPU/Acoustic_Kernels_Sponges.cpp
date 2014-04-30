/*
** num_x_zeropad	If > 0, right side of x sponge is shifted this many cells to the left. The remanining cells on the right side are zeroed.
**			This will happen if actual dimx is not a multiple of 4, in which case num_x_zeropad equals dimx - actual_dimx.
*/
void Compute_Sponges(
	float spongecoeff_x,
	float spongecoeff_y,
	float spongecoeff_z_lo,
	float spongecoeff_z_hi,
	int spongewidth_x,
	int spongewidth_y,
	int spongewidth_z_lo,
	int spongewidth_z_hi,
	int absorbz0,
	int dimx,
	int num_x_zeropad,
	int dimy,
	int num_y_zeropad,
	int dimz,
	int num_z_zeropad,
	float* spgx,
	float* spgy,
	float* spgz
	)
{
	int xstart = dimx - spongewidth_x;
	spgx[0] = spgx[dimx-1-num_x_zeropad] = 0.0f;
	for (int i = 1;  i < spongewidth_x;  ++i)
	{
		spgx[i] = 1.0 - spongecoeff_x * (spongewidth_x - i) * (spongewidth_x - i);
	}
	for (int i = spongewidth_x;  i < xstart;  ++i)
	{
		spgx[i] = 1.0f;
	}
	for (int i = xstart;  i < dimx-1;  ++i)
	{
		spgx[i-num_x_zeropad] = 1.0f - spongecoeff_x * (i - xstart + 1) * (i - xstart + 1);
	}
	if (num_x_zeropad > 0)
	{
		for (int i = dimx-1-num_x_zeropad;  i < dimx;  ++i)
		{
			spgx[i] = 0.0f;
		}
	}
	for (int i = 0;  i < dimx;  ++i) printf("spgx[%d] = %f\n",i,spgx[i]);

	int ystart = dimy - spongewidth_y;
	spgy[0] = spgy[dimy-1-num_y_zeropad] = 0.0f;
	for (int i = 1;  i < spongewidth_y;  ++i)
	{
		spgy[i] = 1.0 - spongecoeff_y * (spongewidth_y - i) * (spongewidth_y - i);
	}
	for (int i = spongewidth_y;  i < ystart;  ++i)
	{
		spgy[i] = 1.0f;
	}
	for (int i = ystart;  i < dimy-1;  ++i)
	{
		spgy[i-num_y_zeropad] = 1.0f - spongecoeff_y * (i - ystart + 1) * (i - ystart + 1);
	}
	if (num_y_zeropad > 0)
	{
		for (int i = dimy-1-num_y_zeropad;  i < dimy;  ++i)
		{
			spgy[i] = 0.0f;
		}
	}
	for (int i = 0;  i < dimy;  ++i) printf("spgy[%d] = %f\n",i,spgy[i]);

	int zstart = dimz - spongewidth_z_hi;
	spgz[0] = absorbz0 ? 0.0f : 1.0f;
	spgz[dimz-1-num_z_zeropad] = 0.0f;
	for (int i = 1;  i < spongewidth_z_lo;  ++i)
	{
		spgz[i] = absorbz0 ? 1.0 - spongecoeff_z_lo * (spongewidth_z_lo - i) * (spongewidth_z_lo - i) : 1.0f;
	}
	for (int i = spongewidth_z_lo;  i < zstart;  ++i)
	{
		spgz[i] = 1.0f;
	}
	for (int i = zstart;  i < dimz-1;  ++i)
	{
		spgz[i-num_z_zeropad] = 1.0f - spongecoeff_z_hi * (i - zstart + 1) * (i - zstart + 1);
	}
	if (num_z_zeropad > 0)
	{
		for (int i = dimz-1-num_z_zeropad;  i < dimz;  ++i)
		{
			spgz[i] = 0.0f;
		}
	}
	for (int i = 0;  i < dimz;  ++i) printf("spgz[%d] = %f\n",i,spgz[i]);
}
