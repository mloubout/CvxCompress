//
// Earth model packing:
//
// Earth model has the following parameters:
//
// Vp		P-wave velocity
// Vs		S-wave velocity
// Density	
// Delta1
// Delta2
// Delta3
// Epsilon1
// Epsilon2
// Gamma1
// Gamma2
// Q
// Dip
// Azimuth
// Rake
//
// int 0: Vp (16 bits), Vs (16 bits)
// int 1: Rake (8 bits), Delta1 (8 bits), Delta2 (8 bits), Delta3 (8 bits)
// int 2: Epsilon1 (8 bits), Epsilon2 (8 bits), Gamma1 (8 bits), Gamma2 (8 bits)
// int 3: Q (8 bits), Density (8 bits), Dip (8 bits), Azimuth (8 bits)
//

__device__ 
void cuUnpack_Q_Density(
	unsigned int Word3,
	float Q_min,
	float Q_range,
	float* Q,
	float Density_min,
	float Density_range,
	float* Density
	)
{
	*Q = Q_min + (float)((Word3 >> 24) & 255) * Q_range;
	*Density = 1000.0f * (Density_min + (float)((Word3 >> 16) & 255) * Density_range);
}

__device__ 
void cuUnpack_Everything( // well, everything except Q.
	unsigned int Word0,
	unsigned int Word1,
	unsigned int Word2,
	unsigned int Word3,
	float Vp_min,
	float Vp_range,
	float* Vp,
	float Vs_min,
	float Vs_range,
	float* Vs,
	float Density_min,
	float Density_range,
	float* Density,
	float Dip_min,
	float Dip_range,
	float* Dip,
	float Azimuth_min,
	float Azimuth_range,
	float* Azimuth,
	float Rake_min,
	float Rake_range,
	float* Rake,
	float Delta1_min,
	float Delta1_range,
	float* Delta1,
	float Delta2_min,
	float Delta2_range,
	float* Delta2,
	float Delta3_min,
	float Delta3_range,
	float* Delta3,
	float Epsilon1_min,
	float Epsilon1_range,
	float* Epsilon1,
	float Epsilon2_min,
	float Epsilon2_range,
	float* Epsilon2,
	float Gamma1_min,
	float Gamma1_range,
	float* Gamma1,
	float Gamma2_min,
	float Gamma2_range,
	float* Gamma2
	)
{
	*Dip = Dip_min + (float)((Word3 >> 8) & 255) * Dip_range;
	*Azimuth = Azimuth_min + (float)(Word3 & 255) * Azimuth_range;
	*Rake = Rake_min + (float)((Word1 >> 24) & 255) * Rake_range;
	*Delta1 = Delta1_min + (float)((Word1 >> 16) & 255) * Delta1_range;
	*Delta2 = Delta2_min + (float)((Word1 >> 8) & 255) * Delta2_range;
	*Delta3 = Delta3_min + (float)(Word1 & 255) * Delta3_range;
	*Epsilon1 = Epsilon1_min + (float)((Word2 >> 24) & 255) * Epsilon1_range;
	*Epsilon2 = Epsilon2_min + (float)((Word2 >> 16) & 255) * Epsilon2_range;
	*Gamma1 = Gamma1_min + (float)((Word2 >> 8) & 255) * Gamma1_range;
	*Gamma2 = Gamma2_min + (float)(Word2 & 255) * Gamma2_range;
	*Density = 1000.0f * (Density_min + (float)((Word3 >> 16) & 255) * Density_range);
	*Vp = Vp_min + (float)((Word0 >> 16) & 65535) * Vp_range;
	*Vs = Vs_min + (float)(Word0 & 65535) * Vs_range;
}

__device__ 
void cuUnpack_And_Compute_On_Kite_CIJs(
	unsigned int Word0,
	unsigned int Word1,
	unsigned int Word2,
	unsigned int Word3,
	float Vp_min,
	float Vp_range,
	float Vs_min,
	float Vs_range,
	float Density_min,
	float Density_range,
	float Dip_min, 
	float Dip_range, 
	float Azimuth_min, 
	float Azimuth_range, 
	float Rake_min, 
	float Rake_range, 
	float Delta1_min,
	float Delta1_range,
	float Delta2_min,
	float Delta2_range,
	float Delta3_min,
	float Delta3_range,
	float Epsilon1_min,
	float Epsilon1_range,
	float Epsilon2_min,
	float Epsilon2_range,
	float Gamma1_min,
	float Gamma1_range,
	float Gamma2_min,
	float Gamma2_range,
	float* Dip,
	float* Azimuth,
	float* Rake,
	float* c11,		// c11 = c33 * (1 + 2 * Epsilon2)
	float* c22,		// c22 = c33 * (1 + 2 * Epsilon1)
	float* c33,		// c33 = bulk modulus = Vp^2 * Density
	float* c44,		// c44 = c66 / (1 + 2 * Gamma2)
	float* c55,		// c55 = Vs^2 * Density
	float* c66,		// c66 = c55 * (1 + 2 * Gamma1)
	float* c12,		// c12 = sqrt((c11 - c66)*(2*Delta3*c11 + c11 - c66)) - c66
	float* c13,		// c13 = sqrt((c33 - c55)*(2*Delta2*c33 + c33 - c55)) - c55
	float* c23		// c23 = sqrt((c33 - c44)*(2*Delta1*c33 + c33 - c44)) - c44
	)
{
	float Vp, Vs, Density, Delta1, Delta2, Delta3, Epsilon1, Epsilon2, Gamma1, Gamma2;
	cuUnpack_Everything(
		Word0,Word1,Word2,Word3,
		Vp_min,Vp_range,&Vp,
		Vs_min,Vs_range,&Vs,
		Density_min,Density_range,&Density,
		Dip_min,Dip_range,Dip,
		Azimuth_min,Azimuth_range,Azimuth,
		Rake_min,Rake_range,Rake,
		Delta1_min,Delta1_range,&Delta1,
		Delta2_min,Delta2_range,&Delta2,
		Delta3_min,Delta3_range,&Delta3,
		Epsilon1_min,Epsilon1_range,&Epsilon1,
		Epsilon2_min,Epsilon2_range,&Epsilon2,
		Gamma1_min,Gamma1_range,&Gamma1,
		Gamma2_min,Gamma2_range,&Gamma2
		);
	*c33 = Vp * Vp * Density;
	*c11 = *c33 * (1.0f + 2.0f * Epsilon2);
	*c22 = *c33 * (1.0f + 2.0f * Epsilon1);
	*c55 = Vs * Vs * Density;
	*c66 = *c55 * (1.0f + 2.0f * Gamma1);
	*c44 = *c66 / (1.0f + 2.0f * Gamma2);
	*c12 = sqrtf((*c11 - *c66) * (2.0f * Delta3 * *c11 + *c11 - *c66)) - *c66;
	*c13 = sqrtf((*c33 - *c55) * (2.0f * Delta2 * *c33 + *c33 - *c55)) - *c55;
	*c23 = sqrtf((*c33 - *c44) * (2.0f * Delta1 * *c33 + *c33 - *c44)) - *c44;
}

__device__ 
void cuUnpack_And_Compute_Rotated_CIJs(
	unsigned int Word0,
	unsigned int Word1,
	unsigned int Word2,
	unsigned int Word3,
	float Vp_min,
	float Vp_range,
	float Vs_min,
	float Vs_range,
	float Density_min,
	float Density_range,
	float Dip_min, 
	float Dip_range, 
	float Azimuth_min, 
	float Azimuth_range, 
	float Rake_min, 
	float Rake_range, 
	float Delta1_min,
	float Delta1_range,
	float Delta2_min,
	float Delta2_range,
	float Delta3_min,
	float Delta3_range,
	float Epsilon1_min,
	float Epsilon1_range,
	float Epsilon2_min,
	float Epsilon2_range,
	float Gamma1_min,
	float Gamma1_range,
	float Gamma2_min,
	float Gamma2_range,
	float* c11,
	float* c12,
	float* c13,
	float* c14,
	float* c15,
	float* c16,
	float* c22,
	float* c23,
	float* c24,
	float* c25,
	float* c26,
	float* c33,
	float* c34,
	float* c35,
	float* c36,
	float* c44,
	float* c45,
	float* c46,
	float* c55,
	float* c56,
	float* c66
	)
{
	float Dip, Azimuth, Rake;
	float _c11, _c22, _c33, _c44, _c55, _c66, _c12, _c13, _c23;
	cuUnpack_And_Compute_On_Kite_CIJs(
		Word0,Word1,Word2,Word3,
		Vp_min,Vp_range,
		Vs_min,Vs_range,
		Density_min,Density_range,
		Dip_min,Dip_range,
		Azimuth_min,Azimuth_range,
		Rake_min,Rake_range,
		Delta1_min,Delta1_range,
		Delta2_min,Delta2_range,
		Delta3_min,Delta3_range,
		Epsilon1_min,Epsilon1_range,
		Epsilon2_min,Epsilon2_range,
		Gamma1_min,Gamma1_range,
		Gamma2_min,Gamma2_range,
		&Dip,&Azimuth,&Rake,
		&_c11,&_c22,&_c33,&_c44,&_c55,&_c66,&_c12,&_c13,&_c23
		);

	float sind, cosd, sina, cosa, sinr, cosr;
	__sincosf(Dip,&sind,&cosd);
	__sincosf(Azimuth,&sina,&cosa);
	__sincosf(Rake,&sinr,&cosr);
	
	//..Euler (Z(2)-Y(1)-Z-order) rotation matrix b: local to global (FD)
	float b11 =  cosa*cosd*cosr - sina*sinr;
	float b12 = -cosa*cosd*sinr - sina*cosr;
	float b13 =  cosa*sind;
	float b21 =  sina*cosd*cosr + cosa*sinr;
	float b22 = -sina*cosd*sinr + cosa*cosr;
	float b23 =  sina*sind;
	float b31 = -sind*cosr;
	float b32 =  sind*sinr;
	float b33 =  cosd;

	//..Bond transformation matrix M:
	float m11 = b11*b11;
	float m12 = b12*b12;
	float m13 = b13*b13;
	float m14 = 2.0f*b12*b13;
	float m15 = 2.0f*b11*b13;
	float m16 = 2.0f*b11*b12;

	float m21 = b21*b21;
	float m22 = b22*b22;
	float m23 = b23*b23;
	float m24 = 2.0f*b22*b23;
	float m25 = 2.0f*b21*b23;
	float m26 = 2.0f*b21*b22;

	float m31 = b31*b31;
	float m32 = b32*b32;
	float m33 = b33*b33;
	float m34 = 2.0f*b32*b33;
	float m35 = 2.0f*b31*b33;
	float m36 = 2.0f*b31*b32;

	float m41 = b21*b31;
	float m42 = b22*b32;
	float m43 = b23*b33;
	float m44 = b22*b33 + b23*b32;
	float m45 = b21*b33 + b23*b31;
	float m46 = b22*b31 + b21*b32;

	float m51 = b11*b31;
	float m52 = b12*b32;
	float m53 = b13*b33;
	float m54 = b12*b33 + b13*b32;
	float m55 = b11*b33 + b13*b31;
	float m56 = b11*b32 + b12*b31;

	float m61 = b11*b21;
	float m62 = b12*b22;
	float m63 = b13*b23;
	float m64 = b13*b22 + b12*b23;
	float m65 = b11*b23 + b13*b21;
	float m66 = b11*b22 + b12*b21;

	//..[M] [CIJ_local] product:
	float mc11 = m11*_c11 + m12*_c12 + m13*_c13;
	float mc12 = m11*_c12 + m12*_c22 + m13*_c23;
	float mc13 = m11*_c13 + m12*_c23 + m13*_c33;
	float mc14 = m14*_c44;
	float mc15 = m15*_c55;
	float mc16 = m16*_c66;

	//..[CIJ_global] = [M] [CIJ_local] [M]^T
	*c11 = mc11*m11 + mc12*m12 + mc13*m13 + mc14*m14 + mc15*m15 + mc16*m16;
	*c12 = mc11*m21 + mc12*m22 + mc13*m23 + mc14*m24 + mc15*m25 + mc16*m26;
	*c13 = mc11*m31 + mc12*m32 + mc13*m33 + mc14*m34 + mc15*m35 + mc16*m36;
	*c14 = mc11*m41 + mc12*m42 + mc13*m43 + mc14*m44 + mc15*m45 + mc16*m46;
	*c15 = mc11*m51 + mc12*m52 + mc13*m53 + mc14*m54 + mc15*m55 + mc16*m56;
	*c16 = mc11*m61 + mc12*m62 + mc13*m63 + mc14*m64 + mc15*m65 + mc16*m66;

	//..[M] [CIJ_local] product:
	float mc21 = m21*_c11 + m22*_c12 + m23*_c13;
	float mc22 = m21*_c12 + m22*_c22 + m23*_c23;
	float mc23 = m21*_c13 + m22*_c23 + m23*_c33;
	float mc24 = m24*_c44;
	float mc25 = m25*_c55;
	float mc26 = m26*_c66;

	//..[CIJ_global] = [M] [CIJ_local] [M]^T
	*c22 = mc21*m21 + mc22*m22 + mc23*m23 + mc24*m24 + mc25*m25 + mc26*m26;
	*c23 = mc21*m31 + mc22*m32 + mc23*m33 + mc24*m34 + mc25*m35 + mc26*m36;
	*c24 = mc21*m41 + mc22*m42 + mc23*m43 + mc24*m44 + mc25*m45 + mc26*m46;
	*c25 = mc21*m51 + mc22*m52 + mc23*m53 + mc24*m54 + mc25*m55 + mc26*m56;
	*c26 = mc21*m61 + mc22*m62 + mc23*m63 + mc24*m64 + mc25*m65 + mc26*m66;

	//..[M] [CIJ_local] product:
	float mc31 = m31*_c11 + m32*_c12 + m33*_c13;
	float mc32 = m31*_c12 + m32*_c22 + m33*_c23;
	float mc33 = m31*_c13 + m32*_c23 + m33*_c33;
	float mc34 = m34*_c44;
	float mc35 = m35*_c55;
	float mc36 = m36*_c66;

	//..[CIJ_global] = [M] [CIJ_local] [M]^T
	*c33 = mc31*m31 + mc32*m32 + mc33*m33 + mc34*m34 + mc35*m35 + mc36*m36;
	*c34 = mc31*m41 + mc32*m42 + mc33*m43 + mc34*m44 + mc35*m45 + mc36*m46;
	*c35 = mc31*m51 + mc32*m52 + mc33*m53 + mc34*m54 + mc35*m55 + mc36*m56;
	*c36 = mc31*m61 + mc32*m62 + mc33*m63 + mc34*m64 + mc35*m65 + mc36*m66;

	//..[M] [CIJ_local] product:
	float mc41 = m41*_c11 + m42*_c12 + m43*_c13;
	float mc42 = m41*_c12 + m42*_c22 + m43*_c23;
	float mc43 = m41*_c13 + m42*_c23 + m43*_c33;
	float mc44 = m44*_c44;
	float mc45 = m45*_c55;
	float mc46 = m46*_c66;

	//..[CIJ_global] = [M] [CIJ_local] [M]^T
	*c44 = mc41*m41 + mc42*m42 + mc43*m43 + mc44*m44 + mc45*m45 + mc46*m46;
	*c45 = mc41*m51 + mc42*m52 + mc43*m53 + mc44*m54 + mc45*m55 + mc46*m56;
	*c46 = mc41*m61 + mc42*m62 + mc43*m63 + mc44*m64 + mc45*m65 + mc46*m66;

	//..[M] [CIJ_local] product:
	float mc51 = m51*_c11 + m52*_c12 + m53*_c13;
	float mc52 = m51*_c12 + m52*_c22 + m53*_c23;
	float mc53 = m51*_c13 + m52*_c23 + m53*_c33;
	float mc54 = m54*_c44;
	float mc55 = m55*_c55;
	float mc56 = m56*_c66;

	//..[CIJ_global] = [M] [CIJ_local] [M]^T
	*c55 = mc51*m51 + mc52*m52 + mc53*m53 + mc54*m54 + mc55*m55 + mc56*m56;
	*c56 = mc51*m61 + mc52*m62 + mc53*m63 + mc54*m64 + mc55*m65 + mc56*m66;

	//..[M] [CIJ_local] product:
	float mc61 = m61*_c11 + m62*_c12 + m63*_c13;
	float mc62 = m61*_c12 + m62*_c22 + m63*_c23;
	float mc63 = m61*_c13 + m62*_c23 + m63*_c33;
	float mc64 = m64*_c44;
	float mc65 = m65*_c55;
	float mc66 = m66*_c66;

	//..[CIJ_global] = [M] [CIJ_local] [M]^T
	*c66 = mc61*m61 + mc62*m62 + mc63*m63 + mc64*m64 + mc65*m65 + mc66*m66;

	/*
	if (_c33 > 63031812.000000f)
	{
		printf("_c11=%f, _c22=%f, _c33=%f, _c44=%f, _c55=%f, _c66=%f, _c12=%f, _c13=%f, _c23=%f\nc16=%f, c15=%f, c14=%f, c13=%f, c12=%f, c11=%f\nc26=%f, c25=%f, c24=%f, c23=%f, c22=%f\nc36=%f, c35=%f, c34=%f, c33=%f\nc46=%f, c45=%f, c44=%f\nc56=%f, c55=%f\nc66=%f\n",
			_c11,_c22,_c33,_c44,_c55,_c66,_c12,_c13,_c23,
			*c16, *c15, *c14, *c13, *c12, *c11,
			*c26, *c25, *c24, *c23, *c22,
			*c36, *c35, *c34, *c33,
			*c46, *c45, *c44,
			*c56, *c55,
			*c66);
		printf("Dip=%f, Azimuth=%f, Rake=%f, sind=%f, cosd=%f, sina=%f, cosa=%f, sinr=%f, cosr=%f\n",
			Dip, Azimuth, Rake,
			sind, cosd, sina, cosa, sinr, cosr
			);
	}
	*/
}

