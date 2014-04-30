#ifndef ATTRIBUTE_TRANSFORMER_HXX
#define ATTRIBUTE_TRANSFORMER_HXX

#include <math.h>
#include <stdio.h>

class Transform_Q
{
public:
	Transform_Q(
		float fq,
		float dt
		);
	virtual ~Transform_Q();

	void Transform(float Q, float& inv_Q);	

	float Get_FQ() {return _fq;}
	float Get_DT() {return _dt;}

private:
	float _fq;
	float _dt;
};

class Transform_Dip_Azm
{
public:
	Transform_Dip_Azm(
		int degrees_Flag,
		int dipxdipy_Flag
		);
	virtual ~Transform_Dip_Azm();
	
	void Transform(float dip_or_dx, float azm_or_dy, float& dip, float& azm);

	int Get_Degrees_Flag() {return _degrees_Flag;}
	int Get_DipxDipy_Flag() {return _dipxdipy_Flag;}

private:
	int _degrees_Flag;
	int _dipxdipy_Flag;
};

class Transform_Eps
{
public:
	Transform_Eps(
		int eta_Flag
		);
	virtual ~Transform_Eps();

	void Transform(float eps_or_eta, float dta, float& eps);

	int Get_ETA_Flag() {return _eta_Flag;}

private:
	int _eta_Flag;
};

#endif

