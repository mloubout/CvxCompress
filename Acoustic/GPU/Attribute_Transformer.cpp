#include "Attribute_Transformer.hxx"

Transform_Q::Transform_Q(
		float fq,
		float dt
		)
{
	_fq = fq;
	_dt = dt;
}

Transform_Q::~Transform_Q()
{
}

void Transform_Q::Transform(float Q, float& inv_Q)
{
       // transform Qatten into a scalar multiplier at each node location: this makes Q proportional to freq,
        // because only one dominant freq is used to represent all frequencies
        // Qatten = 1-(PI/Q)*dt/Tdom = 1-PI/Q*dt*fmax/3 = 1-0.333*PI*fmax*dt/Q (applied to amps, not energy)
        //val = (val > 1e8f) ? 1.0f : 1.0f - 0.333f*M_PI*fmax*dt/val;
	inv_Q = (Q > 1e8f) ? 1.0f : exp(-M_PI*_fq*_dt/Q);
}

Transform_Dip_Azm::Transform_Dip_Azm(
		int degrees_Flag,
		int dipxdipy_Flag
		)
{
	_degrees_Flag = degrees_Flag;
	_dipxdipy_Flag = dipxdipy_Flag;
}

Transform_Dip_Azm::~Transform_Dip_Azm()
{
}

void Transform_Dip_Azm::Transform(float dip_or_dx, float azm_or_dy, float& dip, float& azm)
{
	dip = dip_or_dx;
	azm = azm_or_dy;
	if(_degrees_Flag)
	{
		const float d2r = 0.01745329251994329509f;
		dip *= d2r;
		azm *= d2r;
	}
	if(_dipxdipy_Flag) // input is dipx(in VpIn) & dipy(in EpsIn), so convert to dip,azm
	{
		float tandx = tanf(dip);
		float tandy = tanf(azm);
		dip = sqrtf(tandx*tandx + tandy*tandy);
		azm = atan2f(tandy,tandx);
	}
}

Transform_Eps::Transform_Eps(
		int eta_Flag
		)
{
	_eta_Flag = eta_Flag;
}

Transform_Eps::~Transform_Eps()
{
}

void Transform_Eps::Transform(float eps_or_eta, float dta, float& eps)
{
	eps = eps_or_eta;
	if(_eta_Flag) // input is really eta, so convert to eps
		eps = 0.5f *((1.0f + 2.0f * eps) * (1.0f + 2.0f * dta) - 1.0f);
}

