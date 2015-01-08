#ifndef CVX_ESDRD_MI_TMJ_ELASTIC_INTERPOLATION_HXX
#define CVX_ESDRD_MI_TMJ_ELASTIC_INTERPOLATION_HXX

enum Elastic_Interpolation_t
{
	Point = 1,
	Trilinear = 2,
	Sinc = 3
};

static const char* ToString_Elastic_Interpolation_t(Elastic_Interpolation_t interpolation_method)
{
	switch (interpolation_method)
	{
	case Point:
		return "Point";
	case Trilinear:
		return "Trilinear";
	case Sinc:
		return "Sinc";
	default:
		return "??";
	}
}

#endif

