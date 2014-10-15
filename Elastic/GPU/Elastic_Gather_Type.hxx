#ifndef CVX_ESDRD_MI_TMJ_ELASTIC_GATHER_TYPE_HXX
#define CVX_ESDRD_MI_TMJ_ELASTIC_GATHER_TYPE_HXX

enum Elastic_Gather_Type_t
{
	Common_Shot_Gather = 1,
	Common_Receiver_Gather = 2
};

static char* ToString_Elastic_Gather_Type_t(Elastic_Gather_Type_t gather_type)
{
	switch (gather_type)
	{
	case Common_Shot_Gather:
		return "Common_Shot_Gather";
	case Trilinear:
		return "Common_Receiver_Gather";
	default:
		return "??";
	}
}

#endif

