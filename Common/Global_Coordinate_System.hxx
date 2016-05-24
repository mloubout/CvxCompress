#ifndef CVX_ESDRD_MI_TMJ_GLOBAL_COORDINATE_SYSTEM_HXX
#define CVX_ESDRD_MI_TMJ_GLOBAL_COORDINATE_SYSTEM_HXX

class Global_Coordinate_System
{
public:
	Global_Coordinate_System(const char* path);
	~Global_Coordinate_System();

	// Dump debug information.
	void Dump();

	int Get_NU();
	int Get_NV();
	int Get_NW();

	double Get_DU();
	double Get_DV();
	double Get_DW();

	int Get_NX();
	int Get_NY();
	int Get_NZ();

	double Get_DX();
	double Get_DY();
	double Get_DZ();

	void Convert_Global_To_Normalized_Local(
			double g0,
			double g1,
			double g2,
			double &nlu,
			double &nlv,
			double &nlw
			);

	void Convert_Normalized_Local_To_Global(
			double nlu,
			double nlv,
			double nlw,
			double &g0,
			double &g1,
			double &g2
			);

	void Convert_Local_To_Normalized_Local(
			double lu,
			double lv,
			double lw,
			double &nlu,
			double &nlv,
			double &nlw
			);

	void Convert_Normalized_Local_To_Local(
			double nlu,
			double nlv,
			double nlw,
			double &lu,
			double &lv,
			double &lw
			);

	void Convert_Normalized_Local_To_Local_Index(
			double nlu,
			double nlv,
			double nlw,
			int &ilu,
			int &ilv,
			int &ilw
			);

	void Convert_Local_Index_To_Normalized_Local(
			int ilu,
			int ilv,
			int ilw,
			double &nlu,
			double &nlv,
			double &nlw
			);

	void Convert_Normalized_Local_To_Fractional_Local_Index(
			double nlu,
			double nlv,
			double nlw,
			double &filu,
			double &filv,
			double &filw
			);

	void Convert_Fractional_Local_Index_To_Normalized_Local(
			double filu,
			double filv,
			double filw,
			double &nlu,
			double &nlv,
			double &nlw
			);

	void Convert_Fractional_Local_Index_To_Local_Index(
			double filu,
			double filv,
			double filw,
			int& ilu,
			int& ilv,
			int& ilw
			);

	void Convert_Local_Index_To_Fractional_Local_Index(
			int ilu,
			int ilv,
			int ilw,
			double &filu,
			double &filv,
			double &filw
			);

	void Convert_Global_To_Local(
		double g0, 
		double g1, 
		double g2,
		double &lu, 
		double &lv,
		double &lw
		);

	void Convert_Local_To_Global(
		double lu,
		double lv,
		double lw,
		double &g0,
		double &g1,
		double &g2
		);

	void Convert_Global_To_Local_Index(
		double g0,
		double g1,
		double g2,
		int &ilu,
		int &ilv,
		int &ilw
		);

	void Convert_Local_Index_To_Global(
		int ilu,
		int ilv,
		int ilw,
		double &g0,
		double &g1,
		double &g2
		);

	//
	// Choose desire transpose from U-V-W to X-Y-Z.
	// Provided as a string, the following options are possible:
	//
	// zyx	u-v-w -> z-y-x
	// zxy	u-v-w -> z-x-y
	// yzx	u-v-w -> y-z-x
	// yxz	u-v-w -> y-x-z
	// xzy	u-v-w -> x-z-y
	// xyz	u-v-w -> x-y-z
	//
	// The default transpose is xyz.
	//
	// Return value:
	// true if transpose was set correctly, false otherwise.
	//
	bool Set_Transpose(const char* transpose);

	bool U_Is_Z();
	bool V_Is_Z();
	bool W_Is_Z();

	void Convert_Local_To_Transposed_Fractional_Index(
			double lu,
			double lv,
			double lw,
			double &x,
			double &y,
			double &z
			);

	void Convert_Transposed_Fractional_Index_To_Local(
			double x,
			double y,
			double z,
			double &lu,
			double &lv,
			double &lw
			);

	void Convert_Fractional_Local_Index_To_Transposed_Fractional_Index(
			double filu,
			double filv,
			double filw,
			double &x,
			double &y,
			double &z
			);

	void Convert_Transposed_Fractional_Index_To_Fractional_Local_Index(
			double x,
			double y,
			double z,
			double &filu,
			double &filv,
			double &filw
			);

	void Convert_Local_Index_To_Transposed_Index(
			int ilu,
			int ilv,
			int ilw,
			int &ix,
			int &iy,
			int &iz
			);

	void Convert_Transposed_Index_To_Local_Index(
			int ix,
			int iy,
			int iz,
			int &ilu,
			int &ilv,
			int &ilw
			);

	void Convert_Global_To_Transposed_Fractional_Index(
		double g0,
		double g1,
		double g2,
		double &x,
		double &y,
		double &z
		);

	void Convert_Transposed_Fractional_Index_To_Global(
		double x,
		double y,
		double z,
		double &g0,
		double &g1,
		double &g2
		);

	void Convert_Global_To_Transposed_Index(
		double g0,
		double g1,
		double g2,
		int &ix,
		int &iy,
		int &iz
		);
	
	void Convert_Transposed_Index_To_Global(
		int ix,
		int iy,
		int iz,
		double &g0,
		double &g1,
		double &g2
		);

private:
	int _transpose;
	int _Transpose_String2Idx(const char* str);
	const char* _Transpose_Idx2String(int index);
	bool _Transpose_Is_Valid();

	void _Read(const char* path);

	int _NU, _NV, _NW;
	double _DU, _DV, _DW;

	double _O0, _O1, _O2;
	double _U0, _U1, _U2;
	double _V0, _V1, _V2;
	double _W0, _W1, _W2;

	double _LUMIN, _LUMAX;
	double _LVMIN, _LVMAX;
	double _LWMIN, _LWMAX;
};

#endif

