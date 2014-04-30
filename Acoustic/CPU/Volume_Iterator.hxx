#ifndef VOLUME_ITERATOR_HXX
#define VOLUME_ITERATOR_HXX

#include "Attribute_Reader.hxx"

//
// Class designed to iterate over a voxet style volume.
// Supports sub-volumes.
//

class Volume_Iterator
{
	public:
		Volume_Iterator(
			int log_level,
			int fast_Axis,
			int med_Axis,
			int slow_Axis,
			int dimx,
			int dimy,
			int dimz,
			int x0,
			int x1,
			int y0,
			int y1,
			int z0,
			int z1
			);
		virtual ~Volume_Iterator();

		void Iterate(Attribute_Reader* reader);

	private:
		int _log_level;
		UVW_XYZ_Calculator* _ucalc;
		int _ubeg, _uend, _dimu, _actual_dimu;
		int _vbeg, _vend, _dimv;
		int _wbeg, _wend, _dimw;
};

#endif

