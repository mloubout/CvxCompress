#include <stdio.h>
#include <omp.h>

#include "Volume_Iterator.hxx"

Volume_Iterator::Volume_Iterator(
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
		)
{
	_log_level = log_level;
	_ucalc = new UVW_XYZ_Calculator(fast_Axis,med_Axis,slow_Axis);
	_ucalc->Compute_UVW_From_XYZ(dimx,dimy,dimz,_dimu,_dimv,_dimw);
        _ucalc->Compute_UVW_From_XYZ(x0,y0,z0,_ubeg,_vbeg,_wbeg);
        _ucalc->Compute_UVW_From_XYZ(x1,y1,z1,_uend,_vend,_wend);
        _actual_dimu = _uend - _ubeg + 1;
	if (_log_level >= 4)
	{
		printf("dimu = %d, actual_dimu = %d\n",_dimu,_actual_dimu);
		printf("dimv = %d\n",_dimv);
		printf("dimw = %d\n",_dimw);
		printf("ubeg = %d\n",_ubeg);
		printf("uend = %d\n",_uend);
		printf("vbeg = %d\n",_vbeg);
		printf("vend = %d\n",_vend);
		printf("wbeg = %d\n",_wbeg);
		printf("wend = %d\n",_wend);
	}
}

Volume_Iterator::~Volume_Iterator()
{
}

void Volume_Iterator::Iterate(Attribute_Reader** readers, int num_readers)
{
        if (_ubeg < 0 || _uend >= _dimu || _vbeg < 0 || _vend >= _dimv || _wbeg < 0 || _wend >= _dimw)
        {
                printf("ERROR! Read_Earth_Model_Attributes - Sub volume extends beyond dimensions of original volume!\n");
                printf("u = [%d, %d], _dimu = %d\n",_ubeg,_uend,_dimu);
                printf("v = [%d, %d], _dimv = %d\n",_vbeg,_vend,_dimv);
                printf("w = [%d, %d], _dimw = %d\n",_wbeg,_wend,_dimw);
        }

        // loop over u,v,w.
        // read, convert to right endian-ness and transpose to x-y-z ordering.
        int nwpu = (_wend-_wbeg)/100;
        if (nwpu < 1) nwpu = 1;
	omp_set_num_threads(num_readers);
	int cnt = 0;
#pragma omp parallel for
        for (int w = _wbeg;  w <= _wend;  ++w)
        {
                int curr_thread = omp_get_thread_num();
		Attribute_Reader* reader = readers[curr_thread];
		
                for (int v = _vbeg;  v <= _vend;  ++v)
                {
                        long offset = (long)sizeof(float) * ((long)w * (long)_dimu * (long)_dimv + (long)v * (long)_dimu + (long)_ubeg);
                        // read from files. we do this sequentially, since hitting the filesystem with many parallel requests most likely will trash performance.
			reader->Read(offset, 0, _uend - _ubeg, v - _vbeg, w - _wbeg);
                }
#pragma omp critical
		{
			++cnt;
			if (_log_level >= 3 && (cnt % nwpu) == 0)
			{
				printf("\r%s (%.0f%%)",reader->name,100.0f*(float)cnt/(float)(_wend-_wbeg));
				fflush(stdout);
			}
		}
        }
	if (_log_level >= 3) printf("\n");
}

void Volume_Iterator::Iterate(Attribute_Reader* reader)
{
	int num_threads;
#pragma omp parallel
	{
		num_threads = omp_get_num_threads();
	}

	Attribute_Reader** readers = new Attribute_Reader*[num_threads];
	for (int i = 0;  i < num_threads;  ++i)
	{
		readers[i] = reader->clone();
	}

	Iterate(readers,num_threads);
	reader->Reduce(readers,num_threads);

	for (int i = 0;  i < num_threads;  ++i)
        {
		delete readers[i];
	}

	delete [] readers;
}

