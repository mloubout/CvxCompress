#include <stdio.h>
#include "CvxCompress.hxx"

int main(int argc, char* argv[])
{
	CvxCompress* compressor = new CvxCompress();
	compressor->Run_Module_Tests(false,false);
	return 0;
}
