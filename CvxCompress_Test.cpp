#include <stdio.h>
#include "CvxCompress.hxx"

int main(int argc, char* argv[])
{
  CvxCompress compressor;
  bool passed_tests = compressor.Run_Module_Tests(false,false);
  if (passed_tests) return 0; else return 1;
}
