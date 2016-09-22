#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

std::vector<std::string> split(std::string str)
{
	std::vector<std::string> substrings;
	int mode = 0;
	int start = 0, stop = 0;
	for (int i = 0;  i < str.size();  ++i)
	{
		char c = str[i];
		if (mode == 0)
		{
			// looking for start of next substring
			if (!isspace(c))
			{
				start = i;
				mode = 1;
			}
		}
		else if (mode == 1)
		{
			// looking for end of substring
			if (isspace(c))
			{
				stop = i-1;
				mode = 0;
				std::string sub = str.substr(start,stop-start+1);
				substrings.push_back(sub);
			}
		}
	}
	if (mode == 1)
	{
		std::string sub = str.substr(start,stop-start+1);
		substrings.push_back(sub);
	}
	return substrings;
}

int main(int argc, char** argv)
{
	if (argc != 3)
	{
		printf("\n");
		printf("THIS UTILITY ONLY WORKS FOR THE 2015-2016 GORGON 4DVwxyzt JOB.\n");
		printf("Compiled on 22-Sep-2016\n");
		printf("Usage : %s <input-file> <output-file>\n",argv[0]);
		printf("        <input-file> Original file.\n");
		printf("        <output-file> Converted file.\n");
		printf("\n");
		return -1; 
	}	

	setenv("TZ", "UTC", 1);  // TMJ - We override local timezone, force it to UTC since all field timestamps are in UTC

	FILE* fp = fopen(argv[3],"w");
	assert(fp != 0L);

	std::ifstream ifs(argv[1],std::ifstream::in);
	int isrc = 0;
	while (!ifs.eof())
	{
		++isrc;
		char buf[16385];
		buf[0] = 0;
		ifs.getline(buf,16384);
		if (buf[0] != 0)
		{
			std::vector<std::string> substrings = split((std::string)buf);
			int col = 1;
			for (std::vector<std::string>::iterator it = substrings.begin();  it != substrings.end();  ++it)
			{
				printf("%s ",(*it).c_str());
				++col;
			}
			printf("\n");
			struct tm tm1;
			strptime((substrings[1]+" "+substrings[2]).c_str(),"%d-%b-%y %T",&tm1);
			int year = tm1.tm_year;
			int julian_day = 1 + tm1.tm_yday + (year == 116 ? 365 : 0);
			int hours = tm1.tm_hour;
			int minutes = tm1.tm_min;
			int seconds = tm1.tm_sec;
			fprintf(fp,(substrings[0]+" ").c_str());
			fprintf(fp,"%d %d %d %.2f ",julian_day,hours,minutes,(float)seconds);
			for (int i = 3;  i < substrings.size();  ++i)
			{
				fprintf(fp,(substrings[i]+" ").c_str());
			}
			fprintf(fp," -9.0");  // add source depth.
			fprintf(fp,"\n");
		}
	}

	fclose(fp);

	return 0;
}
