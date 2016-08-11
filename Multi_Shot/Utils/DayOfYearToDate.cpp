#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		printf("\n%s Aug-11 2016\nUsage: %s <day> <year>\n\n",argv[0],argv[0]);
		return -1;
	}

	setenv("TZ", "UTC", 1);
	
	int day = atoi(argv[1]);
	int year = atoi(argv[2]);
	
	struct tm tms;
	tms.tm_sec = 0;
	tms.tm_min = 0;
	tms.tm_hour = 0;
	tms.tm_mday = day;
	tms.tm_mon = 0;
	tms.tm_year = year - 1900;
	tms.tm_isdst = 0;

	time_t ts = mktime(&tms);
	printf("Day %d of Year %d is Date %s",day,year,asctime(localtime(&ts)));

	return 0;
}

