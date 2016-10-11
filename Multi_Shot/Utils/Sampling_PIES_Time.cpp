#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

void Time2Str(struct tm da_time, char* str)
{
	char* p = asctime(&da_time);
        int l = strlen(p);
        memcpy((void*)str,p,l-1);
        str[l-1] = '\0';
}

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		printf("\nUsage : %s <time-beg> <time-end> <time-delta>\n",argv[0]);
		printf("          <time-beg> Start time in DD-MMM-YY HH:MM format.\n");
		printf("          <time-end> End time in DD-MMM-YY HH:MM format.\n");
		printf("          <time-delta> Time delta in minutes.\n");
		printf("\n");
		return -1;
	}

	setenv("TZ", "UTC", 1);  // TMJ - We override local timezone, force it to UTC since all field timestamps are in UTC
	struct tm start_tm;
	strptime(argv[1],"%d-%b-%y %H:%M",&start_tm);
	start_tm.tm_sec = 0;
	struct tm ref_tm;
	ref_tm = start_tm;
	ref_tm.tm_sec = 0;
	ref_tm.tm_min = 0;
	ref_tm.tm_hour = 0;
	ref_tm.tm_mday = 1;
	ref_tm.tm_mon = 0;
	ref_tm.tm_wday = 0;
	ref_tm.tm_yday = 0;
	struct tm end_tm;
	strptime(argv[2],"%d-%b-%y %H:%M",&end_tm);
	end_tm.tm_sec = 0;
	double time_delta = atof(argv[3]);

	char start_str[256];
	Time2Str(start_tm,start_str);
	char end_str[256];
	Time2Str(end_tm,end_str);

	printf("Generating shot times from %s to %s in increments of %.2f minutes.\n",start_str,end_str,time_delta);

	time_t ref_time = mktime(&ref_tm);
	time_t start_time = mktime(&start_tm);
	time_t end_time = mktime(&end_tm);
	
	int ffid = 1;
	for (time_t curr_time = start_time;  curr_time <= end_time;  curr_time = curr_time + time_delta*60.0)
	{
		struct tm curr_tm = *localtime(&curr_time);
		char curr_str[256];
		Time2Str(curr_tm,curr_str);
		double td = difftime(curr_time,ref_time);
		int jday = td / 86400.0;
		td = td - (double)jday * 86400.0;
		int jhour = td / 3600.0;
		td = td - (double)jhour * 3600.0;
		int jmin = td / 60.0;
		td = td - (double)jmin * 60.0;
		printf("%05d %s %d %d %d %.2f\n",ffid,curr_str,jday,jhour,jmin,td);
		++ffid;
	}

	return 0;
}
