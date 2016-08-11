#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <list>
#include <vector>
#include <limits.h>
#include <stdlib.h>

#include <GeomTrace.h>
#include <GeomTraceAuxiliary.h>
#include <ArrND.h>
#include <Common/Voxet.hxx>
#include <Common/Voxet_Property.hxx>
#include <Common/Global_Coordinate_System.hxx>

struct pos
{
	long ffid;
	int seqno;
	int gunseq;
	int il;
	int xl;
	time_t shot_time;
	double shot_time_fractional_usec;
	double x;
	double y;
	double z;
};

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

std::list<pos> read_ascii(
	char* label,
	char* filename,
	int idx_ffid,
	int idx_seqno,
	int idx_gunseq,
	int idx_shot_time,
	int idx_il,
	int idx_xl,
	int idx_X,
	int idx_Y,
	int idx_Z,
	int start_day,
	int start_year
	)
{
	int max_idx = std::max(std::max(idx_X,idx_Y),idx_Z);
	if (idx_ffid > 0) max_idx = std::max(max_idx,idx_ffid);
	if (idx_shot_time > 0) max_idx = std::max(max_idx,idx_shot_time+2);
	std::list<pos> positions;
	std::ifstream ifs(filename,std::ifstream::in);
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
			if (substrings.size() > 0 && substrings.size() < max_idx)
			{
				printf("Error! Too few columns in line (was %d, expected at least %d):\n",substrings.size(),max_idx);
				printf(" ==> %s\n",buf);
				exit(-5);
			}
			pos p;
			p.x = atof(substrings[idx_X-1].c_str());
			p.y = atof(substrings[idx_Y-1].c_str());
			p.z = atof(substrings[idx_Z-1].c_str());
			if (idx_ffid > 0)
			{
				p.ffid = atol(substrings[idx_ffid-1].c_str());
			}
			else
			{
				p.ffid = isrc;
			}
			if (idx_seqno > 0)
			{
				p.seqno = atoi(substrings[idx_seqno-1].c_str());
			}
			else
			{
				p.seqno = 0;
			}
			if (idx_gunseq > 0)
			{
				p.gunseq = atoi(substrings[idx_gunseq-1].c_str());
			}
			else
			{
				p.gunseq = 0;
			}
			if (idx_il > 0)
			{
				p.il = atoi(substrings[idx_il-1].c_str());
			}
			else
			{
				p.il = 0;
			}
			if (idx_xl > 0)
			{
				p.xl = atoi(substrings[idx_xl-1].c_str());
			}
			else
			{
				p.xl = 0;
			}
			if (idx_shot_time > 0)
			{
				struct tm ft;
				ft.tm_mon = 0;
				ft.tm_mday = atoi(substrings[idx_shot_time-1].c_str());
				ft.tm_hour = atoi(substrings[idx_shot_time].c_str());
				ft.tm_min = atoi(substrings[idx_shot_time+1].c_str());
				double fractional_seconds = atof(substrings[idx_shot_time+2].c_str());
				int whole_seconds = (int)trunc(fractional_seconds);
				ft.tm_sec = whole_seconds;
				fractional_seconds -= (double)whole_seconds;
				ft.tm_year = (ft.tm_mday < start_day ? start_year+1 : start_year) - 1900;
				time_t acquisition_time = mktime(&ft);
				p.shot_time = acquisition_time;
				p.shot_time_fractional_usec = fractional_seconds*1e6;
			}
			positions.push_back(p);
			std::string str = (std::string)asctime(localtime(&(p.shot_time)));
			str = str.substr(0,str.size()-1);
			printf("%s (ffid=%d, shot_time=%s) : %.2f %.2f %.2f\n",label,p.ffid,str.c_str(),p.x,p.y,p.z);
		}
	}
	return positions;
}

int main(int argc, char* argv[])
{
	if (argc != 24)
	{
		printf("\n");
		printf("Usage : %s \n",argv[0]);
		printf("        <source-file> <seqno> <gunseq> <source-ffid-column> <source-il-column> <source-xl-column> <source-X-column> <source-Y-column> <source-Z-column> <shot-time-column>\n");
		printf("        <receiver-file> <receiver-ffid_column> <receiver-il-column> <receiver-xl-column> <receiver-X-column> <receiver-Y-column> <receiver-Z-column>\n");
		printf("        <maximum-offset> <voxet-file> <voxet-axis-mapping> <geometry-file>\n");
		printf("        <start-day> <start-year>\n");
		printf("        Set <source-ffid-column> to zero if your input file has no source ffid column.\n");
		printf("        Set <source-il-column> and/or <source-xl-column> if input file has no column for source inline and/or crossline.\n");
		printf("        Set <receiver-ffid-column> to zero if your input file has no receiver ffid column.\n");
		printf("        Set <receiver-il-column> and/or <receiver-xl-column> if input file has no column for receiver inline and/or crossline.\n");
		printf("        Set <maximum-offset> to zero to disable maximum offset filter.\n");
		printf("        <voxet-file> must be full path to a .vo file. It can be a dummy string if <maximum-offset> is zero.\n");
		printf("        <voxet-axis-mapping> is the mapping from U-V-W to X-Y-Z. It can be a dummy string if <maximum-offset> is zero.\n");
		printf("        <start-day> is the first day of the survey, ranging from 1-366.\n");
		printf("        <start-year> is the year the survey started.\n");
		printf("\n");
		return -1;
	}

	setenv("TZ", "UTC", 1);  // TMJ - We override local timezone, force it to UTC since all field timestamps are in UTC
	const int log_level = 2;

	double max_offset = atof(argv[18]);
	if (max_offset <= 0.0)
	{
		printf("Maximum offset filter disabled, all source-receiver pairs will be included.\n");
	}
	else
	{
		printf("Maximum offset filter enabled. Offset capped at %.2lf\n",max_offset);
	}

	int start_day = atoi(argv[22]);
	int start_year = atoi(argv[23]);
	if (start_day < 1 || start_day > 366)
	{
		printf("Error! <start-day> must be >= 1 and <= 366\n");
		return -2;
	}

	struct stat sfs;
	if (stat(argv[1],&sfs) != 0)
	{
		printf("stat(%s) returned error.\n",argv[1]);
		return -2;
	}
	struct stat rfs;
	if (stat(argv[11],&rfs) != 0)
	{
		printf("stat(%s) returned error.\n",argv[11]);
		return -2;
	}

	int seqno = atoi(argv[2]);
	int gunseq = atoi(argv[3]);
	int source_ffid = atoi(argv[4]);
	int source_il = atoi(argv[5]);
	int source_xl = atoi(argv[6]);
	int source_X = atoi(argv[7]);
	int source_Y = atoi(argv[8]);
	int source_Z = atoi(argv[9]);
	int shot_time = atoi(argv[10]);
	if (source_X <= 0 || source_Y <= 0 || source_Z <= 0)
	{
		printf("Error! Source columns must be >= 1.\n");
		return -3;
	}
	if (source_ffid <= 0)
	{
		printf("No source ffid will be read.\n");
		printf("Shots will be assigned monotonically increasing index in the order they appear.\n");
	}
	if (shot_time <= 0)
	{
		printf("No shot time will be read.\n");
		printf("All shots will have zero time.\n");
	}

	int receiver_ffid = atoi(argv[12]);
	int receiver_il = atoi(argv[13]);
	int receiver_xl = atoi(argv[14]);
	int receiver_X = atoi(argv[15]);
	int receiver_Y = atoi(argv[16]);
	int receiver_Z = atoi(argv[17]);
	if (receiver_X <= 0 || receiver_Y <= 0 || receiver_Z <= 0)
	{
		printf("Error! Receiver columns must be >= 1.\n");
		return -4;
	}

	Voxet* voxet = 0L;
	Global_Coordinate_System* gcs = 0L;
	if (max_offset > 0.0)
	{
		voxet = new Voxet(log_level,argv[19]);
		gcs = voxet->Get_Global_Coordinate_System();
		gcs->Set_Transpose(argv[20]);
	}

	std::list<pos> receiver_positions = read_ascii("receiver",argv[11],receiver_ffid,0,0,0,receiver_il,receiver_xl,receiver_X,receiver_Y,receiver_Z,0,0);
	std::list<pos> source_positions = read_ascii("source",argv[1],source_ffid,seqno,gunseq,shot_time,source_il,source_xl,source_X,source_Y,source_Z,start_day,start_year);
	if (max_offset > 0.0)
	{
		int max_rx = 0;
		std::list<pos> filtered_source_positions;
		for (auto const& it1: source_positions)
		{
			double sx,sy,sz;
			gcs->Convert_Global_To_Transposed_Fractional_Index(it1.x,it1.y,it1.z,sx,sy,sz);
			sx *= gcs->Get_DX();
			sy *= gcs->Get_DY();
			sz *= gcs->Get_DZ();
			int num_included_rx = 0;
			for (auto const& it2: receiver_positions)
			{
				double rx,ry,rz;
				gcs->Convert_Global_To_Transposed_Fractional_Index(it2.x,it2.y,it2.z,rx,ry,rz);
				rx *= gcs->Get_DX();
				ry *= gcs->Get_DY();
				rz *= gcs->Get_DZ();
				double delta_x = fabs(sx-rx);
				double delta_y = fabs(sy-ry);
				if (delta_x < max_offset && delta_y < max_offset)
				{
					++num_included_rx;
				}
			}
			if (num_included_rx > 0) filtered_source_positions.push_back(it1);
			if (num_included_rx > max_rx) max_rx = num_included_rx;
		}
		printf("Maximum number of receivers within offset was %d (out of %d possible)\n",max_rx,receiver_positions.size());
		long excluded_shots = source_positions.size() - filtered_source_positions.size();
		double excluded_shots_percentage = 100.0 * (double)excluded_shots / (double)source_positions.size();
		printf("%.2lf%% (%ld) of the shots have no receiver nodes within the maximum offset.\n",excluded_shots_percentage,excluded_shots);

		vector<long> size(2);
		size[0] = filtered_source_positions.size();
		size[1] = max_rx;
		ArrND<GeomTrace> gtarr(size);
		ArrND<GeomTraceAuxiliary> gtauxarr(size);

		int isrc = 0;
		for (auto const& it1: filtered_source_positions)
		{
			printf("Shot time is %s\n",asctime(localtime(&(it1.shot_time))));
			int irec = 0;
			double sx,sy,sz;
			gcs->Convert_Global_To_Transposed_Fractional_Index(it1.x,it1.y,it1.z,sx,sy,sz);
			sx *= gcs->Get_DX();
			sy *= gcs->Get_DY();
			sz *= gcs->Get_DZ();
			for (auto const& it2: receiver_positions)
			{
				double rx,ry,rz;
				gcs->Convert_Global_To_Transposed_Fractional_Index(it2.x,it2.y,it2.z,rx,ry,rz);
				rx *= gcs->Get_DX();
				ry *= gcs->Get_DY();
				rz *= gcs->Get_DZ();
				double delta_x = fabs(sx-rx);
				double delta_y = fabs(sy-ry);
				if (delta_x < max_offset && delta_y < max_offset)
				{
					GeomTrace &gt = gtarr[isrc][irec].dat();
					gt.setSx(it1.x);
					gt.setSy(it1.y);
					gt.setSz(it1.z);
					gt.setRx(it2.x);
					gt.setRy(it2.y);
					gt.setRz(it2.z);
					gt.setLive(true);
					if (shot_time > 0)
					{
						gt.setShotTime(it1.shot_time);
						gt.setShotTimeUSec((int)(it1.shot_time_fractional_usec));
					}
					gt.setSortindex(it1.ffid);
					GeomTraceAuxiliary &gtaux = gtauxarr[isrc][irec].dat();
					gtaux.setSourceFFID(it1.ffid);
					gtaux.setSeqNo(it1.seqno);
					gtaux.setGunSeq(it1.gunseq);
					gtaux.setSourceInline(it1.il);
					gtaux.setSourceXline(it1.xl);
					gtaux.setReceiverFFID(it2.ffid);
					gtaux.setReceiverInline(it2.il);
					gtaux.setReceiverXline(it2.xl);
					++irec;
				}
			}
			for (; irec < max_rx;  ++irec)
			{
				GeomTrace &gt = gtarr[isrc][irec].dat();
				gt.setSx(0.0);
				gt.setSy(0.0);
				gt.setSz(0.0);
				gt.setRx(0.0);
				gt.setRy(0.0);
				gt.setRz(0.0);
				gt.setLive(false);
				gt.setSortindex(0);
				GeomTraceAuxiliary &gtaux = gtauxarr[isrc][irec].dat();
				gtaux.setSourceFFID(0);
				gtaux.setSeqNo(0);
				gtaux.setSourceInline(0);
				gtaux.setSourceXline(0);
				gtaux.setReceiverFFID(0);
				gtaux.setReceiverInline(0);
				gtaux.setReceiverXline(0);
			}
			++ isrc;
		}

		char* resolved_path = new char[PATH_MAX+1];
		realpath(argv[21],resolved_path);
		ArrND<GeomTrace> gtd(size, resolved_path);
		ArrND<GeomTraceAuxiliary> gtauxd(size, ((std::string)(resolved_path)+".ixl").c_str());
		delete [] resolved_path;

		gtd<<gtarr;
		gtauxd<<gtauxarr;

		vector<double> delta(2);
		delta[0] = 1.;
		delta[1] = 1.;
		vector<double> origin(2);
		origin[0] = 0.;
		origin[1] = 0.;
		vector<string> axis(2);
		axis[0] = "src";
		axis[1] = "rcv";
		vector<string> units(2);
		units[0] = "ord";
		units[1] = "ord";
		string format("GeomTrace");
		gtd.writeHed(delta,origin,axis,units,format);
		string formataux("GeomTraceAuxiliary");
		gtauxd.writeHed(delta,origin,axis,units,formataux);
	}
	else
	{
		vector<long> size(2);
		size[0] = source_positions.size();
		size[1] = receiver_positions.size();
		ArrND<GeomTrace> gtarr(size);
		ArrND<GeomTraceAuxiliary> gtauxarr(size);
		printf("Created ArrND object of dimension %d x %d.\n",source_positions.size(),receiver_positions.size());

		int isrc = 0;
		for (auto const& it1: source_positions)
		{
			int irec = 0;
			for (auto const& it2: receiver_positions)
			{
				GeomTrace &gt = gtarr[isrc][irec].dat();
				gt.setSx(it1.x);
				gt.setSy(it1.y);
				gt.setSz(it1.z);
				gt.setRx(it2.x);
				gt.setRy(it2.y);
				gt.setRz(it2.z);
				gt.setLive(true);
				if (shot_time > 0)
				{
					gt.setShotTime(it1.shot_time);
					gt.setShotTimeUSec((int)(it1.shot_time_fractional_usec));
				}
				gt.setSortindex(it1.ffid);
				GeomTraceAuxiliary &gtaux = gtauxarr[isrc][irec].dat();
				gtaux.setSourceFFID(it1.ffid);
				gtaux.setSeqNo(it1.seqno);
				gtaux.setGunSeq(it1.gunseq);
				gtaux.setSourceInline(it1.il);
				gtaux.setSourceXline(it1.xl);
				gtaux.setReceiverFFID(it2.ffid);
				gtaux.setReceiverInline(it2.il);
				gtaux.setReceiverXline(it2.xl);
				++irec;
			}
			++ isrc;
		}

		ArrND<GeomTrace> gtd(size, argv[21]);
		ArrND<GeomTraceAuxiliary> gtauxd(size, ((std::string)(argv[21])+".ixl").c_str());

		gtd<<gtarr;
		gtauxd<<gtauxarr;

		vector<double> delta(2);
		delta[0] = 1.;
		delta[1] = 1.;
		vector<double> origin(2);
		origin[0] = 0.;
		origin[1] = 0.;
		vector<string> axis(2);
		axis[0] = "src";
		axis[1] = "rcv";
		vector<string> units(2);
		units[0] = "ord";
		units[1] = "ord";
		string format("GeomTrace");
		gtd.writeHed(delta,origin,axis,units,format);
		string formataux("GeomTraceAuxiliary");
		gtauxd.writeHed(delta,origin,axis,units,formataux);
	}

	return 0;
}
