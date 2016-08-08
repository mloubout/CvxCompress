#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <GeomTrace.h>
#include <GeomTraceAuxiliary.h>
#include <ArrND.h>
#include <SEGY_File.h>

void Print_Usage(char* cmd)
{
	printf("\nUsage: %s <input-ArrND-geom-file> <receiver-ffid> <output-ArrND-geom-file>\n\n",cmd);
}

void Extract_Receiver(const char* input_path, long rec_ffid, const char* output_path)
{
	std::string geom_path = (std::string)input_path;
	std::string ixl_path = geom_path + ".ixl";

	struct stat st_geom; 
	struct stat st_ixl;
	bool geom_ok = (stat(geom_path.c_str(),&st_geom) == 0) ? true : false;
	bool ixl_ok = (stat(ixl_path.c_str(),&st_ixl) == 0) ? true : false;

	if (!geom_ok)
	{
		printf("ERROR! Unable to open %s for reading.\n",geom_path.c_str());
		return;
	}
	if (!ixl_ok)
	{
		printf("ERROR! This utility requires a .ixl file for each .geom file.\n");
		printf("Version 08022016 and newer of MergeShotAndReceiverASCIIFiles has this feature.\n");
		return;
	}

	ArrND<GeomTrace> gtd(geom_path.c_str());
	ArrND<GeomTraceAuxiliary> gtauxd(ixl_path.c_str());
	vector<long> size = gtd.size();
	int nsrc = size[0];
	int nrcv = size[1];
	vector<long> size1(1);
	size1[0]=size[1];

	ArrND<GeomTrace> gtm(size1);
        ArrND<GeomTraceAuxiliary> gtauxm(size1);

	int num_included = 0;
	for (int isrc = 0; isrc < nsrc; isrc++)
	{
		gtauxm << gtauxd[isrc];
		bool included = false;
		for (int ircv = 0;  ircv < nrcv;  ++ircv)
		{
			GeomTraceAuxiliary &gtaux = gtauxm[ircv].dat();
			if (gtaux.getReceiverFFID() == rec_ffid) included = true;
		}
		if (included)
		{
			++num_included;
		}
	}
	printf("%d/%d shots included\n",num_included,nsrc);

	if (num_included > 0)
	{
		vector<long> size3(2);
		size3[0] = num_included;
		size3[1] = nrcv;
		ArrND<GeomTrace> gtarr(size3);
		ArrND<GeomTraceAuxiliary> gtauxarr(size3);

		int cnt = 0;
		for (int isrc = 0; isrc < nsrc; isrc++)
		{
			gtauxm << gtauxd[isrc];
			bool included = false;
			for (int ircv = 0;  ircv < nrcv;  ++ircv)
			{
				GeomTraceAuxiliary &gtaux = gtauxm[ircv].dat();
				if (gtaux.getReceiverFFID() == rec_ffid) included = true;
			}
			if (included)
			{
				printf("isrc %d is included\n",isrc);
				gtm << gtd[isrc];
				for (int ircv = 0;  ircv < nrcv;  ++ircv)
				{
					GeomTrace &gt = gtm[ircv].dat();
					GeomTraceAuxiliary &gtaux = gtauxm[ircv].dat();
					GeomTrace &gto = gtarr[cnt][ircv].dat();
					GeomTraceAuxiliary &gtauxo = gtauxarr[cnt][ircv].dat();
					gto = gt;
					gtauxo = gtaux;
				}
				/*
				printf("1\n");
				gtarr[cnt] << gtm;
				printf("2\n");
				gtauxarr[cnt] << gtauxm;
				printf("3\n");
				*/
				++cnt;
				printf("%d/%d\n",cnt,num_included);
			}
		}

		ArrND<GeomTrace> gtd2(size3, output_path);
		ArrND<GeomTraceAuxiliary> gtauxd2(size3, ((std::string)(output_path)+".ixl").c_str());

		gtd2<<gtarr;
		gtauxd2<<gtauxarr;

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
		gtd2.writeHed(delta,origin,axis,units,format);
		string formataux("GeomTraceAuxiliary");
		gtauxd2.writeHed(delta,origin,axis,units,formataux);		
	}
}

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		Print_Usage(argv[0]);
		return -1;
	}

	Extract_Receiver(argv[1],atol(argv[2]),argv[3]);

	return 0;
}


