#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <GeomTrace.h>
#include <GeomTraceAuxiliary.h>
#include <ArrND.h>
#include <SEGY_File.h>

void Print_Usage(char* cmd)
{
	printf("\n%s Aug-11-2016\nUsage: %s in=<input-geom-file> rcv_stat=<receiver-ffid> sou_line=<min>:<max> shot_po=<min>:<max> out=<output-geom-file>\n\n",cmd,cmd);
}

std::string input;
std::string output;
bool use_rcv_stat = false;
int rcv_stat = 0;
int sou_line_min = 0;
int sou_line_max = -1;
int shot_po_min = 0;
int shot_po_max = -1;

bool Process_Arguments(int argc, char* argv[])
{
	if (argc > 1)
	{
		int num1, num2;
		char str1[256];
		for (int i = 1;  i < argc;  ++i)
		{
			if (sscanf(argv[i], "in=%s", str1) == 1)
			{
				input = (std::string)str1;
			}
			if (sscanf(argv[i], "out=%s", str1) == 1)
			{
				output = (std::string)str1;
			}
			if (sscanf(argv[i], "rcv_stat=%d", &num1) == 1)
			{
				rcv_stat = num1;
				use_rcv_stat = true;
			}
			if (sscanf(argv[i], "sou_line=%d:%d", &num1, &num2) == 2)
			{
				sou_line_min = num1;
				sou_line_max = num2;
			}
			if (sscanf(argv[i], "shot_po=%d:%d", num1, &num2) == 2)
			{
				shot_po_min = num1;
				shot_po_max = num2;
			}
		}

		printf("Arguments :: ");
		if (input.size() > 0) printf("in=%s ",input.c_str());
		if (output.size() > 0) printf("out=%s ",output.c_str());
		if (use_rcv_stat) printf("rcv_stat=%d ",rcv_stat);
		if (sou_line_min <= sou_line_max) printf("sou_line=%d:%d ",sou_line_min,sou_line_max);
		if (shot_po_min <= shot_po_max) printf("shot_po=%d:%d ",shot_po_min,shot_po_max);
		printf("\n");
	}

	if (input.size() > 0 && output.size() > 0)
		return true;
	else
		return false;
}

bool Is_Included(GeomTraceAuxiliary &gtaux)
{
	if (
		(!use_rcv_stat || rcv_stat == gtaux.getReceiverFFID()) && 
		((sou_line_min > sou_line_max) || (gtaux.getSourceInline() >= sou_line_min && gtaux.getSourceInline() <= sou_line_max)) &&
		((shot_po_min > shot_po_max) || (gtaux.getSourceXline() >= shot_po_min && gtaux.getSourceXline() <= shot_po_max))
		)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Extract_Receiver()
{
	std::string geom_path = (std::string)input;
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
		gtm << gtd[isrc];
		gtauxm << gtauxd[isrc];
		bool included = false;
		for (int ircv = 0;  ircv < nrcv;  ++ircv)
		{
			GeomTrace &gt = gtm[ircv].dat();
			GeomTraceAuxiliary &gtaux = gtauxm[ircv].dat();
			if (gt.isLive())
			{
				//printf("SRC[il=%d,xl=%d,seqno=%d,gunseq=%d] - REC[ffid=%d,il=%d,xl=%d]\n",
				//		gtaux.getSourceInline(),gtaux.getSourceXline(),gtaux.getSeqNo(),gtaux.getGunSeq(),
				//		gtaux.getReceiverFFID(),gtaux.getReceiverInline(),gtaux.getReceiverXline());
				included = included || Is_Included(gtaux);
			}
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
			gtm << gtd[isrc];
			gtauxm << gtauxd[isrc];
			bool included = false;
			for (int ircv = 0;  ircv < nrcv;  ++ircv)
			{
				GeomTrace &gt = gtm[ircv].dat();
				GeomTraceAuxiliary &gtaux = gtauxm[ircv].dat();
				included = included || (gt.isLive() && Is_Included(gtaux));
			}
			if (included)
			{
				printf("isrc %d is included\n",isrc);
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

		char* str = new char[PATH_MAX+1];
		realpath(output.c_str(),str);
		output = (string)str;
		delete [] str;
		ArrND<GeomTrace> gtd2(size3, output.c_str());
		ArrND<GeomTraceAuxiliary> gtauxd2(size3, (output+".ixl").c_str());

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
	if (!Process_Arguments(argc,argv))
	{
		Print_Usage(argv[0]);
		return -1;
	}

	Extract_Receiver();

	return 0;
}


