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
	printf("\nUsage: %s <ArrND-geom-file>\n\n",cmd);
}

void DumpArrND(const char* path)
{
	std::string geom_path = (std::string)path;
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

	ArrND<GeomTrace> gtd(geom_path.c_str());
	vector<long> size = gtd.size();
	int nsrc = size[0];
	int nrcv = size[1];
	vector<long> size1(1);
	size1[0]=size[1];

	if (!ixl_ok)
	{
		ArrND<GeomTrace> gtm(size1);
		for (int isrc = 0; isrc < nsrc; isrc++)
	        {
			printf("START OF SHOT\n");
			gtm << gtd[isrc];
			for (int ircv = 0;  ircv < nrcv;  ++ircv)
			{
				
				GeomTrace &gt = gtm[ircv].dat();
				if (gt.isLive())
				{
					time_t shotTime = gt.getShotTime();
					int usec = gt.getShotTimeUSec();
					char* time_str = strdup(asctime(localtime(&shotTime)));
					time_str[strlen(time_str)-1] = 0;  // remove newline char
					printf("SRC[time=%s,usec=%d,ffid=%d,pos=(%.2f,%.2f,%.2f)] - REC[pos=(%.2f,%.2f,%.2f)]\n",
							time_str,usec,gt.getSortindex(),gt.getSx(),gt.getSy(),gt.getSz(),gt.getRx(),gt.getRy(),gt.getRz());
				}
			}
			printf("END OF SHOT\n\n");
		}
	}
	else
	{
		// geom and ixl
		ArrND<GeomTrace> gtm(size1);
		ArrND<GeomTraceAuxiliary> gtauxm(size1);
		ArrND<GeomTraceAuxiliary> gtauxd(ixl_path.c_str());
		for (int isrc = 0; isrc < nsrc; isrc++)
		{
			printf("START OF SHOT\n");
			gtm << gtd[isrc];
			gtauxm << gtauxd[isrc];
			for (int ircv = 0;  ircv < nrcv;  ++ircv)
			{
				GeomTrace &gt = gtm[ircv].dat();
				GeomTraceAuxiliary &gtaux = gtauxm[ircv].dat();
				if (gt.isLive())
				{
					time_t shotTime = gt.getShotTime();
					int usec = gt.getShotTimeUSec();
					char* time_str = strdup(asctime(localtime(&shotTime)));
					time_str[strlen(time_str)-1] = 0;  // remove newline char
					printf("SRC[time=%s,usec=%d,ffid=%d,il=%d,xl=%d,seqno=%d,gunseq=%d,pos=(%.2f,%.2f,%.2f)] - REC[ffid=%d,il=%d,xl=%d,pos=(%.2f,%.2f,%.2f)]\n",
							time_str,usec,gt.getSortindex(),gtaux.getSourceInline(),gtaux.getSourceXline(),gtaux.getSeqNo(),gtaux.getGunSeq(),gt.getSx(),gt.getSy(),gt.getSz(),
							gtaux.getReceiverFFID(),gtaux.getReceiverInline(),gtaux.getReceiverXline(),gt.getRx(),gt.getRy(),gt.getRz());
				}
			}
			printf("END OF SHOT\n\n");
		}
	}
}

int main(int argc, char* argv[])
{
	if (argc == 1)
	{
		Print_Usage(argv[0]);
		return -1;
	}

	DumpArrND(argv[argc-1]);

	return 0;
}


