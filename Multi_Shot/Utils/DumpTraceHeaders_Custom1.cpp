#include <cstdio>
#include <SEGY_File.h>

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		printf("\n%s Aug-11-2016\nUsage : %s <segy-file>\n\n",argv[0],argv[0]);
		return -1;
	}

	SEGY_File* segy = new SEGY_File(argv[1],false);
	for (int iTrc = 0;  iTrc < segy->Get_Number_Of_Traces();  ++iTrc)
	{
		SEGY_Trace* trc = segy->Get_Trace(iTrc);
		SEGY_Trace_Header* hdr = trc->Get_Trace_Header();
		printf("Trace #%d/%d ::\n",iTrc+1,segy->Get_Number_Of_Traces());
		printf("   SEQ_NO = %d\n",hdr->Get_Custom1_SEQ_NO());
		printf("   GUN_SEQ = %d\n",hdr->Get_Custom1_GUN_SEQ());
		printf("   COMPON = %d\n",hdr->Get_Custom1_COMPON());
		printf("   FFID = %d\n",hdr->Get_Custom1_FFID());
		printf("   OFFSET = %d\n",hdr->Get_Custom1_OFFSET());
		printf("   RCV_ELEV = %.2f\n",(float)hdr->Get_Custom1_RCV_ELEV()/100.0f);
		printf("   DEPTH = %.2f\n",(float)hdr->Get_Custom1_DEPTH()/100.0f);
		printf("   SOU_H2OD = %.2f\n",(float)hdr->Get_Custom1_SOU_H2OD()/100.0f);
		printf("   REC_H2OD = %.2f\n",(float)hdr->Get_Custom1_REC_H2OD()/100.0f);
		printf("   AOFFSET = %.2f\n",hdr->Get_Custom1_AOFFSET());
		printf("   FLAG_VWXYZT = %s\n",hdr->Get_Custom1_FLAG_VWXYZT()==1?"Variable":"Constant");
		printf("   SOU_XD = %.2f\n",hdr->Get_Custom1_SOU_XD());
		printf("   SOU_YD = %.2f\n",hdr->Get_Custom1_SOU_YD());
		printf("   REC_XD = %.2f\n",hdr->Get_Custom1_REC_XD());
		printf("   REC_YD = %.2f\n",hdr->Get_Custom1_REC_YD());
		printf("   RCV_STAT = %d\n",hdr->Get_Custom1_RCV_STAT());
		printf("   YEAR = %d\n",hdr->Get_Custom1_YEAR());
		printf("   DAY_OF_YEAR = %d\n",hdr->Get_Custom1_DAY_OF_YEAR());
		printf("   HOUR_OF_DAY = %d\n",hdr->Get_Custom1_HOUR_OF_DAY());
		printf("   MINUTE_OF_HOUR = %d\n",hdr->Get_Custom1_MINUTE_OF_HOUR());
		printf("   SECOND_OF_MINUTE = %d\n",hdr->Get_Custom1_SECOND_OF_MINUTE());
		printf("   USEC_OF_SECOND = %d\n",hdr->Get_Custom1_USEC_OF_SECOND());
		printf("   SOU_LINE = %d\n",hdr->Get_Custom1_SOU_LINE());
		printf("   SHOT_POINT = %d\n",hdr->Get_Custom1_SHOT_POINT());
		printf("   RCV_LINE = %d\n",hdr->Get_Custom1_RCV_LINE());
		printf("   RCV_POINT = %d\n",hdr->Get_Custom1_RCV_POINT());
		printf("\n");
	}

	delete segy;
	return 0;
}
