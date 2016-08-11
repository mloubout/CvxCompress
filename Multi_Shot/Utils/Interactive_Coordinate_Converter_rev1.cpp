#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <GeomTrace.h>
#include <ArrND.h>

#include <swapbytes.h>
#include <Voxet.hxx>
#include <Global_Coordinate_System.hxx>

// 
// DON'T EDIT THIS FILE ANYMORE, IT HAS BEEN FROZEN
// WORKING COPY IS *_work
// 

int Parse_Command(const char* cmd)
{
	if (strcmp(cmd, "g2l") == 0)
	{
		return 1;
	}
	else if (strcmp(cmd, "g2i") == 0)
        {
                return 2;
        }
	else if (strcmp(cmd, "l2g") == 0)
        {
                return 3;
        }
	else if (strcmp(cmd, "l2i") == 0)
        {
                return 4;
        }
	else if (strcmp(cmd, "i2g") == 0)
        {
                return 5;
        }
	else if (strcmp(cmd, "i2l") == 0)
        {
                return 6;
        }
	else if (strcmp(cmd, "survey") == 0)
	{
		return 7;
	}
	else if (strcmp(cmd, "depth") == 0)
	{
		return 8;
	}
	else
	{
		return -1;
	}
}

void split_string(const char* str, char**& tokens, int& ntokens)
{
	// split string into tokens separated by whitespace
	int len = strlen(str);
	if (len > 0)
	{
		ntokens = 1;
		tokens = new char*[1024];  // yes it will blow up if string has > 1024 tokens.
		char* ptr = new char[len];
		ptr[0] = 0;
		tokens[0] = ptr;
		bool copy = true;
		int j = 0;
		for (int i = 0;  i < len;  ++i)
		{
			if (!isspace(str[i]))
			{
				if (i > 0 && isspace(str[i-1]))
				{
					ptr[j++] = 0;
					tokens[ntokens++] = ptr + j;
				}
				ptr[j++] = str[i];
			}
		}
		ptr[j++] = 0;
	}
	else
	{
		tokens = 0L;
		ntokens = 0;
	}
}

int choose(
	const char* question,
	const char* choice1,
	const char* choice2
	)
{
	int choice = -1;
	char str[4096];
	while (choice < 0)
	{
		printf("%s (%s,%s) ? ",question,choice1,choice2);
		gets(str);
		int len = strlen(str);
		if (len > 0)
		{
			if (len <= strlen(choice1) && strncmp(str,choice1,len) == 0)
			{
				return 0;
			}
			else if (len <= strlen(choice2) && strncmp(str,choice2,len) == 0)
                        {
                                return 1;
                        }
		}
	}
	return -1;  // should never happen
}

int choose(
	const char* question,
	const char* choice1,
	const char* choice2,
	const char* choice3
	)
{
	int choice = -1;
	char str[4096];
	while (choice < 0)
	{
		printf("%s (%s,%s,%s) ? ",question,choice1,choice2,choice3);
		gets(str);
		int len = strlen(str);
		if (len > 0)
		{
			if (len <= strlen(choice1) && strncmp(str,choice1,len) == 0)
			{
				return 0;
			}
			else if (len <= strlen(choice2) && strncmp(str,choice2,len) == 0)
                        {
                                return 1;
                        }
			else if (len <= strlen(choice3) && strncmp(str,choice3,len) == 0)
                        {
                                return 2;
                        }
		}
	}
	return -1; // should never happen
}

int choose_token(
	const char* question,
	char** tokens,
	int ntokens,
	bool allow_skip,
	int* chosen_ones
	)
{
	int opto_count = 0;
	char str[4096];
	while (true)
	{
		for (int i = 0;  i < ntokens;  ++i)
		{
			bool found_it = false;
			for (int j = 0;  j < ntokens && !found_it && chosen_ones[j] != 0;  ++j) if (chosen_ones[j] == i+1) found_it = true;
			if (!found_it) printf("%d : %s\n",i+1,tokens[i]);
		}
		char buf[128];
		sprintf(buf,", %d more ENTER to skip",3-opto_count);
		printf("%s%s (type number in front of your choice%s) ? ",allow_skip?"(optional) ":"",question,allow_skip&&opto_count>0&&opto_count<3?buf:"");
		gets(str);
		if (strlen(str) == 0)
		{
			if (allow_skip) ++opto_count;
			if (opto_count >= 3) return 0;
		}
		else
		{
			opto_count = 0;
			int chosen = atoi(str);
			if (allow_skip && chosen < 0) return 0;
			if (chosen >= 1 && chosen <= ntokens)
			{
				// valid value
				for (int j = 0;  j < ntokens;  ++j)
				{
					if (chosen_ones[j] == 0)
					{
						chosen_ones[j] = chosen;
						break;
					}
				}
				return chosen;
			}
		}
	}
}

void strip_space(char* str)
{
	int len = strlen(str);
	if (len > 0)
	{
		int i0 = 0;
		while (i0 < len && isspace(str[i0])) ++i0;
		if (i0 >= len)
		{
			str[0] = 0;
		}
		else
		{
			int i1 = len-1;
			while (i1 > i0 && isspace(str[i1])) --i1;
			for (int i = i0;  i <= i1;  ++i) str[i-i0] = str[i];
			str[i1-i0+1] = 0;
		}
	}
}

double get_double_value(
	const char* question,
	double min_val,
	double max_val
	)
{
	while (true)
	{
		char buf[128];
		if (min_val > max_val)
			buf[0] = 0;
		else
			sprintf(buf, " (min=%.2f, max=%.2f)", min_val, max_val);
		printf("%s%s ? ",question,buf);
		char str[4096];
		gets(str);
		strip_space(str);
		char* endptr;
		double val = strtod(str, &endptr);
		if (str[0] != 0 && endptr[0] == 0)
		{
			// entire string is a valid double
			if (min_val > max_val || (val >= min_val && val <= max_val)) return val;
		}
	};
}

double get_double_value(
	const char* question
	)
{
	return get_double_value(question,1,0);
}

long get_integer_value(
	const char* question,
	long min_val,
	long max_val
	)
{
	while (true)
	{
		char buf[128];
		if (min_val > max_val)
			buf[0] = 0;
		else
			sprintf(buf, " (min=%ld, max=%ld)", min_val, max_val);
		printf("%s%s ? ",question,buf);
		char str[4096];
		gets(str);
		strip_space(str);
		char* endptr;
		long val = strtol(str, &endptr, 0);
		if (str[0] != 0 && endptr[0] == 0)
		{
			// entire string is a valid long
			if (min_val > max_val || (val >= min_val && val <= max_val)) return val;
		}
	};
}

long get_integer_value(
	const char* question
	)
{
	return get_integer_value(question,1,0);
}

void get_output_filename(
	const char* question,
	char* filename
	)
{
	filename[0] = 0;
	while(true)
	{
		printf("%s ? ",question);
		gets(filename);
		FILE* fp = fopen(filename, "w");
		if (fp != 0L)
		{
			fclose(fp);
			return;
		}
	};
}

void get_input_filename(
	const char* question,
	char* filename
	)
{
	filename[0] = 0;
	while(true)
	{
		printf("%s ? ",question);
		gets(filename);
		FILE* fp = fopen(filename, "r");
		if (fp != 0L)
		{
			fclose(fp);
			return;
		}
	};
}

int get_next_line_from_file_strip_comments(FILE* fp, char* str)
{
	while (true)
	{
		if (feof(fp) != 0) return 0;
		char* s = fgets(str,4096,fp);
		int len = s != 0L ? strlen(s) : 0;
		if (len > 0)
		{
			// strip end of line comment starting with #
			int i = 0;
			for (; i < len && s[i] != '#';  ++i);
			if (i < len)
			{
				printf("Comment ignored: %s",s+i);  // comment is terminated with newline char, don't insert extra one.
				s[i] = 0;
				// update length
				len = strlen(s);
			}

			if (len > 0)
			{
				// strip leading and trailing whitespace
				strip_space(s);
				// update length
				len = strlen(s);
			}
		}
		if (len > 0) return len;
	};
	return 0;
}

int read_source_locations(
		FILE* fp,
		Global_Coordinate_System* gcs,
		int num_expected_tokens,
		int is_global,
		int idx0,
		int idx1,
		int idx2,
		int idxFFID,
		double* srcx,
		double* srcy,
		double* srcz,
		int* FFID
		)
{
	char** tokens;
	int ntokens;
	char str[4096];

	int num_sources = 0;
	rewind(fp);
	bool done = false;
	while (!done)
	{
		str[0] = 0;
		int len = get_next_line_from_file_strip_comments(fp,str);
		if (len <= 0)
		{
			done = true;
		}
		else
		{
			split_string(str, tokens, ntokens);
			double c0 = idx0 > 0 ? atof(tokens[idx0-1]) : 0.0;
			double c1 = idx1 > 0 ? atof(tokens[idx1-1]) : 0.0;
			double c2 = idx2 > 0 ? atof(tokens[idx2-1]) : 0.0;
			int currFFID = (idxFFID > 0) ? (int)round(atof(tokens[idxFFID-1])) : 0;
			if (is_global)
			{
				// convert to local before storing shot locations
				double d0,d1,d2;
				gcs->Convert_Global_To_Transposed_Fractional_Index(c0,c1,c2,d0,d1,d2);
				c0 = d0 * gcs->Get_DX();
				c1 = d1 * gcs->Get_DY();
				c2 = d2 * gcs->Get_DZ();
			}
			if (srcx != 0L) srcx[num_sources] = c0;
			if (srcy != 0L) srcy[num_sources] = c1;
			if (srcz != 0L) srcz[num_sources] = c2;
			if (FFID != 0L) FFID[num_sources] = currFFID;
			++num_sources;
		}
		if (ntokens > 0)
		{
			delete [] tokens[0];
			delete [] tokens;
			tokens = 0L;
			ntokens = 0;
		}
	};
	return num_sources;
}

int Read_Positions_From_ASCII_File(
	Global_Coordinate_System* gcs,
	double*& srcx,
	double*& srcy,
	double*& srcz,
	int*& FFID
	)
{
	int num_src = 0;
	char filename[4096];
	get_input_filename("Please provide full path to file containing source locations: ",filename);
	FILE* fp = fopen(filename, "r");
	if (fp != 0L)
	{
		char str[4096];
		int len = get_next_line_from_file_strip_comments(fp,str);
		if (len > 0)
		{
			char** tokens;
			int ntokens;
			split_string(str, tokens, ntokens);
			if (ntokens > 0)
			{
				for (int i = 0;  i < ntokens;  ++i) if (i == 0) printf("%d:%s",i+1,tokens[i]); else printf(" %d:%s",i+1,tokens[i]);
				printf("\n");

				int* chosen_ones = new int[ntokens];
				for (int i = 0;  i < ntokens;  ++i) chosen_ones[i] = 0;
				int is_global = choose("These coordinates are","local","global");
				int idx0=-1, idx1=-1, idx2=-1, idxFFID=-1;
				int idxZ_valid = 0;
				if (is_global)
				{
					double dx, dy, dz;
					gcs->Convert_Global_To_Transposed_Fractional_Index(1.0,0.0,0.0,dx,dy,dz);
					int u_is_z = dz != 0.0 ? 1 : 0;
					gcs->Convert_Global_To_Transposed_Fractional_Index(0.0,1.0,0.0,dx,dy,dz);
					int v_is_z = dz != 0.0 ? 1 : 0;
					gcs->Convert_Global_To_Transposed_Fractional_Index(0.0,0.0,1.0,dx,dy,dz);
					int w_is_z = dz != 0.0 ? 1 : 0;

					idx0 = choose_token("Which column is the U coordinate",tokens,ntokens,u_is_z,chosen_ones);
					idx1 = choose_token("Which column is the V coordinate",tokens,ntokens,v_is_z,chosen_ones);
					idx2 = choose_token("Which column is the W coordinate",tokens,ntokens,w_is_z,chosen_ones);
					idxFFID = choose_token("Which column is the FFID",tokens,ntokens,true,chosen_ones);
					idxZ_valid = (u_is_z && idx0 > 0) || (v_is_z && idx1 > 0) || (w_is_z && idx2 > 0) ? 1 : 0;
				}
				else
				{
					idx0 = choose_token("Which column is the X coordinate",tokens,ntokens,false,chosen_ones);
					idx1 = choose_token("Which column is the Y coordinate",tokens,ntokens,false,chosen_ones);
					idx2 = choose_token("Which column is the Z coordinate",tokens,ntokens,true,chosen_ones);
					idxFFID = choose_token("Which column is the FFID",tokens,ntokens,true,chosen_ones);
					idxZ_valid = idx2 > 0 ? 1 : 0;
				}
				delete [] chosen_ones;
				delete [] tokens[0];
				delete [] tokens;

				num_src = read_source_locations(fp,gcs,ntokens,is_global,idx0,idx1,idx2,idxFFID,srcx,srcy,srcz,FFID);  // count number of source locations
				srcx = new double[num_src];
				srcy = new double[num_src];
				srcz = new double[num_src];
				FFID = new int[num_src];
				num_src = read_source_locations(fp,gcs,ntokens,is_global,idx0,idx1,idx2,idxFFID,srcx,srcy,srcz,FFID);  // store source locations
				if (idxFFID <= 0)
				{
					// no FFID was provided, put 1..num_src instead
					printf("Warning! Source number will be used in place of FFID since no FFID was provided.\n");
					for (int i = 0;  i < num_src;  ++i) FFID[i] = i+1;
				}

				int what_depth_to_use = 0;
				if (idxZ_valid)
				{
					what_depth_to_use = choose("Do you want original source depth from file or something else","original","constant value","sea floor");
				}
				else
				{
					what_depth_to_use = choose("What source depth do you want","constant value","sea floor") + 1;
				}
				if (what_depth_to_use == 1)
				{
					// constant value
					double source_depth = get_double_value("What constant source depth do you want",gcs->Get_DZ(),gcs->Get_DZ()*(double)(gcs->Get_NZ()-1));
					for (int i = 0;  i < num_src;  ++i) srcz[i] = source_depth;
				}
				else if (what_depth_to_use == 2)
				{
					// sea floor
				}
			}
		}
		fclose(fp);
	}
	return num_src;
}

void Add_Depth(Global_Coordinate_System* gcs)
{
	char str[4096];
	get_input_filename("Please provide full path to file containing source locations: ",str);
	
}

void Create_Survey(Global_Coordinate_System* gcs)
{
	char str[4096];

	int read_source_from_file = !choose("Do you want to read source locations from a file","yes","no");

	int num_src = 0, num_src_recpatch = 0;
	double *srcx = 0L, *srcy = 0L, *srcz = 0L;
	double *srcx_recpatch = 0L, *srcy_recpatch = 0L, *srcz_recpatch = 0L;
	int *FFID = 0L, *FFID_recpatch = 0L;
	if (read_source_from_file)
	{
		num_src = Read_Positions_From_ASCII_File(gcs,srcx,srcy,srcz,FFID);
		int recpatch_separate = !choose("Do you want to generate receiver patch from a different set of source locations","yes","no");
		if (recpatch_separate)
		{
			num_src_recpatch = Read_Positions_From_ASCII_File(gcs,srcx_recpatch,srcy_recpatch,srcz_recpatch,FFID_recpatch);
			if (num_src_recpatch != num_src)
			{
				printf("Error! The two source location files have different number of locations (%d != %d)\n",num_src,num_src_recpatch);
				return;
			}
			else
			{
				for (int i = 0;  i < num_src;  ++i)
				{
					if (FFID[i] != FFID_recpatch[i])
					{
						printf("Error! FFID's don't match in the two source location files.\n");
						return;
					}
				}
			}
		}
		else
		{
			num_src_recpatch = num_src;
			srcx_recpatch = srcx;
			srcy_recpatch = srcy;
			srcz_recpatch = srcz;
			FFID_recpatch = FFID;
		}
	}
	else
	{
		int shape = 0;	
		const char* circular_str = "circular";
		const char* oval_str = "oval";
		const char* rectangular_str = "rectangular";
		while (shape == 0)
		{
			printf("What is the survey shape (%s,%s,%s)? ",circular_str,oval_str,rectangular_str);
			gets(str);
			int len = strlen(str);
			if (len > 0)
			{
				if (len <= strlen(circular_str) && strncmp(str,circular_str,strlen(str)) == 0)
				{
					shape = 1; // circular
				}
				else if (len <= strlen(oval_str) && strncmp(str,oval_str,strlen(str)) == 0)
				{
					shape = 2; // oval
				}
				else if (len <= strlen(rectangular_str) && strncmp(str,rectangular_str,strlen(str)) == 0)
				{
					shape = 3; // rectangular
				}
			}
		}

		double center_off_x, center_off_y;
		printf("Please specify center point for survey shape:\n");
		char* end_ptr = 0L;
		do
		{
			printf("  What is the distance from center of volume along inline axis? ");
			gets(str);
			end_ptr = 0L;
			center_off_x = strtod(str, &end_ptr);
		} while (end_ptr == 0L || end_ptr == str);
		do
		{
			printf("  What is the distance from center of volume along xline axis? ");
			gets(str);
			end_ptr = 0L;
			center_off_y = strtod(str, &end_ptr);
		} while (end_ptr == 0L || end_ptr == str);

		double half_ix = (double)gcs->Get_NX() / 2.0;
		double half_iy = (double)gcs->Get_NY() / 2.0;
		double center_ix = half_ix + center_off_x / gcs->Get_DX();  // center point as a fractional local index
		double center_iy = half_iy + center_off_y / gcs->Get_DY();
		double center_u, center_v, center_w;
		gcs->Convert_Transposed_Fractional_Index_To_Global(center_ix,center_iy,0.0,center_u,center_v,center_w);
		printf("Center point is U=%lf,V=%lf,W=%lf\n",center_u,center_v,center_w);

		double iline_width=0.0, xline_width=0.0;
		if (shape == 1)
		{
			// circular
			while (iline_width <= 0.0)
			{
				printf("Enter diameter of circle? ");
				gets(str);
				double val = atof(str);
				if (val > 0.0)
				{
					iline_width = val;
					xline_width = val;
				}
			};
		}
		else if (shape == 2)
		{
			// oval
			while (iline_width <= 0.0)
			{
				printf("Enter diameter of oval along inline dimension? ");
				gets(str);
				double val = atof(str);
				if (val > 0.0)
				{
					iline_width = val;
				}
			};
			while (xline_width <= 0.0)
			{
				printf("Enter diameter of oval along xline dimension? ");
				gets(str);
				double val = atof(str);
				if (val > 0.0)
				{
					xline_width = val;
				}
			};
		}
		else if (shape == 3)
		{
			// rectangle
			while (iline_width <= 0.0)
			{
				printf("Enter width of rectangle along inline dimension? ");
				gets(str);
				double val = atof(str);
				if (val > 0.0)
				{
					iline_width = val;
				}
			};
			while (xline_width <= 0.0)
			{
				printf("Enter width of rectangle along xline dimension? ");
				gets(str);
				double val = atof(str);
				if (val > 0.0)
				{
					xline_width = val;
				}
			};
		}
	}

	// print the source locations for debug
	FILE* fp2 = fopen("sources.txt", "w");
	for (int iSrc = 0;  iSrc < num_src;  ++iSrc)
	{
		//printf("FFID=%d, srcx=%.2f, srcy=%.2f, srcz=%.2f\n",FFID!=0L?FFID[iSrc]:0,srcx[iSrc],srcy[iSrc],srcz[iSrc]);

		double d0, d1, d2;
		gcs->Convert_Transposed_Fractional_Index_To_Global(srcx[iSrc]/gcs->Get_DX(),srcy[iSrc]/gcs->Get_DY(),srcz[iSrc]/gcs->Get_DZ(),d0,d1,d2);
		fprintf(fp2,"%.6f %.6f %.6f\n",d0,d1,d2);
		if (iSrc > 0 && FFID[iSrc] > FFID[iSrc-1]+1) fprintf(fp2,"\n");
	}
	fclose(fp2);

	// do the receivers now
	printf("\nGenerated %d shots / source locations.\nLet's add receivers.\n\n",num_src);

	int receiver_patch_type = choose("What receiver shape do you want","streamer","centered patch");
	if (receiver_patch_type == 0) // streamer
	{
		long num_streamers = get_integer_value("How many streamers",1,100);
		double streamer_distance = 0.0;
		if (num_streamers > 1)
		{
			streamer_distance = get_double_value("What's the distance between streamers",gcs->Get_DY(),gcs->Get_DY()*(double)(gcs->Get_NY())/2.0);
		}
		double near_trace_gap = get_double_value("The sign of the near trace gap determines whether streamer are behind or ahead of the source looking towards positive X.\nWhat's the near trace gap (horizontal distance from source to nearest receiver)",-gcs->Get_DX()*(double)(gcs->Get_NX())/2.0,gcs->Get_DX()*(double)(gcs->Get_NX())/2.0);
		double GI_sign = near_trace_gap < 0.0 ? -1.0 : 1.0;
		double streamer_length = get_double_value("What's the streamer length",gcs->Get_DX(),gcs->Get_DX()*(double)(gcs->Get_NX())/2.0);
		double streamer_group_interval = get_double_value("What's the group interval",gcs->Get_DX(),gcs->Get_DX()*(double)(gcs->Get_NX())/2.0);
		long nchan_per_streamer = (long)round(streamer_length / streamer_group_interval);
		printf("Each streamer will have %ld channels per component.\n",nchan_per_streamer);
		int clip_rx = choose("Do you want receiver patch to be clipped to model","yes","no") ? 0 : 1;
                double rx_depth = get_double_value("What is the receiver depth",0,gcs->Get_DZ()*(double)gcs->Get_NZ());
		
		// determine maximum number of receivers
		int max_num_rx = 0;
		for (int iSrc = 0;  iSrc < num_src;  ++iSrc)
		{
			int num_rx = 0;
			for (long iStreamer = 0;  iStreamer < num_streamers;  ++iStreamer)
			{
				double y_offset = streamer_distance * (double)iStreamer - (streamer_distance * (double)(num_streamers - 1) / 2.0);
				for (int iChan = 0;  iChan < nchan_per_streamer;  ++iChan)
				{
					double x_offset = GI_sign * streamer_group_interval * (double)iChan + near_trace_gap;
					double lrx = x_offset + srcx[iSrc];
					double lry = y_offset + srcy[iSrc];
					if ( !clip_rx || (
								lrx >= 0.0 && lrx <= gcs->Get_DX() * (double)gcs->Get_NX() &&
								lry >= 0.0 && lry <= gcs->Get_DY() * (double)gcs->Get_NY() 
							 )
					   )
					{
						++num_rx;
					}
				}
			}
			if (num_rx > max_num_rx) max_num_rx = num_rx;
		}
		printf("max_num_rx = %d\n",max_num_rx);
		
		char filename[4096];
		get_output_filename("Please enter full path for geometry file",filename);
		char confirm[4096];
		sprintf(confirm, "Geometry file will be written to %s. Do you want to continue", filename);
		int confirmed = choose(confirm, "yes", "no") == 1 ? 0 : 1;
		
		if (confirmed)
		{
			vector<long> size(2);
			size[0] = num_src;
			size[1] = max_num_rx;
			ArrND<GeomTrace> gtarr(size);

			for (int iSrc = 0;  iSrc < num_src;  ++iSrc)
			{
				int iRcv = 0;
				for (long iStreamer = 0;  iStreamer < num_streamers;  ++iStreamer)
				{
					double y_offset = streamer_distance * (double)iStreamer - (streamer_distance * (double)(num_streamers - 1) / 2.0);
					for (int iChan = 0;  iChan < nchan_per_streamer;  ++iChan)
					{
						double x_offset = GI_sign * streamer_group_interval * (double)iChan + near_trace_gap;
						double lrx = x_offset + srcx[iSrc];
						double lry = y_offset + srcy[iSrc];
						if ( !clip_rx || (
									lrx >= 0.0 && lrx <= gcs->Get_DX() * (double)gcs->Get_NX() &&
									lry >= 0.0 && lry <= gcs->Get_DY() * (double)gcs->Get_NY() 
								 )
						   )
						{
							// convert to global coordinates
							double sd0, sd1, sd2;
							gcs->Convert_Transposed_Fractional_Index_To_Global(srcx[iSrc]/gcs->Get_DX(),srcy[iSrc]/gcs->Get_DY(),srcz[iSrc]/gcs->Get_DZ(),sd0,sd1,sd2);
							double rd0, rd1, rd2;
							gcs->Convert_Transposed_Fractional_Index_To_Global(lrx/gcs->Get_DX(),lry/gcs->Get_DY(),rx_depth/gcs->Get_DZ(),rd0,rd1,rd2);

							GeomTrace &gt = gtarr[iSrc][iRcv].dat();
							gt.setSx(sd0);
							gt.setSy(sd1);
							gt.setSz(sd2);
							gt.setRx(rd0);
							gt.setRy(rd1);
							gt.setRz(rd2);
							gt.setLive(true);
							gt.setSortindex(FFID[iSrc]);
							++iRcv;
						}
					}
				}
				for (;  iRcv < max_num_rx;  ++iRcv)
				{
					GeomTrace &gt = gtarr[iSrc][iRcv].dat();
					gt.setSx(0.0);
					gt.setSy(0.0);
					gt.setSz(0.0);
					gt.setRx(0.0);
					gt.setRy(0.0);
					gt.setRz(0.0);
					gt.setLive(false);
					gt.setSortindex(0);
				}
			}

			ArrND<GeomTrace> gtd(size, filename);

			gtd<<gtarr;

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
		}
	}
	else if (receiver_patch_type == 1) // centered patch
	{
		double iline_width = get_double_value("What is size of patch along inline",0.0,gcs->Get_DX()*(double)gcs->Get_NX());
		double iline_delta = get_double_value("What is the receiver spacing along inline",gcs->Get_DX()/10.0,gcs->Get_DX()*(double)gcs->Get_NX());
		double xline_width = get_double_value("What is size of patch along xline",0.0,gcs->Get_DY()*(double)gcs->Get_NY());
		double xline_delta = get_double_value("What is the receiver spacing along xline",gcs->Get_DY()/10.0,gcs->Get_DY()*(double)gcs->Get_NY());
		int clip_rx = choose("Do you want receiver patch to be clipped to model","yes","no") ? 0 : 1;
		double rx_depth = get_double_value("What is the receiver depth",0,gcs->Get_DZ()*(double)gcs->Get_NZ());
		
		// determine maximum number of receivers
		int max_num_rx = 0;
		for (int iSrc = 0;  iSrc < num_src_recpatch;  ++iSrc)
		{
			int num_rx = 0;
			double min_recx = srcx_recpatch[iSrc] - (iline_width / 2.0);
			double max_recx = min_recx + iline_width;
			double min_recy = srcy_recpatch[iSrc] - (xline_width / 2.0);
			double max_recy = min_recy + xline_width;
			int min_irecx = (int)floor(min_recx / iline_delta);
			int max_irecx = (int)ceil(max_recx / iline_delta);
			int min_irecy = (int)floor(min_recy / xline_delta);
			int max_irecy = (int)ceil(max_recy / xline_delta);
			for (int iy = min_irecy;  iy <= max_irecy;  ++iy)
			{
				for (int ix = min_irecx;  ix <= max_irecx;  ++ix)
				{
					double lrx = (double)ix * iline_delta;
					double lry = (double)iy * xline_delta;
					if (lrx >= min_recx && lrx <= max_recx && lry >= min_recy && lry <= max_recy)
					{
						if (!clip_rx || (lrx >= 0.0 && lrx <= gcs->Get_DX() * (double)gcs->Get_NX() && lry >= 0.0 && lry <= gcs->Get_DY() * (double)gcs->Get_NY()))
						{
							++num_rx;
						}
					}
				}
			}
			if (num_rx > max_num_rx) max_num_rx = num_rx;
		}
		printf("max_num_rx = %d\n",max_num_rx);

		/* OLD CODE - truly centered patch
		int nrx_2 = (int)trunc(iline_width/(2.0*iline_delta));
		int nry_2 = (int)trunc(xline_width/(2.0*xline_delta));
		int max_num_rx = (2 * nrx_2 + 1) * (2 * nry_2 + 1);
		if (clip_rx)
		{
			max_num_rx = 0;
			for (int iSrc = 0;  iSrc < num_src;  ++iSrc)
			{
				int num_rx = 0;
				for (int iy = -nry_2;  iy <= nry_2;  ++iy)
				{
					for (int ix = -nrx_2;  ix <= nrx_2;  ++ix)
					{
						double lrx = srcx[iSrc] + (double)ix * iline_delta;
						double lry = srcy[iSrc] + (double)iy * xline_delta;
						if (
								lrx >= 0.0 && lrx <= gcs->Get_DX() * (double)gcs->Get_NX() &&
								lry >= 0.0 && lry <= gcs->Get_DY() * (double)gcs->Get_NY() 
						   )
						{
							++num_rx;
						}
					}
				}
				if (num_rx > max_num_rx) max_num_rx = num_rx;
				//printf("FFID=%d, srcx=%.2f, srcy=%.2f, srcz=%.2f, num_rx=%d\n",FFID!=0L?FFID[iSrc]:0,srcx[iSrc],srcy[iSrc],srcz[iSrc],num_rx);
			}
		}
		printf("max_num_rx = %d\n",max_num_rx);
		*/
		
		char filename[4096];
		get_output_filename("Please enter full path for geometry file",filename);
		char confirm[4096];
		sprintf(confirm, "Geometry file will be written to %s. Do you want to continue", filename);
		int confirmed = choose(confirm, "yes", "no") == 1 ? 0 : 1;
		
		if (confirmed)
		{
			vector<long> size(2);
			size[0] = num_src;
			size[1] = max_num_rx;
			ArrND<GeomTrace> gtarr(size);

			for (int iSrc = 0;  iSrc < num_src;  ++iSrc)
			{
				int iRcv = 0;
				double min_recx = srcx_recpatch[iSrc] - (iline_width / 2.0);
				double max_recx = min_recx + iline_width;
				double min_recy = srcy_recpatch[iSrc] - (xline_width / 2.0);
				double max_recy = min_recy + xline_width;
				int min_irecx = (int)floor(min_recx / iline_delta);
				int max_irecx = (int)ceil(max_recx / iline_delta);
				int min_irecy = (int)floor(min_recy / xline_delta);
				int max_irecy = (int)ceil(max_recy / xline_delta);
				for (int iy = min_irecy;  iy <= max_irecy;  ++iy)
				{
					for (int ix = min_irecx;  ix <= max_irecx;  ++ix)
					{
						double lrx = (double)ix * iline_delta;
						double lry = (double)iy * xline_delta;
						if (lrx >= min_recx && lrx <= max_recx && lry >= min_recy && lry <= max_recy)
						{
							if (!clip_rx || (lrx >= 0.0 && lrx <= gcs->Get_DX() * (double)gcs->Get_NX() && lry >= 0.0 && lry <= gcs->Get_DY() * (double)gcs->Get_NY()))
							{
								// convert to global coordinates
								double sd0, sd1, sd2;
								gcs->Convert_Transposed_Fractional_Index_To_Global(srcx[iSrc]/gcs->Get_DX(),srcy[iSrc]/gcs->Get_DY(),srcz[iSrc]/gcs->Get_DZ(),sd0,sd1,sd2);
								double rd0, rd1, rd2;
								gcs->Convert_Transposed_Fractional_Index_To_Global(lrx/gcs->Get_DX(),lry/gcs->Get_DY(),rx_depth/gcs->Get_DZ(),rd0,rd1,rd2);

								GeomTrace &gt = gtarr[iSrc][iRcv].dat();
								gt.setSx(sd0);
								gt.setSy(sd1);
								gt.setSz(sd2);
								gt.setRx(rd0);
								gt.setRy(rd1);
								gt.setRz(rd2);
								gt.setLive(true);
								gt.setSortindex(FFID[iSrc]);
								++iRcv;
							}
						}
					}
				}
			}

			/* OLD CODE - truly centered patch
			for (int iSrc = 0;  iSrc < num_src;  ++iSrc)
			{
				int iRcv = 0;
				for (int iy = -nry_2;  iy <= nry_2;  ++iy)
				{
					for (int ix = -nrx_2;  ix <= nrx_2;  ++ix)
					{
						double lrx = srcx[iSrc] + (double)ix * iline_delta;
						double lry = srcy[iSrc] + (double)iy * xline_delta;
						if (!clip_rx || (
									lrx >= 0.0 && lrx <= gcs->Get_DX() * (double)gcs->Get_NX() &&
									lry >= 0.0 && lry <= gcs->Get_DY() * (double)gcs->Get_NY() 
								)
						   )
						{
							// convert to global coordinates
							double sd0, sd1, sd2;
							gcs->Convert_Transposed_Fractional_Index_To_Global(srcx[iSrc]/gcs->Get_DX(),srcy[iSrc]/gcs->Get_DY(),srcz[iSrc]/gcs->Get_DZ(),sd0,sd1,sd2);
							double rd0, rd1, rd2;
							gcs->Convert_Transposed_Fractional_Index_To_Global(lrx/gcs->Get_DX(),lry/gcs->Get_DY(),rx_depth/gcs->Get_DZ(),rd0,rd1,rd2);

							GeomTrace &gt = gtarr[iSrc][iRcv].dat();
							gt.setSx(sd0);
							gt.setSy(sd1);
							gt.setSz(sd2);
							gt.setRx(rd0);
							gt.setRy(rd1);
							gt.setRz(rd2);
							gt.setLive(true);
							gt.setSortindex(FFID[iSrc]);
							++iRcv;
						}
					}
				}
				for (;  iRcv < max_num_rx;  ++iRcv)
				{
					GeomTrace &gt = gtarr[iSrc][iRcv].dat();
					gt.setSx(0.0);
					gt.setSy(0.0);
					gt.setSz(0.0);
					gt.setRx(0.0);
					gt.setRy(0.0);
					gt.setRz(0.0);
					gt.setLive(false);
					gt.setSortindex(0);
				}
			}
			*/

			ArrND<GeomTrace> gtd(size, filename);

			gtd<<gtarr;

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
		}
	}
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		printf("\n%s Aug-11 2016\nUsage : %s <voxet-file>\n\n",argv[0],argv[0]);
		return -1;
	}

	Voxet* voxet = new Voxet(2,argv[1]);
	Global_Coordinate_System* gcs = voxet->Get_Global_Coordinate_System();

	bool done = false;
	while (!done)
	{
		printf("> ");
		char str[4096];
		gets(str);
		char cmd[4096];
		double c0,c1,c2;
		double d0,d1,d2;
		if (strcmp(str, "quit") == 0 || strcmp(str, "q") == 0)
		{
			done = true;
		}
		else if (sscanf(str, "map %s", cmd) == 1)
		{
			if (gcs->Set_Transpose(cmd))
			{
				printf("Mapping set to UVW -> %s\n",cmd);
			}
			else
			{
				printf("Error! Unknown mapping %s.\n",cmd);
			}
		}
		else if (
			(sscanf(str, "%s %lf %lf %lf", cmd, &c0, &c1, &c2) == 4 && Parse_Command(cmd) > 0) ||
			(sscanf(str, "%s %lf,%lf,%lf", cmd, &c0, &c1, &c2) == 4 && Parse_Command(cmd) > 0)
			)
		{
			switch (Parse_Command(cmd))
			{
			case 1: // g2l
				gcs->Convert_Global_To_Transposed_Fractional_Index(c0,c1,c2,d0,d1,d2);
				d0 *= gcs->Get_DX();
				d1 *= gcs->Get_DY();
				d2 *= gcs->Get_DZ();
				printf("Global coordinate %lf,%lf,%lf => Local coordinate %lf,%lf,%lf\n",c0,c1,c2,d0,d1,d2);
				break;
			case 2: // g2i
				gcs->Convert_Global_To_Transposed_Fractional_Index(c0,c1,c2,d0,d1,d2);
				printf("Global coordinate %lf,%lf,%lf => Index %lf,%lf,%lf\n",c0,c1,c2,d0,d1,d2);
				break;
			case 3: // l2g
				gcs->Convert_Transposed_Fractional_Index_To_Global(c0/gcs->Get_DX(),c1/gcs->Get_DY(),c2/gcs->Get_DZ(),d0,d1,d2);
				printf("Local coordinate %lf,%lf,%lf => Global coordinate %lf,%lf,%lf\n",c0,c1,c2,d0,d1,d2);
				break;
			case 4: // l2i
				d0 = c0 / gcs->Get_DX();
				d1 = c1 / gcs->Get_DY();
				d2 = c2 / gcs->Get_DZ();
				printf("Local coordinate %lf,%lf,%lf => Index %lf,%lf,%lf\n",c0,c1,c2,d0,d1,d2);
				break;
			case 5: // i2g
				gcs->Convert_Transposed_Fractional_Index_To_Global(c0,c1,c2,d0,d1,d2);
				printf("Index %lf,%lf,%lf => Global coordinate %lf,%lf,%lf\n",c0,c1,c2,d0,d1,d2);
				break;
			case 6: // i2l
				d0 = c0 * gcs->Get_DX();
				d1 = c1 * gcs->Get_DY();
				d2 = c2 * gcs->Get_DZ();
				printf("Index %lf,%lf,%lf => Local coordinate %lf,%lf,%lf\n",c0,c1,c2,d0,d1,d2);
				break;
			default:
				break;
			}
		}
		else if (
			(sscanf(str, "create %s", cmd) == 1 && Parse_Command(cmd) > 0)
			)
		{
			switch (Parse_Command(cmd))
			{
			case 7: // create survey
				Create_Survey(gcs);
				break;
			default:
				break;
			}
		}
		else if (
			(sscanf(str, "add %s", cmd) == 1 && Parse_Command(cmd) > 0)
			)
		{
			switch (Parse_Command(cmd))
			{
			case 8: // add depth
				Add_Depth(gcs);
				break;
			default:
				break;
			}
		}
		else if (strcmp(cmd, "dim") == 0)
		{
			printf("NU=%d, NV=%d, NW=%d :: NX=%d, NY=%d, NZ=%d\n",gcs->Get_NU(),gcs->Get_NV(),gcs->Get_NW(),gcs->Get_NX(),gcs->Get_NY(),gcs->Get_NZ());
		}
		else
		{
			printf("  Usage : map [xyz | xzy | yxz | yzx | zxy | zyx]\n");
			printf("          maps UVW axes to XYZ axes (choose one of the mappings).\n\n");
			printf("          quit\n");
			printf("          Quit program.\n\n");
			printf("          <cmd> <c0> <c1> <c2>\n");
			printf("          <cmd> is one of g2l, g2i, l2g, l2i, i2g, i2l.\n");
			printf("                cmd indicate desired conversion.\n");
			printf("                g2l means global-2-local, l2-i means local-2-index etc.\n\n");
			printf("          dim\n");
			printf("                Print dimensions of voxet.\n\n");
			printf("          <floor-cmd> <density-attr> <density-threshold> <c0> <c1> <c2>\n");
			printf("                g2floor, l2floor and i2floor computes seafloor depth for a location given in global, local or index coordinates.\n");
			printf("                only the x and y coordinates are used for the seafloor calculation.\n");
			printf("                seafloor is the shallowest cell where <density-attr> is larger than <density-threshold>.\n\n");
			printf("          create survey\n");
			printf("                Create a survey geometry file to be used for forward modeling.\n");
			printf("                This is an interactive command.\n\n");
		}
	} 

	return 0;
}

