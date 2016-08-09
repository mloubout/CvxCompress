//============================================================================
// Name        :

// test.cpp
// Author      : atzv
// Version     :
// Copyright   : 
// Description : Parallel Finite Difference (GPU) using ATM inheriting from 2
//============================================================================
//Modified on 9/24/2014 to assume that source/receiver locations out of Alrt3D/Geom file will be world coordinates.
//Modified on 10/1/2014 to use new ArrNDApp

#include<iostream>
#include<vector>
#include<tr1/array>
#include<ctime>
#include<math.h>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <omp.h>
#include <numa.h>
#include <libgen.h>

#include <ArrND/ArrND.h>
#include <ArrND/ArrNDTask.h>
#include <ArrND/ArrNDAppMastSlaveAtm.h>
#include <Geom/GeomTrace.h>
#include <Geom/GeomTraceAuxiliary.h>

#include <Elastic_Modeling_Job.hxx>
#include <Elastic_Shot.hxx>
#include <Elastic_Propagator.hxx>
#include <Elastic_SEGY_File.hxx>
#include <Voxet.hxx>
#include <Global_Coordinate_System.hxx>
#include <Voxet_Memory_Mapper.hxx>
#include <Variable_Water_Velocity.hxx>

using namespace std;

// Determine number of physical CPU cores.
// Hyper-threaded logical cores are not counted.
// Cache_Size is per core and in whole KB.
int Get_Physical_Core_Count(int& Cache_Size_Per_Core_KB)
{
        FILE* fp = fopen("/proc/cpuinfo","r");
        if (fp != 0L)
        {
                char buf[256];
                int max_phys_id = -1, max_cores = -1, max_cache_size = -1;
                while (fgets(buf, 256, fp) != 0L)
                {
                        int phys_id = -1;
                        if (sscanf(buf, "physical id : %d", &phys_id) == 1)
                        {
                                if (phys_id > max_phys_id) max_phys_id = phys_id;
                        }

                        int cpu_cores = -1;
                        if (sscanf(buf, "cpu cores : %d", &cpu_cores) == 1)
                        {
                                if (cpu_cores > max_cores) max_cores = cpu_cores;
                        }

                        int cache_size = -1;
                        if (sscanf(buf, "cache size : %d", &cache_size) == 1)
                        {
                                if (cache_size > max_cache_size) max_cache_size = cache_size;
                        }
                }
                fclose(fp);
                if (max_phys_id >= 0 && max_cores > 0)
                {
                        if (max_cache_size >= 0)
                        {
                                Cache_Size_Per_Core_KB = max_cache_size / max_cores;
                        }
                        else
                        {
                                Cache_Size_Per_Core_KB = -1;
                        }
                        return (max_phys_id+1) * max_cores;
                }
                else
                {
                        return -1;
                }
        }
        return -1;
}

class FDTask:public ArrNDTask<GeomTrace,GeomTraceAuxiliary>{

protected:
	char* _parmfile;
	int _log_level;
	Voxet_Memory_Mapper* _mapper;
	Variable_Water_Velocity* _Vwxyzt;
	bool debug_wavelet;
	bool debug_slices;
	bool debug_earthmodel_slices;
public:

	FDTask(ArrND<GeomTrace> &in, ArrND<GeomTraceAuxiliary> &in_aux, char* parmfile, int log_level):ArrNDTask<GeomTrace,GeomTraceAuxiliary>(in,in_aux){
		//Read parmfile and assign parameter values; Create shot objects (set source type,amplitues,wavelet);
		//Create SegY file objects; Create Global_Coordinate_System object
		_log_level = log_level;
		_parmfile = parmfile;
		_mapper = new Voxet_Memory_Mapper();
		_Vwxyzt = new Variable_Water_Velocity();

		debug_wavelet=false; //If itask=0, this is set to true; if itask != 0, set to false; Outputs wavelet as ASCII file

		debug_slices=false; //Hardcoded false; Output wavefield slices
		debug_earthmodel_slices=false; //Hardcoded false; Output earth model slices
	};

	~FDTask(){
		delete _Vwxyzt;
		delete _mapper;
	}
	void doit(int itask){


		printf("\n\n\nTask = %d\n",itask);

		//Output wavelet for first shot
		if (itask==0) {
			debug_wavelet=true;
		} else {
			debug_wavelet=false;
		}

		ArrND<GeomTrace> &in=_in;
		ArrND<GeomTraceAuxiliary> &in_aux=_in_aux;

		vector<long> insize = in.size();
		ArrND<GeomTrace> tmp(insize);
		tmp<<in;
		GeomTrace *dat=tmp.datptr();
		ArrND<GeomTraceAuxiliary> tmp_aux(insize);
		tmp_aux<<in_aux;
		GeomTraceAuxiliary *dat_aux=tmp_aux.datptr();

		//All coordinates are local
		double sou_x = dat[0].getSx();
		double sou_y = dat[0].getSy();
		double sou_z = dat[0].getSz();
		cout<<"src="<<sou_x<<" "<<sou_y<<" "<<sou_z<<" "<<endl;

		cout<<"Propagation sub-volume is relative to source"<<endl;

		Elastic_Modeling_Job* job = new Elastic_Modeling_Job(_log_level,_parmfile,_mapper,_Vwxyzt);
		if (job->Is_Valid()){
			if (_Vwxyzt != 0L && _Vwxyzt->Has_Been_Initialized()) _Vwxyzt->Set_Shot_Time(dat[0].getShotTime());
			Voxet* voxet = job->Get_Voxet();
			Global_Coordinate_System* gcs = voxet->Get_Global_Coordinate_System();

			Elastic_Shot* shot = job->Get_Shot_By_Index(0);

			Elastic_SEGY_File* file = shot->Get_SEGY_File_by_Index(0);

			/*
				int file1 = file->Get_File_Index();
				file->Set_File_Index(itask+file1);
			 */

			//Reading in shot id from geom file
			file->Set_File_Index(dat[0].getSortindex());
			file->Set_SeqNo(dat_aux[0].getSeqNo());

			cout<<"File index="<<file->Get_File_Index()<<endl;
			//Convert_Global_To_Transposed_Fractional_Index expects u,v,w as input.
			//ArrND<GeomTrace> file is storing the global u,v,w values
			//gt.sx() is actually returning gt.u and so on.
			//Same is true of receivers

			//Set shot locations
			double x, y, z;
			cout<<"Source in world coord:"<<sou_x<<" "<<sou_y<<" "<<sou_z<<endl;
			gcs->Convert_Global_To_Transposed_Fractional_Index(sou_x,sou_y,sou_z,x,y,z);
			cout<<"Source in fractional index:"<<x<<" "<<y<<" "<<z<<endl;

			shot->Set_Source(x,y,z,dat_aux[0].getSourceInline(),dat_aux[0].getSourceXline());

			//Setting receiver locations
			double* rcv_x = new double[insize[0]]; //This should be the size of the receiver array for one shot
			double* rcv_y = new double[insize[0]];
			double* rcv_z = new double[insize[0]];
			int* rec_il = new int[insize[0]];
			int* rec_xl = new int[insize[0]];
			int* irec = new int[insize[0]];
			int* trcens = new int[insize[0]];
			time_t* acqtime = new time_t[insize[0]];
			int* usec = new int[insize[0]];

			int ilive=0;
			for (int i=0;i<insize[0];i++){
				if (dat[i].isLive()) {
					double rec_x = dat[i].getRx();
					double rec_y = dat[i].getRy();
					double rec_z = dat[i].getRz();
					gcs->Convert_Global_To_Transposed_Fractional_Index(rec_x,rec_y,rec_z,rcv_x[ilive],rcv_y[ilive],rcv_z[ilive]);

					irec[ilive] = dat_aux[i].getReceiverFFID();
					rec_il[ilive] = dat_aux[i].getReceiverInline();
					rec_xl[ilive] = dat_aux[i].getReceiverXline();
					trcens[ilive] = 1;
					acqtime[ilive] = dat[i].getShotTime();
					usec[ilive] = dat[i].getShotTimeUSec();
					//Hardcoded to write out receivers every 100000th location.
					if (i%100000==0){
						printf("%d %12.3f\t%12.3f\t%12.3f\n",i+1,rcv_x[ilive],rcv_y[ilive],rcv_z[ilive]);
					}
					ilive++;
				}
			}

			file->Add_Receiver_Array(ilive,rcv_x,rcv_y,rcv_z,rec_il,rec_xl,trcens,irec,acqtime,usec);
			delete [] rcv_x;
			delete [] rcv_y;
			delete [] rcv_z;
			delete [] rec_xl;
			delete [] rec_il;
			delete [] trcens;
			delete [] irec;
			delete [] acqtime;
			delete [] usec;
			shot->Free_Trace_Resample_Buffers();

			//Recompute the sub volume for the new shot/rec parameters
			job->Compute_Subvolume();
			Elastic_Propagator* prop = new Elastic_Propagator(job);
			prop->Configure();

			//Re-read the earth model for the new subvolume
			prop->Read_Earth_Model();
			if (debug_earthmodel_slices) {
				char str[512];
				sprintf(str,"slices/model_xz_slice_%04d",itask+1);
				job->Write_Earth_Model_XZ_Slice(str,(int)round(shot->Get_Propagation_Source_Y()));
			}
			prop->Propagate_Shot(shot, debug_wavelet, debug_slices);


			//Delete objects - otherwise when objects are re-created all fields are not populated afresh
			delete prop, voxet, gcs, file, shot;
			prop = 0L; voxet = 0L; gcs = 0L;file = 0L; shot = 0L;
			delete job;
			job = 0L;
			cout<<"***Done propagate:"<<itask<<endl;
		}
		return;
	}

};

int main(int argc, char **argv) {
	if (argc != 4 && argc != 5)
        {
		printf("Usage : %s <pathfile> <timeout> <sleeptime> [<numprocesses>]\n",argv[0]);
		printf("        <pathfile> is file created by SeisSpace GUI with all the paths.\n");
		printf("        <timeout> is timeout in seconds before a (presumed) failed shot is restarted.\n");
		printf("        <sleeptime> is sleep time in seconds before idle process checks ATM for more work.\n");
		printf("        <numprocesses> is how many parallel processes are created within a single node.\n");
		printf("                       This argument can be either 1, 2 or 4.\n");
		printf("                       The GPUs are split evenly among the processes.\n");
		printf("                       Using more than one process improves throughput for small shots\n");
		return -1;
	}

	ifstream in_pathfile;
	try{
		in_pathfile.open(argv[1]);
	}
	catch(ifstream::failure & e){
		cout<<"Exception in path file open"<<endl;
	}
	cout<<"Opened input path file -- "<<argv[1]<<endl;

	string comment, geomfile_path, parmfile_path, atmfile_path;
	getline (in_pathfile,comment);
	getline (in_pathfile, atmfile_path);
	cout<<"output file path--"<<atmfile_path<<endl;
	getline (in_pathfile, comment);
	getline (in_pathfile, parmfile_path);
	cout<<"parmfile path --"<<parmfile_path<<endl;
	getline (in_pathfile, comment);
	getline (in_pathfile, geomfile_path);
	cout<<"geometry file path--"<< geomfile_path<<endl;
	//	getline (in_pathfile, comment);
	
	// TMJ set ELAORTHO_CUSTOM_PARMFILE to point to custom parmfile
	/*
	string custom_parmfile_path = parmfile_path + "_custom";
	setenv("ELAORTHO_CUSTOM_PARMFILE",custom_parmfile_path.c_str(),1);
	struct stat st;
	if (stat(custom_parmfile_path.c_str(),&st) == -1)
	{
		printf("Unable to stat(%s). Custom parmfile is REQUIRED for this version of the modeling code.\n",custom_parmfile_path.c_str());
		return -2;
	}
	else
	{
		printf("Changed environment variable ELAORTHO_CUSTOM_PARMFILE to %s.\n",custom_parmfile_path.c_str());
	}
	*/

	in_pathfile.close();
	char parmfile [parmfile_path.size()+1];
	strcpy (parmfile, parmfile_path.c_str());

	int timeOut = 300;
	int sleepTime = 30;
	if (argc>2) {
		timeOut = atoi(argv[2]);
		sleepTime = atoi(argv[3]);
	}
//	cout<<"timeOut = "<<timeOut<<endl;
//	cout<<"timeSleep = "<<sleepTime<<endl;

	int num_processes = argc == 5 ? atoi(argv[4]) : 1;
	if (num_processes != 1 && num_processes != 2 && num_processes != 4)
        {
                printf("Error! <numprocesses> must be either 1, 2 or 4.\n");
                return -2;
        }

	Elastic_Modeling_Job* job = new Elastic_Modeling_Job(4,parmfile_path.c_str(),0L,0L);
	if (job->Is_Valid()){
		if (job->Get_Number_Of_Parallel_Shots() > 0 && job->Get_Number_Of_Parallel_Shots() != num_processes)
		{
			num_processes = job->Get_Number_Of_Parallel_Shots();
			printf("Overriding NUM_PARALLEL_SHOTS from command line argument. New value is %d.\n",num_processes);
		}
	}
	delete job;

	int cu_device_count = 0;
	pid_t child_pid = fork();
	if (child_pid == 0)
	{
		int num_cuda_devices;
                cudaError_t err = cudaGetDeviceCount(&num_cuda_devices);
                return (err == cudaSuccess) ? num_cuda_devices : 0;
	}
	else
	{
		waitpid(child_pid,&cu_device_count,0);
                cu_device_count = cu_device_count >> 8;
	}
        printf("This machine has %d GPUs.\n",cu_device_count);
        if (cu_device_count < num_processes || ((cu_device_count/num_processes)*num_processes) < cu_device_count)
        {
                printf("Error! Number of GPUs must be a multiple of the number of processes.\n");
                return -3;
        }

	int Cache_Size_Per_Core_KB;
	int num_cores = Get_Physical_Core_Count(Cache_Size_Per_Core_KB);
	int threads_per_process = num_cores / num_processes;
	printf("Machine has %d physical cores, using %d openmp threads per process.\n",num_cores,threads_per_process);
	omp_set_num_threads(threads_per_process);

	//
	// num_processes == number of blocks pool of 16 GPUs should be divided into.
	// If num_processes > 1, a new process is created for each block with fork().
	// Each process is bound to the socket that is physically connected to its block of GPUs,
	// and CUDA_VISIBLE_DEVICES is updated to show only the devices the process is allowed to acces.
	//
	int curr_process = 1;
	child_pid = 0;
	if (num_processes > 1)
	{
		bool fork_again = (curr_process < num_processes);
		while (fork_again)
		{
			child_pid = fork();
			if (child_pid == 0)
			{
				// child
				++curr_process;
				fork_again = (curr_process < num_processes);
				sleep(5);
			}
			else
			{
				// parent
				fork_again = false;
			}
		}

		int node_number = ((curr_process-1)*2) / num_processes;
		char node_str[256];
		sprintf(node_str,"%d",node_number);
		printf("PROCESS %d :: numa node %s\n",curr_process,node_str);
		struct bitmask* nodemask = numa_parse_nodestring(node_str);
		numa_bind(nodemask);
		numa_free_nodemask(nodemask);

		int cu_device_count = 16;
		std::string cuda_devices = "";
		for (int iDev = ((curr_process-1)*cu_device_count)/num_processes;  iDev < (curr_process*cu_device_count)/num_processes;  ++iDev)
		{
			char str[256];
			sprintf(str,"%d",iDev);
			if (iDev == ((curr_process-1)*cu_device_count)/num_processes)
			{
				cuda_devices = (std::string)str;
			}
			else
			{
				cuda_devices = cuda_devices + "," + (std::string)str;
			}
		}
		printf("PROCESS %d :: cuda_devices = %s\n",curr_process,cuda_devices.c_str());
		setenv("CUDA_VISIBLE_DEVICES",cuda_devices.c_str(),1);
	}

	//
	// At this point, all processes have been forked off and bound to their respective CPU sockets.
	// CUDA_VISIBLE_DEVICES have been set, so no code change is required in the modeling code.
	//

	time_t start = time(0);

	string ixlfile_path = ((std::string)geomfile_path) + ".ixl";

	ArrND<GeomTrace> in(geomfile_path);
	ArrND<GeomTraceAuxiliary> in_aux(ixlfile_path.c_str());
	vector<long> size1=in.size();
	cout<<"Number of jobs="<<size1[0]<<endl;

	//size of one shot.
	std::vector<long> size2(1);
	size2[0] = size1[1];

	//Sanity check that source id number (stored in gt.sortindex) is unique for each shot
	long nshots = size1[0];
	std::vector<int> shotIDs(size1[0]);

	for (int i =0;i<size1[0];i++) {
		ArrND<GeomTrace> shotGeom(size2);
		shotGeom<<in[i];
		GeomTrace *dat=shotGeom.datptr();
		shotIDs[i] = dat[0].getSortindex();
	}
	std::sort(shotIDs.begin(), shotIDs.end()); //Sort to ascending order
	for (int i=1; i<size1[0]; i++) {
		if (shotIDs[i]==shotIDs[i-1]) {
			printf("Repeated values for shot IDs in the geometry file: shotID[%d]=%d; shotID[%d]=%d\n",
					i-1,shotIDs[i-1],i,shotIDs[i]);
			exit(-1);
		}
	}

	shotIDs.clear();

	ArrND<GeomTrace> inloc(size2);
	ArrND<GeomTraceAuxiliary> inloc_aux(size2);

	FDTask mytask(inloc, inloc_aux, parmfile, 4);
	ArrNDAppMastSlaveAtm<GeomTrace,GeomTraceAuxiliary> ams(in, in_aux, mytask, atmfile_path);
	ams.setTaskTimeOutSecs(timeOut);
	ams.setSleepTime(sleepTime);
	ams.run();

	int child_return_status = 0;
	if (child_pid != 0) waitpid(child_pid,&child_return_status,0);

	time_t end = time(0);
	cout<<"time taken = "<<end-start<<endl;
	return 0;
}
