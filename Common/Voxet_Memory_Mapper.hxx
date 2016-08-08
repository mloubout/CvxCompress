#ifndef CVX_ESDRD_MI_TMJ_VOXET_MEMORY_MAPPER_HXX
#define CVX_ESDRD_MI_TMJ_VOXET_MEMORY_MAPPER_HXX

#include <map>
#include <string>

class Voxet_Memory_Mapper
{
	public:
		Voxet_Memory_Mapper();
		virtual ~Voxet_Memory_Mapper();

		float* Get_Memory_Mapped_File(std::string path);
		size_t Get_Length(std::string path);

	private:
		void _MMAP(std::string path);
		std::map<std::string,int> _fd;
		std::map<std::string,size_t> _lengths;
		std::map<std::string,float*> _memory_mapped_files;
};

#endif
