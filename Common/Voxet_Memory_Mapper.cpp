#include <map>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <Voxet_Memory_Mapper.hxx>

Voxet_Memory_Mapper::Voxet_Memory_Mapper()
{
}

Voxet_Memory_Mapper::~Voxet_Memory_Mapper()
{
	for (std::map<std::string,int>::iterator it = _fd.begin();  it != _fd.end();  ++it)
	{
		std::string path = it->first;
		int fd = it->second;
		size_t length = _lengths[path];
		float* addr = _memory_mapped_files[path];
		munlock(addr,length);
		munmap(addr,length);
		close(fd);
	}
}

void Voxet_Memory_Mapper::_MMAP(std::string path)
{
	if (_memory_mapped_files.find(path) == _memory_mapped_files.end())
	{
		struct stat sb;
		int fd = open(path.c_str(),O_RDONLY);
		assert(fd != -1);
		assert(fstat(fd,&sb) != -1);
		size_t length = sb.st_size;
		float* addr = (float*)mmap(0L,length,PROT_READ,MAP_SHARED,fd,0);
		_fd[path] = fd;
		_lengths[path] = length;
		_memory_mapped_files[path] = addr;
		mlock(addr,length);
	}
}

float* Voxet_Memory_Mapper::Get_Memory_Mapped_File(std::string path)
{
	_MMAP(path);
	return _memory_mapped_files[path];
}

size_t Voxet_Memory_Mapper::Get_Length(std::string path)
{
	_MMAP(path);
	return _lengths[path];
}
