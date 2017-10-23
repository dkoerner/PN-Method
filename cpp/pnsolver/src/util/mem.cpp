#include <util/mem.h>

#include <iostream>





size_t physical_mem_used_by_process()
{
	//MEMORYSTATUSEX memInfo;
	//memInfo.dwLength = sizeof(MEMORYSTATUSEX);
	//GlobalMemoryStatusEx(&memInfo);
	//DWORDLONG totalVirtualMem = memInfo.ullTotalPageFile;
	//DWORDLONG totalPhysMem = memInfo.ullTotalPhys;

	PROCESS_MEMORY_COUNTERS_EX pmc;
	//GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
	bool result = GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
	if( result == false )
		std::cout << "physical_mem_used_by_process: error!\n";
	SIZE_T physMemUsedByMe = pmc.WorkingSetSize;
	return physMemUsedByMe;
}











