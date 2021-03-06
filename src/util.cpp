#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <unistd.h>
#include <ios>
#include <fstream>
#include <string>
#include <math.h>
#include <assert.h>

double randomVal(int range, int min) {
	int newMin = min;
	bool isNeg = false;
	double randDbl = 0.0;
	if (min < 0) {
		newMin = abs(min);
		isNeg = true;
	} 
	if (isNeg) {
		randDbl = (double) (rand() % range - newMin);
	} else {
		randDbl = (double) (rand() % range + newMin);
	}
	return randDbl;
};


///////////////////////////////////////////////////////////////////////////
// 
// process_mem_usage( double &, double &) - takes two doubles by reference
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB
//
// On failure returns 0.0 0.0

void process_mem_usage(double& vm_usage, double& resident_set)
{
	using std::ios_base;
	using std::ifstream;
	using std::string;

	vm_usage = 0.0;
	resident_set = 0.0;

	// 'file' stat seems to give the most reliable results
	ifstream stat_stream("/proc/self/stat", ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	//
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want
	
	unsigned long vsize;
	long rss;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
		>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
		>> utime >> stime >> cutime >> cstime >> priority >> nice 
		>> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about these

	stat_stream.close();

	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage = vsize / 1024.0;
	resident_set = rss * page_size_kb;
}

double to_degrees(double radians) {
	return radians * (180.0 / M_PI);
}

void compareMin(double* mins, double* tmp) {
	if (tmp[0] < mins[0])
		mins[0] = tmp[0];
	
	if (tmp[1] < mins[1]) 
		mins[1] = tmp[1];
	
	if (tmp[2] < mins[2])
		mins[2] = tmp[2];
}

void compareMax(double* maxs, double* tmp) {
	if (tmp[0] > maxs[0])
		maxs[0] = tmp[0];

	if (tmp[1] > maxs[1])
		maxs[1] = tmp[1];

	if (tmp[2] > maxs[2])
		maxs[2] = tmp[2];
}

int unique_count(int arr[], int len, int uniq_arr[])
{
	if (len <= 0) return 0;
	int unique = 1;
	uniq_arr[0] = arr[0];
	int outer, inner, is_unique = 0;
	for (outer = 1; outer < len; ++outer)
	{
		int is_unique = 1;
		// Don't count negative indexes
		if (arr[outer] < 0) {
			continue;
		}
		for (int inner = 0; is_unique && inner < outer; ++inner)
		{
			if (arr[inner] == arr[outer]) is_unique = 0;
		}
		if (is_unique) {
			
			uniq_arr[unique] = arr[outer];
			++unique;
		}
	}
	return unique;
}

int compareFloat(float a, float b, int maxUlps) {
	assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
	int a_int = *(int*)&a;
	// Make a Int lexicographically ordered as a twos-complement int
	if (a_int < 0) 
		a_int = 0x80000000 - a_int;
	// make b int lexicographically ordered as a twos-complement int
	int b_int = *(int*)&b;
	if (b_int < 0)
		b_int = 0x80000000 - b_int;
	int int_diff = abs(a_int - b_int);
	if (int_diff <= maxUlps)
		return 1;
	return 0;
}

