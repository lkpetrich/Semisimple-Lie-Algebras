/*

What's shared by the test suites

*/

#ifndef TEST_SUITES_SHARED
#define TEST_SUITES_SHARED


#include <stdio.h>
#include <string>


void println() {printf("\n");}

inline bool CheckArgMembership(int argc, char **argv, const string &mem)
{
	for (int i=1; i<argc; i++)
	{
		string avstr(argv[i]);
		if (avstr == mem) return true;
	}
	return false;
}


#endif