/*

What's shared by the test suites

*/

#ifndef TEST_SUITES_SHARED
#define TEST_SUITES_SHARED


#include <print>

// Shortcut

void println() {std::println("");}


// Old print:

#include <cstdio>




/*

	Is "mem" in the commend-line arguments?
	
	That's "argv" with length "argc" and the first member skipped over.

*/

inline bool CheckArgMembership(int argc, char **argv, const std::string &mem)
{
	for (int i=1; i<argc; i++)
	{
		std::string avstr(argv[i]);
		if (avstr == mem) return true;
	}
	return false;
}


#endif