#include <stdio.h>

#define dumptype(name,type) printf("%s: %lu\n",name,sizeof(type));

int main()
{
	printf("Dumps the lengths of various C++ integer data types.\n\n");
	dumptype("char",char)
	dumptype("short",short)
	dumptype("int",int)
	dumptype("long",long)
	dumptype("long long",long long)
}