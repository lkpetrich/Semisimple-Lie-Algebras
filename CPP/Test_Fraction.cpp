/*
	Tests the fraction class
*/

#include <stdio.h>
#include <time.h>
using namespace std;

#include "Fraction.h"

typedef Fraction<int> R;

void print(R &frac)
{
	if (frac.get_den() != 1)
		printf("%d/%d",frac.get_num(),frac.get_den());
	else
		printf("%d",frac.get_num());
}

void println() {printf("\n");}

void println(char *text) {printf("%s\n",text);}

void println(R &frac)
{
	print(frac); println();
}

void println(int nconds, bool *conds)
{
	for (int i=0; i<nconds; i++)
		printf("%s ",conds[i] ? "T" : "F");
	printf("\n");
}

template<class N> void SpeedTest(N dummy)
{
	clock_t t0, t1;
	const int ticks = 100;
	typedef long long INT;
	
	N x(13), a(29), b(79);
	t0 = clock();
	const int niter = 1000000000;
	for (int i=0; i<niter; i++)
		x = a*x + b;
	printf("%lld\n",INT(x));
	t1 = clock();
	double dt = double(t1-t0)/double(CLOCKS_PER_SEC);
	double dtr = dt/niter;
	printf("%d %lg %lg\n",niter,dt,dtr);
	
	Fraction<N> xf(13), af(29,3), bf(79);
	t0 = clock();
	const int niterf = 10000000;
	for (int i=0; i<niterf; i++)
		xf = af*xf + bf;
	printf("%lld/%lld\n",INT(xf.get_num()),INT(xf.get_den()));
	t1 = clock();
	double dtf = double(t1-t0)/double(CLOCKS_PER_SEC);
	double dtrf = dtf/niterf;
	printf("%d %lg %lg\n",niterf,dtf,dtrf);
	
	printf("%lg\n",dtrf/dtr);
	printf("\n");
}


int main()
{
	// Correctness test
	R f1, f2(3), f3(6,-8), f4(f3), f5(2,3);
	
	println("Expected: 0 3 -3/4 -3/4 3 -5 -4/3");
	println(f1);
	println(f2);
	println(f3);
	println(f4);
	f4 = f2; println(f4);
	f4 = -5; println(f4);
	f4.set(16,-12); println(f4);
	println("Expected: -2/3 -1/12 9/4 15/4 -9/4 -4");
	f4 = - f5; println(f4);
	f4 = f3 + f5; println(f4);
	f4 = f2 + f3; println(f4);
	f4 = f2 - f3; println(f4);
	f4 = f2*f3; println(f4);
	f4 = f2/f3; println(f4);
	println("Expected: 9/4 15/4 -9/4 -4");
	f4 = f2; f4 += f3; println(f4);
	f4 = f2; f4 -= f3; println(f4);
	f4 = f2; f4 *= f3; println(f4);
	f4 = f2; f4 /= f3; println(f4);
	println("Expected: T F   F T   F F   F T   T F   T T");
	bool conds[2];
	conds[0] = f2 == f2; conds[1] = f2 == f3;
	println(2,conds);
	conds[0] = f2 != f2; conds[1] = f2 != f3;
	println(2,conds);
	conds[0] = f2 < f2; conds[1] = f2 < f3;
	println(2,conds);
	conds[0] = f2 > f2; conds[1] = f2 > f3;
	println(2,conds);
	conds[0] = f2 <= f2; conds[1] = f2 <= f3;
	println(2,conds);
	conds[0] = f2 >= f2; conds[1] = f2 >= f3;
	println(2,conds);
	
	println();
	/*
	// Speed test
	char dum0;
	short dum1;
	int dum2;
	long long dum3;
	SpeedTest(dum0);
	SpeedTest(dum1);
	SpeedTest(dum2);
	SpeedTest(dum3);
	*/
}