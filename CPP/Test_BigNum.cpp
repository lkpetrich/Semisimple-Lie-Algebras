/*
	Tests the included big-integer class
*/

#include <gmpxx.h>
#include "Test_Shared.h"


int main(void) {

	mpz_class n2(2);
	mpz_class n3(3);
	mpz_class ns = n2 + n3;
	
	std::println("{} + {} = {}",
		n2.get_str(), n3.get_str(), ns.get_str());
	
	mpq_class hf(1,2);
	mpq_class td(1,3);
	mpq_class sx(1,6);
	mpq_class fs = hf + td + sx;
	
	std::println("{} + {} + {} = {}",
		hf.get_str(), td.get_str(), sx.get_str(), fs.get_str());
}