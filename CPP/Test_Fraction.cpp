/*
	Tests the fraction class
*/

// #include <time.h>

#include "Fraction.h"
#include "Test_Shared.h"


// The fraction type
using R = Fraction<int>;

// For writing out a fraction
template<>
struct std::formatter<R>
{
	template <typename ParseContext>
	constexpr auto parse(ParseContext &ctxt) {
		return ctxt.begin();
  	}
  	template <typename FormatContext>
	auto format(const R &frac, FormatContext &ctxt) const {
		if (frac.get_den() != 1)
			return std::format_to(ctxt.out(), "{}/{}", frac.get_num(), frac.get_den());
		else
			return std::format_to(ctxt.out(), "{}", frac.get_num());
	}
};

// No success in overriding formatters for built-in types

// boolean to string: T(rue) or F(alse)
inline std::string B2S(bool x) {return x ? "T" : "F";}


int main()
{
	// Correctness test
	R f1, f2(3), f3(6,-8), f4(f3), f5(2,3);
	R fx1, fx2, fx3, fx4, fx5, fx6;
	
	fx1 = f2; fx2 = -5; fx3.set(16,-12);
	std::println("Expected:   0 3 -3/4 -3/4 3 -5 -4/3");
	std::println("Calculated: {} {} {} {} {} {} {}", f1, f2, f3, f4, fx1, fx2, fx3);
	
	fx1 = -f5; fx2 = f3 + f5; fx3 = f2 + f3; fx4 = f2 - f3; fx5 = f2*f3; fx6 = f2/f3;
	std::println("Expected:   -2/3 -1/12 9/4 15/4 -9/4 -4");
	std::println("Calculated: {} {} {} {} {} {}", fx1, fx2, fx3, fx4, fx5, fx6);
	
	fx1 = f2; fx1 += f3; fx2 = f2; fx2 -= f3; fx3 = f2; fx3 *= f3; fx4 = f2; fx4 /= f3;
	std::println("Expected:   9/4 15/4 -9/4 -4");
	std::println("Calculated: {} {} {} {}", fx1, fx2, fx3, fx4);
	
	bool condsets[6][2];
	condsets[0][0] = f2 == f2; condsets[0][1] = f2 == f3;
	condsets[1][0] = f2 != f2; condsets[1][1] = f2 != f3;
	condsets[2][0] = f2 < f2; condsets[2][1] = f2 < f3;
	condsets[3][0] = f2 > f2; condsets[3][1] = f2 > f3;
	condsets[4][0] = f2 <= f2; condsets[4][1] = f2 <= f3;
	condsets[5][0] = f2 >= f2; condsets[5][1] = f2 >= f3;
	
	std::println("Expected:   T F   F T   F F   F T   T F   T T");
	std::print("Calculated: ");
	for (int i=0; i<6; i++) {
		std::print("{} {}   ", B2S(condsets[i][0]), B2S(condsets[i][1]));
	}
	println();
}