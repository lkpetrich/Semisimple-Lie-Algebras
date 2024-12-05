#ifndef FRACTION_TEMPLATE_CLASS
#define FRACTION_TEMPLATE_CLASS
/*

	Template class for doing fractions
	
	Assumes an integer-like class that support all operations on integers,
	and that interprets 0 as its additive identity
	and 1 as its multiplicative identity.
	Assignment (=) makes a new copy, not varying with the old one.
*/

template<typename N> class Fraction
{
private:
	N num, den;
	void Normalize();
	
public:
	// Accessors
	N get_num() const {return num;}
	N get_den() const {return den;}
	
	// Constructors
	Fraction(); // Zero
	Fraction(const N &num_); // Integer
	Fraction(const N &num_, const N &den_); // Numerator, denominator
	Fraction(const Fraction<N> &frac);
	
	// Assignment
	template<typename NI> Fraction<N> &operator = (const Fraction<NI> &val);
	Fraction<N> &set(const N &num_ = N(0), const N &den_ = N(1));
	
	// Arithmetic
	Fraction<N> &operator + ();
	Fraction<N> operator - ();
	
	template<typename NX> friend Fraction<NX> operator + (const Fraction<NX> &a, const Fraction<NX> &b);
	template<typename NX> friend Fraction<NX> operator - (const Fraction<NX> &a, const Fraction<NX> &b);
	template<typename NX> friend Fraction<NX> operator * (const Fraction<NX> &a, const Fraction<NX> &b);
	template<typename NX> friend Fraction<NX> operator / (const Fraction<NX> &a, const Fraction<NX> &b);
	
	Fraction<N> &operator += (const Fraction<N> &a);
	Fraction<N> &operator -= (const Fraction<N> &a);
	Fraction<N> &operator *= (const Fraction<N> &a);
	Fraction<N> &operator /= (const Fraction<N> &a);
	
	// Comparison
	
	template<class NX> friend bool operator == (const Fraction<NX> &a, const Fraction <NX> &b);
	template<class NX> friend bool operator != (const Fraction<NX> &a, const Fraction <NX> &b);
	template<class NX> friend bool operator > (const Fraction<NX> &a, const Fraction <NX> &b);
	template<class NX> friend bool operator < (const Fraction<NX> &a, const Fraction <NX> &b);
	template<class NX> friend bool operator >= (const Fraction<NX> &a, const Fraction <NX> &b);
	template<class NX> friend bool operator <= (const Fraction<NX> &a, const Fraction <NX> &b);
};

// Greatest common denominator, using Euclid's algorithm
// Returns a nonnegative result even if one or both inputs are negative
template<typename N> N GCD(const N &a, const N &b)
{
	// Avoid dividing by zero
	if (a == N(0)) return N(1);
	if (b == N(0)) return N(1);
	
	// Make positive
	N ax = a;
	if (ax < N(0)) ax = - ax;
	N bx = b;
	if (bx < N(0)) bx = - bx;
	
	// Main algorithm
	N cx;
	while (bx != N(0))
	{
		cx = ax % bx;
		ax = bx;
		bx = cx;
	}
	return ax;
}

// Least common multiple
// Also nonnegative even if one or both inputs are negative
template<typename N> N LCM(const N &a, const N &b)
{
	N res = (a/GCD(a,b))*b;
	if (res < N(0)) res = - res;
	return res;
}

// Absolute value
template<typename N> N abs(const N &a)
{
	return Fraction<N>(abs(a.get_num()),a.get_den());
}

// Reduces a fraction to lowest terms,
// makes the denominator nonnegative
template<typename N> void Fraction<N>::Normalize()
{
	// Handle the special cases with zeros first
	if (num == N(0) && den == N(0)) {}
	else if (num == N(0) && den != N(0)) den = N(1);
	else if (num != N(0) && den == N(0)) num = N(1);
	else
	{
		// No zeros
		N div = GCD(num,den);
		if (den < N(0)) div = - div;
		num /= div;
		den /= div;
	}
}

// Constructors

template<typename N> Fraction<N>::Fraction():
	num(0), den(1) {}

template<typename N> Fraction<N>::Fraction(const N &num_):
	num(num_), den(1) {}

template<typename N> Fraction<N>::Fraction(const N &num_, const N &den_):
	num(num_), den(den_) {Normalize();}

template<typename N> Fraction<N>::Fraction(const Fraction<N> &frac):
	num(frac.num), den(frac.den) {}

// Assignment

template<typename N> template<typename NI> Fraction<N> &Fraction<N>::operator = (const Fraction<NI> &val)
	{num = val.num; den = val.den; return *this;}

template<typename N> Fraction<N> &Fraction<N>::set(const N &num_, const N &den_)
	{num = num_; den = den_; Normalize(); return *this;}

// Arithmetic

template<typename N> Fraction<N> &Fraction<N>::operator + ()
	{return *this;}

template<typename N> Fraction<N> Fraction<N>::operator - ()
	{Fraction<N> res(-num,den); return res;}

template<typename N> Fraction<N> operator + (const Fraction<N> &a, const Fraction<N> &b)
{
	N div = GCD(a.den,b.den);
	N num = a.num*(b.den/div) + b.num*(a.den/div);
	N den = (a.den/div)*b.den;
	Fraction<N> res(num,den);
	return res;
}

template<typename N> Fraction<N> operator - (const Fraction<N> &a, const Fraction<N> &b)
{
	N div = GCD(a.den,b.den);
	N num = a.num*(b.den/div) - b.num*(a.den/div);
	N den = (a.den/div)*b.den;
	Fraction<N> res(num,den);
	return res;
}

template<typename N> Fraction<N> operator * (const Fraction<N> &a, const Fraction<N> &b)
{
	N div1 = GCD(a.num,b.den);
	N div2 = GCD(a.den,b.num);
	N num = (a.num/div1)*(b.num/div2);
	N den = (a.den/div2)*(b.den/div1);
	Fraction<N> res(num,den);
	return res;
}

template<typename N> Fraction<N> operator / (const Fraction<N> &a, const Fraction<N> &b)
{
	N div1 = GCD(a.num,b.num);
	N div2 = GCD(a.den,b.den);
	N num = (a.num/div1)*(b.den/div2);
	N den = (a.den/div2)*(b.num/div1);
	Fraction<N> res(num,den);
	return res;
}

template<typename N> Fraction<N> &Fraction<N>::operator += (const Fraction<N> &a)
{
	N div = GCD(den,a.den);
	num = num*(a.den/div) + a.num*(den/div);
	den = (den/div)*a.den;
	Normalize();
	return *this;
}

template<typename N> Fraction<N> &Fraction<N>::operator -= (const Fraction<N> &a)
{
	N div = GCD(den,a.den);
	num = num*(a.den/div) - a.num*(den/div);
	den = (den/div)*a.den;
	Normalize();
	return *this;
}

template<typename N> Fraction<N> &Fraction<N>::operator *= (const Fraction<N> &a)
{
	N div1 = GCD(num,a.den);
	N div2 = GCD(den,a.num);
	num = (num/div1)*(a.num/div2);
	den = (den/div2)*(a.den/div1);
	Normalize();
	return *this;
}

template<typename N> Fraction<N> &Fraction<N>::operator /= (const Fraction<N> &a)
{
	N div1 = GCD(num,a.num);
	N div2 = GCD(den,a.den);
	num = (num/div1)*(a.den/div2);
	den = (den/div2)*(a.num/div1);
	Normalize();
	return *this;
}

// Comparison

template<typename N> bool operator == (const Fraction<N> &a, const Fraction <N> &b)
{
	N div = GCD(a.den,b.den);
	return a.num*(b.den/div) == b.num*(a.den/div);
}

template<typename N> bool operator != (const Fraction<N> &a, const Fraction <N> &b)
{
	N div = GCD(a.den,b.den);
	return a.num*(b.den/div) != b.num*(a.den/div);
}

template<typename N> bool operator > (const Fraction<N> &a, const Fraction <N> &b)
{
	N div = GCD(a.den,b.den);
	return a.num*(b.den/div) > b.num*(a.den/div);
}

template<typename N> bool operator < (const Fraction<N> &a, const Fraction <N> &b)
{
	N div = GCD(a.den,b.den);
	return a.num*(b.den/div) < b.num*(a.den/div);
}

template<typename N> bool operator >= (const Fraction<N> &a, const Fraction <N> &b)
{
	N div = GCD(a.den,b.den);
	return a.num*(b.den/div) >= b.num*(a.den/div);
}

template<typename N> bool operator <= (const Fraction<N> &a, const Fraction <N> &b)
{
	N div = GCD(a.den,b.den);
	return a.num*(b.den/div) <= b.num*(a.den/div);
}


#endif
