#ifndef SEMISIMPLE_LIE_ALGEBRA_CLASS
#define SEMISIMPLE_LIE_ALGEBRA_CLASS
/*
	Class that implements semisimple Lie algebras
*/

#include <vector>
using namespace std;

#include "Fraction.h"
#include "LinearAlgebra.h"

// For most integer data, one has two conflicting concerns about the data-type length.
// Too short: risk of overflow
// Too long: conumption of too much RAM and too many cache misses (thrashing)
// Possible types:
// int8 - char
// int16 - short
// int32 -- int (sometimes long)
// int64 -- long long (sometimes long)


// Integer type used for all persistent data
typedef int LAINT;

// Intermediate-result integer type, for avoiding overflows
typedef int LAXINT;


struct RootConnectionType
{
	LAINT root1, root2, strength; // Uses 1-based indexing
};

struct DynkinDiagramType
{
	vector<LAINT> rtlens; // Root lengths
	vector<RootConnectionType> rtconns; // Root connections
};


// The parameters of the Lie algebra as a single object
// Family and rank
struct LieAlgebraParams {LAINT family, rank;};


// The algebra itself
class LieAlgebra: public LieAlgebraParams
{
	// Made private because these are used internally
	
	// Returns whether the family and rank are valid
	bool MakeDynkin();
	// Returns whether the matrix is singular
	// Also makes associated inverse matrices
	bool MakeMetric(); // Metric matrix
	bool MakeCartan(); // Cartan matirx
	void MakePosRoots(); // Positive roots
	void MakePosRootSum(); // Add them

public:

	// Members are public for convenience
	
	// Overall validity of the algebra
	// Invalid parameters will make it false
	bool IsValid;
	
	// Family: A:1, B:2, C:3, D:4, E:5, F:6, G:7
	DynkinDiagramType dynkin;
	
	// The metric and the Cartan matrix
	// with their inverses and their inverses integerized form.
	// Multiply by the LCM of the denominators to produce
	// an integer-numerator matrix ("num") and an integer denominator ("den").
	Matrix<LAINT> Metric;
	Matrix< Fraction<LAINT> > InverseMetric;
	Matrix<LAINT> InvMetNum;
	LAINT InvMetDen;
	Matrix<LAINT> Cartan;
	Matrix< Fraction<LAINT> > InverseCartan;
	Matrix<LAINT> InvCtnNum;
	LAINT InvCtnDen;
	
	// Positive roots and weights
	// and their sums
	Matrix<LAINT> PosRoots, PosWeights;
	vector<LAINT> PosRootSum, PosWeightSum;
	
	// Call this function to do all the setup
	void Setup();
	
	// Constructors
	// Default arguments are bad values
	LieAlgebra(LAINT family_ = 0, LAINT rank_ = 0)
		{family = family_; rank = rank_; Setup();}
	LieAlgebra(const LieAlgebraParams &Params)
		{family = Params.family; rank = Params.rank; Setup();}
	
	// More with parameter object
	LieAlgebraParams GetParams()
		{LieAlgebraParams Params = *this; return Params;}
	// returns whether the new algebra is valid
	bool UseParams(const LieAlgebraParams &Params)
		{family = Params.family; rank = Params.rank; Setup(); return IsValid;}
};


// For fast sorting
struct LARootHashFunction
{
	size_t operator() (vector<LAINT> &vec);
};

// For caching the algebras

LieAlgebra &GetLieAlgebra(const LieAlgebraParams &Params);
LieAlgebra &GetLieAlgebra(LAINT family, LAINT rank);

void ClearLieAlgebras();

#endif