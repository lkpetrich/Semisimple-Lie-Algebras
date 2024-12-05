#ifndef SEMISIMPLE_LIE_ALGEBRA_CLASS
#define SEMISIMPLE_LIE_ALGEBRA_CLASS
/*
	Class that implements semisimple Lie algebras
*/

#include <vector>

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
using LAINT = int;
using LAINT_VECTOR = std::vector<LAINT>;
using LAINT_MATRIX = Matrix<LAINT>;
using LAINT_MATRIX_ROW = MatrixRow<LAINT>;
using LAINT_MATRIX_ROW_CONST = ConstMatrixRow<LAINT>;
using LAINT_FRACTION = Fraction<LAINT>;
using LAINT_FRACTION_MATRIX = Matrix<LAINT_FRACTION>;

// Intermediate-result integer type, for avoiding overflows
using LAXINT = int;


struct RootConnectionType
{
	LAINT root1, root2, strength; // Uses 1-based indexing
};

struct DynkinDiagramType
{
	LAINT_VECTOR rtlens; // Root lengths
	std::vector<RootConnectionType> rtconns; // Root connections
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
	LAINT_MATRIX Metric;
	LAINT_FRACTION_MATRIX InverseMetric;
	LAINT_MATRIX InvMetNum;
	LAINT InvMetDen;
	LAINT_MATRIX Cartan;
	LAINT_FRACTION_MATRIX InverseCartan;
	LAINT_MATRIX InvCtnNum;
	LAINT InvCtnDen;
	
	// Positive roots and weights
	// and their sums
	LAINT_MATRIX PosRoots, PosWeights;
	LAINT_VECTOR PosRootSum, PosWeightSum;
	
	// Call this function to do all the setup
	void Setup();
	
	// Constructors
	// Default arguments are bad values
	LieAlgebra(LAINT family_ = 0, LAINT rank_ = 0)
		{family = family_; rank = rank_; Setup();}
	LieAlgebra(const LieAlgebraParams &Params)
		{family = Params.family; rank = Params.rank; Setup();}
	
	// More with parameter object
	LieAlgebraParams GetParams() const
		{LieAlgebraParams Params = *this; return Params;}
	// returns whether the new algebra is valid
	bool UseParams(const LieAlgebraParams &Params)
		{family = Params.family; rank = Params.rank; Setup(); return IsValid;}
};


// For fast sorting
struct LARootHashFunction
{
	size_t operator() (LAINT_VECTOR &vec);
};

// For caching the algebras

LieAlgebra &GetLieAlgebra(const LieAlgebraParams &Params);
LieAlgebra &GetLieAlgebra(LAINT family, LAINT rank);

void ClearLieAlgebras();

#endif