#ifndef SEMISIMPLE_LIE_ALGEBRA_REPRESENTATIONS
#define SEMISIMPLE_LIE_ALGEBRA_REPRESENTATIONS
/*
	Functions and classes for semisimple-Lie-algebra representations
	
	Uses the GNU Multple Precision Library to get total degens
	Some of them get huge
*/

#include <gmpxx.h>
#include "LieAlgebra.h"
#include "SmartPointer.h"


// For safety when calculating the total degeneracies / multiplicities of some reps
// typedef long long TDINT;
typedef mpz_class TDINT;

// Total degeneracy / multiplicity of an irrep of algebra la specified as MaxWeights
TDINT TotalDegen(LieAlgebra &la, const LAINT *MaxWeights);

// Conjugate of the irrep
// Returns whether or not self-conjugate
bool RepConjugate(const LieAlgebraParams &AlgParams, const LAINT *MaxWeights, LAINT *ConjgMaxWts);

// Height of the irrep
LAINT RepHeight(const LieAlgebraParams &AlgParams, const LAINT *MaxWeights);

// Conserved-quantity values: sets of (modulus, value)
void RepConserved(const LieAlgebraParams &AlgParams, const LAINT *MaxWeights, vector<LAINT> &Conserved);

// Properties collected
enum {
	REP_REALITY_REAL = 0,
	REP_REALITY_PSEUDOREAL = 1,
	REP_REALITY_COMPLEX = -1
};

struct RepProperties
{
	LieAlgebraParams AlgParams;
	vector<LAINT> MaxWts, ConjgMaxWts;
	int Reality;
	vector<LAINT> Conserved;
	LAINT Height;
	TDINT Size;
	
	void Use(const LieAlgebraParams &AlgParams_, const LAINT *MaxWts_);
};


// Representation data object
struct LieAlgRep
{
	size_t vlen;
	vector<LAINT> Degens;
	Matrix<LAINT> Roots;
	Matrix<LAINT> Weights;
	
	LieAlgRep(): vlen(0) {}
	LieAlgRep(size_t vlen_) {set_vlen(vlen_);}
	
	void set_vlen(size_t vlen_)
		{vlen = vlen_; Degens.clear(); Roots.resize(0,vlen); Weights.resize(0,vlen);}
	
	void AddRoot(LAINT Degen, const LAINT *Root, const LAINT *Weight);
	void AddRoot(LAINT Degen, vector<LAINT> &Root, vector<LAINT> &Weight)
		{AddRoot(Degen, (const LAINT *)&Root[0], (const LAINT *)&Weight[0]);}

	// Sort from max to min in order of root-vector sum, then root-vector elements
	void SortIndices(vector<size_t> &Indxs);
	void Export(LieAlgRep &Rep, vector<size_t> &Indxs);
	void Export(LieAlgRep &Rep);
};

// Data object for building representations
// Uses the rep data object's contents
struct LieAlgRepBuilder: public LieAlgRep
{
	// Indexes the roots, not the weights or the degens
	VectorSet<LAINT,LARootHashFunction> Indexer;
	
	void set_vlen(size_t vlen_)
		{Indexer.set_vlen(vlen_); LieAlgRep::set_vlen(vlen_);}
	size_t get_vlen() {return vlen;}
	
	LieAlgRepBuilder() {}
	LieAlgRepBuilder(size_t vlen_) {set_vlen(vlen_);}
	
	// Adds the root (degen, rt vec, wt vec) if not already present,
	// or its count (degen) if it is.
	size_t AddRootOrCount(LAINT Degen, const LAINT *Root, const LAINT *Weight);
	size_t AddRootOrCount(LAINT Degen, vector<LAINT> &Root, vector<LAINT> &Weight)
		{return AddRootOrCount(Degen, (const LAINT *)&Root[0], (const LAINT *)&Weight[0]);}
	
	void Import(LieAlgRep &Rep, bool Append = false);
};


// Find rep from algebra, highest weight -- directly
bool LieAlgebraRepDirect(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts);

inline bool LieAlgebraRepDirect(LieAlgRep &rep,LieAlgebra &la, vector<LAINT> &MaxWts)
	{return LieAlgebraRepDirect(rep,la,(const LAINT *)&MaxWts[0]);}

// Find orbit from algebra, highest weight
bool LieAlgebraWeylOrbitGeneral(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts);

inline bool LieAlgebraWeylOrbitGeneral(LieAlgRep &rep,LieAlgebra &la, vector<LAINT> &MaxWts)
	{return LieAlgebraWeylOrbitGeneral(rep,la,(const LAINT *)&MaxWts[0]);}

bool LieAlgebraWeylOrbitExplicit(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts);

inline bool LieAlgebraWeylOrbitExplicit(LieAlgRep &rep,LieAlgebra &la, vector<LAINT> &MaxWts)
	{return LieAlgebraWeylOrbitExplicit(rep,la,(const LAINT *)&MaxWts[0]);}

bool LieAlgebraWeylOrbit(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts);

inline bool LieAlgebraWeylOrbit(LieAlgRep &rep,LieAlgebra &la, vector<LAINT> &MaxWts)
	{return LieAlgebraWeylOrbit(rep,la,(const LAINT *)&MaxWts[0]);}

// Find rep orbits from algebra, highest weight
bool LieAlgebraRepWeylOrbits(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts);

inline bool LieAlgebraRepWeylOrbits(LieAlgRep &rep,LieAlgebra &la, vector<LAINT> &MaxWts)
	{return LieAlgebraRepWeylOrbits(rep,la,(const LAINT *)&MaxWts[0]);}

// Find rep from algebra, highest weight -- through orbits
bool LieAlgebraRepByWOs(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts);

inline bool LieAlgebraRepByWOs(LieAlgRep &rep,LieAlgebra &la, vector<LAINT> &MaxWts)
	{return LieAlgebraRepByWOs(rep,la,(const LAINT *)&MaxWts[0]);}

// Find rep from algebra, highest weight -- directly
bool LieAlgebraRepresentation(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts);

inline bool LieAlgebraRepresentation(LieAlgRep &rep,LieAlgebra &la, vector<LAINT> &MaxWts)
	{return LieAlgebraRepresentation(rep,la,(const LAINT *)&MaxWts[0]);}



// For caching the representations and related objects
// What's cached:
enum RepObjType
{
	RO_ORBIT,
	RO_REP_ORBIT,
	RO_REP,
	NUMBER_OF_REP_OBJ_TYPES
};

// Uses reference-counting smart pointers for consistency
typedef RefCounter<LieAlgRep> LieAlgRepPtr;

LieAlgRepPtr &GetRepObject(enum RepObjType rotype, const LieAlgebraParams &AlgParams, const LAINT *MaxWts);
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, const LieAlgebraParams &AlgParams, vector<LAINT> &MaxWts);
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, LAINT family, LAINT rank, const LAINT *MaxWts);
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, LAINT family, LAINT rank, vector<LAINT> &MaxWts);
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, LieAlgebra &la, const LAINT *MaxWts);
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, LieAlgebra &la, vector<LAINT> &MaxWts);

void ClearLieAlgReps();


// Find rep from product algebra -- (semi)simple ones and U(1)'s
struct LieAlgProduct
{
	vector<LieAlgebraParams> ParamList;
	LAINT NumU1s;
	
	LAINT get_rank();
};

// The max-weights vector is the concatenation of the individual max-weights vectors
// and the U(1) factors

void LieAlgProdRepresentation(LieAlgRep &rep, LieAlgProduct &lap, 
	enum RepObjType rotype, const LAINT *MaxWts);

inline void LieAlgProdRepresentation(LieAlgRep &rep,LieAlgProduct &lap,
	enum RepObjType rotype, vector<LAINT> &MaxWts)
	{LieAlgProdRepresentation(rep,lap,rotype,(const LAINT *)&MaxWts[0]);}

// Caching not implmented here, since the calculating an algebra-product rep
// from individual-algebra reps is not as heavyweight as the calculation of the
// individual-algebra reps themselves.
// There's also the problem of caching the U(1) factors.

// Rep handler
// Reducible reps: counted list of maxweight vectors

struct LACntdMaxWtEntry
{
	LAINT Count;
	vector<LAINT> MaxWts;
};

typedef vector<LACntdMaxWtEntry> LACntdMaxWtList;

// Abstract base class for representation handler
// Subclass for single-algebra and product-of-algebras cases
struct RepHandlerBase
{
	// What's the total rank?
	// Implement in subclass
	virtual LAINT get_rank()=0;
	
	// Weight vector to rep
	// Implement in subclass
	// Returned as a smart pointer
	// because the alg-prod reps are not cached
	virtual LieAlgRepPtr GetRepPtr(enum RepObjType rotype, const LAINT *MaxWts) = 0;
	LieAlgRepPtr GetRepPtr(enum RepObjType rotype, vector<LAINT> &MaxWts)
		{return GetRepPtr(rotype, (const LAINT *)&MaxWts[0]);}
	
	// Counted list of weight vectors to rep
	LieAlgRepPtr GetRepPtr(enum RepObjType rotype, LACntdMaxWtList &CWL);
	
	// Extract irreps from a rep with sort data
	// as a counted list of max weights for irreps
	// It will be altered as the extraction goes
	LACntdMaxWtList ExtractWts(enum RepObjType rotype, LieAlgRepBuilder &Bld);
};

struct LASnglRepHandler: public RepHandlerBase
{
	LieAlgebraParams Params;
	LASnglRepHandler() {}
	LASnglRepHandler(LieAlgebraParams &Params_): Params(Params_) {}
	
	LAINT get_rank() {return Params.rank;}
	LieAlgRepPtr GetRepPtr(enum RepObjType rotype, const LAINT *MaxWts)
		{return GetRepObject(rotype,Params,MaxWts);}
};

struct LAProdRepHandler: public RepHandlerBase
{
	LieAlgProduct Params;
	LAProdRepHandler() {}
	LAProdRepHandler(LieAlgProduct &Params_): Params(Params_) {}
	
	LAINT get_rank() {return Params.get_rank();}
	LieAlgRepPtr GetRepPtr(enum RepObjType rotype, const LAINT *MaxWts)
		{LieAlgRepPtr RepPtr; LieAlgProdRepresentation(*RepPtr,Params,rotype,MaxWts);
			return RepPtr;}
};

#endif
