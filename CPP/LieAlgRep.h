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


// For safety when calculating the total degeneracies / multiplicities of reps
using TDINT = mpz_class;

// Total degeneracy / multiplicity of an irrep of algebra la specified as MaxWeights
TDINT TotalDegen(const LieAlgebra &la, const LAINT *MaxWeights);

// Conjugate of the irrep
// Returns whether or not self-conjugate
bool RepConjugate(const LieAlgebraParams &AlgParams, const LAINT *MaxWeights, LAINT *ConjgMaxWts);

// Height of the irrep
LAINT RepHeight(const LieAlgebraParams &AlgParams, const LAINT *MaxWeights);

// Conserved-quantity values: sets of (modulus, value)
void RepConserved(const LieAlgebraParams &AlgParams, const LAINT *MaxWeights, LAINT_VECTOR &Conserved);

// Properties collected
enum {
	REP_REALITY_REAL = 0,
	REP_REALITY_PSEUDOREAL = 1,
	REP_REALITY_COMPLEX = -1
};

struct RepProperties
{
	LieAlgebraParams AlgParams;
	LAINT_VECTOR MaxWts, ConjgMaxWts;
	int Reality;
	LAINT_VECTOR Conserved;
	LAINT Height;
	TDINT Size;
	
	void Use(const LieAlgebraParams &AlgParams_, const LAINT *MaxWts_);
};

using SIZE_T_VECTOR = std::vector<size_t>;

// Representation data object
struct LieAlgRep
{
	size_t vlen;
	LAINT_VECTOR Degens;
	LAINT_MATRIX Roots;
	LAINT_MATRIX Weights;
	
	LieAlgRep(): vlen(0) {}
	LieAlgRep(size_t vlen_) {set_vlen(vlen_);}
	
	void set_vlen(size_t vlen_)
		{vlen = vlen_; Degens.clear(); Roots.resize(0,vlen); Weights.resize(0,vlen);}
	
	void AddRoot(LAINT Degen, const LAINT *Root, const LAINT *Weight);
	void AddRoot(LAINT Degen, LAINT_VECTOR &Root, LAINT_VECTOR &Weight)
		{AddRoot(Degen, &Root[0], &Weight[0]);}

	// Sort from max to min in order of root-vector sum, then root-vector elements
	void SortIndices(SIZE_T_VECTOR &Indxs) const;
	
	void Export(LieAlgRep &Rep, SIZE_T_VECTOR &Indxs) const;
	void Export(LieAlgRep &Rep) const;
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
	size_t AddRootOrCount(LAINT Degen, const LAINT_VECTOR &Root, const LAINT_VECTOR &Weight)
		{return AddRootOrCount(Degen, &Root[0], &Weight[0]);}
	
	void Import(LieAlgRep &Rep, bool Append = false);
};


// Find rep from algebra, highest weight -- directly
bool LieAlgebraRepDirect(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts);

inline bool LieAlgebraRepDirect(LieAlgRep &rep, const LieAlgebra &la, const LAINT_VECTOR &MaxWts)
	{return LieAlgebraRepDirect(rep, la, &MaxWts[0]);}

// Find orbit from algebra, highest weight
bool LieAlgebraWeylOrbitGeneral(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts);

inline bool LieAlgebraWeylOrbitGeneral(LieAlgRep &rep, const LieAlgebra &la, const LAINT_VECTOR &MaxWts)
	{return LieAlgebraWeylOrbitGeneral(rep, la, &MaxWts[0]);}

bool LieAlgebraWeylOrbitExplicit(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts);

inline bool LieAlgebraWeylOrbitExplicit(LieAlgRep &rep, const LieAlgebra &la, const LAINT_VECTOR &MaxWts)
	{return LieAlgebraWeylOrbitExplicit(rep, la, &MaxWts[0]);}

bool LieAlgebraWeylOrbit(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts);

inline bool LieAlgebraWeylOrbit(LieAlgRep &rep, const LieAlgebra &la, const LAINT_VECTOR &MaxWts)
	{return LieAlgebraWeylOrbit(rep, la, &MaxWts[0]);}

// Find rep orbits from algebra, highest weight
bool LieAlgebraRepWeylOrbits(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts);

inline bool LieAlgebraRepWeylOrbits(LieAlgRep &rep, const LieAlgebra &la, const LAINT_VECTOR &MaxWts)
	{return LieAlgebraRepWeylOrbits(rep, la, &MaxWts[0]);}

// Find rep from algebra, highest weight -- through orbits
bool LieAlgebraRepByWOs(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts);

inline bool LieAlgebraRepByWOs(LieAlgRep &rep, const LieAlgebra &la, const LAINT_VECTOR &MaxWts)
	{return LieAlgebraRepByWOs(rep, la, &MaxWts[0]);}

// Find rep from algebra, highest weight -- directly
bool LieAlgebraRepresentation(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts);

inline bool LieAlgebraRepresentation(LieAlgRep &rep, const LieAlgebra &la, const LAINT_VECTOR &MaxWts)
	{return LieAlgebraRepresentation(rep, la, &MaxWts[0]);}



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
using LieAlgRepPtr = RefCounter<LieAlgRep>;

LieAlgRepPtr &GetRepObject(enum RepObjType rotype, const LieAlgebraParams &AlgParams, const LAINT *MaxWts);
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, const LieAlgebraParams &AlgParams, const LAINT_VECTOR &MaxWts);
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, LAINT family, LAINT rank, const LAINT *MaxWts);
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, LAINT family, LAINT rank, const LAINT_VECTOR &MaxWts);
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, const LieAlgebra &la, const LAINT *MaxWts);
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, const LieAlgebra &la, const LAINT_VECTOR &MaxWts);

void ClearLieAlgReps();


// Find rep from product algebra -- (semi)simple ones and U(1)'s
struct LieAlgProduct
{
	std::vector<LieAlgebraParams> ParamList;
	LAINT NumU1s;
	
	LAINT get_rank() const;
};

// The max-weights vector is the concatenation of the individual max-weights vectors
// and the U(1) factors

void LieAlgProdRepresentation(LieAlgRep &rep, const LieAlgProduct &lap, 
	enum RepObjType rotype, const LAINT *MaxWts);

inline void LieAlgProdRepresentation(LieAlgRep &rep, const LieAlgProduct &lap,
	enum RepObjType rotype, const LAINT_VECTOR &MaxWts)
	{LieAlgProdRepresentation(rep, lap, rotype, &MaxWts[0]);}

// Caching not implmented here, since the calculating an algebra-product rep
// from individual-algebra reps is not as heavyweight as the calculation of the
// individual-algebra reps themselves.
// There's also the problem of caching the U(1) factors.

// Rep handler
// Reducible reps: counted list of maxweight vectors

struct LACntdMaxWtEntry
{
	LAINT Count;
	LAINT_VECTOR MaxWts;
};

using LACntdMaxWtList = std::vector<LACntdMaxWtEntry>;

// Abstract base class for representation handler
// Subclass for single-algebra and product-of-algebras cases
struct RepHandlerBase
{
	// What's the total rank?
	// Implement in subclass
	virtual LAINT get_rank() const =0;
	
	// Weight vector to rep
	// Implement in subclass
	// Returned as a smart pointer
	// because the alg-prod reps are not cached
	virtual const LieAlgRepPtr GetRepPtr(enum RepObjType rotype, const LAINT *MaxWts) const = 0;
	const LieAlgRepPtr GetRepPtr(enum RepObjType rotype, const LAINT_VECTOR &MaxWts) const
		{return GetRepPtr(rotype, &MaxWts[0]);}
	
	// Counted list of weight vectors to rep
	const LieAlgRepPtr GetRepPtr(enum RepObjType rotype, const LACntdMaxWtList &CWL) const;
	
	// Extract irreps from a rep with sort data
	// as a counted list of max weights for irreps
	// It will be altered as the extraction goes
	const LACntdMaxWtList ExtractWts(enum RepObjType rotype, LieAlgRepBuilder &Bld) const;
};

struct LASnglRepHandler: public RepHandlerBase
{
	LieAlgebraParams Params;
	LASnglRepHandler() {}
	LASnglRepHandler(const LieAlgebraParams &Params_): Params(Params_) {}
	
	LAINT get_rank() const {return Params.rank;}
	const LieAlgRepPtr GetRepPtr(enum RepObjType rotype, const LAINT *MaxWts) const
		{return GetRepObject(rotype,Params,MaxWts);}
};

struct LAProdRepHandler: public RepHandlerBase
{
	LieAlgProduct Params;
	LAProdRepHandler() {}
	LAProdRepHandler(const LieAlgProduct &Params_): Params(Params_) {}
	
	LAINT get_rank() const {return Params.get_rank();}
	const LieAlgRepPtr GetRepPtr(enum RepObjType rotype, const LAINT *MaxWts) const
		{LieAlgRepPtr RepPtr;
			LieAlgProdRepresentation(*RepPtr,Params,rotype,MaxWts);
			return RepPtr;}
};

#endif
