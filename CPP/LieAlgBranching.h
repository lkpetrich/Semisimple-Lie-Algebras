#ifndef SEMISIMPLE_LIE_ALGEBRA_BRANCHING_RULES
#define SEMISIMPLE_LIE_ALGEBRA_BRANCHING_RULES

/*
	For doing subalgebras and branching rules
*/

#include "LieAlgRep.h"

// Branching-rule object

using BRSUBMAT_ENTRY = LAINT_FRACTION;
using BRSUBMAT_VECTOR = std::vector<BRSUBMAT_ENTRY>;
using BRSUBMAT_MATRIX = Matrix<BRSUBMAT_ENTRY>;
using BRSUBMAT_MATRIX_ROW = MatrixRow<BRSUBMAT_ENTRY>;
using BRSUBMAT_MATRIX_ROW_CONST = ConstMatrixRow<BRSUBMAT_ENTRY>;


struct LABrProjector
{
	// Subalgebra parameters
	LieAlgebraParams Params;
	
	// In the format of my Mathematica and Python versions:
	// (original) * (result subalgebra)
	// Uses fractions for full generality
	BRSUBMAT_MATRIX SubMatrix;
	
	// Integerized forms, for convenience in calculating
	LAINT_MATRIX SubMatNum;
	LAINT SubMatDen;
	
	void clear() {SubMatrix.clear(); SubMatNum.clear();}
	void Setup(); // Does precalculation
};

struct LABrancher
{
	LASnglRepHandler OrigAlg; // the original algebra -- a single one
	LAProdRepHandler ResAlg; // The result algebra -- a product one
	
	// Number can vary
	std::vector<LABrProjector> Projectors;
	// The U(1) factors, specified in two different ways
	// This one is for root demotion
	// As in the Mathematica and Python versions, 1-based and not 0-based
	LAINT_VECTOR U1Indices;
	// This one is for interpreting SO(2) as U(1)
	// Like those two, (result) * (original)
	// Uses fractions for full generality
	BRSUBMAT_MATRIX U1SrcVecs;
	// Start using U1SrcVecs; no effect if already using it
	void U1SrcVecStart()
		{if (U1SrcVecs.get_rows() == 0) U1SrcVecs.resize(0,OrigAlg.Params.rank);}
	
	// Expand the U(1) factors and export
	void ExpandU1s(BRSUBMAT_MATRIX &XpndU1s) const;
		
	// Integerized forms, for convenience in calculating
	LAINT_MATRIX U1SrcVecNums;
	LAINT_VECTOR U1SrcVecDens;
	
	void clear() {ResAlg.Params.ParamList.clear(); Projectors.clear();
		U1Indices.clear(); U1SrcVecs.clear();}
	void Setup(); // Does precalculation
	
	// Evaluates the branching
	LACntdMaxWtList DoBranching(const LieAlgRep &Rep) const;
	LACntdMaxWtList DoBranching(const LAINT *MaxWts) const
		{const LieAlgRepPtr RepPtr = OrigAlg.GetRepPtr(RO_REP,MaxWts);
			const LieAlgRep &Rep = *RepPtr;
			return DoBranching(Rep);}
	LACntdMaxWtList DoBranching(const LAINT_VECTOR &MaxWts) const
		{const LieAlgRepPtr RepPtr = OrigAlg.RepHandlerBase::GetRepPtr(RO_REP,MaxWts);
			const LieAlgRep &Rep = *RepPtr;
			return DoBranching(Rep);}
	LACntdMaxWtList DoBranching(const LACntdMaxWtList &CWL) const
		{const LieAlgRepPtr RepPtr = OrigAlg.RepHandlerBase::GetRepPtr(RO_REP,CWL);
			const LieAlgRep &Rep = *RepPtr;
			return DoBranching(Rep);}
};

// The branchers

// Returns branching to the original algebra
LABrancher SubalgSelf(const LieAlgebraParams &OrigAlgParams);

// Turns root RootNo (1-based) into a U(1) factor
LABrancher MakeRootDemoter(const LieAlgebraParams &OrigAlgParams,
	LAINT RootNo);

// Turns roots RootNos (1-based) into U(1) factors
LABrancher MakeMultiRootDemoter(const LieAlgebraParams &OrigAlgParams,
	const LAINT_VECTOR &RootNos);

// Extends the algebra with a root then removes root RootNo (1-based)
// Does not exist for A(n), B1, C1, D2, or D3
LABrancher MakeExtensionSplitter(const LieAlgebraParams &OrigAlgParams,
	LAINT RootNo);

// Decomposes SU(n1*n2*...) into SU(n1)*SU(n2)*...
// Fundamental rep is outer product of fundamental reps
LABrancher SubalgMultSU(const LAINT_VECTOR &SUOrds);

// Decomposes SO/Sp(n1*n2*...) into SO/Sp(n1)*SO/Sp(n2)*...
// Positive n makes SO, negative n Sp
// Vector rep is outer product of vector reps
LABrancher SubalgMultSOSp(const LAINT_VECTOR &SOSpOrds);

// An version of multiplied-n SU(n)
LABrancher SubalgMultAn(const LAINT_VECTOR &ranks);

// BnCnDn version of multiplied-n SO/Sp(n)
LABrancher SubalgMultBnCnDn(const std::vector<LieAlgebraParams> &paramset);

// SO(even) into one or two SO(odd)'s.
// It does SO(2*rank) -> SO(2*sbrk+1) * SO(2*rank-2*sbrk-1)
// Ranks: sbrk, rank-sbrk-1
// It omits the first one if sbrk = 0
// The second one will always be present: sbrk <= rank-2
//
// The extension splitters do
// even -> even + even
// odd -> odd + even
// Handle SO(2) as a U(1) factor by demoting the first root
LABrancher SubalgSOEvenOdd(LAINT rank, LAINT sbrk);

// SU(n) to SO(n)
LABrancher SubalgSUSO(LAINT n);

// SU(2n) to Sp(2n)
LABrancher SubalgSUSp(LAINT n);

// Height to A1, possible for every algebra
LABrancher SubalgHeightA1(const LieAlgebraParams &OrigAlgParams);

// Vector to irrep, possible for all the four infinite families
// It goes from an algebra of family f to one with type dsttype,
// with the source algebra's fundamental or vector rep
// being mapped onto the destination algebra's irrep with highest weights dstwts.
LABrancher SubalgVector(LAINT family,
	const LieAlgebraParams &DestAlgParams,
	const LAINT *DestMaxWeights);

// Extra subalgebras, for the exceptional algebras individually
// Created from the Mathematic version, but with - turned to _ in the names
enum {
	LA_BR_SUBALG_EXTRA_G2A1,
	LA_BR_SUBALG_EXTRA_B3G2,
	LA_BR_SUBALG_EXTRA_D4A2,
	LA_BR_SUBALG_EXTRA_D4G2,
	LA_BR_SUBALG_EXTRA_F4A1,
	LA_BR_SUBALG_EXTRA_F4A1G2,
	LA_BR_SUBALG_EXTRA_E6A2,
	LA_BR_SUBALG_EXTRA_E6G2,
	LA_BR_SUBALG_EXTRA_E6C4,
	LA_BR_SUBALG_EXTRA_E6F4,
	LA_BR_SUBALG_EXTRA_E6A2G2,
	LA_BR_SUBALG_EXTRA_E7A1_1,
	LA_BR_SUBALG_EXTRA_E7A1_2,
	LA_BR_SUBALG_EXTRA_E7A2,
	LA_BR_SUBALG_EXTRA_E7A1A1,
	LA_BR_SUBALG_EXTRA_E7A1G2,
	LA_BR_SUBALG_EXTRA_E7A1F4,
	LA_BR_SUBALG_EXTRA_E7G2C3,
	LA_BR_SUBALG_EXTRA_E8A1_1,
	LA_BR_SUBALG_EXTRA_E8A1_2,
	LA_BR_SUBALG_EXTRA_E8A1_3,
	LA_BR_SUBALG_EXTRA_E8A1A2,
	LA_BR_SUBALG_EXTRA_E8B2,
	LA_BR_SUBALG_EXTRA_E8G2F4,
	LA_BR_NUMBER_OF_SUBALG_EXTRAS
};

// Call with one of the enums above
LABrancher SubalgExtra(const LAINT SAName);


// These all return a brancher object from their input brancher objects
// and other data
// Indexing is 1-based

// Subalgebra #ix of brancher Brn0 gets branched by Brn1,
// making a combined brancher.
LABrancher ConcatBranchers(const LABrancher &Brn0,
	LAINT ix, const LABrancher &Brn1);

// Subalgebra #ix of brancher Brn0 gets a new family number
// if it's A(1), B(1), or C(1): 1, 2, 3
LABrancher BrancherRenameA1B1C1(const LABrancher &Brn0, LAINT ix, LAINT newfam);

// Subalgebra #ix of brancher Brn0 gets flipped between B(2) and C(2)
LABrancher BrancherRenameB2C2(const LABrancher &Brn0, LAINT ix);

// Subalgebra #ix of brancher Brn0 gets flipped between A(3) and D(3)
LABrancher BrancherRenameA3D3(const LABrancher &Brn0, LAINT ix);

// Subalgebra #ix of brancher Brn0 gets split into two A(1)'s if it is D(2)
LABrancher BrancherSplitD2(const LABrancher &Brn0, LAINT ix);

// Subalgebras #ix1 and #ix2 of brancher Brn0 get joined into D(2)
// if they are A(1)/B(1)/C(1)'s.
LABrancher BrancherJoin2A1(const LABrancher &Brn0, LAINT ix1, LAINT ix2);

// Makes conjugates of subalgebras  of brancher Brn0 with indexes in cjixs;
// different for A(n), n>1, D(n), and E(6). Not a true conjugate for D(2n),
// but exchanged anyway.
LABrancher BrancherConjugate(const LABrancher &Brn0, const LAINT_VECTOR &cjixs);

// Subalgebra #ix gets conjugated with newrts specifiying
// the new roots' location if it is D(4).
// newrts is 1-based with length 4 with the second being 2
LABrancher BrancherConjgD4(const LABrancher &Brn0, LAINT ix, const LAINT_VECTOR &newrts);

// Puts the subalgebras of brancher Brn0 into the order specified in neword.
LABrancher BrancherRearrange(const LABrancher &Brn0, const LAINT_VECTOR &neword);

// Mixes the U(1) factors of brancher Brn0, using matrix u1mix
LABrancher BrancherRearrangeU1s(const LABrancher &Brn0, const BRSUBMAT_MATRIX &u1mix);

#endif