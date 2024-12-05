#ifndef SEMISIMPLE_LIE_ALGEBRA_REPRESENTATION_PRODUCTS
#define SEMISIMPLE_LIE_ALGEBRA_REPRESENTATION_PRODUCTS

/*
	For doing rep-product and rep-power functions
*/

#include "LieAlgRep.h"

// Product of two reps

void DoRepProduct(LieAlgRepBuilder &Bld, size_t TotRank,
	const LieAlgRep &Rep1, const LieAlgRep &Rep2);

template<typename WTVEC1, typename WTVEC2>
LACntdMaxWtList DecomposeRepProduct(const RepHandlerBase &Hdlr, const WTVEC1 &Wt1, const WTVEC2 &Wt2)
{
	const LieAlgRepPtr RepPtr1 = Hdlr.GetRepPtr(RO_REP,Wt1);
	const LieAlgRepPtr RepPtr2 = Hdlr.GetRepPtr(RO_REP,Wt2);
	LieAlgRepBuilder Bld;
	DoRepProduct(Bld, Hdlr.get_rank(), *RepPtr1, *RepPtr2);
	return Hdlr.ExtractWts(RO_REP_ORBIT,Bld);
}

// General power of a rep -- plethysm

// Contains the diagram in form (row lengths from max to min, no zero ones)
// and the multiplicity associated with it
using YoungDiagram = LAINT_VECTOR;

struct YoungDiagramEntry
{
	YoungDiagram YD;
	LAINT Mult;
};

struct LAYDRepBldEntry
{
	LieAlgRepBuilder Bld;
	YoungDiagramEntry YDE;
};

using LAYDRepBldList = std::vector<LAYDRepBldEntry>;

struct LAYDCWLEntry
{
	LACntdMaxWtList CWL;
	YoungDiagramEntry YDE;
};

using LAYDCWLList = std::vector<LAYDCWLEntry>;

void DoRepPower(LAYDRepBldList &BldList, size_t TotRank,
	const LieAlgRep &Rep, LAINT Power);

template <class WTVEC>
LAYDCWLList DecomposeRepPower(const RepHandlerBase &Hdlr, const WTVEC &Wt, LAINT Power)
{
	const LieAlgRepPtr RepPtr = Hdlr.GetRepPtr(RO_REP,Wt);
	LAYDRepBldList BldList;
	DoRepPower(BldList, Hdlr.get_rank(), *RepPtr, Power);
	
	LAYDCWLList CWLList;
	for (LAYDRepBldList::iterator BldIter = BldList.begin(); BldIter != BldList.end(); BldIter++)
	{
		LAYDCWLEntry CWLEntry;
		CWLEntry.YDE = BldIter->YDE;
		CWLEntry.CWL = Hdlr.ExtractWts(RO_REP_ORBIT,BldIter->Bld);
		CWLList.push_back(CWLEntry);
	}
	return CWLList;
}

// (Anti)symmetric power of a rep -- plethysm
// Symm is the symmetry type: +1: symmetric, -1: antisymmetric

void DoRepPwrSym(LieAlgRepBuilder &Bld, size_t TotRank,
	const LieAlgRep &Rep, LAINT Power, LAINT Symm);

template <class WTVEC>
LACntdMaxWtList DecomposeRepPwrSym(const RepHandlerBase &Hdlr, const WTVEC &Wt, LAINT Power, LAINT Symm)
{
	const LieAlgRepPtr RepPtr = Hdlr.GetRepPtr(RO_REP,Wt);
	LieAlgRepBuilder Bld;
	DoRepPwrSym(Bld, Hdlr.get_rank(), *RepPtr, Power, Symm);
	return Hdlr.ExtractWts(RO_REP_ORBIT,Bld);
}

#endif
