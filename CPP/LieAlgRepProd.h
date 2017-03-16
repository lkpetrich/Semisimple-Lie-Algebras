#ifndef SEMISIMPLE_LIE_ALGEBRA_REPRESENTATION_PRODUCTS
#define SEMISIMPLE_LIE_ALGEBRA_REPRESENTATION_PRODUCTS

/*
	For doing rep-product and rep-power functions
*/

#include "LieAlgRep.h"

// Product of two reps

void DoRepProduct(LieAlgRepBuilder &Bld, size_t TotRank,
	LieAlgRep &Rep1, LieAlgRep &Rep2);

template<class WTVEC1, class WTVEC2>
LACntdMaxWtList DecomposeRepProduct(RepHandlerBase &Hdlr, WTVEC1 &Wt1, WTVEC2 &Wt2)
{
	LieAlgRepPtr RepPtr1 = Hdlr.GetRepPtr(RO_REP,Wt1);
	LieAlgRepPtr RepPtr2 = Hdlr.GetRepPtr(RO_REP,Wt2);
	LieAlgRepBuilder Bld;
	DoRepProduct(Bld, Hdlr.get_rank(), *RepPtr1, *RepPtr2);
	return Hdlr.ExtractWts(RO_REP_ORBIT,Bld);
}

// General power of a rep -- plethysm

// Contains the diagram in form (row lengths from max to min, no zero ones)
// and the multiplicity associated with it
typedef vector<LAINT> YoungDiagram;

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

typedef vector<LAYDRepBldEntry> LAYDRepBldList;

struct LAYDCWLEntry
{
	LACntdMaxWtList CWL;
	YoungDiagramEntry YDE;
};

typedef vector<LAYDCWLEntry> LAYDCWLList;

void DoRepPower(LAYDRepBldList &BldList, size_t TotRank,
	LieAlgRep &Rep, LAINT Power);

template <class WTVEC>
LAYDCWLList DecomposeRepPower(RepHandlerBase &Hdlr, WTVEC &Wt, LAINT Power)
{
	LieAlgRepPtr RepPtr = Hdlr.GetRepPtr(RO_REP,Wt);
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
	LieAlgRep &Rep, LAINT Power, LAINT Symm);

template <class WTVEC>
LACntdMaxWtList DecomposeRepPwrSym(RepHandlerBase &Hdlr, WTVEC &Wt, LAINT Power, LAINT Symm)
{
	LieAlgRepPtr RepPtr = Hdlr.GetRepPtr(RO_REP,Wt);
	LieAlgRepBuilder Bld;
	DoRepPwrSym(Bld, Hdlr.get_rank(), *RepPtr, Power, Symm);
	return Hdlr.ExtractWts(RO_REP_ORBIT,Bld);
}

#endif
