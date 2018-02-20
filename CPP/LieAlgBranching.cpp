/*
	Implements the subalgebra and branching-rule functions
*/

#include "LieAlgBranching.h"


void LABrProjector::Setup()
{
	// Integerize the projection matrix
	size_t ir, ic,
		nr = SubMatrix.get_rows(), nc = SubMatrix.get_cols();
	SubMatNum.resize(nr,nc);
	LAINT lcmden = 1;
	for (ir=0; ir<nr; ir++)
		for (ic=0; ic<nc; ic++)
		{
			BrSubMatEntry &val = SubMatrix(ir,ic);
			lcmden = LCM(lcmden,val.get_den());
		}
	SubMatDen = lcmden;
	for (ir=0; ir<nr; ir++)
		for (ic=0; ic<nc; ic++)
		{
			BrSubMatEntry &val = SubMatrix(ir,ic);
			SubMatNum(ir,ic) = val.get_num()*(lcmden/val.get_den());
		}
}


void LABrancher::Setup()
{
	// Prepare the result algebra for assembly
	ResAlg.Params.ParamList.clear();
	ResAlg.Params.NumU1s = U1Indices.size() + U1SrcVecs.get_rows();
	
	for (vector<LABrProjector>::iterator PrjIter = Projectors.begin();
		PrjIter != Projectors.end(); PrjIter++)
	{
		ResAlg.Params.ParamList.push_back(PrjIter->Params);
		PrjIter->Setup();
	}
	
	// Integerize the U1-vector matrix
	size_t ir, ic,
		nr = U1SrcVecs.get_rows(), nc = U1SrcVecs.get_cols();
	U1SrcVecNums.resize(nr,nc);
	U1SrcVecDens.resize(nr);
	for (ir=0; ir<nr; ir++)
	{
		MatrixRow<BrSubMatEntry> U1SrcVecRow(U1SrcVecs,ir);
		MatrixRow<LAINT> U1SrcVecNumRow(U1SrcVecNums,ir);
		LAINT lcmden = 1;
		for (ic=0; ic<nc; ic++)
		{
			BrSubMatEntry &val = U1SrcVecRow[ic];
			lcmden = LCM(lcmden,val.get_den());
		}
		U1SrcVecDens[ir] = lcmden;
		for (ic=0; ic<nc; ic++)
		{
			BrSubMatEntry &val = U1SrcVecRow[ic];
			U1SrcVecNumRow[ic] = val.get_num()*(lcmden/val.get_den());
		}
	}
}


// For assembling irreps from Weyl orbits
// instead of from their entire contents
static bool IsDomWt(LAINT *Wts, size_t n)
{
	if (n == 0) return true;
	LAINT minwt;
	index_min_v(minwt, Wts, n);
	return (minwt >= 0);
}


LACntdMaxWtList LABrancher::DoBranching(LieAlgRep &Rep)
{
	LAINT ResRank = ResAlg.Params.get_rank();
	vector<LAINT> Roots(ResRank), Weights(ResRank);
	LieAlgRepBuilder Bld(ResRank);
	LieAlgebra &OrigLA = GetLieAlgebra(OrigAlg.Params);
	
	for (size_t i=0; i<Rep.Degens.size(); i++)
	{
		LAINT Degen = Rep.Degens[i];
		MatrixRow<LAINT> OrigRoots(Rep.Roots,i);
		
		// Construct the root and weight vectors
		LAINT ColBase = 0;
		
		// The semisimple part
		for (vector<LABrProjector>::iterator PrjIter = Projectors.begin();
			PrjIter != Projectors.end(); PrjIter++)
		{
			mul_vm(&Roots[ColBase],&OrigRoots[0],PrjIter->SubMatNum);
			
			LieAlgebra &PrjLA = GetLieAlgebra(PrjIter->Params);
			LAINT PrjRank = PrjIter->Params.rank;
			LAINT NextColBase = ColBase + PrjRank;
			for (LAINT j=ColBase; j<NextColBase; j++)
				Roots[j] = (PrjLA.InvCtnDen*Roots[j]) /
					(OrigLA.InvCtnDen*PrjIter->SubMatDen);
			
			mul_vm(&Weights[ColBase],&Roots[ColBase],PrjLA.Cartan);
			for (LAINT j=ColBase; j<NextColBase; j++)
				Weights[j] /= PrjLA.InvCtnDen;
			
			// Next one
			ColBase = NextColBase;
		}
		LAINT NumWts = ColBase;
		
		// The U(1) part
		for (vector<LAINT>::iterator U1InIter = U1Indices.begin();
			U1InIter != U1Indices.end(); U1InIter++)
		{
			Weights[ColBase] = Roots[ColBase] = OrigRoots[(*U1InIter)-1];
			ColBase++;
		}
		mul_mv(&Roots[ColBase],U1SrcVecNums,&OrigRoots[0]);
		LAINT NextColBase = ColBase + U1SrcVecNums.get_rows();
		for (LAINT j=ColBase; j<NextColBase; j++)
			Weights[j] = Roots[j];
		
		// Finally!
		if (IsDomWt(&Weights[0],NumWts))
			Bld.AddRootOrCount(Degen,Roots,Weights);
	}
	
	return ResAlg.ExtractWts(RO_REP_ORBIT,Bld);
}


// Makes submats intended for weight space.
// Args: size of original, special column, 
// list of original root for each result root.
// Original root = -1 makes the special column,
// intended for extension splitting.
static void WtSpcSubMat(Matrix<LAINT> &SubMat, LAINT n,
	vector<LAINT> &SpcCol, vector<LAINT> &OrigRts)
{
	LAINT m = OrigRts.size();
	SubMat.resize(n,m);
	SubMat.fill(0);
	
	for (LAINT i=0; i<m; i++)
	{
		LAINT ix = OrigRts[i];
		if (ix == -1)
		{
			for (LAINT j=0; j<n; j++)
				SubMat(j,i) = SpcCol[j];
		}
		else if (ix >= 1 && ix <= n)
			SubMat(ix-1,i) = 1;
	}
}

// Weight space to root space
static void SubMatWtToRoot(Matrix<BrSubMatEntry>& DestSubMat, LieAlgebraParams &OrigParams,
	LieAlgebraParams &SubParams, Matrix<LAINT> &OrigSubMat)
{
	LieAlgebra &OrigLA = GetLieAlgebra(OrigParams);
	LieAlgebra &SubLA = GetLieAlgebra(SubParams);
	size_t nr = OrigSubMat.get_rows();
	size_t nc = OrigSubMat.get_cols();
	
	Matrix<BrSubMatEntry> IntmdSubMat(nr,nc);
	mul_mm(IntmdSubMat, OrigLA.Cartan, OrigSubMat);
	
	DestSubMat.resize(nr,nc);
	mul_mm(DestSubMat, IntmdSubMat, SubLA.InverseCartan);
}


struct RDSegment
{
	LieAlgebraParams Params;
	vector<LAINT> Roots;
	
	void resize() {Roots.resize(Params.rank);}
};

static void SegmentToProjector(LABrProjector &Proj, LieAlgebraParams &OrigParams,
	vector<LAINT> &SpcCol, RDSegment &Segment)
{
	Proj.Params = Segment.Params;
	Matrix<LAINT> IntmdSubMat;
	WtSpcSubMat(IntmdSubMat, OrigParams.rank, SpcCol, Segment.Roots);
	SubMatWtToRoot(Proj.SubMatrix, OrigParams, Segment.Params, IntmdSubMat);
}

static void SegListToBrancher(LABrancher &Brancher, vector<LAINT> &SpcCol,
	vector<RDSegment> &SegList)
{
	for (vector<RDSegment>::iterator SegIter = SegList.begin();
		SegIter != SegList.end(); SegIter++)
	{
		LABrProjector Proj;
		SegmentToProjector(Proj, Brancher.OrigAlg.Params, SpcCol, *SegIter);
		Brancher.Projectors.push_back(Proj);
	}
}

static void UnpackRDSegmentData(vector<RDSegment> &RootSplit, const LAINT *SegmentData)
{
	if (!SegmentData) return;
	RDSegment Segment;
	const LAINT *SD = SegmentData;
	LAINT numsegs = *(SD++);
	for (LAINT iseg=0; iseg<numsegs; iseg++)
	{
		Segment.Params.family = *(SD++);
		Segment.Params.rank = *(SD++);
		Segment.resize();
		for (LAINT i=0; i<Segment.Params.rank; i++)
			Segment.Roots[i] = *(SD++);
		RootSplit.push_back(Segment);
	}
}


static void SplitByDemotedRoot(vector<RDSegment> &RootSplit, RDSegment &RSMem, LAINT RootNo)
{
	// Appends the algebra-split segments to RootSplit
	LAINT family = RSMem.Params.family;
	LAINT rank = RSMem.Params.rank;
	
	// Find the location of the root to split at
	LAINT ix = -1;
	for (LAINT i=0; i<RSMem.Roots.size(); i++)
	{
		if (RSMem.Roots[i] == RootNo)
		{
			ix = i + 1; // 1-based
			break;
		}
	}
	// If not present, then return the input
	if (ix < 0)
	{
		RootSplit.push_back(RSMem);
		return;
	}
	
	// Create a list of
	// (subalgebra type), (which orig roots because which subalgebra roots)
	// The input roots will be substituted in later.
	// 1-based here
	// Append each one to RootSplit, count off as one goes
	LAINT OldNRSCount = RootSplit.size();
	RDSegment Segment;
	const LAINT *SegData = NULL;
	
	switch(family)
	{
	case 1: // A(n)
		if (ix > 1)
		{
			Segment.Params.family = 1;
			Segment.Params.rank = ix-1;
			Segment.resize();
			for (LAINT i=0; i<Segment.Params.rank; i++)
				Segment.Roots[i] = i + 1;
			RootSplit.push_back(Segment);
		}
		if (ix < rank)
		{
			Segment.Params.family = 1;
			Segment.Params.rank = rank-ix;
			Segment.resize();
			for (LAINT i=0; i<Segment.Params.rank; i++)
				Segment.Roots[i] = i + (ix+1);
			RootSplit.push_back(Segment);		
		}
		break;
	case 2: // B(n)
		if (ix > 1)
		{
			Segment.Params.family = 1;
			Segment.Params.rank = ix-1;
			Segment.resize();
			for (LAINT i=0; i<Segment.Params.rank; i++)
				Segment.Roots[i] = i + 1;
			RootSplit.push_back(Segment);
		}
		if (ix < rank)
		{
			Segment.Params.family = 2;
			Segment.Params.rank = rank-ix;
			Segment.resize();
			for (LAINT i=0; i<Segment.Params.rank; i++)
				Segment.Roots[i] = i + (ix+1);
			RootSplit.push_back(Segment);		
		}
		break;
	case 3: // C(n)
		if (ix > 1)
		{
			Segment.Params.family = 1;
			Segment.Params.rank = ix-1;
			Segment.resize();
			for (LAINT i=0; i<Segment.Params.rank; i++)
				Segment.Roots[i] = i + 1;
			RootSplit.push_back(Segment);
		}
		if (ix < rank)
		{
			Segment.Params.family = 3;
			Segment.Params.rank = rank-ix;
			Segment.resize();
			for (LAINT i=0; i<Segment.Params.rank; i++)
				Segment.Roots[i] = i + (ix+1);
			RootSplit.push_back(Segment);		
		}
		break;
	case 4: // D(n)
		if (ix > 1 && ix < rank-1)
		{
			Segment.Params.family = 1;
			Segment.Params.rank = ix-1;
			Segment.resize();
			for (LAINT i=0; i<Segment.Params.rank; i++)
				Segment.Roots[i] = i + 1;
			RootSplit.push_back(Segment);
		}
		if (ix < rank-1)
		{
			Segment.Params.family = 4;
			Segment.Params.rank = rank-ix;
			Segment.resize();
			for (LAINT i=0; i<Segment.Params.rank; i++)
				Segment.Roots[i] = i + (ix+1);
			RootSplit.push_back(Segment);		
		}
		if (ix == rank-1)
		{
			Segment.Params.family = 1;
			Segment.Params.rank = rank-1;
			Segment.resize();
			for (LAINT i=0; i<Segment.Params.rank-1; i++)
				Segment.Roots[i] = i + 1;
			Segment.Roots[Segment.Params.rank-1] = rank;
			RootSplit.push_back(Segment);		
		}
		if (ix == rank)
		{
			Segment.Params.family = 1;
			Segment.Params.rank = rank-1;
			Segment.resize();
			for (LAINT i=0; i<Segment.Params.rank; i++)
				Segment.Roots[i] = i + 1;
			RootSplit.push_back(Segment);		
		}
		break;
	case 5:
		switch(rank)
		{
		case 6: // E(6)
			switch(ix)
			{
			case 1:
				{
				const LAINT Data[] = {1,   4, 5,  5, 4, 3, 2, 6};
				SegData = Data;
				}
				break;
			case 2:
				{
				const LAINT Data[] = {2,   1, 1,  1,   1, 4,  6, 3, 4, 5};
				SegData = Data;
				}
				break;
			case 3:
				{
				const LAINT Data[] = {3,   1, 2,  1, 2,   1, 2,  4, 5,   1, 1,  6};
				SegData = Data;
				}
				break;
			case 4:
				{
				const LAINT Data[] = {2,   1, 4,  1, 2, 3, 6,   1, 1,  5};
				SegData = Data;
				}
				break;
			case 5:
				{
				const LAINT Data[] = {1,   4, 5,  1, 2, 3, 4, 6};
				SegData = Data;
				}
				break;
			case 6:
				{
				const LAINT Data[] = {1,   1, 5,  1, 2, 3, 4, 5};
				SegData = Data;
				}
				break;
			}
			break;
		case 7: // E(7)
			switch(ix)
			{
			case 1:
				{
				const LAINT Data[] = {1,   4, 6,  6, 5, 4, 3, 2, 7};
				SegData = Data;
				}
				break;
			case 2:
				{
				const LAINT Data[] = {2,   1, 1,  1,   1, 5,  7, 3, 4, 5, 6};
				SegData = Data;
				}
				break;
			case 3:
				{
				const LAINT Data[] = {3,   1, 2,  1, 2,   1, 3,  4, 5, 6,   1, 1,  7};
				SegData = Data;
				}
				break;
			case 4:
				{
				const LAINT Data[] = {2,   1, 4,  1, 2, 3, 7,   1, 2,  5, 6};
				SegData = Data;
				}
				break;
			case 5:
				{
				const LAINT Data[] = {2,   4, 5,  1, 2, 3, 4, 7,   1, 1,  6};
				SegData = Data;
				}
				break;
			case 6:
				{
				const LAINT Data[] = {1,   5, 6,  1, 2, 3, 4, 5, 7};
				SegData = Data;
				}
				break;
			case 7:
				{
				const LAINT Data[] = {1,   1, 6,  1, 2, 3, 4, 5, 6};
				SegData = Data;
				}
				break;
			}
			break;
		case 8: // E(8)
			switch(ix)
			{
			case 1:
				{
				const LAINT Data[] = {1,   4, 7,  7, 6, 5, 4, 3, 2, 8};
				SegData = Data;
				}
				break;
			case 2:
				{
				const LAINT Data[] = {2,   1, 1,  1,   1, 6,  8, 3, 4, 5, 6, 7};
				SegData = Data;
				}
				break;
			case 3:
				{
				const LAINT Data[] = {3,   1, 2,  1, 2,   1, 4,  4, 5, 6, 7,   1, 1,  8};
				SegData = Data;
				}
				break;
			case 4:
				{
				const LAINT Data[] = {2,   1, 4,  1, 2, 3, 8,   1, 3,  5, 6, 7};
				SegData = Data;
				}
				break;
			case 5:
				{
				const LAINT Data[] = {2,   4, 5,  1, 2, 3, 4, 8,   1, 2,  6, 7};
				SegData = Data;
				}
				break;
			case 6:
				{
				const LAINT Data[] = {2,   5, 6,  1, 2, 3, 4, 5, 8,   1, 1,  7};
				SegData = Data;
				}
				break;
			case 7:
				{
				const LAINT Data[] = {1,   5, 7,  1, 2, 3, 4, 5, 6, 8};
				SegData = Data;
				}
				break;
			case 8:
				{
				const LAINT Data[] = {1,   1, 7,  1, 2, 3, 4, 5, 6, 7};
				SegData = Data;
				}
				break;
			}
			break;
		}
		break;
	case 6:
		switch(rank)
		{
		case 4: // F(4)
			switch(ix)
			{
			case 1:
				{
				const LAINT Data[] = {1,   3, 3,  4, 3, 2};
				SegData = Data;
				}
				break;
			case 2:
				{
				const LAINT Data[] = {2,   1, 1,  1,   1, 2,  3, 4};
				SegData = Data;
				}
				break;
			case 3:
				{
				const LAINT Data[] = {2,   1, 2,  1, 2,   1, 1,  4};
				SegData = Data;
				}
				break;
			case 4:
				{
				const LAINT Data[] = {1,   2, 3,  1, 2, 3};
				SegData = Data;
				}
				break;
			}
			break;
		}
		break;
	case 7:
		switch(rank)
		{
		case 2: // G(2)
			switch(ix)
			{
			case 1:
				{
				const LAINT Data[] = {1,   1, 1,  2};
				SegData = Data;
				}
				break;
			case 2:
				{
				const LAINT Data[] = {1,   1, 1,  1};
				SegData = Data;
				}
				break;
			}
			break;
		}
		break;
	}
	
	UnpackRDSegmentData(RootSplit, SegData);
	
	// Plug in which actual original roots
	// for the split records that got added
	LAINT NewNRSCount = RootSplit.size();
	for (LAINT i=OldNRSCount; i<NewNRSCount; i++)
	{
		RDSegment &Seg = RootSplit[i];
		for (vector<LAINT>::iterator SRIter = Seg.Roots.begin();
			SRIter != Seg.Roots.end(); SRIter++)
			*SRIter = RSMem.Roots[(*SRIter)-1];
	}
}

LABrancher MakeMultiRootDemoter(const LieAlgebraParams &OrigAlgParams, vector<LAINT> &RootNos)
{
	LABrancher Brancher;
	Brancher.OrigAlg.Params = OrigAlgParams;
	
	LAINT rank = OrigAlgParams.rank;
	for (size_t i=0; i<RootNos.size(); i++)
	{
		if (RootNos[i] < 1) return Brancher;
		if (RootNos[i] > rank) return Brancher;
		for (size_t j=i+1; j<RootNos.size(); j++)
		{
			if (RootNos[j] == RootNos[i]) return Brancher;
		}
	}
	
	vector<RDSegment> RootSplit, NewRootSplit;
	
	// Initially unsplit
	RDSegment Initial;
	Initial.Params = OrigAlgParams;
	Initial.resize();
	for (LAINT i=0; i<Initial.Params.rank; i++)
		Initial.Roots[i] = i + 1; // 1-based, like in Mma and Python
	RootSplit.push_back(Initial);
	
	// Split by each one separately -- easiest to implement
	for (vector<LAINT>::iterator RNIter = RootNos.begin();
		RNIter != RootNos.end(); RNIter++)
	{
		NewRootSplit.clear();
		for (vector<RDSegment>::iterator RSIter = RootSplit.begin();
			RSIter != RootSplit.end(); RSIter++)
			SplitByDemotedRoot(NewRootSplit, *RSIter, *RNIter);
		RootSplit.swap(NewRootSplit);
	}
	
	vector<LAINT> DummyCol;
	SegListToBrancher(Brancher, DummyCol, RootSplit);
	Brancher.U1Indices = RootNos;
	Brancher.Setup();
	return Brancher;
}


LABrancher MakeRootDemoter(const LieAlgebraParams &OrigAlgParams, LAINT RootNo)
{
	vector<LAINT> RootNos;
	RootNos.push_back(RootNo);
	return MakeMultiRootDemoter(OrigAlgParams, RootNos);
}



LABrancher MakeExtensionSplitter(const LieAlgebraParams &OrigAlgParams, LAINT RootNo)
{
	LABrancher Brancher;
	LAINT family = OrigAlgParams.family;
	LAINT rank = OrigAlgParams.rank;
	Brancher.OrigAlg.Params = OrigAlgParams;
	
	if (RootNo < 1 || RootNo > rank) return Brancher;
	//  D2 = A1*A1, D3 = A3, don't use
	if (family == 4 && rank <= 3) return Brancher;
	
	// Indices are for 1 values; they are 1-based
	vector<LAINT> Indices;
	vector<BrSubMatEntry> SpecialRow;
	
	// The general case: starts off empty
	vector<RDSegment> RootSplit;
	// Numbering of roots is 1-based here
	RDSegment Segment;
	// For special-casing for the exceptional algebras
	const LAINT *SegData = NULL;
	
	switch(family)
	{
	case 1: // A(n)
		{
			Segment.Params = OrigAlgParams;
			Segment.resize();
			LAINT ir = 0;
			for (LAINT i=RootNo+1; i<=rank; i++)
				Segment.Roots[ir++] = i;
			Segment.Roots[ir++] = -1;
			for (LAINT i=1; i<=RootNo-1; i++)
				Segment.Roots[ir++] = i;
			RootSplit.push_back(Segment);
		}
		break;
	case 2: // B(n)
		if (RootNo == 1)
		{
			Segment.Params = OrigAlgParams;
			Segment.resize();
			Segment.Roots[0] = -1;
			for (LAINT i=1; i<Segment.Params.rank; i++)
				Segment.Roots[i] = i + 1;
			RootSplit.push_back(Segment);
		}
		else
		{
			Segment.Params.family = 4;
			Segment.Params.rank = RootNo;
			Segment.resize();
			for (LAINT i=0; i<Segment.Params.rank-1; i++)
				Segment.Roots[i] = (RootNo-1) - i;
			Segment.Roots[Segment.Params.rank-1] = -1;
			RootSplit.push_back(Segment);
			if (RootNo < rank)
			{
				Segment.Params.family = 2;
				Segment.Params.rank = rank-RootNo;
				Segment.resize();
				for (LAINT i=0; i<Segment.Params.rank; i++)
					Segment.Roots[i] = i + (RootNo+1);
				RootSplit.push_back(Segment);
			}
		}
		break;
	case 3: // C(n)
		{
			Segment.Params.family = 3;
			Segment.Params.rank = RootNo;
			Segment.resize();
			for (LAINT i=0; i<Segment.Params.rank-1; i++)
				Segment.Roots[i] = (RootNo-1) - i;
			Segment.Roots[Segment.Params.rank-1] = -1;
			RootSplit.push_back(Segment);
			if (RootNo < rank)
			{
				Segment.Params.family = 3;
				Segment.Params.rank = rank-RootNo;
				Segment.resize();
				for (LAINT i=0; i<Segment.Params.rank; i++)
					Segment.Roots[i] = i + (RootNo+1);
				RootSplit.push_back(Segment);
			}
		}
		break;
	case 4: // D(n)
		if (RootNo == 1)
		{
			Segment.Params = OrigAlgParams;
			Segment.resize();
			Segment.Roots[0] = -1;
			for (LAINT i=1; i<Segment.Params.rank; i++)
				Segment.Roots[i] = i + 1;
			RootSplit.push_back(Segment);
		}
		else if (RootNo >= rank-1)
		{
			LAINT RtNX = (2*rank-1) - RootNo;
			Segment.Params = OrigAlgParams;
			Segment.resize();
			Segment.Roots[0] = RtNX;
			for (LAINT i=1; i<Segment.Params.rank-1; i++)
				Segment.Roots[i] = (rank-1) - i;
			Segment.Roots[Segment.Params.rank-1] = -1;
			RootSplit.push_back(Segment);
		}
		else
		{
			Segment.Params.family = 4;
			Segment.Params.rank = RootNo;
			Segment.resize();
			for (LAINT i=0; i<Segment.Params.rank-1; i++)
				Segment.Roots[i] = (RootNo-1) - i;
			Segment.Roots[Segment.Params.rank-1] = -1;
			RootSplit.push_back(Segment);
			if (RootNo < rank)
			{
				Segment.Params.family = 4;
				Segment.Params.rank = rank-RootNo;
				Segment.resize();
				for (LAINT i=0; i<Segment.Params.rank; i++)
					Segment.Roots[i] = i + (RootNo+1);
				RootSplit.push_back(Segment);
			}
		}
		break;
	case 5: // E(6), E(7), E(8)
		switch(rank)
		{
		case 6: // E(6)
			switch(RootNo)
			{
			case 1:
				{
				const LAINT Data[] = {1,   5, 6,  5, 4, 3, 6, -1, 2};
				SegData = Data;
				}
				break;
			case 2:
				{
				const LAINT Data[] = {2,   1, 1,  1,   1, 5,  5, 4, 3, 6, -1};
				SegData = Data;
				}
				break;
			case 3:
				{
				const LAINT Data[] = {3,   1, 2,  1, 2,   1, 2,  4, 5,   1, 2,  6, -1};
				SegData = Data;
				}
				break;
			case 4:
				{
				const LAINT Data[] = {2,   1, 5,  1, 2, 3, 6, -1,   1, 1,  5};
				SegData = Data;
				}
				break;
			case 5:
				{
				const LAINT Data[] = {1,   5, 6,  1, 2, 3, 6, -1, 4};
				SegData = Data;
				}
				break;
			case 6:
				{
				const LAINT Data[] = {2,   1, 5,  1, 2, 3, 4, 5,  1, 1,  -1};
				SegData = Data;
				}
				break;
			}
			break;
		case 7: // E(7)
			switch(RootNo)
			{
			case 1:
				{
				const LAINT Data[] = {2,   1, 1,  -1,   4, 6,  6, 5, 4, 3, 2, 7};
				SegData = Data;
				}
				break;
			case 2:
				{
				const LAINT Data[] = {2,   1, 2,  -1, 1,   1, 5,  6, 5, 4, 3, 7};
				SegData = Data;
				}
				break;
			case 3:
				{
				const LAINT Data[] = {3,   1, 3,  -1, 1, 2,   1, 3,  4, 5, 6,   1, 1,  7};
				SegData = Data;
				}
				break;
			case 4:
				{
				const LAINT Data[] = {2,   1, 5,  -1, 1, 2, 3, 7,   1, 2,  5, 6};
				SegData = Data;
				}
				break;
			case 5:
				{
				const LAINT Data[] = {2,   4, 6,  -1, 1, 2, 3, 4, 7,   1, 1,  6};
				SegData = Data;
				}
				break;
			case 6:
				{
				const LAINT Data[] = {1,   5, 7,  5, 4, 3, 2, 1, -1, 7};
				SegData = Data;
				}
				break;
			case 7:
				{
				const LAINT Data[] = {1,   1, 7,  -1, 1, 2, 3, 4, 5, 6};
				SegData = Data;
				}
				break;
			}
			break;
		case 8: // E(8)
			switch(RootNo)
			{
			case 1:
				{
				const LAINT Data[] = {1,   4, 8,  -1, 7, 6, 5, 4, 3, 2, 8};
				SegData = Data;
				}
				break;
			case 2:
				{
				const LAINT Data[] = {2,   1, 1,  1,   1, 7,  -1, 7, 6, 5, 4, 3, 8};
				SegData = Data;
				}
				break;
			case 3:
				{
				const LAINT Data[] = {3,   1, 2,  1, 2,   1, 5,  4, 5, 6, 7, -1,   1, 1,  8};
				SegData = Data;
				}
				break;
			case 4:
				{
				const LAINT Data[] = {2,   1, 4,  1, 2, 3, 8,   1, 4,  5, 6, 7, -1};
				SegData = Data;
				}
				break;
			case 5:
				{
				const LAINT Data[] = {2,   4, 5,  1, 2, 3, 4, 8,   1, 3,  6, 7, -1};
				SegData = Data;
				}
				break;
			case 6:
				{
				const LAINT Data[] = {2,   5, 6,  1, 2, 3, 4, 5, 8,   1, 2,  7, -1};
				SegData = Data;
				}
				break;
			case 7:
				{
				const LAINT Data[] = {2,   5, 7,  1, 2, 3, 4, 5, 6, 8,   1, 1,  -1};
				SegData = Data;
				}
				break;
			case 8:
				{
				const LAINT Data[] = {1,   1, 8,  1, 2, 3, 4, 5, 6, 7, -1};
				SegData = Data;
				}
				break;
			}
			break;
		}
		break;
	case 6: // F(4)
		switch(rank)
		{
		case 4: // F(4)
			switch(RootNo)
			{
			case 1:
				{
				const LAINT Data[] = {2,   1, 1,  -1,   3, 3,  4, 3, 2};
				SegData = Data;
				}
				break;
			case 2:
				{
				const LAINT Data[] = {2,   1, 2,  -1, 1,   1, 2,  3, 4};
				SegData = Data;
				}
				break;
			case 3:
				{
				const LAINT Data[] = {2,   1, 3,  -1, 1, 2,   1, 1,  4};
				SegData = Data;
				}
				break;
			case 4:
				{
				const LAINT Data[] = {1,   2, 4,  -1, 1, 2, 3};
				SegData = Data;
				}
				break;
			}
			break;
		}
		break;
	case 7: // G(2)
		switch(rank)
		{
		case 2: // G(2)
			switch(RootNo)
			{
			case 1:
				{
				const LAINT Data[] = {2,   1, 1,  -1,   1, 1,  2};
				SegData = Data;
				}
				break;
			case 2:
				{
				const LAINT Data[] = {1,   1, 2,  -1, 1};
				SegData = Data;
				}
				break;
			}
			break;
		}
		break;
	}
	
	UnpackRDSegmentData(RootSplit, SegData);
	
	LieAlgebra LA = GetLieAlgebra(OrigAlgParams);
	
	// Get the metric diagonal and its maximum
	vector<LAINT> MetricDiag(rank);
	for (LAINT i=0; i<rank; i++)
		MetricDiag[i] = LA.Metric(i,i);
	LAINT MDMax = MetricDiag[0];
	for (LAINT i=1; i<rank; i++)
		MDMax = max(MDMax,MetricDiag[i]);
	
	// The maximum positive root (assume sorted)
	MatrixRow<LAINT> MaxPosRoot(LA.PosRoots, LA.PosRoots.get_rows()-1);
	
	// Assemble the special column vector
	vector<LAINT> SpcCol(rank);
	for (LAINT i=0; i<rank; i++)
		SpcCol[i] = - (MaxPosRoot[i]*MetricDiag[i])/MDMax;
	
	SegListToBrancher(Brancher, SpcCol, RootSplit);
	Brancher.Setup();
	return Brancher;
}


// Additional Branching-Rule Generators


// The first two here have an interesting interpretation: 
// the original Lie-group matrices become outer products of the subgroup matrices. 
// Dynkin has proved that these are the only maximal nonsimple
// subalgebra decompositions of the classical Lie algebras.

// Determination of the simple subalgebras can be very difficult, it must be said. 
// This document will thus have only some of the more interesting or easier-to-find ones

// Reduces SU(n) to a product of SU(m)'s where n = product of m's.
LABrancher SubalgMultSU(vector<LAINT> &SUOrds)
{
	LABrancher Brancher;
	LAINT nso = SUOrds.size();
	LAINT ntot = 1;
	for (vector<LAINT>::iterator suoi = SUOrds.begin();
		suoi != SUOrds.end(); suoi++)
			ntot *= (*suoi);
	LAINT ntrnk = ntot - 1;
	
	Brancher.OrigAlg.Params.family = 1;
	Brancher.OrigAlg.Params.rank = ntrnk;
	
	// Using fractions for generality
	Matrix<BrSubMatEntry> VecProj(ntrnk,ntot);
	VecProj.fill(0);
	for (LAINT k=0; k<ntrnk; k++)
	{
		VecProj(k,k) = 1;
		VecProj(k,k+1) = -1;
	}
	Matrix<BrSubMatEntry> SMat0, SMat1;
	for (LAINT kso=0; kso<nso; kso++)
	{
		LABrProjector Proj;
		LAINT smo = SUOrds[kso];
		LAINT smrnk = smo - 1;
		Proj.Params.family = 1;
		Proj.Params.rank = smrnk;
		
		SMat0.resize(ntot,smrnk);
		SMat0.fill(0);
		LAINT sstr = 1;
		for (LAINT k=kso+1; k<nso; k++)
			sstr *= SUOrds[k];
		for (LAINT k=0; k<ntot; k++)
			for (LAINT ks=0; ks<smrnk; ks++)
			{
				LAINT ko = (k / sstr) % smo;
				if (ko == ks) SMat0(k,ks) = 1;
				if (ko == ks+1) SMat0(k,ks) = -1;
			}
		
		SMat1.resize(ntot,smrnk);
		Proj.SubMatrix.resize(ntrnk,smrnk);
		
		LieAlgebra &LA = GetLieAlgebra(Proj.Params);
		mul_mm(SMat1, SMat0, LA.InverseCartan);
		mul_mm(Proj.SubMatrix, VecProj, SMat1);
		Brancher.Projectors.push_back(Proj);
	}
	
	Brancher.Setup();
	return Brancher;
}

// Reduces SO(n) and Sp(n) to a product of SO(m)'s and Sp(m)' s, 
// where n = product of m's. 
// Positive n means SO(n) and negative n means Sp(-n).
// SO(2) is handled as a U(1) factor.

LABrancher SubalgMultSOSp(vector<LAINT> &SOSpOrds)
{
	LABrancher Brancher;
	LAINT nso = SOSpOrds.size();
	LAINT soprod = 1;
	for (vector<LAINT>::iterator sospoi = SOSpOrds.begin();
		sospoi != SOSpOrds.end(); sospoi++)
			soprod *= (*sospoi);
	
	// Construct a projection matrix for roots onto vector rep. 
	// Also get the type. Much like the previous function, 
	// it's designed to turn the vector rep
	// into crossed diagonal lines of 1 and -1
	LAINT ntot, ntrnk;
	Matrix<BrSubMatEntry> vpbase, vecproj;
	if (soprod > 0)
	{
		if (soprod % 2 == 0)
		{
			// SO(even): Dn
			ntot = soprod;
			ntrnk = ntot/2;
			Brancher.OrigAlg.Params.family = 4;
			Brancher.OrigAlg.Params.rank = ntrnk;
			
			vpbase.resize(ntrnk,ntrnk);
			vpbase.fill(0);
			for (LAINT k=0; k<(ntrnk-1); k++)
			{
				vpbase(k,k) = 1;
				vpbase(k+1,k) = -1;
			}
			vpbase(ntrnk-2,ntrnk-1) = vpbase(ntrnk-1,ntrnk-1) = 1;
			vecproj.resize(ntrnk,ntot);
			vecproj.fill(0);
			for (LAINT k=0; k<ntrnk; k++)
				for (LAINT kx=0; kx<ntrnk; kx++)
				{
					vecproj(k,kx) = vpbase(kx,k);
					vecproj(k,(ntot-1)-kx) = - vpbase(kx,k);
				}
		}
		else
		{
			// SO(odd): Bn
			ntot = soprod;
			ntrnk = (ntot-1)/2;
			Brancher.OrigAlg.Params.family = 2;
			Brancher.OrigAlg.Params.rank = ntrnk;
			
			vpbase.resize(ntrnk,ntrnk);
			vpbase.fill(0);
			for (LAINT k=0; k<(ntrnk-1); k++)
			{
				vpbase(k,k) = 1;
				vpbase(k+1,k) = -1;
			}
			vpbase(ntrnk-1,ntrnk-1) = 1;
			vecproj.resize(ntrnk,ntot);
			vecproj.fill(0);
			for (LAINT k=0; k<ntrnk; k++)
				for (LAINT kx=0; kx<ntrnk; kx++)
				{
					vecproj(k,kx) = vpbase(kx,k);
					vecproj(k,(ntot-1)-kx) = - vpbase(kx,k);
				}
		}
	}
	else
	{
		if (soprod % 2 == 0)
		{
			// Sp(even): Cn
			ntot = - soprod;
			ntrnk = ntot/2;
			Brancher.OrigAlg.Params.family = 3;
			Brancher.OrigAlg.Params.rank = ntrnk;
			
			vpbase.resize(ntrnk,ntrnk);
			vpbase.fill(0);
			for (LAINT k=0; k<(ntrnk-1); k++)
			{
				vpbase(k,k) = 1;
				vpbase(k+1,k) = -1;
			}
			vpbase(ntrnk-1,ntrnk-1) = 2;
			vecproj.resize(ntrnk,ntot);
			vecproj.fill(0);
			for (LAINT k=0; k<ntrnk; k++)
				for (LAINT kx=0; kx<ntrnk; kx++)
				{
					vecproj(k,kx) = vpbase(kx,k);
					vecproj(k,(ntot-1)-kx) = - vpbase(kx,k);
				}
		}
		else
		{
			// Sp(odd) does not exist
			return Brancher;
		}
	}
	
	// Now the subalgebra matrices, with SO(2)/D(1) as a U(1) factor
	Matrix<BrSubMatEntry> SMat0, SMat1;
	vector<BrSubMatEntry> U1FacVec(ntrnk);
	Brancher.U1SrcVecStart();
	for (LAINT kso=0; kso<nso; kso++)
	{
		LAINT smo = SOSpOrds[kso];
		LieAlgebraParams Params;
		if (smo > 0)
		{
			if (smo % 2 == 0)
			{
				// SO(even): Dn -- works for root-vector dimension 1
				LAINT smsz = smo;
				LAINT smrnk = smsz/2;
				Params.family = 4;
				Params.rank = smrnk;
				
				SMat0.resize(ntot,smrnk);
				SMat0.fill(0);
				LAINT sstr = 1;
				for (LAINT k=kso+1; k<nso; k++)
					sstr *= SOSpOrds[k];
				if (sstr < 0) sstr *= -1;
				for (LAINT k=0; k<ntot; k++)
				{
					LAINT ko = (k/sstr) % smsz;
					for (LAINT ks=0; ks<smrnk; ks++)
					{
						if (ko == ks) SMat0(k,ks) = BrSubMatEntry(1,2);
						if (ko == ks+1) SMat0(k,ks) = BrSubMatEntry(-1,2);
						if ((ko == smrnk-2) && (ks == smrnk-1)) SMat0(k,ks) = BrSubMatEntry(1,2);
						if (ko == (smsz-1) - ks) SMat0(k,ks) = BrSubMatEntry(-1,2);
						if (ko == (smsz-1) - (ks+1)) SMat0(k,ks) = BrSubMatEntry(1,2);
						if ((ko == smrnk+1) && (ks == smrnk-1)) SMat0(k,ks) = BrSubMatEntry(-1,2);						
					}
				}
			}
			else
			{
				// SO(odd): Bn
				LAINT smsz = smo;
				LAINT smrnk = (smsz-1)/2;
				Params.family = 2;
				Params.rank = smrnk;
				
				SMat0.resize(ntot,smrnk);
				SMat0.fill(0);
				LAINT sstr = 1;
				for (LAINT k=kso+1; k<nso; k++)
					sstr *= SOSpOrds[k];
				if (sstr < 0) sstr *= -1;
				for (LAINT k=0; k<ntot; k++)
				{
					LAINT ko = (k/sstr) % smsz;
					for (LAINT ks=0; ks<smrnk; ks++)
					{
						if (ko == ks) SMat0(k,ks) = (ks == smrnk-1) ? 1 : BrSubMatEntry(1,2);
						if (ko == ks+1) SMat0(k,ks) = BrSubMatEntry(-1,2);
						if (ko == (smsz-1) - ks) SMat0(k,ks) = (ks == smrnk-1) ? -1 : BrSubMatEntry(-1,2);
						if (ko == (smsz-1) - (ks+1)) SMat0(k,ks) = BrSubMatEntry(1,2);
					}
				}			
			}
		}
		else
		{
			if (smo % 2 == 0)
			{
				// Sp(even): Cn
				LAINT smsz = - smo;
				LAINT smrnk = smsz/2;
				Params.family = 3;
				Params.rank = smrnk;
				
				SMat0.resize(ntot,smrnk);
				SMat0.fill(0);
				LAINT sstr = 1;
				for (LAINT k=kso+1; k<nso; k++)
					sstr *= SOSpOrds[k];
				if (sstr < 0) sstr *= -1;
				for (LAINT k=0; k<ntot; k++)
				{
					LAINT ko = (k/sstr) % smsz;
					for (LAINT ks=0; ks<smrnk; ks++)
					{
						if (ko == ks) SMat0(k,ks) = BrSubMatEntry(1,2);
						if (ko == ks+1) SMat0(k,ks) = BrSubMatEntry(-1,2);
						if (ko == (smsz-1) - ks) SMat0(k,ks) = BrSubMatEntry(-1,2);
						if (ko == (smsz-1) - (ks+1)) SMat0(k,ks) = BrSubMatEntry(1,2);
					}
				}			
			}
			else
			{
				Brancher.Projectors.clear();
				return Brancher;
			}
		}
		
		SMat1.resize(ntrnk,Params.rank);
		mul_mm(SMat1, vecproj, SMat0);
		
		if (Params.family == 2 && Params.rank == 0)
		{
			// Skip over SO(1)
		}
		else if (Params.family == 4 && Params.rank == 1)
		{
			// Treat as a U(1) factor
			mulby_sm(SMat1, BrSubMatEntry(1,2));
			
			for (LAINT k=0; k<ntrnk; k++)
				U1FacVec[k] = SMat1(k,0);
			
			Brancher.U1SrcVecs.AppendVector(U1FacVec);
		}
		else
		{
			// All the others are "real" algebras
			LABrProjector Proj;
			Proj.Params = Params;
			Proj.SubMatrix.resize(ntrnk,Proj.Params.rank);
			LieAlgebra &LA = GetLieAlgebra(Proj.Params);
			mul_mm(Proj.SubMatrix, SMat1, LA.InverseCartan);
			Brancher.Projectors.push_back(Proj);
		}
	}
	
	Brancher.Setup();
	return Brancher;
}


// In: length of root-vector lengths
LABrancher SubalgMultAn(vector<LAINT> &ranks)
{
	vector<LAINT> SUOrds;
	for (vector<LAINT>::iterator rkiter = ranks.begin();
		rkiter != ranks.end(); rkiter++)
	{
		LAINT &rank = *rkiter;
		SUOrds.push_back(rank+1);
	}
	return SubalgMultSU(SUOrds);
}


// In: list of algebra types as lists. For B(n),C(n),D(n): (2,n),(3,n),(4,n).
// D(1) or (4,1) is legitimate.
LABrancher SubalgMultBnCnDn(vector<LieAlgebraParams> &paramset)
{
	vector<LAINT> SOSpOrds;
	for (vector<LieAlgebraParams>::iterator psiter = paramset.begin();
		psiter != paramset.end(); psiter++)
	{
		LieAlgebraParams &Params = *psiter;
		switch(Params.family)
		{
		case 2:
			SOSpOrds.push_back(2*Params.rank+1);
			break;
		case 3:
			SOSpOrds.push_back(-2*Params.rank);
			break;
		case 4:
			SOSpOrds.push_back(2*Params.rank);
			break;
		}
	}
	return SubalgMultSOSp(SOSpOrds);
}

// D(rank) -> B(sbrk) * B(rank-sbrk-1)
LABrancher SubalgSOEvenOdd(LAINT rank, LAINT sbrk)
{
	LABrancher Brancher;
	if (rank < 2) return Brancher;
	if (sbrk < 0 || sbrk > rank-2) return Brancher;
	Brancher.OrigAlg.Params.family = 4;
	Brancher.OrigAlg.Params.rank = rank;
	
	if (sbrk > 0)
	{
		LAINT rsrk = sbrk;
		LABrProjector Proj;
		Proj.Params.family = 2;
		Proj.Params.rank = rsrk;
		
		Proj.SubMatrix.resize(rank,rsrk);
		Proj.SubMatrix.fill(0);
		for (LAINT i=2; i<rsrk; i++)
			Proj.SubMatrix(sbrk-i,i-1) = 1;
		if (sbrk > 1)
		{
			for (LAINT i=0; i<rsrk-1; i++)
				Proj.SubMatrix(0,i) = -1;
			Proj.SubMatrix(0,rsrk-1) = 0;
			
			for (LAINT i=1; i<rsrk; i++)
				Proj.SubMatrix(sbrk-1,i) = -1;
			Proj.SubMatrix(sbrk-1,0) = 0;
		}
		else
		{
			Proj.SubMatrix(0,0) = -1;
		}
		Brancher.Projectors.push_back(Proj);
	}
	
	LAINT rsrk = rank - sbrk - 1;
	LABrProjector Proj;
	Proj.Params.family = 2;
	Proj.Params.rank = rsrk;
	
	Proj.SubMatrix.resize(rank,rsrk);
	Proj.SubMatrix.fill(0);
	for (LAINT i=1; i<rsrk; i++)
		Proj.SubMatrix(sbrk+i-1,i-1) = 1;
	if (sbrk > 0)
		for (LAINT i=0; i<rsrk; i++)
			Proj.SubMatrix(sbrk-1,i) = -1;
	Proj.SubMatrix(rank-1,rsrk-1) = Proj.SubMatrix(rank-2,rsrk-1) = 1;
	
	Brancher.Projectors.push_back(Proj);
	
	Brancher.Setup();
	return Brancher;
}

// Turns SU(n) / A(n-1) into SO(n) / D(n/2) or B((n-1)/2)
LABrancher SubalgSUSO(LAINT n)
{
	LABrancher Brancher;
	Brancher.OrigAlg.Params.family = 1;
	Brancher.OrigAlg.Params.rank = n-1;
	
	LABrProjector Proj;
	if ((n % 2) == 0)
	{
		LAINT rsrk = n/2;
		Proj.Params.family = 4;
		Proj.Params.rank = rsrk;
		
		Proj.SubMatrix.resize(n-1,rsrk);
		Proj.SubMatrix.fill(0);
		for (LAINT k=0; k<rsrk-1; k++)
		{
			Proj.SubMatrix(k,k) = 1;
			Proj.SubMatrix(n-k-2,k) = 1;
		}
		Proj.SubMatrix(rsrk-1,rsrk-2) = -1;
		Proj.SubMatrix(rsrk-1,rsrk-1) = 1;
	}
	else
	{
		LAINT rsrk = (n-1)/2;
		Proj.Params.family = 2;
		Proj.Params.rank = rsrk;
		
		Proj.SubMatrix.resize(n-1,rsrk);
		Proj.SubMatrix.fill(0);
		for (LAINT k=0; k<rsrk; k++)
		{
			Proj.SubMatrix(k,k) = 1;
			Proj.SubMatrix(n-k-2,k) = 1;
		}
	}
	
	Brancher.Projectors.push_back(Proj);
	
	Brancher.Setup();
	return Brancher;
}

// Turns SU(2n) / A(2n-1) into Sp(2n) / C(n)
LABrancher SubalgSUSp(LAINT n)
{
	LABrancher Brancher;
	Brancher.OrigAlg.Params.family = 1;
	Brancher.OrigAlg.Params.rank = 2*n-1;
	
	LABrProjector Proj;
	Proj.Params.family = 3;
	Proj.Params.rank = n;
	
	Proj.SubMatrix.resize(2*n-1,n);
	Proj.SubMatrix.fill(0);
	for (LAINT k=0; k<n-1; k++)
	{
		Proj.SubMatrix(k,k) = 1;
		Proj.SubMatrix(2*n-k-2,k) = 1;
	}
	Proj.SubMatrix(n-1,n-1) = 1;
	
	Brancher.Projectors.push_back(Proj);
	
	Brancher.Setup();
	return Brancher;
}

// Takes each root's height (its total value) and finds out what SU(2)/A1 reps
// account for the heights.
LABrancher SubalgHeightA1(const LieAlgebraParams &OrigAlgParams)
{
	LABrancher Brancher;
	Brancher.OrigAlg.Params = OrigAlgParams;
	
	LABrProjector Proj;
	Proj.Params.family = 1;
	Proj.Params.rank = 1;

	Proj.SubMatrix.resize(OrigAlgParams.rank,1);
	Proj.SubMatrix.fill(1);
	
	Brancher.Projectors.push_back(Proj);
	
	Brancher.Setup();
	return Brancher;
}

// Reduces a member of the four infinite families to
// another algebra by taking the source algebra's vector rep
// (SU and Sp fundamental) to the destination algebra's supplied rep.

// Index-sort object
struct SAVSorter
{
	// Data reference
	LieAlgRepPtr RepPtr;
	
	// The sort function in algorthm.h needs this function
	// It must work like < (true for correct order, false if incorrect order or equal)
	bool operator() (LAINT ix1, LAINT ix2);
};

bool SAVSorter::operator() (LAINT ix1, LAINT ix2)
{
	LAINT rank = RepPtr->vlen;
	
	LAINT *root1 = &RepPtr->Roots(ix1,0);
	LAINT *root2 = &RepPtr->Roots(ix2,0);
	
	// Heights first
	LAINT ttl1 = 0;
	LAINT ttl2 = 0;
	sum(ttl1,root1,rank);
	sum(ttl2,root2,rank);
	if (ttl1 != ttl2) return (ttl2 < ttl1);
	
	// Member by member
	return VecLessThan(root2,root1,rank);
}

LABrancher SubalgVector(LAINT family, const LieAlgebraParams &DestAlgParams,
	const LAINT *DestMaxWeights)
{
	// if the inputs are bad, then return an empty brancher
	LABrancher Brancher;
	Brancher.clear();
	
	// The size of the dest rep
	LieAlgebra &DestLA = GetLieAlgebra(DestAlgParams);
	TDINT RepSizeObj = TotalDegen(DestLA, DestMaxWeights);
	
	// It should not be too big for the C++ code
	if (!RepSizeObj.fits_sshort_p()) {return Brancher;}
	
	LAINT ndst = RepSizeObj.get_ui();
	LAINT n = 0;
	
	switch(family)
	{
	case 1:
		// A(n)
		n = ndst - 1;
		if (n < 1) return Brancher;
		break;
	case 2:
		// B(n)
		n = (ndst - 1) >> 1;
		if (ndst != (2*n+1)) return Brancher;
		if (n < 1) return Brancher;
		break;
	case 3:
		// C(n)
		n = ndst >> 1;
		if (ndst != (2*n)) return Brancher;
		if (n < 1) return Brancher;
		break;
	case 4:
		// D(n) -- exclude D(1), D(2)
		n = ndst >> 1;
		if (ndst != (2*n)) return Brancher;
		if (n < 3) return Brancher;
		break;
	default:
		return Brancher;
	}
	
	// Get the destination rep, expand it, sort it, and de-integerize it
	LieAlgRepPtr DestRepPtr = GetRepObject(RO_REP, DestLA, DestMaxWeights);
	LAINT DestRank = DestRepPtr->vlen;
	size_t NumRepBlks = DestRepPtr->Degens.size();
	
	// Create sort indices
	vector<LAINT> SortIxs(NumRepBlks);
	for (int k=0; k<NumRepBlks; k++)
		SortIxs[k] = k;
	
	// Sort them
	SAVSorter Sorter;
	Sorter.RepPtr = DestRepPtr;
	sort(SortIxs.begin(), SortIxs.end(), Sorter);
	
	LAINT RootDen = DestLA.InvCtnDen;
	Matrix<BrSubMatEntry> RepXpnd(ndst,DestRank);
	LAINT ix = 0;
	for (LAINT k=0; k<NumRepBlks; k++)
	{
		LAINT kx = SortIxs[k];
		LAINT degen = DestRepPtr->Degens[kx];
		LAINT *root = &DestRepPtr->Roots(kx,0);
		for (LAINT l=0; l<degen; l++)
		{
			for (LAINT m=0; m<DestRank; m++)
			{
				BrSubMatEntry deirt(root[m],RootDen);
				RepXpnd(ix,m) = deirt;
			}
			ix++;
		}
	}
	
	// Is the destination rep (pseudo)real?
	if (family >= 2)
	{
		for (LAINT k=0; k<ndst; k++)
		{
			Fraction<LAINT> *root1 = &RepXpnd(k,0);
			Fraction<LAINT> *root2 = &RepXpnd(ndst-k-1,0);
			for (LAINT m=0; m<DestRank; m++)
				if (root1[m] != (- root2[m])) return Brancher;
		}
	}
	
	// The inverse of the source rep
	BrSubMatEntry zero(0,1);
	BrSubMatEntry one(1,1);
	BrSubMatEntry half(1,2);
	Matrix<BrSubMatEntry> SrcRepInv(n,ndst);
	SrcRepInv.fill(zero);
	switch(family)
	{
	case 1:
		// A(n)
		for (LAINT k=0; k<n; k++)
		{
			SrcRepInv(k,k) = one;
			SrcRepInv(k,k+1) = - one;
		}
		break;
	case 2:
		// B(n)
		for (LAINT k=0; k<(n-1); k++)
		{
			SrcRepInv(k,k) = half;
			SrcRepInv(k,k+1) = - half;
			SrcRepInv(k,ndst-k-1) = - half;
			SrcRepInv(k,ndst-k-2) = half;			
		}
		SrcRepInv(n-1,n-1) = half;
		SrcRepInv(n-1,ndst-n) = - half;
		break;
	case 3:
		// C(n)
		for (LAINT k=0; k<(n-1); k++)
		{
			SrcRepInv(k,k) = half;
			SrcRepInv(k,k+1) = - half;
			SrcRepInv(k,ndst-k-1) = - half;
			SrcRepInv(k,ndst-k-2) = half;			
		}
		SrcRepInv(n-1,n-1) = one;
		SrcRepInv(n-1,ndst-n) = - one;
		break;
	case 4:
		// D(n)
		for (LAINT k=0; k<(n-2); k++)
		{
			SrcRepInv(k,k) = half;
			SrcRepInv(k,k+1) = - half;
			SrcRepInv(k,ndst-k-1) = - half;
			SrcRepInv(k,ndst-k-2) = half;			
		}
		SrcRepInv(n-2,n-2) = half;
		SrcRepInv(n-2,n-1) = half;
		SrcRepInv(n-2,ndst-n+1) = - half;
		SrcRepInv(n-2,ndst-n) = - half;
		SrcRepInv(n-1,n-2) = half;
		SrcRepInv(n-1,n-1) = - half;
		SrcRepInv(n-1,ndst-n+1) = - half;
		SrcRepInv(n-1,ndst-n) = half;
	}	
	
	// Final setup
	LieAlgebraParams OrigAlgParams;
	OrigAlgParams.family = family;
	OrigAlgParams.rank = n;
	Brancher.OrigAlg.Params = OrigAlgParams;
	
	LABrProjector Proj;
	Proj.Params = DestAlgParams;
	Proj.SubMatrix.resize(n,DestRank);
	
	// Create the projection matrix
	mul_mm(Proj.SubMatrix, SrcRepInv, RepXpnd);
	
	Brancher.Projectors.push_back(Proj);
	
	Brancher.Setup();
	return Brancher;
}

// Exceptional-algebra extras
// Generated from the Mathematica version
static const LAINT SubalgExtraData_G2A1[] = {7, 2, 1, 1, 1, 1, 1};
static const LAINT SubalgExtraData_B3G2[] = {2, 3, 1, 7, 2, 0, 1, 1, 0, 0, 1};
static const LAINT SubalgExtraData_D4A2[] = {4, 4, 1, 1, 2, 1, 0, -1, 1, 1, 0, 1, 0};
static const LAINT SubalgExtraData_D4G2[] = {4, 4, 1, 7, 2, 0, 1, 1, 0, 0, 1, 0, 1};
static const LAINT SubalgExtraData_F4A1[] = {6, 4, 1, 1, 1, 1, 1, 1, 1};
static const LAINT SubalgExtraData_F4A1G2[] = {6, 4, 2, 1, 1, 0, 0, 0, 1, 7, 2, 2, 3, -1, 0, 0, -1, 0, 0};
static const LAINT SubalgExtraData_E6A2[] = {5, 6, 1, 1, 2, 1, 0, 0, 1, -1, 0, 0, 1, 1, 0, 2, -1};
static const LAINT SubalgExtraData_E6G2[] = {5, 6, 1, 7, 2, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, -1, 0};
static const LAINT SubalgExtraData_E6C4[] = {5, 6, 1, 3, 4, 0, 1, 0, 0, 1, 0, 0, 0, -1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1};
static const LAINT SubalgExtraData_E6F4[] = {5, 6, 1, 6, 4, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0};
static const LAINT SubalgExtraData_E6A2G2[] = {5, 6, 2, 1, 2, 1, 0, 0, 1, -1, -1, 1, 0, 0, 1, 0, 0, 7, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};
static const LAINT SubalgExtraData_E7A1_1[] = {5, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
static const LAINT SubalgExtraData_E7A1_2[] = {5, 7, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1};
static const LAINT SubalgExtraData_E7A2[] = {5, 7, 1, 1, 2, 0, 2, 0, 2, 0, -1, 1, 0, 0, 1, 1, 0, 0, -2};
static const LAINT SubalgExtraData_E7A1A1[] = {5, 7, 2, 1, 1, 2, 1, -2, 1, 1, 1, -2, 1, 1, 0, 0, 0, 1, 0, 0, 0};
static const LAINT SubalgExtraData_E7A1G2[] = {5, 7, 2, 1, 1, 0, 0, 0, 0, 0, 0, 1, 7, 2, 0, -1, -2, -2, 0, -1, 1, 1, 0, 1, 0, 1, 1, 2};
static const LAINT SubalgExtraData_E7A1F4[] = {5, 7, 2, 1, 1, 0, 0, 0, 0, 0, 1, 0, 6, 4, 0, 0, 1, 0, 1, 2, 2, 1, -2, -3, -4, -2, 1, 2, 2, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0};
static const LAINT SubalgExtraData_E7G2C3[] = {5, 7, 2, 7, 2, 0, 0, 0, 0, 0, -1, 1, 2, 0, 0, 0, 0, -1, 0, 3, 3, 0, 1, 0, 1, 0, 0, 0, 1, 1, -1, -2, -1, 0, 1, 0, 1, 0, 0, 0, 0, 0};
static const LAINT SubalgExtraData_E8A1_1[] = {5, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
static const LAINT SubalgExtraData_E8A1_2[] = {5, 8, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1};
static const LAINT SubalgExtraData_E8A1_3[] = {5, 8, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1};
static const LAINT SubalgExtraData_E8A1A2[] = {5, 8, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 0, 0, 1, -1, 0, 1, 1, 0, -1, 1, 0, 0, 1, -1, -2};
static const LAINT SubalgExtraData_E8B2[] = {5, 8, 1, 2, 2, 0, 1, 0, -1, 1, 0, 0, 1, 0, 1, -1, 0, 0, 1, 0, -1};
static const LAINT SubalgExtraData_E8G2F4[] = {5, 8, 2, 7, 2, 0, 0, 0, 1, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, -1, 6, 4, 0, 0, 0, -1, -1, -1, -1, 0, 0, -1, -2, -1, 1, 1, 2, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 2, 1};

static const LAINT *SubalgExtraData[] = {
	SubalgExtraData_G2A1,
	SubalgExtraData_B3G2,
	SubalgExtraData_D4A2,
	SubalgExtraData_D4G2,
	SubalgExtraData_F4A1,
	SubalgExtraData_F4A1G2,
	SubalgExtraData_E6A2,
	SubalgExtraData_E6G2,
	SubalgExtraData_E6C4,
	SubalgExtraData_E6F4,
	SubalgExtraData_E6A2G2,
	SubalgExtraData_E7A1_1,
	SubalgExtraData_E7A1_2,
	SubalgExtraData_E7A2,
	SubalgExtraData_E7A1A1,
	SubalgExtraData_E7A1G2,
	SubalgExtraData_E7A1F4,
	SubalgExtraData_E7G2C3,
	SubalgExtraData_E8A1_1,
	SubalgExtraData_E8A1_2,
	SubalgExtraData_E8A1_3,
	SubalgExtraData_E8A1A2,
	SubalgExtraData_E8B2,
	SubalgExtraData_E8G2F4,
};
// End of Mma-generated code

// Use that data
LABrancher SubalgExtra(LAINT SAName)
{
	LABrancher Brancher;
	if (SAName < 0 || SAName >= LA_BR_NUMBER_OF_SUBALG_EXTRAS) return Brancher;
	
	LAINT *SED = (LAINT *)SubalgExtraData[SAName];
	
	Brancher.OrigAlg.Params.family = *(SED++);
	Brancher.OrigAlg.Params.rank = *(SED++);
	LAINT rank = Brancher.OrigAlg.Params.rank;
	
	LAINT NumSubAlgs = *(SED++);
	for (LAINT i=0; i<NumSubAlgs; i++)
	{
		LABrProjector Proj;
		Proj.Params.family = *(SED++);
		Proj.Params.rank = *(SED++);
		LAINT rsrk = Proj.Params.rank;

		Proj.SubMatrix.resize(rank,rsrk);
		for (LAINT ir=0; ir<rank; ir++)
			for (LAINT ic=0; ic<rsrk; ic++)
				Proj.SubMatrix(ir,ic) = *(SED++);
		
		Brancher.Projectors.push_back(Proj);
	}
	
	Brancher.Setup();
	return Brancher;
}


// These all return a brancher object from their input brancher objects
// and other data
// Indexing is 1-based

// Subalgebra #ix of brancher Brn0 gets branched by Brn1,
// making a combined brancher.
LABrancher ConcatBranchers(LABrancher &Brn0, LAINT ix, LABrancher &Brn1)
{
	// Check
	LABrancher Brancher;
	if (ix <= 0 || ix > Brn0.Projectors.size()) return Brancher;
	LABrProjector &Proj01 = Brn0.Projectors[ix-1];
	if (Proj01.Params.family != Brn1.OrigAlg.Params.family) return Brancher;
	if (Proj01.Params.rank != Brn1.OrigAlg.Params.rank) return Brancher;
	
	// Same starting algebra
	Brancher.OrigAlg.Params = Brn0.OrigAlg.Params;
	
	// Do the projection matrices
	// Brn1 ones get run into the appropriate Brn0 one.
	for (LAINT i=0; i<Brn0.Projectors.size(); i++)
	{
		LABrProjector &Proj0 = Brn0.Projectors[i];
		if (i == ix-1)
		{
			for (vector<LABrProjector>::iterator PrjIter = Brn1.Projectors.begin();
				PrjIter != Brn1.Projectors.end(); PrjIter++)
			{
				LABrProjector &Proj1 = *PrjIter;
				LABrProjector Proj;
				Proj.Params = Proj1.Params;
				LAINT rows = Brn0.OrigAlg.Params.rank;
				LAINT cols = Proj1.Params.rank;
				Proj.SubMatrix.resize(rows,cols);
				mul_mm(Proj.SubMatrix, Proj01.SubMatrix, Proj1.SubMatrix);
				Brancher.Projectors.push_back(Proj);
			}
		}
		else
		{
			Brancher.Projectors.push_back(Proj0);
		}
	}
	
	// Do the U(1) factors
	// Need to do a lot of setting up
	// to change from Brn1's root space to Brn0's root space.
	Brancher.U1Indices = Brn0.U1Indices;
	Brancher.U1SrcVecs = Brn0.U1SrcVecs;
	
	LAINT rank0 = Brn0.OrigAlg.Params.rank;
	LAINT rank1 = Brn1.OrigAlg.Params.rank;
	vector<BrSubMatEntry> NewU1Vec(rank0);
	vector<BrSubMatEntry> NewU1VecSrc(rank1);
	
	Brancher.U1SrcVecStart();
	
	for (vector<LAINT>::iterator U1IxIter = Brn1.U1Indices.begin();
		U1IxIter != Brn1.U1Indices.end(); U1IxIter++)
	{
		fill(NewU1VecSrc.begin(), NewU1VecSrc.end(), 0);
		LAINT uix = *U1IxIter;
		NewU1VecSrc[uix-1] = 1;
		
		mul_mv(NewU1Vec, Proj01.SubMatrix, NewU1VecSrc);
		Brancher.U1SrcVecs.AppendVector(NewU1Vec);
	}
	
	for (LAINT i=0; i<Brn1.U1SrcVecs.get_rows(); i++)
	{
		MatrixRow<BrSubMatEntry> U1SV(Brn1.U1SrcVecs,i);
		copy(U1SV.begin(), U1SV.end(), NewU1VecSrc.begin());
		
		mul_mv(NewU1Vec, Proj01.SubMatrix, NewU1VecSrc);
		Brancher.U1SrcVecs.AppendVector(NewU1Vec);
	}
	
	Brancher.Setup();
	return Brancher;
}


// Subalgebra #ix of brancher Brn0 gets a new family number
// if it's A(1), B(1), or C(1): 1, 2, 3
LABrancher BrancherRenameA1B1C1(LABrancher &Brn0, LAINT ix, LAINT newfam)
{
	LABrancher Brancher;
	if (ix <= 0 || ix > Brn0.Projectors.size()) return Brancher;
	if (newfam < 1 || newfam > 3) return Brancher;
	LABrProjector &Proj01 = Brn0.Projectors[ix-1];
	if (Proj01.Params.family < 1 || Proj01.Params.family > 3) return Brancher;
	if (Proj01.Params.rank != 1) return Brancher;
	
	// Same starting algebra, U(1)'s
	Brancher.OrigAlg.Params = Brn0.OrigAlg.Params;
	Brancher.U1Indices = Brn0.U1Indices;
	Brancher.U1SrcVecs = Brn0.U1SrcVecs;
	
	// Do the projection matrices
	for (LAINT i=0; i<Brn0.Projectors.size(); i++)
	{
		LABrProjector &Proj0 = Brn0.Projectors[i];
		if (i == ix-1)
		{
			LABrProjector Proj = Proj0;
			Proj.Params.family = newfam;
			Brancher.Projectors.push_back(Proj);
		}
		else
		{
			Brancher.Projectors.push_back(Proj0);
		}
	}
	
	Brancher.Setup();
	return Brancher;
}


// Subalgebra #ix of brancher Brn0 gets flipped between B(2) and C(2)
LABrancher BrancherRenameB2C2(LABrancher &Brn0, LAINT ix)
{
	LABrancher Brancher;
	if (ix <= 0 || ix > Brn0.Projectors.size()) return Brancher;
	LABrProjector &Proj01 = Brn0.Projectors[ix-1];
	if (Proj01.Params.family != 2 && Proj01.Params.family != 3) return Brancher;
	if (Proj01.Params.rank != 2) return Brancher;
	
	// Same starting algebra, U(1)'s
	Brancher.OrigAlg.Params = Brn0.OrigAlg.Params;
	Brancher.U1Indices = Brn0.U1Indices;
	Brancher.U1SrcVecs = Brn0.U1SrcVecs;
	
	// Do the projection matrices
	for (LAINT i=0; i<Brn0.Projectors.size(); i++)
	{
		LABrProjector &Proj0 = Brn0.Projectors[i];
		if (i == ix-1)
		{
			LABrProjector Proj = Proj0;
			Proj.Params.family = 5 - Proj.Params.family;
			for (LAINT k=0; k<Proj.SubMatrix.get_rows(); k++)
				swap(Proj.SubMatrix(k,0), Proj.SubMatrix(k,1));
			Brancher.Projectors.push_back(Proj);
		}
		else
		{
			Brancher.Projectors.push_back(Proj0);
		}
	}
	
	Brancher.Setup();
	return Brancher;
}


// Subalgebra #ix of brancher Brn0 gets flipped between A(3) and D(3)
LABrancher BrancherRenameA3D3(LABrancher &Brn0, LAINT ix)
{
	LABrancher Brancher;
	if (ix <= 0 || ix > Brn0.Projectors.size()) return Brancher;
	LABrProjector &Proj01 = Brn0.Projectors[ix-1];
	if (Proj01.Params.family != 1 && Proj01.Params.family != 4) return Brancher;
	if (Proj01.Params.rank != 3) return Brancher;
	
	// Same starting algebra, U(1)'s
	Brancher.OrigAlg.Params = Brn0.OrigAlg.Params;
	Brancher.U1Indices = Brn0.U1Indices;
	Brancher.U1SrcVecs = Brn0.U1SrcVecs;
	
	// Do the projection matrices
	for (LAINT i=0; i<Brn0.Projectors.size(); i++)
	{
		LABrProjector &Proj0 = Brn0.Projectors[i];
		if (i == ix-1)
		{
			LABrProjector Proj = Proj0;
			Proj.Params.family = 5 - Proj.Params.family;
			for (LAINT k=0; k<Proj.SubMatrix.get_rows(); k++)
				swap(Proj.SubMatrix(k,0), Proj.SubMatrix(k,1));
			Brancher.Projectors.push_back(Proj);
		}
		else
		{
			Brancher.Projectors.push_back(Proj0);
		}
	}
	
	Brancher.Setup();
	return Brancher;
}


// Subalgebra #ix of brancher Brn0 gets split into two A(1)'s if it is D(2)
LABrancher BrancherSplitD2(LABrancher &Brn0, LAINT ix)
{
	LABrancher Brancher;
	if (ix <= 0 || ix > Brn0.Projectors.size()) return Brancher;
	LABrProjector &Proj01 = Brn0.Projectors[ix-1];
	if (Proj01.Params.family != 4) return Brancher;
	if (Proj01.Params.rank != 2) return Brancher;
	
	// Same starting algebra, U(1)'s
	Brancher.OrigAlg.Params = Brn0.OrigAlg.Params;
	Brancher.U1Indices = Brn0.U1Indices;
	Brancher.U1SrcVecs = Brn0.U1SrcVecs;

	// Do the projection matrices
	for (LAINT i=0; i<Brn0.Projectors.size(); i++)
	{
		LABrProjector &Proj0 = Brn0.Projectors[i];
		if (i == ix-1)
		{
			LABrProjector Proj1, Proj2;
			Proj1.Params.family = 1;
			Proj2.Params.family = 1;
			Proj1.Params.rank = 1;
			Proj2.Params.rank = 1;
			Proj1.SubMatrix.resize(Proj0.SubMatrix.get_rows(),1);
			Proj2.SubMatrix.resize(Proj0.SubMatrix.get_rows(),1);
			for (LAINT k=0; k<Proj0.SubMatrix.get_rows(); k++)
			{
				Proj1.SubMatrix(k,0) = Proj0.SubMatrix(k,0);
				Proj2.SubMatrix(k,0) = Proj0.SubMatrix(k,1);
			}
			Brancher.Projectors.push_back(Proj1);
			Brancher.Projectors.push_back(Proj2);
		}
		else
		{
			Brancher.Projectors.push_back(Proj0);
		}
	}
	
	Brancher.Setup();
	return Brancher;
}


// Subalgebras #ix1 and #ix2 of brancher Brn0 get joined into D(2)
// if they are A(1)/B(1)/C(1)'s.
LABrancher BrancherJoin2A1(LABrancher &Brn0, LAINT ix1, LAINT ix2)
{
	LABrancher Brancher;
	if (ix2 == ix1) return Brancher;
	if (ix1 <= 0 || ix1 > Brn0.Projectors.size()) return Brancher;
	if (ix2 <= 0 || ix2 > Brn0.Projectors.size()) return Brancher;
	LABrProjector &Proj01 = Brn0.Projectors[ix1-1];
	LABrProjector &Proj02 = Brn0.Projectors[ix2-1];
	if (Proj01.Params.rank != 1) return Brancher;
	if (Proj02.Params.rank != 1) return Brancher;
	
	// Same starting algebra, U(1)'s
	Brancher.OrigAlg.Params = Brn0.OrigAlg.Params;
	Brancher.U1Indices = Brn0.U1Indices;
	Brancher.U1SrcVecs = Brn0.U1SrcVecs;
	
	// Do the projection matrices
	LABrProjector Proj1;
	Proj1.Params.family = 4;
	Proj1.Params.rank = 2;
	Proj1.SubMatrix.resize(Proj01.SubMatrix.get_rows(),2);
	for (LAINT i=0; i<Brn0.Projectors.size(); i++)
	{
		LABrProjector &Proj0 = Brn0.Projectors[i];
		if (i == ix1-1)
		{
			for (LAINT k=0; k<Proj0.SubMatrix.get_rows(); k++)
				Proj1.SubMatrix(k,0) = Proj0.SubMatrix(k,0);
		}
		else if (i == ix2-1)
		{
			for (LAINT k=0; k<Proj0.SubMatrix.get_rows(); k++)
				Proj1.SubMatrix(k,1) = Proj0.SubMatrix(k,0);
		}
		else
		{
			Brancher.Projectors.push_back(Proj0);
		}
	}
	Brancher.Projectors.push_back(Proj1);
	
	Brancher.Setup();
	return Brancher;
}


// Makes conjugates of subalgebras  of brancher Brn0 with indexes in cjixs
// different for A(n), n>1, D(n), and E(6). Not a true conjugate for D(2n),
// but exchanged anyway.
LABrancher BrancherConjugate(LABrancher &Brn0, vector<LAINT> &cjixs)
{
	LABrancher Brancher;
	
	// Same starting algebra, U(1)'s
	Brancher.OrigAlg.Params = Brn0.OrigAlg.Params;
	Brancher.U1Indices = Brn0.U1Indices;
	Brancher.U1SrcVecs = Brn0.U1SrcVecs;
	
	// Do the projection matrices
	for (LAINT i=0; i<Brn0.Projectors.size(); i++)
	{
		LABrProjector &Proj0 = Brn0.Projectors[i];
		
		vector<LAINT>::iterator cxiter =
			find(cjixs.begin(), cjixs.end(), i+1);
		if (cxiter != cjixs.end())
		{
			LABrProjector Proj = Proj0;
			
			// Find the conjugate (or pseudo-conjugate / interchange for D(even))
			if (Proj.Params.family == 1)
			{
				// A(n)
				LAINT rnk = Proj.Params.rank;
				LAINT hfrnk = rnk/2;
				for (LAINT k=0; k<Proj.SubMatrix.get_rows(); k++)
				{
					for (LAINT m=0; m<hfrnk; m++)
						swap(Proj.SubMatrix(k,m),Proj.SubMatrix(k,rnk-m-1));
				}
			}
			else if (Proj.Params.family == 4)
			{
				// D(n)
				LAINT rnk = Proj.Params.rank;
				for (LAINT k=0; k<Proj.SubMatrix.get_rows(); k++)
				{
					swap(Proj.SubMatrix(k,rnk-2),Proj.SubMatrix(k,rnk-1));
				}
			}
			else if (Proj.Params.family == 5 && Proj.Params.rank == 6)
			{
				// E(6)
				for (LAINT k=0; k<Proj.SubMatrix.get_rows(); k++)
				{
					swap(Proj.SubMatrix(k,0),Proj.SubMatrix(k,4));
					swap(Proj.SubMatrix(k,1),Proj.SubMatrix(k,3));
				}
			}
			
			Brancher.Projectors.push_back(Proj);
		}
		else
		{
			Brancher.Projectors.push_back(Proj0);
		}
	}
	
	Brancher.Setup();
	return Brancher;
}


// Subalgebra #ix gets conjugated with newrts specifiying
// the new roots' location if it is D(4).
// newrts is 1-based with length 4 with the second being 2
LABrancher BrancherConjgD4(LABrancher &Brn0, LAINT ix, vector<LAINT> &newrts)
{
	LABrancher Brancher;
	if (ix <= 0 || ix > Brn0.Projectors.size()) return Brancher;
	LABrProjector &Proj01 = Brn0.Projectors[ix-1];
	if (Proj01.Params.family != 4) return Brancher;
	if (Proj01.Params.rank != 4) return Brancher;
	
	// Check to see if newrts has all of 1,2,3,4,
	// with the second one being 2
	if (newrts.size() != 4) return Brancher;
	if (newrts[1] != 2) return Brancher;
	vector<LAINT> nrsrt = newrts;
	sort(nrsrt.begin(),nrsrt.end());
	for (LAINT i=0; i<4; i++)
		if (nrsrt[i] != (i+1)) return Brancher;
	
	// Same starting algebra, U(1)'s
	Brancher.OrigAlg.Params = Brn0.OrigAlg.Params;
	Brancher.U1Indices = Brn0.U1Indices;
	Brancher.U1SrcVecs = Brn0.U1SrcVecs;
	
	// Do the projection matrices
	for (LAINT i=0; i<Brn0.Projectors.size(); i++)
	{
		LABrProjector &Proj0 = Brn0.Projectors[i];
		if (i == ix-1)
		{
			LABrProjector Proj = Proj0;
			vector<BrSubMatEntry> MatRow(4);
			for (LAINT k=0; k<Proj.SubMatrix.get_rows(); k++)
			{
				for (LAINT m=0; m<4; m++)
					MatRow[m] = Proj.SubMatrix(k,m);
				for (LAINT m=0; m<4; m++)
					Proj.SubMatrix(k,m) = MatRow[newrts[m]-1];
			}
			Brancher.Projectors.push_back(Proj);
		}
		else
		{
			Brancher.Projectors.push_back(Proj0);
		}
	}
	
	Brancher.Setup();
	return Brancher;
}


// Puts the subalgebras of brancher Brn0 into the order specified in neword.
LABrancher BrancherRearrange(LABrancher &Brn0, vector<LAINT> &neword)
{
	LABrancher Brancher;
	if (neword.size() != Brn0.Projectors.size()) return Brancher;
	vector<LAINT> nwosrt = neword;
	sort(nwosrt.begin(),nwosrt.end());
	for (LAINT i=0; i<Brn0.Projectors.size(); i++)
		if (nwosrt[i] != (i+1)) return Brancher;
	
	// Same starting algebra, U(1)'s
	Brancher.OrigAlg.Params = Brn0.OrigAlg.Params;
	Brancher.U1Indices = Brn0.U1Indices;
	Brancher.U1SrcVecs = Brn0.U1SrcVecs;
	
	for (vector<LAINT>::iterator nwiter = neword.begin();
		nwiter != neword.end(); nwiter++)
	{
		LABrProjector &Proj0 = Brn0.Projectors[(*nwiter)-1];
		Brancher.Projectors.push_back(Proj0);
	}
	
	Brancher.Setup();
	return Brancher;
}


// Returns branching to the original algebra
LABrancher SubalgSelf(const LieAlgebraParams &OrigAlgParams)
{
	LABrancher Brancher;

	Brancher.OrigAlg.Params = OrigAlgParams;
	
	LABrProjector Proj;
	Proj.Params = OrigAlgParams;
	
	// The projection matrix is the identity matrix here
	LAINT rank = OrigAlgParams.rank;
	Proj.SubMatrix.resize(rank);
	Proj.SubMatrix.fill(0);
	for (LAINT k=0; k<rank; k++)
		Proj.SubMatrix(k,k) = 1;
	
	Brancher.Projectors.push_back(Proj);
	
	Brancher.Setup();
	return Brancher;
}



