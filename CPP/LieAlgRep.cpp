/*
	Implements the representation functions
*/

#include "LieAlgRep.h"

using LAXINT_VECTOR = std::vector<LAXINT>;

TDINT TotalDegen(const LieAlgebra &la, const LAINT *MaxWeights)
{
	if (!la.IsValid) return TDINT(0);
	
	// Weyl's celebrated formula
	
	// We need the product of <pr, mr + (1/2)prs> / <pr, (1/2)prs>
	// Calculated in integerized fashion
	// <pr, 2*mw.invctnnum + invctnden*prs> / <pr, invctnden*prs>
	
	LAXINT_VECTOR MaxRoots(la.rank);
	mul_vm(&MaxRoots[0],MaxWeights,la.InvCtnNum);
	mulby_sv(MaxRoots,LAXINT(2));
	
	LAXINT_VECTOR ScaledPRSum(la.rank);
	mul_sv(ScaledPRSum,la.InvCtnDen,la.PosRootSum);
	
	LAXINT_VECTOR ScaledPRWTSum(la.rank);
	add_vv(ScaledPRWTSum,MaxRoots,ScaledPRSum);
	
	LAXINT_VECTOR MetSclPRSum(la.rank);
	mul_mv(MetSclPRSum,la.Metric,ScaledPRSum);
	
	LAXINT_VECTOR MetSclPRWTSum(la.rank);
	mul_mv(MetSclPRWTSum,la.Metric,ScaledPRWTSum);
	
	// The result
	// Fraction<TDINT> ntf = TDINT(1);
	mpq_class ntf(1);
	for (int i=0; i<la.PosRoots.get_rows(); i++)
	{
		LAXINT prnum, prden;
		const LAINT_MATRIX_ROW_CONST PosRoot(la.PosRoots,i);
		mul_vv(prnum,&PosRoot[0],MetSclPRWTSum);
		mul_vv(prden,&PosRoot[0],MetSclPRSum);
		// Fraction<TDINT> pr(prnum,prden);
		// ntf *= pr;
		mpq_class pr(prnum,prden);
		pr.canonicalize();
		ntf *= pr;
	}
	// return ntf.num();
	TDINT res = ntf.get_num();
	return res;
}

// Conjugate of the irrep
// Returns whether or not self-conjugate
bool RepConjugate(const LieAlgebraParams &AlgParams, const LAINT *MaxWeights, LAINT *ConjgMaxWts)
{
	LAINT family = AlgParams.family; LAINT rank = AlgParams.rank;
	if (family == 1)
		std::reverse_copy(MaxWeights, MaxWeights+rank, ConjgMaxWts);
	else if (family == 4 && ((rank % 2) != 0))
	{
		std::copy(MaxWeights, MaxWeights+rank-2, ConjgMaxWts);
		std::reverse_copy(MaxWeights+rank-2, MaxWeights+rank, ConjgMaxWts+rank-2);
	}
	else if (family == 5 && rank == 6)
	{
		std::reverse_copy(MaxWeights, MaxWeights+5, ConjgMaxWts);
		ConjgMaxWts[5] = MaxWeights[5];
	}
	else
		std::copy(MaxWeights, MaxWeights+rank, ConjgMaxWts);
	
	return VecEqual(MaxWeights,ConjgMaxWts,rank);
}

// Height of the irrep
LAINT RepHeight(const LieAlgebraParams &AlgParams, const LAINT *MaxWeights)
{
	LAINT family = AlgParams.family; LAINT rank = AlgParams.rank;
	LAINT h = 0;
	
	if (family == 1)
	{
		for (int k=0; k<rank; k++)
			h += (k+1)*(rank-k)*MaxWeights[k];
	}
	else if (family == 2)
	{
		for (int k=0; k<rank-1; k++)
			h += (k+1)*(2*rank-k)*MaxWeights[k];
		h += (rank*(rank+1)/2)*MaxWeights[rank-1];
	}
	else if (family == 3)
	{
		for (int k=0; k<rank; k++)
			h += (k+1)*(2*rank-k-1)*MaxWeights[k];
	}
	else if (family == 4)
	{
		for (int k=0; k<rank-2; k++)
			h += (k+1)*(2*rank-k-2)*MaxWeights[k];
		h += (rank*(rank-1)/2)*(MaxWeights[rank-2]+MaxWeights[rank-1]);
	}
	else if (family == 5 && rank == 6)
	{
		const LAINT Coeffs[6] = {16, 30, 42, 30, 16, 22};
		mul_vv(h, MaxWeights, Coeffs, rank);
	}
	else if (family == 5 && rank == 7)
	{
		const LAINT Coeffs[7] = {34, 66, 96, 75, 52, 27, 49};
		mul_vv(h, MaxWeights, Coeffs, rank);
	}
	else if (family == 5 && rank == 8)
	{
		const LAINT Coeffs[8] = {92, 182, 270, 220, 168, 114, 58, 136};
		mul_vv(h, MaxWeights, Coeffs, rank);
	}
	else if (family == 6 && rank == 4)
	{
		const LAINT Coeffs[4] = {22, 42, 30, 16};
		mul_vv(h, MaxWeights, Coeffs, rank);
	}
	else if (family == 7 && rank == 2)
	{
		const LAINT Coeffs[2] = {10, 6};
		mul_vv(h, MaxWeights, Coeffs, rank);
	}
	
	return h;
}


static void PushConserved(LAINT_VECTOR &Conserved, LAINT mod, LAINT q)
{
	Conserved.push_back(mod);
	Conserved.push_back(q % mod);
}


// Conserved-quantity values: sets of (modulus, value)
void RepConserved(const LieAlgebraParams &AlgParams, const LAINT *MaxWeights, LAINT_VECTOR &Conserved)
{
	LAINT family = AlgParams.family; LAINT rank = AlgParams.rank;
	Conserved.clear();
	
	if (family == 1)
	{
		LAINT q = 0;
		for (int i=0; i<rank; i++)
			q += (i+1)*MaxWeights[i];
		PushConserved(Conserved,rank+1,q);
	}
	else if (family == 2)
	{
		PushConserved(Conserved,2,MaxWeights[rank-1]);
	}
	else if (family == 3)
	{
		LAINT q = 0;
		for (int i=0; i<rank; i+=2)
			q += MaxWeights[i];
		PushConserved(Conserved,2,q);
	}
	else if (family == 4)
	{
		PushConserved(Conserved,2,MaxWeights[rank-2]+MaxWeights[rank-1]);
		LAINT q = 0;
		for (int i=0; i<rank-2; i+=2)
			q += MaxWeights[i];		
		if ((rank % 2) == 0)
		{
			PushConserved(Conserved,2, q+MaxWeights[rank-2]);
			PushConserved(Conserved,2, q+MaxWeights[rank-1]);
		}
		else
		{
			PushConserved(Conserved,4, 2*q+MaxWeights[rank-2]+3*MaxWeights[rank-1]);
		}
	}
	else if (family == 5 && rank == 6)
	{
		LAINT q = MaxWeights[0] + 2*MaxWeights[1] + MaxWeights[3] + 2*MaxWeights[4];
		PushConserved(Conserved,3,q);
	}
	else if (family == 5 && rank == 7)
	{
		LAINT q = MaxWeights[3] + MaxWeights[5] + MaxWeights[6];
		PushConserved(Conserved,2,q);
	}
}


void RepProperties::Use(const LieAlgebraParams &AlgParams_, const LAINT *MaxWts_)
{
	AlgParams = AlgParams_;
	LAINT rank = AlgParams.rank;
	MaxWts.resize(rank);
	ConjgMaxWts.resize(rank);
	copy(MaxWts_, MaxWts_+rank, MaxWts.begin());
	bool IsReal = RepConjugate(AlgParams,&MaxWts[0],&ConjgMaxWts[0]);
	RepConserved(AlgParams,&MaxWts[0],Conserved);
	Height = RepHeight(AlgParams,&MaxWts[0]);
	if (IsReal)
	{
		if ((Height % 2) == 0)
			Reality = REP_REALITY_REAL;
		else
			Reality = REP_REALITY_PSEUDOREAL;
	}
	else
		Reality = REP_REALITY_COMPLEX;
	Size = TotalDegen(GetLieAlgebra(AlgParams),&MaxWts[0]);
}


// Representation data object

class CompareRootsBySum
{
public:
	LAINT_MATRIX *RtsPtr;
	
	// True if f(ix1) < f(ix2), false otherwise
	bool operator() (size_t ix1, size_t ix2);
};

bool CompareRootsBySum::operator() (size_t ix1, size_t ix2)
{
	LAINT_MATRIX &Rts = *RtsPtr;
	LAINT_MATRIX_ROW R1(Rts,ix1), R2(Rts,ix2);
	size_t n = Rts.get_cols();
	long sres1, sres2;
	sum(sres1,&R1[0],n);
	sum(sres2,&R2[0],n);
	if (sres1 > sres2) return true;
	if (sres1 < sres2) return false;
	return VecLessThan(&R1[0], &R2[0], n);
}

void LieAlgRep::AddRoot(LAINT Degen, const LAINT *Root, const LAINT *Weight)
{
	Degens.push_back(Degen);
	Roots.AppendVector(Root);
	Weights.AppendVector(Weight);
}

void LieAlgRep::SortIndices(SIZE_T_VECTOR &Indxs) const
{	
	// Index sort by root sum; largest to smallest
	size_t n = Degens.size();
	Indxs.resize(n);
	for (size_t i=0; i<n; i++)
		Indxs[i] = i;
	CompareRootsBySum RootCompare;
	RootCompare.RtsPtr = (LAINT_MATRIX *)&Roots;
	sort(Indxs.begin(), Indxs.end(), RootCompare); 
}

void LieAlgRep::Export(LieAlgRep &rep, SIZE_T_VECTOR &Indxs) const
{
	rep.set_vlen(vlen);
	for (auto ix: Indxs)
	{
		LAINT Degen = Degens[ix];
		if (Degen == 0) continue; // Zero means absent - no need to export it
		const LAINT_MATRIX_ROW_CONST Rt(Roots,ix);
		const LAINT_MATRIX_ROW_CONST Wt(Weights,ix);
		rep.Degens.push_back(Degen);
		rep.Roots.AppendVector(&Rt[0]);
		rep.Weights.AppendVector(&Wt[0]);
	}
}

void LieAlgRep::Export(LieAlgRep &Rep) const
{
	SIZE_T_VECTOR Indxs;
	SortIndices(Indxs);
	Export(Rep,Indxs);
}

// Expanded rep object

// Adds the root (degen, rt vec, wt vec) if not already present,
// or its count (degen) if it is.
size_t LieAlgRepBuilder::AddRootOrCount(LAINT Degen, const LAINT *Root, const LAINT *Weight)
{
	std::pair<bool,size_t> ret = Indexer.AppendVector(Root);
	if (ret.first)
	{
		// Add its count
		Degens[ret.second] += Degen;
	}
	else
	{
		// Add it
		AddRoot(Degen,Root,Weight);
	}
	return ret.second;
}

void LieAlgRepBuilder::Import(LieAlgRep &Rep, bool Append)
{
	if (!Append) set_vlen(Rep.vlen);
	for (size_t i=0; i<Rep.Roots.get_rows(); i++)
		AddRootOrCount(Rep.Degens[i],&Rep.Roots(i,0),&Rep.Weights(i,0));
}


// The main calculation of a (semi)simple-algebra irrep

static void RepRootVectors(LieAlgRepBuilder &bld, const LieAlgebra &la, const LAINT *MaxWts)
{
	// Initial root
	LAINT_VECTOR MaxRts(la.rank);
	mul_vm(&MaxRts[0],MaxWts,la.InvCtnNum);
	bld.AddRootOrCount(0, &MaxRts[0], MaxWts);
	
	// Using LAINT instead of bool because in the STL,
	// bool gets turned into packed bits,
	// something that you can't use pointers or refs on
	LAINT_MATRIX WhichWay(1,la.rank);
	WhichWay.fill(true);
	
	LAINT_VECTOR NewRoot(la.rank), NewWeight(la.rank);
	LAINT_VECTOR NewWW(la.rank);
	
	// Find the next root until one cannot find any more
	size_t RtIndx = 0;
	while (RtIndx < bld.Roots.get_rows())
	{
		LAINT_MATRIX_ROW ThisRoot(bld.Roots,RtIndx);
		LAINT_MATRIX_ROW ThisWeight(bld.Weights,RtIndx);
		LAINT_MATRIX_ROW ThisWW(WhichWay,RtIndx);
		
		// Advance in each direction, if possible
		for (int i=0; i<la.rank; i++)
		{
			for (int j=0; j<la.rank; j++)
				NewWW[j] = (j != i);
			if (ThisWW[i])
			{
				for (LAINT j=1; j<=ThisWeight[i]; j++)
				{
					// Calculate roots and weights in parallel,
					// to avoid repeated matrix.vector calculations
					copy(ThisRoot.begin(), ThisRoot.end(), NewRoot.begin());
					NewRoot[i] -= j*la.InvCtnDen;
					copy(ThisWeight.begin(), ThisWeight.end(), NewWeight.begin());
					for (int k=0; k<la.rank; k++)
						NewWeight[k] -= j*la.Cartan(i,k);
					
					// Add the root if possible
					size_t NumRoots = bld.Roots.get_rows();
					size_t NewIndx = bld.AddRootOrCount(0,NewRoot,NewWeight);
					if (NewIndx < NumRoots)
					{
						// Don't add it, and mark off its direction
						// as not to be followed
						LAINT_MATRIX_ROW WWRow(WhichWay,RtIndx);
						for (int k=0; k<la.rank; k++)
							WWRow[k] &= NewWW[k];
					}
					else
					{
						// Add it
						WhichWay.AppendVector(NewWW);
					}
				}
			}
		}
		// All done; get ready to do the next one
		RtIndx++;
	}
}


static void RepRootDegens(LieAlgRepBuilder &bld, const LieAlgebra &la, LieAlgRep &rep)
{
	// Freudenthal's recurrence
	
	// Needs the max-to-min sorting order,
	// because the algorithm uses the highest weight.
	SIZE_T_VECTOR Indxs;
	bld.SortIndices(Indxs);
	
	// Scale up the positive roots and their sum
	LAINT_MATRIX PosRoots = la.PosRoots;
	mulby_sm(PosRoots,la.InvCtnDen);
	LAINT_VECTOR PosRootSum = la.PosRootSum;
	mulby_sv(PosRootSum,la.InvCtnDen);
	
	// The main algorithm
	LAINT_MATRIX_ROW MaxRt(bld.Roots,Indxs[0]);	
	bld.Degens[Indxs[0]] = 1;
	LAINT_VECTOR BkRt(la.rank);
	std::vector<LAXINT> mrpd(la.rank), mrm(la.rank);
	bool IsFirst = true;
	for (auto ix: Indxs)
	{
		if (IsFirst) {IsFirst = false; continue;}
		LAINT_MATRIX_ROW Rt(bld.Roots,ix);
		LAXINT nx = 0;
		for (size_t i=0; i<PosRoots.get_rows(); i++)
		{
			LAINT_MATRIX_ROW ShtRt(PosRoots,i);
			add_vv(BkRt,&Rt[0],&ShtRt[0]);
			while (true)
			{
				std::pair<bool,size_t> ret = bld.Indexer.VectorIndex(&BkRt[0]);
				if (ret.first == false) break;
				LAXINT nbs;
				mul_vmv(nbs,BkRt,la.Metric,&ShtRt[0]);
				nx += bld.Degens[ret.second]*nbs;
				addto_vv(BkRt,&ShtRt[0]);
			}
		}
		add_vv(mrpd, &MaxRt[0], &Rt[0]);
		addto_vv(mrpd, &PosRootSum[0]);
		sub_vv(mrm, &MaxRt[0], &Rt[0]);
		LAXINT nd;
		mul_vmv(nd, mrpd, la.Metric, mrm);
		bld.Degens[ix] = (2*nx)/nd;
	}
	
	// Use the sorted order to export
	bld.Export(rep,Indxs);
}

bool LieAlgebraRepDirect(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts)
{
	if (!la.IsValid) return false;
	for (size_t i=0; i<la.rank; i++)
		if (MaxWts[i] < 0) return false;
	LieAlgRepBuilder bld(la.rank);
	RepRootVectors(bld,la,MaxWts);
	RepRootDegens(bld,la,rep);
	return true;
}


// Weyl orbits

static bool WeylOrbitForDomWt(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts)
{
	if (!la.IsValid) return false;
	
	LieAlgRepBuilder bld(la.rank);
	
	// Initial root
	LAINT_VECTOR MaxRts(la.rank);
	mul_vm(&MaxRts[0],MaxWts,la.InvCtnNum);
	bld.AddRootOrCount(0, &MaxRts[0], MaxWts);
	
	LAINT_VECTOR NewRoot(la.rank), NewWeight(la.rank);
	LAINT_VECTOR MetRoot(la.rank);
	
	// Find the next root until one cannot find any more
	size_t RtIndx = 0;
	while (RtIndx < bld.Roots.get_rows())
	{
		LAINT_MATRIX_ROW ThisRoot(bld.Roots,RtIndx);
		
		mul_mv(MetRoot,la.Metric,&ThisRoot[0]);
		
		// Advance in each direction, if possible
		for (int i=0; i<la.rank; i++)
		{
			LAINT mr = MetRoot[i];
			if (mr > 0)
			{
				// Add this root
				copy(ThisRoot.begin(), ThisRoot.end(), NewRoot.begin());
				NewRoot[i] -= (2*mr)/la.Metric(i,i);
				mul_vm(NewWeight, NewRoot, la.Cartan);
				divby_sv(NewWeight, la.InvCtnDen);
				bld.AddRootOrCount(0,NewRoot,NewWeight);
			}
		}

		// All done; get ready to do the next one
		RtIndx++;
	}
	
	fill(bld.Degens.begin(), bld.Degens.end(), 1);
	bld.Export(rep);
	return true;
}


static bool DomWtFromRoot(LAINT *MaxWts, const LieAlgebra &la, const LAINT *Root)
{
	if (!la.IsValid) return false;
	
	// Set up the root
	LAINT_VECTOR ThisRoot(la.rank), ThisWeight(la.rank);
	LAINT_VECTOR NewRoot(la.rank), NewWeight(la.rank);
	LAINT_VECTOR MetRoot(la.rank);
	
	copy(Root, Root+la.rank, ThisRoot.begin());
	
	while(true)
	{
		// Find the weight. Is it a dominant one?
		mul_vm(ThisWeight, ThisRoot, la.Cartan);
		divby_sv(ThisWeight, la.InvCtnDen);
		LAINT minwt;
		index_min_v(minwt, ThisWeight);
		if (minwt >= 0) break;
		
		mul_mv(MetRoot,la.Metric,ThisRoot);
		
		LAINT mr;
		size_t i = index_min_v(mr, MetRoot);
		ThisRoot[i] -= (2*mr)/la.Metric(i,i);
	}
	
	copy(ThisWeight.begin(), ThisWeight.end(), MaxWts);
	return true;
}


template<typename N> void MxWtToYD(N *YDPtr, const N *MWPtr, size_t n)
{
	YDPtr[n-1] = MWPtr[n-1];
	for (int i=n-2; i>=0; i--)
		YDPtr[i] = YDPtr[i+1] + MWPtr[i];
}

template<typename N> void Accumulate(N *AccSum, const N *Orig, size_t n)
{
	AccSum[0] = Orig[0];
	for (int i=1; i<n; i++)
		AccSum[i] = AccSum[i-1] + Orig[i];
}


template<typename N> void Permutations(Matrix<N> &Mat, const N *Syms, size_t n)
{
	Mat.resize(0,n);
	std::vector<N> SymVec(n);
	std::copy(Syms,Syms+n,SymVec.begin());
	std::sort(SymVec.begin(),SymVec.end());
	Mat.AppendVector(SymVec);
	while(next_permutation(SymVec.begin(),SymVec.end()))
		Mat.AppendVector(SymVec);
}

template<typename N> void Permutations(Matrix<N> &Mat, const std::vector<N> &Syms)
	{Permutations(Mat, &Syms[0], Syms.size());}

template<typename N> void PermsForList(Matrix<N> &Mat, const Matrix<N> &Syms)
{
	Mat.resize(0,Syms.get_cols());
	Matrix<N> SymPerms;
	for (int i=0; i<Syms.get_rows(); i++)
	{
		const ConstMatrixRow<N> SymRow(Syms,i);
		Permutations(SymPerms, &SymRow[0], Syms.get_cols());
		Mat.AppendMatrix(SymPerms);
	}
}

template<typename N> void AddSigns(Matrix<N> &Mat, const N *Syms, size_t n)
{
	Mat.resize(0,n);
	std::vector<N> SymVec(n);
	std::vector<bool> SameAbsVal(n-1);
	
	for (int i=0; i<n; i++)
	{
		N S = Syms[i];
		SymVec[i] = - ((S >= 0) ? S : -S);
	}
	
	for (int i=0; i<n-1; i++)
		SameAbsVal[i] = (SymVec[i] == SymVec[i+1]);
	
	while(true)
	{
		// Add if terms with the same absolute value are not out of order
		// Assumes already-sorted input
		bool AddToMat = true;
		for (int i=0; i<n-1; i++)
		{
			if (SameAbsVal[i] && (SymVec[i] < SymVec[i+1]))
			{
				AddToMat = false;
				break;
			}
		}
		if (AddToMat)
			Mat.AppendVector(SymVec);
		
		// Look for the next sign combination
		bool DidAdvance = false;
		for (int i=n-1; i>=0; i--)
		{
			if (SymVec[i] < 0)
			{
				SymVec[i] *= -1;
				for (int j=i+1; j<n; j++)
					SymVec[j] *= -1;
				DidAdvance = true;
				break;
			}
		}
		if (!DidAdvance) break;
	}
}

template<typename N> void AddSigns(Matrix<N> &Mat, const std::vector<N> &Syms)
	{AddSigns(Mat, &Syms[0], Syms.size());}

template<typename N> void SelectParity(Matrix<N> &SelRts, const Matrix<N> &OrigRts, N Val)
{
	if (Val == 0)
	{
		SelRts = OrigRts;
		return;
	}
	
	SelRts.resize(0,OrigRts.get_cols());
	
	for (int k=0; k<OrigRts.get_rows(); k++)
	{
		const ConstMatrixRow<N> Root(OrigRts,k);
		N SFac = Val;
		for (int i=0; i<Root.size(); i++)
			if (Root[i] < 0) SFac *= -1;
			else if (Root[i] == 0) SFac = 0;
		if (SFac > 0)
			SelRts.AppendVector(&Root[0]);
	}
}

static bool WeylOrbitForDomWtExplicit(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts)
{
	if (!la.IsValid) return false;
	if (la.family > 4) return false;
	
	rep.set_vlen(la.rank);
	LAINT_VECTOR Root(la.rank), Weight(la.rank);
	
	switch(la.family)
	{
	// A(n)
	case 1:
	{
		LAINT_VECTOR MWXtnd(la.rank+1);
		copy(MaxWts,MaxWts+la.rank,MWXtnd.begin());
		MWXtnd[la.rank] = 0;
		LAINT_VECTOR YDXtnd(la.rank+1);
		MxWtToYD(&YDXtnd[0],&MWXtnd[0],la.rank+1);
		LAINT YDSum = 0;
		for (int i=0; i<=la.rank; i++)
			YDSum += YDXtnd[i];
		for (int i=0; i<=la.rank; i++)
			YDXtnd[i] = (la.rank+1)*YDXtnd[i] - YDSum;
		LAINT_MATRIX RtXtnd;
		Permutations(RtXtnd, &YDXtnd[0], la.rank+1);
		for (int k=0; k<RtXtnd.get_rows(); k++)
		{
			Accumulate(&Root[0],&RtXtnd(k,1),la.rank);
			rep.AddRoot(1, Root, Weight);
		}
		break;
	}
	// B(n)
	case 2:
	{
		LAINT_VECTOR MWX(la.rank);
		copy(MaxWts,MaxWts+la.rank,MWX.begin());
		mulby_sv(&MWX[0],2,la.rank-1);
		LAINT_VECTOR YDX(la.rank);
		MxWtToYD(&YDX[0],&MWX[0],la.rank);
		LAINT_MATRIX SgndYD;
		AddSigns(SgndYD,&YDX[0],la.rank);
		LAINT_MATRIX Rts;
		PermsForList(Rts,SgndYD);
		for (int k=0; k<Rts.get_rows(); k++)
		{
			Accumulate(&Root[0],&Rts(k,0),la.rank);
			rep.AddRoot(1, Root, Weight);
		}
		break;
	}
	// C(n)
	case 3:
	{
		LAINT_VECTOR MWX(la.rank);
		copy(MaxWts,MaxWts+la.rank,MWX.begin());
		LAINT_VECTOR YDX(la.rank);
		MxWtToYD(&YDX[0],&MWX[0],la.rank);
		LAINT_MATRIX SgndYD;
		AddSigns(SgndYD,&YDX[0],la.rank);
		LAINT_MATRIX Rts;
		PermsForList(Rts,SgndYD);
		for (int k=0; k<Rts.get_rows(); k++)
		{
			Accumulate(&Root[0],&Rts(k,0),la.rank);
			mulby_sv(&Root[0],2,la.rank-1);
			rep.AddRoot(1, Root, Weight);
		}
		break;
	}
	// D(n)
	case 4:
	{
		LAINT_VECTOR MWX(la.rank);
		copy(MaxWts,MaxWts+la.rank,MWX.begin());
		mulby_sv(&MWX[0],2,la.rank-2);
		LAINT endm2 = MWX[la.rank-2];
		LAINT endm1 = MWX[la.rank-1];
		MWX[la.rank-2] = endm2 + endm1;
		LAINT_VECTOR YDX(la.rank);
		MxWtToYD(&YDX[0],&MWX[0],la.rank-1);
		YDX[la.rank-1] = - endm2 + endm1;
		LAINT_MATRIX SgndYD;
		AddSigns(SgndYD,&YDX[0],la.rank);
		LAINT_MATRIX PrtySgndYD;
		SelectParity(PrtySgndYD,SgndYD,YDX[la.rank-1]);
		LAINT_MATRIX Rts;
		PermsForList(Rts,PrtySgndYD);
		for (int k=0; k<Rts.get_rows(); k++)
		{
			LAINT_MATRIX_ROW Rt(Rts,k);
			Rt[la.rank-2] = Rt[la.rank-2] - Rt[la.rank-1];
			Rt[la.rank-1] *= 2;
			Accumulate(&Root[0],&Rt[0],la.rank);
			if ((la.rank % 2) == 0)
				divby_sv(&Root[la.rank-2],2,2);
			else
				mulby_sv(&Root[0],2,la.rank-2);
			rep.AddRoot(1, Root, Weight);
		}
		break;
	}
	}
	
	mul_mm(rep.Weights,rep.Roots,la.Cartan);
	divby_sm(rep.Weights,la.InvCtnDen);
	return true;
}


bool LieAlgebraWeylOrbitGeneral(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts)
{
	return WeylOrbitForDomWt(rep, la, MaxWts);
}


bool LieAlgebraWeylOrbitExplicit(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts)
{
	return WeylOrbitForDomWtExplicit(rep, la, MaxWts);
}


bool LieAlgebraWeylOrbit(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts)
{
	bool rc;
	rc = LieAlgebraWeylOrbitExplicit(rep, la, MaxWts);
	if (rc) return rc;
	rc = LieAlgebraWeylOrbitGeneral(rep, la, MaxWts);
	return rc;
}


// The main calculation of a (semi)simple-algebra irrep
// Only get dominant weights of Weyl orbits here

static void WeylDomWtRepRootVectors(LieAlgRepBuilder &bld, const LieAlgebra &la, const LAINT *MaxWts)
{
	// Initial root
	LAINT_VECTOR MaxRts(la.rank);
	mul_vm(&MaxRts[0],MaxWts,la.InvCtnNum);
	bld.AddRootOrCount(0, &MaxRts[0], MaxWts);
	
	// Scale up the positive roots and their sum
	LAINT_MATRIX PosRoots = la.PosRoots;
	mulby_sm(PosRoots,la.InvCtnDen);
		
	LAINT_VECTOR NewRoot(la.rank), NewWeight(la.rank);
	
	// Find the next root until one cannot find any more
	size_t RtIndx = 0;
	while (RtIndx < bld.Roots.get_rows())
	{
		LAINT_MATRIX_ROW ThisRoot(bld.Roots,RtIndx);
		LAINT_MATRIX_ROW ThisWeight(bld.Weights,RtIndx);
		
		// Advance in each direction, if possible
		for (int i=0; i<PosRoots.get_rows(); i++)
		{
			LAINT_MATRIX_ROW ThisPosRoot(PosRoots, i);
			
			// New root and its weight
			copy(ThisRoot.begin(), ThisRoot.end(), NewRoot.begin());
			subfm_vv(NewRoot, &ThisPosRoot[0]);
			mul_vm(NewWeight,NewRoot,la.Cartan);
			divby_sv(NewWeight,la.InvCtnDen);
			
			// Is it a dominant weight?
			LAINT minwt;
			index_min_v(minwt, NewWeight);
			if (minwt < 0) continue;
			
			// Add the root if possible
			bld.AddRootOrCount(0,NewRoot,NewWeight);
		}
		// All done; get ready to do the next one
		RtIndx++;
	}
}

static void WeylDomWtRepRootDegens(LieAlgRepBuilder &bld, const LieAlgebra &la, LieAlgRep &rep)
{
	// Freudenthal's recurrence
	
	// Needs the max-to-min sorting order,
	// because the algorithm uses the highest weight.
	SIZE_T_VECTOR Indxs;
	bld.SortIndices(Indxs);
	
	// Scale up the positive roots and their sum
	LAINT_MATRIX PosRoots = la.PosRoots;
	mulby_sm(PosRoots,la.InvCtnDen);
	LAINT_VECTOR PosRootSum = la.PosRootSum;
	mulby_sv(PosRootSum,la.InvCtnDen);
	
	// The main algorithm
	LAINT_MATRIX_ROW MaxRt(bld.Roots,Indxs[0]);	
	bld.Degens[Indxs[0]] = 1;
	LAINT_VECTOR BkRt(la.rank);
	LAINT_VECTOR DomRt(la.rank), DomWt(la.rank);
	LAXINT_VECTOR mrpd(la.rank), mrm(la.rank);
	for (auto rit=Indxs.begin(); rit!= Indxs.end(); rit++)
	{
		if (rit == Indxs.begin()) continue;
		size_t ix = *rit;
		LAINT_MATRIX_ROW Rt(bld.Roots,ix);
		LAXINT nx = 0;
		for (size_t i=0; i<PosRoots.get_rows(); i++)
		{
			LAINT_MATRIX_ROW ShtRt(PosRoots,i);
			add_vv(BkRt,&Rt[0],&ShtRt[0]);
			// Dominant root and weight for it
			DomWtFromRoot(&DomWt[0], la, &BkRt[0]);
			mul_vm(DomRt, DomWt, la.InvCtnNum);
			while (true)
			{
				std::pair<bool,size_t> ret = bld.Indexer.VectorIndex(&DomRt[0]);
				if (ret.first == false) break;
				LAXINT nbs;
				mul_vmv(nbs,BkRt,la.Metric,&ShtRt[0]);
				nx += bld.Degens[ret.second]*nbs;
				addto_vv(BkRt,&ShtRt[0]);
				// Dominant root and weight for it
				DomWtFromRoot(&DomWt[0], la, &BkRt[0]);
				mul_vm(DomRt, DomWt, la.InvCtnNum);
			}
		}
		add_vv(mrpd, &MaxRt[0], &Rt[0]);
		addto_vv(mrpd, &PosRootSum[0]);
		sub_vv(mrm, &MaxRt[0], &Rt[0]);
		LAXINT nd;
		mul_vmv(nd, mrpd, la.Metric, mrm);
		bld.Degens[ix] = (2*nx)/nd;
	}
	
	// Use the sorted order to export
	bld.Export(rep,Indxs);
}

bool LieAlgebraRepWeylOrbits(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts)
{
	if (!la.IsValid) return false;
	for (size_t i=0; i<la.rank; i++)
		if (MaxWts[i] < 0) return false;
	LieAlgRepBuilder bld(la.rank);
	WeylDomWtRepRootVectors(bld,la,MaxWts);
	WeylDomWtRepRootDegens(bld,la,rep);
	return true;
}

bool LieAlgebraRepByWOs(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts)
{
	if (!la.IsValid) return false;
	
	LieAlgRepPtr DWPtr = GetRepObject(RO_REP_ORBIT, la, MaxWts);
	if (DWPtr.IsNull()) return false;
	
	LieAlgRep intrep(la.rank);
	
	// Workspace for roots and weights
	LAINT_VECTOR WODWT(la.rank), Root(la.rank), Weight(la.rank);
	
	size_t NDWS = DWPtr->Degens.size();
	for (int k=0; k<NDWS; k++)
	{
		LAINT Degen = DWPtr->Degens[k];
		LAINT_MATRIX_ROW DWT(DWPtr->Weights,k);
		
		// Find the scale factor for scaling down
		LAINT MWScale = 0;
		for (int i=0; i<la.rank; i++)
		{
			LAINT wt = DWT[i];
			if (wt == 0) continue;
			if (wt < 0) wt *= -1;
			if (MWScale < 1)
				MWScale = wt;
			else
				MWScale = GCD(MWScale,wt);
		}
		if (MWScale < 1) MWScale = 1;
		
		// Scale down
		div_sv(&WODWT[0], MWScale, &DWT[0],la.rank);
		
		LieAlgRepPtr WOPtr = GetRepObject(RO_ORBIT, la, &WODWT[0]);
		if (WOPtr.IsNull()) continue;
		
		size_t NEnts = WOPtr->Degens.size();
		for (int i=0; i< NEnts; i++)
		{
			LAINT_MATRIX_ROW WORoot(WOPtr->Roots,i);
			LAINT_MATRIX_ROW WOWeight(WOPtr->Weights,i);
			
			// Restore the scale
			mul_sv(&Root[0], MWScale, &WORoot[0],la.rank);
			mul_sv(&Weight[0], MWScale, &WOWeight[0],la.rank);
			
			intrep.AddRoot(Degen, &Root[0], &Weight[0]);
		}
	}
	
	intrep.Export(rep);
	return true;
}


bool LieAlgebraRepresentation(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts)
{
	return LieAlgebraRepByWOs(rep, la, MaxWts);
}


// For caching the representations

using InlineGetRep = bool (*)(LieAlgRep &rep, const LieAlgebra &la, const LAINT *MaxWts);

static InlineGetRep RepFuncs[NUMBER_OF_REP_OBJ_TYPES] =
{
	LieAlgebraWeylOrbit,
	LieAlgebraRepWeylOrbits,
	LieAlgebraRepresentation
};

struct LieAlgRepParams
{
	LieAlgebraParams AlgParams;
	LAINT_VECTOR MaxWts;
};

struct LieAlgRepLessThan
{
	bool operator() (const LieAlgRepParams &R1, const LieAlgRepParams &R2) const;
};

bool LieAlgRepLessThan::operator() (const LieAlgRepParams &R1, const LieAlgRepParams &R2) const
{
	if (R1.AlgParams.family < R2.AlgParams.family) return true;
	if (R1.AlgParams.family > R2.AlgParams.family) return false;
	if (R1.AlgParams.rank < R2.AlgParams.rank) return true;
	if (R1.AlgParams.rank > R2.AlgParams.rank) return false;
	return VecLessThan(&R1.MaxWts[0], &R2.MaxWts[0], R1.AlgParams.rank);
}

static std::map<LieAlgRepParams, LieAlgRepPtr, LieAlgRepLessThan> LieAlgRepCaches[NUMBER_OF_REP_OBJ_TYPES];

LieAlgRepPtr InvalidLieAlgRep(NULL);

LieAlgRepPtr &GetRepObject(enum RepObjType rotype, const LieAlgebraParams &AlgParams, const LAINT *MaxWts)
{
	std::map<LieAlgRepParams, LieAlgRepPtr, LieAlgRepLessThan> &LieAlgRepCache = LieAlgRepCaches[rotype];
	InlineGetRep RepFunc = RepFuncs[rotype];
	LieAlgRepParams RepParams;
	RepParams.AlgParams = AlgParams;
	RepParams.MaxWts.resize(RepParams.AlgParams.rank);
	copy(MaxWts, MaxWts+AlgParams.rank, RepParams.MaxWts.begin());
	auto it = LieAlgRepCache.find(RepParams);
	if (it == LieAlgRepCache.end())
	{
		std::pair<LieAlgRepParams, LieAlgRepPtr> rec;
		LieAlgRep &Rep = *rec.second;
		if (!RepFunc(Rep,GetLieAlgebra(AlgParams),MaxWts))
			return InvalidLieAlgRep;
		rec.first = RepParams;
		std::pair<std::map<LieAlgRepParams, LieAlgRepPtr, LieAlgRepLessThan>::iterator, bool> ret =
			LieAlgRepCache.insert(rec);
		return ret.first->second;
	}
	else
		return it->second;
}
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, const LieAlgebraParams &AlgParams, LAINT_VECTOR &MaxWts)
	{return GetRepObject(rotype, AlgParams, &MaxWts[0]);}

LieAlgRepPtr &GetRepObject(enum RepObjType rotype, LAINT family, LAINT rank, const LAINT *MaxWts)
{
	LieAlgebraParams AlgParams;
	AlgParams.family = family; AlgParams.rank = rank;
	return GetRepObject(rotype, AlgParams, MaxWts);
}
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, LAINT family, LAINT rank, LAINT_VECTOR &MaxWts)
{
	LieAlgebraParams AlgParams;
	AlgParams.family = family; AlgParams.rank = rank;
	return GetRepObject(rotype, AlgParams, &MaxWts[0]);
}

LieAlgRepPtr &GetRepObject(enum RepObjType rotype, const LieAlgebra &la, const LAINT *MaxWts)
	{return GetRepObject(rotype, la.GetParams(), MaxWts);}
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, const LieAlgebra &la, const LAINT_VECTOR &MaxWts)
	{return GetRepObject(rotype, la.GetParams(), &MaxWts[0]);}

void ClearLieAlgReps()
{
	for (int rotype=0; rotype<NUMBER_OF_REP_OBJ_TYPES; rotype++)
	{
		std::map<LieAlgRepParams, LieAlgRepPtr, LieAlgRepLessThan> &LieAlgRepCache =
			LieAlgRepCaches[rotype];
		LieAlgRepCache.clear();
	}
}


// Find rep from product algebra -- (semi)simple ones and U(1)'s
LAINT LieAlgProduct::get_rank() const
{
	LAINT rank = 0;
	for (const auto &it: ParamList)
		rank += it.rank;
	rank += NumU1s;
	return rank;
}

void LieAlgProdRepresentation(LieAlgRep &rep, const LieAlgProduct &lap,
	enum RepObjType rotype, const LAINT *MaxWts)
{
	LAINT rank = lap.get_rank();
	rep.set_vlen(rank);
	
	// Create the initial root -- a singlet
	LAINT_VECTOR Rt(rank), Wt(rank);
	fill(Rt.begin(), Rt.end(), 0);
	fill(Wt.begin(), Wt.end(), 0);
	size_t nustrt = rank - lap.NumU1s;
	std::copy(MaxWts+nustrt, MaxWts+rank, &Rt[nustrt]);
	std::copy(MaxWts+nustrt, MaxWts+rank, &Wt[nustrt]);
	rep.AddRoot(1,Rt,Wt);
	
	LieAlgRep NewRep; // Temporary for each run
	size_t stix = 0, ndix; // Start index
	
	for (const auto &Params: lap.ParamList)
	{
		ndix = stix + Params.rank;
		
		LieAlgRepPtr &RepAddedPtr = GetRepObject(rotype, Params, MaxWts+stix);
		LieAlgRep &RepAdded = *RepAddedPtr;
		NewRep.set_vlen(rank);
		for (int i=0; i<rep.Degens.size(); i++)
		{
			// Get the entire root/weight from the existing array
			LAINT Degen = rep.Degens[i];
			LAINT_MATRIX_ROW Rtx(rep.Roots,i);
			LAINT_MATRIX_ROW Wtx(rep.Weights,i);
			std::copy(Rtx.begin(), Rtx.end(), Rt.begin());
			std::copy(Wtx.begin(), Wtx.end(), Wt.begin());
			for (int j=0; j<RepAdded.Degens.size(); j++)
			{
				// Copy in and add the new root
				LAINT NewDegen = Degen*RepAdded.Degens[j];
				LAINT_MATRIX_ROW Rta(RepAdded.Roots,j);
				LAINT_MATRIX_ROW Wta(RepAdded.Weights,j);
				std::copy(Rta.begin(), Rta.end(), &Rt[stix]);
				std::copy(Wta.begin(), Wta.end(), &Wt[stix]);
				NewRep.AddRoot(NewDegen,Rt,Wt);
			}
		}
		std::swap(rep,NewRep); // Put new ones into old ones' location
		stix = ndix;
	}
}


// Abstract base class for representation handler
// Subclass for single-algebra and product-of-algebras cases

// Counted list of weight vectors to rep
const LieAlgRepPtr RepHandlerBase::GetRepPtr(enum RepObjType rotype, const LACntdMaxWtList &CWL) const
{
	LieAlgRepBuilder Bld(get_rank());
	
	for (const auto &Entry: CWL)
	{
		LieAlgRepPtr RepPtr = GetRepPtr(rotype, Entry.MaxWts);
		LieAlgRep &Rep = *RepPtr;
		for (size_t i=0; i<Rep.Degens.size(); i++)
			Bld.AddRootOrCount(Entry.Count*Rep.Degens[i],
				&Rep.Roots(i,0), &Rep.Weights(i,0));
	}
	
	LieAlgRepPtr RepPtr;
	LieAlgRep &Rep = *RepPtr;
	Bld.Export(Rep);
	return RepPtr;
}

// Extract irreps from a rep with sort data
// as a counted list of max weights for irreps
// It will be altered as the extraction goes
const LACntdMaxWtList RepHandlerBase::ExtractWts(enum RepObjType rotype, LieAlgRepBuilder &Bld) const
{
	LAINT rank = get_rank();
	LACntdMaxWtList CWList;
	
	while(true)
	{
		// Search for the root with the highest sum
		// Init these values to avoid "wasn't inited" compiler complaints
		size_t MaxIndx = 0;
		LAINT MaxRtSum = 0;
		bool Found = false;
		for (size_t i=0; i<Bld.Degens.size(); i++)
		{
			// Zero-degen entries: skipped over, not deleted,
			// in order to avoid readjusting the indices in Bld's search tree
			if (Bld.Degens[i] == 0) continue;
			
			LAINT RtSum;
			sum(RtSum,&Bld.Roots(i,0),rank);
			if (Found)
			{
				if (RtSum > MaxRtSum)
				{
					MaxIndx = i;
					MaxRtSum = RtSum;
				}
			}
			else
			{
				Found = true;
				MaxIndx = i;
				MaxRtSum = RtSum;
			}
		}
		// If no nonzero ones were found, then exit
		if (!Found) break;
		
		// Get the degen and max wts
		LACntdMaxWtEntry Entry;
		Entry.Count = Bld.Degens[MaxIndx];
		Entry.MaxWts.resize(rank);
		LAINT_MATRIX_ROW MaxWtRow(Bld.Weights,MaxIndx);
		copy(MaxWtRow.begin(),MaxWtRow.end(),Entry.MaxWts.begin());
		
		// Subtract out the rep
		LieAlgRepPtr RepPtr = GetRepPtr(rotype, Entry.MaxWts);
		LieAlgRep &Rep = *RepPtr;
		for (size_t i=0; i<Rep.Degens.size(); i++)
		{
			Bld.AddRootOrCount(-Entry.Count*Rep.Degens[i],
				&Rep.Roots(i,0), &Rep.Weights(i,0));
		}
		
		// All done
		CWList.push_back(Entry);
	}
	
	return CWList;
}
