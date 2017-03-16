/*
	Implements the representation functions
*/

#include "LieAlgRep.h"


TDINT TotalDegen(LieAlgebra &la, const LAINT *MaxWeights)
{
	if (!la.IsValid) return TDINT(0);
	
	// Weyl's celebrated formula
	
	// We need the product of <pr, mr + (1/2)prs> / <pr, (1/2)prs>
	// Calculated in integerized fashion
	// <pr, 2*mw.invctnnum + invctnden*prs> / <pr, invctnden*prs>
	
	vector<LAXINT> MaxRoots(la.rank);
	mul_vm(&MaxRoots[0],MaxWeights,la.InvCtnNum);
	mulby_sv(MaxRoots,LAXINT(2));
	
	vector<LAXINT> ScaledPRSum(la.rank);
	mul_sv(ScaledPRSum,la.InvCtnDen,la.PosRootSum);
	
	vector<LAXINT> ScaledPRWTSum(la.rank);
	add_vv(ScaledPRWTSum,MaxRoots,ScaledPRSum);
	
	vector<LAXINT> MetSclPRSum(la.rank);
	mul_mv(MetSclPRSum,la.Metric,ScaledPRSum);
	
	vector<LAXINT> MetSclPRWTSum(la.rank);
	mul_mv(MetSclPRWTSum,la.Metric,ScaledPRWTSum);
	
	// The result
	// Fraction<TDINT> ntf = TDINT(1);
	mpq_class ntf(1);
	for (int i=0; i<la.PosRoots.get_rows(); i++)
	{
		LAXINT prnum, prden;
		MatrixRow<LAINT> PosRoot(la.PosRoots,i);
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
		reverse_copy(MaxWeights, MaxWeights+rank, ConjgMaxWts);
	else if (family == 4 && ((rank % 2) != 0))
	{
		copy(MaxWeights, MaxWeights+rank-2, ConjgMaxWts);
		reverse_copy(MaxWeights+rank-2, MaxWeights+rank, ConjgMaxWts+rank-2);
	}
	else if (family == 5 && rank == 6)
	{
		reverse_copy(MaxWeights, MaxWeights+5, ConjgMaxWts);
		ConjgMaxWts[5] = MaxWeights[5];
	}
	else
		copy(MaxWeights, MaxWeights+rank, ConjgMaxWts);
	
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


static void PushConserved(vector<LAINT> &Conserved, LAINT mod, LAINT q)
{
	Conserved.push_back(mod);
	Conserved.push_back(q % mod);
}


// Conserved-quantity values: sets of (modulus, value)
void RepConserved(const LieAlgebraParams &AlgParams, const LAINT *MaxWeights, vector<LAINT> &Conserved)
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
	Matrix<LAINT> *RtsPtr;
	
	// True if f(ix1) < f(ix2), false otherwise
	bool operator() (size_t ix1, size_t ix2);
};

bool CompareRootsBySum::operator() (size_t ix1, size_t ix2)
{
	Matrix<LAINT> &Rts = *RtsPtr;
	MatrixRow<LAINT> R1(Rts,ix1), R2(Rts,ix2);
	size_t n = Rts.get_cols();
	long sres1, sres2;
	sum(sres1,&R1[0],n);
	sum(sres2,&R2[0],n);
	if (sres1 > sres2) return true;
	if (sres1 < sres2) return false;
	return VecLessThan((const LAINT *)&R1[0], (const LAINT *)&R2[0], n);
}

void LieAlgRep::AddRoot(LAINT Degen, const LAINT *Root, const LAINT *Weight)
{
	Degens.push_back(Degen);
	Roots.AppendVector(Root);
	Weights.AppendVector(Weight);
}

void LieAlgRep::SortIndices(vector<size_t> &Indxs)
{	
	// Index sort by root sum; largest to smallest
	size_t n = Degens.size();
	Indxs.resize(n);
	for (size_t i=0; i<n; i++)
		Indxs[i] = i;
	CompareRootsBySum RootCompare;
	RootCompare.RtsPtr = &Roots;
	sort(Indxs.begin(), Indxs.end(), RootCompare); 
}

void LieAlgRep::Export(LieAlgRep &rep, vector<size_t> &Indxs)
{
	rep.set_vlen(vlen);
	for (vector<size_t>::iterator rit=Indxs.begin(); rit!= Indxs.end(); rit++)
	{
		size_t ix = *rit;
		LAINT Degen = Degens[ix];
		if (Degen == 0) continue; // Zero means absent - no need to export it
		MatrixRow<LAINT> Rt(Roots,ix);
		MatrixRow<LAINT> Wt(Weights,ix);
		rep.Degens.push_back(Degen);
		rep.Roots.AppendVector((const LAINT *)&Rt[0]);
		rep.Weights.AppendVector((const LAINT *)&Wt[0]);
	}
}

void LieAlgRep::Export(LieAlgRep &Rep)
{
	vector<size_t> Indxs;
	SortIndices(Indxs);
	Export(Rep,Indxs);
}

// Expanded rep object

// Adds the root (degen, rt vec, wt vec) if not already present,
// or its count (degen) if it is.
size_t LieAlgRepBuilder::AddRootOrCount(LAINT Degen, const LAINT *Root, const LAINT *Weight)
{
	pair<bool,size_t> ret = Indexer.AppendVector(Root);
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

static void RepRootVectors(LieAlgRepBuilder &bld, LieAlgebra &la, const LAINT *MaxWts)
{
	// Initial root
	vector<LAINT> MaxRts(la.rank);
	mul_vm(&MaxRts[0],MaxWts,la.InvCtnNum);
	bld.AddRootOrCount(0, &MaxRts[0], MaxWts);
	
	// Using LAINT instead of bool because in the STL,
	// bool gets turned into packed bits,
	// something that you can't use pointers or refs on
	Matrix<LAINT> WhichWay(1,la.rank);
	WhichWay.fill(true);
	
	vector<LAINT> NewRoot(la.rank), NewWeight(la.rank);
	vector<LAINT> NewWW(la.rank);
	
	// Find the next root until one cannot find any more
	size_t RtIndx = 0;
	while (RtIndx < bld.Roots.get_rows())
	{
		MatrixRow<LAINT> ThisRoot(bld.Roots,RtIndx);
		MatrixRow<LAINT> ThisWeight(bld.Weights,RtIndx);
		MatrixRow<LAINT> ThisWW(WhichWay,RtIndx);
		
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
						MatrixRow<LAINT> WWRow(WhichWay,RtIndx);
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


static void RepRootDegens(LieAlgRepBuilder &bld, LieAlgebra &la, LieAlgRep &rep)
{
	// Freudenthal's recurrence
	
	// Needs the max-to-min sorting order,
	// because the algorithm uses the highest weight.
	vector<size_t> Indxs;
	bld.SortIndices(Indxs);
	
	// Scale up the positive roots and their sum
	Matrix<LAINT> PosRoots = la.PosRoots;
	mulby_sm(PosRoots,la.InvCtnDen);
	vector<LAINT> PosRootSum = la.PosRootSum;
	mulby_sv(PosRootSum,la.InvCtnDen);
	
	// The main algorithm
	MatrixRow<LAINT> MaxRt(bld.Roots,Indxs[0]);	
	bld.Degens[Indxs[0]] = 1;
	vector<LAINT> BkRt(la.rank);
	vector<LAXINT> mrpd(la.rank), mrm(la.rank);
	for (vector<size_t>::iterator rit=Indxs.begin(); rit!= Indxs.end(); rit++)
	{
		if (rit == Indxs.begin()) continue;
		size_t ix = *rit;
		MatrixRow<LAINT> Rt(bld.Roots,ix);
		LAXINT nx = 0;
		for (size_t i=0; i<PosRoots.get_rows(); i++)
		{
			MatrixRow<LAINT> ShtRt(PosRoots,i);
			add_vv(BkRt,&Rt[0],&ShtRt[0]);
			while (true)
			{
				pair<bool,size_t> ret = bld.Indexer.VectorIndex(&BkRt[0]);
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

bool LieAlgebraRepDirect(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts)
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

static bool WeylOrbitForDomWt(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts)
{
	if (!la.IsValid) return false;
	
	LieAlgRepBuilder bld(la.rank);
	
	// Initial root
	vector<LAINT> MaxRts(la.rank);
	mul_vm(&MaxRts[0],MaxWts,la.InvCtnNum);
	bld.AddRootOrCount(0, &MaxRts[0], MaxWts);
	
	vector<LAINT> NewRoot(la.rank), NewWeight(la.rank);
	vector<LAINT> MetRoot(la.rank);
	
	// Find the next root until one cannot find any more
	size_t RtIndx = 0;
	while (RtIndx < bld.Roots.get_rows())
	{
		MatrixRow<LAINT> ThisRoot(bld.Roots,RtIndx);
		
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


static bool DomWtFromRoot(LAINT *MaxWts, LieAlgebra &la, const LAINT *Root)
{
	if (!la.IsValid) return false;
	
	// Set up the root
	vector<LAINT> ThisRoot(la.rank), ThisWeight(la.rank);
	vector<LAINT> NewRoot(la.rank), NewWeight(la.rank);
	vector<LAINT> MetRoot(la.rank);
	
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


template<class N> void MxWtToYD(N *YDPtr, N *MWPtr, size_t n)
{
	YDPtr[n-1] = MWPtr[n-1];
	for (int i=n-2; i>=0; i--)
		YDPtr[i] = YDPtr[i+1] + MWPtr[i];
}

template<class N> void Accumulate(N *AccSum, N *Orig, size_t n)
{
	AccSum[0] = Orig[0];
	for (int i=1; i<n; i++)
		AccSum[i] = AccSum[i-1] + Orig[i];
}


template<class N> void Permutations(Matrix<N> &Mat, const N *Syms, size_t n)
{
	Mat.resize(0,n);
	vector<N> SymVec(n);
	copy(Syms,Syms+n,SymVec.begin());
	sort(SymVec.begin(),SymVec.end());
	Mat.AppendVector(SymVec);
	while(next_permutation(SymVec.begin(),SymVec.end()))
		Mat.AppendVector(SymVec);
}

template<class N> void Permutations(Matrix<N> &Mat, vector<N> &Syms)
	{Permutations(Mat, &Syms[0], Syms.size());}

template<class N> void PermsForList(Matrix<N> &Mat, Matrix<N> &Syms)
{
	Mat.resize(0,Syms.get_cols());
	Matrix<N> SymPerms;
	for (int i=0; i<Syms.get_rows(); i++)
	{
		MatrixRow<N> SymRow(Syms,i);
		Permutations(SymPerms,&SymRow[0],Syms.get_cols());
		Mat.AppendMatrix(SymPerms);
	}
}

template<class N> void AddSigns(Matrix<N> &Mat, const N *Syms, size_t n)
{
	Mat.resize(0,n);
	vector<N> SymVec(n);
	vector<bool> SameAbsVal(n-1);
	
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

template<class N> void AddSigns(Matrix<N> &Mat, vector<N> &Syms)
	{AddSigns(Mat, &Syms[0], Syms.size());}

template<class N> void SelectParity(Matrix<N> &SelRts, Matrix<N> &OrigRts, N Val)
{
	if (Val == 0)
	{
		SelRts = OrigRts;
		return;
	}
	
	SelRts.resize(0,OrigRts.get_cols());
	
	for (int k=0; k<OrigRts.get_rows(); k++)
	{
		MatrixRow<N> Root(OrigRts,k);
		N SFac = Val;
		for (int i=0; i<Root.size(); i++)
			if (Root[i] < 0) SFac *= -1;
			else if (Root[i] == 0) SFac = 0;
		if (SFac > 0)
			SelRts.AppendVector(&Root[0]);
	}
}

static bool WeylOrbitForDomWtExplicit(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts)
{
	if (!la.IsValid) return false;
	if (la.family > 4) return false;
	
	rep.set_vlen(la.rank);
	vector<LAINT> Root(la.rank), Weight(la.rank);
	
	switch(la.family)
	{
	// A(n)
	case 1:
	{
		vector<LAINT> MWXtnd(la.rank+1);
		copy(MaxWts,MaxWts+la.rank,MWXtnd.begin());
		MWXtnd[la.rank] = 0;
		vector<LAINT> YDXtnd(la.rank+1);
		MxWtToYD(&YDXtnd[0],&MWXtnd[0],la.rank+1);
		LAINT YDSum = 0;
		for (int i=0; i<=la.rank; i++)
			YDSum += YDXtnd[i];
		for (int i=0; i<=la.rank; i++)
			YDXtnd[i] = (la.rank+1)*YDXtnd[i] - YDSum;
		Matrix<LAINT> RtXtnd;
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
		vector<LAINT> MWX(la.rank);
		copy(MaxWts,MaxWts+la.rank,MWX.begin());
		mulby_sv(&MWX[0],2,la.rank-1);
		vector<LAINT> YDX(la.rank);
		MxWtToYD(&YDX[0],&MWX[0],la.rank);
		Matrix<LAINT> SgndYD;
		AddSigns(SgndYD,&YDX[0],la.rank);
		Matrix<LAINT> Rts;
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
		vector<LAINT> MWX(la.rank);
		copy(MaxWts,MaxWts+la.rank,MWX.begin());
		vector<LAINT> YDX(la.rank);
		MxWtToYD(&YDX[0],&MWX[0],la.rank);
		Matrix<LAINT> SgndYD;
		AddSigns(SgndYD,&YDX[0],la.rank);
		Matrix<LAINT> Rts;
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
		vector<LAINT> MWX(la.rank);
		copy(MaxWts,MaxWts+la.rank,MWX.begin());
		mulby_sv(&MWX[0],2,la.rank-2);
		LAINT endm2 = MWX[la.rank-2];
		LAINT endm1 = MWX[la.rank-1];
		MWX[la.rank-2] = endm2 + endm1;
		vector<LAINT> YDX(la.rank);
		MxWtToYD(&YDX[0],&MWX[0],la.rank-1);
		YDX[la.rank-1] = - endm2 + endm1;
		Matrix<LAINT> SgndYD;
		AddSigns(SgndYD,&YDX[0],la.rank);
		Matrix<LAINT> PrtySgndYD;
		SelectParity(PrtySgndYD,SgndYD,YDX[la.rank-1]);
		Matrix<LAINT> Rts;
		PermsForList(Rts,PrtySgndYD);
		for (int k=0; k<Rts.get_rows(); k++)
		{
			MatrixRow<LAINT> Rt(Rts,k);
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


bool LieAlgebraWeylOrbitGeneral(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts)
{
	return WeylOrbitForDomWt(rep, la, MaxWts);
}


bool LieAlgebraWeylOrbitExplicit(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts)
{
	return WeylOrbitForDomWtExplicit(rep, la, MaxWts);
}


bool LieAlgebraWeylOrbit(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts)
{
	bool rc;
	rc = LieAlgebraWeylOrbitExplicit(rep, la, MaxWts);
	if (rc) return rc;
	rc = LieAlgebraWeylOrbitGeneral(rep, la, MaxWts);
	return rc;
}


// The main calculation of a (semi)simple-algebra irrep
// Only get dominant weights of Weyl orbits here

static void WeylDomWtRepRootVectors(LieAlgRepBuilder &bld, LieAlgebra &la, const LAINT *MaxWts)
{
	// Initial root
	vector<LAINT> MaxRts(la.rank);
	mul_vm(&MaxRts[0],MaxWts,la.InvCtnNum);
	bld.AddRootOrCount(0, &MaxRts[0], MaxWts);
	
	// Scale up the positive roots and their sum
	Matrix<LAINT> PosRoots = la.PosRoots;
	mulby_sm(PosRoots,la.InvCtnDen);
		
	vector<LAINT> NewRoot(la.rank), NewWeight(la.rank);
	
	// Find the next root until one cannot find any more
	size_t RtIndx = 0;
	while (RtIndx < bld.Roots.get_rows())
	{
		MatrixRow<LAINT> ThisRoot(bld.Roots,RtIndx);
		MatrixRow<LAINT> ThisWeight(bld.Weights,RtIndx);
		
		// Advance in each direction, if possible
		for (int i=0; i<PosRoots.get_rows(); i++)
		{
			MatrixRow<LAINT> ThisPosRoot(PosRoots, i);
			
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

static void WeylDomWtRepRootDegens(LieAlgRepBuilder &bld, LieAlgebra &la, LieAlgRep &rep)
{
	// Freudenthal's recurrence
	
	// Needs the max-to-min sorting order,
	// because the algorithm uses the highest weight.
	vector<size_t> Indxs;
	bld.SortIndices(Indxs);
	
	// Scale up the positive roots and their sum
	Matrix<LAINT> PosRoots = la.PosRoots;
	mulby_sm(PosRoots,la.InvCtnDen);
	vector<LAINT> PosRootSum = la.PosRootSum;
	mulby_sv(PosRootSum,la.InvCtnDen);
	
	// The main algorithm
	MatrixRow<LAINT> MaxRt(bld.Roots,Indxs[0]);	
	bld.Degens[Indxs[0]] = 1;
	vector<LAINT> BkRt(la.rank);
	vector<LAINT> DomRt(la.rank), DomWt(la.rank);
	vector<LAXINT> mrpd(la.rank), mrm(la.rank);
	for (vector<size_t>::iterator rit=Indxs.begin(); rit!= Indxs.end(); rit++)
	{
		if (rit == Indxs.begin()) continue;
		size_t ix = *rit;
		MatrixRow<LAINT> Rt(bld.Roots,ix);
		LAXINT nx = 0;
		for (size_t i=0; i<PosRoots.get_rows(); i++)
		{
			MatrixRow<LAINT> ShtRt(PosRoots,i);
			add_vv(BkRt,&Rt[0],&ShtRt[0]);
			// Dominant root and weight for it
			DomWtFromRoot(&DomWt[0], la, &BkRt[0]);
			mul_vm(DomRt, DomWt, la.InvCtnNum);
			while (true)
			{
				pair<bool,size_t> ret = bld.Indexer.VectorIndex(&DomRt[0]);
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

bool LieAlgebraRepWeylOrbits(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts)
{
	if (!la.IsValid) return false;
	for (size_t i=0; i<la.rank; i++)
		if (MaxWts[i] < 0) return false;
	LieAlgRepBuilder bld(la.rank);
	WeylDomWtRepRootVectors(bld,la,MaxWts);
	WeylDomWtRepRootDegens(bld,la,rep);
	return true;
}

bool LieAlgebraRepByWOs(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts)
{
	if (!la.IsValid) return false;
	
	LieAlgRepPtr DWPtr = GetRepObject(RO_REP_ORBIT, la, MaxWts);
	if (DWPtr.IsNull()) return false;
	
	LieAlgRep intrep(la.rank);
	
	size_t NDWS = DWPtr->Degens.size();
	for (int k=0; k<NDWS; k++)
	{
		LAINT Degen = DWPtr->Degens[k];
		MatrixRow<LAINT> DWT(DWPtr->Weights,k);
		
		LieAlgRepPtr WOPtr = GetRepObject(RO_ORBIT, la, &DWT[0]);
		if (WOPtr.IsNull()) continue;
		
		size_t NEnts = WOPtr->Degens.size();
		for (int i=0; i< NEnts; i++)
		{
			MatrixRow<LAINT> Root(WOPtr->Roots,i);
			MatrixRow<LAINT> Weight(WOPtr->Weights,i);
			intrep.AddRoot(Degen, &Root[0], &Weight[0]);
		}
	}
	
	intrep.Export(rep);
	return true;
}


bool LieAlgebraRepresentation(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts)
{
	return LieAlgebraRepByWOs(rep, la, MaxWts);
}


// For caching the representations

typedef bool (* InlineGetRep)(LieAlgRep &rep, LieAlgebra &la, const LAINT *MaxWts);

static InlineGetRep RepFuncs[NUMBER_OF_REP_OBJ_TYPES] =
{
	LieAlgebraWeylOrbit,
	LieAlgebraRepWeylOrbits,
	LieAlgebraRepresentation
};

struct LieAlgRepParams
{
	LieAlgebraParams AlgParams;
	vector<LAINT> MaxWts;
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
	return VecLessThan((const LAINT *)&R1.MaxWts[0], (const LAINT *)&R2.MaxWts[0], R1.AlgParams.rank);
}

static map<LieAlgRepParams, LieAlgRepPtr, LieAlgRepLessThan> LieAlgRepCaches[NUMBER_OF_REP_OBJ_TYPES];

LieAlgRepPtr InvalidLieAlgRep(NULL);

LieAlgRepPtr &GetRepObject(enum RepObjType rotype, const LieAlgebraParams &AlgParams, const LAINT *MaxWts)
{
	map<LieAlgRepParams, LieAlgRepPtr, LieAlgRepLessThan> &LieAlgRepCache = LieAlgRepCaches[rotype];
	InlineGetRep RepFunc = RepFuncs[rotype];
	LieAlgRepParams RepParams;
	RepParams.AlgParams = AlgParams;
	RepParams.MaxWts.resize(RepParams.AlgParams.rank);
	copy(MaxWts, MaxWts+AlgParams.rank, RepParams.MaxWts.begin());
	map<LieAlgRepParams, LieAlgRepPtr, LieAlgRepLessThan>::iterator it =
		LieAlgRepCache.find(RepParams);
	if (it == LieAlgRepCache.end())
	{
		pair<LieAlgRepParams, LieAlgRepPtr> rec;
		LieAlgRep &Rep = *rec.second;
		if (!RepFunc(Rep,GetLieAlgebra(AlgParams),MaxWts))
			return InvalidLieAlgRep;
		rec.first = RepParams;
		pair<map<LieAlgRepParams, LieAlgRepPtr, LieAlgRepLessThan>::iterator, bool> ret =
			LieAlgRepCache.insert(rec);
		return ret.first->second;
	}
	else
		return it->second;
}
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, const LieAlgebraParams &AlgParams, vector<LAINT> &MaxWts)
	{return GetRepObject(rotype, AlgParams, (const LAINT *)&MaxWts[0]);}

LieAlgRepPtr &GetRepObject(enum RepObjType rotype, LAINT family, LAINT rank, const LAINT *MaxWts)
{
	LieAlgebraParams AlgParams;
	AlgParams.family = family; AlgParams.rank = rank;
	return GetRepObject(rotype, AlgParams, MaxWts);
}
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, LAINT family, LAINT rank, vector<LAINT> &MaxWts)
{
	LieAlgebraParams AlgParams;
	AlgParams.family = family; AlgParams.rank = rank;
	return GetRepObject(rotype, AlgParams, (const LAINT *)&MaxWts[0]);
}

LieAlgRepPtr &GetRepObject(enum RepObjType rotype, LieAlgebra &la, const LAINT *MaxWts)
	{return GetRepObject(rotype, la.GetParams(),MaxWts);}
LieAlgRepPtr &GetRepObject(enum RepObjType rotype, LieAlgebra &la, LAINT rank, vector<LAINT> &MaxWts)
	{return GetRepObject(rotype, la.GetParams(),(const LAINT *)&MaxWts[0]);}

void ClearLieAlgReps()
{
	for (int rotype=0; rotype<NUMBER_OF_REP_OBJ_TYPES; rotype++)
	{
		map<LieAlgRepParams, LieAlgRepPtr, LieAlgRepLessThan> &LieAlgRepCache = LieAlgRepCaches[rotype];
		LieAlgRepCache.clear();
	}
}


// Find rep from product algebra -- (semi)simple ones and U(1)'s
LAINT LieAlgProduct::get_rank()
{
	LAINT rank = 0;
	for (vector<LieAlgebraParams>::iterator it = ParamList.begin();
		it != ParamList.end(); it++)
		rank += it->rank;
	rank += NumU1s;
	return rank;
}

void LieAlgProdRepresentation(LieAlgRep &rep, LieAlgProduct &lap,
	enum RepObjType rotype, const LAINT *MaxWts)
{
	LAINT rank = lap.get_rank();
	rep.set_vlen(rank);
	
	// Create the initial root -- a singlet
	vector<LAINT> Rt(rank), Wt(rank);
	fill(Rt.begin(), Rt.end(), 0);
	fill(Wt.begin(), Wt.end(), 0);
	size_t nustrt = rank - lap.NumU1s;
	copy(MaxWts+nustrt, MaxWts+rank, &Rt[nustrt]);
	copy(MaxWts+nustrt, MaxWts+rank, &Wt[nustrt]);
	rep.AddRoot(1,Rt,Wt);
	
	LieAlgRep NewRep; // Temporary for each run
	size_t stix = 0, ndix; // Start index
	
	for (vector<LieAlgebraParams>::iterator it = lap.ParamList.begin();
		it != lap.ParamList.end(); it++)
	{
		LieAlgebraParams &Params = *it; 
		ndix = stix + Params.rank;
		
		LieAlgRepPtr &RepAddedPtr = GetRepObject(rotype, Params,MaxWts+stix);
		LieAlgRep &RepAdded = *RepAddedPtr;
		NewRep.set_vlen(rank);
		for (int i=0; i<rep.Degens.size(); i++)
		{
			// Get the entire root/weight from the existing array
			LAINT Degen = rep.Degens[i];
			MatrixRow<LAINT> Rtx(rep.Roots,i);
			MatrixRow<LAINT> Wtx(rep.Weights,i);
			copy(Rtx.begin(), Rtx.end(), Rt.begin());
			copy(Wtx.begin(), Wtx.end(), Wt.begin());
			for (int j=0; j<RepAdded.Degens.size(); j++)
			{
				// Copy in and add the new root
				LAINT NewDegen = Degen*RepAdded.Degens[j];
				MatrixRow<LAINT> Rta(RepAdded.Roots,j);
				MatrixRow<LAINT> Wta(RepAdded.Weights,j);
				copy(Rta.begin(), Rta.end(), &Rt[stix]);
				copy(Wta.begin(), Wta.end(), &Wt[stix]);
				NewRep.AddRoot(NewDegen,Rt,Wt);
			}
		}
		swap(rep,NewRep); // Put new ones into old ones' location
		stix = ndix;
	}
}


// Abstract base class for representation handler
// Subclass for single-algebra and product-of-algebras cases

// Counted list of weight vectors to rep
LieAlgRepPtr RepHandlerBase::GetRepPtr(enum RepObjType rotype, LACntdMaxWtList &CWL)
{
	LieAlgRepBuilder Bld(get_rank());
	
	for (LACntdMaxWtList::iterator eniter = CWL.begin();
		eniter != CWL.end(); eniter++)
	{
		LACntdMaxWtEntry &Entry = *eniter;
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
LACntdMaxWtList RepHandlerBase::ExtractWts(enum RepObjType rotype, LieAlgRepBuilder &Bld)
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
		MatrixRow<LAINT> MaxWtRow(Bld.Weights,MaxIndx);
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
