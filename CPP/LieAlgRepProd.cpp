/*
	Implements the rep-product and rep-power functions
*/

#include "LieAlgRepProd.h"


// For assembling irreps from Weyl orbits
// instead of from their entire contents
static bool IsDomWt(LAINT *Wts, size_t n)
{
	if (n == 0) return true;
	LAINT minwt;
	index_min_v(minwt, Wts, n);
	return (minwt >= 0);
}


void DoRepProduct(LieAlgRepBuilder &Bld, size_t TotRank,
	LieAlgRep &Rep1, LieAlgRep &Rep2)
{
	LAINT rank = Rep1.Roots.get_cols();
	
	Bld.set_vlen(rank);
	vector<LAINT> Roots(rank), Weights(rank);
	
	for (size_t i1=0; i1<Rep1.Degens.size(); i1++)
		for (size_t i2=0; i2<Rep2.Degens.size(); i2++)
		{
			LAINT Degen = Rep1.Degens[i1]*Rep2.Degens[i2];
			add_vv(Roots, &Rep1.Roots(i1,0), &Rep2.Roots(i2,0));
			add_vv(Weights, &Rep1.Weights(i1,0), &Rep2.Weights(i2,0));
			if (IsDomWt(&Weights[0],TotRank))
				Bld.AddRootOrCount(Degen,Roots,Weights);
		}
}


// For iterating over all the possible index combinations.
// With the condition that each one >= the one previous.
struct IndexSetIterator
{
	size_t MaxIndex;
	vector<size_t> Indices;
	map<size_t,LAINT> IndexCount;
	
	// Creates and sets it up for a run
	// Could be broken up if necessary
	IndexSetIterator(int Power): Indices(Power)
		{fill(Indices.begin(), Indices.end(), 0); Count();}
	
	// Go the the next index; return whether doing so was successful
	// Failing to do so will end the run
	bool Next();
	
	// Count up the indices
	void Count();
};

bool IndexSetIterator::Next()
{
	int NumIxs = Indices.size();
	int i;
	for (i = NumIxs-1; i>=0; i--)
	{
		if ((++Indices[i]) < MaxIndex) break;
	}
	if (i < 0) return false;
	size_t NextIndex = Indices[i];
	for (i=i+1; i<NumIxs; i++)
		Indices[i] = NextIndex;
	
	Count();
	return true;
}

void IndexSetIterator::Count()
{
	IndexCount.clear();
	for (vector<size_t>::iterator IItr = Indices.begin();
		IItr != Indices.end(); IItr++)
	{
		size_t &ix = *IItr;
		map<size_t,LAINT>::iterator CItr = IndexCount.find(ix);
		if (CItr == IndexCount.end())
		{
			IndexCount[ix] = 1;
		}
		else
		{
			CItr->second += 1;
		}
	}
}

void DoRepPwrSym(LieAlgRepBuilder &Bld, size_t TotRank,
	LieAlgRep &Rep, LAINT Power, LAINT Symm)
{
	LAINT rank = Rep.Roots.get_cols();
	vector<LAINT> Roots(rank), Weights(rank);
	if (Power < 1)
	{
		// Zero: scalar rep
		fill(Roots.begin(),Roots.end(),0);
		fill(Weights.begin(),Weights.end(),0);
		Bld.set_vlen(rank);
		if (IsDomWt(&Weights[0],TotRank))
			Bld.AddRootOrCount(1,Roots,Weights);
		return;
	}
	
	Bld.set_vlen(rank);
	LAXINT Symx = (Symm >= 0) ? 1 : -1;
	
	IndexSetIterator Iter(Power);
	Iter.MaxIndex = Rep.Degens.size();
	do
	{
		LAXINT Degen = 1;
		for (map<size_t,LAINT>::iterator CItr = Iter.IndexCount.begin();
			CItr != Iter.IndexCount.end(); CItr++)
		{
			size_t ix = CItr->first;
			LAINT dgn = Rep.Degens[ix];
			LAINT mult = CItr->second;
			LAXINT IDgn = 1;
			for (LAINT i=0; i<mult; i++)
				IDgn = ((dgn+Symx*i)*IDgn)/(i+1);
			Degen *= IDgn;
		}
		
		// Assemble the root and weight vectors
		fill(Roots.begin(),Roots.end(),0);
		fill(Weights.begin(),Weights.end(),0);
		for (map<size_t,LAINT>::iterator CItr = Iter.IndexCount.begin();
			CItr != Iter.IndexCount.end(); CItr++)
		{
			size_t ix = CItr->first;
			LAINT mult = CItr->second;
			muladdto_vsv(Roots,mult,&Rep.Roots(ix,0));
			muladdto_vsv(Weights,mult,&Rep.Weights(ix,0));
		}
		
		// Add it into the target rep
		if (IsDomWt(&Weights[0],TotRank))
			Bld.AddRootOrCount(Degen,Roots,Weights);
	}
	while(Iter.Next());
}

// Need to construct tensor powers of vectors in SL(n) ~ SU(n) for the plethysms
// the rep-power decompositions by symmetry type

struct YDE_Extended: public YoungDiagramEntry
{
	// Extra stuff for convenience in calculating
	LAINT YDLen;
	LAINT MultNomFac;
};

// Nesting matrix:
// Row indexing: by starting YD
// Column indexing: by concatenation of nested YD's

struct NestingEntry
{
	vector<LAINT> SubYD_Indices;
	vector<LAINT> SubYD_Lens;
	LAINT MultNomFac;
};

typedef vector<NestingEntry> NestingList;

struct YoungDiagramTensorPower
{
	vector<YDE_Extended> YDList;
	Matrix<LAINT> Kostka;
	Matrix<NestingList> Nesting;
	
	// Starting from the YD's, calculate everything else
	// YDLen and MultNomFac in the YDList elements
	// Kostka and Nesting in overall object
	void FillOut();
};

typedef vector<YoungDiagramTensorPower> YoungDiagramTensorPowerList;

static YoungDiagramTensorPowerList TPList;

static YoungDiagramTensorPower InvalidTPEntry;

struct YDIndex
{
	LAINT Power, Index;
};

static map<YoungDiagram,YDIndex> YDIndices;

static void YDTensAddVector(vector<YoungDiagram> &NewYDs, YoungDiagram &YD)
{
	NewYDs.clear();
	YoungDiagram NewYD(YD.size());
	
	// Top one
	copy(YD.begin(),YD.end(),NewYD.begin());
	NewYD[0] += 1;
	NewYDs.push_back(NewYD);
	
	// Middle ones:
	for (size_t i=1; i<YD.size(); i++)
	{
		if (YD[i] < YD[i-1])
		{
			copy(YD.begin(),YD.end(),NewYD.begin());
			NewYD[i] += 1;
			NewYDs.push_back(NewYD);
		}
	}
	
	// Bottom one:
	copy(YD.begin(),YD.end(),NewYD.begin());
	NewYD.push_back(1);
	NewYDs.push_back(NewYD);
}

void YDTensSetAddVector(vector<YDE_Extended> &NewYDList, vector<YDE_Extended> &YDList)
{
	vector<YoungDiagram> NewYDs;
	NewYDList.clear();
	for (vector<YDE_Extended>::iterator YDEIter = YDList.begin();
		YDEIter != YDList.end(); YDEIter++)
	{
		YDE_Extended &YDEX = *YDEIter;
		YDTensAddVector(NewYDs,YDEX.YD);
		for (vector<YoungDiagram>::iterator NewYDIter = NewYDs.begin();
			NewYDIter != NewYDs.end(); NewYDIter++)
		{
			YoungDiagram &YD = *NewYDIter;
			map<YoungDiagram,YDIndex>::iterator YDIFind = YDIndices.find(YD);
			if (YDIFind == YDIndices.end())
			{
				// Add the diagram to the indexer
				YDIndex YIX;
				sum(YIX.Power,YD);
				YIX.Index = NewYDList.size();
				YDIndices[YD] = YIX;
				YDE_Extended NewYDEX;
				NewYDEX.YD = YD;
				NewYDEX.Mult = YDEX.Mult;
				NewYDList.push_back(NewYDEX);
			}
			else
			{
				// Bump up the multiplicity appropriately
				YDIndex &YIX = YDIFind->second;
				NewYDList[YIX.Index].Mult += YDEX.Mult;
			}
		}
	}
}

typedef long long FACTORIAL_TYPE;

static FACTORIAL_TYPE Factorial(LAINT n)
{
	long long fac = 1;
	for (LAINT i=1; i<=n; i++)
		fac *= i;
	return fac;
}

static LAINT MultNomFactorial(YoungDiagram &YD)
{
	map<LAINT,LAINT> LenCounts;
	LAINT TotLen = 0;
	for (YoungDiagram::iterator Iter = YD.begin(); Iter != YD.end(); Iter++)
	{
		LAINT Row = *Iter;
		map<LAINT,LAINT>::iterator LCFind = LenCounts.find(Row);
		if (LCFind == LenCounts.end())
		{
			// Start
			LenCounts[Row] = 1;
		}
		else
		{
			// Bump up
			LCFind->second += 1;
		}
		TotLen += 1;
	}
	FACTORIAL_TYPE MNFNum = Factorial(TotLen);
	FACTORIAL_TYPE MNFDen = 1;
	for (map<LAINT,LAINT>::iterator LCIter = LenCounts.begin();
		LCIter != LenCounts.end(); LCIter++)
	{
		LAINT NumRow = LCIter->second;
		MNFDen *= Factorial(NumRow);
	}
	return LAINT(MNFNum/MNFDen);
}

static void KostkaMatrix(Matrix<LAINT> &Kostka, vector<YDE_Extended> &YDList)
{
	size_t NumYDs = YDList.size();
	Matrix<LAXINT> LookBack(NumYDs);
	LookBack.fill(0);
	vector<LAXINT> DenFac(NumYDs);
	YoungDiagram LBYD;
	
	for (size_t i=0; i<NumYDs; i++)
	{
		YoungDiagram &YD = YDList[i].YD;
		LAINT YDLen = YD.size();
		// Do lookback; the += takes care of the multiplicities automatically.
		for (LAINT j=0; j<(YDLen-1); j++)
			for (LAINT jx=j+1; jx<YDLen; jx++)
				for (LAINT k=0; k<=YD[jx]; k++)
				{
					LBYD = YD;
					LBYD[j] += k;
					LBYD[jx] -= k;
					LAXINT dgd = LBYD[j] - LBYD[jx];
					sort(LBYD.begin(),LBYD.end());
					reverse(LBYD.begin(),LBYD.end());
					if (LBYD.back() == 0) LBYD.pop_back();
					YDIndex &YDIX = YDIndices[LBYD];
					LookBack(i,YDIX.Index) += 2*dgd;
				}
		// Denominator factor
		LAXINT DFSum = 0;
		for (LAINT j=0; j<YDLen; j++)
		{
			LAINT &d = YD[j];
			DFSum += d*(d - 2*(j+1));
		}
		DenFac[i] = DFSum;
	}
	
	Kostka.resize(NumYDs);
	Kostka.fill(0);
	for (size_t i=0; i<NumYDs; i++)
	{
		MatrixRow<LAINT> KKRow(Kostka,i);
		KKRow[i] = 1;
		for (size_t j = i+1; j<NumYDs; j++)
		{
			LAXINT KKNum = 0;
			MatrixRow<LAXINT> LBRow(LookBack,j);
			for (size_t k=i; k<j; k++)
				KKNum += LBRow[k]*KKRow[k];
			if (KKNum != 0)
				KKRow[j] = KKNum/(DenFac[i] - DenFac[j]);
		}
	}
}

// Get whichever list of YD's will be necessary
// Assumes that all the lists before the recently-calculated one
// are present in TPList
static vector<YDE_Extended> &NestingGetYDList(LAINT Power,
	vector<YDE_Extended> &YDList)
{
	if (Power > TPList.size()) return YDList;
	else return TPList[Power-1].YDList;
}

static void NestingMatrix(Matrix< vector<NestingEntry> > &Nesting,
	vector<YDE_Extended> &YDList)
{
	size_t NumYDs = YDList.size();
	Nesting.resize(NumYDs);
	vector<LAINT> SubdgrmIndices, NewSDI;
	YoungDiagram ConcatSubdgrms;
	for (size_t i=0; i<NumYDs; i++)
	{
		YoungDiagram &YD = YDList[i].YD;
		LAINT YDLen = YD.size();
		// Assemble a list of all possible indices of subdiagrams.
		// The subdiagrams have size (row length) in YD.
		// They are stored back-to-back in the list
		SubdgrmIndices.clear();
		for (LAINT ir=0; ir<YDLen; ir++)
		{
			size_t NumIndices = NestingGetYDList(YD[ir],YDList).size();
			NewSDI.clear();
			size_t NumOldEntries = (ir > 0) ? SubdgrmIndices.size()/ir : 1;
			// Existing-list loop
			for (size_t j0=0; j0<NumOldEntries; j0++)
			{
				// New indices
				for (LAINT j1=0; j1<NumIndices; j1++)
				{
					// Copy over then append
					for (LAINT k=0; k<ir; k++)
						NewSDI.push_back(SubdgrmIndices[ir*j0+k]);
					NewSDI.push_back(j1);
				}
			}
			swap(SubdgrmIndices,NewSDI); // Put new ones into old ones' location
		}
		
		size_t NumIndexSets = SubdgrmIndices.size()/YDLen;
		for (size_t ix=0; ix<NumIndexSets; ix++)
		{
			NestingEntry NE;
			for (LAINT ir=0; ir<YDLen; ir++)
				NE.SubYD_Indices.push_back(SubdgrmIndices[YDLen*ix+ir]);
			NE.MultNomFac = 1;
			
			ConcatSubdgrms.clear();
			for (LAINT ir=0; ir<YDLen; ir++)
			{
				vector<YDE_Extended> &SubYDList = NestingGetYDList(YD[ir],YDList);
				YDE_Extended &SubYDEntry = SubYDList[NE.SubYD_Indices[ir]];
				for (YoungDiagram::iterator YDIter = SubYDEntry.YD.begin();
					YDIter != SubYDEntry.YD.end(); YDIter++)
					ConcatSubdgrms.push_back(*YDIter);
				NE.SubYD_Lens.push_back(SubYDEntry.YDLen);
				NE.MultNomFac *= SubYDEntry.MultNomFac;
			}
			
			sort(ConcatSubdgrms.begin(),ConcatSubdgrms.end());
			reverse(ConcatSubdgrms.begin(),ConcatSubdgrms.end());
			YDIndex &YDIX = YDIndices[ConcatSubdgrms];
			Nesting(i,YDIX.Index).push_back(NE);
		}
	}
}


void YoungDiagramTensorPower::FillOut()
{
	for (vector<YDE_Extended>::iterator YDIter = YDList.begin();
		YDIter != YDList.end(); YDIter++)
	{
		YoungDiagram &YD = YDIter->YD;
		YDIter->YDLen = YD.size();
		YDIter->MultNomFac = MultNomFactorial(YD);
	}
	
	KostkaMatrix(Kostka,YDList);
	NestingMatrix(Nesting,YDList);
};

static YoungDiagramTensorPower FirstTensorPower()
{
	// For a vector: power = 1
	YoungDiagramTensorPower TPEntry;
	YDE_Extended YDEntry;
	YDEntry.YD.push_back(1);
	YDEntry.Mult = 1;
	TPEntry.YDList.push_back(YDEntry);
	TPEntry.FillOut();
	
	return TPEntry;
}

// Next entry
static YoungDiagramTensorPower NextTensorPower(YoungDiagramTensorPower &TPEntry)
{
	// Multiply tensor by vector -- easiest
	YoungDiagramTensorPower NextTPEntry;
	YDTensSetAddVector(NextTPEntry.YDList,TPEntry.YDList);
	NextTPEntry.FillOut();
	
	return NextTPEntry;
}

// Which entry
static YoungDiagramTensorPower &GetTensorPower(LAINT p)
{
	if (p <= 0) return InvalidTPEntry;
	
	if (TPList.size() < 1)
		TPList.push_back(FirstTensorPower());
	
	size_t tpsz = TPList.size();
	for (LAINT pi = tpsz; pi < p; pi++)
		TPList.push_back(NextTensorPower(TPList.back()));
	
	return TPList[p-1];
}


struct IndexMult
{
	size_t Index;
	LAINT Mult;
};

bool IndexMultLessThan(const IndexMult &IM1, const IndexMult &IM2)
{
	// First sort by reverse order in multiplicities
	if (IM1.Mult > IM2.Mult) return true;
	if (IM1.Mult < IM2.Mult) return false;
	// Sort in forward order by indices
	return (IM1.Index < IM2.Index);
}


void DoRepPower(LAYDRepBldList &BldList, size_t TotRank,
	LieAlgRep &Rep, LAINT Power)
{
	LAINT rank = Rep.Roots.get_cols();
	vector<LAINT> Roots(rank), Weights(rank);
	if (Power < 1)
	{
		// Zero: scalar rep
		fill(Roots.begin(),Roots.end(),0);
		fill(Weights.begin(),Weights.end(),0);
		BldList.resize(1);
		BldList[0].Bld.set_vlen(rank);
		if (IsDomWt(&Weights[0],TotRank))
			BldList[0].Bld.AddRootOrCount(1,Roots,Weights);
		return;
	}
	
	YoungDiagramTensorPower &TPEntry = GetTensorPower(Power);
	vector<YDE_Extended> &YDList = TPEntry.YDList;
	LAINT NumYDs = YDList.size();
	BldList.resize(NumYDs);
	for (LAINT i=0; i<NumYDs; i++)
	{
		YDE_Extended &YDEntry = YDList[i];
		LAYDRepBldEntry &BldEntry = BldList[i];
		BldEntry.Bld.set_vlen(rank);
		BldEntry.YDE = YDEntry;
	}
	
	// Index multiplicities -> Young diagrams
	// Also get degens for indices in the YD order
	vector<IndexMult> IndexMultList;
	YoungDiagram MultYD;
	vector<LAINT> RepIxDegens;
	// Intermediates: orbit and rep multiplicities
	// Orbit directly from nesting, rep from orbit and Kostka matrix
	vector<LAXINT> OrbDegens(NumYDs), RepPwrDegens(NumYDs);
	
	IndexSetIterator Iter(Power);
	Iter.MaxIndex = Rep.Degens.size();
	do
	{
		// Find the YD from the index multiplicities
		// Be sure that the indices are sorted along with their multiplicities
		// when getting those multiplicities into YD order
		// Get the rep degens with those indices
		IndexMultList.clear();
		for (map<size_t,LAINT>::iterator CItr = Iter.IndexCount.begin();
			CItr != Iter.IndexCount.end(); CItr++)
		{
			size_t ix = CItr->first;
			LAINT mult = CItr->second;
			IndexMult IM;
			IM.Index = ix;
			IM.Mult = mult;
			IndexMultList.push_back(IM);
		}
		sort(IndexMultList.begin(),IndexMultList.end(),IndexMultLessThan);
		MultYD.clear();
		RepIxDegens.clear();
		for (vector<IndexMult>::iterator IMLIter = IndexMultList.begin();
			IMLIter != IndexMultList.end(); IMLIter++)
		{
			MultYD.push_back(IMLIter->Mult);
			RepIxDegens.push_back(Rep.Degens[IMLIter->Index]);
		}
		YDIndex &YDIX = YDIndices[MultYD];
		LAINT MultIndex = YDIX.Index;
		MatrixRow<NestingList> NestingRow(TPEntry.Nesting,MultIndex);
		
		// Find the orb mults
		for (LAINT i=0; i<NumYDs; i++)
		{
			NestingList &NList = NestingRow[i];
			LAXINT DegenSum = 0;
			for (NestingList::iterator NLIter = NList.begin();
				NLIter != NList.end(); NLIter++)
			{
				NestingEntry &NE = *NLIter;
				LAINT DegenIndSum = NE.MultNomFac;
				for (LAINT j=0; j<NE.SubYD_Lens.size(); j++)
				{
					LAINT SYLen = NE.SubYD_Lens[j];
					LAINT Degen = RepIxDegens[j];
					for (LAINT k=0; k<SYLen; k++)
						DegenIndSum = ((Degen-k)*DegenIndSum)/(k+1);
				}
				DegenSum += DegenIndSum;
			}
			OrbDegens[i] = DegenSum;
		}
		
		// Find the rep mults
		mul_mv(RepPwrDegens,TPEntry.Kostka,OrbDegens);
		
		// Assemble the root and weight vectors
		fill(Roots.begin(),Roots.end(),0);
		fill(Weights.begin(),Weights.end(),0);
		for (map<size_t,LAINT>::iterator CItr = Iter.IndexCount.begin();
			CItr != Iter.IndexCount.end(); CItr++)
		{
			size_t ix = CItr->first;
			LAINT mult = CItr->second;
			muladdto_vsv(Roots,mult,&Rep.Roots(ix,0));
			muladdto_vsv(Weights,mult,&Rep.Weights(ix,0));
		}
		
		// Add them
		for (LAINT i=0; i<NumYDs; i++)
			if (IsDomWt(&Weights[0],TotRank))
				BldList[i].Bld.AddRootOrCount(RepPwrDegens[i],Roots,Weights);
	}
	while(Iter.Next());
}
