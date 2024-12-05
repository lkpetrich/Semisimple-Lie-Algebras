/*
	Implements the Lie-algebra class
*/

#include <map>
#include "LieAlgebra.h"


// Returns whether the family and rank are valid
bool LieAlgebra::MakeDynkin()
{
	// Declarations
	int i;
	RootConnectionType rtconn;
	
	// n = rank
	switch(family)
	{
	case 1:
	{
		// A(n) or SU(n+1); n >= 1
		if (rank < 1) return false;
		for (i=1; i<=rank; i++)
			dynkin.rtlens.push_back(1);
		for (i=1; i<=rank-1; i++)
		{
			rtconn.root1 = i;
			rtconn.root2 = i+1;
			rtconn.strength = 1;
			dynkin.rtconns.push_back(rtconn);
		}
		return true;
	}
	case 2:
	{
		// B(n) or SO(2n+1); n >= 1
		if (rank < 1) return false;
		else if (rank == 1)
		{
			dynkin.rtlens.push_back(1);
			return true;
		}
		for (i=1; i<=rank-1; i++)
			dynkin.rtlens.push_back(2);
		dynkin.rtlens.push_back(1);
		for (i=1; i<=rank-2; i++)
		{
			rtconn.root1 = i;
			rtconn.root2 = i+1;
			rtconn.strength = 2;
			dynkin.rtconns.push_back(rtconn);
		}
		rtconn.root1 = rank-1;
		rtconn.root2 = rank;
		rtconn.strength = 2;
		dynkin.rtconns.push_back(rtconn);
		return true;
	}
	case 3:
	{
		// C(n) or Sp(2n); n >= 1
		if (rank < 1) return false;
		else if (rank == 1)
		{
			dynkin.rtlens.push_back(1);
			return true;
		}
		for (i=1; i<=rank-1; i++)
			dynkin.rtlens.push_back(1);
		dynkin.rtlens.push_back(2);
		for (i=1; i<=rank-2; i++)
		{
			rtconn.root1 = i;
			rtconn.root2 = i+1;
			rtconn.strength = 1;
			dynkin.rtconns.push_back(rtconn);
		}
		rtconn.root1 = rank-1;
		rtconn.root2 = rank;
		rtconn.strength = 2;
		dynkin.rtconns.push_back(rtconn);
		return true;
	}
	case 4:
	{
		// D(n) or SO(2n); n >= 2
		if (rank < 2) return false;
		else if (rank == 2)
		{
			dynkin.rtlens.push_back(1);
			dynkin.rtlens.push_back(1);
			return true;
		}
		for (i=1; i<=rank; i++)
			dynkin.rtlens.push_back(1);
		for (i=1; i<=rank-2; i++)
		{
			rtconn.root1 = i;
			rtconn.root2 = i+1;
			rtconn.strength = 1;
			dynkin.rtconns.push_back(rtconn);
		}
		rtconn.root1 = rank-2;
		rtconn.root2 = rank;
		rtconn.strength = 1;
		dynkin.rtconns.push_back(rtconn);
		return true;
	}
	case 5:
	{
		// E(n); n = 6,7,8
		if (rank < 6 || rank > 8) return false;
		for (i=1; i<=rank; i++)
			dynkin.rtlens.push_back(1);
		for (i=1; i<=rank-2; i++)
		{
			rtconn.root1 = i;
			rtconn.root2 = i+1;
			rtconn.strength = 1;
			dynkin.rtconns.push_back(rtconn);
		}
		rtconn.root1 = 3;
		rtconn.root2 = rank;
		rtconn.strength = 1;
		dynkin.rtconns.push_back(rtconn);
		return true;
	}
	case 6:
	{
		// F(n); n = 4
		if (rank != 4) return false;
		dynkin.rtlens.push_back(2);
		dynkin.rtlens.push_back(2);
		dynkin.rtlens.push_back(1);
		dynkin.rtlens.push_back(1);
		rtconn.root1 = 1;
		rtconn.root2 = 2;
		rtconn.strength = 2;
		dynkin.rtconns.push_back(rtconn);
		rtconn.root1 = 2;
		rtconn.root2 = 3;
		rtconn.strength = 2;
		dynkin.rtconns.push_back(rtconn);
		rtconn.root1 = 3;
		rtconn.root2 = 4;
		rtconn.strength = 1;
		dynkin.rtconns.push_back(rtconn);
		return true;
	}
	case 7:
	{
		// G(n); n = 2
		if (rank != 2) return false;
		dynkin.rtlens.push_back(3);
		dynkin.rtlens.push_back(1);
		rtconn.root1 = 1;
		rtconn.root2 = 2;
		rtconn.strength = 3;
		dynkin.rtconns.push_back(rtconn);
		return true;
	}
	default:
		return false;
	}
}

// Set up the metric and the Cartan matrix

// Returns whether the matrix is nonsingular
template<typename N> bool SetupMatrices(const Matrix<N> &Mat,
	Matrix< Fraction<N> > &InvMat,
	Matrix<N> &InvMatNum,
	N &InvMatDen)
{
	// Singular?
	if (!MatrixInverse(InvMat,Mat)) return false;
	
	// Find the least common multiple of the denominators
	InvMatDen = N(1);
	for (int i=0; i<InvMat.get_rows(); i++)
		for (int j=0; j<InvMat.get_cols(); j++)
		{
			N den = InvMat(i,j).get_den();
			InvMatDen = LCM(InvMatDen,den);
		}
	
	// Multiply the numerators by the LCM, then divide by the denominators
	InvMatNum.resize(InvMat.get_rows(),InvMat.get_cols());
	for (int i=0; i<InvMat.get_rows(); i++)
		for (int j=0; j<InvMat.get_cols(); j++)
		{
			Fraction<N> &val = InvMat(i,j);
			InvMatNum(i,j) = val.get_num()*(InvMatDen/val.get_den());
		}
	
	return true;
}

bool LieAlgebra::MakeMetric()
{
	int nc = 0;
	Metric.resize(rank);
	Metric.fill(0);
	
	for (int k=0; k<dynkin.rtlens.size(); k++)
		Metric(k,k) = LAINT(2)*dynkin.rtlens[k];
	
	for (int k=0; k<dynkin.rtconns.size(); k++)
	{
		RootConnectionType &rtconn = dynkin.rtconns[k];
		LAINT i = rtconn.root1 - LAINT(1);
		LAINT j = rtconn.root2 - LAINT(1);
		LAINT &str = rtconn.strength;
		Metric(i,j) = - str;
		Metric(j,i) = - str;
	}
	
	return SetupMatrices(Metric, InverseMetric, InvMetNum, InvMetDen);
}

bool LieAlgebra::MakeCartan()
{
	Cartan.resize(rank);
	Cartan.fill(0);
	
	for (int k=0; k<dynkin.rtlens.size(); k++)
		Cartan(k,k) = LAINT(2);
	
	for (int k=0; k<dynkin.rtconns.size(); k++)
	{
		RootConnectionType &rtconn = dynkin.rtconns[k];
		LAINT i = rtconn.root1 - LAINT(1);
		LAINT j = rtconn.root2 - LAINT(1);
		LAINT &str = rtconn.strength;
		Cartan(i,j) = - str/dynkin.rtlens[j];
		Cartan(j,i) = - str/dynkin.rtlens[i];
	}
	
	return SetupMatrices(Cartan, InverseCartan, InvCtnNum, InvCtnDen);
}

// Positive roots
void LieAlgebra::MakePosRoots()
{
	VectorSet<LAINT,LARootHashFunction> Indexer(rank);
	
	PosRoots.resize(rank,rank);
	PosRoots.fill(LAINT(0));
	for (int i=0; i<rank; i++)
		PosRoots(i,i) = LAINT(1);
	
	PosWeights = Cartan;
	
	// Using LAINT instead of bool because in the STL,
	// bool gets turned into packed bits,
	// something that you can't use pointers or refs on
	LAINT_MATRIX WhichWay(rank,rank);
	WhichWay.fill(true);
	
	LAINT_VECTOR NewRoot(rank), NewWeight(rank), NwRtRed(rank);
	LAINT_VECTOR NewWW(rank);
	
	for (int i=0; i<rank; i++)
		Indexer.AppendVector(&PosRoots(i,0));
	
	// Find the next root until one cannot find any more
	size_t RtIndx = 0;
	while (RtIndx < PosRoots.get_rows())
	{
		MatrixRow<LAINT> ThisRoot(PosRoots,RtIndx);
		MatrixRow<LAINT> ThisWeight(PosWeights,RtIndx);
		MatrixRow<LAINT> ThisWW(WhichWay,RtIndx);
		
		// Advance in each direction, if possible
		for (int i=0; i<rank; i++)
		{
			for (int j=0; j<rank; j++)
				NewWW[j] = (j != i);
			if (ThisWW[i])
			{
				for (LAINT j=1; j<=-ThisWeight[i]; j++)
				{
					// Calculate roots and weights in parallel,
					// to avoid repeated matrix.vector calculations
					copy(ThisRoot.begin(), ThisRoot.end(), NewRoot.begin());
					NewRoot[i] += j;
					copy(ThisWeight.begin(), ThisWeight.end(), NewWeight.begin());
					for (int k=0; k<rank; k++)
						NewWeight[k] += j*Cartan(i,k);
					
					// Avoid integer multiples of previous root vectors
					LAINT NRLen = 0;
					for (int k=0; k<rank; k++)
						NRLen += NewRoot[k];
					
					bool RootOK = true;
					for (LAINT NRDiv=2; NRDiv<NRLen; NRDiv++)
					{
						// Divides the whole length?
						if ((NRLen % NRDiv) != 0) continue;
						
						// Divides each member?
						bool DVI = true;
						for (int k=0; k<rank; k++)
						{
							if ((NewRoot[k] % NRDiv) != 0)
							{
								DVI = false;
								break;
							}
						}
						if (!DVI) continue;
						
						// Divide and check
						for (int k=0; k<rank; k++)
							NwRtRed[k] = NewRoot[k]/NRDiv;
						if (Indexer.CheckVector(NwRtRed))
						{
							RootOK = false;
							break;
						}
					}
					if (!RootOK) continue;
					
					// Add the root if possible
					std::pair<bool,size_t> ret = Indexer.AppendVector(NewRoot);
					if (ret.first)
					{
						// Don't add it, and mark off its direction
						// as not to be followed
						MatrixRow<LAINT> WWRow(WhichWay,RtIndx);
						for (int k=0; k<rank; k++)
							WWRow[k] &= NewWW[k];
					}
					else
					{
						// Add it
						PosRoots.AppendVector(NewRoot);
						PosWeights.AppendVector(NewWeight);
						WhichWay.AppendVector(NewWW);
					}
				}
			}
		}
		// All done; get ready to do the next one
		RtIndx++;
	}
}

void LieAlgebra::MakePosRootSum()
{
	PosRootSum.resize(rank);
	PosWeightSum.resize(rank);
	for (int i=0; i<PosRoots.get_rows(); i++)
		for (int j=0; j<rank; j++)
		{
			PosRootSum[j] += PosRoots(i,j);
			PosWeightSum[j] += PosWeights(i,j);
		}
}

// General setup
void LieAlgebra::Setup()
{
	int nc = 0;
	if (!(IsValid = MakeDynkin())) return;
	if (!(IsValid = MakeMetric())) return;
	if (!(IsValid = MakeCartan())) return;
	MakePosRoots();
	MakePosRootSum();
}


// For fast sorting
size_t LARootHashFunction::operator() (LAINT_VECTOR &vec)
{
	size_t hashcode = 0;
	for (LAINT_VECTOR::iterator it = vec.begin(); it != vec.end(); it++)
		hashcode = 17*hashcode + (*it);
	
	return hashcode;
}


// For caching the algebras

struct LieAlgebraLessThan
{
	bool operator() (const LieAlgebraParams &L1, const LieAlgebraParams &L2) const;
};

bool LieAlgebraLessThan::operator() (const LieAlgebraParams &L1, const LieAlgebraParams &L2) const
{
	if (L1.family < L2.family) return true;
	if (L1.family > L2.family) return false;
	return (L1.rank < L2.rank);
}


typedef std::map<LieAlgebraParams, LieAlgebra, LieAlgebraLessThan> LieAlgebraCacheType;
static LieAlgebraCacheType LieAlgebraCache;
static LieAlgebra InvalidLieAlgebra;


LieAlgebra &GetLieAlgebra(const LieAlgebraParams &Params)
{
	auto iter = LieAlgebraCache.find(Params);
	if (iter == LieAlgebraCache.end())
	{
		std::pair<LieAlgebraParams, LieAlgebra> rec;
		rec.first = Params;
		if (!rec.second.UseParams(Params)) return InvalidLieAlgebra;
		std::pair<LieAlgebraCacheType::iterator, bool> ret =
			LieAlgebraCache.insert(rec);
		return ret.first->second;
	}
	else
		return iter->second;
}

LieAlgebra &GetLieAlgebra(LAINT family, LAINT rank)
{
	LieAlgebraParams Params;
	Params.family = family; Params.rank = rank;
	return GetLieAlgebra(Params);
}

void ClearLieAlgebras() {LieAlgebraCache.clear();}