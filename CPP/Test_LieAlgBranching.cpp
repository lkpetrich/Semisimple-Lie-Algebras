/*
	Tests the Lie-algebra subalgebras and branching rules
*/

#include "LieAlgBranching.h"
#include "Test_Shared.h"


// For writing out a Lie-algebra parameter set
template<>
struct std::formatter<LieAlgebraParams>
{
	template <typename ParseContext>
	constexpr auto parse(ParseContext &ctxt) {
		return ctxt.begin();
  	}
  	template <typename FormatContext>
	auto format(const LieAlgebraParams &Params, FormatContext &ctxt) const {
		return std::format_to(ctxt.out(), "({},{})", Params.family, Params.rank);
	}
};


// Zero out an array
template<typename CType> void Zeros(CType &Ctnr)
	{std::fill(Ctnr.begin(), Ctnr.end(), 0);}


void DumpCWL(const LACntdMaxWtList &CWL)
{
	for (const auto &Entry: CWL)
	{
		std::print("{:3}   ", Entry.Count);
		for (const auto wt: Entry.MaxWts)
			std::print(" {:3}", wt);
		println();
	}
}


// Mathematica can do this with DownValues[brancher]
void DumpBrancher(const LABrancher &Brancher)
{
	std::println("Brancher Dump");
	const LieAlgebraParams &OrigParams = Brancher.OrigAlg.Params;
	std::println("Orig algebra (family, rank) = {}", OrigParams);
	std::println("Result algebra: (family, rank) sets, number of U(1) factors");
	const LieAlgProduct &ResParams = Brancher.ResAlg.Params;
	for (const auto &Param: ResParams.ParamList)
			std::println("{}", Param);
	std::println("{}", ResParams.NumU1s);
	
	std::println("Projectors");
	LAINT ir, ic, nr, nc;
	for (const auto &Proj: Brancher.Projectors)
	{
		std::println("Params = {}", Proj.Params);
		nr = Proj.SubMatrix.get_rows(); nc = Proj.SubMatrix.get_cols();
		std::println("Submatrix size = {}, {}", nr, nc);
		for (ir=0; ir<nr; ir++)
		{
			for (ic=0; ic<nc; ic++)
			{
				const BrSubMatEntry &val = Proj.SubMatrix(ir, ic);
				std::print("   {:3} {:3}", val.get_num(), val.get_den());
			}
			std::println("");
		}
		nr = Proj.SubMatNum.get_rows(); nc = Proj.SubMatNum.get_cols();
		std::println("Numerator size = {}, {}", nr, nc);
		for (ir=0; ir<nr; ir++)
		{
			for (ic=0; ic<nc; ic++)
			{
				const LAINT &val = Proj.SubMatNum(ir, ic);
				std::print(" {:3}", val);
			}
			std::println("");
		}
		std::println("Denominator = {}", Proj.SubMatDen);
	}
	
	std::println("U(1) Indices");
	for (const auto u1i: Brancher.U1Indices)
		std::println("{}", u1i);
	
	nr = Brancher.U1SrcVecs.get_rows(); nc = Brancher.U1SrcVecs.get_cols();
	std::println("U(1) Vectors, size = {} {}", nr, nc);
	for (ir=0; ir<nr; ir++)
	{
		for (ic=0; ic<nc; ic++)
		{
			const BrSubMatEntry &val = Brancher.U1SrcVecs(ir, ic);
			std::print(" {}/{}", val.get_num(), val.get_den());
		}
		std::println("");
	}
	
	println();
}

// A stripped-down version
void DumpBrancherSubalgebras(const LABrancher &Brancher)
{
	std::println("Brancher Subalgebra Dump");
	const LieAlgebraParams &OrigParams = Brancher.OrigAlg.Params;
	std::println("Orig algebra (family, rank) = {}", OrigParams);
	std::println("Result algebra: (family, rank) sets, number of U(1) factors");
	const LieAlgProduct &ResParams = Brancher.ResAlg.Params;
	for (const auto Proj: ResParams.ParamList)
			std::println("{}", Proj);
	std::println("{}", ResParams.NumU1s);

	println();
}


void DumpRootRelated(LABrancher (*BranchMaker)(const LieAlgebraParams &, LAINT),
	const LieAlgebraParams &Params, LAINT NumMaxWts, const LAINT *MaxWtSet)
{
	std::println("Algebra: {}", Params); println();
	for (LAINT RootNo = 1; RootNo <= Params.rank; RootNo++)
	{
		std::println("Root number = {}", RootNo); println();
		const LABrancher Brancher = BranchMaker(Params, RootNo);
		// DumpBrancher(Brancher);
		std::print("Result algebra = ");
		const LieAlgProduct &ResParams = Brancher.ResAlg.Params;
		for (const auto &Param: ResParams.ParamList)
			std::print("{}  ", Param);
		std::println("{}", ResParams.NumU1s);

		for (LAINT iwt=0; iwt<NumMaxWts; iwt++)
		{
			const LAINT *MaxWts = MaxWtSet + Params.rank*iwt;
			std::print("Max Wts =");
			for (LAINT i=0; i<Params.rank; i++) std::print(" {:3}", MaxWts[i]);
			std::println("");
			LACntdMaxWtList CWL = Brancher.DoBranching(MaxWts);
			DumpCWL(CWL);
		}
		println();
	}
}


void DumpWeights(const LAINT_VECTOR &MaxWts)
{
	for (const auto mw: MaxWts)
		std::print(" {:3}", mw);
	std::println("");
}


void DumpBranchingSU(const LABrancher &Brancher)
{
	LAINT_VECTOR MaxWts;
	LACntdMaxWtList CWL;
	LAINT rank = Brancher.OrigAlg.Params.rank;
	MaxWts.resize(rank);
	
	Zeros(MaxWts);
	MaxWts[0] = 1;
	DumpWeights(MaxWts);
	CWL = Brancher.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	Zeros(MaxWts);
	MaxWts[rank-1] = 1;
	DumpWeights(MaxWts);
	CWL = Brancher.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	Zeros(MaxWts);
	if (rank > 1)
	{
		MaxWts[0] = 1;
		MaxWts[rank-1] = 1;
	}
	else
		MaxWts[0] = 2;
	DumpWeights(MaxWts);
	CWL = Brancher.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
}

void DumpBranchingSOSp(const LABrancher &Brancher, LAINT NumSpinors = 0)
{
	LAINT_VECTOR MaxWts;
	LACntdMaxWtList CWL;
	LAINT rank = Brancher.OrigAlg.Params.rank;
	MaxWts.resize(rank);
	
	Zeros(MaxWts);
	MaxWts[0] = 1;
	DumpWeights(MaxWts);
	CWL = Brancher.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	Zeros(MaxWts);
	MaxWts[0] = 2;
	DumpWeights(MaxWts);
	CWL = Brancher.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	Zeros(MaxWts);
	if (rank > 1)
		MaxWts[1] = 1;
	else
		MaxWts[0] = 2;
	DumpWeights(MaxWts);
	CWL = Brancher.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	if (NumSpinors >= 2)
	{
		Zeros(MaxWts);
		MaxWts[rank-2] = 1;
		DumpWeights(MaxWts);
		CWL = Brancher.DoBranching(MaxWts);
		DumpCWL(CWL);
		println();
	}
	
	if (NumSpinors >= 1)
	{
		Zeros(MaxWts);
		MaxWts[rank-1] = 1;
		DumpWeights(MaxWts);
		CWL = Brancher.DoBranching(MaxWts);
		DumpCWL(CWL);
		println();
	}
	
	println();
}


void DumpConjugates(const LABrancher &Br0, const LAINT *MaxWts)
{
	LACntdMaxWtList CWL = Br0.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	LAINT_VECTOR cjgixs;
	LABrancher Br00 = BrancherConjugate(Br0, cjgixs);
	CWL = Br00.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	cjgixs.push_back(2);
	LABrancher Br01 = BrancherConjugate(Br0, cjgixs);
	CWL = Br01.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	cjgixs[0] = 1;
	LABrancher Br10 = BrancherConjugate(Br0, cjgixs);
	CWL = Br10.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	cjgixs.push_back(2);
	LABrancher Br11 = BrancherConjugate(Br0, cjgixs);
	CWL = Br11.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
}


// Should be operator==

bool operator==(const LieAlgebraParams &Params1, const LieAlgebraParams &Params2)
{
	if (Params1.family != Params2.family) return false;
	if (Params1.rank != Params2.rank) return false;
	return true;
}

template<typename N> bool operator==(const Matrix<N> &M1, const Matrix<N> &M2)
{
	size_t n1 = M1.get_rows();
	if (M2.get_rows() != n1) return false;
	
	size_t n2 = M1.get_cols();
	if (M2.get_cols() != n2) return false;
	
	for (size_t i1=0; i1<n1; i1++)
		for (size_t i2=0; i2<n2; i2++)
			if (M1(i1, i2) != M2(i1, i2)) return false;
	
	return true;
}


int main(int argc, char **argv)
{
	if (argc <= 1)
	{
		std::println("Tests what's in LieAlgBranching.h");
		std::println("Command-line args:");
		std::println("   demote -- root demotions");
		std::println("   multdmt -- multiple root demotions");
		std::println("   extsplit -- extension splits");
		std::println("   voprods -- vector-outer-product reps");
		std::println("   clssxtra -- classical extra");
		std::println("   excpxtra -- exceptional extra");
		std::println("   rearr -- rearranging, concatenating, etc.");
		return 0;
	}
	
	const bool RootDemotions = CheckArgMembership(argc, argv, "demote");
	const bool MultiRootDemote = CheckArgMembership(argc, argv, "multdmt");
	const bool ExtensionSplits = CheckArgMembership(argc, argv, "extsplit");
	const bool VectorOuterProds = CheckArgMembership(argc, argv, "voprods");
	const bool ClassicalExtra = CheckArgMembership(argc, argv, "clssxtra");
	const bool ExceptionalExtra = CheckArgMembership(argc, argv, "excpxtra");
	const bool Rearranging = CheckArgMembership(argc, argv, "rearr");
	
	if (RootDemotions)
	{
		{
			const LieAlgebraParams Params = {1,6};
			const LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,1,0,0,0,0, 1,0,0,0,0,1};
			DumpRootRelated(MakeRootDemoter, Params, 3, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {2,6};
			const LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,0,0,0,1};
			DumpRootRelated(MakeRootDemoter, Params, 3, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {3,6};
			const LAINT MaxWtSet[] = {1,0,0,0,0,0, 2,0,0,0,0,0};
			DumpRootRelated(MakeRootDemoter, Params, 2, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {4,6};
			const LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1};
			DumpRootRelated(MakeRootDemoter, Params, 4, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {7,2};
			const LAINT MaxWtSet[] = {1,0, 0,1};
			DumpRootRelated(MakeRootDemoter, Params, 2, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {6,4};
			const LAINT MaxWtSet[] = {1,0,0,0, 0,0,0,1};
			DumpRootRelated(MakeRootDemoter, Params, 2, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {5,6};
			const LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1};
			DumpRootRelated(MakeRootDemoter, Params, 3, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {5,7};
			const LAINT MaxWtSet[] = {1,0,0,0,0,0,0, 0,0,0,0,0,1,0};
			DumpRootRelated(MakeRootDemoter, Params, 2, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {5,8};
			const LAINT MaxWtSet[] = {0,0,0,0,0,0,1,0};
			DumpRootRelated(MakeRootDemoter, Params, 1, MaxWtSet);
		}
	}
	
	if (MultiRootDemote)
	{
		const LieAlgebraParams ParamList[] = {{1,6}, {2,6}, {3,6}, {4,6},
			{7,2}, {6,4}, {5,6}, {5,7}, {5,8}};
		size_t NumParams = sizeof(ParamList)/sizeof(LieAlgebraParams);
		for (size_t ip=0; ip<NumParams; ip++)
		{
			const LieAlgebraParams &Params = ParamList[ip];
			std::println("Params = {}", Params);
			LAINT_VECTOR RootNos(1);
			for (LAINT ir=1; ir<=Params.rank; ir++)
			{
				RootNos[0] = ir;
				std::println(" Root = {}", ir);
				const LABrancher Brancher = MakeMultiRootDemoter(Params, RootNos);
				const LABrancher Brancher0 = MakeRootDemoter(Params, ir);
				bool AreTheSame = true;
				if (Brancher.OrigAlg.Params != Brancher0.OrigAlg.Params)
					AreTheSame = false;
				if (Brancher.Projectors.size() == Brancher0.Projectors.size())
				{
					for (LAINT i=0; i<Brancher.Projectors.size(); i++)
					{
						const LABrProjector &Proj = Brancher.Projectors[i];
						const LABrProjector &Proj0 = Brancher0.Projectors[i];
						if (Proj.Params != Proj0.Params)
							AreTheSame = false;
						if (Proj.SubMatrix != Proj0.SubMatrix)
							AreTheSame = false;
					}
				}
				else
					AreTheSame = false;
				if (Brancher.U1Indices.size() == Brancher0.U1Indices.size())
				{
					for (size_t i=0; i<Brancher.U1Indices.size(); i++)
						if (Brancher.U1Indices[i] != Brancher0.U1Indices[i])
							AreTheSame = false;
				}
				else
					AreTheSame = false;
				
				if (Brancher.U1SrcVecs != Brancher0.U1SrcVecs)
					AreTheSame = false;
				
				if (AreTheSame)
					std::println("Same");
				else
				{
					std::println("*** Not the Same ***");
					std::println("*** Multi ***");
					DumpBrancher(Brancher);
					std::println("*** Orig ***");
					DumpBrancher(Brancher0);
					std::println("*** End ***");
				}
			}
		}
		println();
		
		std::println("SU(5) fundamental rep");
		{
			const LieAlgebraParams Params = {1,4};
			const LAINT_VECTOR RootNos = {1,2,3,4};
			const LABrancher Brancher = MakeMultiRootDemoter(Params, RootNos);
			
			const LAINT Wts[] = {1,0,0,0};
			const LACntdMaxWtList CWL = Brancher.DoBranching(Wts);
		
			std::println("All U(1)'s");
			DumpCWL(CWL);
		}
		{
			const LieAlgebraParams Params = {1,4};
			const LAINT_VECTOR RootNos = {3};
			const LABrancher Brancher = MakeMultiRootDemoter(Params, RootNos);
			
			const LAINT Wts[] = {1,0,0,0};
			const LACntdMaxWtList CWL = Brancher.DoBranching(Wts);
		
			std::println("Standard Model");
			DumpCWL(CWL);
		}
		println();
		std::println("SO(10) vector rep");
		{
			const LieAlgebraParams Params = {4,5};
			const LAINT_VECTOR RootNos = {1,2,3,4,5};
			const LABrancher Brancher = MakeMultiRootDemoter(Params, RootNos);
			
			const LAINT Wts[] = {1,0,0,0,0};
			const LACntdMaxWtList CWL = Brancher.DoBranching(Wts);
		
			std::println("All U(1)'s");
			DumpCWL(CWL);
		}
		{
			const LieAlgebraParams Params = {4,5};
			const LAINT_VECTOR RootNos = {3,4};
			const LABrancher Brancher = MakeMultiRootDemoter(Params, RootNos);
			
			const LAINT Wts[] = {1,0,0,0,0};
			const LACntdMaxWtList CWL = Brancher.DoBranching(Wts);
		
			std::println("Standard Model");
			DumpCWL(CWL);
		}
	}
	
	if (ExtensionSplits)
	{
		{
			const LieAlgebraParams Params = {1,6};
			const LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,1,0,0,0,0, 1,0,0,0,0,1};
			DumpRootRelated(MakeExtensionSplitter, Params, 3, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {2,6};
			const LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,0,0,0,1};
			DumpRootRelated(MakeExtensionSplitter, Params, 3, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {3,6};
			const LAINT MaxWtSet[] = {1,0,0,0,0,0, 2,0,0,0,0,0};
			DumpRootRelated(MakeExtensionSplitter, Params, 2, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {4,6};
			const LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1};
			DumpRootRelated(MakeExtensionSplitter, Params, 4, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {7,2};
			const LAINT MaxWtSet[] = {1,0, 0,1};
			DumpRootRelated(MakeExtensionSplitter, Params, 2, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {6,4};
			const LAINT MaxWtSet[] = {1,0,0,0, 0,0,0,1};
			DumpRootRelated(MakeExtensionSplitter, Params, 2, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {5,6};
			const LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1};
			DumpRootRelated(MakeExtensionSplitter, Params, 3, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {5,7};
			const LAINT MaxWtSet[] = {1,0,0,0,0,0,0, 0,0,0,0,0,1,0};
			DumpRootRelated(MakeExtensionSplitter, Params, 2, MaxWtSet);
		}
		{
			const LieAlgebraParams Params = {5,8};
			const LAINT MaxWtSet[] = {0,0,0,0,0,0,1,0};
			DumpRootRelated(MakeExtensionSplitter, Params, 1, MaxWtSet);
		}
	}
	
	if (VectorOuterProds)
	{	
		{
			std::println("SU(30) -> SU(5)*SU(3)*SU(2)"); println();
			const LAINT_VECTOR SUOrds = {5,3,2};
			const LABrancher Brancher = SubalgMultSU(SUOrds);
			DumpBranchingSU(Brancher);
		}
		
		{
			std::println("SO(30) -> SO(6)*SO(5)"); println();
			const LAINT_VECTOR SOSpOrds = {6,5};
			const LABrancher Brancher = SubalgMultSOSp(SOSpOrds);
			DumpBranchingSOSp(Brancher);
		}

		{
			std::print("Sp(30) -> SO(5)*Sp(6)");
			const LAINT_VECTOR SOSpOrds = {5,-6};
			const LABrancher Brancher = SubalgMultSOSp(SOSpOrds);
			DumpBranchingSOSp(Brancher);
		}

		{
			std::println("SO(35) -> SO(5)*SO(7)"); println();
			const LAINT_VECTOR SOSpOrds = {5,7};
			const LABrancher Brancher = SubalgMultSOSp(SOSpOrds);
			DumpBranchingSOSp(Brancher);
		}

		{
			std::println("SO(8) -> SO(2)*SO(2)*SO(2)"); println();
			const LAINT_VECTOR SOSpOrds = {2,2,2};
			const LABrancher Brancher = SubalgMultSOSp(SOSpOrds);
			DumpBranchingSOSp(Brancher);
		}

		{
			std::println("Sp(16) -> Sp(2)*SO(2)*SO(4)"); println();
			const LAINT_VECTOR SOSpOrds = {-2,2,4};
			const LABrancher Brancher = SubalgMultSOSp(SOSpOrds);
			DumpBranchingSOSp(Brancher);
		}
	}
	
	if (ClassicalExtra)
	{
		const LAINT eornk = 5;
		for (LAINT esrk=0; esrk <= eornk-2; esrk++)
		{
			std::println("{}: SO({}) -> SO({})*SO({})",
				esrk, 2*eornk, 2*esrk+1 , 2*eornk-2*esrk-1); println();
			const LABrancher Brancher = SubalgSOEvenOdd(eornk, esrk);
			DumpBranchingSOSp(Brancher,2);
		}
		
		{
			std::println("SU(10) -> SO(10)"); println();
			const LABrancher Brancher = SubalgSUSO(10);
			DumpBranchingSU(Brancher);
		}
		
		{
			std::println("SU(11) -> SO(11)"); println();
			const LABrancher Brancher = SubalgSUSO(11);
			DumpBranchingSU(Brancher);
		}
		
		{
			std::println("SU(10) -> Sp(10)"); println();
			const LABrancher Brancher = SubalgSUSp(5);
			DumpBranchingSU(Brancher);
		}
		
		{
			const LieAlgebraParams Params = {1,8};
			std::println("Height-A1 original: {}", Params); println();
			const LABrancher Brancher = SubalgHeightA1(Params);
			LAINT_VECTOR MaxWts(Params.rank);
			for (LAINT ix=0; ix<Params.rank; ix++)
			{
				std::print("{:3}:  ", ix);
				Zeros(MaxWts);
				MaxWts[ix] = 1;
				DumpWeights(MaxWts);
				const LACntdMaxWtList CWL = Brancher.DoBranching(MaxWts);
				DumpCWL(CWL);
				println();
			}
		}
		
		{
			std::println("Vector-to-irrep: A5 -> A2 {{2,0}}"); println();
			const LieAlgebraParams DestParams = {1,2};
			const LAINT DestWts[2] = {2,0};
			const LABrancher Brancher = SubalgVector(1, DestParams, DestWts);
			DumpBranchingSU(Brancher);
		}
		
		{
			std::println("Vector-to-irrep: B10 -> B3 {{0,1,0}}"); println();
			const LieAlgebraParams DestParams = {2,3};
			const LAINT DestWts[3] = {0,1,0};
			const LABrancher Brancher = SubalgVector(2, DestParams, DestWts);
			DumpBranchingSOSp(Brancher,1);
		}
		
		{
			std::println("Vector-to-irrep: C10 -> B2 {{0,3}}"); println();
			const LieAlgebraParams DestParams = {2,2};
			const LAINT DestWts[2] = {0,3};
			const LABrancher Brancher = SubalgVector(3, DestParams, DestWts);
			DumpBranchingSOSp(Brancher,0);
		}
		
		{
			std::println("Vector-to-irrep: D5 -> B2 {{0,2}}"); println();
			const LieAlgebraParams DestParams = {2,2};
			const LAINT DestWts[2] = {0,2};
			const LABrancher Brancher = SubalgVector(4, DestParams, DestWts);
			DumpBranchingSOSp(Brancher,2);
		}
	}
	
	if (ExceptionalExtra)
	{
		{
			std::println("B3 / SO(7) -> G2"); println();
			const LABrancher Brancher = SubalgExtra(LA_BR_SUBALG_EXTRA_B3G2);
			LAINT_VECTOR MaxWts(3);
			for (LAINT ix=0; ix<3; ix++)
			{
				std::print("{:3}:  ", ix);
				Zeros(MaxWts);
				MaxWts[ix] = 1;
				DumpWeights(MaxWts);
				const LACntdMaxWtList CWL = Brancher.DoBranching(MaxWts);
				DumpCWL(CWL);
				println();
			}
		}
		
		{
			std::println("E8 -> G2*F4"); println();
			LABrancher Brancher = SubalgExtra(LA_BR_SUBALG_EXTRA_E8G2F4);
			const LAINT MaxWts[] = {0,0,0,0, 0,0,1,0};
			const LACntdMaxWtList CWL = Brancher.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();	
		}
	}
	
	if (Rearranging)
	{
		{
			std::println("SO(10) spinor breakdown"); println();
			
			const LAINT MaxWts[] = {0,0,0,1,0};
			const LieAlgebraParams Params0 = {4,5};
			const LABrancher Br0 = MakeRootDemoter(Params0,5);
			const LACntdMaxWtList CWL0 = Br0.DoBranching(MaxWts);
			DumpCWL(CWL0);
			println();
			
			const LieAlgebraParams Params1 = {1,4};
			const LABrancher Br1 = MakeRootDemoter(Params1,3);
			const LABrancher Br01 = ConcatBranchers(Br0,1, Br1);
			const LACntdMaxWtList CWL1 = Br01.DoBranching(MaxWts);
			DumpCWL(CWL1);
			println();
			
			LAINT_VECTOR mrts;
			mrts.push_back(5);
			mrts.push_back(3);
			const LABrancher Br2 = MakeMultiRootDemoter(Params0, mrts);
			const LACntdMaxWtList CWL2 = Br2.DoBranching(MaxWts);
			DumpCWL(CWL2);
			println();
			
			const LABrancher Br3 = MakeExtensionSplitter(Params0,3);
			const LACntdMaxWtList CWL3 = Br3.DoBranching(MaxWts);
			DumpCWL(CWL3);
			println();
			
			LAINT_VECTOR sos;
			sos.push_back(2);
			sos.push_back(2);
			const LABrancher Br31 = SubalgMultSOSp(sos);
			const LABrancher Br311 = ConcatBranchers(Br3,2, Br31);
			const LACntdMaxWtList CWL31 = Br311.DoBranching(MaxWts);
			DumpCWL(CWL31);
			println();
		}
		
		{
			std::println("A(1) - B(1) - C(1) relabeling"); println();
			LAINT_VECTOR mrts;
			mrts.push_back(-2);
			mrts.push_back(-2);
			const LABrancher Br0 = SubalgMultSOSp(mrts);
			DumpBrancherSubalgebras(Br0);
			const LABrancher Br1 = BrancherRenameA1B1C1(Br0,1,1);
			DumpBrancherSubalgebras(Br1);
			const LABrancher Br2 = BrancherRenameA1B1C1(Br1,2,2);
			DumpBrancherSubalgebras(Br2);
			const LABrancher Br3 = BrancherRenameA1B1C1(Br2,1,3);
			DumpBrancherSubalgebras(Br3);
			const LAINT MaxWts[2] = {2,3};
			LACntdMaxWtList CWL;
			CWL = Br0.DoBranching(MaxWts);
			DumpCWL(CWL);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			CWL = Br2.DoBranching(MaxWts);
			DumpCWL(CWL);
			CWL = Br3.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
		}
		
		{
			std::println("B(2) - C(2) relabeling"); println();
			const LieAlgebraParams Params = {2,4};
			const LABrancher Br0 = MakeExtensionSplitter(Params,2);
			const LABrancher Br1 = BrancherRenameB2C2(Br0,2);
			const LABrancher Br2 = BrancherRenameB2C2(Br1,2);
			const LAINT MaxWts[4] = {1,0,0,0};
			LACntdMaxWtList CWL;
			CWL = Br0.DoBranching(MaxWts);
			DumpCWL(CWL);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			CWL = Br2.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
		}
		
		{
			std::println("A(3) - D(3) relabeling"); println();
			const LieAlgebraParams Params = {4,4};
			const LABrancher Br0 = MakeRootDemoter(Params,3);
			const LABrancher Br1 = BrancherRenameA3D3(Br0,1);
			const LABrancher Br2 = BrancherRenameA3D3(Br1,1);
			const LAINT MaxWts[4] = {1,0,0,0};
			LACntdMaxWtList CWL;
			CWL = Br0.DoBranching(MaxWts);
			DumpCWL(CWL);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			CWL = Br2.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
		}
		
		{
			std::println("D(2) - A(1)^2 split and join"); println();
			const LieAlgebraParams Params = {4,4};
			const LABrancher Br0 = MakeExtensionSplitter(Params,2);
			const LABrancher Br1 = BrancherSplitD2(Br0,1);
			const LABrancher Br2 = BrancherSplitD2(Br1,3);
			const LABrancher Br3 = BrancherJoin2A1(Br2,1,3);
			const LABrancher Br4 = BrancherJoin2A1(Br3,1,2);
			const LAINT MaxWts[4] = {1,0,0,0};
			LACntdMaxWtList CWL;
			CWL = Br0.DoBranching(MaxWts);
			DumpCWL(CWL);
			CWL = Br2.DoBranching(MaxWts);
			DumpCWL(CWL);
			CWL = Br4.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
		}
		
		{
			std::println("Conjugate reps"); println();
			
			std::println("SU(9) -> SU(4)*SU(5)*U(1)");
			const LieAlgebraParams Params1 = {1,8};
			const LABrancher Br1 = MakeRootDemoter(Params1,4);
			const LAINT MaxWts1[] = {1,0,0,0,0,0,0,0};
			DumpConjugates(Br1, MaxWts1);
			
			std::println("SO(16) -> SU(3)*SO(10)*U(1)");
			const LieAlgebraParams Params2 = {4,8};
			const LABrancher Br2 = MakeRootDemoter(Params2,3);
			const LAINT MaxWts2[] = {0,0,0,0,0,0,0,1};
			DumpConjugates(Br2, MaxWts2);
			
			std::println("SO(17) -> SU(3)*SO(11)*U(1)");
			const LieAlgebraParams Params3 = {2,8};
			const LABrancher Br3 = MakeRootDemoter(Params3,3);
			const LAINT MaxWts3[] = {0,0,0,0,0,0,0,1};
			DumpConjugates(Br3, MaxWts3);
			
			std::println("E(8) -> E(6)*SU(2)*U(1)");
			const LieAlgebraParams Params4 = {5,8};
			const LABrancher Br4 = MakeRootDemoter(Params4,6);
			const LAINT MaxWts4[] = {0,0,0,0,0,0,1,0};
			DumpConjugates(Br4, MaxWts4);
		}
		
		{
			std::println("SO(8) vector-spinor interchangeability"); println();
			const LieAlgebraParams Params = {4,4};
			LABrancher Br0 = SubalgSelf(Params);
			const LAINT MaxWts[] = {0,0,1,2};
			LAINT_VECTOR newrts;
			for (LAINT i=1; i<=4;i++) newrts.push_back(i);
			
			LACntdMaxWtList CWL = Br0.DoBranching(MaxWts);
			DumpCWL(CWL);
			
			LABrancher Br1 = BrancherConjgD4(Br0,1, newrts);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			
			newrts[0] = 1; newrts[2] = 4; newrts[3] = 3;
			Br1 = BrancherConjgD4(Br0,1, newrts);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			
			newrts[0] = 3; newrts[2] = 1; newrts[3] = 4;
			Br1 = BrancherConjgD4(Br0,1, newrts);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			
			newrts[0] = 3; newrts[2] = 4; newrts[3] = 1;
			Br1 = BrancherConjgD4(Br0,1, newrts);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			
			newrts[0] = 4; newrts[2] = 1; newrts[3] = 3;
			Br1 = BrancherConjgD4(Br0,1, newrts);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			
			newrts[0] = 4; newrts[2] = 3; newrts[3] = 1;
			Br1 = BrancherConjgD4(Br0,1, newrts);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			
			println();
		}
		
		{
			std::println("SO(10) -> SO(6)*SU(2)*SU(2) rearrangement"); println();
			const LieAlgebraParams Params = {4,5};
			LABrancher Br00 = MakeExtensionSplitter(Params,3);
			LABrancher Br0 = BrancherSplitD2(Br00,2);
			const LAINT MaxWts[] = {0,0,0,0,1};
			LAINT_VECTOR neword;
			for (LAINT i=1; i<=3;i++) neword.push_back(i);
			
			LACntdMaxWtList CWL = Br0.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
			
			LABrancher Br1 = BrancherRearrange(Br0, neword);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
			
			neword[0] = 1; neword[1] = 3; neword[2] = 2;
			Br1 = BrancherRearrange(Br0, neword);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
			
			neword[0] = 2; neword[1] = 1; neword[2] = 3;
			Br1 = BrancherRearrange(Br0, neword);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
			
			neword[0] = 2; neword[1] = 3; neword[2] = 1;
			Br1 = BrancherRearrange(Br0, neword);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
			
			neword[0] = 3; neword[1] = 1; neword[2] = 2;
			Br1 = BrancherRearrange(Br0, neword);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
			
			neword[0] = 3; neword[1] = 2; neword[2] = 1;
			Br1 = BrancherRearrange(Br0, neword);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
		}
	}
}