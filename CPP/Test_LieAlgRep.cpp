/*
	Tests the Lie-algebra representation finding
*/

#include "LieAlgRep.h"
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


void MakeLowMaxWts(LAINT_MATRIX &mws, size_t n)
{
	mws.resize(0, n);
	std::vector<LAINT> mw(n);
	
	fill(mw.begin(), mw.end(),0);
	mws.AppendVector(mw);
	
	for (size_t i=0; i<n; i++)
	{
		fill(mw.begin(), mw.end(),0);
		mw[i] += 1;
		mws.AppendVector(mw);
	}
	
	for (size_t i=0; i<n; i++)
		for (size_t j=i; j<n; j++)
		{
			fill(mw.begin(), mw.end(),0);
			mw[i] += 1;
			mw[j] += 1;
			mws.AppendVector(mw);
		}
}

void DumpTotalDegens(const LieAlgebraParams &Params)
{
	std::println("Lie algebra: {}", Params);
	LieAlgebra &la = GetLieAlgebra(Params);
	
	LAINT_MATRIX mws;
	MakeLowMaxWts(mws, Params.rank);
	
	for (size_t i=0; i<mws.get_rows(); i++)
	{
		MatrixRow<LAINT> mw(mws, i);
		
		for (size_t j=0; j<Params.rank; j++)
			std::print("{} ", mw[j]);
		
		TDINT res = TotalDegen(la,&mw[0]);
		std::println("  {}" , res.get_str());
	}
	
	println();
}

const size_t LARPWtSize = 20;
struct LARepParams
{
	LieAlgebraParams LAP;
	LAINT Wt[LARPWtSize];
};

void DumpParams(const LARepParams &Params)
{
	std::print("Params: {}  ", Params.LAP);
	for (int i=0; i<Params.LAP.rank; i++)
		std::print(" {:3}", Params.Wt[i]);
	println();
}

void DumpRepObject(const LieAlgRep &Rep)
{
	size_t rows = Rep.Roots.get_rows();
	for (int i=0; i<rows; i++)
	{
		LAINT Degen = Rep.Degens[i];
		const LAINT_MATRIX_ROW_CONST Root(Rep.Roots, i), Weight(Rep.Weights, i);
		std::print("{:3}   ", Degen);
		for (int j=0; j<Rep.Roots.get_cols(); j++)
			std::print("{:3} ", Root[j]);
		std::print("  ");
		for (int j=0; j<Rep.Weights.get_cols(); j++)
			std::print("{:3} ", Weight[j]);
		println();
	}
}

void DumpCWL(const LACntdMaxWtList &CWL)
{
	for (const auto &Entry: CWL)
	{
		std::print("{:3}   ", Entry.Count);
		for (auto wt: Entry.MaxWts)
			std::print(" {:3}", wt);
		println();
	}
}

void CompareReps(const LieAlgRep &Rep1, const LieAlgRep &Rep2)
{
	LieAlgRep SR1, SR2;
	Rep1.Export(SR1);
	Rep2.Export(SR2);
	
	size_t kn = SR1.Degens.size();
	size_t knx = SR2.Degens.size();
	size_t n = SR1.vlen;
	size_t nx = SR2.vlen;
	if (knx != kn) goto dump;
	if (nx != n) goto dump;
	
	for (int k=0; k < kn; k++)
	{
		if (SR2.Degens[k] != SR1.Degens[k]) goto dump;
	
		for (int j = 0; j < n; j++)
		{
			if (SR2.Roots(k, j) != SR1.Roots(k, j)) goto dump;
			if (SR2.Weights(k, j) != SR1.Weights(k, j)) goto dump;
		}
	}
	
	std::println("Reps equal");
	return;
	dump:
	
	std::println("Rep 1");
	DumpRepObject(SR1);
	std::println("Rep 2");
	DumpRepObject(SR2);
}

void DumpIndivRep(const LARepParams &Params)
{
	std::print("Params: {}  ", Params.LAP);
	for (int i=0; i<Params.LAP.rank; i++)
		std::print(" {}", Params.Wt[i]);
	println();
	
	LieAlgebra &la = GetLieAlgebra(Params.LAP);
	clock_t t0 = clock();
	LieAlgRep Rep0;
	LieAlgebraRepDirect(Rep0, la, Params.Wt);
	clock_t t1 = clock();
	std::println("Time: {}", (t1-t0)/double(CLOCKS_PER_SEC));
	LieAlgRepPtr &RepPtr = GetRepObject(RO_REP, Params.LAP, Params.Wt);
	LieAlgRep &Rep = *RepPtr;
	clock_t t2 = clock();
	std::println("Time: {}", (t2-t1)/double(CLOCKS_PER_SEC));
	CompareReps(Rep0, Rep);
	/*
	DumpRepObject(Rep);
	*/
	long ntot; sum(ntot, Rep.Degens);
	std::println("Total = {}", ntot);
	RepProperties Properties;
	Properties.Use(Params.LAP, Params.Wt);
	std::println("True total = {}\n", Properties.Size.get_str());
	switch(Properties.Reality)
	{
	case REP_REALITY_REAL:
		std::println("Real");
		break;
	case REP_REALITY_PSEUDOREAL:
		std::println("Pseudoreal");
		break;
	case REP_REALITY_COMPLEX:
		std::println("Complex");
		break;
	}
	std::print("Conserved:");
	for (const auto val: Properties.Conserved)
		std::print(" {}", val);
	std::println("");
	std::println("Height: {}", Properties.Height);
	
	println();
}


int main(int argc, char **argv)
{
	if (argc <= 1)
	{
		std::println("Tests what's in LieAlgRep.h");
		std::println("Command-line args:");
		std::println("   totals -- total sizes");
		std::println("   indivs -- individual reps");
		std::println("   vecspn -- vector and spinor reps of SO(n)");
		std::println("   algprod -- algebra-product algebras");
		std::println("   hdlr -- rep handler: creation and decomposition");
		std::println("   weylorb -- Weyl orbits");
		return 0;
	}

	const bool TotalTest = CheckArgMembership(argc, argv, "totals");
	const bool IndivTest = CheckArgMembership(argc, argv, "indivs");
	const bool VecSpnTest = CheckArgMembership(argc, argv, "vecspn");
	const bool AlgProdTest = CheckArgMembership(argc, argv, "algprod");
	const bool HandlerTest = CheckArgMembership(argc, argv, "hdlr");
	const bool WeylOrbitTest = CheckArgMembership(argc, argv, "weylorb");

	if (TotalTest)
	{
		// Algebras to verify
		const LieAlgebraParams ParamList[] =
			{{1,1}, {1,5}, {2,5}, {3,5}, {4,5}, {7,2}, {6,4}, {5,6}, {5,7}, {5,8}};
		const int NumLAs = sizeof(ParamList)/sizeof(LieAlgebraParams);
		
		for (int i=0; i<NumLAs; i++)
			DumpTotalDegens(ParamList[i]);
	}
	
	if (IndivTest)
	{
		const LARepParams ParamList[] = 
		{
		{{1,1},{0}}, {{1,1},{1}}, {{1,1},{2}}, {{1,1},{3}}, {{1,1},{4}},
		{{1,2},{0,0}}, {{1,2},{1,0}}, {{1,2},{0,1}}, {{1,2},{2,0}}, {{1,2},{0,2}}, {{1,2},{1,1}}, {{1,2},{3,0}}, {{1,2},{0,3}}, {{1,2},{2,1}}, {{1,2},{1,2}},
		{{2,2},{0,0}}, {{2,2},{1,0}}, {{2,2},{0,1}}, {{2,2},{2,0}}, {{2,2},{0,2}}, {{2,2},{1,1}}, {{2,2},{3,0}}, {{2,2},{0,3}}, {{2,2},{2,1}}, {{2,2},{1,2}},
		{{3,2},{0,0}}, {{3,2},{1,0}}, {{3,2},{0,1}}, {{3,2},{2,0}}, {{3,2},{0,2}}, {{3,2},{1,1}}, {{3,2},{3,0}}, {{3,2},{0,3}}, {{3,2},{2,1}}, {{3,2},{1,2}},
		{{4,2},{0,0}}, {{4,2},{1,0}}, {{4,2},{0,1}}, {{4,2},{2,0}}, {{4,2},{0,2}}, {{4,2},{1,1}}, {{4,2},{3,0}}, {{4,2},{0,3}}, {{4,2},{2,1}}, {{4,2},{1,2}},
		{{7,2},{0,0}}, {{7,2},{1,0}}, {{7,2},{0,1}}, {{7,2},{2,0}}, {{7,2},{0,2}}, {{7,2},{1,1}}, {{7,2},{3,0}}, {{7,2},{0,3}}, {{7,2},{2,1}}, {{7,2},{1,2}},
		{{1,4},{0,0,0,0}}, {{1,4},{1,0,0,0}}, {{1,4},{0,1,0,0}}, {{1,4},{0,0,1,0}}, {{1,4},{0,0,0,1}},
		{{2,4},{0,0,0,0}}, {{2,4},{1,0,0,0}}, {{2,4},{0,1,0,0}}, {{2,4},{0,0,1,0}}, {{2,4},{0,0,0,1}},
		{{3,4},{0,0,0,0}}, {{3,4},{1,0,0,0}}, {{3,4},{0,1,0,0}}, {{3,4},{0,0,1,0}}, {{3,4},{0,0,0,1}},
		{{4,4},{0,0,0,0}}, {{4,4},{1,0,0,0}}, {{4,4},{0,1,0,0}}, {{4,4},{0,0,1,0}}, {{4,4},{0,0,0,1}},
		{{6,4},{0,0,0,0}}, {{6,4},{1,0,0,0}}, {{6,4},{0,1,0,0}}, {{6,4},{0,0,1,0}}, {{6,4},{0,0,0,1}},
		{{5,6},{0,0,0,0,0,0}}, {{5,6},{1,0,0,0,0,0}}, {{5,6},{0,0,0,0,1,0}}, {{5,6},{0,0,0,0,0,1}},
		{{5,7},{0,0,0,0,0,0,0}}, {{5,7},{1,0,0,0,0,0,0}}, {{5,7},{0,0,0,0,0,1,0}},
		{{5,8},{0,0,0,0,0,0,0,0}}, {{5,8},{0,0,0,0,0,0,1,0}}, {{5,8},{1,0,0,0,0,0,0,0}},
		{{5,8},{0,0,0,0,0,0,2,0}}, {{5,8},{0,0,0,0,0,1,0,0}}, {{5,8},{0,0,0,0,0,0,0,1}},
		};
		const int NumReps = sizeof(ParamList)/sizeof(LARepParams);
		
		for (int i=0; i<NumReps; i++)
			DumpIndivRep(ParamList[i]);
	}
	
	if (VecSpnTest)
	{
		// Look for Bott periodicities
		// Properties of spinors for n space dimensions should repeat for n + k*8
		// That's n + k*4 for rank
		LARepParams ZeroParams;
		memset(&ZeroParams,0, sizeof(ZeroParams));
		LARepParams Params;
		
		std::println("-- Vectors --"); println();
		for (size_t i=1; i<=16; i++)
		{
			std::println("SO({})",2*i); println();
			
			if (i == 1)
			{
			}
			else if (i == 2)
			{
				Params = ZeroParams;
				Params.LAP.family = 4;
				Params.LAP.rank = i;
				Params.Wt[0] = 1;
				Params.Wt[1] = 1;
				DumpIndivRep(Params);
			}
			else
			{
				Params = ZeroParams;
				Params.LAP.family = 4;
				Params.LAP.rank = i;
				Params.Wt[0] = 1;
				DumpIndivRep(Params);
			}
			
			std::println("SO({})",2*i+1); println();
			
			if (i == 1)
			{
				Params = ZeroParams;
				Params.LAP.family = 2;
				Params.LAP.rank = i;
				Params.Wt[0] = 2;
				DumpIndivRep(Params);
			}
			else
			{
				Params = ZeroParams;
				Params.LAP.family = 2;
				Params.LAP.rank = i;
				Params.Wt[0] = 1;
				DumpIndivRep(Params);
			}
		}
		std::println("-- Adjoints --"); println();
		for (size_t i=1; i<=16; i++)
		{
			std::println("SO({})",2*i); println();
			
			if (i == 1)
			{
			}
			else if (i == 2)
			{
				Params = ZeroParams;
				Params.LAP.family = 4;
				Params.LAP.rank = i;
				Params.Wt[0] = 2;
				Params.Wt[1] = 0;
				DumpIndivRep(Params);
				Params = ZeroParams;
				Params.LAP.family = 4;
				Params.LAP.rank = i;
				Params.Wt[0] = 0;
				Params.Wt[1] = 2;
				DumpIndivRep(Params);
			}
			else if (i == 3)
			{
				Params = ZeroParams;
				Params.LAP.family = 4;
				Params.LAP.rank = i;
				Params.Wt[1] = 1;
				Params.Wt[2] = 1;
				DumpIndivRep(Params);
			}
			else
			{
				Params = ZeroParams;
				Params.LAP.family = 4;
				Params.LAP.rank = i;
				Params.Wt[1] = 1;
				DumpIndivRep(Params);
			}
			
			std::println("SO({})",2*i+1); println();
			
			if (i == 1)
			{
				Params = ZeroParams;
				Params.LAP.family = 2;
				Params.LAP.rank = i;
				Params.Wt[0] = 2;
				DumpIndivRep(Params);
			}
			else if (i == 2)
			{
				Params = ZeroParams;
				Params.LAP.family = 2;
				Params.LAP.rank = i;
				Params.Wt[1] = 2;
				DumpIndivRep(Params);
			}
			else
			{
				Params = ZeroParams;
				Params.LAP.family = 2;
				Params.LAP.rank = i;
				Params.Wt[1] = 1;
				DumpIndivRep(Params);
			}
		}
		std::println("-- Spinors --"); println();
		for (size_t i=1; i<=16; i++)
		{
			const std::string EvenType[4] = {"Real", "Complex", "Pseudoreal", "Complex"};
			std::println("SO({}) -- {}",2*i, EvenType[i%4]); println();
			
			if (i == 1)
			{
			}
			else
			{
				Params = ZeroParams;
				Params.LAP.family = 4;
				Params.LAP.rank = i;
				Params.Wt[i-2] = 1;
				DumpIndivRep(Params);
				Params = ZeroParams;
				Params.LAP.family = 4;
				Params.LAP.rank = i;
				Params.Wt[i-1] = 1;
				DumpIndivRep(Params);
			}
			
			const std::string OddType[4] = {"Real", "Pseudoreal", "Pseudoreal", "Real"};
			std::println("SO({}) -- {}",2*i+1, OddType[i%4]); println();
			
			{
				Params = ZeroParams;
				Params.LAP.family = 2;
				Params.LAP.rank = i;
				Params.Wt[i-1] = 1;
				DumpIndivRep(Params);
			}
		}
	}
	
	if (AlgProdTest)
	{
		std::println("-- Product --");
		std::println("Alg: {{{{1,3}}, {{1,2}}}}");
		std::println("Rep: {{1,0,1, 1,1, -1,2,-3}}");
		println();
		LieAlgProduct LAP;
		const LieAlgebraParams AlgParamList[] = {{1,3},{1,2}};
		const int NumAlgs = sizeof(AlgParamList)/sizeof(LieAlgebraParams);
		for (int i=0; i<NumAlgs; i++)
			LAP.ParamList.push_back(AlgParamList[i]);
		LAP.NumU1s = 3;
		const LAINT MaxWts[] = {1,0,1, 1,1, -1,2,-3};
		
		LieAlgRep Rep;
		LieAlgProdRepresentation(Rep, LAP, RO_REP, MaxWts);
		DumpRepObject(Rep);
	}
	
	if (HandlerTest)
	{
		// Create sums of irreps and then decompose them
		std::println("-- Rep-handler test --"); println();
		LieAlgRep Rep;
		LACntdMaxWtList CWL;
		
		LieAlgebraParams p1 = {1,1};
		LAINT wt11[] = {1};
		LAINT wt12[] = {3};
		
		LASnglRepHandler RH1(p1);
		LieAlgRepPtr rp11 = RH1.GetRepPtr(RO_REP, wt11);
		LieAlgRepPtr rp12 = RH1.GetRepPtr(RO_REP, wt12);
		LieAlgRepBuilder bld1(p1.rank);
		bld1.Import(*rp11, true);
		bld1.Import(*rp12, true);
		bld1.Import(*rp11, true);
		bld1.Export(Rep);
		DumpRepObject(Rep);
		CWL = RH1.ExtractWts(RO_REP, bld1);
		DumpCWL(CWL);
		println();

		LieAlgebraParams p2 = {1,2};
		LAINT wt21[] = {1,0};
		LAINT wt22[] = {0,2};
		
		LASnglRepHandler RH2(p2);
		LieAlgRepPtr rp21 = RH2.GetRepPtr(RO_REP, wt21);
		LieAlgRepPtr rp22 = RH2.GetRepPtr(RO_REP, wt22);
		LieAlgRepBuilder bld2(p2.rank);
		bld2.Import(*rp21, true);
		bld2.Import(*rp22, true);
		bld2.Import(*rp21, true);
		bld2.Export(Rep);
		DumpRepObject(Rep);
		CWL = RH2.ExtractWts(RO_REP, bld2);
		DumpCWL(CWL);
		println();
		
		LieAlgProduct p3;
		LieAlgebraParams p31 = {1,2}, p32 = {1,1};
		p3.ParamList.push_back(p31);
		p3.ParamList.push_back(p32);
		p3.NumU1s = 1;
		LAINT wt31[] = {1,0,1,2};
		LAINT wt32[] = {0,2,3,2};
		LAINT wt33[] = {1,0,3,2};
		LAINT wt34[] = {0,2,1,2};
		LAINT wt35[] = {1,0,1,-3};
		
		LAProdRepHandler RH3(p3);
		LieAlgRepPtr rp31 = RH3.GetRepPtr(RO_REP, wt31);
		LieAlgRepPtr rp32 = RH3.GetRepPtr(RO_REP, wt32);
		LieAlgRepPtr rp33 = RH3.GetRepPtr(RO_REP, wt33);
		LieAlgRepPtr rp34 = RH3.GetRepPtr(RO_REP, wt34);
		LieAlgRepPtr rp35 = RH3.GetRepPtr(RO_REP, wt35);
		LieAlgRepBuilder bld3(p3.get_rank());
		bld3.Import(*rp31, true);
		bld3.Import(*rp32, true);
		bld3.Import(*rp33, true);
		bld3.Import(*rp34, true);
		bld3.Import(*rp35, true);
		bld3.Import(*rp32, true);
		bld3.Import(*rp33, true);
		bld3.Import(*rp34, true);
		bld3.Import(*rp35, true);
		bld3.Import(*rp33, true);
		bld3.Import(*rp34, true);
		bld3.Import(*rp35, true);
		bld3.Import(*rp34, true);
		bld3.Import(*rp35, true);
		bld3.Import(*rp35, true);
		bld3.Export(Rep);
		DumpRepObject(Rep);
		CWL = RH3.ExtractWts(RO_REP, bld3);
		DumpCWL(CWL);
	}
	
	if (WeylOrbitTest)
	{
		std::println("-- Weyl-orbit tests --"); println();
		
		std::println("Explicit orbits:"); println();

		const LARepParams XpParamList[] = 
		{
		{{1,4},{2,1,0,0}}, {{2,4},{2,1,0,0}}, {{3,4},{2,1,0,0}}, {{4,4},{2,1,0,0}},
		{{1,4},{4,3,2,1}}, {{2,4},{4,3,2,1}}, {{3,4},{4,3,2,1}}, {{4,4},{4,3,2,1}},
		};
		const int XpNumReps = sizeof(XpParamList)/sizeof(LARepParams);
		
		for (int i=0; i<XpNumReps; i++)
		{
			const LARepParams &Params = XpParamList[i];
			DumpParams(Params);
			LieAlgebra &la = GetLieAlgebra(Params.LAP);
			LieAlgRep worep;
			LieAlgebraWeylOrbitGeneral(worep, la, Params.Wt);
			LieAlgRep worepx;
			LieAlgebraWeylOrbitExplicit(worepx, la, Params.Wt);
			CompareReps(worep, worepx);
		}
		
		println();
		
		std::println("Rep Weyl orbits:"); println();
		
		const LARepParams RWOParamList[] = 
		{
		{{1,4},{1,0,0,0}}, {{1,4},{1,0,0,1}}, {{1,4},{2,1,0,0}},
		{{2,4},{2,0,0,1}}, {{3,4},{1,0,0,1}}, {{4,4},{0,0,1,2}}, 
		{{7,2},{2,1}}, {{6,4},{1,0,0,1}}, {{5,6},{1,0,0,0,0,1}},
		{{5,7},{2,0,0,0,0,0,0}}, {{5,8},{0,0,0,0,0,0,2,0}}
		};
		const int RWONumReps = sizeof(RWOParamList)/sizeof(LARepParams);
		
		for (int i=0; i<RWONumReps; i++)
		{
			const LARepParams &Params = RWOParamList[i];
			DumpParams(Params);
			LieAlgebra &la = GetLieAlgebra(Params.LAP);
			LieAlgRep rwrep;
			LieAlgebraRepWeylOrbits(rwrep, la, Params.Wt);
			DumpRepObject(rwrep);
		}
		
		std::println("Reps from orbits:"); println();
		const LARepParams ORParamList[] = 
		{
		{{1,4},{1,0,1,0}}, {{2,4},{1,1,0,1}}, {{3,4},{1,2,0,0}},
		{{4,4},{0,1,1,0}}, {{7,2},{3,2}}, {{6,4}, {1,1,0,0}},
		{{5,6},{2,0,0,0,0,1}}, {{5,7},{2,0,0,0,0,0,0}},
		{{5,8},{0,0,0,0,0,0,2,0}}
		};
		const int ORNumReps = sizeof(ORParamList)/sizeof(LARepParams);
		
		for (int i=0; i<ORNumReps; i++)
		{
			const LARepParams &Params = ORParamList[i];
			DumpParams(Params);
			LieAlgebra &la = GetLieAlgebra(Params.LAP);
			LieAlgRep orrep, drrep;
			LieAlgebraRepByWOs(orrep, la, Params.Wt);
			LieAlgebraRepDirect(drrep, la, Params.Wt);
			CompareReps(orrep, drrep);
		}
	}
}
