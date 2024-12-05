/*
	Tests the Lie-algebra class
*/

#include "LieAlgebra.h"
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


void DumpDesc(const LieAlgebra &la)
{
	std::println("Lie algebra: {}", (LieAlgebraParams &)la);
}

void DumpDynkin(const LieAlgebra &la)
{
	const DynkinDiagramType &dd = la.dynkin;
	std::println("Dynkin diagram: root lengths, root connections");
	for (int i=0; i<dd.rtlens.size(); i++)
		std::println("{}: {}", i+1, dd.rtlens[i]);
	for (int i=0; i<dd.rtconns.size(); i++)
	{
		const RootConnectionType &rtconn = dd.rtconns[i];
		std::println("{}: {} {} {}", i, rtconn.root1, rtconn.root2, rtconn.strength);
	}
	println();
}

template<typename N> void DumpIntMatrix(const Matrix<N> &Mat)
{
	for (int ir=0; ir<Mat.get_rows(); ir++)
	{
		for (int ic=0; ic<Mat.get_cols(); ic++)
			std::print("{:6}", Mat(ir, ic));
		println();
	}
}

template<typename N> void DumpFracMatrix(const Matrix< Fraction<N> > &Mat)
{
	for (int ir=0; ir<Mat.get_rows(); ir++)
	{
		for (int ic=0; ic<Mat.get_cols(); ic++)
			std::print("{:4}/{}", Mat(ir, ic).get_num(), Mat(ir, ic).get_den());
		println();
	}
}


int main(int argc, char **argv)
{
	if (argc <= 1)
	{
		std::println("Tests what's in LieAlgebra.h");
		std::println("Command-line args:");
		std::println("   params -- parameter-validity-test results");
		std::println("   dynkin -- Dynkin diagrams");
		std::println("   matrices -- metric and Cartan matrices");
		std::println("   posroots -- positive roots");
		return 0;
	}

	const bool ParamsTest = CheckArgMembership(argc, argv, "params");
	const bool DynkinTest = CheckArgMembership(argc, argv, "dynkin");
	const bool MatrixTest = CheckArgMembership(argc, argv, "matrices");
	const bool PosrootTest = CheckArgMembership(argc, argv, "posroots");
	
	// Algebras to verify
	const LieAlgebraParams ParamList[] = {
		// Good parameter values, for doing a good sample of all the Lie algebras
		// The first and a later classical one for all 4 classical ones,
		// then the 5 exceptional ones
		{1,1}, {1,5}, {2,1}, {2,5}, {3,1}, {3,5}, {4,2}, {4,5},
		{7,2}, {6,4}, {5,6}, {5,7}, {5,8},
		// Bad parameter values, for testing bad-parameter rejection
		{0,1}, {1,0}, {2,0}, {3,0}, {4,1},
		{5,5}, {5,9}, {6,3}, {6,5}, {7,1}, {7,3},
		{8,1}
	};
	const int NumLAs = sizeof(ParamList)/sizeof(LieAlgebraParams);
	
	if (ParamsTest)
	{
		for (int i=0; i<NumLAs; i++)
		{
			std::print("{}: ", i);
			const LieAlgebra &la = GetLieAlgebra(ParamList[i]);
			if (la.IsValid)
				std::print("Good: ");
			else
				std::print("Bad: ");
			std::println("{}", ParamList[i]);
		}
		println();
	}
	
	if (DynkinTest)
	{
		for (int i=0; i<NumLAs; i++)
		{
			const LieAlgebra &la = GetLieAlgebra(ParamList[i]);
			if (!la.IsValid) continue;
			DumpDesc(la);
			DumpDynkin(la);
		}
	}
	
	if (MatrixTest)
	{
		for (int i=0; i<NumLAs; i++)
		{
			const LieAlgebra &la = GetLieAlgebra(ParamList[i]);
			if (!la.IsValid) continue;
			DumpDesc(la);
			std::println("Metric");
			DumpIntMatrix(la.Metric);
			DumpFracMatrix(la.InverseMetric);
			DumpIntMatrix(la.InvMetNum);
			std::println("{}", la.InvMetDen);
			println();
			std::println("Cartan");
			DumpIntMatrix(la.Cartan);
			DumpFracMatrix(la.InverseCartan);
			DumpIntMatrix(la.InvCtnNum);
			std::println("{}", la.InvCtnDen);
			println();
		}
	}
	
	if (PosrootTest)
	{
		for (int i=0; i<NumLAs; i++)
		{
			const LieAlgebra &la = GetLieAlgebra(ParamList[i]);
			if (!la.IsValid) continue;
			DumpDesc(la);
			std::println("Positive Roots");
			DumpIntMatrix(la.PosRoots);
			std::println("Positive Weights");
			DumpIntMatrix(la.PosWeights);
		}	
	}
}
