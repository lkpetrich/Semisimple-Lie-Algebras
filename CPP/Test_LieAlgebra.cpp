/*
	Tests the Lie-algebra class
*/

#include "LieAlgebra.h"
#include "Test_Shared.h"


void DumpDesc(LieAlgebra &la)
{
	printf("Lie Algebra family, rank: %d %d\n",la.family,la.rank);
}

void DumpDynkin(LieAlgebra &la)
{
	DynkinDiagramType &dd = la.dynkin;
	printf("Dynkin diagram: root lengths, root connections\n");
	for (int i=0; i<dd.rtlens.size(); i++)
		printf("%d: %d\n", i+1, dd.rtlens[i]);
	for (int i=0; i<dd.rtconns.size(); i++)
	{
		RootConnectionType &rtconn = dd.rtconns[i];
		printf("%d: %d %d %d\n",i,rtconn.root1,rtconn.root2,rtconn.strength);
	}
	printf("\n");
}

template<class N> void DumpIntMatrix(Matrix<N> &Mat)
{
	for (int ir=0; ir<Mat.get_rows(); ir++)
	{
		for (int ic=0; ic<Mat.get_cols(); ic++)
			printf("%6d",Mat(ir,ic));
		println();
	}
}

template<class N> void DumpFracMatrix(Matrix< Fraction<N> > &Mat)
{
	for (int ir=0; ir<Mat.get_rows(); ir++)
	{
		for (int ic=0; ic<Mat.get_cols(); ic++)
			printf("%4d/%1d",Mat(ir,ic).get_num(),Mat(ir,ic).get_den());
		println();
	}
}


int main(int argc, char **argv)
{
	if (argc <= 1)
	{
		printf("Tests what's in LieAlgebra.h\n");
		printf("Command-line args:\n");
		printf("   params -- parameter-validity-test results\n");
		printf("   dynkin -- Dynkin diagrams\n");
		printf("   matrices -- metric and Cartan matrices\n");
		printf("   posroots -- positive roots\n");
		return 0;
	}

	const bool ParamsTest = CheckArgMembership(argc,argv,"params");
	const bool DynkinTest = CheckArgMembership(argc,argv,"dynkin");
	const bool MatrixTest = CheckArgMembership(argc,argv,"matrices");
	const bool PosrootTest = CheckArgMembership(argc,argv,"posroots");
	
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
	for (int i=0; i<NumLAs; i++)
	{
		LieAlgebra &la = GetLieAlgebra(ParamList[i]);
		if (!la.IsValid) \
			if (ParamsTest)
				printf("Bad LA index, params: %d: %d %d\n", i, ParamList[i].family, ParamList[i].rank);
	}
	if (ParamsTest) println();
	
	if (DynkinTest)
	{
		for (int i=0; i<NumLAs; i++)
		{
			LieAlgebra &la = GetLieAlgebra(ParamList[i]);
			if (!la.IsValid) continue;
			DumpDesc(la);
			DumpDynkin(la);
		}
	}
	
	if (MatrixTest)
	{
		for (int i=0; i<NumLAs; i++)
		{
			LieAlgebra &la = GetLieAlgebra(ParamList[i]);
			if (!la.IsValid) continue;
			DumpDesc(la);
			printf("Metric\n");
			DumpIntMatrix(la.Metric);
			DumpFracMatrix(la.InverseMetric);
			DumpIntMatrix(la.InvMetNum);
			printf("%d\n\n",la.InvMetDen);
			printf("Cartan\n");
			DumpIntMatrix(la.Cartan);
			DumpFracMatrix(la.InverseCartan);
			DumpIntMatrix(la.InvCtnNum);
			printf("%d\n\n",la.InvCtnDen);
		}
	}
	
	if (PosrootTest)
	{
		for (int i=0; i<NumLAs; i++)
		{
			LieAlgebra &la = GetLieAlgebra(ParamList[i]);
			if (!la.IsValid) continue;
			DumpDesc(la);
			printf("Positive Roots\n");
			DumpIntMatrix(la.PosRoots);
			printf("Positive Weights\n");
			DumpIntMatrix(la.PosWeights);
		}	
	}
}
