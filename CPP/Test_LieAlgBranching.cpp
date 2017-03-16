/*
	Tests the Lie-algebra subalgebras and branching rules
*/

#include "LieAlgBranching.h"
#include "Test_Shared.h"


void DumpCWL(LACntdMaxWtList &CWL)
{
	for (LACntdMaxWtList::iterator wcei = CWL.begin(); wcei != CWL.end(); wcei++)
	{
		LACntdMaxWtEntry &Entry = *wcei;
		printf("%3d   ",Entry.Count);
		for (vector<LAINT>::iterator wti = Entry.MaxWts.begin(); wti != Entry.MaxWts.end(); wti++)
			printf(" %3d",*wti);
		println();
	}
}


// Mathematica can do this with DownValues[brancher]
void DumpBrancher(LABrancher &Brancher)
{
	printf("Brancher Dump\n");
	LieAlgebraParams &OrigParams = Brancher.OrigAlg.Params;
	printf("Orig algebra family, rank = %d %d\n", OrigParams.family, OrigParams.rank);
	printf("Result algebra: (family, rank) sets, number of U(1) factors\n");
	LieAlgProduct &ResParams = Brancher.ResAlg.Params;
	for (vector<LieAlgebraParams>::iterator PIter = ResParams.ParamList.begin();
		PIter != ResParams.ParamList.end(); PIter++)
			printf("%d %d\n",PIter->family,PIter->rank);
	printf("%d\n",ResParams.NumU1s);
	
	printf("Projectors\n");
	LAINT ir, ic, nr, nc;
	for (vector<LABrProjector>::iterator PIter = Brancher.Projectors.begin();
		PIter != Brancher.Projectors.end(); PIter++)
	{
		printf("Params = %d %d\n",PIter->Params.family,PIter->Params.rank);
		nr = PIter->SubMatrix.get_rows(); nc = PIter->SubMatrix.get_cols();
		printf("Submatrix size = %d, %d\n", nr, nc);
		for (ir=0; ir<nr; ir++)
		{
			for (ic=0; ic<nc; ic++)
			{
				BrSubMatEntry &val = PIter->SubMatrix(ir,ic);
				printf("   %3d %3d",val.get_num(),val.get_den());
			}
			printf("\n");
		}
		nr = PIter->SubMatNum.get_rows(); nc = PIter->SubMatNum.get_cols();
		printf("Numerator size = %d, %d\n", nr, nc);
		for (ir=0; ir<nr; ir++)
		{
			for (ic=0; ic<nc; ic++)
			{
				LAINT &val = PIter->SubMatNum(ir,ic);
				printf(" %3d",val);
			}
			printf("\n");
		}
		printf("Denominator = %d\n", PIter->SubMatDen);
	}
	
	printf("U(1) Indices\n");
	for (vector<LAINT>::iterator uiiter = Brancher.U1Indices.begin();
		uiiter != Brancher.U1Indices.end(); uiiter++)
	printf("%d\n",*uiiter);
	
	nr = Brancher.U1SrcVecs.get_rows(); nc = Brancher.U1SrcVecs.get_cols();
	printf("U(1) Vectors, size = %d %d\n", nr, nc);
	for (ir=0; ir<nr; ir++)
	{
		for (ic=0; ic<nc; ic++)
		{
			BrSubMatEntry &val = Brancher.U1SrcVecs(ir,ic);
			printf(" %3d/%3d",val.get_num(),val.get_den());
		}
		printf("\n");
	}
	
	println();
}

// A stripped-down version
void DumpBrancherSubalgebras(LABrancher &Brancher)
{
	printf("Brancher Subalgebra Dump\n");
	LieAlgebraParams &OrigParams = Brancher.OrigAlg.Params;
	printf("Orig algebra family, rank = %d %d\n", OrigParams.family, OrigParams.rank);
	printf("Result algebra: (family, rank) sets, number of U(1) factors\n");
	LieAlgProduct &ResParams = Brancher.ResAlg.Params;
	for (vector<LieAlgebraParams>::iterator PIter = ResParams.ParamList.begin();
		PIter != ResParams.ParamList.end(); PIter++)
			printf("%d %d\n",PIter->family,PIter->rank);
	printf("%d\n",ResParams.NumU1s);

	println();
}


void DumpRootRelated(LABrancher (*BranchMaker)(const LieAlgebraParams &, LAINT),
	LieAlgebraParams &Params, LAINT NumMaxWts, const LAINT *MaxWtSet)
{
	printf("Algebra family, rank = %d %d\n\n",Params.family,Params.rank);
	for (LAINT RootNo = 1; RootNo <= Params.rank; RootNo++)
	{
		printf("Root number = %d\n\n",RootNo);
		LABrancher Brancher = BranchMaker(Params, RootNo);
		// DumpBrancher(Brancher);
		printf("Result algebra = ");
		LieAlgProduct &ResParams = Brancher.ResAlg.Params;
		for (vector<LieAlgebraParams>::iterator PIter = ResParams.ParamList.begin();
			PIter != ResParams.ParamList.end(); PIter++)
				printf("%d %d   ",PIter->family,PIter->rank);
		printf("%d\n",ResParams.NumU1s);

		for (LAINT iwt=0; iwt<NumMaxWts; iwt++)
		{
			const LAINT *MaxWts = MaxWtSet + Params.rank*iwt;
			printf("Max Wts =");
			for (LAINT i=0; i<Params.rank; i++) printf(" %3d",MaxWts[i]);
			printf("\n");
			LACntdMaxWtList CWL = Brancher.DoBranching(MaxWts);
			DumpCWL(CWL);
		}
		println();
	}
}


void DumpWeights(vector<LAINT> &MaxWts)
{
	for (vector<LAINT>::iterator mwiter = MaxWts.begin();
		mwiter != MaxWts.end(); mwiter++)
			printf(" %3d",*mwiter);
	printf("\n");
}

void DumpBranchingSU(LABrancher &Brancher)
{
	vector<LAINT> MaxWts;
	LACntdMaxWtList CWL;
	LAINT rank = Brancher.OrigAlg.Params.rank;
	MaxWts.resize(rank);
	
	fill(MaxWts.begin(),MaxWts.end(),0);
	MaxWts[0] = 1;
	DumpWeights(MaxWts);
	CWL = Brancher.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	fill(MaxWts.begin(),MaxWts.end(),0);
	MaxWts[rank-1] = 1;
	DumpWeights(MaxWts);
	CWL = Brancher.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	fill(MaxWts.begin(),MaxWts.end(),0);
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

void DumpBranchingSOSp(LABrancher &Brancher, LAINT NumSpinors = 0)
{
	vector<LAINT> MaxWts;
	LACntdMaxWtList CWL;
	LAINT rank = Brancher.OrigAlg.Params.rank;
	MaxWts.resize(rank);
	
	fill(MaxWts.begin(),MaxWts.end(),0);
	MaxWts[0] = 1;
	DumpWeights(MaxWts);
	CWL = Brancher.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	fill(MaxWts.begin(),MaxWts.end(),0);
	MaxWts[0] = 2;
	DumpWeights(MaxWts);
	CWL = Brancher.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	fill(MaxWts.begin(),MaxWts.end(),0);
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
		fill(MaxWts.begin(),MaxWts.end(),0);
		MaxWts[rank-2] = 1;
		DumpWeights(MaxWts);
		CWL = Brancher.DoBranching(MaxWts);
		DumpCWL(CWL);
		println();
	}
	
	if (NumSpinors >= 1)
	{
		fill(MaxWts.begin(),MaxWts.end(),0);
		MaxWts[rank-1] = 1;
		DumpWeights(MaxWts);
		CWL = Brancher.DoBranching(MaxWts);
		DumpCWL(CWL);
		println();
	}
	
	println();
}


void DumpConjugates(LABrancher &Br0, const LAINT *MaxWts)
{
	LACntdMaxWtList CWL = Br0.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	vector<LAINT> cjgixs;
	LABrancher Br00 = BrancherConjugate(Br0,cjgixs);
	CWL = Br00.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	cjgixs.push_back(2);
	LABrancher Br01 = BrancherConjugate(Br0,cjgixs);
	CWL = Br01.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	cjgixs[0] = 1;
	LABrancher Br10 = BrancherConjugate(Br0,cjgixs);
	CWL = Br10.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
	
	cjgixs.push_back(2);
	LABrancher Br11 = BrancherConjugate(Br0,cjgixs);
	CWL = Br11.DoBranching(MaxWts);
	DumpCWL(CWL);
	println();
}


bool AreEqual(LieAlgebraParams &Params1, LieAlgebraParams &Params2)
{
	if (Params1.family != Params2.family) return false;
	if (Params1.rank != Params2.rank) return false;
	return true;
}

template<class N> bool AreEqual(Matrix<N> &M1, Matrix<N> &M2)
{
	size_t n1 = M1.get_rows();
	if (M2.get_rows() != n1) return false;
	
	size_t n2 = M1.get_cols();
	if (M2.get_cols() != n2) return false;
	
	for (size_t i1=0; i1<n1; i1++)
		for (size_t i2=0; i2<n2; i2++)
			if (M1(i1,i2) != M2(i1,i2)) return false;
	
	return true;
}


int main(int argc, char **argv)
{
	if (argc <= 1)
	{
		printf("Tests what's in LieAlgBranching.h\n");
		printf("Command-line args:\n");
		printf("   demote -- root demotions\n");
		printf("   multdmt -- multiple root demotions\n");
		printf("   extsplit -- extension splits\n");
		printf("   voprods -- vector-outer-product reps\n");
		printf("   clssxtra -- classical extra\n");
		printf("   excpxtra -- exceptional extra\n");
		printf("   rearr -- rearranging, concatenating, etc.\n");
		return 0;
	}
	
	const bool RootDemotions = CheckArgMembership(argc,argv,"demote");
	const bool MultiRootDemote = CheckArgMembership(argc,argv,"multdmt");
	const bool ExtensionSplits = CheckArgMembership(argc,argv,"extsplit");
	const bool VectorOuterProds = CheckArgMembership(argc,argv,"voprods");
	const bool ClassicalExtra = CheckArgMembership(argc,argv,"clssxtra");
	const bool ExceptionalExtra = CheckArgMembership(argc,argv,"excpxtra");
	const bool Rearranging = CheckArgMembership(argc,argv,"rearr");
	
	if (RootDemotions)
	{
		{
			LieAlgebraParams Params = {1,6};
			LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,1,0,0,0,0, 1,0,0,0,0,1};
			DumpRootRelated(MakeRootDemoter, Params, 3, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {2,6};
			LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,0,0,0,1};
			DumpRootRelated(MakeRootDemoter, Params, 3, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {3,6};
			LAINT MaxWtSet[] = {1,0,0,0,0,0, 2,0,0,0,0,0};
			DumpRootRelated(MakeRootDemoter, Params, 2, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {4,6};
			LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1};
			DumpRootRelated(MakeRootDemoter, Params, 4, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {7,2};
			LAINT MaxWtSet[] = {1,0, 0,1};
			DumpRootRelated(MakeRootDemoter, Params, 2, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {6,4};
			LAINT MaxWtSet[] = {1,0,0,0, 0,0,0,1};
			DumpRootRelated(MakeRootDemoter, Params, 2, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {5,6};
			LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1};
			DumpRootRelated(MakeRootDemoter, Params, 3, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {5,7};
			LAINT MaxWtSet[] = {1,0,0,0,0,0,0, 0,0,0,0,0,1,0};
			DumpRootRelated(MakeRootDemoter, Params, 2, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {5,8};
			LAINT MaxWtSet[] = {0,0,0,0,0,0,1,0};
			DumpRootRelated(MakeRootDemoter, Params, 1, MaxWtSet);
		}
	}
	
	if (MultiRootDemote)
	{
		LieAlgebraParams ParamList[] = {{1,6}, {2,6}, {3,6}, {4,6},
			{7,2}, {6,4}, {5,6}, {5,7}, {5,8}};
		size_t NumParams = sizeof(ParamList)/sizeof(LieAlgebraParams);
		for (size_t ip=0; ip<NumParams; ip++)
		{
			LieAlgebraParams &Params = ParamList[ip];
			printf("Params = %d %d\n",Params.family,Params.rank);
			vector<LAINT> RootNos(1);
			for (LAINT ir=1; ir<=Params.rank; ir++)
			{
				RootNos[0] = ir;
				printf(" Root = %d\n",ir);
				LABrancher Brancher = MakeMultiRootDemoter(Params,RootNos);
				LABrancher Brancher0 = MakeRootDemoter(Params,ir);
				bool AreTheSame = true;
				if (!AreEqual(Brancher.OrigAlg.Params, Brancher0.OrigAlg.Params))
					AreTheSame = false;
				if (Brancher.Projectors.size() == Brancher0.Projectors.size())
				{
					for (LAINT i=0; i<Brancher.Projectors.size(); i++)
					{
						LABrProjector &Proj = Brancher.Projectors[i];
						LABrProjector &Proj0 = Brancher0.Projectors[i];
						if (!AreEqual(Proj.Params,Proj0.Params))
							AreTheSame = false;
						if (!AreEqual(Proj.SubMatrix,Proj0.SubMatrix))
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
				
				if (!AreEqual(Brancher.U1SrcVecs,Brancher0.U1SrcVecs))
					AreTheSame = false;
				
				if (AreTheSame)
					printf("Same\n");
				else
				{
					printf("*** Not the Same ***\n");
					printf("*** Multi ***\n");
					DumpBrancher(Brancher);
					printf("*** Orig ***\n");
					DumpBrancher(Brancher0);
					printf("*** End ***\n");
				}
			}
		}
		println();
		
		printf("SU(5) fundamental rep\n");
		{
			LieAlgebraParams Params = {1,4};
			vector<LAINT> RootNos;
			RootNos.push_back(1);
			RootNos.push_back(2);
			RootNos.push_back(3);
			RootNos.push_back(4);
			LABrancher Brancher = MakeMultiRootDemoter(Params,RootNos);
			
			LAINT Wts[] = {1,0,0,0};
			LACntdMaxWtList CWL = Brancher.DoBranching(Wts);
		
			printf("All U(1)'s\n");
			DumpCWL(CWL);
		}
		{
			LieAlgebraParams Params = {1,4};
			vector<LAINT> RootNos;
			RootNos.push_back(3);
			LABrancher Brancher = MakeMultiRootDemoter(Params,RootNos);
			
			LAINT Wts[] = {1,0,0,0};
			LACntdMaxWtList CWL = Brancher.DoBranching(Wts);
		
			printf("Standard Model\n");
			DumpCWL(CWL);
		}
	}
	
	if (ExtensionSplits)
	{
		{
			LieAlgebraParams Params = {1,6};
			LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,1,0,0,0,0, 1,0,0,0,0,1};
			DumpRootRelated(MakeExtensionSplitter, Params, 3, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {2,6};
			LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,0,0,0,1};
			DumpRootRelated(MakeExtensionSplitter, Params, 3, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {3,6};
			LAINT MaxWtSet[] = {1,0,0,0,0,0, 2,0,0,0,0,0};
			DumpRootRelated(MakeExtensionSplitter, Params, 2, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {4,6};
			LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1};
			DumpRootRelated(MakeExtensionSplitter, Params, 4, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {7,2};
			LAINT MaxWtSet[] = {1,0, 0,1};
			DumpRootRelated(MakeExtensionSplitter, Params, 2, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {6,4};
			LAINT MaxWtSet[] = {1,0,0,0, 0,0,0,1};
			DumpRootRelated(MakeExtensionSplitter, Params, 2, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {5,6};
			LAINT MaxWtSet[] = {1,0,0,0,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1};
			DumpRootRelated(MakeExtensionSplitter, Params, 3, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {5,7};
			LAINT MaxWtSet[] = {1,0,0,0,0,0,0, 0,0,0,0,0,1,0};
			DumpRootRelated(MakeExtensionSplitter, Params, 2, MaxWtSet);
		}
		{
			LieAlgebraParams Params = {5,8};
			LAINT MaxWtSet[] = {0,0,0,0,0,0,1,0};
			DumpRootRelated(MakeExtensionSplitter, Params, 1, MaxWtSet);
		}
	}
	
	if (VectorOuterProds)
	{	
		{
			printf("SU(30) -> SU(5)*SU(3)*SU(2)\n\n");
			
			vector<LAINT> SUOrds;
			SUOrds.push_back(5);
			SUOrds.push_back(3);
			SUOrds.push_back(2);
			LABrancher Brancher = SubalgMultSU(SUOrds);
			DumpBranchingSU(Brancher);
		}
		
		{
			printf("SO(30) -> SO(6)*SO(5)\n\n");
			vector<LAINT> SOSpOrds;
			SOSpOrds.push_back(6);
			SOSpOrds.push_back(5);
			LABrancher Brancher = SubalgMultSOSp(SOSpOrds);
			DumpBranchingSOSp(Brancher);
		}

		{
			printf("Sp(30) -> SO(5)*Sp(6)\n\n");
			vector<LAINT> SOSpOrds;
			SOSpOrds.push_back(5);
			SOSpOrds.push_back(-6);
			LABrancher Brancher = SubalgMultSOSp(SOSpOrds);
			DumpBranchingSOSp(Brancher);
		}

		{
			printf("SO(35) -> SO(5)*SO(7)\n\n");
			vector<LAINT> SOSpOrds;
			SOSpOrds.push_back(5);
			SOSpOrds.push_back(7);
			LABrancher Brancher = SubalgMultSOSp(SOSpOrds);
			DumpBranchingSOSp(Brancher);
		}

		{
			printf("SO(8) -> SO(2)*SO(2)*SO(2)\n\n");
			vector<LAINT> SOSpOrds;
			SOSpOrds.push_back(2);
			SOSpOrds.push_back(2);
			SOSpOrds.push_back(2);
			LABrancher Brancher = SubalgMultSOSp(SOSpOrds);
			DumpBranchingSOSp(Brancher);
		}

		{
			printf("Sp(16) -> Sp(2)*SO(2)*SO(4)\n\n");
			vector<LAINT> SOSpOrds;
			SOSpOrds.push_back(-2);
			SOSpOrds.push_back(2);
			SOSpOrds.push_back(4);
			LABrancher Brancher = SubalgMultSOSp(SOSpOrds);
			DumpBranchingSOSp(Brancher);
		}
	}
	
	if (ClassicalExtra)
	{
		const LAINT eornk = 5;
		for (LAINT esrk=0; esrk <= eornk-2; esrk++)
		{
			printf("%d: SO(%d) -> SO(%d)*SO(%d)\n\n",esrk,2*eornk,2*esrk+1 , 2*eornk- 2*esrk -1);
			LABrancher Brancher = SubalgSOEvenOdd(eornk,esrk);
			DumpBranchingSOSp(Brancher,2);
		}
		
		{
			printf("SU(10) -> SO(10)\n\n");
			LABrancher Brancher = SubalgSUSO(10);
			DumpBranchingSU(Brancher);
		}
		
		{
			printf("SU(11) -> SO(11)\n\n");
			LABrancher Brancher = SubalgSUSO(11);
			DumpBranchingSU(Brancher);
		}
		
		{
			printf("SU(10) -> Sp(10)\n\n");
			LABrancher Brancher = SubalgSUSp(5);
			DumpBranchingSU(Brancher);
		}
		
		LieAlgebraParams Params = {1,8};
		printf("Height-A1 original: family, rank = %d, %d\n\n",Params.family,Params.rank);
		LABrancher Brancher = SubalgHeightA1(Params);
		vector<LAINT> MaxWts(Params.rank);
		for (LAINT ix=0; ix<Params.rank; ix++)
		{
			printf("%3d:  ",ix);
			fill(MaxWts.begin(),MaxWts.end(),0);
			MaxWts[ix] = 1;
			DumpWeights(MaxWts);
			LACntdMaxWtList CWL = Brancher.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
		}
	}
	
	if (ExceptionalExtra)
	{
		{
			printf("B3 / SO(7) -> G2\n\n");
			LABrancher Brancher = SubalgExtra(LA_BR_SUBALG_EXTRA_B3G2);
			vector<LAINT> MaxWts(3);
			for (LAINT ix=0; ix<3; ix++)
			{
				printf("%3d:  ",ix);
				fill(MaxWts.begin(),MaxWts.end(),0);
				MaxWts[ix] = 1;
				DumpWeights(MaxWts);
				LACntdMaxWtList CWL = Brancher.DoBranching(MaxWts);
				DumpCWL(CWL);
				println();
			}
		}
		
		{
			printf("E8 -> G2*F4\n\n");
			LABrancher Brancher = SubalgExtra(LA_BR_SUBALG_EXTRA_E8G2F4);
			LAINT MaxWts[] = {0,0,0,0, 0,0,1,0};
			LACntdMaxWtList CWL = Brancher.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();	
		}
	}
	
	if (Rearranging)
	{
		{
			printf("SO(10) spinor breakdown\n\n");
			
			const LAINT MaxWts[] = {0,0,0,1,0};
			const LieAlgebraParams Params0 = {4,5};
			LABrancher Br0 = MakeRootDemoter(Params0,5);
			LACntdMaxWtList CWL0 = Br0.DoBranching(MaxWts);
			DumpCWL(CWL0);
			println();
			
			const LieAlgebraParams Params1 = {1,4};
			LABrancher Br1 = MakeRootDemoter(Params1,3);
			LABrancher Br01 = ConcatBranchers(Br0,1,Br1);
			LACntdMaxWtList CWL1 = Br01.DoBranching(MaxWts);
			DumpCWL(CWL1);
			println();
			
			vector<LAINT> mrts;
			mrts.push_back(5);
			mrts.push_back(3);
			LABrancher Br2 = MakeMultiRootDemoter(Params0,mrts);
			LACntdMaxWtList CWL2 = Br2.DoBranching(MaxWts);
			DumpCWL(CWL2);
			println();
			
			LABrancher Br3 = MakeExtensionSplitter(Params0,3);
			LACntdMaxWtList CWL3 = Br3.DoBranching(MaxWts);
			DumpCWL(CWL3);
			println();
			
			vector<LAINT> sos;
			sos.push_back(2);
			sos.push_back(2);
			LABrancher Br31 = SubalgMultSOSp(sos);
			LABrancher Br311 = ConcatBranchers(Br3,2,Br31);
			LACntdMaxWtList CWL31 = Br311.DoBranching(MaxWts);
			DumpCWL(CWL31);
			println();
		}
		
		{
			printf("A(1) - B(1) - C(1) relabeling\n\n");
			vector<LAINT> mrts;
			mrts.push_back(-2);
			mrts.push_back(-2);
			LABrancher Br0 = SubalgMultSOSp(mrts);
			DumpBrancherSubalgebras(Br0);
			LABrancher Br1 = BrancherRenameA1B1C1(Br0,1,1);
			DumpBrancherSubalgebras(Br1);
			LABrancher Br2 = BrancherRenameA1B1C1(Br1,2,2);
			DumpBrancherSubalgebras(Br2);
			LABrancher Br3 = BrancherRenameA1B1C1(Br2,1,3);
			DumpBrancherSubalgebras(Br3);
			LAINT MaxWts[2] = {2,3};
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
			printf("B(2) - C(2) relabeling\n\n");
			const LieAlgebraParams Params = {2,4};
			LABrancher Br0 = MakeExtensionSplitter(Params,2);
			LABrancher Br1 = BrancherRenameB2C2(Br0,2);
			LABrancher Br2 = BrancherRenameB2C2(Br1,2);
			LAINT MaxWts[4] = {1,0,0,0};
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
			printf("A(3) - D(3) relabeling\n\n");
			const LieAlgebraParams Params = {4,4};
			LABrancher Br0 = MakeRootDemoter(Params,3);
			LABrancher Br1 = BrancherRenameA3D3(Br0,1);
			LABrancher Br2 = BrancherRenameA3D3(Br1,1);
			LAINT MaxWts[4] = {1,0,0,0};
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
			printf("D(2) - A(1)^2 split and join\n\n");
			const LieAlgebraParams Params = {4,4};
			LABrancher Br0 = MakeExtensionSplitter(Params,2);
			LABrancher Br1 = BrancherSplitD2(Br0,1);
			LABrancher Br2 = BrancherSplitD2(Br1,3);
			LABrancher Br3 = BrancherJoin2A1(Br2,1,3);
			LABrancher Br4 = BrancherJoin2A1(Br3,1,2);
			LAINT MaxWts[4] = {1,0,0,0};
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
			printf("Conjugate reps\n\n");
			
			printf("SU(9) -> SU(4)*SU(5)*U(1)\n");
			const LieAlgebraParams Params1 = {1,8};
			LABrancher Br1 = MakeRootDemoter(Params1,4);
			const LAINT MaxWts1[] = {1,0,0,0,0,0,0,0};
			DumpConjugates(Br1,MaxWts1);
			
			printf("SO(16) -> SU(3)*SO(10)*U(1)\n");
			const LieAlgebraParams Params2 = {4,8};
			LABrancher Br2 = MakeRootDemoter(Params2,3);
			const LAINT MaxWts2[] = {0,0,0,0,0,0,0,1};
			DumpConjugates(Br2,MaxWts2);
			
			printf("SO(17) -> SU(3)*SO(11)*U(1)\n");
			const LieAlgebraParams Params3 = {2,8};
			LABrancher Br3 = MakeRootDemoter(Params3,3);
			const LAINT MaxWts3[] = {0,0,0,0,0,0,0,1};
			DumpConjugates(Br3,MaxWts3);
			
			printf("E(8) -> E(6)*SU(2)*U(1)\n");
			const LieAlgebraParams Params4 = {5,8};
			LABrancher Br4 = MakeRootDemoter(Params4,6);
			const LAINT MaxWts4[] = {0,0,0,0,0,0,1,0};
			DumpConjugates(Br4,MaxWts4);
		}
		
		{
			printf("SO(8) vector-spinor interchangeability\n\n");
			const LieAlgebraParams Params = {4,4};
			LABrancher Br0 = SubalgSelf(Params);
			const LAINT MaxWts[] = {0,0,1,2};
			vector<LAINT> newrts;
			for (LAINT i=1; i<=4;i++) newrts.push_back(i);
			
			LACntdMaxWtList CWL = Br0.DoBranching(MaxWts);
			DumpCWL(CWL);
			
			LABrancher Br1 = BrancherConjgD4(Br0,1,newrts);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			
			newrts[0] = 1; newrts[2] = 4; newrts[3] = 3;
			Br1 = BrancherConjgD4(Br0,1,newrts);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			
			newrts[0] = 3; newrts[2] = 1; newrts[3] = 4;
			Br1 = BrancherConjgD4(Br0,1,newrts);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			
			newrts[0] = 3; newrts[2] = 4; newrts[3] = 1;
			Br1 = BrancherConjgD4(Br0,1,newrts);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			
			newrts[0] = 4; newrts[2] = 1; newrts[3] = 3;
			Br1 = BrancherConjgD4(Br0,1,newrts);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			
			newrts[0] = 4; newrts[2] = 3; newrts[3] = 1;
			Br1 = BrancherConjgD4(Br0,1,newrts);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			
			println();
		}
		
		{
			printf("SO(10) -> SO(6)*SU(2)*SU(2) rearrangement\n\n");
			const LieAlgebraParams Params = {4,5};
			LABrancher Br00 = MakeExtensionSplitter(Params,3);
			LABrancher Br0 = BrancherSplitD2(Br00,2);
			const LAINT MaxWts[] = {0,0,0,0,1};
			vector<LAINT> neword;
			for (LAINT i=1; i<=3;i++) neword.push_back(i);
			
			LACntdMaxWtList CWL = Br0.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
			
			LABrancher Br1 = BrancherRearrange(Br0,neword);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
			
			neword[0] = 1; neword[1] = 3; neword[2] = 2;
			Br1 = BrancherRearrange(Br0,neword);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
			
			neword[0] = 2; neword[1] = 1; neword[2] = 3;
			Br1 = BrancherRearrange(Br0,neword);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
			
			neword[0] = 2; neword[1] = 3; neword[2] = 1;
			Br1 = BrancherRearrange(Br0,neword);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
			
			neword[0] = 3; neword[1] = 1; neword[2] = 2;
			Br1 = BrancherRearrange(Br0,neword);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
			
			neword[0] = 3; neword[1] = 2; neword[2] = 1;
			Br1 = BrancherRearrange(Br0,neword);
			CWL = Br1.DoBranching(MaxWts);
			DumpCWL(CWL);
			println();
		}
	}
}