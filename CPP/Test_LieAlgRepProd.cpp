/*
	Tests the Lie-algebra rep-product finding
*/

#include "LieAlgRepProd.h"
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

void DumpCWLList(LAYDCWLList &CWLList)
{
	for (LAYDCWLList::iterator CWLIter = CWLList.begin(); CWLIter != CWLList.end(); CWLIter++)
	{
		printf("Mult = %d\n",CWLIter->YDE.Mult);
		printf("YD =");
		for (YoungDiagram::iterator YDIter = CWLIter->YDE.YD.begin(); YDIter != CWLIter->YDE.YD.end(); YDIter++)
			printf(" %3d",*YDIter);
		printf("\n");
		DumpCWL(CWLIter->CWL);
	}
}

int main(int argc, char **argv)
{
	if (argc <= 1)
	{
		printf("Tests what's in LieAlgRepProd.h\n");
		printf("Command-line args:\n");
		printf("   prods -- rep products\n");
		printf("   powers -- rep powers (plethysms) \n");
		printf("   pwrsyms -- (anti)symmetric rep powers \n");
		return 0;
	}

	const bool ProdTest = CheckArgMembership(argc,argv,"prods");
	const bool PowerTest = CheckArgMembership(argc,argv,"powers");
	const bool PwrSymTest = CheckArgMembership(argc,argv,"pwrsyms");
	
	LieAlgebraParams p1 = {1,4};
	LASnglRepHandler h1(p1);
	LAINT wt11[] = {1,0,0,0};
	LAINT wt12[] = {1,0,0,1};
	
	// For a light quark: SU(3) color, SU(3) flavor, SU(2)/SO(3) spin
	LieAlgProduct p2;
	LieAlgebraParams pp2[] = {{1,2},{1,2},{1,1}};
	for (int i=0; i<3; i++)
		p2.ParamList.push_back(pp2[i]);
	p2.NumU1s = 0;
	LAProdRepHandler h2(p2);
	LAINT wt2[] = {1,0,1,0,1};
	LAINT wtx2[] = {0,1,0,1,1};
	
	if (ProdTest)
	{
		LACntdMaxWtList CWL;
		
		/*
			Mathematica equivalent:
			MakeLieAlgebra[su3, {1, 2}];
			MatrixForm @ DecomposeRepProduct[su3, {1, 2}, {3, 1}]
		*/
		LieAlgebraParams p0 = {1,2};
		LASnglRepHandler h0(p0);
		LAINT wt01[] = {1,2};
		LAINT wt02[] = {3,1};
		
		printf("Product of SU(3) {1,2} and {3,1}\n");
		CWL = DecomposeRepProduct(h0,wt01,wt02);
		DumpCWL(CWL);
		
		printf("Quark and Antiquark: Color, Light Flavors, Spin\n");
		CWL = DecomposeRepProduct(h2,wt2,wtx2);
		DumpCWL(CWL);
		
		println();
	}
	
	if (PowerTest)
	{
		/*
			Mathematica equivalent:
			MakeLieAlgebra[la, {1, 4}];
			MatrixForm[Table[MatrixForm /@ DecomposeRepPower[la, {1, 0, 0, 0}, p], {p, 1, 4}]]
			MatrixForm[Table[MatrixForm /@ DecomposeRepPower[la, {1, 0, 0, 1}, p], {p, 1, 4}]]
			Three Quarks:
			MakeLieAlgebra[su2, {1, 1}]; MakeLieAlgebra[su3, {1, 2}]; 
			MatrixForm /@ DecomposeAlgProdRepPower[{su3, su3, su2}, {{1, 0}, {1, 0}, {1}}, 3]
		*/
		LAYDCWLList CWLList;
		
		printf("Powers of SU(5) {1,0,0,0} decomposed by Young-diagram symmetry type\n");
		for (int p=0; p<=4; p++)
		{
			printf("Power = %d\n",p);
			CWLList = DecomposeRepPower(h1,wt11,p);
			DumpCWLList(CWLList);
		}
		printf("Powers of SU(5) {1,0,0,1} decomposed by Young-diagram symmetry type\n");
		for (int p=0; p<=4; p++)
		{
			printf("Power = %d\n",p);
			CWLList = DecomposeRepPower(h1,wt12,p);
			DumpCWLList(CWLList);
		}
		printf("Three Quarks\n");
		CWLList = DecomposeRepPower(h2,wt2,3);
		DumpCWLList(CWLList);
		println();
	}
	
	if (PwrSymTest)
	{		
		/*
			Mathematica equivalent:
			MakeLieAlgebra[la, {1, 4}];
			MatrixForm[Table[MatrixForm @ DecomposeRepPwrSym[la, {1, 0, 0, 0}, p, s], {p, 1, 4}, {s, 1, -1, -2}]]
			MatrixForm[Table[MatrixForm @ DecomposeRepPwrSym[la, {1, 0, 0, 1}, p, s], {p, 1, 4}, {s, 1, -1, -2}]]
			Three Quarks:
			MakeLieAlgebra[su2, {1, 1}]; MakeLieAlgebra[su3, {1, 2}]; 
			MatrixForm /@ Table[DecomposeAlgProdRepPwrSym[{su3, su3, su2}, {{1, 0}, {1, 0}, {1}}, 3, s], {s, 1, -1, -2}]
		*/
		LACntdMaxWtList CWL;
		
		printf("Powers of SU(5) {1,0,0,0} -- (anti)symmetric\n");
		for (int p=0; p<=4; p++)
		{
			printf("Power = %d\n",p);
			for (int s=1; s>=-1; s-=2)
			{
				printf("Symm = %d\n",s);
				CWL = DecomposeRepPwrSym(h1,wt11,p,s);
				DumpCWL(CWL);
			}
		}
		printf("Powers of SU(5) {1,0,0,1} -- (anti)symmetric\n");
		for (int p=0; p<=4; p++)
		{
			printf("Power = %d\n",p);
			for (int s=1; s>=-1; s-=2)
			{
				printf("Symm = %d\n",s);
				CWL = DecomposeRepPwrSym(h1,wt12,p,s);
				DumpCWL(CWL);
			}
		}
		printf("Three Quarks\n");
		for (int s=1; s>=-1; s-=2)
		{
			printf("Symm = %d\n",s);
			CWL = DecomposeRepPwrSym(h2,wt2,3,s);
			DumpCWL(CWL);
		}
		println();
	}
}
