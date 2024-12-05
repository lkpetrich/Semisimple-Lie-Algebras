/*
	Tests the Lie-algebra rep-product finding
*/

#include "LieAlgRepProd.h"
#include "Test_Shared.h"


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

void DumpCWLList(const LAYDCWLList &CWLList)
{
	for (const auto &ListEntry: CWLList)
	{
		std::println("Mult = {}", ListEntry.YDE.Mult);
		std::print("YD =");
		for (const auto YD: ListEntry.YDE.YD)
			std::print(" {:3}", YD);
		println();
		DumpCWL(ListEntry.CWL);
	}
}

int main(int argc, char **argv)
{
	if (argc <= 1)
	{
		std::println("Tests what's in LieAlgRepProd.h");
		std::println("Command-line args:");
		std::println("   prods -- rep products");
		std::println("   powers -- rep powers (plethysms) ");
		std::println("   pwrsyms -- (anti)symmetric rep powers ");
		return 0;
	}

	const bool ProdTest = CheckArgMembership(argc, argv, "prods");
	const bool PowerTest = CheckArgMembership(argc, argv, "powers");
	const bool PwrSymTest = CheckArgMembership(argc, argv, "pwrsyms");
	
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
		
		LieAlgebraParams p0 = {1,2};
		LASnglRepHandler h0(p0);
		LAINT wt01[] = {1,2};
		LAINT wt02[] = {3,1};
		
		std::println("Product of SU(3) {{1,2}} and {{3,1}}");
		CWL = DecomposeRepProduct(h0, wt01, wt02);
		DumpCWL(CWL);
		
		std::println("Quark and Antiquark: Color, Light Flavors, Spin");
		CWL = DecomposeRepProduct(h2, wt2, wtx2);
		DumpCWL(CWL);
		
		println();
	}
	
	if (PowerTest)
	{
		LAYDCWLList CWLList;
		
		std::println("Powers of SU(5) {{1,0,0,0}} decomposed by Young-diagram symmetry type");
		for (int p=0; p<=4; p++)
		{
			std::println("Power = {}", p);
			CWLList = DecomposeRepPower(h1, wt11, p);
			DumpCWLList(CWLList);
		}
		std::println("Powers of SU(5) {{1,0,0,1}} decomposed by Young-diagram symmetry type");
		for (int p=0; p<=4; p++)
		{
			std::println("Power = {}", p);
			CWLList = DecomposeRepPower(h1, wt12, p);
			DumpCWLList(CWLList);
		}
		std::println("Three Quarks");
		CWLList = DecomposeRepPower(h2, wt2,3);
		DumpCWLList(CWLList);
		println();
	}
	
	if (PwrSymTest)
	{
		LACntdMaxWtList CWL;
		
		std::println("Powers of SU(5) {{1,0,0,0}} -- (anti)symmetric");
		for (int p=0; p<=4; p++)
		{
			std::println("Power = {}", p);
			for (int s=1; s>=-1; s-=2)
			{
				std::println("Symm = {}", s);
				CWL = DecomposeRepPwrSym(h1, wt11, p, s);
				DumpCWL(CWL);
			}
		}
		std::println("Powers of SU(5) {{1,0,0,1}} -- (anti)symmetric");
		for (int p=0; p<=4; p++)
		{
			std::println("Power = {}", p);
			for (int s=1; s>=-1; s-=2)
			{
				std::println("Symm = {}", s);
				CWL = DecomposeRepPwrSym(h1, wt12, p, s);
				DumpCWL(CWL);
			}
		}
		std::println("Three Quarks");
		for (int s=1; s>=-1; s-=2)
		{
			std::println("Symm = {}", s);
			CWL = DecomposeRepPwrSym(h2, wt2,3, s);
			DumpCWL(CWL);
		}
		println();
	}
}
