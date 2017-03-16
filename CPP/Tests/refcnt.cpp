// Tests a simple reference counter

#include <stdio.h>
#include "../SmartPointer.h"


// Label
int OverallLabel = 0;

// Test class
class Test
{
	int Lbl; // for labeling each instance
	
public:
	Test();
	~Test();
	
	void Show(const char *Prefix); // For showing the class's label
};

Test::Test()
{
	Lbl = OverallLabel++;
	Show("Creating");
}

Test::~Test()
{
	Show("Destroying");
}

void Test::Show(const char *Prefix)
{
	printf("%s: %d\n",Prefix,Lbl);
}

RefCounter<Test> EmitTest()
{
	RefCounter<Test> TRC;
	Test &T = *TRC;
	T.Show("In EmitTest");
	return TRC;
}

int main()
{
	{
		RefCounter<Test> TRC0(NULL);
		Test &T = *TRC0;
		T.Show("In block");
	}

	RefCounter<Test> TRC = EmitTest();
	
	Test &T = *TRC;
	T.Show("In main");
	TRC->Show("In main again");
	
	return 0;
}