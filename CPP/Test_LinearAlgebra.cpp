/*
	Tests the linear-algebra stuff
*/

#include "Fraction.h"
#include "LinearAlgebra.h"
#include "Test_Shared.h"


template<class N> void DumpScalar(N s)
{
	printf(" %3d\n",s);
}

template<class N, class ND> void MakeVector(vector<N> &v, const ND *dv, int n)
{
	v.clear();
	for (int i=0; i<n; i++)
		v.push_back(dv[i]);
}

template<class N> void DumpVector(vector<N> &v)
{
	for (int i=0; i<v.size(); i++)
		printf(" %3d",v[i]);
	printf("\n");
}

inline const char *TrueFalse(bool b) {return b ? "T" : "F";}

template<class N, class ND> void MakeMatrix(Matrix<N> &M, const ND *dv, int nr, int nc)
{
	M.resize(nr,nc);
	for (int ir=0; ir<nr; ir++)
		for (int ic=0; ic<nc; ic++)
			M(ir,ic) = *(dv++);
}

template<class N> void DumpMatrix(Matrix<N> &M)
{
	for (int ir=0; ir<M.get_rows(); ir++)
	{
		printf("   ");
		for (int ic=0; ic<M.get_cols(); ic++)
			printf(" %3d",M(ir,ic));
	}
	printf("\n");
}


struct HashFunction
{
	char operator() (vector<int> &vec);
};

char HashFunction::operator() (vector<int> &vec)
{
	size_t hashcode = 0;
	for (vector<int>::iterator it = vec.begin(); it != vec.end(); it++)
		hashcode += 13666*hashcode + (*it);
	return hashcode;
}

void VectorSetTest(int rows, int cols, const int *data)
{
	VectorSet<int,HashFunction> VS(cols);
	printf("%4d%4d\n",rows,cols);
	for (int ir=0; ir<rows; ir++)
	{
		pair<bool,size_t> ret = VS.AppendVector(data + ir*cols);
		printf("%4lu\n",ret.second);
	}
	vector<size_t> indxs;
	VS.ExportIndices(indxs);
	printf("%4lu -- ",indxs.size());
	for (vector<size_t>::iterator it=indxs.begin(); it!=indxs.end(); it++)
		printf("%4lu",(*it));
	println();
	Matrix<int> mat;
	VS.ExportMatrix(mat);
	for (int i=0; i<mat.get_rows(); i++)
	{
		for (int j=0; j<mat.get_cols(); j++)
			printf("%4d",mat(i,j));
		printf("\n");
	}
	println();
}


void MatrixInverseTest(int dim, const int *data)
{
	Matrix<int> mat(dim);
	int ix = 0;
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			mat(i,j) = data[ix++];
	for (int i=0; i<dim; i++)
		for (int j=0; j<dim; j++)
			printf("%d %d: %d\n",i,j,mat(i,j));
	Matrix< Fraction<int> > invmat;
	if (MatrixInverse(invmat,mat))
	{
		for (int i=0; i<dim; i++)
			for (int j=0; j<dim; j++)
			{
				Fraction<int> &val = invmat(i,j);
				printf("%d %d: %d/%d\n",i,j,val.get_num(),val.get_den());
			}
		println();
	}
	else
	{
		printf("Singular matrix\n");
	}
}

int main(int argc, char **argv)
{
	if (argc <= 1)
	{
		printf("Tests what's in LinearAlgebra.h\n");
		printf("Command-line args:\n");
		printf("   ops -- elementary operations\n");
		printf("   vset -- vector set: list with uniqueness enforced\n");
		printf("   inverse -- matrix inverse\n");
		return 0;
	}
	const bool OpsTest = CheckArgMembership(argc,argv,"ops");
	const bool VSetTest = CheckArgMembership(argc,argv,"vset");
	const bool InverseTest = CheckArgMembership(argc,argv,"inverse");
	
	if (OpsTest)
	{
		int sres;
		vector<int> vres(4);
		const int dv0[] = {1,3,-2,-4};
		vector<int> v0;
		MakeVector(v0,dv0,4);
		const int dv1[] = {3,4,2,1};
		vector<int> v1;
		MakeVector(v1,dv1,4);
		
		printf("sum\n  -2\n");
		sum(sres,v0);
		DumpScalar(sres);
		
		printf("add_vv\n   4   7   0  -3\n");
		add_vv(vres,v0,v1);
		DumpVector(vres);
		
		printf("sub_vv\n  -2  -1  -4  -5\n");
		sub_vv(vres,v0,v1);
		DumpVector(vres);
		
		printf("addto_vv\n   4   7   0  -3\n");
		copy(v0.begin(),v0.end(),vres.begin());
		addto_vv(vres,v1);
		DumpVector(vres);
		
		printf("subfm_vv\n  -2  -1  -4  -5\n");
		copy(v0.begin(),v0.end(),vres.begin());
		subfm_vv(vres,v1);
		DumpVector(vres);
		
		printf("mul_sv\n  -3  -9   6  12\n");
		mul_sv(vres,-3,v0);
		DumpVector(vres);
		
		printf("mulby_sv\n  -3  -9   6  12\n");
		copy(v0.begin(),v0.end(),vres.begin());
		mulby_sv(vres,-3);
		DumpVector(vres);
		
		printf("mul_vv\n   7\n");
		mul_vv(sres,v0,v1);
		DumpScalar(sres);
		
		printf("muladdto_vsv\n  -8  -9  -8  -7\n");
		copy(v0.begin(),v0.end(),vres.begin());
		muladdto_vsv(vres,-3,v1);
		DumpVector(vres);
		
		printf("Vector comparisons: {1,3,-2,-4} to itself and to {3,4,2,1}\n");
		printf("==  %s %s\n",TrueFalse(VecEqual(v0,v0)),TrueFalse(VecEqual(v0,v1)));
		printf("!=  %s %s\n",TrueFalse(VecNotEqual(v0,v0)),TrueFalse(VecNotEqual(v0,v1)));
		printf(">   %s %s\n",TrueFalse(VecGreaterThan(v0,v0)),TrueFalse(VecGreaterThan(v0,v1)));
		printf("<   %s %s\n",TrueFalse(VecLessThan(v0,v0)),TrueFalse(VecLessThan(v0,v1)));
		printf(">=  %s %s\n",TrueFalse(VecGreaterEqual(v0,v0)),TrueFalse(VecGreaterEqual(v0,v1)));
		printf("<=  %s %s\n",TrueFalse(VecLessEqual(v0,v0)),TrueFalse(VecLessEqual(v0,v1)));
		
		Matrix<int> mres;
		Matrix<int> m0;
		const int dm0[] = {1,2,-3,4, 2,-3,1,5, 3,2,-3,1, 4,4,2,3};
		MakeMatrix(m0,dm0,4,4);
		Matrix<int>m1;
		const int dm1[] = {2,-3,-3, 1,-1,2, 3,-3,-4, 2,-5,4};
		MakeMatrix(m1,dm1,4,3);
		
		printf("mul_sm\n     -3  -6   9 -12     -6   9  -3 -15     -9  -6   9  -3    -12 -12  -6  -9\n");
		mres.resize(4,4);
		mul_sm(mres,-3,m0);
		DumpMatrix(mres);
		
		printf("mulby_sm\n     -3  -6   9 -12     -6   9  -3 -15     -9  -6   9  -3    -12 -12  -6  -9\n");
		mres = m0;
		mulby_sm(mres,-3);
		DumpMatrix(mres);
		
		printf("mul_vm\n -15 -27  -2   5\n");
		mul_vm(vres,v0,m0);
		DumpVector(vres);
		
		printf("mul_mv\n  -3 -29  11   0\n");
		mul_mv(vres,m0,v0);
		DumpVector(vres);
		
		printf("mul_vmv\n -152\n");
		mul_vmv(sres,v0,m0,v1);
		DumpScalar(sres);
		
		printf("mul_mm\n      3 -16  29     14 -31   4      1  -7  11     24 -37   0\n");
		mres.resize(4,3);
		mul_mm(mres,m0,m1);
		DumpMatrix(mres);
	}

	if (VSetTest)
	{
		const int row0 = 6;
		const int col0 = 3;
		const int data0[row0*col0] = {1,3,2, 2,5,4, 1,3,2, 1,-5,-4, 2,5,4, 1,-5,1};
		VectorSetTest(row0,col0,data0);
	}
	
	if (InverseTest)
	{
		const int dim0 = 1;
		const int data0[dim0*dim0] = {1};
		MatrixInverseTest(dim0,data0);
		const int dim1 = 1;
		const int data1[dim1*dim1] = {3};
		MatrixInverseTest(dim1,data1);
		const int dim2 = 2;
		const int data2[dim2*dim2] = {0,1, 2,3};
		MatrixInverseTest(dim2,data2);
		const int dim3 = 3;
		const int data3[dim3*dim3] = {1,1,3, 1,1,2, 3,2,1};
		MatrixInverseTest(dim3,data3);
		const int dim = 4;
		const int data[dim*dim] = {2,-4,-4,-1, -1,-4,0,3, 1,0,-1,1, -2,3,1,-4};
		MatrixInverseTest(dim,data);
	}
}
