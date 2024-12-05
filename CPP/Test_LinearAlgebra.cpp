/*
	Tests the linear-algebra stuff
*/

#include "Fraction.h"
#include "LinearAlgebra.h"
#include "Test_Shared.h"


// The fraction type
using R = Fraction<int>;

// For writing out a fraction
template<>
struct std::formatter<R>
{
	template <typename ParseContext>
	constexpr auto parse(ParseContext &ctxt) {
		return ctxt.begin();
  	}
  	template <typename FormatContext>
	auto format(const R &frac, FormatContext &ctxt) const {
		if (frac.get_den() != 1)
			return std::format_to(ctxt.out(), "{}/{}", frac.get_num(), frac.get_den());
		else
			return std::format_to(ctxt.out(), "{}", frac.get_num());
	}
};


template<typename N>
void DumpScalar(N s) {std::print(" {}", s);}

template<typename Ctnr>
requires ArrayLike<Ctnr>
void DumpVector(Ctnr &v)
{
	for (int i=0; i<v.size(); i++) DumpScalar(v[i]);
	println();
}

// boolean to string: T(rue) or F(alse)
inline std::string B2S(bool x) {return x ? "T" : "F";}

inline const char *TrueFalse(bool b) {return b ? "T" : "F";}

template<typename N, typename ND>
void MakeMatrix(Matrix<N> &M, const ND *dv, int nr, int nc)
{
	M.resize(nr, nc);
	for (int ir=0; ir<nr; ir++)
		for (int ic=0; ic<nc; ic++)
			M(ir, ic) = *(dv++);
}

template<typename N>
void DumpMatrix(Matrix<N> &M)
{
	for (int ir=0; ir<M.get_rows(); ir++)
	{
		std::print("   ");
		for (int ic=0; ic<M.get_cols(); ic++)
			DumpScalar(M(ir, ic));
	}
	println();
}


struct HashFunction
{
	char operator() (const std::vector<int> &vec);
};

char HashFunction::operator() (const std::vector<int> &vec)
{
	size_t hashcode = 0;
	for (auto entry: vec)
		hashcode += 13666*hashcode + entry;
	return hashcode;
}

void VectorSetTest(int rows, int cols, const int *data)
{
	std::println("Input size: {} {}", rows, cols);
	std::println("Data:");
	for (int ir=0; ir<rows; ir++) {
		for (int ic=0; ic<cols; ic++)
			DumpScalar(data[ir*cols+ic]);
		println();
	}
			
	std::println("Index of data vector in list");
	
	VectorSet<int, HashFunction> VS(cols);
	for (int ir=0; ir<rows; ir++)
	{
		std::pair<bool, size_t> ret = VS.AppendVector(data + ir*cols);
		std::println("{:4}", ret.second);
	}
	
	std::vector<size_t> indxs;
	VS.ExportIndices(indxs);
	std::println("Number of entries: {}", indxs.size());
	
	std::println("Index values:");
	for (auto entry: indxs) DumpScalar(entry);
	println();
	
	std::println("The added-in vectors");
	Matrix<int> mat;
	VS.ExportMatrix(mat);
	for (int i=0; i<mat.get_rows(); i++)
	{
		for (int j=0; j<mat.get_cols(); j++)
			DumpScalar(mat(i, j));
		println();
	}
	println();
}


void MatrixInverseTest(int dim, const int *data)
{
	Matrix<int> mat(dim);
	MakeMatrix(mat, data, dim, dim);
	DumpMatrix(mat);
	
	Matrix<R> invmat;
	if (MatrixInverse(invmat, mat)) {DumpMatrix(invmat);}
	else {std::print("Singular matrix\n");}
}


int main(int argc, char **argv)
{
	if (argc <= 1)
	{
		std::println("Tests what's in LinearAlgebra.h");
		std::println("Command-line args:");
		std::println("   ops -- elementary operations");
		std::println("   vset -- vector set: list with uniqueness enforced");
		std::println("   inverse -- matrix inverse");
		return 0;
	}
	const bool OpsTest = CheckArgMembership(argc, argv, "ops");
	const bool VSetTest = CheckArgMembership(argc, argv, "vset");
	const bool InverseTest = CheckArgMembership(argc, argv, "inverse");
	
	if (OpsTest)
	{
		int sres;
		std::vector<int> vres(4);
		const std::vector<int> v0 = {1,3,-2,-4};
		const std::vector<int> v1 = {3,4,2,1};
		
		std::println("sum");
		DumpScalar(-2); println();
		sum(sres, v0);
		DumpScalar(sres); println();
		
		std::println("add_vv");
		const std::vector<int> vcadd = {4,7,0,-3};
		DumpVector(vcadd);
		add_vv(vres, v0, v1);
		DumpVector(vres);
		
		std::println("sub_vv");
		const std::vector<int> vcsub = {-2,-1,-4,-5};
		DumpVector(vcsub);
		sub_vv(vres, v0, v1);
		DumpVector(vres);
		
		std::println("addto_vv");
		const std::vector<int> vcaddto = {4,7,0,-3};
		DumpVector(vcaddto);
		copy(v0.begin(), v0.end(), vres.begin());
		addto_vv(vres, v1);
		DumpVector(vres);
		
		std::println("subfm_vv");
		const std::vector<int> vcsubfm = {-2,-1,-4,-5};
		DumpVector(vcsubfm);
		copy(v0.begin(), v0.end(), vres.begin());
		subfm_vv(vres, v1);
		DumpVector(vres);
		
		std::println("mul_sv");
		const std::vector<int> vcmul = {-3,-9,6,12};
		DumpVector(vcmul);
		mul_sv(vres,-3, v0);
		DumpVector(vres);
		
		std::println("mulby_sv");
		const std::vector<int> vcmulby = {-3,-9,6,12};
		DumpVector(vcmulby);
		copy(v0.begin(), v0.end(), vres.begin());
		mulby_sv(vres,-3);
		DumpVector(vres);
		
		std::println("mul_vv");
		DumpScalar(7); println();
		mul_vv(sres, v0, v1);
		DumpScalar(sres); println();
		
		std::println("muladdto_vsv");
		const std::vector<int> vcmladto = {-8,-9,-8,-7};
		DumpVector(vcmladto);
		copy(v0.begin(), v0.end(), vres.begin());
		muladdto_vsv(vres,-3, v1);
		DumpVector(vres);
		
		std::println("Vector comparisons: {{1,3,-2,-4}} to itself and to {{3,4,2,1}}");
		std::println("==  {} {}", B2S(VecEqual(v0, v0)), B2S(VecEqual(v0, v1)));
		std::println("!=  {} {}", B2S(VecNotEqual(v0, v0)), B2S(VecNotEqual(v0, v1)));
		std::println(">   {} {}", B2S(VecGreaterThan(v0, v0)), B2S(VecGreaterThan(v0, v1)));
		std::println("<   {} {}", B2S(VecLessThan(v0, v0)), B2S(VecLessThan(v0, v1)));
		std::println(">=  {} {}", B2S(VecGreaterEqual(v0, v0)), B2S(VecGreaterEqual(v0, v1)));
		std::println("<=  {} {}", B2S(VecLessEqual(v0, v0)), B2S(VecLessEqual(v0, v1)));
		
		Matrix<int> mres;
		
		Matrix<int> m0;
		const int dm0[] = {1,2,-3,4, 2,-3,1,5, 3,2,-3,1, 4,4,2,3};
		MakeMatrix(m0, dm0,4,4);
		Matrix<int>m1;
		const int dm1[] = {2,-3,-3, 1,-1,2, 3,-3,-4, 2,-5,4};
		MakeMatrix(m1, dm1,4,3);
		
		std::println("mul_sm");
		Matrix<int> mrm;
		const int dmrm[] = {-3,-6,9,-12, -6,9,-3,-15, -9,-6,9,-3, -12,-12,-6,-9};
		MakeMatrix(mrm, dmrm,4,4);
		DumpMatrix(mrm);
		mres.resize(4,4);
		mul_sm(mres,-3, m0);
		DumpMatrix(mres);

		std::println("mulby_sm");
		Matrix<int> mrmb;
		const int dmrmb[] = {-3,-6,9,-12, -6,9,-3,-15, -9,-6,9,-3, -12,-12,-6,-9};
		MakeMatrix(mrmb, dmrmb,4,4);
		DumpMatrix(mrmb);
		mres = m0;
		mulby_sm(mres,-3);
		DumpMatrix(mres);
		
		std::println("mul_vm");
		const std::vector<int> vrvm = {-15,-27,-2,5};
		DumpVector(vrvm);
		mul_vm(vres, v0, m0);
		DumpVector(vres);
		
		std::println("mul_mv");
		const std::vector<int> vrmv = {-3,-29,11,0};
		DumpVector(vrmv);
		mul_mv(vres, m0, v0);
		DumpVector(vres);
		
		std::println("mul_vmv");
		DumpScalar(-152); println();
		mul_vmv(sres, v0, m0, v1);
		DumpScalar(sres); println();
		
		std::println("mul_mm");
		Matrix<int> mrmm;
		const int dmrmm[] = {3,-16,29, 14,-31,4, 1,-7,11, 24,-37,0};
		MakeMatrix(mrmm, dmrmm,4,3);
		DumpMatrix(mrmm);
		mres.resize(4,3);
		mul_mm(mres, m0, m1);
		DumpMatrix(mres);
	}

	if (VSetTest)
	{
		const int row0 = 6;
		const int col0 = 3;
		const int data0[row0*col0] = {1,3,2, 2,5,4, 1,3,2, 1,-5,-4, 2,5,4, 1,-5,1};
		VectorSetTest(row0, col0, data0);
	}
	
	if (InverseTest)
	{
		const int dim0 = 1;
		const int data0[dim0*dim0] = {1};
		MatrixInverseTest(dim0, data0);
		const int dim1 = 1;
		const int data1[dim1*dim1] = {3};
		MatrixInverseTest(dim1, data1);
		const int dim2 = 2;
		const int data2[dim2*dim2] = {0,1, 2,3};
		MatrixInverseTest(dim2, data2);
		const int dim3 = 3;
		const int data3[dim3*dim3] = {1,1,3, 1,1,2, 3,2,1};
		MatrixInverseTest(dim3, data3);
		const int dim = 4;
		const int data[dim*dim] = {2,-4,-4,-1, -1,-4,0,3, 1,0,-1,1, -2,3,1,-4};
		MatrixInverseTest(dim, data);
	}
}
