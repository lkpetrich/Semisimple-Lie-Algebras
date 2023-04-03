/*
	Symmetrized Angular-Momentum Invariants
	
	Calculating with C++ for speed
	
	A command-line program
	
	Args:
	Maximum value of 2j: nonnegative integer
	Power: positive integer
	Is symmetric: 0 or 1
	Method: 0, 1, 2, 3, 4
	
	Prints out results for twice angular momentum in integer steps from 0
	For all but the last method,
		the power is the power of matrices
	M(i,j) = e(i,j,k)*L(k)
	composed from the antisymmetric symbol e and generators L:
	M(i,j)*M(j,k)*...*M(l,i)
	
	The different choices are different stages in code development,
	included here to illustrate run times.
	
	The last method makes a product of 2*power operators L:
	(L.L)*(L.L)*...*(L.L)
	In effect, powers of (L.L)
	
	The generators themselves are implemented as matrices
	for QM angular-momentum states m = -j, -j+1, ..., j-1, j
	
	The output is, for each value of 2*j
	The invariant value
	The difference between the sum of products of generators, and
	(invariant value) * (identity matrix)
	
	The invariant values are afterward collected.
	
	The calculation is entirely numerical.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>


// Doesn't have any vetting of inputs

// For iterating through index values:
// (len) of them, each with range (base)

class Digits {

	int base, len;
	std::vector<int> digits;

public:
	// Digit base, then number of digits
	Digits(int _base, int _len): base(_base), len(_len), digits(_len) {Reset();}
	
	void Reset() {std::fill(digits.begin(),digits.end(),0);}
	
	auto GetBase() const {return base;}
	auto GetLen() const {return len;}
	
	auto operator[](int ix) const {return digits[ix];}
	
	// Returns true if advanced, false if rolled over
	bool Next() {
		int ix = 0;
		while(ix < len) {
			if ((++digits[ix]) >= base) {
				digits[ix++] = 0;
			} else {
				return true;
			}
		}
		return false;
	}
};


void Dump(const Digits &dgts) {
	for (int i=0; i<dgts.GetLen(); i++)
		printf(" %d",dgts[i]);
	printf("\n");
}


// For iterating through permutations: 0, 1, ..., (len-1)

class Permute {

	int len;
	std::vector<int> pmt;

public:
	Permute(int _len): len(_len), pmt(_len) {Reset();}
	
	void Reset() {for(int ix=0; ix<len; ix++) pmt[ix]=ix;}
	
	auto GetLen() const {return len;}
	
	auto operator[](int ix) const {return pmt[ix];}
	
	// Returns true if advanced, false if rolled over
	bool Next() {return next_permutation(pmt.begin(), pmt.end());}
};


void Dump(const Permute &pmt) {
	for (int i=0; i<pmt.GetLen(); i++)
		printf(" %d",pmt[i]);
	printf("\n");
}

// Like the above, but does pairs of the same value

class DoublePermute {

	int len;
	std::vector<int> pmt;

public:
	DoublePermute(int _len): len(_len), pmt(2*_len) {Reset();}
	
	void Reset() {for(int ix=0; ix<len; ix++) {pmt[2*ix]=ix; pmt[2*ix+1]=ix;}}
	
	auto GetLen() const {return len;}
	
	auto GetExtLen() const {return 2*len;}
	
	auto operator[](int ix) const {return pmt[ix];}
	
	// Returns true if advanced, false if rolled over
	bool Next() {return next_permutation(pmt.begin(), pmt.end());}
};


// Multidimensional arrays: 2, 3, 4 dimensions

template<typename T> class Array2 {

	int dim1, dim2;
	std::vector<T> data;
	
public:
	Array2(int dim1_, int dim2_):
		dim1(dim1_), dim2(dim2_),
		data(dim1*dim2) {}
	
	Array2() {} // Empty constructor
	
	auto GetDim1() const {return dim1;}
	auto GetDim2() const {return dim2;}
	
	void fill(const T &value) {std::fill(data.begin(),data.end(),value);}
	
	auto &at(int ix1, int ix2)
		{return data[dim2*ix1+ix2];}
	
	void CopyIn(const Array2<T> &arr)
		{std::copy(arr.data.begin(), arr.data.end(), data.begin());}
};

template<typename T> class Array3 {

	int dim1, dim2, dim3;
	std::vector<T> data;
	
public:
	Array3(int dim1_, int dim2_, int dim3_):
		dim1(dim1_), dim2(dim2_), dim3(dim3_),
		data(dim1*dim2*dim3) {}
		
	Array3() {} // Empty constructor

	auto GetDim1() const {return dim1;}
	auto GetDim2() const {return dim2;}
	auto GetDim3() const {return dim3;}
	
	void fill(const T &value) {std::fill(data.begin(),data.end(),value);}
	
	auto &at(int ix1, int ix2, int ix3)
		{return data[dim3*(dim2*ix1+ix2)+ix3];} 
	
	void CopyIn(const Array3<T> &arr)
		{std::copy(arr.data.begin(), arr.data.end(), data.begin());}

};

template<typename T> class Array4 {

	int dim1, dim2, dim3, dim4;
	std::vector<T> data;
	
public:
	Array4(int dim1_, int dim2_, int dim3_, int dim4_):
		dim1(dim1_), dim2(dim2_), dim3(dim3_), dim4(dim4_),
		data(dim1*dim2*dim3*dim4) {}
		
	Array4() {} // Empty constructor

	auto GetDim1() const {return dim1;}
	auto GetDim2() const {return dim2;}
	auto GetDim3() const {return dim3;}
	auto GetDim4() const {return dim4;}
	
	void fill(const T &value) {std::fill(data.begin(),data.end(),value);}
	
	auto &at(int ix1, int ix2, int ix3, int ix4)
		{return data[dim4*(dim3*(dim2*ix1+ix2)+ix3)+ix4];} 
	
	void CopyIn(const Array4<T> &arr)
		{std::copy(arr.data.begin(), arr.data.end(), data.begin());}
};

void Dump(Array2<double> &arr) {
	for (int i1 =0; i1 < arr.GetDim1(); i1++) {
		for (int i2=0; i2 < arr.GetDim2(); i2++) {
			printf(" %lg",arr.at(i1,i2));
		}
		printf("\n");
	}
}

void Dump(Array3<double> &arr) {
	for (int i1 =0; i1 < arr.GetDim1(); i1++) {
		for (int i2=0; i2 < arr.GetDim2(); i2++) {
			for (int i3=0; i3<arr.GetDim3(); i3++) {
				printf(" %lg",arr.at(i1,i2,i3));
			}
			printf("\n");
		}
		printf("\n");
	}
}

void Dump(Array4<double> &arr) {
	for (int i1 =0; i1 < arr.GetDim1(); i1++) {
		for (int i2=0; i2 < arr.GetDim2(); i2++) {
			for (int i3=0; i3<arr.GetDim3(); i3++) {
				for (int i4=0; i4<arr.GetDim4(); i4++) {
					printf(" %lg",arr.at(i1,i2,i3,i4));
				}
				printf("\n");
			}
			printf("\n");
		}
		printf("---\n\n");
	}
}

// Moved out here if I want it repeated with different parameters
// Needs 2*j, power, whether symmetrized

// Operator matrix - each component
// Zero multiplied by anything is zero, so using that as a shortcut


// For recognizing when to skip over matrix multiplications
template<typename T> bool IsMatrixNonzero(Array2<T> &mat) {
	for (int i1=0; i1<mat.GetDim1(); i1++)
		for (int i2=0; i2<mat.GetDim2(); i2++)
			if (mat.at(i1,i2) != 0)
				return true;
	return false;
}

// Operator matrix - each component
// Zero multiplied by anything is zero, so using that as a shortcut
struct OperMatrix {
	bool IsNonzero;
	Array2<double> data;
	
	OperMatrix(int Dim1, int Dim2): IsNonzero(false), data(Dim1,Dim2) {data.fill(0);}
	
	OperMatrix(): IsNonzero(false) {} // Empty constructor
	
	bool CheckNonzero();
};

bool OperMatrix::CheckNonzero() {
	IsNonzero = false;
	for (int i1=0; i1<data.GetDim1(); i1++)
		for (int i2=0; i2<data.GetDim2(); i2++) {
			if (data.at(i1,i2) != 0) {
				IsNonzero = true;
				goto Escape;
			}
		}
	Escape:
	return IsNonzero;
}

// Using sparse arrays for the algebra operators to speed them up
template<int N, typename T> struct SparseArrayEntry {int ixs[N]; T val;};
template<int N, typename T> using SparseArray = std::vector< SparseArrayEntry<N,T> >;


// Dense data matrix
using DataMat = Array2<double>;

// Use sparse arrays for the operators:
using OperEnt = SparseArrayEntry<2,double>;
using OperMat = SparseArray<2,double>;


void AddTo(DataMat &SumMat, DataMat &Mat, int dim) {
	for (int i1=0; i1<dim; i1++)
		for (int i2=0; i2<dim; i2++)
			SumMat.at(i1,i2) += Mat.at(i1,i2);
}

void ScaleByNumPermute(DataMat &Mat, int NumPermute) {
	if (NumPermute != 1) {
		// Divide by the number of permutations
		double nprecip = 1.0/NumPermute;
		
		for (int i1=0; i1<Mat.GetDim1(); i1++)
			for (int i2=0; i2<Mat.GetDim2(); i2++)
				Mat.at(i1,i2) *= nprecip;
	}
}

double FindDiagonal(DataMat &Mat, int dim) {

	
	double DiagVal = 0;
	for (int i=0; i<dim; i++)
		DiagVal += Mat.at(i,i);
	DiagVal /= dim;
	
	for (int i=0; i<dim; i++)
		Mat.at(i,i) -= DiagVal;
	
	Dump(Mat);
		
	return DiagVal;
}

double FindInvariant0(int twoj, int power, bool is_symm) {

	// The operator vectors
	Array3<double> OpVec(3,twoj+1,twoj+1);
	OpVec.fill(0);
	
	for (int k=0; k<twoj; k++) {
		OpVec.at(0, k+1, k) = OpVec.at(1, k, k+1) = sqrt(1.*(k+1)*(twoj-k));
	}
	for (int k=0; k<=twoj; k++) {
		OpVec.at(2, k, k) = - 0.5*twoj + k;
	}
	
	Array3<double> Commutator(3,3,3);
	Commutator.fill(0);
	Commutator.at(0,0,2) = - 1;
	Commutator.at(0,2,1) = 0.5;
	Commutator.at(1,1,2) = 1;
	Commutator.at(1,2,0) = - 0.5;
	Commutator.at(2,0,0) = 1;
	Commutator.at(2,1,1) = - 1;
	
	Array4<double> OpMat(3,3,twoj+1,twoj+1);
	for (int i1=0; i1<3; i1++)
		for (int i2=0; i2<3; i2++)
			for (int i3=0; i3<=twoj; i3++)
				for (int i4=0; i4<=twoj; i4++) {
					double sum = 0;
					for (int j=0; j<3; j++)
						sum += Commutator.at(i1,i2,j) * OpVec.at(j,i3,i4);
					OpMat.at(i1,i2,i3,i4) = sum;
				}
	
	DataMat ProdMat(twoj+1,twoj+1), NewProdMat(twoj+1,twoj+1);
	DataMat ItmdSumProdMat(twoj+1,twoj+1), SumProdMat(twoj+1,twoj+1);
	
	int NumPermute = 0;
	SumProdMat.fill(0);
	Digits OpIndices(3,power);
	Permute OpLocs(power);

	do {
		NumPermute++;
		ItmdSumProdMat.fill(0);
		do {
			for (int iop = 0; iop < power; iop++) {
				// Locations in the operator matrices
				auto ixp = OpLocs[iop];
				auto ixpnxt = (ixp + 1) % power;
				auto ix = OpIndices[ixp];
				auto ixnxt = OpIndices[ixpnxt];
				
				if (iop > 0) {
					// Multiply by operator
					for (int i1=0; i1<=twoj; i1++)
						for (int i2=0; i2<=twoj; i2++) {
							double sum = 0;
							for (int j=0; j<=twoj; j++)
								sum += ProdMat.at(i1,j) * OpMat.at(ix,ixnxt,j,i2);
							NewProdMat.at(i1,i2) = sum;
						}
					
					// Copy in the new product
					for (int i1=0; i1<=twoj; i1++)
						for (int i2=0; i2<=twoj; i2++)
							ProdMat.at(i1,i2) = NewProdMat.at(i1,i2);
				} else {
					// Copy in the operator
					for (int i1=0; i1<=twoj; i1++)
						for (int i2=0; i2<=twoj; i2++)
							ProdMat.at(i1,i2) = OpMat.at(ix,ixnxt,i1,i2);
				}
			}
			// Add in the operator product
			AddTo(ItmdSumProdMat, ProdMat, twoj+1);
		
		} while (OpIndices.Next());
		// Add in the intermediate sum
		AddTo(SumProdMat, ItmdSumProdMat, twoj+1);
		
		if (!is_symm) break;
		
	} while (OpLocs.Next());
	
	ScaleByNumPermute(SumProdMat, NumPermute);
	return FindDiagonal(SumProdMat, twoj+1);
}


// With vectors and matrices of operator matrix blocks

double FindInvariant1(int twoj, int power, bool is_symm) {

	// The operator vectors
	std::vector<OperMatrix> OpVec(3);
	for (int i=0; i<3; i++)
		OpVec[i] = OperMatrix(twoj+1,twoj+1);
	
	for (int k=0; k<twoj; k++) {
		OpVec[0].data.at(k+1, k) = OpVec[1].data.at(k, k+1) = sqrt(1.*(k+1)*(twoj-k));
	}
	for (int k=0; k<=twoj; k++) {
		OpVec[2].data.at(k, k) = - 0.5*twoj + k;
	}
	
	Array3<double> Commutator(3,3,3);
	Commutator.fill(0);
	Commutator.at(0,0,2) = - 1;
	Commutator.at(0,2,1) = 0.5;
	Commutator.at(1,1,2) = 1;
	Commutator.at(1,2,0) = - 0.5;
	Commutator.at(2,0,0) = 1;
	Commutator.at(2,1,1) = - 1;
	
	// The operator matrices
	Array2<OperMatrix> OpMat(3,3);
	for (int i1=0; i1<3; i1++)
		for (int i2=0; i2<3; i2++) {
			Array2<double> OpData(twoj+1,twoj+1);
			OpData.fill(0);
			OpMat.at(i1,i2).data = OpData;
		}
	
	for (int i1=0; i1<3; i1++)
		for (int i2=0; i2<3; i2++) {
			OperMatrix &OpCpnt = OpMat.at(i1,i2);
			for (int i3=0; i3<=twoj; i3++)
				for (int i4=0; i4<=twoj; i4++) {
					double sum = 0;
					for (int j=0; j<3; j++)
						sum += Commutator.at(i1,i2,j) * OpVec[j].data.at(i3,i4);
					OpCpnt.data.at(i3,i4) = sum;
				}
			OpCpnt.CheckNonzero();
		}
	
	OperMatrix ProdMat(twoj+1,twoj+1);
	DataMat NewProdMat(twoj+1,twoj+1);
	DataMat ItmdSumProdMat(twoj+1,twoj+1), SumProdMat(twoj+1,twoj+1);
	
	int NumPermute = 0;
	SumProdMat.fill(0);
	Digits OpIndices(3,power);
	Permute OpLocs(power);

	do {
		NumPermute++;
		ItmdSumProdMat.fill(0);
		do {
			for (int iop = 0; iop < power; iop++) {
				// Locations in the operator matrices
				auto ixp = OpLocs[iop];
				auto ixpnxt = (ixp + 1) % power;
				auto ix = OpIndices[ixp];
				auto ixnxt = OpIndices[ixpnxt];
				OperMatrix &OpCpnt = OpMat.at(ix,ixnxt);
				
				if (iop > 0) {
					// Multiply by operator
					for (int i1=0; i1<=twoj; i1++)
						for (int i2=0; i2<=twoj; i2++) {
							double sum = 0;
							for (int j=0; j<=twoj; j++)
								sum += ProdMat.data.at(i1,j) * OpCpnt.data.at(j,i2);
							NewProdMat.at(i1,i2) = sum;
						}
						
					// Copy in the new product
					ProdMat.data.CopyIn(NewProdMat);
					ProdMat.CheckNonzero();
				} else {
					// Copy in the operator
					ProdMat.IsNonzero = true;
					ProdMat.data.CopyIn(OpCpnt.data);
				}
			}
			// Add in the operator product
			AddTo(ItmdSumProdMat, ProdMat.data, twoj+1);
		
		} while (OpIndices.Next());
		// Add in the intermediate sum
		AddTo(SumProdMat, ItmdSumProdMat, twoj+1);
		
		if (!is_symm) break;
		
	} while (OpLocs.Next());
	
	ScaleByNumPermute(SumProdMat, NumPermute);
	return FindDiagonal(SumProdMat, twoj+1);
}

// Skip over multiplication of all-zero matrices

double FindInvariant2(int twoj, int power, bool is_symm) {

	// The operator vectors
	std::vector<OperMatrix> OpVec(3);
	for (int i=0; i<3; i++)
		OpVec[i] = OperMatrix(twoj+1,twoj+1);
	
	for (int k=0; k<twoj; k++) {
		OpVec[0].data.at(k+1, k) = OpVec[1].data.at(k, k+1) = sqrt(1.*(k+1)*(twoj-k));
	}
	for (int k=0; k<=twoj; k++) {
		OpVec[2].data.at(k, k) = - 0.5*twoj + k;
	}
	
	Array3<double> Commutator(3,3,3);
	Commutator.fill(0);
	Commutator.at(0,0,2) = - 1;
	Commutator.at(0,2,1) = 0.5;
	Commutator.at(1,1,2) = 1;
	Commutator.at(1,2,0) = - 0.5;
	Commutator.at(2,0,0) = 1;
	Commutator.at(2,1,1) = - 1;
	
	// The operator matrices
	Array2<OperMatrix> OpMat(3,3);
	for (int i1=0; i1<3; i1++)
		for (int i2=0; i2<3; i2++) {
			Array2<double> OpData(twoj+1,twoj+1);
			OpData.fill(0);
			OpMat.at(i1,i2).data = OpData;
		}
	
	for (int i1=0; i1<3; i1++)
		for (int i2=0; i2<3; i2++) {
			OperMatrix &OpCpnt = OpMat.at(i1,i2);
			for (int i3=0; i3<=twoj; i3++)
				for (int i4=0; i4<=twoj; i4++) {
					double sum = 0;
					for (int j=0; j<3; j++)
						sum += Commutator.at(i1,i2,j) * OpVec[j].data.at(i3,i4);
					OpCpnt.data.at(i3,i4) = sum;
				}
			OpCpnt.CheckNonzero();
		}
	
	OperMatrix ProdMat(twoj+1,twoj+1);
	DataMat NewProdMat(twoj+1,twoj+1);
	DataMat ItmdSumProdMat(twoj+1,twoj+1), SumProdMat(twoj+1,twoj+1);
	
	int NumPermute = 0;
	SumProdMat.fill(0);
	Digits OpIndices(3,power);
	Permute OpLocs(power);

	do {
		NumPermute++;
		ItmdSumProdMat.fill(0);
		
		do {
			
			// Check on whether any of the operators is zero.
			bool AllNonzero = true;
			
			for (int iop = 0; iop < power; iop++) {
				// Locations in the operator matrices
				auto ixp = OpLocs[iop];
				auto ixpnxt = (ixp + 1) % power;
				auto ix = OpIndices[ixp];
				auto ixnxt = OpIndices[ixpnxt];
				OperMatrix &OpCpnt = OpMat.at(ix,ixnxt);
				
				if (!OpCpnt.IsNonzero) {
					AllNonzero = false;
					break;
				}
			}
			
			if (!AllNonzero) continue;
			
			// All nonzero; continue
			for (int iop = 0; iop < power; iop++) {
				// Locations in the operator matrices
				auto ixp = OpLocs[iop];
				auto ixpnxt = (ixp + 1) % power;
				auto ix = OpIndices[ixp];
				auto ixnxt = OpIndices[ixpnxt];
				OperMatrix &OpCpnt = OpMat.at(ix,ixnxt);
				
				if (iop > 0) {
					// Multiply by operator
					for (int i1=0; i1<=twoj; i1++)
						for (int i2=0; i2<=twoj; i2++) {
							double sum = 0;
							for (int j=0; j<=twoj; j++)
								sum += ProdMat.data.at(i1,j) * OpCpnt.data.at(j,i2);
							NewProdMat.at(i1,i2) = sum;
						}
						
					// Copy in the new product
					ProdMat.data.CopyIn(NewProdMat);
					// Non need to continue if one gets zero
					ProdMat.CheckNonzero();
					if (!ProdMat.IsNonzero) break;
				} else {
					// Copy in the operator
					ProdMat.IsNonzero = true;
					ProdMat.data.CopyIn(OpCpnt.data);
				}
			}
			// Add in the operator product
			AddTo(ItmdSumProdMat, ProdMat.data, twoj+1);
		
		} while (OpIndices.Next());
		// Add in the intermediate sum
		AddTo(SumProdMat, ItmdSumProdMat, twoj+1);
		
		if (!is_symm) break;
		
	} while (OpLocs.Next());
	
	ScaleByNumPermute(SumProdMat, NumPermute);
	return FindDiagonal(SumProdMat, twoj+1);
}

// Use sparse arrays for the operators

double FindInvariant3(int twoj, int power, bool is_symm) {

	// The operator vectors
	std::vector<OperMat> OpVec(3);
	for (int k=0; k<twoj; k++) {
		double val = sqrt(1.*(k+1)*(twoj-k));
		OperEnt ent0{k+1,k,val};
		OperEnt ent1{k,k+1,val};
		OpVec[0].push_back(ent0);
		OpVec[1].push_back(ent1);
	}
	for (int k=0; k<=twoj; k++) {
		double val = - 0.5*twoj + k;
		OperEnt ent2{k,k,val};
		OpVec[2].push_back(ent2);
	}
	
	Array3<double> Commutator(3,3,3);
	Commutator.fill(0);
	Commutator.at(0,0,2) = - 1;
	Commutator.at(0,2,1) = 0.5;
	Commutator.at(1,1,2) = 1;
	Commutator.at(1,2,0) = - 0.5;
	Commutator.at(2,0,0) = 1;
	Commutator.at(2,1,1) = - 1;
	
	// The operator matrices
	Array2<OperMat> OpMat(3,3);
	DataMat OMTemp(twoj+1,twoj+1);
	for (int i1=0; i1<3; i1++)
		for (int i2=0; i2<3; i2++) {
			OMTemp.fill(0);
			auto &OpCpnt = OpMat.at(i1,i2);
			for (int j=0; j<3; j++) {
				auto cmt = Commutator.at(i1,i2,j);
				for (auto ent: OpVec[j]) {
					OMTemp.at(ent.ixs[0],ent.ixs[1]) += cmt * ent.val;
				}
			}
			auto &OME = OpMat.at(i1,i2);
			for (int i3=0; i3<=twoj; i3++)
				for (int i4=0; i4<=twoj; i4++) {
					double val = OMTemp.at(i3,i4);
					if (val != 0) {
						OperEnt ent{i3,i4,val};
						OME.push_back(ent);
					}
				}
			}
	
	DataMat ProdMat(twoj+1,twoj+1), NewProdMat(twoj+1,twoj+1);
	DataMat ItmdSumProdMat(twoj+1,twoj+1), SumProdMat(twoj+1,twoj+1);
	
	int NumPermute = 0;
	SumProdMat.fill(0);
	Digits OpIndices(3,power);
	Permute OpLocs(power);

	do {
		NumPermute++;
		ItmdSumProdMat.fill(0);
		
		do {
			// Check on whether any of the operators is zero.
			// If it is, then then product will also be zero,
			// and there is no need to waste CPU cycles
			// on that calculation.
			bool AllNonzero = true;
			
			for (int iop = 0; iop < power; iop++) {
				// Locations in the operator matrices
				auto ixp = OpLocs[iop];
				auto ixpnxt = (ixp + 1) % power;
				auto ix = OpIndices[ixp];
				auto ixnxt = OpIndices[ixpnxt];
				auto &OpCpnt = OpMat.at(ix,ixnxt);
				
				if (OpCpnt.empty()) {
					AllNonzero = false;
					break;
				}
			}
			
			if (!AllNonzero) continue;
			
			for (int iop = 0; iop < power; iop++) {
				// Locations in the operator matrices
				auto ixp = OpLocs[iop];
				auto ixpnxt = (ixp + 1) % power;
				auto ix = OpIndices[ixp];
				auto ixnxt = OpIndices[ixpnxt];
				auto &OpCpnt = OpMat.at(ix,ixnxt);
				
				if (iop > 0) {
					// Multiply by operator
					NewProdMat.fill(0);
					for (int i=0; i<=twoj; i++)
						for (auto ent: OpCpnt) {
							NewProdMat.at(i,ent.ixs[1]) +=
								ProdMat.at(i,ent.ixs[0]) * ent.val;
						}
					
					// Copy in the new product
					ProdMat.CopyIn(NewProdMat);
					// No need to continue if one gets zero
					if (!IsMatrixNonzero(ProdMat)) break;
				} else {
					// Copy in the operator
					ProdMat.fill(0);
					for (auto ent: OpCpnt) {
						ProdMat.at(ent.ixs[0],ent.ixs[1]) = ent.val;
					}					
				}
			}
			// Add in the operator product
			AddTo(ItmdSumProdMat, ProdMat, twoj+1);
		
		} while (OpIndices.Next());
		// Add in the intermediate sum
		AddTo(SumProdMat, ItmdSumProdMat, twoj+1);
		
		if (!is_symm) break;
		
	} while (OpLocs.Next());
	
	ScaleByNumPermute(SumProdMat, NumPermute);
	return FindDiagonal(SumProdMat, twoj+1);
}

// Use vectors instead of matrices

double FindInvariant4(int twoj, int power, bool is_symm) {
	
	// The operator vectors
	// First the plain operator,
	// then the operator multiplied by the inverse of the algebra metric
	// Do that to avoid using complex numbers
	// 4 times the calculation for the same results
	std::vector<OperMat> OpVec(6);
	for (int k=0; k<twoj; k++) {
		double val = sqrt(1.*(k+1)*(twoj-k));
		OperEnt ent0{k+1,k,val};
		OperEnt ent1{k,k+1,val};
		OpVec[0].push_back(ent0);
		OpVec[1].push_back(ent1);
		ent0.val *= 0.5;
		ent1.val *= 0.5;
		OpVec[3].push_back(ent1);
		OpVec[4].push_back(ent0);
		/*
		double val = sqrt(1.*(k+1)*(twoj-k))/sqrt(2.);
		OperEnt ent0{k+1,k,val};
		OperEnt ent1{k,k+1,val};
		OpVec[0].push_back(ent0);
		OpVec[0].push_back(ent1);
		ent1.val *= -1;
		OpVec[1].push_back(ent0);
		OpVec[1].push_back(ent1);
		ent1.val *= -1;
		OpVec[3].push_back(ent0);
		OpVec[3].push_back(ent1);
		ent0.val *= -1;
		OpVec[4].push_back(ent0);
		OpVec[4].push_back(ent1);
		
		double val = sqrt(1.*(k+1)*(twoj-k));
		OperEnt ent0{k+1,k,val};
		OperEnt ent1{k,k+1,val};
		OpVec[0].push_back(ent0);
		OpVec[0].push_back(ent1);
		ent1.val *= -1;
		OpVec[1].push_back(ent0);
		OpVec[1].push_back(ent1);
		ent0.val *= 0.5;
		ent1.val *= 0.5;
		ent1.val *= -1;
		OpVec[3].push_back(ent0);
		OpVec[3].push_back(ent1);
		ent0.val *= -1;
		OpVec[4].push_back(ent0);
		OpVec[4].push_back(ent1);
		*/
	}
	for (int k=0; k<=twoj; k++) {
		double val = - 0.5*twoj + k;
		OperEnt ent2{k,k,val};
		OpVec[2].push_back(ent2);
		OpVec[5].push_back(ent2);
	}
	
	DataMat ProdMat(twoj+1,twoj+1), NewProdMat(twoj+1,twoj+1);
	DataMat ItmdSumProdMat(twoj+1,twoj+1), SumProdMat(twoj+1,twoj+1);
	
	int NumPermute = 0;
	SumProdMat.fill(0);
	Digits OpIndices(3,power);
	std::vector<int> OpIndExtend(2*power);
	DoublePermute OpLocs(power);
	std::vector<int> OpLocExtend(2*power);
	std::vector<int> OpLocSecond(power);
	
	do {
		bool InOrder = true;
		int CurrentOne = -1;
		for (int i=0; i<2*power; i++)
			if (OpLocs[i] > CurrentOne) {
				if (OpLocs[i] > CurrentOne + 1) {
					InOrder = false;
					break;
				} else {
					CurrentOne = OpLocs[i];
				}
			}
		
		if (!InOrder) continue;
		
		// Works by overwriting the first one's index
		for (int i=0; i<2*power; i++)
			OpLocSecond[OpLocs[i]] = i;
		
		// Copy out of the permutation iterator
		for (int i=0; i<2*power; i++)
			OpLocExtend[i] = OpLocs[i];
		
		// Bump up the second one of each pair
		for (int i=0; i<power; i++)
			OpLocExtend[OpLocSecond[i]] += power;
		
		// A permutation to use
		NumPermute++;
		ItmdSumProdMat.fill(0);
		
		do {
			for (int i=0; i<power; i++) {
				OpIndExtend[i] = OpIndices[i];
				OpIndExtend[i+power] = OpIndices[i] + 3;
			}
			
			// Need to loop over the whole length: (power) + extended (power)
			for (int iop = 0; iop < 2*power; iop++) {
				auto &OpCpnt = OpVec[OpIndExtend[OpLocExtend[iop]]];
				
				if (iop > 0) {
					// Multiply by operator
					NewProdMat.fill(0);
					for (int i=0; i<=twoj; i++)
						for (auto ent: OpCpnt) {
							NewProdMat.at(i,ent.ixs[1]) +=
								ProdMat.at(i,ent.ixs[0]) * ent.val;
						}
					
					// Copy in the new product
					ProdMat.CopyIn(NewProdMat);
					// No need to continue if one gets zero
					if (!IsMatrixNonzero(ProdMat)) break;
				} else {
					// Copy in the operator
					ProdMat.fill(0);
					for (auto ent: OpCpnt) {
						ProdMat.at(ent.ixs[0],ent.ixs[1]) = ent.val;
					}					
				}
				
			}			
			// Add in the operator product
			AddTo(ItmdSumProdMat, ProdMat, twoj+1);
			
		} while (OpIndices.Next());
		// Add in the intermediate sum
		AddTo(SumProdMat, ItmdSumProdMat, twoj+1);
				
		if (!is_symm) break;
		
	} while (OpLocs.Next());
	
	ScaleByNumPermute(SumProdMat, NumPermute);
	return FindDiagonal(SumProdMat, twoj+1);
}

const int NUMBER_OF_METHODS = 5;

double (*FindInvariantFuncs[NUMBER_OF_METHODS])(int, int, bool) =
{FindInvariant0, FindInvariant1, FindInvariant2, FindInvariant3, FindInvariant4};

int main(int argc, char **argv) {

	if (argc <= 4) {
		printf("Need: maximum two*j, power, is symmetric (0,1), method (0,1,2,3,4) \n");
		return 0;
	}
	
	const auto twoj_max = atoi(argv[1]);
	printf("Maximum two*j = %d\n",twoj_max);
	
	const auto power = atoi(argv[2]);
	printf("Power = %d\n",power);
	
	const auto is_symm_input = atoi(argv[3]);
	auto is_symm = is_symm_input != 0;
	printf("Is Symmetric = %d %s\n",is_symm_input,is_symm ? "True" : "False");
	
	const auto method = atoi(argv[4]);
	printf("Method = %d\n",method);
	if (method < 0 || method >= NUMBER_OF_METHODS) {
		printf("Method not supported\n");
		return 0;
	}
	
	printf("\n");
	
	std::vector<double> InvarList;
	
	printf("Invariant values:\n");
	for (int twoj=0; twoj<=twoj_max; twoj++) {
		
		auto FindInvariant = FindInvariantFuncs[method];
		
		double Invar = FindInvariant(twoj, power, is_symm);
		InvarList.push_back(Invar);
		printf("%d - %.16lg\n", twoj, Invar);
	}
	
	printf("\n");
	for (double Invar: InvarList) {printf(" %.16g",Invar);}
	printf("\n");
		
	return 0;
}