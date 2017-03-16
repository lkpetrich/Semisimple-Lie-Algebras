#ifndef LINEAR_ALGEBRA_TEMPLATE_CLASSES
#define LINEAR_ALGEBRA_TEMPLATE_CLASSES

/*
	Template class for doing a matrix and
	template functions for doing linear-algebra operations
	
	Uses STL vector for data storage, vectors
*/

#include <vector>
#include <map>
#include <algorithm>
using namespace std;

// Vector elementary operations
// One must allocate space for the results in advance


template<class N>
inline size_t index_min_v(N &vmin, N *v, size_t n)
{
	size_t ix = 0;
	vmin = v[ix];
	for (size_t i=1; i<n; i++)
	{
		N &val = v[i];
		if (val < vmin)
		{
			ix = i;
			vmin = val;
		}
	}
	return ix;
}

template<class N>
inline size_t index_min_v(N &vmin, vector<N> &v)
{
	return index_min_v(vmin, &v[0], v.size());
}

template<class N>
inline size_t index_max_v(N &vmax, N *v, size_t n)
{
	size_t ix = 0;
	vmax = v[ix];
	for (size_t i=1; i<n; i++)
	{
		N &val = v[i];
		if (val > vmax)
		{
			ix = i;
			vmax = val;
		}
	}
	return ix;
}

template<class N>
inline size_t index_max_v(N &vmax, vector<N> &v)
{
	return index_max_v(vmax, &v[0], v.size());
}


template<class N, class NV>
inline void sum(N &sres, const NV *v, size_t n)
{
	sres = N(0);
	for (size_t i=0; i<n; i++)
		sres += v[i];
}
template<class N, class NV>
inline void sum(N &sres, vector<NV> &v)
	{sum(sres,&v[0],v.size());}


template<class N, class NV1, class NV2>
inline void add_vv(N *vres, const NV1 *v1, const NV2 *v2, size_t n)
{
	for (size_t i=0; i<n; i++)
		vres[i] = v1[i] + v2[i];
}

template<class N, class NV1, class NV2>
inline void add_vv(vector<N> &vres, const NV1 *v1, const NV2 *v2)
	{add_vv(&vres[0],v1,v2,vres.size());}

template<class N, class NV1, class NV2>
inline void add_vv(N *vres, vector<NV1> &v1, const NV2 *v2)
	{add_vv(vres,&v1[0],v2,v1.size());}

template<class N, class NV1, class NV2>
inline void add_vv(vector<N> &vres, vector<NV1> &v1, const NV2 *v2)
	{add_vv(&vres[0],&v1[0],v2,vres.size());}

template<class N, class NV1, class NV2>
inline void add_vv(N *vres, const NV1 *v1, vector<NV2> &v2)
	{add_vv(vres,v1,&v2[0],v2.size());}

template<class N, class NV1, class NV2>
inline void add_vv(vector<N> &vres, const NV1 *v1, vector<NV2> &v2)
	{add_vv(&vres[0],v1,&v2[0],vres.size());}

template<class N, class NV1, class NV2>
inline void add_vv(N *vres, vector<NV1> &v1, vector<NV2> &v2)
	{add_vv(vres,&v1[0],&v2[0],v1.size());}

template<class N, class NV1, class NV2>
inline void add_vv(vector<N> &vres, vector<NV1> &v1, vector<NV2> &v2)
	{add_vv(&vres[0],&v1[0],&v2[0],vres.size());}


template<class N, class NV1, class NV2>
inline void sub_vv(N *vres, const NV1 *v1, const NV2 *v2, size_t n)
{
	for (size_t i=0; i<n; i++)
		vres[i] = v1[i] - v2[i];
}

template<class N, class NV1, class NV2>
inline void sub_vv(vector<N> &vres, const NV1 *v1, const NV2 *v2)
	{sub_vv(&vres[0],v1,v2,vres.size());}

template<class N, class NV1, class NV2>
inline void sub_vv(N *vres, vector<NV1> &v1, const NV2 *v2)
	{sub_vv(vres,&v1[0],v2,v1.size());}

template<class N, class NV1, class NV2>
inline void sub_vv(vector<N> &vres, vector<NV1> &v1, const NV2 *v2)
	{sub_vv(&vres[0],&v1[0],v2,vres.size());}

template<class N, class NV1, class NV2>
inline void sub_vv(N *vres, const NV1 *v1, vector<NV2> &v2)
	{sub_vv(vres,v1,&v2[0],v2.size());}

template<class N, class NV1, class NV2>
inline void sub_vv(vector<N> &vres, const NV1 *v1, vector<NV2> &v2)
	{sub_vv(&vres[0],v1,&v2[0],vres.size());}

template<class N, class NV1, class NV2>
inline void sub_vv(N *vres, vector<NV1> &v1, vector<NV2> &v2)
	{sub_vv(vres,&v1[0],&v2[0],v1.size());}

template<class N, class NV1, class NV2>
inline void sub_vv(vector<N> &vres, vector<NV1> &v1, vector<NV2> &v2)
	{sub_vv(&vres[0],&v1[0],&v2[0],vres.size());}

template<class N, class NV>
inline void addto_vv(N *vres, const NV *v, size_t n)
{
	for (size_t i=0; i<n; i++)
		vres[i] += v[i];
}

template<class N, class NV>
inline void addto_vv(vector<N> &vres, const NV *v)
	{addto_vv(&vres[0],v,vres.size());}

template<class N, class NV>
inline void addto_vv(N *vres, vector<NV> &v)
	{addto_vv(vres,&v[0],v.size());}

template<class N, class NV>
inline void addto_vv(vector<N> &vres, vector<NV> &v)
	{addto_vv(&vres[0],&v[0],vres.size());}

template<class N, class NV>
inline void subfm_vv(N *vres, const NV *v, size_t n)
{
	for (size_t i=0; i<n; i++)
		vres[i] -= v[i];
}

template<class N, class NV>
inline void subfm_vv(vector<N> &vres, const NV *v)
	{subfm_vv(&vres[0],v,vres.size());}

template<class N, class NV>
inline void subfm_vv(N *vres, vector<NV> &v)
	{subfm_vv(vres,&v[0],v.size());}

template<class N, class NV>
inline void subfm_vv(vector<N> &vres, vector<NV> &v)
	{subfm_vv(&vres[0],&v[0],vres.size());}

template<class N, class NS, class NV>
inline void mul_sv(N *vres, const NS &s, const NV *v, size_t n)
{
	for (size_t i=0; i<n; i++)
		vres[i] = s*v[i];
}

template<class N, class NS, class NV>
inline void mul_sv(vector<N> &vres, const NS &s, const NV *v)
	{mul_sv(&vres[0],v,vres.size());}

template<class N, class NS, class NV>
inline void mul_sv(N *vres, const NS &s, vector<NV> &v)
	{mul_sv(vres,s,&v[0],v.size());}

template<class N, class NS, class NV>
inline void mul_sv(vector<N> &vres, const NS &s, vector<NV> &v)
	{mul_sv(&vres[0],s,&v[0],vres.size());}

template<class N, class NS>
inline void mulby_sv(N *vres, const NS &s, size_t n)
{
	for (size_t i=0; i<n; i++)
		vres[i] *= s;
}

template<class N, class NS>
inline void mulby_sv(vector<N> &vres, const NS &s)
	{mulby_sv(&vres[0],s,vres.size());}

template<class N, class NS, class NV>
inline void div_sv(N *vres, const NS &s, const NV *v, size_t n)
{
	for (size_t i=0; i<n; i++)
		vres[i] = v[i]/s;
}

template<class N, class NS, class NV>
inline void div_sv(vector<N> &vres, const NS &s, const NV *v)
	{div_sv(&vres[0],v,vres.size());}

template<class N, class NS, class NV>
inline void div_sv(N *vres, const NS &s, vector<NV> &v)
	{div_sv(vres,s,&v[0],v.size());}

template<class N, class NS, class NV>
inline void div_sv(vector<N> &vres, const NS &s, vector<NV> &v)
	{div_sv(&vres[0],s,&v[0],vres.size());}

template<class N, class NS>
inline void divby_sv(N *vres, const NS &s, size_t n)
{
	for (size_t i=0; i<n; i++)
		vres[i] /= s;
}

template<class N, class NS>
inline void divby_sv(vector<N> &vres, const NS &s)
	{divby_sv(&vres[0],s,vres.size());}

template<class N, class NV1, class NV2>
inline void mul_vv(N &sres, const NV1 *v1, const NV2 *v2, size_t n)
{
	sres = N(0);
	for (size_t i=0; i<n; i++)
		sres += v1[i]*v2[i];
}

template<class N, class NV1, class NV2>
inline void mul_vv(N &sres, vector<NV1> &v1, const NV2 *v2)
	{mul_vv(sres,&v1[0],v2,v1.size());}

template<class N, class NV1, class NV2>
inline void mul_vv(N &sres, const NV1 *v1, vector<NV2> &v2)
	{mul_vv(sres,v1,&v2[0],v2.size());}

template<class N, class NV1, class NV2>
inline void mul_vv(N &sres, vector<NV1> &v1, vector<NV2> &v2)
	{mul_vv(sres,&v1[0],&v2[0],v1.size());}

template<class N, class NS, class NV>
inline void muladdto_vsv(N *vres, NS s, const NV *v, size_t n)
{
	for (size_t i=0; i<n; i++)
		vres[i] += s*v[i];
}

template<class N, class NS, class NV>
inline void muladdto_vsv(vector<N> &vres, NS s, const NV *v)
	{muladdto_vsv(&vres[0],s,v,vres.size());}

template<class N, class NS, class NV>
inline void muladdto_vsv(N *vres, NS s, vector<NV> &v)
	{muladdto_vsv(vres,s,&v[0],v.size());}

template<class N, class NS, class NV>
inline void muladdto_vsv(vector<N> &vres, NS s, vector<NV> &v)
	{muladdto_vsv(&vres[0],s,&v[0],vres.size());}


// vec1 == vec2 ?
template<class N1, class N2> inline bool VecEqual(const N1 *vec1, const N2 *vec2, size_t vlen)
{
	for (size_t i=0; i<vlen; i++)
	{
		if (vec1[i] != vec2[i]) return false;
	}
	return true;
}

template<class N1, class N2> inline bool VecEqual(vector<N1> &vec1, const N2 *vec2)
	{return VecEqual(&vec1[0],vec2,vec1.size());}

template<class N1, class N2> inline bool VecEqual(const N1 *vec1, vector<N2> &vec2)
	{return VecEqual(vec1,&vec2[0],vec2.size());}

template<class N1, class N2> inline bool VecEqual(vector<N1> &vec1, vector<N2> &vec2)
	{return VecEqual(&vec1[0],&vec2[0],vec1.size());}

// vec1 != vec2 ?
template<class N1, class N2> inline bool VecNotEqual(const N1 *vec1, const N2 *vec2, size_t vlen)
{
	for (size_t i=0; i<vlen; i++)
	{
		if (vec1[i] != vec2[i]) return true;
	}
	return false;
}

template<class N1, class N2> inline bool VecNotEqual(vector<N1> &vec1, const N2 *vec2)
	{return VecNotEqual(&vec1[0],vec2,vec1.size());}

template<class N1, class N2> inline bool VecNotEqual(const N1 *vec1, vector<N2> &vec2)
	{return VecNotEqual(vec1,&vec2[0],vec2.size());}

template<class N1, class N2> inline bool VecNotEqual(vector<N1> &vec1, vector<N2> &vec2)
	{return VecNotEqual(&vec1[0],&vec2[0],vec1.size());}

// vec1 > vec2 ?
template<class N1, class N2> inline bool VecGreaterThan(const N1 *vec1, const N2 *vec2, size_t vlen)
{
	for (size_t i=0; i<vlen; i++)
	{
		if (vec1[i] > vec2[i]) return true;
		if (vec1[i] < vec2[i]) return false;
	}
	return false;
}

template<class N1, class N2> inline bool VecGreaterThan(vector<N1> &vec1, const N2 *vec2)
	{return VecGreaterThan(&vec1[0],vec2,vec1.size());}

template<class N1, class N2> inline bool VecGreaterThan(const N1 *vec1, vector<N2> &vec2)
	{return VecGreaterThan(vec1,&vec2[0],vec2.size());}

template<class N1, class N2> inline bool VecGreaterThan(vector<N1> &vec1, vector<N2> &vec2)
	{return VecGreaterThan(&vec1[0],&vec2[0],vec1.size());}

// vec1 < vec2 ?
template<class N1, class N2> inline bool VecLessThan(const N1 *vec1, const N2 *vec2, size_t vlen)
{
	for (size_t i=0; i<vlen; i++)
	{
		if (vec1[i] > vec2[i]) return false;
		if (vec1[i] < vec2[i]) return true;
	}
	return false;
}

template<class N1, class N2> inline bool VecLessThan(vector<N1> &vec1, const N2 *vec2)
	{return VecLessThan(&vec1[0],vec2,vec1.size());}

template<class N1, class N2> inline bool VecLessThan(const N1 *vec1, vector<N2> &vec2)
	{return VecLessThan(vec1,&vec2[0],vec2.size());}

template<class N1, class N2> inline bool VecLessThan(vector<N1> &vec1, vector<N2> &vec2)
	{return VecLessThan(&vec1[0],&vec2[0],vec1.size());}

// vec1 >= vec2 ?
template<class N1, class N2> inline bool VecGreaterEqual(const N1 *vec1, const N2 *vec2, size_t vlen)
{
	for (size_t i=0; i<vlen; i++)
	{
		if (vec1[i] > vec2[i]) return true;
		if (vec1[i] < vec2[i]) return false;
	}
	return true;
}

template<class N1, class N2> inline bool VecGreaterEqual(vector<N1> &vec1, const N2 *vec2)
	{return VecGreaterEqual(&vec1[0],vec2,vec1.size());}

template<class N1, class N2> inline bool VecGreaterEqual(const N1 *vec1, vector<N2> &vec2)
	{return VecGreaterEqual(vec1,&vec2[0],vec2.size());}

template<class N1, class N2> inline bool VecGreaterEqual(vector<N1> &vec1, vector<N2> &vec2)
	{return VecGreaterEqual(&vec1[0],&vec2[0],vec1.size());}

// vec1 <= vec2 ?
template<class N1, class N2> inline bool VecLessEqual(const N1 *vec1, const N2 *vec2, size_t vlen)
{
	for (size_t i=0; i<vlen; i++)
	{
		if (vec1[i] > vec2[i]) return false;
		if (vec1[i] < vec2[i]) return true;
	}
	return true;
}

template<class N1, class N2> inline bool VecLessEqual(vector<N1> &vec1, const N2 *vec2)
	{return VecLessEqual(&vec1[0],vec2,vec1.size());}

template<class N1, class N2> inline bool VecLessEqual(const N1 *vec1, vector<N2> &vec2)
	{return VecLessEqual(vec1,&vec2[0],vec2.size());}

template<class N1, class N2> inline bool VecLessEqual(vector<N1> &vec1, vector<N2> &vec2)
	{return VecLessEqual(&vec1[0],&vec2[0],vec1.size());}


// Does a rectangular matrix
template<class N> class Matrix
{
	vector<N> contents;
	size_t rows, cols;
	
public:
	Matrix(): rows(0), cols(0) {}
	Matrix(size_t sqsize) {resize(sqsize);}
	Matrix(size_t rows_, size_t cols_) {resize(rows_,cols_);}
	
	// Copying: deep copy, so it can be changed independently of the original
	Matrix(const Matrix<N> &Mat);
	Matrix<N>& operator= (const Matrix<N> &Mat);
	
	size_t get_rows() {return rows;}
	size_t get_cols() {return cols;}
	void resize(size_t rows_, size_t cols_);
	void resize(size_t sqsize) {resize(sqsize,sqsize);}
	void clear() {resize(0,0);}
	
	bool IsSquare() {return rows == cols;}
	
	// operator[] accepts only one argument
	N &operator() (size_t irow, size_t icol)
		{return contents[cols*irow + icol];}
	
	void fill(const N &fillval)
		{::fill(contents.begin(),contents.end(),fillval);}
	void swap(Matrix<N> &Mat)
		{contents.swap(Mat.contents); ::swap(rows,Mat.rows); ::swap(cols,Mat.cols);}
	
	void AppendVector(const N *vec);
	void AppendVector(const vector<N> &vec)
		{AppendVector(&vec[0]);}
	
	void AppendMatrix(Matrix &mat);
};

template<class N> void Matrix<N>::resize(size_t rows_, size_t cols_)
{
	rows = rows_; cols = cols_;
	contents.resize(rows*cols);
}

template<class N> Matrix<N>::Matrix(const Matrix<N> &Mat)
{
	resize(Mat.rows,Mat.cols);
	copy(Mat.contents.begin(),Mat.contents.end(),contents.begin());
}

template<class N> Matrix<N>& Matrix<N>::operator= (const Matrix<N> &Mat)
{
	resize(Mat.rows,Mat.cols);
	copy(Mat.contents.begin(),Mat.contents.end(),contents.begin());
	return *this;
}

template<class N0, class N> void matcopy(Matrix<N> &dest, Matrix<N0> &src)
{
	size_t nr = src.get_rows();
	size_t nc = src.get_cols();
	dest.resize(nr,nc);
	for (size_t ir=0; ir<nr; ir++)
		for (size_t ic=0; ic<nc; ic++)
			dest(ir,ic) = src(ir,ic);
}

template<class N0, class N> void transpose(Matrix<N> &dest, Matrix<N0> &src)
{
	size_t nr = src.get_rows();
	size_t nc = src.get_cols();
	dest.resize(nc,nr);
	for (size_t ir=0; ir<nr; ir++)
		for (size_t ic=0; ic<nc; ic++)
			dest(ic,ir) = src(ir,ic);
}

template<class N> void Matrix<N>::AppendVector(const N *vec)
{
	for (size_t i=0; i<cols; i++)
		contents.push_back(vec[i]);
	rows++;
}

template<class N> void Matrix<N>::AppendMatrix(Matrix &mat)
{
	for (size_t i = 0; i < mat.rows; i++)
		AppendVector(&mat(i,0));
}


template<class N> void swap(Matrix<N> &M1, Matrix<N> &M2)
	{M1.swap(M2);}


// Proxy class for a row
// Does not get invalidated by a matrix-contents reallocation
template<class N> class MatrixRow
{
	// Avoid stale-pointer bugs by pointing to the matrix object
	// instead of to its contents
	Matrix<N> *MatPtr;
	size_t irow, cols;
public:
	MatrixRow(Matrix<N> &Mat, size_t irow):
		MatPtr(&Mat), irow(irow), cols(Mat.get_cols()) {}
	
	// Which row, how many columns
	size_t get_irow() {return irow;}
	size_t get_cols() {return cols;}
	
	N &operator[] (size_t icol)
		{Matrix<N> &Mat = *MatPtr; return Mat(irow,icol);}
	
	// Imitate STL container objects
	size_t size() {return cols;}
	N *begin() {return &(*this)[0];}
	N *end() {return &(*this)[size()];}
	N *rbegin() {return begin()-1;}
	N *rend() {return end()-1;}
	N &front() {return (*this)[0];}
	N &back() {return (*this)[size()-1];}
	
	// Which row of the original matrix?
	size_t row() {return irow;}
};


// Matrix elementary operations
// One must allocate space for the results in advance
// Assumes the "right" size for all the results

template<class N, class NS, class NM>
inline void mul_sm(Matrix<N> &mres, const NS &s, Matrix<NM> &m)
{
	for (size_t i=0; i<m.get_rows(); i++)
		for (size_t j=0; j<m.get_cols(); j++)
			mres(i,j) = s*m(i,j);
}

template<class N, class NS>
inline void mulby_sm(Matrix<N> &mres, const NS &s)
{
	for (size_t i=0; i<mres.get_rows(); i++)
		for (size_t j=0; j<mres.get_cols(); j++)
			mres(i,j) *= s;
}

template<class N, class NS, class NM>
inline void div_sm(Matrix<N> &mres, const NS &s, Matrix<NM> &m)
{
	for (size_t i=0; i<m.get_rows(); i++)
		for (size_t j=0; j<m.get_cols(); j++)
			mres(i,j) = s*m(i,j);
}

template<class N, class NS>
inline void divby_sm(Matrix<N> &mres, const NS &s)
{
	for (size_t i=0; i<mres.get_rows(); i++)
		for (size_t j=0; j<mres.get_cols(); j++)
			mres(i,j) /= s;
}

template<class N, class NM, class NV>
inline void mul_mv(N *vres, Matrix<NM> &m, const NV *v)
{
	for (size_t i=0; i<m.get_rows(); i++)
	{
		N sum = N(0);
		for (size_t j=0; j<m.get_cols(); j++)
			sum += m(i,j)*v[j];
		vres[i] = sum;
	}

}

template<class N, class NM, class NV>
inline void mul_mv(vector<N> &vres, Matrix<NM> &m, const NV *v)
	{mul_mv(&vres[0],m,v);}

template<class N, class NM, class NV>
inline void mul_mv(N *vres, Matrix<NM> &m, vector<NV> &v)
	{mul_mv(vres,m,&v[0]);}

template<class N, class NM, class NV>
inline void mul_mv(vector<N> &vres, Matrix<NM> &m, vector<NV> &v)
	{mul_mv(&vres[0],m,&v[0]);}

template<class N, class NV, class NM>
inline void mul_vm(N *vres, const NV *v, Matrix<NM> &m)
{
	for (size_t i=0; i<m.get_cols(); i++)
	{
		N sum = N(0);
		for (size_t j=0; j<m.get_rows(); j++)
			sum += m(j,i)*v[j];
		vres[i] = sum;
	}
}

template<class N, class NV, class NM>
inline void mul_vm(vector<N> &vres, const NV *v, Matrix<NM> &m)
	{mul_vm(&vres[0],v,m);}

template<class N, class NV, class NM>
inline void mul_vm(N *vres, vector<NV> &v, Matrix<NM> &m)
	{mul_vm(vres,&v[0],m);}

template<class N, class NV, class NM>
inline void mul_vm(vector<N> &vres, vector<NV> &v, Matrix<NM> &m)
	{mul_vm(&vres[0],&v[0],m);}

template<class N, class NV1, class NM, class NV2>
inline void mul_vmv(N &sres, const NV1 *v1, Matrix<NM> &m, const NV2 *v2)
{
	sres = N(0);
	for (size_t i=0; i<m.get_rows(); i++)
		for (size_t j=0; j<m.get_cols(); j++)
			sres += v1[i]*m(i,j)*v2[j];
}

template<class N, class NV1, class NM, class NV2>
inline void mul_vmv(N &sres, vector<NV1> &v1, Matrix<NM> &m, const NV2 *v2)
	{mul_vmv(sres,&v1[0],m,v2);}

template<class N, class NV1, class NM, class NV2>
inline void mul_vmv(N &sres, const NV1 *v1, Matrix<NM> &m, vector<NV2> &v2)
	{mul_vmv(sres,v1,m,&v2[0]);}

template<class N, class NV1, class NM, class NV2>
inline void mul_vmv(N &sres, vector<NV1> &v1, Matrix<NM> &m, vector<NV2> &v2)
	{mul_vmv(sres,&v1[0],m,&v2[0]);}

template<class N, class NM1, class NM2>
inline void mul_mm(Matrix<N> &mres, Matrix<NM1> &m1, Matrix<NM2> &m2)
{
	for (size_t i=0; i<m1.get_rows(); i++)
		for (size_t j=0; j<m2.get_cols(); j++)
		{
			N sum = N(0);
			for (int k=0; k<m1.get_cols(); k++)
				sum += m1(i,k)*m2(k,j);
			mres(i,j) = sum;
		}
}


// Vector-set indexer
// Rather kludgy, but it avoids having to write one's own red-black tree
// It uses STL map, which is usually implemented with a red-black tree.
// Template parameters:
// Type of data and hashcode function.
// That function takes an arg of data vector and returns a size_t integer.
// It's a class that contains an instance method operator() with an argument:
// vector<N> &vec

template<class N, class HashFunc>
struct VSRecord
{
	size_t hashcode; // For fast compare
	vector<N> data;
	
	void resize(size_t vlen) {data.resize(vlen); hashcode = 0;}
	void GetData(const N *vec);
};

template<class N, class HashFunc>
void VSRecord<N,HashFunc>::GetData(const N *vec)
{
	copy(vec,vec+data.size(),data.begin());
	HashFunc hf;
	hashcode = hf(data);
}

template<class N, class HashFunc>
struct VSRecLessThan
{
	bool operator() (const VSRecord<N,HashFunc> &rec1,
		const VSRecord<N,HashFunc> &rec2) const;
};

template<class N, class HashFunc>
bool VSRecLessThan<N,HashFunc>::operator() (const VSRecord<N,HashFunc> &rec1,
	const VSRecord<N,HashFunc> &rec2) const
{
	// Does not compare the index,
	// because a new vector may equal an old vector,
	// and we don't know in advance what the index is
	if (rec1.hashcode < rec2.hashcode) return true;
	if (rec1.hashcode > rec2.hashcode) return false;
	size_t sz1 = rec1.data.size();
	size_t sz2 = rec2.data.size();
	if (sz1 < sz2) return true;
	if (sz1 > sz2) return false;
	return VecLessThan(&rec1.data[0],&rec2.data[0],sz1);
}


template<class N, class HashFunc>
class VectorSet
{
	size_t vlen;
	VSRecord<N,HashFunc> rec;
	typedef map< VSRecord<N,HashFunc>, size_t, VSRecLessThan<N,HashFunc> >
		VSMap;
	typedef typename VSMap::iterator VSIter;
	VSMap SortedRecs;
	
public:
	// vlen: how long the vectors are
	VectorSet(): vlen(0) {}
	VectorSet(size_t vlen_): vlen(vlen_) {rec.resize(vlen);}
	size_t get_vlen() {return vlen;}
	void set_vlen(size_t vlen_)
		{vlen = vlen_; rec.resize(vlen);}
	
	void Reset() {SortedRecs.clear();}
	
	size_t size() {return SortedRecs.size();}
	
	// Checks on whether it is present
	bool CheckVector(const N *vec);
	bool CheckVector(const vector<N> &vec)
		{return CheckVector(&vec[0]);}
	
	// Returns the index of the vector of the set,
	// and appends the vector if it is not found in the set.
	// Also returns whether the vector was found in the set.
	pair<bool,size_t> AppendVector(const N *vec);
	pair<bool,size_t> AppendVector(const vector<N> &vec)
		{return AppendVector(&vec[0]);}
	
	// Like CheckVector, but it also returns the index
	// as part of the pair
	pair<bool,size_t> VectorIndex(const N *vec);
	pair<bool,size_t> VectorIndex(const vector<N> &vec)
		{return VectorIndex(&vec[0]);}
	
	// Exports in sorted order
	void ExportIndices(vector<size_t> &indxs);
	void ExportMatrix(Matrix<N> &mat);
};

template<class N, class HashFunc>
bool VectorSet<N,HashFunc>::CheckVector(const N *vec)
{
	rec.GetData(vec);
	VSIter it = SortedRecs.find(rec);
	
	return (it != SortedRecs.end());
}

template<class N, class HashFunc>
pair<bool,size_t> VectorSet<N,HashFunc>::AppendVector(const N *vec)
{
	rec.GetData(vec);
	VSIter it = SortedRecs.find(rec);
	
	pair<bool,size_t> ret;
	if (it == SortedRecs.end())
	{
		ret.first = false;
		ret.second = SortedRecs.size();
		SortedRecs[rec] = ret.second;
	}
	else
	{
		ret.first = true;
		ret.second = it->second;
	}
	return ret;
}

template<class N, class HashFunc>
pair<bool,size_t> VectorSet<N,HashFunc>::VectorIndex(const N *vec)
{
	rec.GetData(vec);
	VSIter it = SortedRecs.find(rec);
	
	pair<bool,size_t> ret;
	if (it == SortedRecs.end())
	{
		ret.first = false;
		ret.second = size_t(-1);
	}
	else
	{
		ret.first = true;
		ret.second = it->second;
	}
	return ret;
}

template<class N, class HashFunc>
void VectorSet<N,HashFunc>::ExportIndices(vector<size_t> &indxs)
{
	indxs.clear();
	VSIter it;
	
	for (it = SortedRecs.begin(); it != SortedRecs.end(); it++)
	{
		indxs.push_back(it->second);
	}
}

template<class N, class HashFunc>
void VectorSet<N,HashFunc>::ExportMatrix(Matrix<N> &mat)
{
	mat.resize(0,vlen);
	VSIter it;
	
	for (it = SortedRecs.begin(); it != SortedRecs.end(); it++)
	{
		mat.AppendVector(it->first.data);
	}
}


// Inverse of square matrix implemented as a template function
// Result's members are in class N, input's in class N0
// Returns whether the inversion could be performed:
// false for singular, non-square
template<class N, class N0> bool MatrixInverse(Matrix<N> &invmat, Matrix<N0> &mat)
{
	int nc = 0;
	if (!mat.IsSquare()) return false;
	
	size_t n = mat.get_rows();
	Matrix<N> workmat(n,2*n);
	workmat.fill(0);
	
	for (size_t i=0; i<n; i++)
	{
		// Initial matrix in first block
		for (size_t j=0; j<n; j++)
			workmat(i,j) = mat(i,j);
		// Identity matrix in second block
		workmat(i,n+i) = 1;
	}
	
	// Index vector for pivoting
	vector<size_t> ixp(n);
	for (size_t i=0; i<n; i++) ixp[i] = i;
	
	// Do forward substitution
	for (size_t icol=0; icol<n; icol++)
	{
		if (workmat(ixp[icol],icol) == N(0))
		{
			bool PivotFound = false;
			size_t ipvt;
			for (size_t i=(icol+1); i<n; i++)
				if (workmat(ixp[i],icol) != N(0))
				{
					PivotFound = true;
					ipvt = i;
					break;
				}
			// Singular?
			if (!PivotFound) return false;
			swap(ixp[ipvt],ixp[icol]);
		}
		// Make diagonal 1
		N dgvrecip = N(1)/workmat(ixp[icol],icol);
		for (size_t i=0; i<2*n; i++)
			workmat(ixp[icol],i) *= dgvrecip;
		// Forward substitute
		for (size_t i=icol+1; i<n; i++)
		{
			N elimval = workmat(ixp[i],icol);
			for (size_t j=icol; j<2*n; j++)
				workmat(ixp[i],j) -= elimval*workmat(ixp[icol],j);
		}
	}
	
	// Do back substitution
	for (size_t iclr=0; iclr<n; iclr++)
	{
		size_t icol = (n-1) - iclr;
		for (size_t i=0; i<icol; i++)
		{
			N elimval = workmat(ixp[i],icol);
			for (size_t j=0; j<2*n; j++)
				workmat(ixp[i],j) -= elimval*workmat(ixp[icol],j);
		}
	}
	
	// Done!
	invmat.resize(n);
	for (size_t i=0; i<n; i++)
		for (size_t j=0; j<n; j++)
			invmat(i,j) = workmat(ixp[i],n+j);
	
	return true;
}

#endif
