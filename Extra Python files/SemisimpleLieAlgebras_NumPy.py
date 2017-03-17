#!/usr/bin/env python
#
# Does calculations for semisimple Lie algebras
# Port of Mathematica notebook Semisimple Lie Algebras.nb
#
# Port of vector and array calculations to NumPy (Numerical Python)
# as much as possible.
#
# Uses tuples as much as possible instead of lists,
# because they are immutable and hashable
#
# Algebras specified here as (family number, Cartan-subalgebra dimension)
# where family number is 1,2,3,4,5,6,7 for A,B,C,D,E,F,G
#
# To get a Lie algebra, do
# GetLieAlgebra(latype)
# For valid latype, it will return its Lie-algebra object
# For invalid latype, it returns None
# It also caches each algebra object it creates
#
# latype is the algebra type, a list with members
# (family, rank or Cartan dimension)
# family is 1 thru 7: A,B,C,D,E,F,G
# Thus, E6 is (5,6)
#
# The cache is the global variable LieAlgebraCache, a dict.
# It is keyed by algebra type, and to clear it, set it to {}
#
# Members of the Lie-algebra object:
# name -- name as letter-number with SU/SO/Sp version where it exists
# special -- dict with key being type, value being max-weight vector
#    The type is "fundamental", "adjoint", "vector", "spinor", etc.
# dynkin -- Dynkin diagram with format
#   list of (root lengths, root connections)
#   The root connections have format (1st root, 2nd root, strength)
#   Has 1-based indexing
# metric -- for root.metric.root
# invmet -- inverse of metric
# imetnum -- numerator of integerized invmet
# imetden -- denominator of integerized invmet
# ctnmat -- Cartan matrix
# invctn -- inverse of Cartan matrix
# ictnnum -- numerator of integerized invctn
# ictnden -- denominator of integerized invctn
# invctn -- inverse of Cartan matrix
# posroots -- positive roots
# posrootsum -- sum of positive roots
#
# The integerized versions are present so that
# calculations can be done with integer arithmetic as much as possible
#
# str(algebra object) returns a description:
# algebra type in number-list and letter-number form
# its dimension
#
# The dimension can be returned by the member function dimension()
#
# To calculate the degeneracies, root vectors, and weight vectors of irreps
# (irreducible representations), take the Lie-algebra type
# and the highest-weight vector maxwts and do
# GetRep(latype, maxwts)
#
# Original code for getting a rep:
# GetRepDirect(latype, maxwts)
#
# For a Weyl orbit, emits the root vectors and weight vectors.
# maxwts is the dominant weights of the orbit,
# analogous to the highest weights of an irrep.
# GetOrbit(latype, maxwts)
#
# An irrep's Weyl-orbit content, in the format of GetRep(latype, maxwts)
# GetRepOrbits(latype, maxwts)
#
# Overall list of properties of a rep
# RepProperties(latype, maxwts)
# Call with a type and the rep's max-weight vector
# Returns a dict with these keys and values:
#
# type: LA type
# maxwts: max-weights vector
# mwconjg: its conjugate
# isselfconjg: whether those two are equal
# dimension: total degeneracy or dimension of the irrep
# height: max - min of sum(values in each root vector)
# reality: real(0), pseudoreal(1), or complex(-1)
# conserved: conserved quantities
#   A tuple of (modulus/divisor, conserved qty)
# Conserved quantities are the same for all the members of a representation,
# and they are modulo additive for product reps
#
# Casimir invariant:
# CasimirInvariant(latype, maxwts)
# Representation index:
# RepIndex(latype, maxwts)
#
# 
# For a list of irreps, specified as a list of max-weight vectors.
# It produces the same kind of output as GetRep().
# GetRepList(latype, mwlist)
# 
# For a counted list of irreps, specified as a list of
# (count, max-weight vector). It produces the same kind of output as GetRep().
# GetRepCntdList(latype, mwclist)
# 
# For a general rep argument, specified with "rptype". The type is "sngl"
# for a single irrep, "list" for a list of irreps, and "cntd"
# for a counted list of irreps.
# GetRepXtnd(latype, rptype, maxwts)
# 
# These are for representations for algebra products.
# Their algebra arg is (la1, la2, ...), while their irreps are
# (mw1, mw2, ...,u1, u2, ...). mw1 is the max weights for algebra la1,
# mw2 is the max weights for algebra la2, and u1, u2,
# are additional U(1) factors. The root vectors are (r1, r2, ..., u1, u2, ...)
# and the weight vectors are (w1, w2, ..., u1, u2, ...),
# where r1 and w1 are for la1 and mw2, etc.
# GetAlgProdRep(latype, maxwts)
# 
# For a list of algebra-product irreps. Output is in the format
# of GetAlgProdRep().
# GetAlgProdRepList(latype, mwlist)
# 
# For a counted list of algebra-product irreps. Output is in the format
# of GetAlgProdRep().
# GetAlgProdRepCntdList(latype, mwclist)
# 
# For a general rep argument, specified with "rptype".
# The type is "sngl"  for a single irrep, "list" for a list of irreps,
# and "cntd" for a counted list of irreps.
# GetAlgProdRepXtnd(latlist, rptype, maxwts)
# 
# Adds up the degeneracies to find the total degeneracy (multiplicity)
# TotalDegenOfExpRep(x) 
# 
# Uses Weyl's celebrated formula; it's much faster than calculating
# an irrep's root/weight vectors
# TotalDegen(latype, maxwts) 
# Generalized:
# TotalDegenXtnd(latype, type, maxwts)
# AlgProdTotalDegen(latlist, maxwts)
# AlgProdTotalDegenXtnd(latlist, rptype, maxwts)
#
# Get the "true" forms, with fractional, non-integerized root values
# GetTrueRep(latype, rep)
# GetTrueAlgProdRep(latlist, rep)
#
#
# For decomposing a product of two irreps into its irrep content,
# use their highest weights here. Makes a counted list of irreps,
# a list of (degeneracy, irrep highest weight).
# DecomposeRepProduct(latype, maxwt1, maxwt2)
# 
# For more general types of highest weights, where  the types are:
# "sngl" is single one, "list" is a list, and "cntd" is a counted list.
# DecomposeRepProductXtnd(latype, rptype1, maxwt1, rptype2, maxwt2)
# 
# Uses an algebra product and highest weights consisting of a list
# of one for each algebra and U(1) factors
# DecomposeAlgProdRepProduct(latlist, maxwt1, maxwt2)
# 
# For more general types of reps, as before
# DecomposeAlgProdRepProductXtnd(latlist, rptype1, maxwt1, rptype2, maxwt2)
# 
# Alternates: these do rep products only on single irreps of each algebra,
# instead of for all of them together. They give the same results,
# though with different performance.
# DecomposeRepProductXtndAlt(latype, rptype1, maxwt1, rptype2, maxwt2)
# DecomposeAlgProdRepProductAlt(latlist, maxwt1, maxwt2)
# DecomposeAlgProdRepProductXtndAlt(latlist, rptype1, maxwt1, rptype2, maxwt2)
# 
# DecomposeRepProdList(latype, maxwtlist) works on a list of irreps' highest weights.
# 
# 
# For decomposing a rep power into parts with different symmetries
# for general positive-integer powers pwr.
# The symmetric part is the first part
# and the antisymmetric part is the last part,
# with (pwr-2) mixed parts in between.
# DecomposeRepPower(latype, maxwts, pwr)
#
# Uses the highest-weight type system: "sngl", "list", "cntd" (counted list)
# DecomposeRepPowerXtnd(latype, rptype, maxwts, pwr)
# 
# For algebra products
# DecomposeAlgProdRepPower(latlist, maxwts, pwr)
# 
# For algebra products and weight types
# DecomposeAlgProdRepPowerXtnd(latlist, rptype, maxwts, pwr)
# 
# Pure symmetric and antisymmetric cases: sym = +1 and -1
# These are analogous to DecomposeRepListPower(latype, mwlist, pwr), etc.
# 
# DecomposeRepPwrSym(latype, maxwts, pwr, sym)
# DecomposeRepPwrSymXtnd(latype, rptype, maxwts, pwr, sym)
# DecomposeAlgProdRepPwrSym(latlist, maxwts, pwr, sym)
# DecomposeAlgProdRepPwrSymXtnd(latlist, rptype, maxwts, pwr, sym)
#
#
# One can find details on the symmetry types with the function
# GetTensorPowerYDX(pwr)
# It returns an object with members
#   ydcnts
#   kostka
#   nesting
#
# ydcnts is the diagram list, and its entries are
#   (multiplicity in the power expansion, the diagram itself,
#   the diagram's length, a multinomial multiplicity factor)
# The overall object has these convenience functions for extracting
# each sort of entry:
# cnts(), yds(), lens(), mnfs()
#
# Overall, one might do something like
# GetTensorPowerYDX(pwr).cnts() or GetTensorPowerYDX(pwr).yds()
#
# kostka is the Kostka matrix is what Weyl orbits for each irrep of A(n),
# using the order of the diagram list.
#
# nesting is the rep-nesting matrix:
# each rep's numbers becomes a rep and one iterates over them.
# Row: index of starting rep, column: index of nested reps concatenated.
# At each row and column, a list of entries:
# Indices of subreps, multiplicities of subreps,
#   multinomial multiplicity factor.
# The reps here are A(n) irreps or Young diagrams
#
#
# Subalgebras and branching rules
#
# One first creates a branching-rule object
# Its members:
# latype: Original-algebra type
# stsms: list for each subalgebra of
#    (type, projection matrix, numerator of integerized projmat,
#        its denominator)
# u1s: Indices of U(1)-factor roots (1-based)
#    For no U(1) factors, use ()
#
# SubAlgebras(): the subalgebra types
# SubAlgebraNmes(): their names
#
# DoBranching(maxwts): method for doing the branching
# Call with max-weight vector for the rep to break down
# Returns a counted list, a tuple of
#   (degen, tuple of subalgebras' maxweight vectors with U(1) factors)
#
# DoBranchingXtnd:(rptype, maxwts):
# Like above, but with type "sngl" for a single irrep,
# "list" for a list of them, and "cntd" for a counted list of them
#
# Available branching-rule generators. Each one returns
# a branching-rule object
#
# Root Demotion
# MakeRootDemoter(latype,m)
# Call with original-algebra type and root to demote
# to a U(1) factor (1-based)
#
# ListRootDemotions(latype)
#    does what it says for all the roots
#
# Multiple Root Demotion
# MakeMultiRootDemoter(latype,dmrts)
# Like previous, but with a list of roots to demote
# to a U(1) factor
#
# Extension Splitting
# MakeExtensionSplitter(latype,m)
# Call with original-algebra type and
# which root to split the extended diagram at (1-based)
# Not implemented for A(n), B1, C2, D2, or D3
#
# ListExtensionSplits(latype)
#    does what it says for all the roots
#
# Additional ones:
#
# In these two, the matrices of the original groups get decomposed
# into outer products of the subgroup matrices,
# all in the vector representation.
#
# SubalgMultSU(suords)
# reduces a SU(n) to a product of SU(m)'s; suords is a list of those m's,
# and n = product(m's).
#
# SubalgMultSOSp(sospords)
# reduces a SO(n) or Sp(n) to a product of SO(m)'s and/or Sp(m's);
# sospords is a list of those m's, and n = product(m's).
# Positive n means SO(n) and negative n means Sp(-n).
# SO(2) is handled as a U(1) factor.
#
# Forms of the previous two for A,B,C,D-style designation.
#
# SubalgMultAn(ords)
# takes list of subalgebra root-vector lengths. The 1 in (1,n) is assumed.
# SubalgMultBCDn(stypes)
# takes list of subalgebra types that are Bn, Cn, and Dn: (2,n), (3,n), (4,n).
# SO(2) / D(1) / (4,1) is legitimate here; it becomes a U(1) factor.
#
# SubalgSOEvenOdd(n,m)
# breaks SO(2n) -> SO(2m+1) + SO(2n-2m-1) / D(n) -> B(m) + B(n-m-1)
# For m = 0, only does SO(2n-1) / B(n-1)
# The extension splitter does
# even -> even + even / D -> D + D
# odd -> odd + even / B -> B + D
# Use the root demoter and demote root #1
# to turn SO(n) into SO(2) * SO(n-2),
# that is, B(n'), D(n') -> B(n'-1), D(n'-1) + U(1) factor
#
# SubalgSUSO(n) -- turns SU(n) / A(n-1) into SO(n) / D(n/2) or B((n-1)/2)
# SubalgSUSp(n) -- turns SU(2n) / A(2n-1) into Sp(2n) / C(n)
#
# SubalgHeightA1(latype)
# reduces any algebra to A1, with the rep height becoming
# the largest highest weight.
#
# Some extra ones named individually.
#
# Mentioned by John Baez in "The Octonions":
#
# B3G2 -- B3/SO(7) to G2 -- G2 is the isomorphism group of the octonions,
# and one gets a construction of G2 from it that' s manifestly
# a subgroup of SO(7) -- 14 7D antisymmetric real matrices.
# 
# D4G2 -- D4/SO(8) to G2 (not maximal, but included for completeness)
#
# The others entioned by Slansky :
# G2A1 -- G2 to A1/SU(2)
# etc.
# Full list returned by SubalgExtraData.keys()
#
# None of these branchings have U(1) factors
#
# SubalgExtra(saname)
# Call with one of these character-string names
#
# These all return a brancher object from their input brancher objects
# and other data
# Indexing is 1-based
#
# ConcatBranchers(ld0, ix, ld1)
# subalgebra #ix of brancher ld0 gets branched by ld1,
# making a combined brancher.
#
# BrancherRenameA1B1C1(ld0, ix, newfam)
# subalgebra #ix of brancher ld0 gets a new family number
# if it's A(1), B(1), or C(1): 1, 2, 3
#
# BrancherRenameB2C2(ld0, ix)
# subalgebra #ix of brancher ld0 gets flipped between B(2) and C(2)
#
# BrancherRenameA3D3(ld0, ix)
# subalgebra #ix of brancher ld0 gets flipped between A(3) and D(3)
#
# BrancherSplitD2(ld0, ix)
# subalgebra #ix of brancher ld0 gets split into two A(1)'s if it is D(2)
#
# BrancherJoin2A1(ld0, ix1, ix2)
# subalgebras #ix1 and #ix2 of brancher ld0 get joined into D(2)
# if they are A(1)/B(1)/C(1)'s.
#
# BrancherConjugate(ld0, cjixs)
# makes conjugates of subalgebras  of brancher ld0 with indexes in cjixs;
# different for A(n), n>1, D(n), and E(6). Not a true conjugate for D(2n),
# but exchanged anyway.
#
# BrancherConjgD4(ld0, ix, newrts)
# subalgebra #ix gets conjugated with newrts specifiying
# the new roots' location if it is D(4).
# newrts is 1-based with length 4 with the second being 2
#
# BrancherRearrange[ld0, neword)
# puts the subalgebras of brancher ld0 into the order specified in neword.
#
# SubalgSelf(latype)
# returns branching to the original algebra
#
#
# References:
#
# Semi-Simple Lie Algebras and Their Representations, by Robert N. Cahn
# http://phyweb.lbl.gov/~rncahn/www/liealgebras/book.html
#
# Group Theory for Unified Model Building, by R. Slansky
# http://www-spires.slac.stanford.edu/spires/find/hep/www?j = PRPLC, 79, 1
#
# This code uses the root-numbering conventions of Robert N. Cahn's book
#
#
# Kostka-matrix algorithm from
# Determinantal Expression and Recursion for Jack Polynomials,
# by L. Lapointe, A. Lascoux, J. Morse
#
# Erroneous inputs are handled by raising exceptions
#
# Bignums are part of Python integers
# Rational numbers are the Fraction class
from fractions import Fraction, gcd

# Numerical Python
import numpy as np

# Utilities

# Convert an array from a list of lists to a tuple of tuples (immutable lists)
def MatrixTuple(x): return tuple(map(tuple,x))
def zeros1(n): return n*[0]
def zeros2(m,n): return [zeros1(n) for k in xrange(m)]
def identmat(n): return [[1 if j == i else 0 for j in xrange(n)] for i in xrange(n)]

# Create NumPy arrays
NPType = np.int32;
def MakeNP(x): return np.array(x, dtype=NPType)
def npzeros1(n): return np.zeros(n, dtype=NPType)
def npzeros2(m,n): return np.zeros((m,n), dtype=NPType)


# Various vector and matrix operations.
# addto is +=
# arg types: s = scalar, v = vector, m = matrix
# Uses integer 0 and 1

def transpose(mat):
	n1 = len(mat)
	n2 = len(mat[0])
	res = zeros2(n2,n1)
	for i1 in xrange(n1):
		for i2 in xrange(n2):
			res[i2][i1] = mat[i1][i2]
	return res

def add_vv(vec1,vec2):
	n = len(vec1)
	res = zeros1(n)
	for i in xrange(n):
		res[i] = vec1[i] + vec2[i]
	return res

def addto_vv(vec1,vec2):
	n = len(vec1)
	for i in xrange(n):
		vec1[i] += vec2[i]
	return vec1

def add_vvv(vec1,vec2,vec3):
	n = len(vec1)
	res = zeros1(n)
	for i in xrange(n):
		res[i] = vec1[i] + vec2[i] + vec3[i]
	return res

def sub_vv(vec1,vec2):
	n = len(vec1)
	res = zeros1(n)
	for i in xrange(n):
		res[i] = vec1[i] - vec2[i]
	return res

def subfm_vv(vec1,vec2):
	n = len(vec1)
	for i in xrange(n):
		vec1[i] -= vec2[i]
	return vec1

def mul_sv(scl,vec):
	return [scl*mem for mem in vec]

def mulby_sv(scl,vec):
	for k in xrange(len(vec)):
		vec[k] *= scl
	return vec

def div_sv(scl,vec):
	return [mem/scl for mem in vec]

def divby_sv(scl,vec):
	for k in xrange(len(vec)):
		vec[k] /= scl
	return vec

def mul_sm(scl,mat):
	return [mul_sv(scl,mem) for mem in mat]

def div_sm(scl,vec):
	return [div_sv(scl,mem) for mem in mat]

def mul_vv(vec,vecx):
	n = len(vec)
	ttl = 0
	for i in xrange(n):
		ttl += vec[i]*vecx[i]
	return ttl

def mul_vm(vec,mat):
	nx = len(vec)
	n = len(mat[0])
	res = zeros1(n)
	for i in xrange(n):
		ttl = 0
		for j in xrange(nx):
			ttl += vec[j]*mat[j][i]
		res[i] = ttl
	return res

def mul_mv(mat,vec):
	return [mul_vv(mtrow,vec) for mtrow in mat]

def mul_vmv(vec1,mat,vec2):
	n1 = len(vec1)
	n2 = len(vec2)
	res = 0
	for i1 in xrange(n1):
		vval = vec1[i1]
		mtrow = mat[i1]
		rsi = 0
		for i2 in xrange(n2):
			rsi += mtrow[i2]*vec2[i2]
		res += vval*rsi
	return res
	
def mul_mm(mat1,mat2):
	n1 = len(mat1)
	nx = len(mat1[0])
	n2 = len(mat2[0])
	res = zeros2(n1,n2)
	for i1 in xrange(n1):
		mtrow = mat1[i1]
		for i2 in xrange(n2):
			ttl = 0
			for j in xrange(nx):
				ttl += mtrow[j]*mat2[j][i2]
			res[i1][i2] = ttl
	return res

def muladdto_vsv(vec1,scl,vec2):
	n = len(vec1)
	for i in xrange(n):
		vec1[i] += scl*vec2[i]
	return vec1

# Do inverse with the Gauss-Jordan algorithm
# Starts with an integer matrix, and uses rational numbers
def MatrixInverse(mat):
	# Set up the work matrix: originally (original,identity)
	# Transform into (identity,inverse)
	n = len(mat)
	workmat = [(2*n)*[Fraction(0)] for k in xrange(n)]
	for i in xrange(n):
		mrow = mat[i]
		wmrow = workmat[i]
		for j in xrange(n):
			wmrow[j] += int(mrow[j])
	for k in xrange(n):
		workmat[k][n+k] += 1
	
	# Do forward substitution
	for icol in xrange(n):
		# Necessary to exchange rows
		# to bring a nonzero value into position?
		# Return None if singular
		if workmat[icol][icol] == 0:
			ipvt = None
			for i in xrange(icol+1,n):
				if workmat[i][icol] != 0:
					ipvt = i
					break
			if ipvt == None: return None
			temp = workmat[icol]
			workmat[icol] = workmat[ipvt]
			workmat[ipvt] = temp
		# Make diagonal 1:
		wmicol = workmat[icol]
		dgvrecip = 1/wmicol[icol]
		for i in xrange(icol,2*n):
			wmicol[i] *= dgvrecip
		# Forward substitute:
		for i in xrange(icol+1,n):
			wmi = workmat[i]
			elimval = wmi[icol]
			for j in xrange(icol,2*n):
				wmi[j] -= elimval*wmicol[j]
	
	# Do back substitution
	for icol in xrange(n-1,0,-1):
		wmicol = workmat[icol]
		for i in xrange(icol):
			wmi = workmat[i]
			elimval = wmi[icol]
			for j in xrange(icol,2*n):
				wmi[j] -= elimval*wmicol[j]
	
	# Done!
	return [[workmat[i][n+j] for j in xrange(n)] for i in xrange(n)]

# Find shared denominator of rational-number vectors and matrices;
# turn them into integer vectors and matrices
# Returns (integerized vector/matrix, shared denominator)
def lcm(a,b): return 1 if a*b == 0 else (a/gcd(a,b))*b

def SharedDen_Vector(vec):
	den = 1
	for mem in vec: den = lcm(den,mem.denominator)
	return ([int(den*mem) for mem in vec], den)

def SharedDen_Matrix(mat):
	matexp = [SharedDen_Vector(vec) for vec in mat]
	den = 1
	for vec,vdn in matexp: den = lcm(den,vdn)
	return ([mul_sv(den/vdn,vec) for vec,vdn in matexp], den)


# Lie-Algebra Setup:

# Algebra metric: pure integers
def LA_Metric(dynkin):
	rtwts = dynkin[0]
	n = len(rtwts)
	mat = npzeros2(n,n)
	for k,str in enumerate(rtwts):
		mat[k,k] = 2*str
	for conn in dynkin[1]:
		i = conn[0]-1
		j = conn[1]-1
		str = conn[2]
		mat[i,j] = - str
		mat[j,i] = - str
	return mat

# Cartan matrix: pure integers
def LA_Cartan(dynkin):
	rtwts = dynkin[0]
	n = len(rtwts)
	mat = npzeros2(n,n)
	for k in xrange(n):
		mat[k,k] = 2
	for conn in dynkin[1]:
		i = conn[0]-1
		j = conn[1]-1
		str = conn[2]
		mat[i,j] = - str/rtwts[j]
		mat[j,i] = - str/rtwts[i]
	return mat

# Sort by root values 
def RootWtSortFunc(rx1,rx2):
	r1 = rx1[0]; r2 = rx2[0];
	l1 = sum(r1); l2 = sum(r2)
	return cmp(tuple(r1),tuple(r2)) if l1 == l2 else cmp(l2,l1)

# Find positive roots from the Cartan Matrix
# They all have integer values
# Returns tuple of values of (root, weight)
def LA_PositiveRoots(ctnmat):
	# Initial positive roots
	n = len(ctnmat)
	PosRoots = zip(map(MakeNP,identmat(n)),ctnmat)
	
	# Test for whether we've found a root, and if so, return its index
	PRTest = {}
	for k,rwt in enumerate(PosRoots):
		rtkey = tuple(rwt[0])
		PRTest[rtkey] = k
	
	# Find the next root until one cannot find any more
	WhichWay = [np.ones(n,dtype=np.bool_) for k in xrange(n)]
	RtIndx = 0
	while RtIndx < len(PosRoots):
		rwt = PosRoots[RtIndx]
		ThisRoot = rwt[0]
		ThisWeight = rwt[1]
		ctpd = list(rwt[1])
		# Advance in each direction, if possible
		for i in xrange(n):
			ThisWW = WhichWay[RtIndx]
			NewWW = np.array([j != i for j in xrange(n)], dtype=np.bool_)
			if ThisWW[i]:
				for j in xrange(1,-ctpd[i]+1):
					# Calculate roots and weights in parallel,
					# to avoid repeated matrix.vector calculations
					NewRoot = np.copy(ThisRoot)
					NewRoot[i] += j
					NewWeight = np.copy(ThisWeight)
					NewWeight += j*ctnmat[i]
					# Avoid integer multiples of previous root vectors
					NRLen = np.sum(NewRoot)
					RootOK = True
					for NRDiv in xrange(2,NRLen):
						if (NRLen % NRDiv) != 0: continue
						DVI = True
						for rcmp in NewRoot:
							if rcmp % NRDiv != 0:
								DVI = False
								break
						if DVI:
							DivRoot = NewRoot/NRDiv
							if tuple(DivRoot) in PRTest:
								RootOK = False
								break
					if RootOK:
						# Add the root if possible
						rtkey = tuple(NewRoot)
						if rtkey in PRTest:
							WhichWay[RtIndx] = \
								np.logical_and(ThisWW, NewWW)
						else:
							PRTest[rtkey] = len(PosRoots)
							PosRoots.append((NewRoot,NewWeight))
							WhichWay.append(NewWW)
		RtIndx += 1
	
	PosRoots.sort(RootWtSortFunc)
	PosRoots.reverse()
	return tuple(PosRoots)

def AlgName(latype):
	return "%s%d" % (chr(ord('A')+(latype[0]-1)),latype[1])

def InvalidLATypeError(latype):
	return TypeError("Invalid Lie-Algebra type: %s" % str(latype))

# The Lie-algebra class itself
# If the algebra type is invalid, the constructor throws an exception.
#
# Members of the Lie-algebra object:
# name -- name as letter-number with SU/SO/Sp version where it exists
# special -- dict with key being type, value being max-weight vector
#    The type is "fundamental", "adjoint", "vector", "spinor", etc.
# dynkin -- Dynkin diagram with format : 
#   list of (root lengths, root connections)
#   The root connections have format (1st root, 2nd root, strength)
#   Has 1-based indexing
# metric -- for root.metric.root
# invmet -- inverse of metric
# imetnum -- numerator of integerized invmet
# imetden -- denominator of integerized invmet
# ctnmat -- Cartan matrix
# invctn -- inverse of Cartan matrix
# ictnnum -- numerator of integerized invctn
# ictnden -- denominator of integerized invctn
# invctn -- inverse of Cartan matrix
# posroots -- positive roots
# posrootsum -- sum of positive roots
#
# The integerized versions are present so that
# calculations can be done with integer arithmetic as much as possible
#
# __str__: has name, dimension
# __repr__: how to create a Lie-algebra object
# dimension() : the algebra's dimension
#
class LieAlgebra:
	def __init__(self,latype):
		if len(latype) != 2: raise InvalidLATypeError(latype)
		family = latype[0]
		n = latype[1]
		self.latype = latype
		self.special = {}
		if type(family) != type(0): raise InvalidLATypeError(latype)
		if type(n) != type(0): raise InvalidLATypeError(latype)
		if n < 1: raise InvalidLATypeError(latype)
		
		# Special irreps and Dynkin diagrams
		nmsfx = ""
		self.special["singlet"] = tuple(zeros1(n))
		if family == 1:
			# A(n)
			nmsfx = "SU(%d)" % (n+1)
			spcl = zeros1(n)
			if n == 1:
				spcl[0] = 2
			else:
				spcl[0] = 1; spcl[-1] = 1
			self.special["adjoint"] = tuple(spcl)
			spcl = zeros1(n)
			spcl[0] = 1
			self.special["vector"] = tuple(spcl)
			self.special["fundamental"] = tuple(spcl)
			spcl = zeros1(n)
			spcl[-1] = 1
			self.special["vector-mirror"] = tuple(spcl)
			self.special["fundamental-mirror"] = tuple(spcl)
			self.dynkin = (tuple(n*[1]), \
				tuple([(k, k+1, 1) for k in xrange(1,n)]))
		elif family == 2:
			# B(n)
			nmsfx = "SO(%d)" % (2*n+1)
			spcl = zeros1(n)
			if n == 1:
				spcl[0] = 2
			elif n == 2:
				spcl[1] = 2
			else:
				spcl[1] = 1
			self.special["adjoint"] = tuple(spcl)
			spcl = zeros1(n)
			if n == 1:
				spcl[0] = 2
			else:
				spcl[0] = 1
			self.special["vector"] = tuple(spcl)
			if n % 2 == 1 and n > 1: self.special["fundamental-2"] = tuple(spcl)
			spcl = zeros1(n)
			spcl[-1] = 1
			self.special["spinor"] = tuple(spcl)
			if n % 2 == 0 or n == 1: self.special["fundamental"] = tuple(spcl)
			if n % 2 == 1 and n > 1: self.special["fundamental-1"] = tuple(spcl)
			if n > 1:
				self.dynkin = (tuple(((n-1)*[2]) + [1]), \
					tuple([(k, k+1, 2) for k in xrange(1,n-1)] + [(n-1,n,2)]))
			else:
				self.dynkin = ((1,),())
		elif family == 3:
			# C(n)
			nmsfx = "Sp(%d)" % (2*n)
			spcl = zeros1(n)
			spcl[0] = 2
			self.special["adjoint"] = tuple(spcl)
			spcl = zeros1(n)
			spcl[0] = 1
			self.special["vector"] = tuple(spcl)
			self.special["fundamental"] = tuple(spcl)
			if n > 1:
				self.dynkin = (tuple(((n-1)*[1]) + [2]), \
					tuple([(k, k+1, 1) for k in xrange(1,n-1)] + [(n-1,n,2)]))
			else:
				self.dynkin = ((1,),())
		elif family == 4:
			if n < 2: raise InvalidLATypeError(latype)
			# D(n)
			nmsfx = "SO(%d)" % (2*n)
			spcl = zeros1(n)
			if n > 2:
				spcl = zeros1(n)
				if n == 3:
					spcl[1] = 1
					spcl[2] = 1
				else:
					spcl[1] = 1
				self.special["adjoint"] = tuple(spcl)
			else:
				spcl = zeros1(n)
				spcl[0] = 2
				self.special["adjoint-1"] = tuple(spcl)
				spcl = zeros1(n)
				spcl[1] = 2
				self.special["adjoint-2"] = tuple(spcl)
			spcl = zeros1(n)
			if n == 2:
				spcl[0] = 1
				spcl[1] = 1
			else:
				spcl[0] = 1
			self.special["vector"] = tuple(spcl)
			spcl = zeros1(n)
			spcl[-2] = 1
			self.special["spinor-1"] = tuple(spcl)
			if n % 2 == 0: self.special["fundamental-1"] = tuple(spcl)
			if n % 2 == 1: self.special["fundamental"] = tuple(spcl)
			spcl = zeros1(n)
			spcl[-1] = 1
			self.special["spinor-2"] = tuple(spcl)
			if n % 2 == 0: self.special["fundamental-2"] = tuple(spcl)
			if n % 2 == 1: self.special["fundamental-mirror"] = tuple(spcl)
			if n > 2:
				self.dynkin = (tuple(n*[1]), \
					tuple([(k, k+1, 1) for k in xrange(1,n-1)] + [(n-2,n,1)]))
			else:
				self.dynkin = ((1,1), ())
		elif family == 5:
			if n < 5 or n > 8: raise InvalidLATypeError(latype)
			if n == 6:
				# E6
				spcl = zeros1(n)
				spcl[-1] = 1
				self.special["adjoint"] = tuple(spcl)
				spcl = zeros1(n)
				spcl[-1] = 1
				self.special["fundamental"] = tuple(spcl)
				spcl = zeros1(n)
				spcl[-2] = 1
				self.special["fundamental-mirror"] = tuple(spcl)
			elif n == 7:
				# E7
				spcl = zeros1(n)
				spcl[0] = 1
				self.special["adjoint"] = tuple(spcl)
				spcl = zeros1(n)
				spcl[-2] = 1
				self.special["fundamental"] = tuple(spcl)
			elif n == 8:
				# E8
				spcl = zeros1(n)
				spcl[-2] = 1
				self.special["adjoint"] = tuple(spcl)
				self.special["fundamental"] = tuple(spcl)
			self.dynkin = (tuple(n*[1]), \
				tuple([(k, k+1, 1) for k in xrange(1,n-1)] + [(3,n,1)]))
		elif family == 6:
			if n != 4: raise InvalidLATypeError(latype)
			# F4
			spcl = zeros1(n)
			spcl[0] = 1
			self.special["adjoint"] = tuple(spcl)
			spcl = zeros1(n)
			spcl[-1] = 1
			self.special["fundamental"] = tuple(spcl)
			self.dynkin = ((2,2,1,1), ((1,2,2),(2,3,2),(3,4,1)))
		elif family == 7:
			if n != 2: raise InvalidLATypeError(latype)
			# G2
			spcl = zeros1(n)
			spcl[0] = 1
			self.special["adjoint"] = tuple(spcl)
			spcl = zeros1(n)
			spcl[-1] = 1
			self.special["fundamental"] = tuple(spcl)
			self.dynkin = ((3,1), ((1,2,3),))
		
		self.name = "%s %s" % (AlgName(self.latype), nmsfx)
		
		self.metric = LA_Metric(self.dynkin)
		self.invmet = MatrixTuple(MatrixInverse(self.metric))
		self.imetnum, self.imetden = SharedDen_Matrix(self.invmet)
		self.imetnum = MakeNP(self.imetnum)
		self.ctnmat = LA_Cartan(self.dynkin)
		self.invctn = MatrixTuple(MatrixInverse(self.ctnmat))
		self.ictnnum, self.ictnden = SharedDen_Matrix(self.invctn)
		self.ictnnum = MakeNP(self.ictnnum)
		
		self.posroots = LA_PositiveRoots(self.ctnmat);
		prsum = npzeros1(n)
		for rt in self.posroots:
			prsum += rt[0]
		self.posrootsum = prsum
	
	# Dimension of the algebra
	def dimension(self):
		return 2*len(self.posroots) + len(self.ctnmat)
	
	# Self-description:
	def __str__(self):
		return "Lie algebra: %s %s dim=%d" % \
			(self.latype, self.name, self.dimension())
	
	def __repr__(self):
		return "LieAlgebra(%d,%d)" % (self.latype[0], self.latype[1])

# Global cache of algebra values
LieAlgebraCache = {}

def GetLieAlgebra(latype):
	tptpl = tuple(latype)
	if tptpl not in LieAlgebraCache:
		LieAlgebraCache[tptpl] = LieAlgebra(tptpl)
	return LieAlgebraCache[tptpl]


# Root and weight interconversion
def WeightToRoot(la, wt):
	return np.dot(wt,la.ictnnum)
	
def RootToWeight(la, rt):
	return np.dot(rt,la.ctnmat)/la.ictnden

# Valid highest-weight vector?
def InvalidMaxWtsError(latype, maxwts):
	return TypeError("For type %s, invalid max weights: %s" % \
		(str(latype), str(maxwts)))

# Make NumPy array
# Replace CheckLARep with making a NumPy vector
# for a rep's highest weight or an orbit's dominant weight
def MakeMWVec(latype, maxwts_):
	la = GetLieAlgebra(latype)

	n = latype[1]
	if len(maxwts_) != n: raise InvalidMaxWtsError(latype, maxwts_)
	
	maxwts = MakeNP(map(int,maxwts_))
	for w in maxwts:
		if w < 0:
			raise InvalidMaxWtsError(latype, maxwts)
	
	return maxwts

# From an algebra and an irrep's max weights
# Since roots can have fractional values, multiply them
# by the shared denominator of the inverse of the Cartan matrix
# Returns tuple of values of (root, weight)
def RepRootVectors(latype, maxwts_):
	la = GetLieAlgebra(latype)
	maxwts = MakeMWVec(latype, maxwts_)
	
	# Initial root
	n = latype[1]
	ctnmat = la.ctnmat
	den = la.ictnden
	rtint = WeightToRoot(la,maxwts)
	RepRoots = [(rtint,maxwts)]
	
	# Test for whether we've found a root, and if so, return its index
	RRTest = {}
	rtkey = tuple(rtint)
	RRTest[rtkey] = 0
	
	# Find the next root until one cannot find any more
	WhichWay = [np.ones(n,dtype=np.bool_)]
	RtIndx = 0
	while RtIndx < len(RepRoots):
		root = RepRoots[RtIndx]
		ThisRootInt = root[0]
		ThisWeight = root[1]
		ctpd = list(root[1])
		# Advance in each direction, if possible
		for i in xrange(n):
			ThisWW = WhichWay[RtIndx]
			NewWW = np.array([j != i for j in xrange(n)], dtype=np.bool_)
			if ThisWW[i]:
				for j in xrange(1,ctpd[i]+1):
					# Calculate roots and weights in parallel,
					# to avoid repeated matrix.vector calculations
					NewRootInt = np.copy(ThisRootInt)
					NewRootInt[i] -= j*den
					NewWeight = np.copy(ThisWeight)
					NewWeight -= j*ctnmat[i]
					# Add the root if possible
					rtkey = tuple(NewRootInt)
					if rtkey in RRTest:
						WhichWay[RtIndx] = \
								np.logical_and(ThisWW, NewWW)
					else:
						RRTest[rtkey] = len(RepRoots)
						RepRoots.append((NewRootInt,NewWeight))
						WhichWay.append(NewWW)
		RtIndx += 1
		
	RepRoots.sort(RootWtSortFunc)
	return tuple(RepRoots)

# To find the root degeneracies, use Freudenthal's recurrence
# Needs only the rep-root values, not the weights
# The rep roots had been scaled up by the shared denominator of the inverse Cartan matrix
# The positive roots have been scaled likewise to integerize the calculations for speed
def RepRootDegens(latype, RepRoots):
	la = GetLieAlgebra(latype)
	
	den = la.ictnden
	posroots = [den*rt[0] for rt in la.posroots]
	posrootsum = den*la.posrootsum
	
	# Assumes RepRoots sorted upward by sum, like with RootWtSortFunc
	# Set up quick search for known ones
	rrix = {}
	for i,rt in enumerate(RepRoots):
		rtkey = tuple(rt)
		rrix[rtkey] = i
	# Top one
	maxrt = RepRoots[0]
	nrr = len(RepRoots)
	degens = npzeros1(nrr)
	degens[0] = 1
	# Next ones
	for k in xrange(1,nrr):
		rt = RepRoots[k]
		nx = 0
		for shtrt in posroots:
			bkrt = rt + shtrt
			brtpl = tuple(bkrt)
			while brtpl in rrix:
				i = rrix[brtpl]
				nx += degens[i]*np.dot(bkrt,np.dot(la.metric,shtrt))
				bkrt += shtrt
				brtpl = tuple(bkrt)
		mrpd = maxrt + rt + posrootsum
		mrm = maxrt - rt
		degens[k] = int(2*nx/np.dot(mrpd,np.dot(la.metric,mrm)))
	return tuple(degens)

# Returns list of (degeneracy, root value, weight value)
def RepRootVectorsDegens(latype, maxwts):
	reproots = RepRootVectors(latype, maxwts)
	degens = RepRootDegens(latype, [rt[0] for rt in reproots])
	return tuple([(degens[k],rt[0],rt[1]) for k,rt in enumerate(reproots)])

# Short version.
# Call with the Lie-algebra type and the max-weight vector
# It returns a tuple of
# (degeneracy, root vector * inverse-Cartan denominator, weight vector)
# The root vectors are scaled so as to integerize them for fast calculation
def GetRepDirect(latype, maxwts):
	return RepRootVectorsDegens(latype, maxwts)

# Total degeneracy of explicit representation
def TotalDegenOfExpRep(rep):
	return sum([rt[0] for rt in rep])

# Total degeneracy: Weyl's celebrated formula
def TotalDegen(latype, maxwts_):
	la = GetLieAlgebra(latype)
	maxwts = MakeMWVec(latype, maxwts_)
	
	posroots = [rt[0] for rt in la.posroots]
	twodelta = la.posrootsum
	maxrtsdel = 2*np.dot(maxwts,la.invctn) + twodelta
	n = Fraction(1,1)
	for rt in posroots:
		ra = np.dot(rt,la.metric)
		n *= Fraction(int(np.dot(ra,maxrtsdel)), int(np.dot(ra,twodelta)))
	return int(n)


# Weyl orbits

# From the orbit's dominant weight. Step down using the simple roots
# Dominant weight specified like an irrep's highest weight
def WeylOrbitForDomWt(latype, maxwts_):
	la = GetLieAlgebra(latype)
	maxwts = MakeMWVec(latype, maxwts_)
	n = latype[1]
	metric = la.metric
	
	# Initial root
	rtint = WeightToRoot(la,maxwts)
	oldrts = set([tuple(rtint)])
	rtlist = [rtint]
	
	# Next roots
	while len(oldrts) > 0:
		newrts = set()
		for rt in oldrts:
			mtrt = np.dot(metric,rt)
			for k in xrange(n):
				mr = mtrt[k]
				if mr > 0:
					nwrt = np.copy(rt)
					nwrt[k] -= (2*mr)/metric[k,k]
					newrts.add(tuple(nwrt))
		rtlist += map(MakeNP,newrts)
		oldrts = newrts
	
	return tuple([(rt,RootToWeight(la,rt)) for rt in rtlist])


def DomWtFromRoot(latype, root):
	la = GetLieAlgebra(latype)
	n = latype[1]
	metric = la.metric
	
	rt = np.copy(root)
	wt = RootToWeight(la,rt)
	while np.amin(wt) < 0:
		mtrt = np.dot(metric,rt)
		for k, mr in enumerate(mtrt):
			if k == 0:
				ix = k
				mtrtmin = mr
			else:
				if mr < mtrtmin:
					ix = k
					mtrtmin = mr
		mr = mtrt[ix]
		if mr > 0: break
		rt[ix] -= (2*mr)/metric[ix,ix]
		wt = RootToWeight(la,rt)
	
	return (rt, wt)

# Explicit expressions for Weyl orbits 

def accumulate(lst):
	return np.cumsum(MakeNP(lst),-1)

def MxWtToYD(maxwts):
	return MakeNP(list(reversed(accumulate(list(reversed(maxwts))))))

def AddSignsToRoot(root):
	rt0 = np.copy(root)
	n = len(rt0)
	sgnroot = [rt0]
	for k in xrange(n):
		if root[k] == 0: continue
		newsgnroot = []
		for rt in sgnroot:
			newsgnroot.append(rt)
			rtx = np.copy(rt)
			rtx[k] *= -1
			newsgnroot.append(rtx)
		sgnroot = newsgnroot
	# Select ones that are unique to within sorting;
	# assumes that the input is already sorted
	for k in xrange(n-1):
		if abs(rt0[k]) == abs(rt0[k+1]):
			newsgnroot = []
			for rt in sgnroot:
				if rt[k] >= rt[k+1]:
					newsgnroot.append(rt)
			sgnroot = newsgnroot
	
	return tuple(sgnroot)

# Implemented as an iterator
# To get list of them, do list(permutations(list to permute))
def permutations(lst):
	x = list(lst)
	n = len(x)
	x.sort()
	yield list(x)
	while True:
		WasFound = False
		for k in xrange(n-2,-1,-1):
			if x[k] < x[k+1]:
				ix = k
				WasFound = True
				break
		if not WasFound: break
		for l in xrange(n-1,ix,-1):
			if x[ix] < x[l]:
				iy = l
				break
		xs = x[ix]
		x[ix] = x[iy]
		x[iy] = xs
		x = x[:ix+1] + list(reversed(x[ix+1:]))
		yield list(x)

def flatten(lsts): return reduce(lambda a,b: a+b, map(list,lsts), [])

def PermuteRootsInList(rlst):
	return flatten(map(permutations,rlst))


# Has all but the exceptional algebras: G2, F4, E6, E7, E8

def WeylOrbitForDomWtExplicit(latype, maxwts_):
	la = GetLieAlgebra(latype)
	maxwts = MakeMWVec(latype, maxwts_)
	family,n = latype
	
	if family == 1:
		mwx = npzeros1(n+1)
		mwx[:n] = maxwts
		mwx = MxWtToYD(mwx)
		mwxtot = np.sum(mwx)
		mwx = (n+1)*mwx - mwxtot
		mwxl = accumulate(list(permutations(mwx)))
		orbrts = mwxl[:,:-1]
	
	elif family == 2:
		mwx = np.copy(maxwts)
		mwx[:-1] *= 2
		mwx = MxWtToYD(mwx)
		mwxl = AddSignsToRoot(mwx)
		orbrts = accumulate(PermuteRootsInList(mwxl))
	
	elif family == 3:
		mwx = np.copy(maxwts)
		mwx = MxWtToYD(mwx)
		mwxl = AddSignsToRoot(mwx)
		mwxl = accumulate(PermuteRootsInList(mwxl))
		mwxl[:,:-1] *= 2
		orbrts = mwxl
	
	elif family == 4:
		mwx = np.copy(maxwts)
		mwxm2 = mwx[-2]
		mwxm1 = mwx[-1]
		mwx[:-2] *= 2
		mwx[-2] = mwxm2 + mwxm1
		mwx[-1] = -mwxm2 + mwxm1
		mwx[:-1] = MxWtToYD(mwx[:-1])
		mwxl = AddSignsToRoot(mwx)
		if mwx[-1] != 0:
			mwxl = [rt for rt in mwxl if np.prod(rt)*mwx[-1] > 0]
		mwxl = MakeNP(PermuteRootsInList(mwxl))
		mlm2 = mwxl[:,-2]
		mlm1 = mwxl[:,-1]
		mwxl[:,-2] = mlm2 - mlm1
		mwxl[:,-1] = 2*mlm1
		mwxl = accumulate(mwxl)
		if (n % 2) == 0:
			mwxl[:,-2:] /= 2
		else:
			mwxl[:,:-2] *= 2
		orbrts = mwxl
	
	else:
		orbrts = []
	
	orbwts = RootToWeight(la,orbrts)
	orblen = len(orbrts)
	return tuple(((orbrts[k],orbwts[k]) for k in xrange(orblen)))


# Caching of Weyl orbits. Uses a global cache
WeylOrbitCache = {}

def GetOrbit(latype, maxwts):
	CacheKey = (tuple(latype), tuple(maxwts))
	if CacheKey not in WeylOrbitCache:
		if latype[0] <= 4:
			orb = WeylOrbitForDomWtExplicit(latype, maxwts)
		else:
			orb = WeylOrbitForDomWt(latype, maxwts)
		WeylOrbitCache[CacheKey] = orb
	return WeylOrbitCache[CacheKey]


# Statistics -- returns dom root, dom wt, length of orbit,
# and length of roots

def WeylOrbitStats(latype, maxwts):
	orb = GetOrbit(latype, maxwts)
	
	la = GetLieAlgebra(latype)
	maxroots = WeightToRoot(la, maxwts)
	return (MakeNP(maxwts), maxroots, len(orb), \
		np.dot(maxroots,np.dot(la.metric,maxroots)))

# Note: the orbit-stabilizer theorem constrains the number of members of
# Weyl-group orbits:
#
# (order of orbit of x in G) * (order of stabilizer of x in G) =
# (order of G) --
# the stabilizer subgroup of G for x has all the elements of G
# that keep x fixed


# Irreps from Weyl Orbits

# Dominant weights of Weyl orbits of a rep

def WeylDomWtRepRootVectors(latype, maxwts_):
	la = GetLieAlgebra(latype)
	maxwts = MakeMWVec(latype, maxwts_)
	
	# Initial root
	n = latype[1]
	den = la.ictnden
	rtint = WeightToRoot(la,maxwts)
	RepRoots = [(rtint,maxwts)]
	
	# Positive roots of the algebra
	posroots = [den*rt[0] for rt in la.posroots]
	
	# Test for whether we've found a root, and if so, return its index
	RRTest = {}
	rtkey = tuple(rtint)
	RRTest[rtkey] = 0
	
	# Find the next root until one cannot find any more
	RtIndx = 0
	while RtIndx < len(RepRoots):
		root = RepRoots[RtIndx]
		ThisRootInt = root[0]
		# Advance in each direction, if possible
		for psrt in posroots:
			NewRootInt = ThisRootInt - psrt
			NewWeight = RootToWeight(la,NewRootInt)
			if np.amin(NewWeight) < 0: continue
			rtkey = tuple(NewRootInt)
			if rtkey not in RRTest:
				RRTest[rtkey] = len(RepRoots)
				RepRoots.append((NewRootInt,NewWeight))
		RtIndx += 1
	
	RepRoots.sort(RootWtSortFunc)
	return tuple(RepRoots)

def WeylDomWtRepRootDegens(latype, RepRoots):
	la = GetLieAlgebra(latype)
	
	den = la.ictnden
	posroots = [den*rt[0] for rt in la.posroots]
	posrootsum = den*la.posrootsum
	
	# Assumes RepRoots sorted upward by sum, like with RootWtSortFunc
	# Set up quick search for known ones
	rrix = {}
	for i,rt in enumerate(RepRoots):
		rtkey = tuple(rt)
		rrix[rtkey] = i
	# Top one
	maxrt = RepRoots[0]
	nrr = len(RepRoots)
	degens = npzeros1(nrr)
	degens[0] = 1
	# Next ones
	for k in xrange(1,nrr):
		rt = RepRoots[k]
		nx = 0
		for shtrt in posroots:
			bkrt = rt + shtrt
			bkdmrw = DomWtFromRoot(latype, bkrt)
			brtpl = tuple(bkdmrw[0])
			while brtpl in rrix:
				i = rrix[brtpl]
				nx += degens[i]*np.dot(bkrt,np.dot(la.metric,shtrt))
				bkrt += shtrt
				bkdmrw = DomWtFromRoot(latype, bkrt)
				brtpl = tuple(bkdmrw[0])
		mrpd = maxrt + rt + posrootsum
		mrm = maxrt - rt
		degens[k] = int(2*nx/np.dot(mrpd,np.dot(la.metric,mrm)))
	return tuple(degens)

def WeylDomWtRepRootVectorsDegens(latype, maxwts):
	reproots = WeylDomWtRepRootVectors(latype, maxwts)
	degens = WeylDomWtRepRootDegens(latype, [rt[0] for rt in reproots])
	return tuple([(degens[k],rt[0],rt[1]) for k,rt in enumerate(reproots)])


# Caching of irreps' Weyl orbits. Uses a global cache
RepresentationWeylOrbitCache = {}

def GetRepOrbits(latype, maxwts):
	CacheKey = (tuple(latype), tuple(maxwts))
	if CacheKey not in RepresentationWeylOrbitCache:
		RepresentationWeylOrbitCache[CacheKey] = \
			WeylDomWtRepRootVectorsDegens(latype, maxwts)
	return RepresentationWeylOrbitCache[CacheKey]

# Orbits input as list of (multiplicity, dominant root, dominant weight)

def ExpandWeylOrbits(latype, orbits):
	exporbs = [[(orb[0],wr[0],wr[1]) for wr in GetOrbit(latype,orb[2])] \
		for orb in orbits]
	return tuple(flatten(exporbs))

# Expand the rep with the Weyl orbits

def GetRepOrbitsExpanded(latype, maxwts):
	return ExpandWeylOrbits(latype, GetRepOrbits(latype, maxwts))


# Caching of representations. Uses a global cache 
RepresentationCache = {}

def GetRep(latype, maxwts):
	CacheKey = (tuple(latype), tuple(maxwts))
	if CacheKey not in RepresentationCache:
		RepresentationCache[CacheKey] = GetRepOrbitsExpanded(latype, maxwts)
	return RepresentationCache[CacheKey]

# Compare to GetRepDirect(latype, maxwts) -- the old code


# Properties of irreps

# Finds the max-weight vector of the conjugate irrep, 
# the one with the root - vector/weight - vector signs reversed
def RepConjugate(latype, maxwts_):
	la = GetLieAlgebra(latype)
	maxwts = MakeMWVec(latype, maxwts_)
	family,n = latype
	
	wt = list(maxwts)
	if family == 1:
		# A(n)
		wt.reverse()
	elif family == 4 and (n % 2) == 1:
		# D(n) for odd n; for even n, self-conjugate
		wt0 = wt[:-2]
		wt1 = wt[-2:]
		wt1.reverse()
		wt = wt0 + wt1
	elif family == 5 and n == 6:
		# E6
		wt0 = wt[:-1]
		wt1 = wt[-1:]
		wt0.reverse()
		wt = wt0 + wt1
	
	return tuple(wt)

def RepIsSelfConjugate(latype, maxwts):
	return tuple(maxwts) == tuple(RepConjugate(latype,maxwts))

# The height of a rep is max - min of  all the root - vector component sums
# Multiplier for max weight; height is linear in the max weight
#
# If a rep is self-conjugate,
# then if its height is even, it is real,
# while if its height is odd, it is pseudoreal
def RepHeightMult(latype):
	la = GetLieAlgebra(latype)
	family,n = latype
	
	if family == 1: # A(n)
		mult = [(k+1)*(n-k) for k in xrange(n)]
	elif family == 2: # B(n)
		mult = [(k+1)*(2*n-k) for k in xrange(n-1)] + [n*(n+1)/2]
	elif family == 3: # C(n)
		mult = [(k+1)*(2*n-k-1) for k in xrange(n)]
	elif family == 4: # D(n)
		mult = [(k+1)*(2*n-k-2) for k in xrange(n-2)] + [n*(n-1)/2, n*(n-1)/2]
	elif family == 5: # E6, E7, E8
		if n == 6:
			mult = [16, 30, 42, 30, 16, 22]
		elif n == 7:
			mult = [34, 66, 96, 75, 52, 27, 49]
		elif n == 8:
			mult = [92, 182, 270, 220, 168, 114, 58, 136]
		else:
			mult = []
	elif family == 6: # F4
		if n == 4:
			mult = [22, 42, 30, 16];
		else:
			mult = []
	elif family == 7: # G2
		if n == 2:
			mult = [10, 6]
		else:
			mult = []
	else:
		mult = []
	
	return MakeNP(mult)

def RepHeight(latype, maxwts): return int(np.dot(RepHeightMult(latype),MakeNP(maxwts)))

# The reality of a representation:
# Real (self-conjugate, even height): 0
# Pseudoreal (self-conjugate, odd height) = 1
# Complex (unequal to conjugate) = -1
# Args: algebra type, max-weights vector
def RepReality(latype, maxwts_):
	la = GetLieAlgebra(latype)
	maxwts = MakeMWVec(latype, maxwts_)
	family,n = latype
	if not RepIsSelfConjugate(latype,maxwts): return -1
	
	if family == 1 and (n % 2) == 1:
		# A(n), n odd
		psr = np.dot(maxwts,np.arange(1,n+1))
	elif family == 2 and ((n-1) % 4) < 2:
		# B(n), n = 4k+1,4k+2
		if ((n-1) % 4) < 2:
			psr = maxwts[-1]
		else:
			psr = 0
	elif family == 3:
		# C(n)
		psr = np.dot(maxwts,np.arange(1,n+1))
	elif family == 4:
		# D(n), n = 4k+2,4k+3
		if ((n-2) % 4) < 2:
			psr = maxwts[-2] + maxwts[-1]
		else:
			psr = 0
	elif family == 5 and n == 7:
		# E7
		psr = maxwts[4] + maxwts[6] + maxwts[7],
	else:
		# G2, F4, E6, E8
		psr = 0
	
	return (psr % 2)

# Weight congruency classes / conserved quantities
#
# Conserved quantities in weights in a rep : 
# for a vector c such that (Cartan matrix) . c = 0 mod p, 
# then all weights w will have w.c = same mod p. 
# Rep products will have highest weights related by
# (w - w1 - w2).c mod p = 0
# Arg: LA type
def RepConservModMult(latype):
	la = GetLieAlgebra(latype)
	family,n = latype
	
	if family == 1: # A(n)
		cvms = ((n+1, tuple([(k+1) for k in xrange(n)])),)
	elif family == 2: # B(n)
		cvms = ((2, tuple(zeros1(n-1) + [1])),)
	elif family == 3: # C(n)
		cvms = ((2, tuple([(k+1) % 2 for k in xrange(n)])),)
	elif family == 4: # D(n)
		cvms0 = (2, tuple(zeros1(n-2) + [1,1]))
		if (n % 2) == 0:
			cvms = (cvms0, \
				(2, tuple([(k+1) % 2 for k in xrange(n-2)] + [1,0])), \
				(2, tuple([(k+1) % 2 for k in xrange(n-2)] + [0,1])) \
				)
		else:
			cvms = (cvms0, \
				(4, tuple([2*((k+1) % 2) for k in xrange(n-2)] + [1,3])) \
				)
	elif family == 5:
		if n == 6: # E6
			cvms = ((3, (1,2,0,1,2,0)),)
		elif n == 7: # E7
			cvms = ((2, (0,0,0,1,0,1,1)),)
		else: # E8
			cvms = ()
	else: # F4, G2
		cvms = ()

	return cvms

# The modulus values, the values that the raw conserved values are divided by
def RepConservModuli(latype):
	return tuple([cvm[0] for cvm in RepConservModMult(latype)])

# List of (modulus, conserved qty)
# Args: algebra type, max-weights vector
def RepConserv(latype, maxwts):
	return tuple([(cvm[0], mul_vv(cvm[1],maxwts) % cvm[0]) \
		for cvm in RepConservModMult(latype)])

# Overall list of properties, as a dict:
# type: LA type
# maxwts: max-weights vector
# mwconjg: its conjugate
# isselfconjg: whether those two are equal
# dimension: total degeneracy or dimension of the irrep
# height: max - min of sum(values in each root vector)
# reality: real(0), pseudoreal(1), or complex(-1)
# conserved: conserved quantities
# Args: algebra type, max-weights vector
def RepProperties(latype, maxwts):
	props = {}
	props["type"] = latype
	mwts = tuple(maxwts)
	props["maxwts"] = mwts
	mwcjg = RepConjugate(latype,mwts)
	props["mwconjg"] = mwcjg
	sc = mwts == mwcjg
	props["isselfconjg"] = sc
	ht = RepHeight(latype,mwts)
	props["dimension"] = TotalDegen(latype,mwts)
	props["height"] = ht
	props["reality"] = (ht % 2) if sc else -1
	props["conserved"] = RepConserv(latype,mwts)
	return props

# Casimir invariant (2nd-order one)
# Args: algebra type, max-weights vector
def CasimirInvariant(latype, maxwts_):
	la = GetLieAlgebra(latype)
	maxwts = MakeMWVec(latype, maxwts_)
	
	posroots = [rtx[0] for rtx in la.posroots]
	maxroots = np.dot(maxwts,la.invctn)
	return np.dot(maxroots,np.dot(la.metric,(maxroots+la.posrootsum)))

# Representation index
# Args: algebra type, max-weights vector
def RepIndex(latype, maxwts):
	la = GetLieAlgebra(latype)

	td = TotalDegen(latype,maxwts)
	ci = CasimirInvariant(latype,maxwts)
	adjdim = la.dimension()
	return Fraction(int(td*ci),adjdim)


# Reducible representations: express as sum of irreps,
# a list of irrep highest weights

def RepSortFunc(rt1,rt2): return cmp(tuple(rt2[1]),tuple(rt1[1]))

def TrimRepZeros(rep):
	return tuple([tuple(r) for r in sorted(rep,RepSortFunc) if r[0] != 0])

def GetRepListGeneral(latype, mwlist, repfunc):
	# Find the combined representation.
	# Would like degens, roots, and weights, since we'll need all three.
	# Use index function rppwrix to speed up lookups in the power-rep list;
	# index with the weight vector. Also pull out
	# the root multiplicities for convenience
	rpix = {}
	for maxwts in mwlist:
		indrep = repfunc(latype,maxwts)
		for rprt in indrep:
			rrkey = tuple(rprt[2])
			if rrkey not in rpix:
				rpix[rrkey] = list(rprt)
			else:
				rpix[rrkey][0] += rprt[0]
	
	return TrimRepZeros(rpix.values())

def GetRepCntdListGeneral(latype, mwclist, repfunc):
	# Find the combined representation for a counted list: (count, maxwts).
	# Would like degens, roots, and weights, since we'll need all three.
	# Use index function rppwrix to speed up lookups in the power-rep list;
	# index with the weight vector. Also pull out
	# the root multiplicities for convenience
	rpix = {}
	for cnt,maxwts in mwclist:
		indrep = repfunc(latype,maxwts)
		indrep = [tuple([cnt*ir[0]] + list(ir[1:])) for ir in indrep]
		for rprt in indrep:
			rrkey = tuple(rprt[2])
			if rrkey not in rpix:
				rpix[rrkey] = list(rprt)
			else:
				rpix[rrkey][0] += rprt[0]
	
	return TrimRepZeros(rpix.values())

def GetRepList(latype,mwlist):
	return GetRepListGeneral(latype, mwlist, GetRep)

def GetRepCntdList(latype,mwclist):
	return GetRepCntdListGeneral(latype, mwclist, GetRep)

# Takes a rep and turns the roots into appropriate fractional versions

def mulroot(rootscale, rootvector):
	return tuple([rootscale*int(r) for r in rootvector])

def GetTrueRep(latype, rep):
	la = GetLieAlgebra(latype)
	rootscale = Fraction(1,la.ictnden)
	return tuple([(r[0], mulroot(rootscale,r[1]), r[2]) for r in rep])

def mulroots(rootscales, rootvectors):
	nrs = len(rootscales)
	return tuple([mulroot(rs,rootvectors[k]) for k,rs in enumerate(rootscales)] + \
		list(rootvectors[nrs:]))

def GetTrueAlgProdRep(latlist, rep):
	rootscales = []
	for latype in latlist:
		la = GetLieAlgebra(latype)
		rootscale = Fraction(1,la.ictnden)
		rootscales.append(rootscale)
	return tuple([(r[0], mulroots(rootscales, r[1]), r[2]) for r in rep])

# Irrep of product of Lie algebras: outer product of irreps of each of them.
# Follows the format of the individual irreps:
#
# list of
# (degen, list of root-vector values, list of max-weight-vector values)
#
# U(1) factors are implemented by adding integers or rational numbers
# after the max-weight vectors in maxwts.
# Every root in a rep will have these U(1) factors.
# It's possible for a rep to contain only U(1) factors.

def GetAlgProdRepGeneral(latlist, maxwts, repfunc):
	rep = [[1, [], []]]
	for i,latype in enumerate(latlist):
		indmwts = maxwts[i]
		indrep = repfunc(latype,indmwts)
		newrep = []
		for ptrt in rep:
			for rt in indrep:
				nwrt = [ptrt[0]*rt[0], ptrt[1]+[rt[1]], ptrt[2]+[rt[2]]]
				newrep.append(nwrt)
		rep = newrep
	u1s = maxwts[len(latlist):]
	for r in rep:
		for k in xrange(1,2+1):
			r[k] += u1s
			r[k] = tuple(r[k])
	
	return MatrixTuple(rep)

def GetAlgProdRep(latlist, maxwts):
	return GetAlgProdRepGeneral(latlist, maxwts, GetRep)

def GetAlgProdRepOrbits(latlist, maxwts):
	return GetAlgProdRepGeneral(latlist, maxwts, GetRepOrbits)

def GetAlgProdRepList(latype,mwlist):
	return GetRepListGeneral(latype, mwlist, GetAlgProdRep)

def GetAlgProdRepCntdList(latype,mwclist):
	return GetRepCntdListGeneral(latype, mwclist, GetAlgProdRep)

# Select extended: type is "sngl" -- single irrep, 
# "list" -- list of irreps, "cntd" -- counted list of irreps

def GetRepXtnd(latype, rptype, maxwts):
	if rptype == "sngl": return GetRep(latype,maxwts)
	elif rptype == "list": return GetRepList(latype,maxwts)
	elif rptype == "cntd": return GetRepCntdList(latype,maxwts)

def GetAlgProdRepXtnd(latlist, rptype, maxwts):
	if rptype == "sngl": return GetAlgProdRep(latlist,maxwts)
	elif rptype == "list": return GetAlgProdRepList(latlist,maxwts)
	elif rptype == "cntd": return GetAlgProdRepCntdList(latlist,maxwts)

# Extensions. Type is "sngl" -- single irrep, "list" -- list of irreps,
# "cntd" -- counted list of irreps. The algebra products are
# (la1, ..., lan), with the corresponding irreps (wt1, ..., wtn, u1, ... um)
# where the la's are the individual algebras, the wt's the individual weights,
# and the u's U(1) factors.

def TotalDegenXtnd(latype, rptype, maxwts):
	if rptype == "sngl":
		return TotalDegen(latype,maxwts)
	elif rptype == "list":
		n = 0
		for mw in maxwts:
			n += TotalDegen(latype,mw)
		return n
	elif rptype == "cntd":
		n = 0
		for mw in maxwts:
			n += mw[0]*TotalDegen(latype,mw[1])
		return n

def AlgProdTotalDegen(latlist, maxwts):
	n = 1
	for i,lat in enumerate(latlist):
		n *= TotalDegen(lat,maxwts[i])
	return n

def AlgProdTotalDegenXtnd(latlist, rptype, maxwts):
	if rptype == "sngl":
		return AlgProdTotalDegen(latlist,maxwts)
	elif rptype == "list":
		n = 0
		for mw in maxwts:
			n += AlgProdTotalDegen(latlist,mw)
		return n
	elif rptype == "cntd":
		n = 0
		for mw in maxwts:
			n += mw[0]*AlgProdTotalDegen(latlist,mw[1])
		return n

# Much simpler than Mathematica
# Indexed rep is a dict from weight vector
# to list of multiplicity/degeneracy, root vector, weight vector
# Not a tuple, since we need to adjust the mult

def DoRepIndexing(rep):
	rpix = {}
	for r in rep: rpix[tuple(r[2])] = list(r)
	return rpix

def AddEntryToIndexedRep(rpix, rpen, mult=1):
	rrkey = tuple((tuple(r) if type(r) == np.ndarray else r for r in rpen[2]))
	if rrkey not in rpix:
		# Add a new entry
		rpls = list(rpen)
		rpls[0] *= mult
		rpix[rrkey] = rpls
	else:
		# Adjust the multiplicity of an existing entry
		rpix[rrkey][0] += mult*rpen[0]

def IndexedRepToList(rpix):
	rep = rpix.values()
	return tuple([tuple(r) for r in sorted(rep,RepSortFunc) if r[0] != 0])

# Extracts irreps from an indexed rep, altering it as it goes.
# Args: algebra, indexed rep, maximum-root function, irrep-finding function
# Alters rpix as it goes.
def ExtractRepIrrepsGeneral(la, rpix, maxrootfunc, repfunc):
	irreplist = []
	while True:
		# Clear out zero ones first
		wts = rpix.keys()
		for wt in wts:
			if rpix[wt][0] == 0: del rpix[wt]
		if len(rpix) <= 0: break
		
		maxroot = None
		for wt in rpix:
			maxroot = rpix[wt] if maxroot == None \
				else maxrootfunc(la,maxroot,rpix[wt])
		
		nwt = maxroot[0]
		wts = maxroot[2]
		irreplist.append((nwt,wts))
		rep = repfunc(la,wts)
		
		for rpen in rep:
			AddEntryToIndexedRep(rpix,rpen,-nwt)
	
	return tuple(irreplist)

def MaxRepRoot(la,en,enx):
	rt = en[1]; rtx = enx[1]
	maxdiff = sum(sub_vv(rtx,rt))
	return enx if maxdiff > 0 else en

def ExtractRepIrreps(la,rpix):
	return ExtractRepIrrepsGeneral(la,rpix,MaxRepRoot,GetRep)

def ExtractRepOrbitIrreps(la,rpix):
	return ExtractRepIrrepsGeneral(la,rpix,MaxRepRoot,GetRepOrbits)

def MaxAlgProdRepRoot(la,en,enx):
	rt = en[1]; rtx = enx[1]
	lalen = len(la)
	maxdiff = sum([sum(sub_vv(rtx[i],rt[i])) for i in xrange(lalen)])
	rtlen = len(rt)
	if rtlen > lalen:
		maxdiff += sum(sub_vv(rtx[lalen:],rt[lalen:]))
	return enx if maxdiff > 0 else en

def ExtractAlgProdRepIrreps(la,rpix):
	return ExtractRepIrrepsGeneral(la,rpix,MaxAlgProdRepRoot,GetAlgProdRep)

def ExtractAlgProdRepOrbitIrreps(la,rpix):
	return ExtractRepIrrepsGeneral(la,rpix,MaxAlgProdRepRoot,GetAlgProdRepOrbits)


# Find product reps as counted lists of irreps, lists of (degen, maxwts)

def AddRootRecords(la,rt1,rt2):
	degen = rt1[0]*rt2[0]
	root = rt1[1] + rt2[1]
	weight = rt1[2] + rt2[2]
	return [degen, root, weight]

def AddRootVectors(la,v1,v2):
	lalen = len(la)
	vsum = [v1[i]+v2[i] for i in xrange(lalen)]
	vsum += add_vv(v1[lalen:],v2[lalen:])
	return tuple(vsum)

def AddAlgProdRootRecords(la,rt1,rt2):
	degen = rt1[0]*rt2[0]
	root = AddRootVectors(la,rt1[1],rt2[1])
	weight = AddRootVectors(la,rt1[2],rt2[2])
	return [degen, root, weight]

# UE functions - whether or not to use an entry. The args are
# Length of the SS-algebra-product part,
# Weight vector -- SS-algebra weights or list of weights + U(1) factors.
#
# The length argument is ignored for a single weight vector or for
# using all the weights (UEAll). The "Orbits" functions are for using
# only the dominant weights of Weyl orbits, a speedup that also saves
# memory.

def UEAll(nsa, wts): return True

def UEOrbits(nsa, wts): return min(wts) >= 0

def UEAlgProdOrbits(nsa, wts):
	return all([min(wt) >= 0 for wt in wts[:nsa]]) \
		if nsa > 0 else True


# From an algebra and a pair of irreps
def DecomposeRepProductRoots(latype, rep1, rep2, addfunc, maxfunc, \
	repfunc, useentry, nsa):
	rpix = {}
	
	# The product rep
	for rt1 in rep1:
		for rt2 in rep2:
			rtc = addfunc(latype,rt1,rt2)
			if useentry(nsa, rtc[2]):
				AddEntryToIndexedRep(rpix,rtc)
	
	return ExtractRepIrrepsGeneral(latype, rpix, maxfunc, repfunc)

def DecomposeRepProductXtnd(latype, rptype1, maxwts1, rptype2, maxwts2):
	return DecomposeRepProductRoots(latype, \
		GetRepXtnd(latype, rptype1, maxwts1), \
		GetRepXtnd(latype, rptype2, maxwts2), \
		AddRootRecords, MaxRepRoot, \
		GetRepOrbits, UEOrbits, 0)

# Weyl - orbit dom wts -- all entries : GetRep, UEAll, 0

def DecomposeRepProduct(latype, maxwts1, maxwts2):
	return DecomposeRepProductXtnd(latype, "sngl", maxwts1, "sngl", maxwts2)

def ExpandXtndWeightList(rptype, maxwts):
	if rptype == "sngl":
		return ((1, maxwts),)
	elif rptype == "list":
		return tuple([(1,mw) for mw in maxwts])
	elif rptype == "cntd":
		return maxwts

def DecomposeRepProductXtndAlt(latype, rptype1, maxwts1, rptype2, maxwts2):
	mwcl1 = ExpandXtndWeightList(rptype1, maxwts1)
	mwcl2 = ExpandXtndWeightList(rptype2, maxwts2)
	rsno = {}
	for mw1 in mwcl1:
		for mw2 in mwcl2:
			mres = DecomposeRepProduct(latype, mw1[1], mw2[1])
			for mr in mres:
				nwdg = mw1[0]*mw2[0]*mr[0]
				nwwt = mr[1]
				if nwwt not in rsno:
					rsno[nwwt] = 0
				rsno[nwwt] += nwdg
	rsres = [(rsno[wt],wt) for wt in rsno if rsno[wt] != 0]
	rsres.sort(RepSortFunc)
	return tuple(rsres)

def DecomposeAlgProdRepProductXtnd(latlist, rptype1, maxwts1, rptype2, maxwts2):
	return DecomposeRepProductRoots(latlist, \
		GetAlgProdRepXtnd(latlist, rptype1, maxwts1), \
		GetAlgProdRepXtnd(latlist, rptype2, maxwts2), \
		AddAlgProdRootRecords, MaxAlgProdRepRoot, \
		GetAlgProdRepOrbits, UEAlgProdOrbits, len(latlist))

# Weyl - orbit dom wts -- all entries : GetAlgProdRep, UEAll, 0

def DecomposeAlgProdRepProduct(latlist, maxwts1, maxwts2):
	return DecomposeAlgProdRepProductXtnd(latlist, "sngl", maxwts1, "sngl", maxwts2)

def DecomposeAlgProdRepProductXtndAlt(latlist, rptype1, maxwt1, rptype2, maxwt2):
	mwcl1 = ExpandXtndWeightList(rptype1, maxwt1)
	mwcl2 = ExpandXtndWeightList(rptype2, maxwt2)
	ltlen = len(latlist)
	rsno = {}
	for mw1 in mwcl1:
		for mw2 in mwcl2:
			mres = [[1, []]]
			for i,latype in enumerate(latlist):
				mw1i = mw1[1][i]
				mw2i = mw2[1][i]
				mwri = DecomposeRepProduct(latype,mw1i,mw2i)
				newmres = []
				for mr in mres:
					for mi in mwri:
						npd = [mr[0]*mi[0], mr[1] + [mi[1]]]
						newmres.append(npd)
				mres = newmres
			u1s1 = mw1[1][ltlen:]
			u1s2 = mw2[1][ltlen:]
			u1s = add_vv(u1s1,u1s2)
			mres = [(mr[0], tuple(mr[1]+u1s)) for mr in mres]
			for mr in mres:
				nwdg = mw1[0]*mw2[0]*mr[0]
				nwwt = mr[1]
				if nwwt not in rsno:
					rsno[nwwt] = 0
				rsno[nwwt] += nwdg
	rsres = [(rsno[wt],wt) for wt in rsno if rsno[wt] != 0]
	rsres.sort(RepSortFunc)
	return tuple(rsres)

def DecomposeAlgProdRepProductAlt(latlist, maxwt1, maxwt2):
	return DecomposeAlgProdRepProductXtndAlt(latlist, "sngl", maxwt1, "sngl", maxwt2)

# Like the previous rep-product decomposer,
# but it uses a list of irrep max-weight vectors instead of two of them
def DecomposeRepProdList(latype, maxwtlist):
	res = ((1,tuple(zeros1(latype[1]))),)
	for maxwts in maxwtlist:
		res = DecomposeRepProductXtnd(latype, "cntd", res, "sngl", maxwts)
	return res


# Irrep Powers

# This is separate from Irrep Products, because it has a lot of stuff
# that has to be calculated for it 

# Tensor-product A(n)/SU(n+1) reps calculated with Young diagrams.
# B(n)/SO(2n+1), D(n)/SO(2n) have traces factored out
# and C(n)/Sp(2n) have sort-of-traces with the symplectic metrix (J = - J^T)

# This product has additional precalculated info --
# Kostka matrices (orbits to irreps) and nesting of Young diagrams
# inside of other Young diagrams

# Data format of irrep-power list members:
# Young=diagram info, Kostka matrix, rep-nesting matrix

# Young-diagram info:
# List of multiplicity, Young diagram, length of YD,
# multinomial multiplicity factor for it

# Rep nesting: each rep's numbers becomes a rep and one iterates over them.
# Row: index of starting rep, column: index of nested reps concatenated
# At each row and column, a list of entries:
# Indices of subreps, multiplicities of subreps,
# multinomial multiplicity factor

# Reset by doing: TensorPowersYDX = []

TensorPowersYDX = []

def YDTensAddVector(ylst):
	ylnl = []
	ylnew = list(ylst)
	ylnew[0] += 1
	ylnl.append(tuple(ylnew))
	for i in xrange(len(ylst)-1):
		if ylst[i] > ylst[i+1]:
			ylnew = list(ylst)
			ylnew[i+1] += 1
			ylnl.append(tuple(ylnew))
	ylnew = list(ylst)
	ylnew.append(1)
	ylnl.append(tuple(ylnew))
	return tuple(ylnl)

def memcount(lst):
	mc = {}
	for mem in lst:
		if mem not in mc: mc[mem] = 0
		mc[mem] += 1
	return mc

def prodlist(lst): return reduce(lambda x,y: x*y, lst, 1)

def factorial(n):
	if n < 0: return 0
	return prodlist(xrange(1,int(n)+1))

def multnomfac(dgrm):
	cnts = memcount(dgrm).values()
	return int(factorial(sum(cnts))/prodlist(map(factorial,cnts)))

# Find next tensor-product diagram list from a tensor-product list of
# (multiplicity, contracted Young diagram)
# For inputs in a standard ordering of same-box-number Young diagrams,
# it produces outputs in that order.
def YDTensSetAddVector(dgrmlist):
	nxdgs = []
	dgcnt = {}
	for num,dg,ln,mnf in dgrmlist:
		ndgs = YDTensAddVector(dg)
		for ndg in ndgs:
			if ndg not in dgcnt:
				nxdgs.append(ndg)
				dgcnt[ndg] = 0
			dgcnt[ndg] += num
	nxdgmlst = [(dgcnt[ndg],ndg,len(ndg),multnomfac(ndg)) for ndg in nxdgs]
	return tuple(nxdgmlst)


def MakeYDSetIndexing(dgrms):
	dgix = {}
	for i,dg in enumerate(dgrms):
		dgix[dg] = i
	return dgix

def YDKostkaMatrix(dgrms):
	# Index the diagrams
	dgix = MakeYDSetIndexing(dgrms)
	# Find "C" lookback matrix and "g" denominator factor
	# Multiply by 2 to avoid fractions
	dglen = len(dgrms)
	lookback = zeros2(dglen,dglen)
	denfac = zeros1(dglen)
	for i,dg in enumerate(dgrms):
		dgl = len(dg)
		# Do lookback; the += takes care of the multiplicities automatically.
		for j in xrange(dgl-1):
			for jx in xrange(j+1,dgl):
				for k in xrange(1,dg[jx]+1):
					dgx = list(dg)
					dgx[j] += k
					dgx[jx] -= k
					dgd = dgx[j] - dgx[jx]
					dgx.sort()
					dgx.reverse()
					if dgx[-1] == 0: del dgx[-1]
					dgx = tuple(dgx)
					ix = dgix[dgx]
					lookback[i][ix] += 2*dgd
		# Denominator factor
		dfsum = 0
		for j,d in enumerate(dg):
			dfsum += d*(d - 2*(j+1))
		denfac[i] = dfsum
	kostka = zeros2(dglen,dglen)
	for i in xrange(dglen):
		kkrow = kostka[i]
		kkrow[i] = 1
		for j in xrange(i+1,dglen):
			kknum = 0
			lbrow = lookback[j]
			for k in xrange(i,j):
				kknum += lbrow[k]*kkrow[k]
			if kknum != 0:
				kkrow[j] = kknum/(denfac[i] - denfac[j])
	return MatrixTuple(kostka)

# Imitation of Mathematica Tuples function for a list of lists:
def mlstsel(mltlst):
	res = [[]]
	for lst in mltlst:
		newres = []
		for rsmem in res:
			for mem in lst:
				newres.append(rsmem + [mem])
		res = newres
	return res

def YDNesting(dgrmlist):
	dgrms = [dgl[1] for dgl in dgrmlist]
	dlens = [dgl[2] for dgl in dgrmlist]
	dmnfs = [dgl[3] for dgl in dgrmlist]
	dglen = len(dgrms)
	dgix = MakeYDSetIndexing(dgrms)
	nesting = [[[] for j in xrange(dglen)] for i in xrange(dglen)]
	for i,dgrm in enumerate(dgrms):
		subdgrms = [TensorPowersYDX[dm-1].yds() for dm in dgrm] \
			if len(dgrm) > 1 else [dgrms]
		sublens = [TensorPowersYDX[dm-1].lens() for dm in dgrm] \
			if len(dgrm) > 1 else [dlens]
		submnfs = [TensorPowersYDX[dm-1].mnfs() for dm in dgrm] \
			if len(dgrm) > 1 else [dmnfs]
		dgxp = mlstsel([range(len(dm)) for dm in subdgrms])
		for dgxpm in dgxp:
			insdg = []
			for l,dgxpmx in enumerate(dgxpm):
				insdg += subdgrms[l][dgxpmx]
			insdg.sort()
			insdg.reverse()
			insdg = tuple(insdg)
			j = dgix[insdg]
			lns = tuple([sublens[l][dgxpmx] for l,dgxpmx in enumerate(dgxpm)])
			mnf = prodlist([submnfs[l][dgxpmx] for l,dgxpmx in enumerate(dgxpm)])
			dxplst = (tuple(dgxpm), lns, mnf)
			nesting[i][j].append(dxplst)
	nesting = tuple([MatrixTuple(nx) for nx in nesting])
	return nesting


class TnsrPwrYDX_Cls:
	def __init__(self,ydcnts_):
		self.ydcnts = ydcnts_
		self.kostka = YDKostkaMatrix(self.yds())
		self.nesting = YDNesting(self.ydcnts)
	
	# "Length" value for it
	def __len__(self):
		return len(self.ydcnts)
	
	# Multiplicities:
	def cnts(self):
		return tuple((yc[0] for yc in self.ydcnts))
	
	# Representation Diagrams:
	def yds(self):
		return tuple((yc[1] for yc in self.ydcnts))
	
	# Their lengths:
	def lens(self):
		return tuple((yc[2] for yc in self.ydcnts))
	
	# Their multinomial factors:
	def mnfs(self):
		return tuple((yc[3] for yc in self.ydcnts))
	
	# Display
	def __repr__(self):
		return "TnsrPwrYDX_Cls(%s)" % repr(self.ydcnts)


def GetTensorPowerYDX(n):
	if n <= 0:
		return TnsrPwrYDX_Cls(())
	tplen = len(TensorPowersYDX)
	if tplen < 1:
		dgrm = (1,)
		dgrmlist = ((1,dgrm,len(dgrm),multnomfac(dgrm)),)
		TensorPowersYDX.append(TnsrPwrYDX_Cls(dgrmlist))
		tplen = len(TensorPowersYDX)
	
	for k in xrange(tplen+1,n+1):
		tnsp = TensorPowersYDX[-1]
		dgrmlist = tnsp.ydcnts
		dgrmlist = YDTensSetAddVector(dgrmlist)
		TensorPowersYDX.append(TnsrPwrYDX_Cls(dgrmlist))
	
	return TensorPowersYDX[n-1]

# Args: tensor-power-ydx value, indexer of rep/digram,
# list of indexes, multiplicity list (gets indexed).
# Returns: list of (multiplicity,index), multiplicities
# for product-rep calculation

def TensorPowerRootMultiplicities(TensPwrYDX,dgix,idxs,mults):
	mcdict = memcount(idxs)
	mc = zip(mcdict.values(),mcdict.keys())
	mc.sort(lambda a,b: cmp(b[0],a[0]))
	mc = tuple(mc)
	repyd = tuple([mcm[0] for mcm in mc])
	ixms = tuple([mults[mcm[1]] for mcm in mc])
	mspwrs = [[im] for im in ixms]
	kostka = TensPwrYDX.kostka
	nesting = TensPwrYDX.nesting[dgix[repyd]]
	mpres = []
	for nstmem in nesting:
		mprsind = 0
		for nstmmi in nstmem:
			pwrs = nstmmi[1]
			mncf = nstmmi[2]
			mprsii = mncf
			for k,pwr in enumerate(pwrs):
				mspwr = mspwrs[k]
				mspl = len(mspwr)
				if pwr > mspl:
					for l in xrange(mspl+1,pwr+1):
						mspwr.append(mspwr[-1]*(ixms[k]-l+1)/l)
				mprsii *= mspwr[pwr-1]
			mprsind += mprsii
		mpres.append(mprsind)
	return (mc, tuple(mul_mv(kostka,mpres)))

# Generator of index sets. The indexes are semi-ascending.
# Python generator - use where one would expect an iterator
# Generator expressions: () instead of list-comprehension []
# Can use function's ()

def MakeIndexVectors(vlen,ntot):
	ixvec = zeros1(vlen)
	while True:
		yield ixvec # Returned values
		didadvance = False
		for k in xrange(vlen-1,-1,-1):
			ixvec[k] += 1
			if ixvec[k] < ntot:
				ixvec[k+1:] = (vlen-k-1)*[ixvec[k]]
				didadvance = True
				break
		if not didadvance: break

# This is for decomposing a representation decomposed as roots,
# for the purpose of handling reducible
# as well as irreducible representations.

def DecomposeRepRootsPower(lax, rep, pwr, \
		mulfunc, muladdtofunc, tuplefunc, maxfunc, \
		repfunc, useentry, nsa):
	
	# Zero : scalar rep
	if pwr <= 0:
		return (((1,tuplefunc(lax,mulfunc(lax,0,rep[0][2]))),),)
	
	# Get tensor-product data for power 
	tpydx = GetTensorPowerYDX(pwr)
	dgrms = tpydx.yds()
	dglen = len(tpydx)
	dgix = MakeYDSetIndexing(dgrms)
		
	# Use index function rppwrix to speed up lookups in the power-rep list;
	# index with the weight vector.
	# Also pull out the root multiplicities for convenience
	rtdgns = [r[0] for r in rep]
	rppwrdcmp = [{} for k in xrange(dglen)]
	
	for ixvec in MakeIndexVectors(pwr,len(rtdgns)):
		tpdgrm,tpmlts = \
			TensorPowerRootMultiplicities(tpydx,dgix,ixvec,rtdgns)
		rtsum = [0, None, None]
		
		inited = False
		for ixmult in tpdgrm:
			rpi = rep[ixmult[1]]
			for k in xrange(1,2+1):
				if inited:
					muladdtofunc(rtsum[k],lax,ixmult[0],rpi[k])
				else:
					rtsum[k] = mulfunc(lax,ixmult[0],rpi[k])
			inited = True
		for k in xrange(1,2+1):
			rtsum[k] = tuplefunc(lax,rtsum[k])
			
		for i in xrange(dglen):
			rtsum[0] = tpmlts[i]
			if useentry(nsa, rtsum[2]):
				AddEntryToIndexedRep(rppwrdcmp[i],rtsum)
	
	return tuple([ExtractRepIrrepsGeneral(lax,rppwrdcmp[i],maxfunc,repfunc) \
		for i in xrange(dglen)])

def RPMul(la,mlt,rt): return mlt*rt

def RPMulAddTo(rtsum,la,mlt,rt): rtsum += mlt*rt

def RPTuple(la,vec): return tuple(vec)

def DecomposeRepPowerXtnd(latype, rptype, maxwts, pwr):
	# The individual reps
	rep = GetRepXtnd(latype,rptype,maxwts)
	if rep == None: return None
	return DecomposeRepRootsPower(latype, rep, pwr, \
		RPMul, RPMulAddTo, RPTuple, \
		MaxRepRoot, \
		GetRepOrbits, UEOrbits, 0)

# Weyl - orbit dom wts -- all entries : GetRep, UEAll, 0 

def DecomposeRepPower(latype, maxwts, pwr):
	return DecomposeRepPowerXtnd(latype, "sngl", maxwts, pwr)

def AlgProdRPMul(lalist,mlt,rt):
	lalen = len(lalist)
	res = [mlt*rt[i] for i in xrange(lalen)]
	res += mul_sv(mlt,rt[lalen:])
	return res

def AlgProdRPMulAddTo(rtsum,lalist,mlt,rt):
	lalen = len(lalist)
	for i in xrange(lalen):
		rtsum[i] += mlt*rt[i]
	for i in xrange(lalen,len(rtsum)):
		rtsum[i] += mlt*rt[i]

def AlgProdRPTuple(lalist,rt):
	# Necessary to make a root vector hashable,
	# which is necessary for counting them up
	lalen = len(lalist)
	res = [tuple(r) for r in rt[:lalen]]
	res += rt[lalen:]
	return tuple(res)

def DecomposeAlgProdRepPowerXtnd(latlist, rptype, maxwts, pwr):
	# The individual reps
	rep = GetAlgProdRepXtnd(latlist,rptype,maxwts)
	if rep == None: return None
	return DecomposeRepRootsPower(latlist, rep, pwr, \
		AlgProdRPMul, AlgProdRPMulAddTo, AlgProdRPTuple, \
		MaxAlgProdRepRoot, \
		GetAlgProdRepOrbits, UEAlgProdOrbits, len(latlist))

# Weyl - orbit dom wts -- all entries : GetAlgProdRep, UEAll, 0

def DecomposeAlgProdRepPower(latlist, maxwts, pwr):
	return DecomposeAlgProdRepPowerXtnd(latlist, "sngl", maxwts, pwr)

# Pure symmetric and antisymmetric cases: sym = +1 and -1

def DecomposeRepRootsPwrSym(lax, rep, pwr, sym, \
		mulfunc, muladdtofunc, tuplefunc, maxfunc, \
		repfunc, useentry, nsa):
	
	# Zero : scalar rep
	if pwr <= 0:
		return ((1,tuplefunc(lax,mulfunc(lax,0,rep[0][2]))),)
	
	# Fix the symmetry factor
	symx = 1 if sym >= 0 else -1
	
	# Use index function rppwrix to speed up lookups in the power-rep list;
	# index with the weight vector.
	# Also pull out the root multiplicities for convenience
	rtdgns = [r[0] for r in rep]
	rppwrdcmp = {}
	
	for ixvec in MakeIndexVectors(pwr,len(rtdgns)):
		tpdgdict = memcount(ixvec)
		tpdgrm = zip(tpdgdict.values(),tpdgdict.keys())
		tpdgrm.sort(lambda a,b: cmp(b[0],a[0]))
		tpdgrm = tuple(tpdgrm)
		tpmlts = 1
		for dgi in tpdgrm:
			dnm = dgi[0]
			nmi = rtdgns[dgi[1]]
			for j in xrange(dnm):
				tpmlts *= (nmi + symx*j)
				tpmlts /= (j+1)
		tpmlts = int(tpmlts)
		
		rtsum = [0, None, None]
		
		inited = False
		for ixmult in tpdgrm:
			rpi = rep[ixmult[1]]
			for k in xrange(1,2+1):
				if inited:
					muladdtofunc(rtsum[k],lax,ixmult[0],rpi[k])
				else:
					rtsum[k] = mulfunc(lax,ixmult[0],rpi[k])
			inited = True
		for k in xrange(1,2+1):
			rtsum[k] = tuplefunc(lax,rtsum[k])
		
		rtsum[0] = tpmlts
		if useentry(nsa, rtsum[2]):
			AddEntryToIndexedRep(rppwrdcmp,rtsum)
	
	return ExtractRepIrrepsGeneral(lax,rppwrdcmp,maxfunc,repfunc)

def DecomposeRepPwrSymXtnd(latype, rptype, maxwts, pwr, sym):
	# The individual reps
	rep = GetRepXtnd(latype,rptype,maxwts)
	if rep == None: return None
	return DecomposeRepRootsPwrSym(latype, rep, pwr, sym, \
		RPMul, RPMulAddTo, RPTuple, \
		MaxRepRoot, \
		GetRepOrbits, UEOrbits, 0)

# Weyl - orbit dom wts -- all entries : GetRep, UEAll, 0 

def DecomposeRepPwrSym(latype, maxwts, pwr, sym):
	return DecomposeRepPwrSymXtnd(latype, "sngl", maxwts, pwr, sym)

def DecomposeAlgProdRepPwrSymXtnd(latlist, rptype, maxwts, pwr, sym):
	# The individual reps
	rep = GetAlgProdRepXtnd(latlist,rptype,maxwts)
	if rep == None: return None
	return DecomposeRepRootsPwrSym(latlist, rep, pwr, sym, \
		AlgProdRPMul, AlgProdRPMulAddTo, AlgProdRPTuple, \
		MaxAlgProdRepRoot, \
		GetAlgProdRepOrbits, UEAlgProdOrbits, len(latlist))

# Weyl - orbit dom wts -- all entries : GetAlgProdRep, UEAll, 0

def DecomposeAlgProdRepPwrSym(latlist, maxwts, pwr, sym):
	return DecomposeAlgProdRepPwrSymXtnd(latlist, "sngl", maxwts, pwr, sym)

# Find Young-diagram products using the Littlewood-Richardson rule.
# Need all of them for each power.
# These are for doing products of rep-power decompositions,
# and also for doing decompositions of powers of sums of reps.

def YDClipEmpty(dgrm):
	newdgrm = dgrm
	for k, drow in enumerate(dgrm):
		if len(drow) == 0:
			newdgrm = dgrm[:k]
			break
	return newdgrm

def YDAddOne(dgrm, token, start):
	# token must be absent from the original diagram
	olddgrm = dgrm + [[]]
	outlist = []
	for k in xrange(start,len(olddgrm)):
		# Explicit copy to avoid inadvertently adding to the lists
		# inside a diagram list.
		# Unlike Mathematica, Python does not have automatic copy-on-write
		newdgrm = [list(dgrow[:]) for dgrow in olddgrm]
		newdgrm[k].append(token)
		ndlen = len(newdgrm[k])
		if k > 0:
			apnd = ndlen <= len(newdgrm[k-1])
		else:
			apnd = True
		if apnd:
			for l in xrange(k):
				if newdgrm[l][ndlen-1] == token:
					apnd = False
					break
		if apnd:
			outlist.append((newdgrm,k))
	
	outlist = [(YDClipEmpty(outmem[0]), outmem[1]) for outmem in outlist]
	return outlist

def YDAddSet(dgrm,token,number):
	# token must be absent from the original diagram
	outlist = [(dgrm,0)]
	oldlist = outlist[:]
	for k in xrange(number):
		outlist = []
		for oldmem in oldlist:
			outlist += YDAddOne(oldmem[0],token,oldmem[1])
		oldlist = outlist[:]
	
	if len(outlist) > 0:
		outlist = [outmem[0] for outmem in outlist]
	return outlist
	
# Works with expanded diagrams; returns a set of expanded diagrams
# with boxes labeled as a result of the algorithm.
# 0 is first diagram and 1 2 3 are the rows of the second diagram.
def YDProduct(dgrm1,dgrm2):
	# Be sure to keep all the tokens distinct;
	# 0 for the original diagram is distinct from
	# 1, 2, ... for the added one.
	dgrmx = YDClipEmpty([zeros1(dgrow) for dgrow in dgrm1])
	
	# Be sure that the the original diagram is all-list
	oldlist = [dgrmx]
	for k,dgrow2 in enumerate(dgrm2):
		outlist = []
		for oldmem in oldlist:
			outlist += YDAddSet(oldmem, k+1, dgrow2)
		oldlist = outlist[:]

	outlist = []
	for dgx in oldlist:
		cnt = zeros1(len(dgrm2))
		apnd = True
		for dgrow in dgx:
			for x in reversed(dgrow):
				if x > 0:
					if x > 1:
						if cnt[x-1] >= cnt[x-2]:
							apnd = False
							break
					cnt[x-1] += 1
			if not apnd: break
		if apnd: outlist.append(dgx)
	
	# Don't need to be lists anymore
	return tuple([tuple([len(dgrow) for dgrow in dgx]) for dgx in outlist])
	
# Output: list of matrices of length-k diagrams * length-(pwr-k) diagrams.
# The products are list of indices of the diagrams.
# Repeated indices are too rare to be worthy of the counted-index approach.

def YDPowerProducts(pwr):
	dgrms = GetTensorPowerYDX(pwr).yds()
	dgix = MakeYDSetIndexing(dgrms)
	
	return tuple([ tuple([ tuple([ tuple([ \
		dgix[dgrm] for dgrm in YDProduct(dgrm1,dgrm2)]) \
		for dgrm2 in GetTensorPowerYDX(pwr-k).yds()]) \
		for dgrm1 in GetTensorPowerYDX(k).yds()]) \
		for k in xrange(1,pwr)])

YDPowerProductList = []

def GetYDPowerProds(n):
	llen = len(YDPowerProductList)
	if llen < n:
		for k in xrange(llen+1,n+1):
			YDPowerProductList.append(YDPowerProducts(k))
	
	return YDPowerProductList[n-1]

# Not yet translated into Python:
# Products of rep-power decompositions
# Decompositions of powers of sums of reps
# The first is for verification of the general applicability
# of the Littlewood-Richardson rule,
# while the second is slower than doing a rep power on a rep list,
# from timing tests in Mma.


# Subalgebras and branching rules

# Doing the Branching
#
# Members:
# latype: Original-algebra type
# stsms: list for each subalgebra of
#    (type, projection matrix, numerator of integerized projmat, its denominator)
# u1s: Indices of U(1)-factor roots (1-based)
#
# SubAlgebras(): the subalgebra types
# SubAlgebraNmes(): their names
#
# DoBranching: method for doing the branching
# Arg: max-weights vector of original irrep
class SubAlgebraBrancher:
	def __init__(self,latype,stsms,u1s):
		self.latype = tuple(latype)
		self.stsms = tuple(map(ExpandSTSM,stsms))
		self.u1s = tuple(u1s)
		
	def verify(self):
		if GetLieAlgebra(self.latype) == None: return False
		for subtype in self.subtypes:
			if GetLieAlgebra(subtype) == None: return False
		return True
	
	def SubAlgebras(self): return tuple([s[0] for s in self.stsms])
	def SubAlgebraNames(self): return ','.join(map(AlgName,self.SubAlgebras()))
	
	# Self-description
	def __str__(self):
		return "Lie subalgebra brancher: type=%s %s subtypes=%s %s #u1s=%d" % \
			(str(self.latype), AlgName(self.latype), \
				str(self.SubAlgebras()), self.SubAlgebraNames(), len(self.u1s))
	
	def __repr__(self): return str(self)

def ExpandSTSM(stsm):
	subtype = tuple(stsm[0])
	submat = MatrixTuple(stsm[1])
	submat_num, submat_den = SharedDen_Matrix(mul_sm(Fraction(1,1),submat))
	submat_num = MakeNP(submat_num)
	return (subtype, submat, submat_num, submat_den)

# Doing the actual branching

def FindMaxRtDcmp(rdlist):
	if len(rdlist) < 0: return None
	maxrd = rdlist[0]
	ns = len(maxrd[0])
	for rd in rdlist[1:]:
		mindiff = min([min(sub_vv(rd[0][k],maxrd[0][k])) for k in xrange(ns)])
		if mindiff >= 0:
			maxrd = rd		
	return maxrd

# Args: branching object, list of representation roots
# Returns a counted list, a tuple of
#   (degen, tuple of subalgebras' maxweight vectors with U(1) factors)
def SAB_DoBranchingRoots(ld,rep):
	# Denominator for source algebra's roots
	la = GetLieAlgebra(ld.latype)
	laden = la.ictnden
	# Subalgebras and their root-vector denominators
	styps = [stsm[0] for stsm in ld.stsms]
	slas = [GetLieAlgebra(styp) for styp in styps]
	slactns = [sla.ctnmat for sla in slas]
	sladens = [sla.ictnden for sla in slas]
	# Matrix numerators and denominators
	smnums = [stsm[2] for stsm in ld.stsms]
	smdens = [stsm[3] for stsm in ld.stsms]
	
	rpix = {}
	
	for rtx in rep:
		degen = rtx[0]
		root = rtx[1]
		subroots = []
		subweights = []
		for k in xrange(len(ld.stsms)):
			srt = np.dot(root,smnums[k])
			srt = [(sladens[k]*r)/(laden*smdens[k]) for r in srt]
			swt = np.dot(srt,slactns[k])
			swt = [r/sladens[k] for r in swt]
			subroots.append(tuple(srt))
			subweights.append(tuple(swt))
		# Handle both the cases of a vector to multiply the roots and a root index
		# Root indices are 1-based
		u1vals = []
		for u1fac in ld.u1s:
			if type(u1fac) == type(0):
				u1val = root[u1fac-1]
			else:
				u1val = mul_vv(root,u1fac)
			u1vals.append(u1val*Fraction(1,laden))
		subroots += u1vals
		subweights += u1vals
		rtdcmp = (degen, tuple(subroots), tuple(subweights))
		if UEAlgProdOrbits(len(styps), rtdcmp[2]):
			AddEntryToIndexedRep(rpix,rtdcmp)
	
	return ExtractAlgProdRepOrbitIrreps(styps, rpix)

def SAB_DoBranchingXtnd(ld, rptype, maxwts):
	return SAB_DoBranchingRoots(ld, GetRepXtnd(ld.latype, rptype, maxwts))

def SAB_DoBranching(ld, maxwts):
	return SAB_DoBranchingXtnd(ld, "sngl", maxwts)

SubAlgebraBrancher.DoBranchingRoots = SAB_DoBranchingRoots
SubAlgebraBrancher.DoBranchingXtnd = SAB_DoBranchingXtnd
SubAlgebraBrancher.DoBranching = SAB_DoBranching


# Makes submats intended for weight space.
# Args: size of original, special column, 
# list of original root for each result root.
# Original root = -1 makes the special column,
# intended for extension splitting.
def WtSpcSubmat(n, spccol, origrts):
	m = len(origrts)
	smat = zeros2(n,m)
	for i,ix in enumerate(origrts):
		if ix == -1:
			for j in xrange(n):
				smat[j][i] = spccol[j]
		elif ix >= 1 and ix <= n:
			smat[ix-1][i] = 1
	return smat

# To be compatible with Fraction():
def IntegerMatrix(m): return [[int(mi) for mi in mr] for mr in m]

# Weight space to root space
def SubmatWtToRoot(latype, sublt, smat):
	la = GetLieAlgebra(latype)
	cmt = IntegerMatrix(la.ctnmat)
	sla = GetLieAlgebra(sublt)
	slc = sla.invctn
	return mul_mm(cmt, mul_mm(smat,slc))

def WSRtSubmat(latype, spccol, sgmt):
	sublt = sgmt[0]
	origrts = sgmt[1]
	smat = WtSpcSubmat(latype[1], spccol, origrts)
	return (sublt, SubmatWtToRoot(latype, sublt, smat))

def WSRtSbmtList(latype, spccol, sgmtlist):
	return [WSRtSubmat(latype, spccol, sgmt) for sgmt in sgmtlist]

# Multiple Root Demotion:
# Makes U(1) factors from a list of roots and splits the algebra.
# Args: type and which roots to split at (1-based)

def SplitByDemotedRoot(algrts,drt):
	latype = algrts[0]
	family = latype[0]
	n = latype[1]
	rts = algrts[1]
	if drt not in rts: return (algrts,)
	
	# Find the location of the root to split at
	for i, rt in enumerate(rts):
		if rt == drt:
			ix = i+1
			break
	
	# Create a list of
	# (subalgebra type), (which orig roots because which subalgebra roots)
	# The input roots will be substituted in later.
	# 1-based here
	if family == 1: # A(n)
		rsdt = []
		if ix > 1:
			rd = ((1,ix-1), range(1,ix))
			rsdt.append(rd)
		if ix < n:
			rd = ((1,n-ix), range(ix+1,n+1))
			rsdt.append(rd)	
	elif family == 2: # B(n)
		rsdt = []
		if ix > 1:
			rd = ((1,ix-1), range(1,ix))
			rsdt.append(rd)
		if ix < n:
			rd = ((2,n-ix), range(ix+1,n+1))
			rsdt.append(rd)	
	elif family == 3: # C(n)
		rsdt = []
		if ix > 1:
			rd = ((1,ix-1), range(1,ix))
			rsdt.append(rd)
		if ix < n:
			rd = ((3,n-ix), range(ix+1,n+1))
			rsdt.append(rd)	
	elif family == 4: # D(n)
		rsdt = []
		if ix > 1 and ix < n-1:
			rd = ((1,ix-1), range(1,ix))
			rsdt.append(rd)
		if ix < n-1:
			rd = ((4,n-ix), range(ix+1,n+1))
			rsdt.append(rd)
		if ix == n-1:
			rd = ((1,n-1), range(1,n-1) + [n])
			rsdt.append(rd)
		if ix == n:
			rd = ((1,n-1), range(1,n))
			rsdt.append(rd)
	elif family == 5:
		if n == 6: # E6
			if ix == 1:
				rsdt = (((4, 5), (5, 4, 3, 2, 6)), )
			elif ix == 2:
				rsdt = (((1, 1), (1, )), ((1, 4), (6, 3, 4, 5)))
			elif ix == 3:
				rsdt = (((1, 2), (1, 2)), ((1, 2), (4, 5)), ((1, 1), (6, )))
			elif ix == 4:
				rsdt = (((1, 4), (1, 2, 3, 6)), ((1, 1), (5, )))
			elif ix == 5:
				rsdt = (((4, 5), (1, 2, 3, 4, 6)), )
			elif ix == 6:
				rsdt = (((1, 5), (1, 2, 3, 4, 5)), )
		elif n == 7: # E7
			if ix == 1:
				rsdt = (((4, 6), (6, 5, 4, 3, 2, 7)), )
			elif ix == 2:
				rsdt = (((1, 1), (1, )), ((1, 5), (7, 3, 4, 5, 6)))
			elif ix == 3:
				rsdt = (((1, 2), (1, 2)), ((1, 3), (4, 5, 6)), ((1, 1), (7, )))
			elif ix == 4:
				rsdt = (((1, 4), (1, 2, 3, 7)), ((1, 2), (5, 6)))
			elif ix == 5:
				rsdt = (((4, 5), (1, 2, 3, 4, 7)), ((1, 1), (6, )))
			elif ix == 6:
				rsdt = (((5, 6), (1, 2, 3, 4, 5, 7)), )
			elif ix == 7:
				rsdt = (((1, 6), (1, 2, 3, 4, 5, 6)), )
		elif n == 8: # E8
			if ix == 1:
				rsdt = (((4, 7), (7, 6, 5, 4, 3, 2, 8)), )
			elif ix == 2:
				rsdt = (((1, 1), (1, )), ((1, 6), (8, 3, 4, 5, 6, 7)))
			elif ix == 3:
				rsdt = (((1, 2), (1, 2)), ((1, 4), (4, 5, 6, 7)), ((1, 1), (8, )))
			elif ix == 4:
				rsdt = (((1, 4), (1, 2, 3, 8)), ((1, 3), (5, 6, 7)))
			elif ix == 5:
				rsdt = (((4, 5), (1, 2, 3, 4, 8)), ((1, 2), (6, 7)))
			elif ix == 6:
				rsdt = (((5, 6), (1, 2, 3, 4, 5, 8)), ((1, 1), (7, )))
			elif ix == 7:
				rsdt = (((5, 7), (1, 2, 3, 4, 5, 6, 8)), )
			elif ix == 8:
				rsdt = (((1, 7), (1, 2, 3, 4, 5, 6, 7)), )
	elif family == 6:
		if n == 4: # F4
			if ix == 1:
				rsdt = (((3, 3), (4, 3, 2)), )
			elif ix == 2:
				rsdt = (((1, 1), (1, )), ((1, 2), (3, 4)))
			elif ix == 3:
				rsdt = (((1, 2), (1, 2)), ((1, 1), (4, )))
			elif ix == 4:
				rsdt = (((2, 3), (1, 2, 3)), )
	elif family == 7:
		if n == 2: # G2
			if ix == 1:
				rsdt = (((1, 1), (2, )), )
			elif ix == 2:
				rsdt = (((1, 1), (1, )), )

	# Plug in which actual original roots and return it
	orrsdt = []
	for rd in rsdt:
		orrd = (rd[0], [rts[i-1] for i in rd[1]])
		orrsdt.append(orrd)
	return orrsdt

def MakeMultiRootDemoter(latype,dmrts):
	n = latype[1]
	if min(dmrts) < 1: return None
	if max(dmrts) > n: return None
	if len(dmrts) > len(set(dmrts)): return None
	splt = ((latype,range(1,n+1)),)
	
	# Split by each one separately -- easiest to implement
	for drt in dmrts:
		newsplt = []
		for slmem in splt:
			newsplt += SplitByDemotedRoot(slmem,drt)
		splt = newsplt
	
	sdlist = WSRtSbmtList(latype, None, splt)
	return SubAlgebraBrancher(latype, sdlist, dmrts)

# Root Demotion:
# Makes a U(1) factor from a root and splits the algebra.
# Args: type and which root to split at (1-based)
def MakeRootDemoter(latype,m):
	return MakeMultiRootDemoter(latype,(m,))

def ListRootDemotions(latype):
	for m in xrange(latype[1]):
		ld = MakeRootDemoter(latype,m)
		print m, ld.SubAlgebras(), ld.SubAlgebraNames()


# Extension Splitting
# Extension of the Dynkin diagram with the adjoint-connected root
# and splitting of the algebra by removing a root.
# Args: type and which root to split at (1-based)
def MakeExtensionSplitter(latype,m):
	family = latype[0]
	n = latype[1]
	if family == 4 and n <= 3: return None # D2 = A1*A1, D3 = A3, don't use
	if m < 1 or m > n: return None
	u1s = ()
	
	splt = []
	if family == 1: # A(n)
		splt.append((latype, range(m+1,n+1)+[-1]+range(1,m)))
	elif family == 2: # B(n)
		if m == 1:
			splt.append((latype, [-1]+range(2,n+1)))
		else:
			splt.append(((4,m), range(m-1,0,-1)+[-1]))
			if m < n:
				splt.append(((2,n-m), range(m+1,n+1)))
	elif family == 3: # C(n)
		splt.append(((3,m), range(m-1,0,-1)+[-1]))
		if m < n:
			splt.append(((3,n-m), range(m+1,n+1)))
	elif family == 4: # D(n)
		if m == 1:
			splt.append((latype, [-1]+range(2,n+1)))
		elif m >= n-1:
			mx = (2*n-1) - m
			splt.append((latype, [mx]+range(n-2,0,-1)+[-1]))
		else:
			splt.append(((4,m), range(m-1,0,-1)+[-1]))
			if m < n:
				splt.append(((4,n-m), range(m+1,n+1)))
	elif family == 5:
		if n == 6: # E6
			if m == 1:
				splt = (((5, 6), (5, 4, 3, 6, -1, 2)), )
			elif m == 2:
				splt = (((1, 1), (1, )), ((1, 5), (5, 4, 3, 6, -1)))
			elif m == 3:
				splt = (((1, 2), (1, 2)), ((1, 2), (4, 5)), ((1, 2), (6, -1)))
			elif m == 4:
				splt = (((1, 5), (1, 2, 3, 6, -1)), ((1, 1), (5, )))
			elif m == 5:
				splt = (((5, 6), (1, 2, 3, 6, -1, 4)), )
			elif m == 6:
				splt = (((1, 5), (1, 2, 3, 4, 5)), ((1, 1), (-1, )))
		elif n == 7: # E7
			if m == 1:
				splt = (((1, 1), (-1, )), ((4, 6), (6, 5, 4, 3, 2, 7)))
			elif m == 2:
				splt = (((1, 2), (-1, 1)), ((1, 5), (6, 5, 4, 3, 7)))
			elif m == 3:
				splt = (((1, 3), (-1, 1, 2)), ((1, 3), (4, 5, 6)), ((1, 1), (7, )))
			elif m == 4:
				splt = (((1, 5), (-1, 1, 2, 3, 7)), ((1, 2), (5, 6)))
			elif m == 5:
				splt = (((4, 6), (-1, 1, 2, 3, 4, 7)), ((1, 1), (6, )))
			elif m == 6:
				splt = (((5, 7), (5, 4, 3, 2, 1, -1, 7)), )
			elif m == 7:
				splt = (((1, 7), (-1, 1, 2, 3, 4, 5, 6)), )
		elif n == 8: # E8
			if m == 1:
				splt = (((4, 8), (-1, 7, 6, 5, 4, 3, 2, 8)), )
			elif m == 2:
				splt = (((1, 1), (1, )), ((1, 7), (-1, 7, 6, 5, 4, 3, 8)))
			elif m == 3:
				splt = (((1, 2), (1, 2)), ((1, 5), (4, 5, 6, 7, -1)), ((1, 1), (8, )))
			elif m == 4:
				splt = (((1, 4), (1, 2, 3, 8)), ((1, 4), (5, 6, 7, -1)))
			elif m == 5:
				splt = (((4, 5), (1, 2, 3, 4, 8)), ((1, 3), (6, 7, -1)))
			elif m == 6:
				splt = (((5, 6), (1, 2, 3, 4, 5, 8)), ((1, 2), (7, -1)))
			elif m == 7:
				splt = (((5, 7), (1, 2, 3, 4, 5, 6, 8)), ((1, 1), (-1, )))
			elif m == 8:
				splt = (((1, 8), (1, 2, 3, 4, 5, 6, 7, -1)), )
	elif family == 6:
		if n == 4: # F4
			if m == 1:
				splt = (((1, 1), (-1, )), ((3, 3), (4, 3, 2)))
			elif m == 2:
				splt = (((1, 2), (-1, 1)), ((1, 2), (3, 4)))
			elif m == 3:
				splt = (((1, 3), (2, 1, -1)), ((1, 1), (4, )))
			elif m == 4:
				splt = (((2, 4), (-1, 1, 2, 3)), )
	elif family == 7:
		if n == 2: # G2
			if m == 1:
				splt = (((1, 1), (-1, )), ((1, 1), (2, )))
			elif m == 2:
				splt = (((1, 2), (-1, 1)), )
	
	if len(splt) == 0: return None
	
	la = GetLieAlgebra(latype)
	metdiag = [la.metric[i][i] for i in xrange(n)]
	mdmax = max(metdiag)
	psrtmax = la.posroots[-1][0]
	spccol = [- psrtmax[i]*metdiag[i]/mdmax for i in xrange(n)]
	
	sdlist = WSRtSbmtList(latype, spccol, splt)
	return SubAlgebraBrancher(latype, sdlist, u1s)

def ListExtensionSplits(latype):
	for m in xrange(latype[1]):
		ld = MakeExtensionSplitter(latype,m)
		print m, ld.SubAlgebras(), ld.SubAlgebraNames()


# Additional Branching-Rule Generators


# The first two here have an interesting interpretation: 
# the original Lie-group matrices become outer products of the subgroup matrices. 
# Dynkin has proved that these are the only maximal nonsimple
# subalgebra decompositions of the classical Lie algebras.

# Determination of the simple subalgebras can be very difficult, it must be said. 
# This document will thus have only some of the more interesting or easier-to-find ones

# Reduces SU(n) to a product of SU(m)'s where n = product of m's.
def SubalgMultSU(suords):
	nso = len(suords)
	soprod = 1
	for k in xrange(nso): soprod *= suords[k]
	n = soprod
	nt = n - 1
	latype = (1,nt)
	vecproj = zeros2(nt,n)
	for k in xrange(nt):
		vecproj[k][k] = 1
		vecproj[k][k+1] = -1
	stsms = []
	for kso in xrange(nso):
		sm = suords[kso]
		smt = sm - 1
		stype = (1, smt)
		smat = zeros2(n,smt)
		sstr = 1
		for k in xrange(kso+1,nso): sstr *= suords[k]
		for k in xrange(n):
			for ks in xrange(smt):
				ko = (k / sstr) % sm
				if ko == ks: smat[k][ks] = 1
				if ko == ks+1: smat[k][ks] = -1
		la = GetLieAlgebra(stype)
		smat = mul_mm(smat,la.invctn)
		smat = mul_mm(vecproj,smat)
		stsms.append((stype,smat))

	u1s = ()
	return SubAlgebraBrancher(latype,stsms,u1s)

# Reduces SO(n) and Sp(n) to a product of SO(m)'s and Sp(m)' s, 
# where n = product of m's. 
# Positive n means SO(n) and negative n means Sp(-n).
# SO(2) is handled as a U(1) factor.

def SubalgMultSOSp(sospords):
	nso = len(sospords)
	soprod = 1
	for k in xrange(nso): soprod *= sospords[k]
	# Construct a projection matrix for roots onto vector rep. 
	# Also get the type. Much like the previous function, 
	# it's designed to turn the vector rep
	# into crossed diagonal lines of 1 and -1
	if soprod > 0:
		if soprod % 2 == 0:
			# SO(even): Dn
			n = soprod
			nt = n/2
			latype = (4,nt)
			vpbase = zeros2(nt,nt)
			for k in xrange(nt-1):
				vpbase[k][k] = 1
				vpbase[k+1][k] = -1
			vpbase[nt-2][nt-1] = vpbase[nt-1][nt-1] = 1
			vecproj = zeros2(nt,n)
			for k in xrange(nt):
				for kx in xrange(nt):
					vecproj[k][kx] = vpbase[kx][k]
					vecproj[k][(n-1)-kx] = - vpbase[kx][k]
		else:
			# SO(odd): Bn
			n = soprod
			nt = (n-1)/2
			latype = (2,nt)
			vpbase = zeros2(nt,nt)
			for k in xrange(nt-1):
				vpbase[k][k] = 1
				vpbase[k+1][k] = -1
			vpbase[nt-1][nt-1] = 1
			vecproj = zeros2(nt,n)
			for k in xrange(nt):
				for kx in xrange(nt):
					vecproj[k][kx] = vpbase[kx][k]
					vecproj[k][(n-1)-kx] = - vpbase[kx][k]
	else:
		if soprod % 2 == 0:
			# Sp(even): Cn
			n = - soprod
			nt = n/2
			latype = (3,nt)
			vpbase = zeros2(nt,nt)
			for k in xrange(nt-1):
				vpbase[k][k] = 1
				vpbase[k+1][k] = -1
			vpbase[nt-1][nt-1] = 2
			vecproj = zeros2(nt,n)
			for k in xrange(nt):
				for kx in xrange(nt):
					vecproj[k][kx] = vpbase[kx][k]
					vecproj[k][(n-1)-kx] = - vpbase[kx][k]
		else:
			# Sp(odd) does not exist
			return None
	
	# Now the subalgebra matrices, with SO(2)/D(1) as a U(1) factor 
	stsms = []
	u1s = []
	for kso, sm in enumerate(sospords):
		if sm > 0:
			if sm % 2 == 0:
				# SO(even): Dn -- works for root-vector dimension 1
				sma = sm
				smt = sma/2
				stype = (4,smt)
				smat = zeros2(n,smt)
				sstr = 1
				for k in xrange(kso+1,nso): sstr *= sospords[k]
				if sstr < 0: sstr *= -1
				for k in xrange(n):
					ko = (k/sstr) % sma
					for ks in xrange(smt):
						if ko == ks: smat[k][ks] = Fraction(1,2)
						if ko == ks+1: smat[k][ks] = -Fraction(1,2)
						if ko == smt-2 and ks == smt-1: smat[k][ks] = Fraction(1,2)
						if ko == (sma-1) - ks: smat[k][ks] = -Fraction(1,2)
						if ko == (sma-1) - (ks+1): smat[k][ks] = Fraction(1,2)
						if ko == smt+1 and ks == smt-1: smat[k][ks] = -Fraction(1,2)
			else:
				# SO(odd): Bn
				sma = sm
				smt = (sma-1)/2
				stype = (2,smt)
				smat = zeros2(n,smt)
				sstr = 1
				for k in xrange(kso+1,nso): sstr *= sospords[k]
				if sstr < 0: sstr *= -1
				for k in xrange(n):
					ko = (k/sstr) % sma
					for ks in xrange(smt):
						if ko == ks: smat[k][ks] = (1 if ks==smt-1 else Fraction(1,2))
						if ko == ks+1: smat[k][ks] = -Fraction(1,2)
						if ko == (sma-1) - ks: smat[k][ks] = -(1 if ks==smt-1 else Fraction(1,2))
						if ko == (sma-1) - (ks+1): smat[k][ks] = Fraction(1,2)
		else:
			if sm % 2 == 0:		
				# Sp(even): Cn
				sma = -sm
				smt = sma/2
				stype = (3,smt)
				smat = zeros2(n,smt)
				sstr = 1
				for k in xrange(kso+1,nso): sstr *= sospords[k]
				if sstr < 0: sstr *= -1
				for k in xrange(n):
					ko = (k/sstr) % sma
					for ks in xrange(smt):
						if ko == ks: smat[k][ks] = Fraction(1,2)
						if ko == ks+1: smat[k][ks] = -Fraction(1,2)
						if ko == (sma-1) - ks: smat[k][ks] = -Fraction(1,2)
						if ko == (sma-1) - (ks+1): smat[k][ks] = Fraction(1,2)
			else:
				# Sp(odd) does not exist
				return None
		smat = mul_mm(vecproj,smat)
		if stype == (2,0):
			# Skip over SO(1)
			pass
		elif stype == (4,1):
			# Treat as a U(1) factor
			u1f = tuple(mul_sv(Fraction(1,2),[sme[0] for sme in smat]))
			u1s.append(u1f)
		else:
			# All the others are "real" algebras
			la = GetLieAlgebra(stype)
			smat = mul_mm(smat,la.invctn)
			stsms.append((stype,smat))
	
	return SubAlgebraBrancher(latype,stsms,u1s)

# Forms of the previous two for A,B,C,D-style designation.

# In: list of root-vector lengths
def SubalgMultAn(ords):
	return SubalgMultSU([ord+1 for ord in ords])

# In: list of algebra types as lists. For B(n),C(n),D(n): (2,n),(3,n),(4,n).
# D(1) or (4,1) is legitimate.
def SubalgMultBCDn(stypes):
	ords = []
	for stype in stypes:
		family = stype[0]; n = stype[1]
		if family == 2:
			ord = 2*n+1
		elif family == 3:
			ord = -2*n
		elif family == 4:
			ord = 2*n
		else:
			return None
		ords.append(ord)
	
	return SubalgMultSOSp(ords)

# This function does even SO -> odd SO or odd SO + odd SO; 
# It does SO(2n) -> SO(2m+1) * SO(2n-2m-1)
# The extension splitter does
# even -> even + even
# odd -> odd + even
# Use the root demoter and demote root #1
# to turn SO(n) into SO(n-2) * SO(2) (the U(1) factor)
#
# For m = 0, it will make only one odd SO; 
# for m >= 1 and m <= n-2, it will make two odd SO's.
def SubalgSOEvenOdd(n,m):
	if n < 2: return None
	if m < 0 or m > n-2: return None
	latype = (4,n)
	stsms = []
	u1s = ()
	if m > 0:
		k = m
		subtype = (2,k)
		smat = zeros2(n,k)
		for i in xrange(2,k): smat[m-i][i-1] = 1
		if m > 1:
			smat[0] = ((k-1)*[-1]) + [0]
			smat[m-1] = [0] + ((k-1)*[-1])
		else:
			smat[0] = [-1]
		stsms.append((subtype,smat))
	k = n - m - 1
	subtype = (2,k)
	smat = zeros2(n,k)
	for i in xrange(1,k): smat[m+i-1][i-1] = 1
	if m > 0:
		smat[m-1] = k*[-1]
	smat[n-1][k-1] = smat[n-2][k-1] = 1
	stsms.append((subtype,smat))
	return SubAlgebraBrancher(latype,stsms,u1s)

# Turns SU(n) / A(n-1) into SO(n) / D(n/2) or B((n-1)/2)
def SubalgSUSO(n):
	if n % 2 == 0:
		# Even n
		nx = n/2
		stype = (4,nx)
		smat = zeros2(n-1,nx)
		for k in xrange(nx-1):
			smat[k][k] = 1
			smat[n-k-2][k] = 1
		smat[nx-1][nx-2] = -1
		smat[nx-1][nx-1] = 1
	else:
		# Odd n
		nx = (n-1)/2
		stype = (2,nx)
		smat = zeros2(n-1,nx)
		for k in xrange(nx):
			smat[k][k] = 1
			smat[n-k-2][k] = 1
	stsms = ((stype, smat),)
	return SubAlgebraBrancher((1,n-1),stsms,())

# Turns SU(2n) / A(2n-1) into Sp(2n) / C(n)
def SubalgSUSp(n):
	smat = zeros2(2*n-1,n)
	for k in xrange(n-1):
		smat[k][k] = 1
		smat[2*n-k-2][k] = 1
	smat[n-1][n-1] = 1
	stsms = (((3,n), smat),)
	return SubAlgebraBrancher((1,2*n-1),stsms,())

# Reduces any algebra to A1, with the rep height becoming the largest highest weight.
# Maximal for B(n), C(n), G2, F4, E7, E8.
# There are some solutions for some subalgebra-matrix values different from 1.
def SubalgHeightA1(latype):
	stsms = (((1,1),latype[1]*[(1,)]),)
	u1s = ()
	return SubAlgebraBrancher(latype,stsms,u1s)

# Some extra ones named individually.
# (Parent algebra) (subalgebras)
#
# Mentioned by John Baez in "The Octonions" :
#
# B3G2 -- B3/SO(7) to G2 -- G2 is the isomorphism group of the octonions,
# and one gets a construction of G2 from it that' s manifestly
# a subgroup of SO(7) -- 14 7D antisymmetric real matrices.
# 
# D4G2 -- D4/SO (8) to G2
#
# The others are mentioned by Slansky.
#
# None of these branchings have U(1) factors

SubalgExtraData = {}

SubalgExtraData["G2A1"] = ((7, 2), (((1, 1), ((1, ), (1, ))), ))

SubalgExtraData["B3G2"] = ((2, 3), (((7, 2), ((0, 1), (1, 0), (0, 1))), ))

SubalgExtraData["D4A2"] = ((4, 4), (((1, 2), ((1, 0), (-1, 1), (1, 0), (1, 0))), ))

SubalgExtraData["D4G2"] = ((4, 4), (((7, 2), ((0, 1), (1, 0), (0, 1), (0, 1))), ))

SubalgExtraData["F4A1"] = ((6, 4), (((1, 1), ((1, ), (1, ), (1, ), (1, ))), ))

SubalgExtraData["F4A1G2"] = ((6, 4), (((1, 1), ((0, ), (0, ), (0, ), (1, ))), ((7, 2), ((2, 3), (-1, 0), (0, -1), (0, 0)))))

SubalgExtraData["E6A2"] = ((5, 6), (((1, 2), ((1, 0), (0, 1), (-1, 0), (0, 1), (1, 0), (2, -1))), ))

SubalgExtraData["E6G2"] = ((5, 6), (((7, 2), ((0, 1), (1, 0), (0, 1), (1, 0), (0, 1), (-1, 0))), ))

SubalgExtraData["E6C4"] = ((5, 6), (((3, 4), ((0, 1, 0, 0), (1, 0, 0, 0), (-1, 0, 1, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 0, 1))), ))

SubalgExtraData["E6F4"] = ((5, 6), (((6, 4), ((0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (1, 0, 0, 0))), ))

SubalgExtraData["E6A2G2"] = ((5, 6), (((1, 2), ((1, 0), (0, 1), (-1, -1), (1, 0), (0, 1), (0, 0))), ((7, 2), ((0, 0), (0, 0), (0, 1), (0, 0), (0, 0), (1, 0)))))

SubalgExtraData["E7A1-1"] = ((5, 7), (((1, 1), ((1, ), (1, ), (1, ), (1, ), (1, ), (1, ), (1, ))), ))

SubalgExtraData["E7A1-2"] = ((5, 7), (((1, 1), ((1, ), (1, ), (0, ), (1, ), (1, ), (1, ), (1, ))), ))

SubalgExtraData["E7A2"] = ((5, 7), (((1, 2), ((0, 2), (0, 2), (0, -1), (1, 0), (0, 1), (1, 0), (0, -2))), ))

SubalgExtraData["E7A1A1"] = ((5, 7), (((1, 1), ((2, ), (1, ), (-2, ), (1, ), (1, ), (1, ), (-2, ))), ((1, 1), ((0, ), (0, ), (0, ), (1, ), (0, ), (0, ), (0, )))))

SubalgExtraData["E7A1G2"] = ((5, 7), (((1, 1), ((0, ), (0, ), (0, ), (0, ), (0, ), (0, ), (1, ))), ((7, 2), ((0, -1), (-2, -2), (0, -1), (1, 1), (0, 1), (0, 1), (1, 2)))))

SubalgExtraData["E7A1F4"] = ((5, 7), (((1, 1), ((0, ), (0, ), (0, ), (0, ), (0, ), (1, ), (0, ))), ((6, 4), ((0, 0, 1, 0), (1, 2, 2, 1), (-2, -3, -4, -2), (1, 2, 2, 1), (0, 0, 1, 0), (0, 0, 0, 1), (1, 0, 0, 0)))))

SubalgExtraData["E7G2C3"] = ((5, 7), (((7, 2), ((0, 0), (0, 0), (0, -1), (1, 2), (0, 0), (0, 0), (-1, 0))), ((3, 3), ((0, 1, 0), (1, 0, 0), (0, 1, 1), (-1, -2, -1), (0, 1, 0), (1, 0, 0), (0, 0, 0)))))

SubalgExtraData["E8A1-1"] = ((5, 8), (((1, 1), ((1, ), (1, ), (1, ), (1, ), (1, ), (1, ), (1, ), (1, ))), ))

SubalgExtraData["E8A1-2"] = ((5, 8), (((1, 1), ((1, ), (1, ), (0, ), (1, ), (1, ), (1, ), (1, ), (1, ))), ))

SubalgExtraData["E8A1-3"] = ((5, 8), (((1, 1), ((1, ), (1, ), (0, ), (1, ), (0, ), (1, ), (1, ), (1, ))), ))

SubalgExtraData["E8B2"] = ((5, 8), (((2, 2), ((0, 1), (0, -1), (1, 0), (0, 1), (0, 1), (-1, 0), (0, 1), (0, -1))), ))

SubalgExtraData["E8A1A2"] = ((5, 8), (((1, 1), ((0, ), (0, ), (0, ), (0, ), (0, ), (0, ), (0, ), (1, ))), ((1, 2), ((1, 0), (0, 1), (-1, 0), (1, 1), (0, -1), (1, 0), (0, 1), (-1, -2)))))

SubalgExtraData["E8G2F4"] = ((5, 8), (((7, 2), ((0, 0), (0, 1), (1, 1), (-1, -1), (0, 0), (0, 0), (0, 0), (0, -1))), ((6, 4), ((0, 0, 0, -1), (-1, -1, -1, 0), (0, -1, -2, -1), (1, 1, 2, 1), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (0, 1, 2, 1)))))

# Slansky's book includes several for the classical algebras that I could not
# fit into any pattern, at least not just yet

# Call with one of these character-string names
def SubalgExtra(saname):
	if saname not in SubalgExtraData: return None
	data = SubalgExtraData[saname]
	latype = data[0]
	stsms = data[1]
	u1s = ()
	return SubAlgebraBrancher(latype,stsms,u1s)


# These all make brancher object ld from their input brancher objects
# and other data
# Indexing is 1-based
def ConcatBranchers(ld0, ix, ld1):
	if ix <= 0 or ix > len(ld0.stsms): return
	if ld0.stsms[ix-1][0] != ld1.latype: return
	submat0 = ld0.stsms[ix-1][1]
	
	stsms = []
	for i,stsm0 in enumerate(ld0.stsms):
		if i == (ix-1):
			for stsms1 in ld1.stsms:
				stsms.append((stsms1[0],mul_mm(submat0,stsms1[1])))
		else:
			stsms.append(stsm0)
	
	newu1s = list(ld0.u1s)
	for u1fac in ld1.u1s:
		if type(u1fac) == type(0):
			u1vec = [1 if (i == u1fac) else 0 for i in xrange(1,ld1.latype[1]+1)]
		else:
			u1vec = u1fac
		newu1s.append(mul_mv(submat0,u1vec))
	
	return SubAlgebraBrancher(ld0.latype, stsms, newu1s)
	

def BrancherRenameA1B1C1(ld0, ix, newfam):
	if ix <= 0 or ix > len(ld0.stsms): return
	if newfam < 1 or newfam > 3: return
	
	stsms = []
	for i,stsm0 in enumerate(ld0.stsms):
		if i == (ix-1):
			latype = stsm0[0]
			if latype[0] < 1 or latype[0] > 3: return
			if latype[1] != 1: return
			latype = (newfam,1)
			smat = stsm0[1]
			stsms.append((latype,smat))
		else:
			stsms.append(stsm0)
	
	return SubAlgebraBrancher(ld0.latype, stsms, ld0.u1s)

def BrancherRenameB2C2(ld0, ix):
	if ix <= 0 or ix > len(ld0.stsms): return
	
	stsms = []
	for i,stsm0 in enumerate(ld0.stsms):
		if i == (ix-1):
			latype = stsm0[0]
			if latype[0] != 2 and latype[0] != 3: return
			if latype[1] != 2: return
			latype = (5-latype[0],2)
			smat = []
			for smrow in stsm0[1]:
				nwsmr = (smrow[1],smrow[0])
				smat.append(nwsmr)
			stsms.append((latype,smat))
		else:
			stsms.append(stsm0)
	
	return SubAlgebraBrancher(ld0.latype, stsms, ld0.u1s)


def BrancherRenameA3D3(ld0, ix):
	if ix <= 0 or ix > len(ld0.stsms): return
	
	stsms = []
	for i,stsm0 in enumerate(ld0.stsms):
		if i == (ix-1):
			latype = stsm0[0]
			if latype[0] != 1 and latype[0] != 4: return
			if latype[1] != 3: return
			latype = (5-latype[0],3)
			smat = []
			for smrow in stsm0[1]:
				nwsmr = (smrow[1],smrow[0],smrow[2])
				smat.append(nwsmr)
			stsms.append((latype,smat))
		else:
			stsms.append(stsm0)
	
	return SubAlgebraBrancher(ld0.latype, stsms, ld0.u1s)


def BrancherSplitD2(ld0, ix):
	if ix <= 0 or ix > len(ld0.stsms): return
	
	stsms = []
	for i,stsm0 in enumerate(ld0.stsms):
		if i == (ix-1):
			latype = stsm0[0]
			if latype[0] != 4: return
			if latype[1] != 2: return
			latype = (1,1)
			smat1 = []
			smat2 = []
			for smrow in stsm0[1]:
				smat1.append((smrow[0],))
				smat2.append((smrow[1],))
			stsms.append((latype,smat1))
			stsms.append((latype,smat2))
		else:
			stsms.append(stsm0)
	
	return SubAlgebraBrancher(ld0.latype, stsms, ld0.u1s)


def BrancherJoin2A1(ld0, ix1, ix2):
	if ix1 == ix2: return
	if ix1 <= 0 or ix1 > len(ld0.stsms): return
	if ix2 <= 0 or ix2 > len(ld0.stsms): return
	
	stsms = []
	smat1 = []
	smat2 = []
	for i,stsm0 in enumerate(ld0.stsms):
		if i == (ix1-1):
			latype = stsm0[0]
			if latype[1] != 1: return
			for smrow in stsm0[1]:
				smat1.append(smrow[0])
		elif i == (ix2-1):
			latype = stsm0[0]
			if latype[1] != 1: return
			for smrow in stsm0[1]:
				smat2.append(smrow[0])
		else:
			stsms.append(stsm0)
	smatc = zip(smat1,smat2)
	stsms.append(((4,2),smatc))
	
	return SubAlgebraBrancher(ld0.latype, stsms, ld0.u1s)

def BrancherConjugate(ld0, cjixs):
	stsms = []
	for i,stsm0 in enumerate(ld0.stsms):
		if (i+1) in cjixs:
			latype = stsm0[0]
			smat = []
			for smrow in stsm0[1]:
				if latype[0] == 1: # A(n)
					nwsmr = list(smrow)
					nwsmr.reverse()
				elif latype[0] == 4: # D(n)
					nwsmr = list(smrow[:-2]) + [smrow[-1],smrow[-2]]
				elif latype[0] == 5 and latype[1] == 6: # E(6)
					nwsmr = list(smrow[:-1])
					nwsmr.reverse()
					nwsmr.append(smrow[-1])
				else:
					nwsmr = smrow
				smat.append(nwsmr)
			stsms.append((latype,smat))
		else:
			stsms.append(stsm0)
	
	return SubAlgebraBrancher(ld0.latype, stsms, ld0.u1s)


def BrancherConjgD4(ld0, ix, newrts):
	if ix <= 0 or ix > len(ld0.stsms): return
	if len(newrts) != 4: return
	if newrts[1] != 2: return
	nwrsrt = list(newrts)
	nwrsrt.sort()
	for i, nr in enumerate(nwrsrt):
		if nr != (i+1): return
	
	stsms = []
	for i,stsm0 in enumerate(ld0.stsms):
		if i == (ix-1):
			latype = stsm0[0]
			if latype[0] != 4: return
			if latype[1] != 4: return
			smat = []
			for smrow in stsm0[1]:
				nwsmr = []
				for nwix in newrts:
					nwsmr.append(smrow[nwix-1])
				smat.append(nwsmr)
			stsms.append((latype,smat))
		else:
			stsms.append(stsm0)
	
	return SubAlgebraBrancher(ld0.latype, stsms, ld0.u1s)


def BrancherRearrange(ld0, neword):
	if len(neword) != len(ld0.stsms): return
	nwrsrt = list(neword)
	nwrsrt.sort()
	for i, nr in enumerate(nwrsrt):
		if nr != (i+1): return
	
	stsms = [ld0.stsms[ix-1] for ix in neword]
	
	return SubAlgebraBrancher(ld0.latype, stsms, ld0.u1s)


def SubalgSelf(latype):
	return SubAlgebraBrancher(latype, ((latype, identmat(latype[1])),), ())
