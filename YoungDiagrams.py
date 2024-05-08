#!python3
# 
# For handling SU(n) multiplets with Young diagrams.
# SU(n) is the Lie algebra A(n-1)
# SU(2) ~ SO(3) ~ Sp(4) -- 3D angular momentum
# SU(4) ~ SO(6)
# 
# A SU(n) irreducible representation (irrep) is defined by
# a highest-weight vector, which is a a set of (n-1) weights,
# each of which is a nonnegative integer.
# For SU(2), the highest weight is 2*(angular momentum).
# 
# There are four data types that this notebook works with.
# 
# * Expanded Young diagrams: list of lists of values, default Null
# * Contracted Young diagrams: list of row lengths
# * Variable-length weight vectors
# * Fixed-length weight vectors: length n-1 for SU(n) / A(n-1)
# 
# The code avoids trailing empty lists for the first one,
# and trailing zeros for the second and third ones.
# 
#
# ** Data Conversion **
# converts between all four data types
#
# Going to fixed-length weights:
# n = (weight dimension) + 1, for SU(n)
# Going to expanded Young:
# func = labeling function, with args (row #),(col #), counted from 0
#
# ContractYoung(ydgm) -- expanded to contracted Young
# ExpandYoung(ydgm,func) -- contracted to expanded Young
# CtYoungToVarWts(yrls) -- contracted Young to variable-length weights
# VarWtsToCtYoung(vwts) -- variable-length weights to contracted Young
# FixWeights(vwts,n) -- variable-length to fixed-length weights
# VarWeights(fwts) -- fixed-length to variable-length weights
# XpYoungToVarWts(ydgm) -- expanded Young to variable-length weights
# VarWtsToXpYoung(vwts,func) --  variable-length weights to expanded Young
# CtYoungToFxdWts(yrls,n) -- contracted Young to fixed-length weights
# FxdWtsToCtYoung(fwts) -- fixed-length weights to contracted Young
# XpYoungToFxdWts(ydgm,n) -- expanded Young to fixed-length weights
# FxdWtsToXpYoung(fwts,func) -- fixed-length weights to expanded Young
#
#
# ** Young-Diagram Operations **
# does relabeling
# and transposing of Young diagrams.
# 
# RelabelYoung(ydgm,func) -- relabels an expanded diagram
# applying labeling function func with args row, column
#
# Transposes a Young diagram to find its dual
# TransposeXpYoung(ydgm) -- expanded Young
# TransposeCtYoung(yrls) -- contracted Young
# FastTpCtYoung(yrls)  -- contracted young, fast
#
#
# ** Degeneracies: total dimensions for SU(n), SO(n), and Sp(n) **
#
# For number of weights w,
# SU(n): n = w + 1, SO(2n): n = 2w, SO(2n+1): n = 2w+1, Sp(2n): n = 2w
#
# From fixed weights:
# DegenSUnFxd(fwts) -- SU(n)
# DegenSO2nFxd(fwts) -- SO(2n)
# DegenSOnp1Fxd(fwts) -- SO(2n+1)
# DegenSpnFxd(fwts) -- Sp(n), usually written Sp(2n) with n = w
#
# For contracted Young diagrams:
# DegenSUnCtYoung(yrls,n)
#
# Polynomial coefficients, since Python doesn't do computer algebra
# (list of nmi, den) -- for n, product of (n + nmi) divided by den
# DegenFacsSUnCtYoung(yrls) -- SU(n)
# DegenFacsSOnCtYoung(yrls) -- SO(n)
# DegenFacsSpnCtYoung(yrls) -- Sp(n)
# 
# Uses the polynomial-coefficient lists
# DegenEval(facs, n)
#
#
# ** SU(n) Representation Products **
# Implements the Littlewood-Richardson rule
# for decomposing product representations into irreps.
# 
# SUnProductXpYoung(dgrm1,dgrm2) takes two expanded diagrams and
# creates a list of expanded product diagrams that are labeled
# by the algorithm.
# 0 = first diagram and 1,2,3,... = rows of second diagram.
# 
# SUnProductCtYoung(yrls1,yrls2) with contracted diagrams
# SUnProductCtYoungCntd(yrls1,yrls2) with a counted result:
# list of (count, diagram)
#
# SUnProductFxd(fwt1,fwt2) with fixed weights
# SUnProductFxdCntd(fwt1,fwt2) with a counted result
#
#
# ** SU(n) Vector-Representation Powers **
# These are useful
# for decomposing by symmetry.
# 
# VPFindMultDgrmList(p) decomposes power p of the vector or
# fundametnal representation of SU(n). It makes a list with form
# (multiplicity, contracted Young diagram), and they have the order
# (most symmetric) to (most antisymmetric).
#
# VPFindDgrmList(p) finds only the diagrams for power p.
#
# VPFindMultDgrmListSeq(p) finds a list of VPFindMultDgrmList()
# for powers from 1 to p.
#
# VPFindDgrmListSeq(p) finds a list of VPFindDgrmList()
# for powers from 1 to p.
#
#
# ** SU(n) Weyl-Orbit Dimensions ** 
#
# SUnWeylOrbitDimCtYoung(yrls,n) -- contrated Young
# SUnWeylOrbitDimFxd(fwts) -- fixed-length weights
# 
# ** SU(n) Representation W-Orbit Multiplicities **
# finds the multiplicity of each root in an irrep.
# The rep roots are grouped by Weyl orbits,
# which can be labeled by their highest or dominant weights.
# These multiplicities are the Kostka numbers, which come in Kostka matrices.
# 
# KostkaMatrix(p) finds the Kostka matrix for power p. Its rows and columns
# have the order of the Young diagrams found with VPFindDgrmList[p].
# Rows: indexed by irrep diagrams. Columns: indexed by Weyl-orbit diagrams.
#
# KostkaMatrices(p) finds a list of Kostka matrices for powers from 1 to p.
#
#
# ** SU(n) W-Orbit Products and Powers ** -- NOT VERIFIED
#
# SUnWeylOrbitProductCtYoung(yrls1,yrls2) takes two contracted Young diagrams,
# and finds a list of product Weyl orbits
# as a list of (multiplicity, contracted Young)
#
#
# ** SU(n) Branching Rules ** -- NOT VERIFIED
#
# SUnWeylOrbitBranching[n1,n2,maxwts) args n1 of SU(n1), n2 of SU(n2),
# maxwts of original algebra, SU(n1+n2).
# Makes list of SU(n1) highest weight, SU(n2) highest weight, U(1) factor.
# U(1) factor is multiplied by (n1+n2) to make it an integer.
# Handle SU(n) -> SU(n-1)*U(1) by treating it as SU(n-1)*SU(1)*U(1).
#
#
# ** Young-Diagram Nesting **
#
# YDDgrmNesting(p) finds the YD nesting matrix for power p.
# Rows and columns are indexed by the YD's for power p. For the row YD,
# take its row lengths and find the YD's for each of the length.
# Take one each and concatenate them to form a YD that is the column YD.
# Add to the matrix's entry at that location the indices of each of these YD's.
#
# YDDgrmNestings(p) finds a list of YD nesting matrices for powers from 1 to p.
#
#
# ** Schur Functions and Related Ones **
#
# *** Symmetric-Group Characters ***
#
# For checking some Schur-polynomial formulas.
# SymmGroupCharacters(p) up to power p.
# A list of matrices like KostkaMatrices.
# Rows are irreps, columns are conjugacy classes.
#
# *** Symmetric Functions ***
# *** Schur Functions ***
#
# Omitted because they are symbolic -- they are in the Mma version.
#
# *** Plethysms ***
#
# In Lie algebras, representation-power decompositions by symmetry type.
# The symmetry type is the plethysmer
# and the starting rep is the plethysmed one.
#
# SymYDsOfPower(pwr, sym, nmax) list of plethysms
# of (anti)symmetric YD's/Schur functions of power sums with power pwr.
# The YD's are all {p} to {nmax} for sym >=0  (symmetric)
# and all {1^p} to {1^nmax} for sym < 0 (antisymmetric).
# GenYDsOfPower(pwr, nmax) like above, but for general YD's,
# arranged like in VPFindDgrmListSeq(nmax)
#
# SelYDPlethysm(pthsel, pthnmax, sel, nmax) makes
# Rows of plethysmed ones designated with sel and nmax
# Columns of plethysmers designed with pthsel and pthmax
# The designators: selectors of which ones to emit,
# and maximum order / size / power
# The selectors:
# Sym -- the YD {nmax}
# SymAll -- the YD's {p} up to p = nmax
# Ats -- the YD {1^nmax}
# AtsAll-- the YD's {1^p} up to p = nmax
# Gen -- all the YD's with size nmax
# GenAll -- all the YD's with sizes up to nmax
# All these are one-level lists of counted lists of YD's
# and are {}'d or flattened as necessary
# Ordering: same as in VPFindDgrmListSeq[]
#
#
# Skipping the graphics, though one could do ASCII art.
# 
# Source for multiplicity and product algorithms:
# Particle Data Group:
# http://pdg.lbl.gov > Reviews > Mathematical Tools > SU(n) multiplets and Young diagrams
# 
# Also see
# http://en.wikipedia.org/wiki/Littlewood-Richardson_rule
# http://en.wikipedia.org/wiki/Kostka_number
# 
# Kostka-matrix algorithm from
# Determinantal Expression and Recursion for Jack Polynomials,
# by L. Lapointe, A. Lascoux, J. Morse
#
# Doing tuples as much as possible


###
### Data Conversion
###

# Expanded Young < - > contracted Young

def ContractYoung(ydgm): return tuple((len(yrow) for yrow in ydgm))

# The function's row and col args are 0-based, not 1-based
def ExpandYoung(yrls,func=(lambda r,c: None)): 
	return tuple(( tuple((func(r,c) for c in range(yrls[r]))) \
		for r in range(len(yrls))))

# Contracted Young < - > variable weights

def CtYoungToVarWts(yrls):
	yrlsx = list(yrls) + [0]
	return tuple((yrlsx[i]-yrlsx[i+1] for i in range(len(yrls))))

def VarWtsToCtYoung(vwts):
	yrls = len(vwts)*[0]
	for i in range(len(vwts)):
		for j in range(i+1):
			yrls[j] += vwts[i]
	return tuple(yrls)

# Variable-length weights < - > fixed-length weights

def FixWeights(vwts,n):
	wl = len(vwts)
	if wl < n:
		return tuple(list(vwts) + (n - 1 - wl)*[0])
	elif wl == n:
		return tuple(vwts[:-1])
	else:
		return None

def VarWeights(fwts):
	if len(fwts) > 0:
		if fwts[-1] == 0:
			return VarWeights(fwts[:-1])
		else:
			return tuple(fwts)
	else:
		return ()

# Expanded Young < - > variable-length weights

def XpYoungToVarWts(ydgm):
	return CtYoungToVarWts(ContractYoung(ydgm))

def VarWtsToXpYoung(vwts,func=(lambda r,c: None)):
	return ExpandYoung(VarWtsToCtYoung(vwts),func)

# Contracted Young < - > fixed-length weights

def CtYoungToFxdWts(yrls,n):
	return FixWeights(CtYoungToVarWts(yrls),n)

def FxdWtsToCtYoung(fwts):
	return VarWtsToCtYoung(VarWeights(fwts))

# Expanded Young < - > fixed-length weights

def XpYoungToFxdWts(ydgm,n):
	return FixWeights(CtYoungToVarWts(ContractYoung(ydgm)),n)

def FxdWtsToXpYoung(fwts,func=(lambda r,c: None)):
	return ExpandYoung(VarWtsToCtYoung(VarWeights(fwts)),func)

###
### Young-Diagram Operations
###

def RelabelYoung(ydgm,func=(lambda r,c: None)):
	return tuple(( tuple((func(r,c) for c in range(len(ydgm[r])))) \
		for r in range(len(ydgm))))

def TransposeXpYoung(ydgm):
	n1 = len(ydgm)
	if n1 == 0: return ()
	n2 = max(map(len,ydgm))
	newydgm = [[] for k in range(n2)]
	for i1 in range(n1):
		for i2 in range(len(ydgm[i1])):
			newydgm[i2].append(ydgm[i1][i2])
	return tuple(map(tuple,newydgm))

def TransposeCtYoung(yrls):
	return ContractYoung(TransposeXpYoung(ExpandYoung(yrls)))

# Fast version for the contracted case

def FastTpCtYoung(yrls):
	nyrs = len(yrls)
	yrlsx = list(yrls) + [0]
	yrlsx.reverse()
	yrdf = [yrlsx[i+1] - yrlsx[i] for i in range(nyrs)]
	tpyr = []
	for i in range(nyrs):
		tpyr += yrdf[i] * [nyrs-i]
	return tuple(tpyr)

###
### SU(n) SO(2n) SO(2n+1) Sp(2n) Degeneracy: Total Dimensions
###

# Fixed-length weight vectors

def DegenSUnFxd(fwts):
	wlen = len(fwts)
	dnum = 1
	dden = 1
	for k in range(wlen):
		dfac = k + 1
		for l in range(wlen-k):
			dnum *= sum(fwts[l:l+k+1]) + dfac
			dden *= dfac
	return dnum//dden

def DegenSO2nFxd(fwts):
	wlen = len(fwts)
	dnum = DegenSUnFxd(fwts[:-1])
	dden = 1
	for k in range(wlen-1):
		dfac = k + 1
		dnum *= sum(fwts[-k-2:-2]) + fwts[-1] + dfac
		dden *= dfac
		for l in range(1,wlen-k-1):
			dfac = 2*k + l + 2
			dnum *= sum(fwts[-k-l-2:-k-2]) + 2*sum(fwts[-k-2:-2]) + \
				fwts[-2] + fwts[-1] + dfac
			dden *= dfac
	return dnum//dden

def DegenSOnp1Fxd(fwts):
	wlen = len(fwts)
	dnum = DegenSUnFxd(fwts)
	dden = 1
	for k in range(1,wlen):
		for l in range(wlen-k):
			dfac = 2*k + l + 1
			dnum *= sum(fwts[-k-l-1:-k-1]) + 2*sum(fwts[-k-1:-1]) + \
				fwts[-1] + dfac
			dden *= dfac
	return dnum//dden

def DegenSpnFxd(fwts):
	wlen = len(fwts)
	dnum = DegenSUnFxd(fwts)
	dden = 1
	for k in range(wlen):
		for l in range(1,wlen-k):
			dfac = 2*k + l + 2
			dnum *= sum(fwts[-k-l-1:-k-1]) + 2*sum(fwts[-k-1:-1]) + \
				2*fwts[-1] + dfac
			dden *= dfac
	return dnum//dden

# Contracted Young diagrams

def HookCtYoung(yrls):
	dual = TransposeCtYoung(yrls)
	hook = 1
	for i, yr in enumerate(yrls):
		for j in range(yr):
			hook *= 1 + (yrls[i] - i - 1) + (dual[j] - j - 1)
	return hook

def DegenSUnCtYoung(yrls,n):
	numer = 1
	denom = HookCtYoung(yrls)
	for i, yr in enumerate(yrls):
		for j in range(yr):
			numer *= n - i + j
	return numer//denom

def DegenFacsSUnCtYoung(yrls):
	numlist = []
	denom = HookCtYoung(yrls)
	for i, yr in enumerate(yrls):
		for j in range(yr):
			numlist.append(-i+j)
	numlist.sort()
	return (tuple(numlist), denom)

def DegenFacsSOSpnCtYoung(yrls,q):
	yrllen = len(yrls)
	yrx = lambda i: yrls[i] if i >= 0 and i < yrllen else 0
	numlist = []
	denlist = []
	denom = HookCtYoung(yrls)
	for i in range(yrllen):
		for j in range(yrls[i]):
			numlist.append(-2*i-j+q-2+yrx(i)+yrx(i+j+q))
		for j in range(0,yrllen-i-yrx(i)-1):
			numlist.append(-2*i-j+q-2+yrx(i+yrx(i)+j+q))
			denlist.append(-2*i-j+q-2)
	# Remove all members of denlist from numlist
	numlist.sort()
	for d in denlist:
		numlist.remove(d)
	return (tuple(numlist), denom)
	

def DegenFacsSOnCtYoung(yrls): return DegenFacsSOSpnCtYoung(yrls,0)
def DegenFacsSpnCtYoung(yrls): return DegenFacsSOSpnCtYoung(yrls,1)

def DegenEval(facs, n):
	numlist, denom = facs
	numer = 1
	for nf in numlist:
		numer *= n + nf
	return numer//denom

###
### SU(n) Representation Products
###

# Product of reps; uses the Littlewood-Richardson rule.

# Auxiliary functions. These work with expanded diagrams
# They assume lists, not tuples
# Keep in mind 0-based and not 1-based indexing

def ClipEmpty(dgrm):
	newdgrm = dgrm
	for k, drow in enumerate(dgrm):
		if len(drow) == 0:
			newdgrm = dgrm[:k]
			break
	return newdgrm

def SUnAddOneToYoung(dgrm, token, start):
	# token must be absent from the original diagram
	olddgrm = dgrm + [[]]
	outlist = []
	for k in range(start,len(olddgrm)):
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
			for l in range(k):
				if newdgrm[l][ndlen-1] == token:
					apnd = False
					break
		if apnd:
			outlist.append((newdgrm,k))
	
	outlist = [(ClipEmpty(outmem[0]), outmem[1]) for outmem in outlist]
	return outlist

def SUnAddSetToYoung(dgrm,token,number):
	# token must be absent from the original diagram
	outlist = [(dgrm,0)]
	oldlist = outlist[:]
	for k in range(number):
		outlist = []
		for oldmem in oldlist:
			outlist += SUnAddOneToYoung(oldmem[0],token,oldmem[1])
		oldlist = outlist[:]
	
	if len(outlist) > 0:
		outlist = [outmem[0] for outmem in outlist]
	return outlist
	
# Works with expanded diagrams; returns a set of expanded diagrams
# with boxes labeled as a result of the algorithm.
# 0 is first diagram and 1 2 3 are the rows of the second diagram.
def SUnProductXpYoung(dgrm1,dgrm2):
	# Be sure to keep all the tokens distinct;
	# 0 for the original diagram is distinct from
	# 1, 2, ... for the added one.
	fxddgrm = RelabelYoung(ClipEmpty(dgrm1), lambda r,c: 0)
	
	# Be sure that the the original diagram is all-list
	oldlist = [[list(fxddgrow) for fxddgrow in fxddgrm]]
	for k,dgrow2 in enumerate(dgrm2):
		outlist = []
		for oldmem in oldlist:
			outlist += SUnAddSetToYoung(oldmem, k+1, len(dgrow2))
		oldlist = outlist[:]

	outlist = []
	for dgx in oldlist:
		cnt = len(dgrm2)*[0]
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
	return tuple((tuple((tuple(dgrow) for dgrow in dgx)) for dgx in outlist))

# Like the above, but with contracted diagrams
def SUnProductCtYoung(yrls1,yrls2):
	dgrm1 = ExpandYoung(yrls1)
	dgrm2 = ExpandYoung(yrls2)
	prod = SUnProductXpYoung(dgrm1,dgrm2)
	return tuple((ContractYoung(dgrm) for dgrm in prod))

def CountedSet(lst):
	cnt = {}
	for mem in lst:
		if mem not in cnt: cnt[mem] = 0
		cnt[mem] += 1
	return tuple(((cnt[mem],mem) for mem in sorted(cnt.keys())))

def SUnProductCtYoungCntd(yrls1,yrls2):
	return CountedSet(SUnProductCtYoung(yrls1,yrls2))

def SUnProductFxd(fwt1,fwt2):
	nf = len(fwt1)
	if len(fwt2) != nf: return None
	ydgm1 = FxdWtsToXpYoung(fwt1)
	ydgm2 = FxdWtsToXpYoung(fwt2)
	res = SUnProductXpYoung(ydgm1,ydgm2)
	return tuple((XpYoungToFxdWts(dgrm,nf+1) for dgrm in res if len(dgrm) <= nf+1))

def SUnProductFxdCntd(fwt1,fwt2):
	return CountedSet(SUnProductFxd(fwt1,fwt2))

###
### SU(n) Vector-Representation Powers
###

def VPYoungAddOne(ylst):
	ylnl = []
	ylnew = list(ylst)
	ylnew[0] += 1
	ylnl.append(tuple(ylnew))
	for i in range(len(ylst)-1):
		if ylst[i] > ylst[i+1]:
			ylnew = list(ylst)
			ylnew[i+1] += 1
			ylnl.append(tuple(ylnew))
	ylnew = list(ylst)
	ylnew.append(1)
	ylnl.append(tuple(ylnew))
	return tuple(ylnl)

def VPYoungSubtractOne(ylst):
	ylnl = []
	for i in range(len(ylst)-1):
		if ylst[i] > ylst[i+1]:
			ylnew = list(ylst)
			ylnew[i] -= 1
			ylnl.append(tuple(ylnew))
	ylnew = list(ylst)
	ylnew[-1] -= 1
	if ylnew[-1] == 0: del ylnew[-1]
	ylnl.append(tuple(ylnew))
	return tuple(ylnl)
	
# Find next tensor-product diagram list from a tensor-product list of
# (multiplicity, contracted Young diagram)
# For inputs in a standard ordering of same-box-number Young diagrams,
# it produces outputs in that order.
def VPFindNextDgrmList(dgrmlist):
	nxdgs = []
	dgcnt = {}
	for num,dg in dgrmlist:
		ndgs = VPYoungAddOne(dg)
		for ndg in ndgs:
			if ndg not in dgcnt:
				nxdgs.append(ndg)
				dgcnt[ndg] = 0
			dgcnt[ndg] += num
	nxdgmlst = [(dgcnt[ndg],ndg) for ndg in nxdgs]
	return tuple(nxdgmlst)

# Finds all those up to power p.
VPMultDgrmListCache = []

def VPFindMultDgrmList(p):
	if p == 0: return ((1,()),)
	
	if p > len(VPMultDgrmListCache):
		if len(VPMultDgrmListCache) < 1:
			VPMultDgrmListCache.append(((1,(1,)),))
		
		numrem = p - len(VPMultDgrmListCache)
		for i in range(numrem):
			nxtlst = VPFindNextDgrmList(VPMultDgrmListCache[-1])
			VPMultDgrmListCache.append(nxtlst)
	
	return VPMultDgrmListCache[p-1]

def VPFindMultDgrmListSeq(p):
	return tuple((VPFindMultDgrmList(ip) for ip in range(1,p+1)))

def VPFindDgrmList(p):
	return tuple((ndg[1] for ndg in VPFindMultDgrmList(p)))

def VPFindDgrmListSeq(p):
	return tuple((VPFindDgrmList(ip) for ip in range(1,p+1)))

###
### SU(n) Weyl-Orbit Dimensions
###

# Denominator for count of permutations
def CountedSetFCTPRod(yrls):
	cnt = {}
	for yrl in yrls:
		if yrl not in cnt: cnt[yrl] = 0
		cnt[yrl] += 1
	prod = 1
	for n in cnt.values():
		for k in range(1,n+1,1):
			prod *= k
	return prod

# Dimension of Weyl orbit: contracted Young and fixed-length weights

def SUnWeylOrbitDimCtYoung(yrls,n):
	prod = 1
	for k in range(len(yrls)):
		prod *= (n-k)
	return prod//CountedSetFCTPRod(yrls)

def SUnWeylOrbitDimFxd(fwts):
	return SUnWeylOrbitDimCtYoung(FxdWtsToCtYoung(fwts),len(fwts)+1)

###
### SU(n) Representation W-Orbit Multiplicities
###

def KostkaMatrix(p):
	dgrms = VPFindDgrmList(p)
	# Index the diagrams
	dglen = len(dgrms)
	dgix = {}
	for i,dg in enumerate(dgrms):
		dgix[dg] = i
	# Find "C" lookback matrix and "g" denominator factor
	# Multiply by 2 to avoid fractions
	lookback = [dglen*[0] for i in range(dglen)]
	denfac = dglen*[0]
	for i,dg in enumerate(dgrms):
		dgl = len(dg)
		# Do lookback; the += takes care of the multiplicities automatically.
		for j in range(dgl-1):
			for jx in range(j+1,dgl):
				for k in range(1,dg[jx]+1):
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
	kostka = [dglen*[0] for i in range(dglen)]
	for i in range(dglen):
		kkrow = kostka[i]
		kkrow[i] = 1
		for j in range(i+1,dglen):
			kknum = 0
			lbrow = lookback[j]
			for k in range(i,j):
				kknum += lbrow[k]*kkrow[k]
			if kknum != 0:
				kkrow[j] = kknum/(denfac[i] - denfac[j])
	return tuple(map(tuple,kostka))

def KostkaMatrices(p):
	return tuple((KostkaMatrix(ip) for ip in range(1,p+1)))

###
### SU(n) W-Orbit Products
###

# All possible splits of a list. Can have repeats in case we want to count them. 

def ListSplits(lst):
	n = len(lst)
	slsts = []
	# Go through all the possible selection combinations. 
	# Use an algorithm that can easily be ported. Start off with all zeros
	dgts = n*[0]
	while True:
		slst1 = []
		slst2 = []
		for k,dgt in enumerate(dgts):
			if dgt != 0:
				slst1.append(lst[k])
			else:
				slst2.append(lst[k])
		slst12 = (tuple(slst1),tuple(slst2))
		slsts.append(slst12)
		#  Move to next one, and quit when one has done all of them.
		pos = 0
		while pos < n:
			if dgts[pos] == 1:
				dgts[pos] = 0
				pos += 1
			else:
				dgts[pos] = 1
				break
			if pos >= n: break
		if pos >= n: break
	return tuple(slsts)

# Like above, but unique splits

def UniqueListSplits(lst):
	n = len(lst)
	slsts = set()
	# Go through all the possible selection combinations. 
	# Use an algorithm that can easily be ported. Start off with all zeros
	dgts = n*[0]
	while True:
		slst1 = []
		slst2 = []
		for k,dgt in enumerate(dgts):
			if dgt != 0:
				slst1.append(lst[k])
			else:
				slst2.append(lst[k])
		slst12 = (tuple(slst1),tuple(slst2))
		slsts.add(slst12)
		#  Move to next one, and quit when one has done all of them.
		pos = 0
		while pos < n:
			if dgts[pos] == 1:
				dgts[pos] = 0
				pos += 1
			else:
				dgts[pos] = 1
				break
			if pos >= n: break
		if pos >= n: break
	slsts = list(slsts)
	slsts.sort()
	return tuple(slsts)

# Like above, but finding counts as one goes

def CountedListSplits(lst):
	n = len(lst)
	slstcnts = {}
	# Go through all the possible selection combinations. 
	# Use an algorithm that can easily be ported. Start off with all zeros
	dgts = n*[0]
	while True:
		slst1 = []
		slst2 = []
		for k,dgt in enumerate(dgts):
			if dgt != 0:
				slst1.append(lst[k])
			else:
				slst2.append(lst[k])
		slst12 = (tuple(slst1),tuple(slst2))
		if slst12 not in slstcnts: slstcnts[slst12] = 0
		slstcnts[slst12] += 1
		#  Move to next one, and quit when one has done all of them.
		pos = 0
		while pos < n:
			if dgts[pos] == 1:
				dgts[pos] = 0
				pos += 1
			else:
				dgts[pos] = 1
				break
			if pos >= n: break
		if pos >= n: break
	
	slsts = slstcnts.keys()
	slsts.sort()
	slsts = [(slstcnts[slst],slst[0],slst[1]) for slst in slsts]
	return tuple(slsts)

# Implementation of a Mathematica function
def permutations(lst):
	n = len(lst)
	# The zero and one cases are very simple
	if n <= 1: return set([lst])
	# For more, do it recursively. Stack depth: about n
	plst = set()
	for k in range(n):
		ps = permutations(lst[:k] + lst[k+1:])
		lmlst = lst[k:k+1]
		plst = plst.union(map(lambda lm: lmlst+lm, ps))
	return plst

def SUnWeylOrbitProductCtYoung(yrls1,yrls2):
	ys1 = tuple(yrls1)
	ys2 = tuple(yrls2)
	n1 = len(ys1)
	n2 = len(ys2)
	n0 = min(n1,n2)
	lss1 = UniqueListSplits(ys1)
	lss2 = UniqueListSplits(ys2)
	lss1 = filter(lambda lm: len(lm[0]) <= n0, lss1)
	lss2 = filter(lambda lm: len(lm[0]) <= n0, lss2)
	lsrs = set()
	for ls11,ls12 in lss1:
		ps = permutations(ls11)
		for ls21,ls22 in lss2:
			if len(ls21) == len(ls11):
				px = ls21
				for p in ps:
					pz = list(zip(p,px))
					pz.sort()
					lsrs.add((tuple(pz),ls12,ls22))
	
	lsrres = {}
	for lsr0,lsr1,lsr2 in lsrs:
		lsr0s = [r0+r1 for r0,r1 in lsr0]
		lsrnew = lsr0s
		lsrnew += list(lsr1)
		lsrnew += list(lsr2)
		lsrnew.sort()
		lsrnew.reverse()
		lsrnew = tuple(lsrnew)
		num = CountedSetFCTPRod(lsrnew)
		num /= CountedSetFCTPRod(lsr0s)
		num /= CountedSetFCTPRod(lsr1)
		num /= CountedSetFCTPRod(lsr2)
		num = int(num)
		if lsrnew not in lsrres:
			lsrres[lsrnew] = 0
		lsrres[lsrnew] += num
	lsrlst = lsrres.keys()
	lsrout = [(lsrres[lsr], lsr) for lsr in sorted(lsrlst,reverse=True)]
	return tuple(lsrout)

# Test Weyl-orbit products
def SUnWeylOrbitProductCtYoungVerify(yrls1,yrls2,n):
	prod = SUnWeylOrbitProductCtYoung(yrls1,yrls2)
	sums1 = [SUnWeylOrbitDimCtYoung(yr,n) for yr in (yrls1,yrls2)]
	sums2 = [num*SUnWeylOrbitDimCtYoung(yr,n) for num,yr in prod]
	return sums1[0]*sums1[1] - sum(sums2)

# Decomposes the square into symmetric and antisymmetric parts

def SUnWeylOrbitSquare(yrls):
	prod = SUnWeylOrbitProductCtYoung(yrls,yrls)
	yr2 = tuple((2*yr for yr in yrls))
	sym = []
	ats = []
	for pe in prod:
		if pe[1] == yr2:
			sym.append(pe)
		else:
			pex = (pe[0]/2,pe[1])
			sym.append(pex)
			ats.append(pex)
	return (tuple(sym),tuple(ats))

###
### SU(n) Branching Rules
###

# Weyl-orbit branching rules for demoting a root of SU (n1+n2)
# and splitting it into SU(n1)*SU(n2)*U(1).
#
# The SU(n) to SU(n-1)*U(1) case is implemented
# by handling a split into SU(n-1)*SU(1)*U(1).
#
# The weight arg and outputs are appropriate fixed weights.
# Output is (SU(n1) weight) (SU(n2) weight) U(1) factor.

def SUnWeylOrbitBranching(n1,n2,maxwts):
	lspl = UniqueListSplits(FxdWtsToCtYoung(maxwts))
	lspl = [(CtYoungToFxdWts(lm[0],n1), CtYoungToFxdWts(lm[1],n2), \
		n2*sum(lm[0]) - n1*sum(lm[1])) for lm in lspl]
	lspl = filter(lambda lm: lm[0] != None and lm[1] != None, lspl)
	return tuple(lspl)

###
### Young-Diagram Nesting
###

# Imitation of Mathematica Tuples function for a list of lists:
def multlstsel(multlst):
	res = [[]]
	for lst in multlst:
		newres = []
		for rsmem in res:
			for mem in lst:
				newres.append(rsmem + [mem])
		res = newres
	return res

def YDDgrmNesting(p):
	dgrms = VPFindDgrmList(p)
	dglen = len(dgrms)
	dgix = {}
	for i,dgrm in enumerate(dgrms):
		dgix[dgrm] = i
	
	nesting = [[[] for j in range(dglen)] for i in range(dglen)]
	for i,dgrm in enumerate(dgrms):
		subdgrms = list(map(VPFindDgrmList,dgrm))
		dgxp = multlstsel([range(len(dm)) for dm in subdgrms])	
		for dgxpm in dgxp:
			insdg = []
			for l,dgxpmx in enumerate(dgxpm):
				insdg += subdgrms[l][dgxpmx]
			insdg.sort()
			insdg.reverse()
			insdg = tuple(insdg)
			j = dgix[insdg]
			nesting[i][j].append(tuple(dgxpm))
	return tuple(map(lambda m: tuple(map(tuple,m)) ,nesting))

def YDDgrmNestings(p):
	return tuple((YDDgrmNesting(ip) for ip in range(1,p+1)))

###
### Schur Functions and Related Ones
###

### Symmetric-Group Characters

def trimzeros(x): return tuple((xm for xm in x if xm > 0))

# If "verify" is set on, then it will remove every row in a conjugacy class
# and compare the results, instead of removing only one
def SymmGroupCharacters(p, verify=False):
	vpdset = VPFindDgrmListSeq(p)
	sgchset = []
	IsCorrect = True
	# Index the diagrams -- which power, where in the power.
	vpdix = {}
	for i,vpds in enumerate(vpdset):
		for j,vpd in enumerate(vpds):
			vpdix[vpd] = (i,j)
	# The empty diagram
	emptyix = (-1,-1)
	vpdix[()] = emptyix
	
	# Use the Murnaghan-Nakayama rule
	for i,vpds in enumerate(vpdset):
		# Rows: irreps. Find all the hook removals for this irrep.
		# Find all the hook removals for this irrep. 
		# The dehooked set dhset is indexed by hook length. 
		# Each member of sgch is a list of (power, dgram index, hook sign)
		sgch = []
		for j,vpdirp in enumerate(vpds):
			vptrirp = FastTpCtYoung(vpdirp)
			dhset = [[] for ix in range(i+1)]
			
			for k in range(len(vpdirp)):
				for l in range(vpdirp[k]):
					# Hook sign and length
					hksgn = (-1)**(vptrirp[l]-k-1)
					hklen = (vpdirp[k]-l-1) + (vptrirp[l]-k-1) + 1
					# The top and left of the dehooked diagram
					dehooked = list(vpdirp)
					for m in range(k,vptrirp[l]):
						dehooked[m] = min(dehooked[m],l)
					# The bottom and right of the dehooked diagram
					dhkxtra = [vpdirp[m]-l-1 for m in range(k+1,vptrirp[l])]
					# Slide the bottom right into the top left
					for m in range(len(dhkxtra)):
						dehooked[k+m] += dhkxtra[m]
					# Finally...
					dehooked = trimzeros(dehooked)
					dhsval = tuple(list(vpdix[dehooked]) + [hksgn])
					dhset[hklen-1].append(dhsval)
			
			# Calculate the characters using this info
			sgchrow = []
			for k,vpdcls in enumerate(vpds):
				vdimax = len(vpdcls) if verify else min(len(vpdcls),1)
				# Make a list of results for all the deletions
				chtrial = []
				for l in range(vdimax):
					rwlen = vpdcls[l]
					vpclred = tuple(list(vpdcls[:l]) + list(vpdcls[l+1:]))
					vpcrix = vpdix[vpclred]
					# Sum over all the appropriate-length hooks
					dhs = dhset[rwlen-1]
					chtsm = 0
					for dhkix in dhs:
						if dhkix[0] == -1:
							dvl = 1
						else:
							dvl = sgchset[dhkix[0]][dhkix[1]][vpcrix[1]]
						chtsm += dhkix[2]*dvl
					chtrial.append(chtsm)
				
				if verify and len(chtrial) > 1:
					AllEqual = True
					for l in range(len(chtrial)-1):
						if chtrial[l] != chtrial[l+1]:
							AllEqual = False
							IsCorrect = False
							break
					if not AllEqual:
						print(i, j, k, "--", chtrial)
				
				chtrval = chtrial[0] if len(chtrial) > 0 else 0
				sgchrow.append(chtrval)
			sgch.append(tuple(sgchrow))
		sgchset.append(tuple(sgch))
	sgchset = tuple(sgchset)
	if verify: sgchset = (sgchset,IsCorrect)
	return sgchset

### Plethysms

# Power sum to sum of Schur functions, expressed to YD's.

def PowerSumToYDs(n):
	return tuple((((-1)**k, tuple([n-k] + k*[1])) for k in range(n)))

# Caching of YD products:

CacheYDProduct = {}

def GetYDProduct(yd1,yd2):
	prodkey = [yd1,yd2]
	prodkey.sort()
	prodkey = tuple(prodkey)
	if prodkey not in CacheYDProduct:
		CacheYDProduct[prodkey] = SUnProductCtYoungCntd(yd1,yd2)
	return CacheYDProduct[prodkey]

# Operations on counted lists of YD's

def CountedYDListZero(): return ()

def CountedYDListScalMult(sclr,cydl):
	return tuple(((sclr*cyd[0],cyd[1]) for cyd in cydl)) if sclr != 0 \
		else CountedYDListZero()

def CountedYDListScalDiv(sclr,cydl):
	return tuple(((cyd[0]/sclr,cyd[1]) for cyd in cydl)) if sclr != 0 \
		else CountedYDListZero()

def CountedListMake(vpdix):
	yds = vpdix.keys()
	return tuple(((vpdix[yd],yd) for yd in sorted(yds,reverse=True) if vpdix[yd] != 0))

def CountedYDListAddList(cydllist):
	vpdix = {}
	for cydl in cydllist:
		for cyd in cydl:
			if cyd[1] not in vpdix:
				vpdix[cyd[1]] = 0
			vpdix[cyd[1]] += cyd[0]
	return CountedListMake(vpdix)

def CountedYDListAdd(cydl1,cydl2): return CountedYDListAddList((cydl1,cydl2))

def CountedYDListMult(cydl1,cydl2):
	vpdix = {}
	for cyd1 in cydl1:
		for cyd2 in cydl2:
			cydp = GetYDProduct(cyd1[1],cyd2[1])
			cydp = CountedYDListScalMult(cyd1[0]*cyd2[0],cydp)
			for cyd in cydp:
				if cyd[1] not in vpdix:
					vpdix[cyd[1]] = 0
				vpdix[cyd[1]] += cyd[0]
	return CountedListMake(vpdix)

# Functions of (anti)symmetric Young diagrams,
# found from list of functions of powers
def SymYDsFromPowers(sym,PwrFuncList):
	# The operations
	ZeroOp = CountedYDListZero
	ScalMultOp = CountedYDListScalMult
	ScalDivOp = CountedYDListScalDiv
	AddOp = CountedYDListAdd
	AddListOp = CountedYDListAddList
	MultOp = CountedYDListMult
	
	nmax = len(PwrFuncList)
	splst = []
	for n in range(1,nmax+1):
		spdlist = []
		for p in range(1,n):
			prod = MultOp(PwrFuncList[p-1],splst[n-p-1])
			if sym < 0:
				prod = ScalMultOp((-1)**(p-1),prod)
			spdlist.append(prod)
		prod = PwrFuncList[n-1]
		if sym < 0:
			prod = ScalMultOp((-1)**(n-1),prod)
		spdlist.append(prod)
		spd = AddListOp(spdlist)
		spd = ScalDivOp(n,spd)
		splst.append(spd)
	return tuple(splst)

# The YD's corresponding to the Schur expansion
# of a (anti)symmetric Schur function of a power of a vector --
# plethysm of (anti)symmetric Schur function on a power sum.
# Args are that power, the symmetry, and the maximum symmetric-diagram order.
# Coefficients are all +1 0 -1
def SymYDsOfPower(pwr,sym,nmax):
	PFL = [PowerSumToYDs(pwr*p) for p in range(1,nmax+1)]
	return SymYDsFromPowers(sym,PFL)


# Chen's algorithm for these ones. For symmetric ones, construct all the YD's
# with rim hooks with length pwr and with at least one box in the top row.
# The multiplier is the product of the signs of all these rim hooks.
# For antisymmetric ones, likewise but left column instead of top row.
#
# Y.M. Chen, "Combinatorial Algorithms for Plethysm", PhD Thesis,
# University of California of San Diego, 1982

def AddRHToDgrm(dgrm,hlen):
	res = []
	itdgrm = list(dgrm)
	rhil = 0
	sgn = 1
	for i in range(len(dgrm)):
		if i == 0:
			itdgrm[i] += 1
			rhil += 1
		else:
			itdgrm[i] = dgrm[i-1] + 1
			rhil += (itdgrm[i] - dgrm[i])
		if rhil > hlen: break
		nwdgrm = list(itdgrm)
		nwdgrm[0] += (hlen-rhil)
		nwsgn = sgn
		sgn *= -1
		res.append((nwsgn,tuple(nwdgrm)))
	return res

def AddRHToDcls(dcls,hlen):
	res = [CountedYDListScalMult(dcl[0],AddRHToDgrm(dcl[1],hlen)) for dcl in dcls]
	return CountedYDListAddList(res)

def RHSymYDsOfPower(pwr,sym,nmax):
	res = []
	if nmax >= 1:
		res.append(PowerSumToYDs(pwr))
		for i in range(1,nmax):
			res.append(AddRHToDcls(res[-1],pwr))
	
	if sym < 0:
		resnew = []
		for k,rs in enumerate(res):
			rsnw = [(r[0],FastTpCtYoung(r[1])) for r in rs]
			rsnw = list(CountedYDListScalMult((-1)**((pwr-1)*(k+1)), rsnw))
			rsnw.reverse()
			resnew.append(tuple(rsnw))
		res = resnew
	
	return tuple(res)

# Find general YD's of a power -- plethysm of general Schur/YD on a power sum

# Add vectors to make a matrix, and accumulate its inverse while doing so
# Reject vectors that are linearly dependent on those already added

def GCD(a,b):
	# Return zero, make positive
	if a == 0: return 0
	ax = a if a > 0 else -a
	if b == 0: return 0
	bx = b if b > 0 else -b
	
	# Euclid's algorithm
	while bx > 0:
		cx = ax % bx
		ax = bx
		bx = cx
	return ax

def GCDList(nlst):
	d = 0
	for n in nlst:
		if n != 0:
			if d == 0:
				d = n
			else:
				d = GCD(d,n)
	return d

def IndexFirstNonzero(vec):
	ix = -1
	for i,x in enumerate(vec):
		if x != 0:
			ix = i
			break
	return ix

# Returns True if the new vector is linearly independent
# and is successfully added, otherwise returns False
def AddVectorToMatrixIfLinIndep(mat,mtixs,vec):
	# Add the appropriate identity-matrix row at the end
	nvc = list(vec)
	vlen = len(vec)
	idm = vlen*[0]
	idm[len(mat)] = 1
	nvc += idm
	
	# Subtract out the previous ones,
	# scale if necessary to keep it all-integer
	newmat = map(list,mat)
	for i,nmln in enumerate(newmat):
		ix = mtixs[i]
		x0 = nmln[ix]
		x1 = nvc[ix]
		xgcd = GCD(x0,x1)
		if xgcd != 0:
			x0 /= xgcd
			x1 /= xgcd
			for j in range(2*vlen):
				nvc[j] = x0*nvc[j] - x1*nmln[j]
		
	# Accept it?
	ix = IndexFirstNonzero(nvc[:vlen])
	if ix < 0: return False
	
	# Do backsubstitution
	for i,nmln in enumerate(newmat):
		x0 = nmln[ix]
		x1 = nvc[ix]
		xgcd = GCD(x0,x1)
		if xgcd != 0:
			x0 /= xgcd
			x1 /= xgcd
			for j in range(2*vlen):
				nmln[j] = x1*nmln[j] - x0*nvc[j]
	
	# Append the new row
	for i,nmln in enumerate(newmat):
		mat[i] = nmln
	mat.append(nvc)
	mtixs.append(ix)
	
	# Normalize by dividing by the gcd
	for mln in mat:
		g = GCDList(mln)
		for i in range(len(mln)):
			mln[i] /= g
		
	return True

def FindScaledInverse(mat,mtixs,vlen):
	midg = vlen*[0]
	miscld = [vlen*[0] for i in range(vlen)]
	for i in range(vlen):
		midg[mtixs[i]] = mat[i][mtixs[i]]
		for j in range(vlen):
			miscld[mtixs[i]][j] = mat[i][vlen+j]
	
	midg = tuple(midg)
	miscld = tuple(map(tuple,miscld))
	return (midg,miscld)


# First, symmetric*general products - symmetric to general coefficients.
# Assume symmetric is first in the product irreps found.

def YDSymGenProdsToGenYDs(pwr):
	dgrmlist = VPFindDgrmListSeq(pwr)
	dgrms = dgrmlist[pwr-1]
	ndgs = len(dgrms)
	
	dgix = {}
	for i,dgrm in enumerate(dgrms):
		dgix[dgrm] = i
	bkdgix = {}
	for i in range(pwr):
		for j,dgrm in enumerate(dgrmlist[i]):
			bkdgix[dgrm] = (i,j)
	
	prods = []
	symctrbs = []
	xfrmat = []
	xfmtixs = []
	for p,dgrms in enumerate(dgrmlist[:-1]):
		for dgrm in dgrms:
			prod = GetYDProduct((pwr-p-1,),dgrm)
			prdx = ndgs*[0]
			for dgx in prod:
				prdx[dgix[dgx[1]]] = dgx[0]
			smcf = prdx[0]
			prdx = prdx[1:]
			rc = AddVectorToMatrixIfLinIndep(xfrmat,xfmtixs,prdx)
			if rc:
				prods.append((pwr-p-1,bkdgix[dgrm]))
				symctrbs.append(smcf)
	
	invmat = FindScaledInverse(xfrmat,xfmtixs,len(xfrmat[0])//2) \
		if len(xfrmat) > 0 else ()
	return (prods,symctrbs,invmat)


# Symmetric YD's to General YD's
def SymYDsToGenYDsFunct(SymFuncList):
	# The operations
	ZeroOp = CountedYDListZero
	ScalMultOp = CountedYDListScalMult
	ScalDivOp = CountedYDListScalDiv
	AddOp = CountedYDListAdd
	AddListOp = CountedYDListAddList
	MultOp = CountedYDListMult
	
	pwr = len(SymFuncList)
	dgrmlist = []
	for p in range(pwr):
		yspg = YDSymGenProdsToGenYDs(p+1)
		prds = yspg[0]
		scts = yspg[1]
		imts = yspg[2]
		pdvs = []
		for i,prd in enumerate(prds):
			pdv0 = SymFuncList[prd[0]-1]
			pdv1 = dgrmlist[prd[1][0]][prd[1][1]]
			pdv2 = MultOp(pdv0,pdv1)
			pdv3 = ScalMultOp(-scts[i],SymFuncList[p])
			pdv = AddOp(pdv2,pdv3)
			pdvs.append(pdv)
		dgrms = [SymFuncList[p]]
		if len(imts) > 0:
			imtsml = imts[0]
			imtsmt = imts[1]
			for i in range(len(prds)):
				imd = imtsml[i]
				imt = imtsmt[i]
				pdvlst = [ScalMultOp(imt[j],pdv) for j,pdv in enumerate(pdvs)]
				pdvres = ScalDivOp(imd,AddListOp(pdvlst))
				dgrms.append(pdvres)
		dgrmlist.append(dgrms)
	return dgrmlist

# Power Functions to General YD's *)
def GenYDsFromPowers(PwrFuncList):
	return SymYDsToGenYDsFunct(SymYDsFromPowers(1,PwrFuncList))

# General YD's of powers: in Schur functions, S(a,x^p) = sum of S(b,x)'s
def GenYDsOfPower(pwr,nmax):
	PFL = [PowerSumToYDs(pwr*p) for p in range(1,nmax+1)]
	return GenYDsFromPowers(PFL)

# Select which one desired:
#
# last symmetric -- "Sym"
# all symmetric -- "SymAll"
# last antisymmetric -- "Ats"
# all antisymmetric -- "AtsAll"
# last set of general -- "Gen"
# all sets of general -- "GenAll"
#
# All results as single-level lists,
# making or flattening result lists if necessary
#
# Necessary for plethysm code to control how much one wants to calculate.

def SelYDsFunction(sel,PwrFuncList):
	if sel == "Sym":
		res = SymYDsFromPowers(1,PwrFuncList)
		res = [res[-1]]
	elif sel == "SymAll":
		res = SymYDsFromPowers(1,PwrFuncList)
	elif sel == "Ats":
		res = SymYDsFromPowers(-1,PwrFuncList)
		res = [res[-1]]
	elif sel == "AtsAll":
		res = SymYDsFromPowers(-1,PwrFuncList)
	elif sel == "Gen":
		res = GenYDsFromPowers(PwrFuncList)
		res = res[-1]
	elif sel == "GenAll":
		res = GenYDsFromPowers(PwrFuncList)
		rsx = []
		for r in res: rsx += r
		res = rsx
	else:
		res = []
	return tuple(res)

def SelYDsOfPower(sel,pwr,nmax):
	PFL = [PowerSumToYDs(pwr*p) for p in range(1,nmax+1)]
	return SelYDsFunction(sel,PFL)
	
# Selection and max order of plathysmer,
# selection and max order of plethysmed one --
# different plethysmed one in each row
# and different plethysmer in each column.
# In each cell is counted list of YD's

def SelYDPlethysm(pthsel,pthnmax,sel,nmax):
	PFLArr = [SelYDsOfPower(sel,p,nmax) for p in range(1,pthnmax+1)]
	d1 = len(PFLArr); d2 = len(PFLArr[0])
	PFLATr = [[PFLArr[j][i] for j in range(d1)] for i in range(d2)]
	return tuple((SelYDsFunction(pthsel,PFL) for PFL in PFLATr))

# Analogue to
#
# S(a,x) = sum over p of (N(p)/N)*Xsym(a,p) *
# product over i of (power sum for power i)^b(i)
# -- multiplicity p (i) of row length i
#
# implemented in in MultiSchurPowerSum here
#
# S(a)[S(b,x)] = sum over p of (N(p)/N)*Xsym (a,p) *
# product over i of (S(b)[power sum for power i])^b(i) --
# multiplicity p (i) of row length i
#
# Power sum -> plethysm of power sum

	
