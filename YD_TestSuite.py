#!python3
#
# Test suite for the Young-diagram code
#

from YoungDiagrams import *
from sys import exit, argv

if len(argv) <= 1:
	print("Need at least one type:")
	print("op: Young-Diagram Operations")
	print("dg: Total Degeneracies: Multiplicities")
	print("pd: Products of Young diagrams for SU(n) reps")
	print("vp: Vector-power verification")
	print("km: Kostka-matrix verification: Weyl orbits in SU(n) reps")
	print("sc: Symmetric-group characters")
	print("pl: Plethysms: symmetrized powers")
	exit()

args = argv[1:]

def ShowYRLS(yrls):
	print(yrls)
	for row in yrls:
		for k in range(row):
			print("# ",end='')
		print()
	print()

def ShowYDGM(ydgm):
	for row in ydgm: print(row)
	print()

def CSV(lst):
	return ','.join(map(str,lst))

def CSVX(*lst):
	return ','.join(map(str,lst))

# Verify rep products

def SUnProductCtYoungVerify(yrls1,yrls2,n):
	prod = SUnProductCtYoung(yrls1,yrls2)
	origdims = tuple(((yrls, DegenSUnCtYoung(yrls,n)) for yrls in (yrls1,yrls2)))
	proddims = tuple(((yrls, DegenSUnCtYoung(yrls,n)) for yrls in prod))
	origlen = 1
	for mem in origdims: origlen *= mem[1]
	prodlen = 0
	for mem in proddims: prodlen += mem[1]
	return ((origdims,proddims),(origlen,prodlen))

def SUnProductFxdVerify(fwts1,fwts2):
	prod = SUnProductFxd(fwts1,fwts2)
	if prod == None: return
	origdims = tuple(((fwts, DegenSUnFxd(fwts)) for fwts in (fwts1,fwts2)))
	proddims = tuple(((fwts, DegenSUnFxd(fwts)) for fwts in prod))
	origlen = 1
	for mem in origdims: origlen *= mem[1]
	prodlen = 0
	for mem in proddims: prodlen += mem[1]
	return ((origdims,proddims),(origlen,prodlen))

# SU(n) Vector-Representation Powers
# Verify that they are in the proper order -- add to earlier row length,
# subtract that amount from later row length
# should give an earlier diagram, and only an earlier one.

def VPDgrmListOrderVerify(p):
	dglst = VPFindDgrmList(p)
	# Set up indexing
	dgidx = {}
	for i,dg in enumerate(dglst):
		dgidx[dg] = i
	# Shoot backwards
	for i,dg in enumerate(dglst):
		for j1 in range(len(dg)-1):
			for j2 in range(j1+1,len(dg)):
				for k in range(1,dg[j2]+1):
					dgb = list(dg)
					dgb[j1] += k
					dgb[j2] -= k
					dgb.sort()
					dgb.reverse()
					if dgb[-1] == 0: del dgb[-1]
					dgb = tuple(dgb)
					ix = dgidx[dgb]
					if ix >= i: return False
	return True

# Should make a list of True's
def VPOrderVerify(p):
	return [VPDgrmListOrderVerify(ip) for ip in range(1,p+1)]

# Verify total sizes for each level
def VPDgrmListSizeVerify(p,n):
	dglst = VPFindMultDgrmList(p)
	total = 0
	for num,dg in dglst:
		total += num*DegenSUnCtYoung(dg,n)
	total -= n**p
	return int(total)

# Should make a list of 0' s
def VPSizeVerify(p,n):
	return [VPDgrmListSizeVerify(ip,n) for ip in range(1,p+1)]

# For Kostka matrices: decomposition of SU(n) irreps into Weyl orbits
def KostkaMatrixVerify(p,n):
	diffs = []
	for ip in range(1,p+1):
		vpds = VPFindDgrmList(ip)
		kostka = KostkaMatrix(ip)
		tgtlens = [DegenSUnCtYoung(vpd,n) for vpd in vpds]
		indlens = [SUnWeylOrbitDimCtYoung(vpd,n) for vpd in vpds]
		vpdlen = len(vpds)
		diff = vpdlen*[0]
		for j in range(vpdlen):
			kkrow = kostka[j]
			iksum = 0
			for k in range(vpdlen):
				iksum += kkrow[k]*indlens[k]
			diff[j] = tgtlens[j] - iksum
		diffs.append(tuple(map(int,diff)))
	return tuple(diffs)

if "op" in args:

	print("Young-Diagram Operations")
	print()
	
	yrls = (4,2,1)
	ydgm = ExpandYoung(yrls, CSVX)
	ywtvr = CtYoungToVarWts(yrls)
	ywtfx = CtYoungToFxdWts(yrls,5)
	
	ShowYRLS(yrls)
	ShowYRLS(ContractYoung(ydgm))
	ShowYRLS(VarWtsToCtYoung(ywtvr))
	ShowYRLS(FxdWtsToCtYoung(ywtfx))
	
	ShowYDGM(ydgm)
	ShowYDGM(FxdWtsToXpYoung(ywtfx, CSVX))
	ShowYDGM(ydgm)
	
	ShowYRLS(ywtvr)
	ShowYRLS(XpYoungToVarWts(ydgm))
	ShowYRLS(VarWeights(ywtfx))
	
	ShowYRLS(ywtfx)
	ShowYRLS(XpYoungToFxdWts(ydgm,5))
	ShowYRLS(FixWeights(ywtvr,5))

	ShowYDGM(RelabelYoung(((1,2,3,4),(5,6),(7,)), CSVX))
	
	ShowYDGM(TransposeXpYoung(ydgm))
	ShowYRLS(TransposeCtYoung(yrls))
	ShowYRLS(FastTpCtYoung(yrls))
	print()

if "dg" in args:

	print("Degeneracies: Representation Sizes")
	print()
	
	print("Weight vector, then SU(n), SO(2n), SO(2n+1), Sp(2n)")
	zero = 4*[0]
	for k in range(len(zero)):
		wts = zero.copy()
		wts[k] = 2
		print(wts)
		print(DegenSUnFxd(wts),end=" ")
		print(DegenSO2nFxd(wts),end=" ")
		print(DegenSOnp1Fxd(wts),end=" ")
		print(DegenSpnFxd(wts),end="\n\n")
	
	print("General-size formulas")
	ydref = (4,2,1)
	wtref = CtYoungToFxdWts(ydref,6)
	print(ydref, " ", wtref)
	print("SU(n)")
	print(DegenSUnFxd(CtYoungToFxdWts(ydref,6)),end=" ")
	print(DegenSUnCtYoung(ydref,6),end=" ")
	print(DegenEval(DegenFacsSUnCtYoung(ydref),6))
	print("SO(2n)")
	print(DegenSO2nFxd(wtref),end=" ")
	print(DegenEval(DegenFacsSOnCtYoung(ydref),10))
	print("SO(2n+1)")
	print(DegenSOnp1Fxd(wtref),end=" ")
	print(DegenEval(DegenFacsSOnCtYoung(ydref),11))
	print("Sp(2n)")
	print(DegenSpnFxd(wtref),end=" ")
	print(DegenEval(DegenFacsSpnCtYoung(ydref),10))
	print()
	
if "pd" in args:
	print("Rep products: sizes are in last pair and should be equal")
	print(SUnProductCtYoungVerify( (2,), (1,1), 5))
	print(SUnProductFxdVerify( (1,0), (0,1) ))
	print()

if "vp" in args:
	print("Vector-power verify: should make lists of True's and 0's")
	print(VPOrderVerify(4))
	print(VPSizeVerify(4,6))
	print()

if "km" in args:
	print("Kostka-matrix verify: should make a list of lists of 0's")
	print(KostkaMatrixVerify(4,6))
	print()

if "sc" in args:
	print("Symmetric-group characters with class and irrep YD's")
	pmax = 4
	ydss = VPFindDgrmListSeq(pmax)
	sgcs = SymmGroupCharacters(pmax)
	for p in range(1,pmax+1):
		print(p)
		print(ydss[p-1])
		print(sgcs[p-1])
	print()
	
if "pl" in args:
	print("Plethysms: symmetrized powers")
	ps = SelYDPlethysm("GenAll",2,"GenAll",2)
	for r, psr in enumerate(ps):
		for c, psc in enumerate(psr):
			print(r, c, "-", psc)