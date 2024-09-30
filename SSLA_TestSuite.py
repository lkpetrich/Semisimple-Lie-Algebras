#!python3
#
# Suite of tests of semisimple-Lie-algebra routines
#

TestList = """
vm - Vector, matrix arithmetic
la - Lie algebras
rpwo - Weyl orbits
rpir - Irreps (irreducible representations)
rppt - Irrep properties
rpcb - Rep combinations
rppd - Rep products
rppw - Rep powers
sbdm - Subalgebras: root demotion
sbxs - Subalgebras: extension splitting
sbmp - Subalgebras: matrix outer products
sbeo - Subalgebras: even to two odd SO
sbmr - Subalgebras: matrix reduction: SU to SO, Sp
sbht - Subalgebras: height ones
sbvc - Subalgebras: vector ones
sbxt - Subalgebras: extra ones
brct - Branchers: concatenation
brcj - Branchers: conjugation
brrn - Branchers: renaming
brra - Branchers: rearranging
brse - Branchers: self

all - do all of these tests
"""

from SemisimpleLieAlgebras import *
# from YoungDiagrams import *

from sys import argv, exit

from time import thread_time

class Timer:
	def __init__(self):
		self.reset()
	
	def reset(self):
		self.time = thread_time()
	
	def elapsed(self):
		return (thread_time() - self.time)


# Highest weights for the fundamental and adjoint irrep of E8
E8FA = (0,0,0,0,0,0,1,0)

# List with only that one, smallest one at size 248
E8FAList = (E8FA,)

# Identity matrix as tuples	
def idmt(n): return tuple(map(tuple,identmat(n)))


def pes(*args): print(*args, end=" ")

def ShowSubalgBranching(br, wtlist):
	pes(br.latype)
	pes("--")
	for alg in br.SubAlgebras(): pes(alg)
	pes("--")
	for u1 in br.u1s: pes(u1)
	print()
		
	for wts in wtlist:
		print("-",wts)
		sawts = br.DoBranching(wts)
		for wt in sawts: print(wt)
	print()


def DisplayMembers(xlst):
	for x in xlst: print(x)

def DisplayNumberedMembers(xlst):
	for n, x in enumerate(xlst): print(n+1, "-", x)


Tests = argv[1:]
if len(Tests) == 0:
	print(TestList)
	exit()

if "all" in Tests:
	Tests = ("vm" "la" \
		"rpwo" "rpir" "rppt" "rpcb" "rppd" "rppw" \
		"sbdm" "sbxs" "sbmp" "sbeo" "sbmr" "sbht" "sbvc" "sbxt" \
		"brct" "brcj" "brrn" "brra" "brse")


if "vm" in Tests:
	print("Vector, matrix arithmetic")
	print()
	mat = [[1,2,3],[4,5,6],[7,8,9]]
	print("matrix",mat)
	print("transpose", transpose(mat))
	vec1 = [1,2,3]
	vec2 = [4,5,6]
	vec3 = [7,8,9]
	print("vec1,2,3", vec1, vec2, vec3)
	print("add_vv", add_vv(vec1,vec2))
	print("add_vvv", add_vvv(vec1,vec2,vec3))
	print("sub_vv", sub_vv(vec1,vec2))
	vecx = vec1.copy()
	print()


if "la" in Tests:
	print("Lie algebras")
	print()
	
	def DumpLieAlgebra(la):
		print("Type:",la.latype)
		print("Name:",la.name)
		print("Special irreps:",la.special)
		print("Dynkin:",la.dynkin)
		print("Metric:",la.metric)
		print("Inverse metric:",la.invmet)
		print("Inverse metric numerators:",la.imetnum)
		print("Inverse metric denominator:",la.imetden)
		print("Cartan matrix:",la.ctnmat)
		print("Inverse Cartan matrix:",la.invctn)
		print("Inverse Cartan matrix numerators:",la.ictnnum)
		print("Inverse Cartan matrix denominator:",la.ictnden)
		print("posroots:",la.posroots)
		print("posrootsum:",la.posrootsum)
		print()
	
	DumpLieAlgebra(LieAlgebra((1,2)))
	DumpLieAlgebra(LieAlgebra((2,2)))
	DumpLieAlgebra(LieAlgebra((3,2)))
	DumpLieAlgebra(LieAlgebra((4,2)))
	DumpLieAlgebra(LieAlgebra((7,2)))
	
	def AlgDimension(latype): return GetLieAlgebra(latype).dimension()
	
	print("Check on algebra dimensions")
	print([AlgDimension((1,n)) - n*(n+2) for n in range(1,10+1)])
	print([AlgDimension((2,n)) - n*(2*n+1) for n in range(1,10+1)])
	print([AlgDimension((3,n)) - n*(2*n+1) for n in range(1,10+1)])
	print([AlgDimension((4,n)) - n*(2*n-1) for n in range(2,10+1)])
	xcdim = [AlgDimension(latype) for latype in ((7,2),(6,4),(5,6),(5,7),(5,8))]
	xcval = [14, 52, 78, 133, 248]
	print(sub_vv(xcdim,xcval))
	print()


if "rpwo" in Tests:
	print("Weyl orbits")
	print()
	
	print("Exceptional algebras' Weyl-group orders relative to subalgebras' ones")
	def ShowWGOrigSub(laorig, lasub):
		ordorig = WeylGroupOrder(laorig)
		ordsub = WeylGroupOrder(lasub)
		print(laorig, lasub, (ordorig/ordsub))
	ShowWGOrigSub( (7,2), (1,2) )
	ShowWGOrigSub( (6,4), (2,4) )
	for n in range(6,8+1):
		ShowWGOrigSub( (5,n), (1,n-1) )
		ShowWGOrigSub( (5,n), (4,n-1) )
	
	def OrbitTimingVerify(latype, maxwts):
		CheckLARep(latype, maxwts)
		
		print(latype, maxwts)
		
		tmr = Timer()
		
		tmr.reset()		
		wodw = WeylOrbitForDomWt(latype, maxwts)
		print(tmr.elapsed())
		
		tmr.reset()
		woxp = WeylOrbitForDomWtExplicit(latype, maxwts)
		print(tmr.elapsed())
		
		print(wodw == woxp)
	
	wts2 = (11,19)
	for n in (1,2,3,4,7):
		OrbitTimingVerify( (n,2), wts2)
	
	wts4 = (11,19,29,41)
	for n in (1,2,3,4,6):
		OrbitTimingVerify( (n,4), wts4)


if "rpir" in Tests:
	print("Irreps (irreducible representations)")
	print()
		
	def IrrepTimingVerify(latype, maxwts):
		CheckLARep(latype, maxwts)
		
		print(latype, maxwts)
		
		tmr = Timer()
		
		tmr.reset()		
		irdr = GetRepDirect(latype, maxwts)
		print(tmr.elapsed())
		
		WeylOrbitCache = {}
		tmr.reset()
		irxp = GetRepOrbitsExpanded(latype, maxwts)
		print(tmr.elapsed())
		
		print(irdr == irxp)

	
	wts2 = (2,3)
	for n in (1,2,3,4,7):
		IrrepTimingVerify( (n,2), wts2)
	
	wts4 = (1,0,1,2)
	for n in (1,2,3,4,6):
		IrrepTimingVerify( (n,4), wts4)
	
if "rppt" in Tests:
	print("Irrep properties")
	print()
	
	def DegenTest(latype, maxwts):
		dgn = TotalDegen(latype, maxwts)
		dgnxp = TotalDegenOfExpRep(GetRep(latype, maxwts))
		print(latype, maxwts, dgnxp == dgn, dgn)
	
	print("Degeneracy Tests: from formula ?= from irrep")
	
	wts2 = (2,3)
	for n in (1,2,3,4,7):
		DegenTest( (n,2), wts2)
	
	wts4 = (1,0,1,2)
	for n in (1,2,3,4,6):
		DegenTest( (n,4), wts4)
	
	print()
	
	print("Find conjugates and test them")
	
	def ConjgTest(latype):
		family,n = latype
		wts0 = tuple(range(1,n+1))
		wts1 = RepConjugate(latype,wts0)
		print(latype,end="")
		if wts0 == wts1:
			print(" Self-conjugate")
		else:
			print(" Pair")
			print(wts0)
			print(wts1)
	
	# Exceptional algebras
	excs = ((7,2),(6,4),(5,6),(5,7),(5,8))
	
	for f in range(1,4+1):
		for n in range(4,5+1):
			latype = (f,n)
			ConjgTest(latype)
	for latype in excs:
		ConjgTest(latype)
	
	print()
			
	print("Basic-Root Features")
	
	def BasicIrrepFeatures(latype):
		family,n = latype
		CnsvMods = RepConservModuli(latype)
		print("LA, mods:", latype, CnsvMods)
		wts0 = n*[0]
		dgns = []
		rlts = []
		csvs = []
		for k in range(n):
			wts = wts0.copy()
			wts[k] = 1
			dgn = TotalDegen(latype,wts)
			dgns.append(dgn)
			rlt = RepReality(latype,wts)
			if rlt == 0: rlt = "R"
			elif rlt == 1: rlt = "H"
			elif rlt == -1: rlt = "C"
			rlts.append(rlt)
			csv = RepConserv(latype,wts)
			csvs.append(csv)
		print("Degens:", dgns)
		print("Reality:", ' '.join(rlts))
		print("Conserv:", csvs)
		print()

	for f in range(1,4+1):
		for n in range(1,4+1):
			if f == 4 and n == 1: continue
			latype = (f,n)
			BasicIrrepFeatures(latype)
	for latype in excs:
		BasicIrrepFeatures(latype)
	
	print()

if "rpcb" in Tests:
	print("Rep combinations")
	
	def RepsFromOrbits(latype, orbits):
		rpix = {}
		for orbit in orbits:
			AddEntryToIndexedRep(rpix, orbit)
		
		return ExtractRepOrbitIrreps(latype, rpix)
	
	def MakeCountedList(rptype, x):
		rptc = rptype.casefold()
		if rptc == "sngl": return ((1,x))
		elif rptc == "list": return tuple( ( (1,xi) for xi in x) )
		elif rptc == "cntd": return x
	
	print("Rep Base")
	print(GetRep((1,2), (0,1)))
	print("Single")
	print(GetRepXtnd((1,2), "Sngl", (2,0)))
	print("List")
	print(GetRepXtnd((1,2), "List", ((2,0),(0,1))))
	print("Counted list")
	print(GetRepXtnd((1,2), "Cntd", ((1,(2,0)),(2,(0,1)))))
	print()
	
	print("Orbits in a rep")
	print(GetRepOrbits((1,2), (1,1)))
	print("Reps from orbits")
	print(RepsFromOrbits((1,2),((2,(1,1),(1,1)),(7,(0,0),(0,0)))))
	print("Make counted lists")
	print(MakeCountedList("sngl","a"))
	print(MakeCountedList("list",("a","b")))
	print(MakeCountedList("cntd",((1,"a"),(2,"b"),(3,"c"))))
	print()
	
	print("Irrep of algebra product")
	print(GetAlgProdRep( ((1,1),(1,2)), ((1,),(1,0)) ))
	print("Xtnd sngl")
	print(GetAlgProdRepXtnd( ((1,1),(1,2)), "sngl", ((1,),(1,0)) ))
	print("Xtnd list")
	print(GetAlgProdRepXtnd( ((1,1),(1,2)), "list", (((1,),(1,0)),) ))
	print("Xtnd cntd")
	print(GetAlgProdRepXtnd( ((1,1),(1,2)), "cntd", ((1,((1,),(1,0))),) ))

if "rppd" in Tests:
	print("Rep products")
	
	print(DecomposeRepProduct((1,2), (1,0), (1,0)))
	print(DecomposeRepProduct((1,2), (1,0), (0,1)))
	print(DecomposeRepProduct((1,2), (1,0), (1,1)))
	print(DecomposeRepProductXtnd((1,2), \
		"sngl", (1,0), "list", ((1,1),(0,0))))
	print(DecomposeRepProductXtnd((1,2), \
		"sngl", (1,0), "cntd", ((1,(1,1)),(2,(0,0)))))
	print()
	
	print("Rep-list products")
	print(DecomposeRepProdList((1,2), ()))
	print(DecomposeRepProdList((1,2), ((1,0),)))
	print(DecomposeRepProdList((1,2), ((1,0),(1,0),(1,0))))
	print()
	
	print("Algebra-product rep products")
	print(DecomposeAlgProdRepProduct( ((1,1),(1,2)), ((1,),(1,0)), ((1,), (1,0)) ))
	print(DecomposeAlgProdRepProductXtnd( ((1,1),(1,2)), \
		"sngl", ((1,),(1,0)), "cntd", ( (1,((2,),(1,1))), (2,((0,),(0,0))) ) ))
	print(DecomposeAlgProdRepProduct( (), (1,2,3), (4,5,6) ))
	print()
	
	print("Alg-prod rep-list products")
	print(DecomposeAlgProdRepProdList( ((1,1),(1,2)), ( ((1,),(1,0)), ) ))
	print(DecomposeAlgProdRepProdList( ((1,1),(1,2)), \
		( ((1,),(1,0)), ((1,),(1,0)), ((1,),(1,0)) ) ))
	
if "rppw" in Tests:
	print("Rep powers")
	
	print(GetTensorPowerYDX(0))
	DisplayMembers(DecomposeRepPower((1,1), (1,), 0))
	print(GetTensorPowerYDX(1))
	DisplayMembers(DecomposeRepPower((1,1), (1,), 0))
	print(GetTensorPowerYDX(2))
	DisplayMembers(DecomposeRepPower((1,1), (1,), 2))
	DisplayMembers(DecomposeRepPower((1,2), (1,0), 2))
	print(GetTensorPowerYDX(3))
	DisplayMembers(DecomposeRepPower((1,2), (1,0), 3))
	print()
	
	print("Symmetric and antisymmetric")
	print(DecomposeRepPwrSym((1,2), (1,0), 3, 1))
	print(DecomposeRepPwrSym((1,2), (1,0), 3, -1))
	print()
	
	print("Extended, with symmetric and antisymmeric")
	print(GetTensorPowerYDX(4))
	DisplayMembers(DecomposeRepPowerXtnd((1,2), "sngl", (1,0), 4))
	print(DecomposeRepPwrSymXtnd((1,2), "sngl", (1,0), 4, 1))
	print(DecomposeRepPwrSymXtnd((1,2), "sngl", (1,0), 4, -1))
	print()
	
	print("Algebra-product rep powers]")
	DisplayMembers(DecomposeAlgProdRepPower(((1,1),(1,2)), ((1,),(1,0)), 3))
	DisplayMembers(DecomposeAlgProdRepPowerXtnd(((1,1),(1,2)), "sngl", ((1,),(1,0)), 3))
	print()
	
	print("Symmetric and antisymmetric")
	print(DecomposeAlgProdRepPwrSym(((1,1),(1,2)), ((1,),(1,0)), 3, 1))
	print(DecomposeAlgProdRepPwrSymXtnd(((1,1),(1,2)), "sngl", ((1,),(1,0)), 3, 1))
	print(DecomposeAlgProdRepPwrSym(((1,1),(1,2)), ((1,),(1,0)), 3, -1))
	print(DecomposeAlgProdRepPwrSymXtnd(((1,1),(1,2)), "sngl", ((1,),(1,0)), 3, -1))
	print()
	
	print("Empty product")
	DisplayMembers(DecomposeAlgProdRepPower((), (1,2,3), 3))
	print()
	

if "sbdm" in Tests:
	print("Subalgebras: root demotion")
	for latype in ((1,5), (2,5), (3,5), (4,5), (5,6), (5,7), (5,8), (6,4), (7,2)):
		print(latype)
		DisplayNumberedMembers(ListRootDemotions(latype))
	
	print()
	
	def ShowRootDemotions(latype, dmroot, wtlist):
		pes(dmroot, "--")
		ShowSubalgBranching(MakeRootDemoter(latype, dmroot), wtlist)
	
	ShowRootDemotions((1,1), 1, ((k,) for k in range(4+1)))
	
	ShowRootDemotions((1,4), 3, \
		((1,0,0,0), (0,0,0,1), (0,0,1,0), (0,1,0,0), (1,0,0,1)))
	
	ShowRootDemotions((4,5), 5, \
		((0,0,0,1,0), (0,0,0,0,1), (1,0,0,0,0), (0,1,0,0,0)))
	
	ShowRootDemotions((5,6), 5, \
		((1,0,0,0,0,0), (0,0,0,0,1,0), (0,0,0,0,0,1)))
	
	
if "sbxs" in Tests:
	print("Subalgebras: extension splitting")
	for latype in ((1,6), (2,6), (3,6), (4,6), (5,6), (5,7), (5,8), (6,4), (7,2)):
		print(latype)
		DisplayNumberedMembers(ListExtensionSplits(latype))
	
	print()
	
	def ShowExtensionSplits(latype, xsroot, wtlist):
		pes(xsroot, "--")
		ShowSubalgBranching(MakeExtensionSplitter(latype, xsroot), wtlist)
	
	ShowExtensionSplits((2,7), 4, \
		((1,0,0,0,0,0,0), (0,0,0,0,0,0,1), (0,1,0,0,0,0,0)))
	
	ShowExtensionSplits((4,7), 4, \
		((1,0,0,0,0,0,0), (0,0,0,0,0,1,0), (0,0,0,0,0,0,1), (0,1,0,0,0,0,0)))
	
	ShowExtensionSplits((5,8), 6, E8FAList)
	
	ShowExtensionSplits((5,8), 5, E8FAList)
	
	ShowExtensionSplits((5,8), 4, E8FAList)
	
	
if "sbmp" in Tests:
	print("Subalgebras: matrix outer products")
	print()
	
	# Using a single weight vector instead of a list of them
	
	def ShowSubalgMultAn(ords, wts):
		pes(ords, "--")
		ShowSubalgBranching(SubalgMultAn(ords), (wts,))

	def ShowSubalgMultBCDn(stypes, wts):
		pes(stypes, "--")
		ShowSubalgBranching(SubalgMultBCDn(stypes), (wts,))
	
	ShowSubalgMultAn((3,1), (1,0,0,0,0,0,0))
	ShowSubalgMultAn((2,2), (1,0,0,0,0,0,0,0))
	
	ShowSubalgMultBCDn(((4,4),(4,1)),(1,0,0,0,0,0,0,0))
	ShowSubalgMultBCDn(((4,4),(3,1)),(1,0,0,0,0,0,0,0))
	ShowSubalgMultBCDn(((3,4),(4,1)),(1,0,0,0,0,0,0,0))
	ShowSubalgMultBCDn(((3,4),(3,1)),(1,0,0,0,0,0,0,0))
	
	ShowSubalgMultBCDn(((4,3),(2,1)),(1,0,0,0,0,0,0,0,0))
	ShowSubalgMultBCDn(((3,3),(2,1)),(1,0,0,0,0,0,0,0,0))
	
	ShowSubalgMultBCDn(((2,2),(2,1)),(1,0,0,0,0,0,0))


if "sbeo" in Tests:
	print("Subalgebras: even to two odd SO")
	print()
	
	def ShowSubalgEvenOdd(n, k, wtlist):
		pes(n, k, "--")
		ShowSubalgBranching(SubalgSOEvenOdd(n, k), wtlist)
	
	for k in range(0,3+1):
		ShowSubalgEvenOdd(5,k, \
			((1,0,0,0,0), (0,0,0,1,0), (0,0,0,0,1), (0,1,0,0,0)))


if "sbmr" in Tests:
	print("Subalgebras: matrix reduction: SU to SO, Sp")
	print()
	
	def ShowSubalgSUSO(n, wtlist):
		pes(n, "--")
		ShowSubalgBranching(SubalgSUSO(n), wtlist)
	
	def ShowSubalgSUSp(n, wtlist):
		pes(n, "--")
		ShowSubalgBranching(SubalgSUSp(n), wtlist)

	ShowSubalgSUSO(2, ((k,) for k in range(0,4+1)))
	
	ShowSubalgSUSO(8, \
		((1,0,0,0,0,0,0), (0,0,0,0,0,0,1), (1,0,0,0,0,0,1)))
	
	ShowSubalgSUSO(9, \
		((1,0,0,0,0,0,0,0), (0,0,0,0,0,0,0,1), (1,0,0,0,0,0,0,1)))
	
	ShowSubalgSUSp(4, \
		((1,0,0,0,0,0,0), (0,0,0,0,0,0,1), (1,0,0,0,0,0,1)))


if "sbht" in Tests:
	print("Subalgebras: height ones")
	print()
	
	def ShowSubalgHeight(latype, wtlist):
		ShowSubalgBranching(SubalgHeightA1(latype), wtlist)
	
	ShowSubalgHeight((1,5), \
		((1,0,0,0,0),(0,0,0,0,1),(1,0,0,0,1)))
	
	ShowSubalgHeight((2,5), \
		((1,0,0,0,0),(0,0,0,0,1),(0,1,0,0,0)))
	
	ShowSubalgHeight((3,5), \
		((1,0,0,0,0),(2,0,0,0,0)))
	
	ShowSubalgHeight((4,5), \
		((1,0,0,0,0),(0,0,0,1,0),(0,0,0,0,1),(0,1,0,0,0)))
	
	ShowSubalgHeight((7,2), \
		((0,1),(1,0)))
	
	ShowSubalgHeight((6,4), \
		((0,0,0,1),(1,0,0,0)))
	
	ShowSubalgHeight((5,6), \
		((1,0,0,0,0,0),(0,0,0,0,1,0),(0,0,0,0,0,1)))
	
	ShowSubalgHeight((5,7), \
		((0,0,0,0,0,1,0),(1,0,0,0,0,0,0)))
	
	ShowSubalgHeight((5,8), \
		((0,0,0,0,0,0,1,0),))

	
if "sbvc" in Tests:
	print("Subalgebras: vector ones")
	print()
		
	def ShowSubalgVector(family, dsttype, dstwts, wtlist):
		pes(dstwts,"--")
		ShowSubalgBranching(SubalgVector(family, dsttype, dstwts), wtlist)
	
	ShowSubalgVector(1,(1,2),(2,0), \
		((1,0,0,0,0), (0,0,1,0,0), (1,0,0,0,1)))
		
	ShowSubalgVector(4,(2,2),(0,2), \
		((1,0,0,0,0), (0,0,0,0,1), (0,1,0,0,0)))	
		
	ShowSubalgVector(2,(1,3),(1,0,1), \
		((1,0,0,0,0,0,0), (0,1,0,0,0,0,0)))
		
	ShowSubalgVector(4,(2,2),(2,0), \
		((1,0,0,0,0,0,0), (0,1,0,0,0,0,0)))
		
	ShowSubalgVector(4,(3,3),(0,1,0), \
		((1,0,0,0,0,0,0), (0,1,0,0,0,0,0)))
		
	ShowSubalgVector(4,(7,2),(1,0), \
		((1,0,0,0,0,0,0), (0,1,0,0,0,0,0)))
		
	ShowSubalgVector(3,(2,2),(1,1), \
		((1,0,0,0,0,0,0,0), (0,1,0,0,0,0,0,0)))
		
	ShowSubalgVector(4,(2,4),(0,0,0,1), \
		((1,0,0,0,0,0,0,0), (0,1,0,0,0,0,0,0)))


if "sbxt" in Tests:
	print("Subalgebras: extra ones")
	print()

	def ShowSubalgExtra(saname, wtlist):
		pes(saname,"--")
		ShowSubalgBranching(SubalgExtra(saname), wtlist)
	
	def ShowSubalgVector(family, dsttype, dstwts, wtlist):
		pes(dstwts,"--")
		ShowSubalgBranching(SubalgVector(family, dsttype, dstwts), wtlist)
	
	def ShowExtensionSplits(latype, xsroot, wtlist):
		pes(xsroot, "--")
		ShowSubalgBranching(MakeExtensionSplitter(latype, xsroot), wtlist)
	
	ShowSubalgExtra("E8G2F4",E8FAList)
	
	e6ws = ((1,0,0,0,0,0),(0,0,0,0,1,0),(0,0,0,0,0,1))
	
	ShowSubalgExtra("E6F4",e6ws)
	
	ShowSubalgExtra("D4G2",idmt(4))
	
	ShowSubalgExtra("B3G2",idmt(3))
	
	ShowSubalgVector(2, (7,2), (0,1), idmt(3))
	
	ShowSubalgExtra("D4A2",idmt(4))
	
	ShowSubalgVector(4, (1,2), (1,1), idmt(4))
	
	ShowExtensionSplits((5,6), 3, e6ws)
	
	ShowExtensionSplits((6,4), 1, ((1,0,0,0),(0,0,0,1)))
	
	ShowExtensionSplits((7,2), 2, idmt(2))


if "brct" in Tests:
	print("Branchers: concatenation")
	print()
	
	ba64 = MakeRootDemoter((1,6), 4)
	ba32 = MakeRootDemoter((1,3), 2)
	ba22 = MakeRootDemoter((1,2), 2)
	
	wtlist = ((1,0,0,0,0,0),(0,0,0,0,0,1))
	
	ShowSubalgBranching(ba64, wtlist)
	
	ba6432 = ConcatBranchers(ba64,1,ba32)
	
	ShowSubalgBranching(ba6432, wtlist)
	
	ba6422 = ConcatBranchers(ba64,2,ba22)
	
	ShowSubalgBranching(ba6422, wtlist)
	
	ba643222 = ConcatBranchers(ba6432,3,ba22)
	
	ShowSubalgBranching(ba643222, wtlist)
	

if "brcj" in Tests:
	print("Branchers: conjugation")
	print()
	
	ba153 = MakeRootDemoter((1,5),3)
	wtls15 = ((1,0,0,0,0), (0,0,0,0,1))
	
	ShowSubalgBranching(ba153, wtls15)
	
	ba1531 = BrancherConjugate(ba153, (1,))
	
	ShowSubalgBranching(ba1531, wtls15)
	
	ba1532 = BrancherConjugate(ba1531, (1,2))
	
	ShowSubalgBranching(ba1532, wtls15)
	
	ba244 = MakeExtensionSplitter((2,4), 4)
	
	ShowSubalgBranching(ba244, idmt(4))
	
	ba2441 = BrancherConjugateD4(ba244, 1, (2,3,1))
	
	ShowSubalgBranching(ba2441, idmt(4))
	

if "brrn" in Tests:
	print("Branchers: renaming")
	print()
	
	ba721 = MakeExtensionSplitter((7,2), 1)
	
	print("Rename A1, B1, C1: SU(2), SO(3), Sp(2)")
	print(ba721.SubAlgebras())
	
	ba7211 = BrancherRenameABC1(ba721, 1, 2)
	
	print("1: 1 to 2")
	print(ba7211.SubAlgebras())
	
	ba7212 = BrancherRenameABC1(ba7211, 2, 3)
	
	print("2: 1 to 3")
	print(ba7212.SubAlgebras())
	
	print()
	
	print("Rename B2, C2: SO(5), Sp(4)")
	baeo52 = SubalgSOEvenOdd(5,2)
	baeowts = ((1,0,0,0,0),(0,0,0,0,1))
	ShowSubalgBranching(baeo52, baeowts)
	print()
	
	baeo521 = BrancherRenameBC2(baeo52,1)
	ShowSubalgBranching(baeo521, baeowts)
	print()
	
	baeo522 = BrancherRenameBC2(baeo521,2)
	ShowSubalgBranching(baeo522, baeowts)
	print()
	
	baeo523 = BrancherRenameBC2(baeo522,1)
	ShowSubalgBranching(baeo523, baeowts)
	print()
	
	print("Rename A3, D3: SU(4), SO(6)")
	ba463 = MakeExtensionSplitter((4,6),3)
	ba463wts = ((1,0,0,0,0,0),(0,0,0,0,1,0),(0,0,0,0,0,1))
	ShowSubalgBranching(ba463, ba463wts)
	print()
	
	ba4631 = BrancherRenameAD3(ba463,1)
	ShowSubalgBranching(ba4631, ba463wts)
	print()
	
	ba4632 = BrancherRenameAD3(ba4631,2)
	ShowSubalgBranching(ba4632, ba463wts)
	print()
	
	ba4633 = BrancherRenameAD3(ba4632,1)
	ShowSubalgBranching(ba4633, ba463wts)
	print()
	
	print("Split D2: A1,A1, join A1, A1: D2")
	ba442 = MakeExtensionSplitter((4,4),2)
	ba442wts = ((1,0,0,0),(0,0,1,0),(0,0,0,1))
	ShowSubalgBranching(ba442, ba442wts)
	print()
	
	ba4421 = BrancherSplitD2(ba442, 1)
	ShowSubalgBranching(ba4421, ba442wts)
	print()	
	
	ba4422 = BrancherSplitD2(ba4421, 3)
	ShowSubalgBranching(ba4422, ba442wts)
	print()	
	
	ba4423 = BrancherJoin2A1(ba4422, 2,3)
	ShowSubalgBranching(ba4423, ba442wts)
	print()	
	
	ba4424 = BrancherJoin2A1(ba4423, 1,2)
	ShowSubalgBranching(ba4424, ba442wts)
	print()	
	

if "brra" in Tests:
	print("Branchers: rearranging")
	print()
	
	ba511 = SubalgSOEvenOdd(5,1)
	ba51wts = ((1,0,0,0,0),(0,0,0,1,0),(0,0,0,0,1),(0,1,0,0,0))
	ShowSubalgBranching(ba511, ba51wts)
	print()
	
	ba5111 = BrancherRearrange(ba511, (2,1))
	ShowSubalgBranching(ba5111, ba51wts)
	print()
	
	ba514 = SubalgSOEvenOdd(5,3)
	ShowSubalgBranching(ba5111, ba51wts)
	print()

if "brse" in Tests:
	print("Branchers: self")
	print()
	
	base = SubalgSelf((1,3))
	print(base.latype)
	for stsm in base.stsms:
		print(stsm)
	print(base.u1s)
