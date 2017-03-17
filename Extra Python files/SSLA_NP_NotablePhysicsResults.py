#!/usr/bin/env python
#
# Does a dump of some interesting results

from SemisimpleLieAlgebras_NumPy import *

def DumpRep(latype, maxwts):
	rep = GetTrueRep(latype, GetRep(latype,maxwts))
	print "Degen, Roots, Weights"
	for rtx in rep:
		outlst = []
		outlst.append("%3d" % rtx[0])
		
		outlst.append("   ")
		
		for r in rtx[1]:
			outlst.append("%5s" % r)
		
		outlst.append("   ")
		
		for w in rtx[2]:
			outlst.append("%3d" % w)
		
		print ''.join(outlst)

reality = {0: "r", 1: "p", -1: "c"}

def DumpLblRepList(latype,lmwlist):
	for label, mwts in lmwlist:
		print label
		DumpRep(latype,mwts)

def DumpRepListProperties(latype, lblmwtlist):
	print "Wts, Cjg, Dim, Ht, Real, Cons"
	for label, maxwts in lblmwtlist:
		props = RepProperties(latype, maxwts)
		outlst = []
		for w in props["maxwts"]:
			outlst.append("%3d" % w)
		outlst.append("   ")
		for w in props["mwconjg"]:
			outlst.append("%3d" % w)
		outlst.append("%7d" % props["dimension"])
		outlst.append("%4d" % props["height"])
		outlst.append(" %1s" % reality[props["reality"]])
		outlst.append("   ")
		for cons in props["conserved"]:
			outlst.append("(%d,%d) " % cons)
		outlst.append("   %s" % label)
		print ''.join(outlst)

def DumpRepProdGeneric(maxwtlist,prod):
	print "Input maxwt vectors:"
	outlst = []
	for maxwts in maxwtlist:
		outlst.append("   ")
		for w in maxwts:
			outlst.append("%3d" % w)
	print ''.join(outlst)
	
	print "Mult, Wts"
	for p in prod:
		outlst = []
		outlst.append("%3d" % p[0])
		outlst.append("   ")
		for w in p[1]:
			outlst.append("%3d" % w)
		print ''.join(outlst)

def DumpRepProdTypesGeneric(maxwts,cmbntypes,prod):
	print "Input maxwt vectors:"
	outlst = []
	outlst.append("   ")
	for w in maxwts:
		outlst.append("%3d" % w)
	print ''.join(outlst)
	
	print "Mult, Wts"
	cbprd = zip(cmbntypes,prod)
	for tp,ps in cbprd:
		print tp
		for p in ps:
			outlst = []
			outlst.append("%3d" % p[0])
			outlst.append("   ")
			for w in p[1]:
				outlst.append("%3d" % w)
			print ''.join(outlst)

def DumpRepProduct(latype, maxwt1, maxwt2):
	DumpRepProdGeneric((maxwt1,maxwt2), \
		DecomposeRepProduct(latype, maxwt1, maxwt2))

def DumpRepProdList(latype, maxwtlst):
	DumpRepProdGeneric(maxwtlst, \
		DecomposeRepProdList(latype, maxwtlst))

def DumpRepSquare(latype, maxwts):
	DumpRepProdTypesGeneric(maxwts, ("Symm","Anti"), \
		DecomposeRepPower(latype, maxwts, 2))

def DumpRepCube(latype, maxwts):
	DumpRepProdTypesGeneric(maxwts, ("Symm","Mixed (2*)","Anti"), \
		DecomposeRepPower(latype, maxwts, 3))

def DumpBranching(brancher, maxwts):
	print "Original max-wts vector:"
	outlst = []
	for w in maxwts:
		outlst.append("%3d" % w)
	print ''.join(outlst)
	print "Mult, Subalg wts, U(1)'s"
	br = brancher.DoBranching(maxwts)
	lalen = len(brancher.stsms)
	for b in br:
		outlst = []
		outlst.append("%3d" % b[0])
		outlst.append("  .   ")
		for wt in b[1][:lalen]:
			for w in wt:
				outlst.append("%3d" % w)
			outlst.append("   ")
		outlst.append(". ")
		for uf in b[1][lalen:]:
			outlst.append("%5s" % uf)
		print ''.join(outlst)


def PrintHeader(text):
	print "****"
	print "****"
	print "**** " + text
	print "****"
	print "****"
	print

def PrintSubheader(text):
	print "**"
	print "** " + text
	print "**"
	print


PrintHeader("Lie Algebras: Notable physics results")

print """
This Python file applies my Python file "SemisimpleLieAlgebras.py"
to various symmetry groups important in physics.

The ones covered here are space-time geometry, particle symmetries,
and some extras.

The space-time ones are 3D, 4D, 8D, 9D, 10D, and 11D.
Though all results will be calculated for Euclidean space
(metric one-sign definite), they can carry over into Minkowskian space
(one sign flipped) and other pseudo-Euclidean spaces by analytic continuation.

The 3D one is associated with quantum-mechanical angular momentum,
and the 4D one with familiar space-time.

The 8D and 9D ones are 10D and 11D ones without the light cone,
the time dimension and one space dimension, meaning that they are for
the transverse degrees of freedom for massless particles.
Their counterpart in 4D is 2D, with symmetry group SO(2) / U(1): helicity.

The 10D one is for string theory.
There are five types of supersymmetric strings (superstrings),
and they are quantum-mechanically consistent only in 10D spacetime,
and some of them offer plausible Grand Unified Theories.

The 11D one is for M-theory: a poorly-understood theory that unites
all the five 10D superstrings in one theory that lives in 11D spacetime.


Turning to  particle symmetries, I consider gauge and flavor ones.

The simplest particle symmetry is U(1), associated with adding a complex phase
with a particle without mixing it with other kinds of particles.
Electromagnetism is a U(1) gauge symmetry, and flavor conservation laws,
like of baryon and lepton number, have U(1) symmetries associated with them.

The next sort is SU(2), which acts like 3D angular momentum,
thus "isospin" from "isotopic spin".
The unbroken Standard Model has a SU(2) symmetry called "weak isospin",
for up-like quarks vs. down-like ones and neutrinos vs. electron-like leptons.
Plain isospin is an approximate flavor symmetry between up and down quarks,
which are much less massive than hadrons.
It is respected by the strong / QCD interaction,
but not by the electromagnetic or weak interactions.

The next one is SU(3), the QCD (quantum chromodynamics) gauge symmetry,
and an approximate flavor symmetry for the three lightest quarks:
up, down, and strange.
The flavor one is somewhat broken by the strange quark's mass, so I illustrate
its breaking into SU(2) * U(1) -- isospin and hypercharge
(average electric charge of an isospin multiplet).

This is followed by some GUT ones: Georgi-Glashow SU(5),
Pati-Salam SU(4) * SU(2) * SU(2) or SO(6) * SO(4),
and SO(10), SU(3)^3, E(6), and E(8) models. The latter one
is part of the E(8) * E(8) gauge symmetry in the HE heterotic superstring,
one of the five kinds.


The mathemtical curiosities include G(2) being a subgroup of SO(7) and SO(8),
mentioned by John Baez in "The Octonions".

Representation properties:
Wts: maximum-weight vector
Cjg: conjugate of it
Dim: dimension of rep
Ht: its height
Real: its reality (r: real, p: pseudoreal, c: complex),
Cons: conserved quantities: list of (modulus/divisor, value)
"""
print

PrintHeader("Space-Time Geometry")

PrintSubheader("3-Space Angular Momentum: SO(3) ~ SU(2)")

spc3type = (1,1)
spc3replist = tuple([(twoj,) for twoj in xrange(0,2*2+1,1)])
spc3lbllist = ("Scalar; ang mom = 0", "Spinor; ang mom = 1/2", \
	"Vector; Adjoint; ang mom = 1", "Vector-Spinor; ang mom = 3/2", \
	"Symmetric Traceless 2-tensor; ang mom = 2")
spc3list = zip(spc3lbllist,spc3replist)

DumpRepListProperties(spc3type,spc3list)
print
print "An antisymmetric 2-tensor becomes a vector"
print "by multiplying by the 3D antisymmetric symbol."
print

DumpLblRepList(spc3type,spc3list)
print

print "Spinor-Spinor: 1/2 1/2"
print "Gives 1 0"
DumpRepProduct(spc3type,spc3replist[1],spc3replist[1])
print "Decomposed by symmetry"
DumpRepSquare(spc3type,spc3replist[1])
print

print "Vector-Spinor: 1 1/2"
print "Gives 3/2 1/2"
DumpRepProduct(spc3type,spc3replist[2],spc3replist[1])
print

print "Vector-Vector: 1 1"
print "Gives 2 1 0"
DumpRepProduct(spc3type,spc3replist[2],spc3replist[2])
print "Decomposed by symmetry"
DumpRepSquare(spc3type,spc3replist[2])
print

print "Three spinors (quarks in a baryon): 1/2 1/2 1/2"
print "Gives one 3/2 and two 1/2's"
DumpRepProdList(spc3type,3*[spc3replist[1]])
print "Decomposed by symmetry"
DumpRepCube(spc3type,spc3replist[1])
print

PrintSubheader("4-Space-Time Lorentz Group: SO(4) ~ SO(3) * SO(3)")

spc4type = (4,2)
spc4replist = ((0, 0), (1, 0), (0, 1), (1, 1), (2, 1), (1, 2), \
	(2, 0), (0, 2), (2, 2), (4, 0), (0, 4))
spc4lbllist = ("Scalar (0,0)", "Spinor 1 (1/2,0)", "Spinor 2 (0,1/2)", \
	"Vector (1/2,1/2)", "Vector-Spinor 1 (1,1/2)", "Vector-Spinor 2 (1/2,1)",
	"AS 2-Tensor 1 (1,0)", "AS 2-Tensor 2 (0,1)", "Symm Tls 2-Tensor (1,1)",
	"Part AS 4-Tensor 1 (2,0)", "Part AS 4-Tensor 2 (0,2)")
spc4list = zip(spc4lbllist,spc4replist)

DumpRepListProperties(spc4type,spc4list)
print

print "A Dirac spinor is a set of Spinors 1,2, and likewise for a Dirac vector-spinor"
print "Antisymmetric 2-tensors 1,2 are opposite duality eigenstates,"
print "as are their AS 4-tensor counterparts."
print "Duality is multiplying a pair of indices by the 4D antisymmetric symbol."
print "An AS 2-tensor splits into 1 and 2; this rep is the adjoint rep."
print "The Faraday tensor is an AS 2-tensor; it splits into E+i*B and E-i*B."
print "The Weyl curvature tensor is a Part AS 4-tensor that also splits."

DumpLblRepList(spc4type,spc4list)
print

print "Vector * Vector = symmetric traceless 2-tensor +"
print "antisymmetric 2-tensor + scalar"
DumpRepSquare(spc4type,spc4replist[3])
print

print "Vector * Spinor = vector-spinor + spinor"
DumpRepProduct(spc4type,spc4replist[3],spc4replist[2])
print

print "Get the Riemann tensor : symmetric product of antisymmetric 2-tensors." 
print "It contains the Weyl tensor (Part AS 4-tensor) and the Ricci tensor"
print "(traceless symmetric 2-tensor + scalar), the scalar being the Ricci scalar."
print

print "AS Tensor * AS Tensor, same duality : Weyl-like tensor, AS 2-tensor, scalar"
DumpRepSquare(spc4type,spc4replist[6])
print

print "AS Tensor * AS Tensor, opposite duality : Symmetric traceless 2-tensor"
DumpRepProduct(spc4type,spc4replist[6],spc4replist[7])
print

PrintSubheader("10-D Space-Time without the Light Cone: SO(8)")

spc8type = (4,4)
spc8replist = ((0, 0, 0, 0), (1, 0, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), \
	(0, 1, 0, 0), (2, 0, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1))
spc8lbllist = ("Scalar", "Vector", "Spinor 1", "Spinor 2", "Adjoint; AS 2-tensor", \
	"Symm tls 2-tensor", "Vector-spinor 1", "Vector-spinor 2")
spc8list = zip(spc8lbllist,spc8replist)

DumpRepListProperties(spc8type,spc8list)
print
print "Notice that the spinors are interchangeable with the vector."
print "SO(8) / D4 is the only Lie algebra that has that property."
print
print "Supergravity multiplet:"
print "- Bosonic: 64"
print "- - Symmetric traceless 2-tensor: 35"
print "- - Antisymmetric 2-tensor: 28"
print "- - Scalar: 1"
print "- Fermionic: 64"
print "- - Vector-spinor: 56"
print "- - Spinor: 8"
print

print "Vector * Vector = symmetric traceless 2-tensor +"
print "antisymmetric 2-tensor + scalar"
DumpRepSquare(spc8type,spc8replist[1])
print

print "Vector * Spinor = vector-spinor + spinor"
DumpRepProduct(spc8type,spc8replist[1],spc8replist[2])
print

print "Compactification : breakdown of 10D into 6D * 4D."
print "The latter has the light cone, making it effectively 2D."
print "Results are multiplets in the 6 dimensions - SO(6)"
print "and the helicity in 4D - U(1)"
spc8d62 = MakeRootDemoter(spc8type,1)
print spc8d62
print

print "Vector -> 6D vector + 6D scalar with helicities +-1"
DumpBranching(spc8d62,spc8replist[1])
print

print "Spinor -> two 6D spinors, each with helicity +1/2 or -1/2"
DumpBranching(spc8d62,spc8replist[2])
print

print "Antisymmetric 2-tensor -> 6D AS 2-tensor +"
print "6D vectors with helicities +-1 + a scalar with zero helicity"
DumpBranching(spc8d62,spc8replist[4])
print

print "Symmetric traceless 2-tensor -> 6D symm tcl 2-tensor +"
print "6D vectors with helicities +-1 + scalars with helicities 0, +-2"
DumpBranching(spc8d62,spc8replist[5])

print "Vector-spinor -> 6D vector-spinors with helicities +-1/2 +"
print "6D spinors with helicities +-3/2, +-1/2"
DumpBranching(spc8d62,spc8replist[6])
print

PrintSubheader("11-D Space-Time without the Light Cone: SO(9)")

spc9type = (2,4)
spc9replist = ((0, 0, 0, 0), (1, 0, 0, 0), (0, 0, 0, 1), (0, 1, 0, 0), \
	(2, 0, 0, 0), (0, 0, 1, 0), (1, 0, 0, 1))
spc9lbllist = ("Scalar", "Vector", "Spinor", "Adjoint; AS 2-tensor", \
	"Symm tls 2-tensor", "AS 3-tensor", "Vector-spinor")
spc9list = zip(spc9lbllist,spc9replist)

DumpRepListProperties(spc9type,spc9list)
print

print "Supergravity multiplet:"
print "- Bosonic: 128"
print "- - Symmetric traceless 2-tensor: 44"
print "- - Antisymmetric 3-tensor: 84"
print "- Fermionic: 128"
print "- - Vector-Spinor: 128"

print "Vector * Vector = symmetric traceless 2-tensor +"
print "antisymmetric 2-tensor + scalar"
DumpRepSquare(spc9type,spc9replist[1])
print

print "Vector * Adjoint: includes vector, antisymmetric 3-tensor"
DumpRepProduct(spc9type,spc9replist[1],spc9replist[3])
print

print "Vector * Spinor = vector-spinor + spinor"
DumpRepProduct(spc9type,spc9replist[1],spc9replist[2])
print

print "Go from M-theory to string theory: 9D -> 8D + 1D"
spc9d8 = MakeExtensionSplitter(spc9type,4)
print spc9d8
print

print "Vector -> vector + scalar"
DumpBranching(spc9d8,spc9replist[1])
print

print "Spinor -> two spinors"
DumpBranching(spc9d8,spc9replist[2])
print

print "Antisymmetric 2-tensor -> AS tensor + vector"
DumpBranching(spc9d8,spc9replist[3])
print

print "Symmetric traceless 2-tensor -> Symm tcl 2-tensor + vector + scalar"
DumpBranching(spc9d8,spc9replist[4])
print

print "Antisymmetric 3-tensor -> AS 3-tensor + AS 2-tensor"
DumpBranching(spc9d8,spc9replist[5])
print

PrintSubheader("10-D Space-Time: SO(10)")

spc10type = (4,5)
spc10replist = ((0, 0, 0, 0, 0), (1, 0, 0, 0, 0), (0, 0, 0, 1, 0), (0, 0, 0, 0, 1), \
	(0, 1, 0, 0, 0), (2, 0, 0, 0, 0), (1, 0, 0, 1, 0), (1, 0, 0, 0, 1))
spc10lbllist = ("Scalar", "Vector", "Spinor 1", "Spinor 2", "Adjoint; AS 2-tensor", \
	"Symm tls 2-tensor", "Vector-spinor 1", "Vector-spinor 2")
spc10list = zip(spc10lbllist,spc10replist)

DumpRepListProperties(spc10type,spc10list)
print

print "Vector * Vector = symmetric traceless 2-tensor +"
print "antisymmetric 2-tensor + scalar"
DumpRepSquare(spc10type,spc10replist[1])
print

print "Vector * Spinor = vector-spinor + spinor"
DumpRepProduct(spc10type,spc10replist[1],spc10replist[2])
print

print "Compactification: breakdown of 10D into 6D + 4D"
spc10d64 = MakeExtensionSplitter(spc10type,3)
print spc10d64
print

print "Vector -> 6D vector + 6D scalar, 4D vector"
DumpBranching(spc10d64,spc10replist[1])
print

print "Spinor -> Two 6D spinors, each paired with a 4D spinor"
DumpBranching(spc10d64,spc10replist[2])
print

print "Antisymmetric 2-tensor -> 6D AS 2-tensor + 4D AS 2-tensor +"
print "6D vector * 4D vector"
DumpBranching(spc10d64,spc10replist[4])
print

print "Symmetric traceless 2-tensor -> 6D and 4D Symm tcl 2-tensors +"
print "6D vector * 4D vector + scalar"
DumpBranching(spc10d64,spc10replist[5])
print

print "Vector-spinor -> 6D vector-spinor * 4D spinor + 6D spinor * 4D vector-spinor +"
print "6D spinor * 4D spinor"
DumpBranching(spc10d64,spc10replist[6])
print

PrintSubheader("11-D Space-Time: SO(11)")

spc11type = (2,5)
spc11replist = ((0, 0, 0, 0, 0), (1, 0, 0, 0, 0), (0, 0, 0, 0, 1), (0, 1, 0, 0, 0), \
	(2, 0, 0, 0, 0), (0, 0, 1, 0, 0), (1, 0, 0, 0, 1))
spc11lbllist = ("Scalar", "Vector", "Spinor", "Adjoint; AS 2-tensor", \
	"Symm tls 2-tensor", "AS 3-tensor", "Vector-spinor")
spc11list = zip(spc11lbllist,spc11replist)

DumpRepListProperties(spc11type,spc11list)
print

print "Vector * Vector = symmetric traceless 2-tensor +"
print "antisymmetric 2-tensor + scalar"
DumpRepSquare(spc11type,spc11replist[1])
print

print "Vector * Adjoint: includes vector, antisymmetric 3-tensor"
DumpRepProduct(spc11type,spc11replist[1],spc11replist[3])
print

print "Vector * Spinor = vector-spinor + spinor"
DumpRepProduct(spc11type,spc11replist[1],spc11replist[2])
print

print "Go from M-theory to string theory: 11D -> 10D + 1D"
spc11d8 = MakeExtensionSplitter(spc11type,5)
print spc11d8
print

print "Vector -> vector + scalar"
DumpBranching(spc11d8,spc11replist[1])
print

print "Spinor -> two spinors"
DumpBranching(spc11d8,spc11replist[2])
print

print "Antisymmetric 2-tensor -> AS tensor + vector"
DumpBranching(spc11d8,spc11replist[3])
print

print "Symmetric traceless 2-tensor -> Symm tcl 2-tensor + vector + scalar"
DumpBranching(spc11d8,spc11replist[4])
print

print "Antisymmetric 3-tensor -> AS 3-tensor + AS 2-tensor"
DumpBranching(spc11d8,spc11replist[5])
print

PrintHeader("Gauge and Flavor Symmetries")

PrintSubheader("QCD Gauge Symmetry and Light-Quark Flavor Symmetry: SU(3)")

print "They are treated together here because they have the same symmetry group."
print "The QCD symmetry is exact, and it is for three \"color\" states."
print "The flavor symmetry is approximate, and it is for up, down, and strange."
gfssu3type = (1,2)
gfssu3replist = ((0, 0), (1, 0), (0, 1), (1, 1), (3, 0), (0, 3))
gfssu3lbllist = ("Scalar: color, meson flavor",
	"Fundamental: quark color, flavor", \
	"Cjg fund.: a-quark color, flavor", \
	"Adjoint: gluon color, meson, spin-1/2 (a-)baryon flavor", \
	"Spin-3/2 baryon flavor", \
	"Spin-3/2 a-baryon flavor")
gfssu3list = zip(gfssu3lbllist,gfssu3replist)

DumpRepListProperties(gfssu3type,gfssu3list)
print
DumpLblRepList(gfssu3type,gfssu3list)
print "All of these states have degen 1, except for the neutral one in the adjoint,"
print "with degen 2. That corresponds to the Cartan subalgebra,"
print "used for designating states with eigenvalues."
print "In flavors, it makes the lambda-0 and sigma-0 (anti)baryons."
print

print "Mesons"
print "Inputs: 3 (1,0), 3* (0,1)"
print "Gives 8 (1,1): flavor; 1 (0,0): color, flavor"
DumpRepProduct(gfssu3type,gfssu3replist[1],gfssu3replist[2])
print

print "Baryons"
print "Inputs: 3 of 3 (1,0)"
print "Gives 10 (3,0): spin-3/2 flavor, two of 8 (1,1): spin-1/2 flavor,"
print "and 1 (0,0): color"
DumpRepCube(gfssu3type,gfssu3replist[1])
print

print "Antibaryons"
print "Inputs: 3 of 3* (0,1)"
print "Gives 10 (0,3): spin-3/2 flavor, two of 8 (1,1): spin-1/2 flavor,"
print "and 1 (0,0): color"
DumpRepCube(gfssu3type,gfssu3replist[2])
print


print "Decomposition of spin and flavor together."
print "They are, in order, symmetric, mixed and, antisymmetric."
print "The observed baryon states are the symmetric ones,"
print "which was contrary to what one would expect of a spin-1/2 particle."
print "One expects Fermi-Dirac statistics - antisymmetry -"
print "rather than Bose-Einstein statistics - symmetry."
print
print "Quarks: spin 1/2 (1,), flavor 3 (1,0)"
print
spinsu2type = (1,1)
spinrep = (1,)
cmbntype = (spinsu2type, gfssu3type)
cmbnrep = (spinrep, gfssu3replist[1])
daprres = DecomposeAlgProdRepPower(cmbntype, cmbnrep, 3)
print "Symmetric"
for rpr in daprres[0]: print rpr
print
print "Mixed"
for rpr in daprres[1]: print rpr
print
print "Antisymmetric"
for rpr in daprres[2]: print rpr
print

print "Decomposition of spin, flavor, and color"
print "together with Fermi-Dirac statistics - antisymmetry."
print "Note that when the last part is antisymmetry (0,0)"
print "the rest of it is what one finds in the observed baryons."
cmbntype = (spinsu2type, gfssu3type, gfssu3type)
cmbnrep = (spinrep, gfssu3replist[1], gfssu3replist[1])
daprres = DecomposeAlgProdRepPwrSym(cmbntype, cmbnrep, 3, -1)
for rpr in daprres: print rpr
print

PrintSubheader("Separation of Up/Down from Strange Quarks: SU(3) -> SU(2) * U(1)")

print "Decompose three-quark flavor SU(3) into:"
print "Isospin (up and down: SU(2), listed as \"angular momentum\" I)"
print "Hypercharge (up/down and strange : U(1), listed as 1/2 the U(1) factor)"
print "Electric charge is (isospin eigenvalue: -I to +I) + (hypercharge)"
fssu3d2 = MakeRootDemoter(gfssu3type,2)
print fssu3d2
print

print "Quark: u,d: (1/2,1/6); s: (0,-1/3)"
DumpBranching(fssu3d2,gfssu3replist[1])
print

print "A-quark: u*,d*: (1/2,-1/6); s*: (0,1/3)"
DumpBranching(fssu3d2,gfssu3replist[2])
print

print "Quark/a-quark or 3-(a-)quark octet"
print "   Pion - sigma - a-sigma - (1,0)"
print "   Kaon - xi - a-nucleon - (1/2,-1/2)"
print "   Kaon - nucleon - a-xi - (1/2,1/2)"
print "   Eta - lambda - a-lambda - (0,0)"
DumpBranching(fssu3d2,gfssu3replist[3])
print

print "Quark decuplet:"
print "   Delta - (3/2,1/2)"
print "   Sigma-x - (1,0)"
print "   Xi-x - (1/2,-1/2)"
print "   Omega - (0,-1)"
DumpBranching(fssu3d2,gfssu3replist[4])
print

print "A-quark decuplet:"
print "   A-delta - (3/2,-1/2)"
print "   A-sigma-x - (1,0)"
print "   A-xi-x - (1/2,1/2)"
print "   A-omega - (0,1)"
DumpBranching(fssu3d2,gfssu3replist[5])
print

PrintSubheader("Georgi-Glashow: SU(5)")

print """
The unbroken (Minimal Supersymmetric) Standard Model has symmetry
SU(3) * SU(2) * U(1), corresponding to QCD (color),
weak isospin, and weak hypercharge.
The elementary fermions (ElFm) and Higgses have left-handed and
right-handed parts, each one the antiparticle of the other.
A multiplet's electric-charge values: (-I, -I+1, ... I-1, I) + WHC,
where I is the WIS "angular momentum".
The multiplets are labeled with
(QCD multiplicity, WIS multiplicity: 2I+1, WHC: (-5/3)*(U(1) factor)):
Gauge: gluon, W, B
   g - SU(3) QCD - (8,1,0) - gluon
   W - SU(2) WIS - (1,3,0)
   B - U(1) WHC - (1,1,0)
Higgs:
   Hu - (1,2,1/2) / Hu* - (1,2,-1/2) - up Higgs
   Hd - (1,2,-1/2) / Hd* - (1,2,1/2) - down Higgs
Elementary fermions:
   Q - (3,2,1/6) / Q* - (3*,2,-1/6) - quark
   U - (3,1,2/3) / U* - (3*,1,-2/3) - up
   D - (3,1,-1/3) / D* - (3*,1,1/3) - down
   L - (1,2,-1/2) / L* - (1,2,1/2) - lepton
   N - (1,1,0) / N* - (1,1,0) - neutrino
   E - (1,1,-1) / E* - (1,1,1) - electron
Extras:
   HQ, HQ' - (3,1,-1/3) / HQ*, HQ'* - (3*,1,1/3) - \"Higgs quark\"
   LQ - (3,1,-5/6) / LQ* (3*,1,5/6) - gauge leptoquark
"""

gssu5type = (1,4)
gssu5replist = \
	((0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (1, 0, 0, 1))
gssu5lbllist = ("Scalar: 1 L and R ElFm", "5: R ElFm", "10 L ElFm", \
	"10* R ElFm", "5* L ElFm", "Adjoint: 24 Gauge")
gssu5list = zip(gssu5lbllist,gssu5replist)

DumpRepListProperties(gssu5type,gssu5list)
print

gssu5d321 = MakeRootDemoter(gssu5type,3)
print gssu5d321
print

print "1 L of ElFm: N*, 1 R of ElFm: N -- right-handed neutrinos"
print

print "5 R of ElFm: D, L*"
print "5 R of Higgs: HQ, Hd*"
print "5 L of Higgs: HQ'*, Hu"
DumpBranching(gssu5d321,gssu5replist[1])
print

print "10 L of ElFm: U*, Q, E*"
DumpBranching(gssu5d321,gssu5replist[2])
print

print "10* L of ElFm: U, Q*, E"
DumpBranching(gssu5d321,gssu5replist[3])
print

print "5* L of ElFm: D*, L"
print "5* L of Higgs: HQ*, Hd"
print "5* R of Higgs: HQ', Hu*"
DumpBranching(gssu5d321,gssu5replist[4])
print

print "24 of Gauge: g, LQ*, LQ, W, B"
DumpBranching(gssu5d321,gssu5replist[5])
print

PrintSubheader("Pati-Salam: SU(4) * SU(2) * SU(2)")

print "Alternative to Georgi-Glashow;"
print "can be expressed as SO(6) * SO(4): color-left-right."
print "The SU(4) breaks down into SU(3)*U(1) and the second SU(2) to U(1)."
print "The SU(3) makes QCD, the first SU(2) weak isospin,"
print "and the weak hypercharge comes from the U(1) factors:"
print "(2/3)*(first) + (second)."
print

gssu4type = (1,3)
gssu4replist = \
	((0, 0, 0), (1, 0, 0), (0, 0, 1), (1, 0, 1), (0, 1, 0))
gssu4lbllist = ("Scalar: Higgs", "Fund: ElFm L,R", "Cjg Fund: ElFm L,R", "Adjoint: Gauge", "Vector")
gssu4list = zip(gssu4lbllist,gssu4replist)

DumpRepListProperties(gssu4type,gssu4list)
print

gssu4d3 = MakeRootDemoter(gssu4type,3)
print gssu4d3
print

print "Higgs particles are (1,2,2) -- L: Hu, Hd -- R: Hu*, Hd*"
print

print "Elementary fermions:"
print "(4,2,1) L: (Q,L)"
print "(4,1,2) R: ((U,D),(N,E))"
DumpBranching(gssu4d3,gssu4replist[1])
print

print "Elementary fermions:"
print "(4*,2,1) R: (Q,L) all *"
print "(4*,1,2) L: ((D,U),(E,N)) all *"
DumpBranching(gssu4d3,gssu4replist[2])
print

print "Gauge: (15,1,1) + (1,3,1) + (1,1,3)"
print "The (15,1,1) decomposes into the gluons, some quark-like particles,"
print "and a B-like particle."
print "The (1,3,1) becomes the W, while the (1,1,3) decomposes into a B-like particle"
print "and a charged version."
DumpBranching(gssu4d3,gssu4replist[3])
print

print "Extra from SO(10) : (6,1,1) -- extra \"Higgs quark\" here also."
DumpBranching(gssu4d3,gssu4replist[4])
print

PrintSubheader("SO(10) to Georgi-Glashow")

print "Decomposes into SU(5) * U(1)."
print "The U(1) factor combined with hypercharge makes (B-L):"
print "(baryon number) - (lepton number)."
print "B-L = (4/5)*(- (U(1) factor) + (hypercharge))."

gsso10type = (4,5)
gsso10replist = \
	((0, 0, 0, 0, 0), (1, 0, 0, 0, 0), (0, 0, 0, 1, 0), (0, 0, 0, 0, 1), (0, 1, 0, 0, 0))
gsso10lbllist = ("Scalar", "Vector: Higgs", "Spinor 1: ElFm R", "Spinor 2: ElFm L", \
	"Adjoint: Gauge")
gsso10list = zip(gsso10lbllist,gsso10replist)

DumpRepListProperties(gsso10type,gsso10list)
print

gsso10d5 = MakeRootDemoter(gsso10type,5)
print gssu4d3
print

print "10 Higgs: 5, 5* for both L and R"
DumpBranching(gsso10d5,gsso10replist[1])
print

print "16 L ElFm: 10, 5*, 1"
DumpBranching(gsso10d5,gsso10replist[2])
print

print "16* R ElFm: 5, 10*, 1"
DumpBranching(gsso10d5,gsso10replist[3])
print

print "45 Gauge: 10*, 24, 10, 1"
print "The 24 is the SU(5) gauge multiplet, while the 1, 10, and 10* are extras."
DumpBranching(gsso10d5,gsso10replist[4])
print

PrintSubheader("SO(10) to Pati-Salam")

gsso10d64 = MakeExtensionSplitter(gsso10type,3)
print gsso10d64
print

print "10 Higgs: (6,1,1) + (1,2,2)"
DumpBranching(gsso10d64,gsso10replist[1])
print

print "16 L ElFm: (4*,1,2) + (4,2,1)"
DumpBranching(gsso10d64,gsso10replist[2])
print

print "16* R ElFm (4,1,2) + (4*,2,1)"
DumpBranching(gsso10d64,gsso10replist[3])
print

print "45 Gauge: (6,2,2) + (15,1,1) + (1,3,1) + (1,1,3)"
print "The last three are the Pati-Salam gauge multiplets."
DumpBranching(gsso10d64,gsso10replist[4])
print

print "**"
print "** E(6) to SO(10)"
print "**"
print
print "E(6) offers the possibility of unifying the Higgs with the elementary fermions."
print "That could explain the generation hierarchy, but the other Higgses must be made"
print "very massive by symmetry breaking."
print

gse6type = (5,6)
gse6replist = \
	((0, 0, 0, 0, 0, 0), (1, 0, 0, 0, 0, 0), (0, 0, 0, 0, 1, 0), (0, 0, 0, 0, 0, 1))
gse6lbllist = ("Scalar", "Fund 1: ElFm R", "Fund 2: ElFm L", "Adjoint: Gauge")
gse6list = zip(gse6lbllist,gse6replist)

DumpRepListProperties(gse6type,gse6list)
print

print "Cube of the fundamental rep"
print "Note that the symmetric part contains a scalar"
print "Thus making it fit with SUSY-style interactions"
print
DumpRepCube(gse6type,gse6replist[1])
print

gse6d10 = MakeRootDemoter(gse6type,5)
print gse6d10
print

print "27 R ElFm, Higgs: 16* + 10 + 1"
DumpBranching(gse6d10,gse6replist[1])
print

print "27* L ElFm, Higgs: 10 + 16 + 1"
DumpBranching(gse6d10,gse6replist[2])
print

print "78 Gauge: 45 + 16 + 16* + 1"
DumpBranching(gse6d10,gse6replist[3])
print

PrintSubheader("E(6) to Trinification")

print "That's three SU(3)'s."
print

gse6d330 = MakeExtensionSplitter(gse6type,3)
# To make the highest weights look more symmetric
gse6d33 = BrancherConjugate(gse6d330,(1,))
print gse6d33
print

print "Notice how the 27 and the 27* are mirror images of each other,"
print "and notice also their cyclicity."
print

print "27 R ElFm, Higgs: (3*,1,3) + (3,3*,1) + (1,3,3*)"
DumpBranching(gse6d33,gse6replist[1])
print

print "27* L ElFm, Higgs: (3,1,3*) + (1,3*,3) + (3*,3,1)"
DumpBranching(gse6d33,gse6replist[2])
print

print "78 Gauge: (3*,3*,3*) + (3,3,3) + (1,8,1) + (8,1,1) + (1,1,8)"
print "The trinification adjoint with some nicely-symmetric extras."
DumpBranching(gse6d33,gse6replist[3])
print

PrintSubheader("E(6) to SU(6)")

gse6d6 = MakeExtensionSplitter(gse6type,6)
print gse6d6
print

print "27 R ElFm, Higgs: (6,2) + (15*,1)"
DumpBranching(gse6d6,gse6replist[1])
print

print "27* L ElFm, Higgs: (6*,2) + (15,1)"
DumpBranching(gse6d6,gse6replist[2])
print

print "78 Gauge: (35,1) + (1,3) + (20,2)"
DumpBranching(gse6d6,gse6replist[3])
print

PrintSubheader("E(8) Splits")

print "This is part of the HE heterotic string;"
print "it has an E(8)*E(8) gauge symmetry."
print "One E(8) representation can contain *all* of the Standard Model's particles."
print "The spin-0 SM particles result from compactification,"
print "reduction of 10 space-time dimensions into 4 large and 6 small dimensions."
print
print "E(8) has the curious property of being the only Lie algebra"
print "whose fundamental and adjoint representations are the same."

gse8type = (5,8)
gse8replist = ((0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 1, 0))
gse8lbllist = ("Scalar", "Fund/Adj")
gse8list = zip(gse8lbllist,gse8replist)

DumpRepListProperties(gse8type,gse8list)
print

gse8d63 = MakeExtensionSplitter(gse8type,6)
print gse8d63
print

print "(27*,3) + (27,3*) + (78,1) + (1,8)"
print "Elementary fermions / Higgs, gauge, extra"
DumpBranching(gse8d63,gse8replist[1])
print

gse8d104 = MakeExtensionSplitter(gse8type,5)
print gse8d104
print

print "(16,4) + (45,1) + (16*,4*) + (10,6)"
print "Elementary fermions, gauge, Higgs"
DumpBranching(gse8d104,gse8replist[1])
print

gse8d55 = MakeExtensionSplitter(gse8type,4)
print gse8d55
print

print "(1,24) + (10,5*) + (5*,10*) + (24,0) + (5,10) + (10*,5)"
print "Georgi-Glashow + big extra"
DumpBranching(gse8d55,gse8replist[1])
print

PrintHeader("Mathematical Curiosities")

print "John Baez, in his treatise on octonions,"
print "http://math.ucr.edu/home/baez/octonions/"
print "mentions that G(2) is a subalgebra of SO(7)/B(3), which is in turn"
print "a subalgebra of SO(8)/D(4). This section will demonstrate those relations."
print

# Adds scalar to beginning
def sclim(n):
	return [n*[0]] + identmat(n)

jbg2type = (7,2)
jbso7type = (2,3)
jbso8type = (4,4)
jbg2replist = sclim(2)
jbso7replist = sclim(3)
jbso8replist = sclim(4)
jbg2lbllist = ("Scalar", "Adjoint", "Fundamental")
jbso7lbllist = ("Scalar", "Vector", "Adjoint", "Spinor")
jbso8lbllist = ("Scalar", "Vector", "Adjoint", "Spinor-1", "Spinor-2")
jbg2list = zip(jbg2lbllist,jbg2replist)
jbso7list = zip(jbso7lbllist,jbso7replist)
jbso8list = zip(jbso8lbllist,jbso8replist)

DumpRepListProperties(jbg2type,jbg2list)
print

DumpRepListProperties(jbso7type,jbso7list)
print

DumpRepListProperties(jbso8type,jbso8list)
print

print "Cube of the fundamental rep of G2"
print "Note that the antisymmetric part contains a scalar"
print
DumpRepCube(jbg2type,jbg2replist[2])
print

print "SO(7) to G(2):"
jbb3g2 = SubalgExtra("B3G2")
print jbb3g2
print

print "Vector -> fund"
DumpBranching(jbb3g2,jbso7replist[1])
print

print "Adjoint -> adjoint + fund"
DumpBranching(jbb3g2,jbso7replist[2])
print

print "Spinor -> fund + scalar"
DumpBranching(jbb3g2,jbso7replist[3])
print

print "SO(8) to SO(7):"
jbd4b3 = SubalgSOEvenOdd(4,0)
print jbd4b3
print

print "Vector -> vector + scalar"
DumpBranching(jbd4b3,jbso8replist[1])
print

print "Adjoint -> adjoint + vector"
DumpBranching(jbd4b3,jbso8replist[2])
print

print "Spinor 1 -> spinor"
DumpBranching(jbd4b3,jbso8replist[3])
print

print "Spinor 2 -> spinor"
DumpBranching(jbd4b3,jbso8replist[4])
print

print "SO(8) to G(2):"
jbd4g2 = SubalgExtra("D4G2")
print jbd4g2
print

print "Vector -> fund + scalar"
DumpBranching(jbd4g2,jbso8replist[1])
print

print "Adjoint -> adjoint + fund"
DumpBranching(jbd4g2,jbso8replist[2])
print

print "Spinor 1 -> fund + scalar"
DumpBranching(jbd4g2,jbso8replist[3])
print

print "Spinor 2 -> fund + scalar"
DumpBranching(jbd4g2,jbso8replist[4])
print

print
print "G(2) itself breaks down into some other algebras,"
print "namely, SU(2)/SO(3)/A(1) and SU(3)/A(2)."
print

print "Demote the first root to a U(1); what remains is SU(2)"
jbg2rd1 = MakeRootDemoter(jbg2type,1)
print jbg2rd1
print

print "G(2) Adjoint"
DumpBranching(jbg2rd1,jbg2replist[1])
print

print "G(2) Fundamental"
DumpBranching(jbg2rd1,jbg2replist[2])
print

print "Demote the second root to a U(1); what remains is SU(2)"
jbg2rd2 = MakeRootDemoter(jbg2type,2)
print jbg2rd2
print

print "G(2) Adjoint"
DumpBranching(jbg2rd2,jbg2replist[1])
print

print "G(2) Fundamental"
DumpBranching(jbg2rd2,jbg2replist[2])
print

print "Turn the whole algebra into a SU(2)"
jbg2ht = SubalgHeightA1(jbg2type)
print jbg2ht
print

print "G(2) Adjoint"
DumpBranching(jbg2ht,jbg2replist[1])
print

print "G(2) Fundamental"
DumpBranching(jbg2ht,jbg2replist[2])
print

print "Split at the first root, turning it into SU(2) * SU(2)"
jbg2sp1 = MakeExtensionSplitter(jbg2type,1)
print jbg2sp1
print

print "G(2) Adjoint"
DumpBranching(jbg2sp1,jbg2replist[1])
print

print "G(2) Fundamental"
DumpBranching(jbg2sp1,jbg2replist[2])
print

print "Split at the second root, turning it into SU(3)"
jbg2sp2 = MakeExtensionSplitter(jbg2type,2)
print jbg2sp2
print

print "G(2) Adjoint -> 8 + 3 + 3*"
DumpBranching(jbg2sp2,jbg2replist[1])
print

print "G(2) Fundamental -> 3 + 3* + 1"
DumpBranching(jbg2sp2,jbg2replist[2])
print


print "In 1949, Giulio Racah attempted to analyze"
print "the f electrons of lanthanides and actinides."
print "He created a cascade of algebras: SU(7) -> SO(7) -> G(2) -> SU(2)."
print "I will show here what his work is equivalent to:"
print "finding the height subalgebra of SU(7)."
print

gr1 = SubalgSUSO(7)
gr2 = SubalgExtra("B3G2")
gr3 = SubalgExtra("G2A1")
gr12 = ConcatBranchers(gr1,1,gr2)
gr123 = ConcatBranchers(gr12,1,gr3)
grh = SubalgHeightA1((1,6))

def DumpStsms(ld):
	print ld.latype
	for stsm in ld.stsms:
		print stsm[0]
		print stsm[1]

# Necessary to avoid comparing numpy arrays, since they are not hashable
def ReduceStsms(stsms):
	return tuple((tuple(stsm[:2]) for stsm in stsms))

print "Original algebra, sets of subalgebra and projection matrix:"
print "Giulio Racah's chain"
DumpStsms(gr123)
print "Height subalgebra"
DumpStsms(grh)
print "Equal?", ReduceStsms(gr123.stsms) == ReduceStsms(grh.stsms)
print

grgh = SubalgHeightA1((7,2))
print "G2A1 is a height subalgebra:"
print "G2A1"
DumpStsms(gr3)
print "Height subalgebra"
DumpStsms(grgh)
print "Equal?", ReduceStsms(gr3.stsms) == ReduceStsms(grgh.stsms)
print