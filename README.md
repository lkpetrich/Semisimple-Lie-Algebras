# Semisimple Lie Algebras

What are semisimple Lie algebras?

These mathematical structures are important in physics for expressing the internal symmetries of space-time and elementary particles and the like. They have applications like getting Standard-Model particle multiplets from Grand-Unified-Theory particle multiplets.

I wrote this software because it was hard to find anything online that does certain things that I wanted to do. More specifically, decomposing product representations into irreducible ones and calculating branching rules: what representations of subalgebras does a representation of an algebra produce? It has applications like getting Standard-Model particle multiplets from Grand-Unified-Theory particle multiplets.

I have written this software in Mathematica, Python, and C++, and all three versions are feature-parallel except for Mathematica's graphics code. This software does:
- Setting up the algebras themselves
- Finding irreducible representations
- Decomposing product and representations
- Finding subalgebra representations: branching rules
- Worked examples of much of the code

## Mathematica
- `Semisimple Lie Algebras.nb`
- `Semisimple Lie Algebras - Notable Physics Results.nb`
- `Semisimple Lie Algebras - Test Suite.nb`

## Python
- `SemisimpleLieAlgebras.py`
- `SSLA_NotablePhysicsResults.py`

## Notable Physics Results

The "notable physics results" are a sort of "greatest hits" of physics results that use Lie algebras, including results for different space-time dimensions, gauge symmetries, and flavor symmetries. Also included are some mathematical curiosities. More specifically:
- Scalars, spinors, vectors, vector-spinors, and tensors in 4D, 10D, 11D, and also transverse dimensions 2D, 8D, 9D
- Quark color and light-flavor symmetries; split of strange quarks from up and down ones
- Grand unified theories: SU(5), SO(10), E6, E8, ...
- Mathematical curiosities: John Baez and G2, Giulio Racah and G2, ...

Also included is some software for working with "Young diagrams", some nice graphical representations of SU(n) representations. Can do rep products (Littlewood-Richardson), and decomposition of reps into orbits (Kostka).
- `Young Diagrams for SU(n) Multiplets.nb`
- `YoungDiagrams.py`

## C++ version in directory CPP

It uses the GNU Multiple Precision library for the total-degeneracy size, since some total degeneracies can be surprisingly large, surpsising since they are for small highest-weight values. For E8, all zero but one gives degeneracies ranging from 248 to 6,899,079,264.

In "Makefile", set variable GMPDIR to wherever the GNU Multiple Precision library resides on your system
Build with "make"

### Files
- `Fraction.h` - Template class and functions for fractions. Use an integer class in it or some other class where Euclid's algorithm is meaningful, like an integer-coefficient polynomial ring.
- `Test_Fraction.cpp` - for testing it
- `LinearAlgebra.h` -- Template classes for matrix and vector set, template functions for vector and matrix operations: elementary ones, matrix inversion.
- `Test_LinearAlgebra.cpp`
- `LieAlgebra.h LieAlgebra.cpp` - the Lie algebra itself: the Dynkin diagram, metric and Cartan matrices, and positive roots
- `Test_LieAlgebra.cpp` - for testing it
- `LieAlgRep.h LieAlgRep.cpp` - Lie-algebra representations: irreducible ones, the total sizes of irreps (what uses GMP), reps of product algebras, generic rep handler for finding which irrep are present in some rep.
- `Test_LieAlgRep.cpp` - for testing it
- `LieAlgRepProd.h LieAlgRepProd.cpp` - for products and powers of representations. The powers (plethysms) are broken down by symmetry type, like completely (anti)symmetric or a type designated by a Young diagram.
- `Test_ LieAlgRepProd.h.cpp` - for testing it

### Compile output
- Intermediate files `*.o` - don't really need to keep them
- `libSSLA.so` - shared library for all the Lie-algebra, rep, and product/power code

### Test programs
- `Test_Fraction`
- `Test_LinearAlgebra`
- `Test_LieAlgebra`
- `Test_LieAlgRep`
- `Test_LieAlgRepProd`
All of them but Test_Fraction do tests specified in command-line args. To find out which ones are available, run them without command-line args.

Its subdirectory Tests contains
- `refcnt.cpp` -- test of reference-counted pointer
- `typelens.cpp` -- lengths of integer types

# Young Diagrams

Young diagrams are a graphical technique that can be convenient for finding features of representations of Lie algebras of kinds of matrices: SU(n), SO(n), Sp(n) for size parameter n. These diagrams consist of boxes that extend rightward and downward from a top-left box. For representations:

Vector:  
[]  
Symmetric 2-tensor:  
[] []  
Antisymmetric 2-tensor:  
[]  
[]  
Mixed 3-tensor:  
[] []  
[]  

Notebook: `Young Diagrams.nb`

In it:
- Young-diagram operations: conversion to and from highest-weight vectors, and transpose
- Degeneracy: representation size
- Repreentation products for SU(n): Littlewood-Richardson rule
- Vector-Representation Powers for SU(n)
- Weyl Orbits for SU(n)
- Branching Rules: SU(n) to SO(n) and Sp(n) -- good for decomposing a general tensor rep into traceless ones, with SU(n) into SO(n)
- Young-Diagram Nesting
- Symmetric Functions and Related Ones
  - Symmetric and Alternating Group Characters
  - Symmetric Polynomials
  - Schur Polynomials - a kind of symmetric polynomial related to SU(n) algebra characters
  - Plethysms - Symmetrized Powers

# Lie-Algebra Matrices

Calculates explicit matrices for some representations of several Lie algebras, and also various associated quantities.

Notebook: `Lie-Algebra Matrices.nb`

The algebras: SU(n), SO(n), Sp(2n), G2, F4, E6, E7

The representations: vector (fundamental) ones for SU(n) and Sp(2n), vector and spinor ones for SO(2n), fundamental ones for G2, F4, E6, E7

Does interrelationships:
- SO(2),U(1) -- SO(1,1),GL(1,R+)
- SO(3),SU(2),SL(1,H),Sp(2) -- SO(2,1),SU(1,1),SL(2,R),Sp(2,R)
- SO(4),SU(2)xSU(2) -- SO(3,1),SL(2,C) -- SO(2,2),SU(1,1)xSU(1,1) -- SO(2,H),SU(2)xSU(1,1)
- SO(5),Sp(4) -- SO(4,1),Sp(2,2) -- SO(3,2),Sp(4,R)
- SO(6),SU(4) -- SO(5,1),SL(2,H) -- SO(4,2),SU(2,2) -- SO(3,3),SL(4,R) -- SO(3,H),SU(3,1)

Notes:
- R: over real numbers, C: over complex numbers, H: over quaternions (Pauli matrices)
- SL(n,H) is SU*(2n) and SO(n,H) is SO*(2n)
- GL(n,R) ~ SL(n,R) x GL(1,R) and U(n) ~ SU(n) x U(1)
- GL(n,C) ~ GL(n,R) x GL(n,R) and likewise for SL(n,C)

For the exceptional algebras, the notebook composes their matrices from other algebras' matrices.
- G2: octonion automorphism, SU(3)
- F4: SO(9)
- E6: F4, SU(3) x SU(3) x SU(3)
- E7: SU(8)
- E8: not done

Also:
- Finds combinations of an algebra's matrices that form the Cartan-Weyl basis for that algebra
- Finds some invariants for some of the algebra matrices, like 2-index and 3-index ones.

# Lie-Algebra Automorphisms

An automorphism is a transformation on an algebra that yields the same algebra. Inner automorphisms are generated by members of that algebra, while outer automorphisms are all others.

Notebook: `Lie-Algebra Automorphisms.nb`

Finds automorphism Lie algebras for these algebras:
- GL(n,R) - inner: SL(n,R), outer: scaling of GL(1,R) part
- SO(n), G2 - inner only
- Euclid-Poincaré algebras: inner: rotation part, outer: scaling of translation part

# Small Lie Algebras

Discusses Lie algebras with 1, 2, 3, 4, and 5 generators. Most of these algebras are not semisimple, though the 3-generator ones include two simple ones: SO(3) and SO(2,1)

Notebook: `Small Lie Algebras.nb`

It calculates:
- Commutator series: from algebra G, series G(i) starting with G(1) = [G,G]
  - Derived series: G(i+1) = [G(i),G(i)]
  - Lower central series: G(i+1) = [G,G(i)]
- Algebra metric
- Some invariants
- Some automorphisms
- Some ideals: subalgebra H of algebra G satisfying [G,H] = H

# Higher-Dimensional Supersymmetry and Supergravity

Supersymmetry and supergravity in at least four space-time dimensions, all the way up to eleven of them. Ten-dimensional supergravity theories are a part of supersymmetric-string ("superstring") theories, and the single eleven-dimensional theory is a part of a proposed unification of string theory called "M theory".

Notebook: `Higher-Dimensional Supersymmetry and Supergravity.nb`

Has the multiplet structure for every SUGRA theory, expressed in the format used in notebook Semisimple Lie Algebras.nb

Also has differences in multiplet content between theories, and generation of some theories by multiplet products.

It also has some attempts at defining SUSY gauge theories between four and ten dimensions, where such theories are most naturally defined.

# Virasoro and Super-Virasoro Infinite-Dimensional Algebras

Infinite-dimensional Lie algebras important in string theory.

Notebook: `Virasoro and Super-Virasoro algebras.nb`

The Witt-Virasoro algebra is this infinite-dimensional one, with generators L(n) satisfying commutator relation
```math
[L_n,L_m] = L_n \cdot L_m - L_m \cdot L_n = (n - m) L_{n+m} + (c_0 + c_1 n^2) \delta_{n+m,0}
```
The second term in the right expression is a "central charge".

The supersymmetric extension of this algebra has fermionic generators G(r), satisfying an anticommutation relation. By comparison, the L(n) are bosonic, like generators of Lie algebras in general.
```math
[L_n,G_r] = L_n \cdot G_r - G_r \cdot L_n = (n/2 - r) G_{n+r}
```
```math
\{G_r,G_s\} = G_r \cdot G_s + G_s \cdot G_r = 2 L_{r+s} + (c_0 + 4 c_1 r^2) \delta_{r+s,0}
```
Note the plus in the middle expression, the expansion of the anticommutator in the left expression.

This notebook is for verifying and deriving those expressions, with only partial success in doing that derivation.

# Extra files

Extra Mma Notebooks -- contains files for doing such extras as explicit constructions of various representations and attempts to find higher-power invariants.

Extra Python Files -- contains versions of the main Python files that use [Numerical Python (NumPy)](http://www.numpy.org/) Though their performance is disappointing, I am including them for the benefit of anyone who might want to experiment with them.
