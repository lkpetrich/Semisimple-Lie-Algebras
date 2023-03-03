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
- Semisimple Lie Algebras.nb
- Semisimple Lie Algebras - Notable Physics Results.nb
- Semisimple Lie Algebras - Test Suite.nb

## Python
- SemisimpleLieAlgebras.py
- SSLA_NotablePhysicsResults.py

## Notable Physics Results

The "notable physics results" are a sort of "greatest hits" of physics results that use Lie algebras, including results for different space-time dimensions, gauge symmetries, and flavor symmetries. Also included are some mathematical curiosities. More specifically:
- Scalars, spinors, vectors, vector-spinors, and tensors in 4D, 10D, 11D, and also transverse dimensions 2D, 8D, 9D
- Quark color and light-flavor symmetries; split of strange quarks from up and down ones
- Grand unified theories: SU(5), SO(10), E6, E8, ...
- Mathematical curiosities: John Baez and G2, Giulio Racah and G2, ...

Also included is some software for working with "Young diagrams", some nice graphical representations of SU(n) representations. Can do rep products (Littlewood-Richardson), and decomposition of reps into orbits (Kostka).
- Young Diagrams for SU(n) Multiplets.nb
- YoungDiagrams.py

## C++ version in directory CPP

It uses the GNU Multiple Precision library for the total-degeneracy size, since some total degeneracies can be surprisingly large, surpsising since they are for small highest-weight values. For E8, all zero but one gives degeneracies ranging from 248 to 6,899,079,264.

In "Makefile", set variable GMPDIR to wherever the GNU Multiple Precision library resides on your system
Build with "make"

### Files
- Fraction.h - Template class and functions for fractions. Use an integer class in it or some other class where Euclid's algorithm is meaningful, like an integer-coefficient polynomial ring.
- Test_Fraction.cpp - for testing it
- LinearAlgebra.h -- Template classes for matrix and vector set, template functions for vector and matrix operations: elementary ones, matrix inversion.
- Test_LinearAlgebra.cpp
- LieAlgebra.h LieAlgebra.cpp - the Lie algebra itself: the Dynkin diagram, metric and Cartan matrices, and positive roots
- Test_LieAlgebra.cpp - for testing it
- LieAlgRep.h LieAlgRep.cpp - Lie-algebra representations: irreducible ones, the total sizes of irreps (what uses GMP), reps of product algebras, generic rep handler for finding which irrep are present in some rep.
- Test_LieAlgRep.cpp - for testing it
- LieAlgRepProd.h LieAlgRepProd.cpp - for products and powers of representations. The powers (plethysms) are broken down by symmetry type, like completely (anti)symmetric or a type designated by a Young diagram.
- Test_ LieAlgRepProd.h.cpp - for testing it

### Compile output
- Intermediate files *.o - don't really need to keep them
- libSSLA.so - shared library for all the Lie-algebra, rep, and product/power code

### Test programs
- Test_Fraction
- Test_LinearAlgebra
- Test_LieAlgebra
- Test_LieAlgRep
- Test_LieAlgRepProd
All of them but Test_Fraction do tests specified in command-line args. To find out which ones are available, run them without command-line args.

Its subdirectory Tests contains
- refcnt.cpp -- test of reference-counted pointer
- typelens.cpp -- lengths of integer types

## Extra files

Extra Mma Notebooks -- contains files for doing such extras as explicit constructions of various representations and attempts to find higher-power invariants.

Extra Python Files -- contains versions of the main Python files that use Numerical Python (NumPy): http://www.numpy.org/ Though their performance is disappointing, I am including them for the benefit of anyone who might want to experiment with them.