I have written some files working for semisimple Lie algebras.

Setting up the algebras themselves
Finding irreducible representations
Decomposing product representations
Finding subalgebra representations: branching rules
Worked examples of much of the code


Mathematica notebooks:
Semisimple Lie Algebras.nb
Semisimple Lie Algebras - Notable Physics Results.nb
Semisimple Lie Algebras - Test Suite.nb
Young Diagrams for SU(n) Multiplets.nb

Python source files:
SemisimpleLieAlgebras.py
SSLA_NotablePhysicsResults.py

SemisimpleLieAlgebras_NumPy.py
SSLA_NP_NotablePhysicsResults.py



What they are:


* Semisimple Lie Algebras.nb
* SemisimpleLieAlgebras.py

Many functions for working with such Lie algebras, including functions for finding product representations and branching rules.

The second Python one uses Numerical Python (NumPy), at 


* Semisimple Lie Algebras - Notable Physics Results.nb
* SSLA_NotablePhysicsResults.py

A sort of "greatest hits" of physics results that use Lie algebras, including results for different space-time dimensions, gauge symmetries, and flavor symmetries. It also has some mathematical curiosities.

The Mathematica one depends on "Semisimple Lie Algebras.nb". Open and run it before running "Semisimple Lie Algebras - Notable Physics Results.nb".

The Python one depends on "SemisimpleLieAlgebras.py". It will automatically be used if it's in the same directory as "SSLA_NotablePhysicsResults.py" or in Python's search path. Run it at the command line; it needs no args. It will make its output at the command line; you may redirect it to a file or a text editor if you so choose.


* Young Diagrams for SU(n) Multiplets.nb
* YoungDiagrams.py

Useful for representations with general n. Can do rep products (Littlewood-Richardson), and decomposition of reps into orbits (Kostka).


C++ version in directory CPP

It uses the GNU Multiple Precision library for 

In "Makefile", set variable GMPDIR to wherever the GNU Multiple Precision library resides on your system
Build with "make"

The files:
Fraction.h - Template class and functions for fractions. Use an integer class in it or some other class where Euclid's algorithm is meaningful, like an integer-coefficient polynomial ring.

Test_Fraction.cpp - for testing it

LinearAlgebra.h -- Template classes for matrix and vector set, template functions for vector and matrix operations: elementary ones, matrix inversion.

Test_LinearAlgebra.cpp
LieAlgebra.h LieAlgebra.cpp - the Lie algebra itself: the Dynkin diagram, metric and Cartan matrices, and positive roots

Test_LieAlgebra.cpp - for testing it

LieAlgRep.h LieAlgRep.cpp - Lie-algebra representations: irreducible ones, the total sizes of irreps (what uses GMP), reps of product algebras, generic rep handler for finding which irrep are present in some rep.

Test_LieAlgRep.cpp - for testing it

LieAlgRepProd.h LieAlgRepProd.cpp - for products and powers of representations. The powers (plethysms) are broken down by symmetry type, like completely (anti)symmetric or a type designated by a Young diagram.

Test_ LieAlgRepProd.h.cpp - for testing it

Output:

Intermediate files *.o - don't really need to keep them

libSSLA.so - shared library for all the Lie-algebra, rep, and product/power code

Test programs:Test_FractionTest_LinearAlgebraTest_LieAlgebraTest_LieAlgRepTest_LieAlgRepProd

All of them but Test_Fraction do tests specified in command-line args. To find out which ones are available, run them without command-line args.

Its subdirectory Tests contains
refcnt.cpp -- test of reference-counted pointer
typelens.cpp -- lengths of integer types


Extra Mma Notebooks -- contains various extra sorts of stuff, like constructions of various representations and attempts to find higher-power invariants.

Extra Python Files -- contains versions of the main Python files that use Numerical Python (NumPy): http://www.numpy.org/
* SemisimpleLieAlgebras_NumPy.py
* SSLA_NP_NotablePhysicsResults.py
Though their performance is disappointing, I am including them for the benefit of anyone who might want to experiment with them.


Licensing: all this software is licensed under the MIT license, as stated in the file LICENSE.txt You can do whatever you want with it as long as you credit me with my code and as long as you accept that I make no legally binding guarantee that it is going to work.


-- Loren Petrich
