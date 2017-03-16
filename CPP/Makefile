# For building C++ semi-simple Lie-algebra code

# It uses the GNU Multiple Precision library (GMP)
# residing at directory GMPDIR
# It's sw here, the default root directory of OSX package manager fink

CC = g++ -O3

GMPDIR = /sw
GMPINC = $(GMPDIR)/include
GMPLIB = $(GMPDIR)/lib

all: libSSLA.so test

clean:
	/bin/rm *.o

libSSLA.so: LieAlgebra.o LieAlgRep.o LieAlgRepProd.o LieAlgBranching.o
	$(CC) -shared -o libSSLA.so LieAlgebra.o LieAlgRep.o \
		LieAlgRepProd.o LieAlgBranching.o \
		-L$(GMPLIB) -lgmpxx -lgmp

LieAlgebra.o: LieAlgebra.cpp \
	LieAlgebra.h LinearAlgebra.h Fraction.h
	$(CC) -c LieAlgebra.cpp

LieAlgRep.o: LieAlgRep.cpp LieAlgRep.h SmartPointer.h \
	LieAlgebra.h LinearAlgebra.h Fraction.h
	$(CC) -c LieAlgRep.cpp -I$(GMPINC)

LieAlgRepProd.o: LieAlgRepProd.cpp LieAlgRepProd.h LieAlgRep.h SmartPointer.h \
	LieAlgebra.h LinearAlgebra.h Fraction.h
	$(CC) -c LieAlgRepProd.cpp -I$(GMPINC)

LieAlgBranching.o: LieAlgBranching.cpp LieAlgBranching.h LieAlgRep.h SmartPointer.h \
	LieAlgebra.h LinearAlgebra.h Fraction.h
	$(CC) -c LieAlgBranching.cpp -I$(GMPINC)

# Representation, Rep Product, Rep Power, Branching

test: Test_Fraction Test_LinearAlgebra Test_LieAlgebra \
	Test_LieAlgRep Test_LieAlgRepProd Test_LieAlgBranching

Test_Fraction: Test_Fraction.cpp Fraction.h
	$(CC) -o Test_Fraction Test_Fraction.cpp

Test_LinearAlgebra: Test_LinearAlgebra.cpp LinearAlgebra.h Fraction.h
	$(CC) -o Test_LinearAlgebra Test_LinearAlgebra.cpp

Test_LieAlgebra: Test_LieAlgebra.cpp libSSLA.so
	$(CC) -o Test_LieAlgebra Test_LieAlgebra.cpp -L. -lSSLA

Test_LieAlgRep: Test_LieAlgRep.cpp libSSLA.so
	$(CC) -o Test_LieAlgRep Test_LieAlgRep.cpp -L. -lSSLA \
		-I$(GMPINC) -L$(GMPLIB) -lgmpxx -lgmp

Test_LieAlgRepProd: Test_LieAlgRepProd.cpp libSSLA.so
	$(CC) -o Test_LieAlgRepProd Test_LieAlgRepProd.cpp -L. -lSSLA \
		-I$(GMPINC) -L$(GMPLIB) -lgmpxx -lgmp

Test_LieAlgBranching: Test_LieAlgBranching.cpp libSSLA.so
	$(CC) -o Test_LieAlgBranching Test_LieAlgBranching.cpp -L. -lSSLA \
		-I$(GMPINC) -L$(GMPLIB) -lgmpxx -lgmp
