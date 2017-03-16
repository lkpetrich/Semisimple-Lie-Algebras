#ifndef REFERENCE_COUNTING_SMART_POINTER
#define REFERENCE_COUNTING_SMART_POINTER

/*
	Reference-counting smart pointer
*/

template<class P> class RefCounter
{
	// Uses pointers to both the object and the count,
	// so that the count can be passed around alongside the object
	P *ObjPtr;
	int *CountPtr;
	
	// Add and subtract a reference
	void AddRef() {++(*CountPtr);}
	void SubRef() {if (--(*CountPtr) <= 0) {delete ObjPtr; delete CountPtr;}}
	
public:

	// Uses the default constructor
	RefCounter():
		ObjPtr(new P(0)), CountPtr(new int) {*CountPtr = 1;}

	// Gets an object that was allocated with "new"
	// Its pointer may then safely go out of scope
	// If it is a null pointer, the default constructor is used instead
	RefCounter(P *OrigPtr):
		ObjPtr(OrigPtr), CountPtr(new int)
			{if (!ObjPtr) ObjPtr = new P; *CountPtr = 1;}
	
	// Copy constructor
	RefCounter(const RefCounter<P> &SRCOrig):
		ObjPtr(SRCOrig.ObjPtr), CountPtr(SRCOrig.CountPtr) {AddRef();}
	
	// Assignment
	RefCounter<P> &operator = (const RefCounter<P> &SRCOrig)
	{
		// Don't do anything further if assigning to itself
		if (ObjPtr == SRCOrig.ObjPtr) return *this;
		// Will no longer refer to its original object`
		SubRef();
		// Refers to its arg's object
		ObjPtr = SRCOrig.ObjPtr;
		CountPtr = SRCOrig.CountPtr;
		AddRef();
		// Pass along for repeated assignment
		return *this;
	}
	
	// Destructor
	~RefCounter() {SubRef();}
	
	// Returns object by reference
	P &operator * () {return *ObjPtr;}
	const P &operator * () const {return *ObjPtr;}
	
	// For accessing members of the object
	P *operator -> () {return ObjPtr;}
	const P *operator -> () const {return ObjPtr;}
	
	// How many?
	int RefCount() const {return *CountPtr;}
	
	// Same pointer?
	bool operator == (RefCounter &RC) {return (ObjPtr == RC.ObjPtr);}
	
	// Test for the null pointer being present
	bool IsNull() {return ObjPtr == 0;}
};

#endif
