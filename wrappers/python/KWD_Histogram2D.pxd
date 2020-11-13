cdef extern from "KWD_Histogram2D.cpp":
    pass

# Declare the class with cdef
cdef extern from "KWD_Histogram2D.h" namespace "KWD":
    cdef cppclass Histogram2D:
        Histogram2D() except +
        void add(int i, int j, double _weight)
        int size()
        double balance()
        void normalize()

    cdef cppclass Solver:
        Solver() except +
        #double distance(Histogram2D A, Histogram2D B, int L)
