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
        double distance(const Histogram2D& A, const Histogram2D& B, int L)
