from libcpp.string cimport string
from libcpp.vector cimport vector

# Declare the class with cdef
cdef extern from "KWD_Histogram2D.h" namespace "KWD":
    cdef cppclass Histogram2D:
        Histogram2D() except +
        Histogram2D(int n, int* x, int* y, double* w) except +
        void add(int i, int j, double _weight)
        void update(int i, int j, double _weight)
        int size()
        double balance()
        void normalize()

    cdef cppclass Solver:
        Solver() except +
        double compareApprox(int _n, int* _Xs, int* _Ys, double* _W1, double* _W2, int _L)
        vector[double] compareApprox(int _n, int _m, int* _Xs, int* _Ys, double* _W1, double* _Ws, int _L)
        # double* compareApprox(int _n, int _m, int* _Xs, int* _Ys, double* _Ws, int _L)        
        double distance(const Histogram2D& A, const Histogram2D& B, int L)
        double column_generation(const Histogram2D& A, const Histogram2D& B, int L)
        double dense(const Histogram2D& A, const Histogram2D& B)
        double runtime()
        void setStrParam(string param, string value)
        void setDblParam(string param, float value)
        void getStrParam(string param)
        void getDblParam(string param)