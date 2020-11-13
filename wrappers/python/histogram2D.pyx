# distutils: language = c++

from KWD_Histogram2D cimport Histogram2D as PyHistogram2D

from KWD_Histogram2D cimport Solver as PySolver

cdef class Histogram2D:
    cdef PyHistogram2D mu

    def __cinit__(self):
        self.mu = PyHistogram2D()

    def add(self, i, j, w):
        return self.mu.add(i, j, w)

    def size(self):
        return self.mu.size()

    def balance(self):
        return self.mu.balance()

    def normalize(self):
        return self.mu.normalize()

cdef class Solver:
    cdef PySolver m

    def __cinit__(self):
        self.m = PySolver()

    #def distance(self, A, B, L):
    #    return self.m.distance(A, B, L)
