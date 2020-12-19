# distutils: language = c++

from KWD_Histogram2D cimport Histogram2D as PyHistogram2D

from KWD_Histogram2D cimport Solver as PySolver

import numpy as np


cdef class Histogram2D:
    cdef PyHistogram2D mu

    def __cinit__(self):
        self.mu = PyHistogram2D()

    def __cinit__(self, n, X, Y, W):
        if not X.flags['C_CONTIGUOUS']:
            X = np.ascontiguousarray(X) # Makes a contiguous copy of the numpy array.
        cdef int[::1] Xmv = X

        if not Y.flags['C_CONTIGUOUS']:
            Y = np.ascontiguousarray(Y) # Makes a contiguous copy of the numpy array.
        cdef int[::1] Ymv = Y

        if not W.flags['C_CONTIGUOUS']:
            W = np.ascontiguousarray(W) # Makes a contiguous copy of the numpy array.
        cdef double[::1] Wmv = W

        self.mu = PyHistogram2D(n, &Xmv[0], &Ymv[0], &Wmv[0])

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

    def distance(self, Histogram2D A, Histogram2D B, L):
        return self.m.distance(A.mu, B.mu, L)

    def dense(self, Histogram2D A, Histogram2D B):
        return self.m.dense(A.mu, B.mu)

    def runtime(self):
        return self.m.runtime()
