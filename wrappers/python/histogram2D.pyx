# distutils: language = c++

from KWD_Histogram2D cimport Histogram2D as PyHistogram2D


cdef class Histogram2D:
    cdef PyHistogram2D mu

    def __cinit__(self):
        self.mu = PyHistogram2D()

    def add(self, i, j, w):
        return self.mu.add(i, j, w)

    def size(self):
        return self.mu.size()