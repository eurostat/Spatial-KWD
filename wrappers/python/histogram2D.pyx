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
            X = np.ascontiguousarray(X, dtype=int) # Makes a contiguous copy of the numpy array.
        cdef int[::1] Xmv = X

        if not Y.flags['C_CONTIGUOUS']:
            Y = np.ascontiguousarray(Y, dtype=int) # Makes a contiguous copy of the numpy array.
        cdef int[::1] Ymv = Y

        if not W.flags['C_CONTIGUOUS']:
            W = np.ascontiguousarray(W, dtype=float) # Makes a contiguous copy of the numpy array.
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
        self.m = PySolver()  # TODO: change the log frequency
        
    def compareExact(self, n, X, Y, W1, W2):
        if not X.flags['C_CONTIGUOUS']:
            X = np.ascontiguousarray(X, dtype=np.int32) # Makes a contiguous copy of the numpy array.
        cdef int[::1] Xmv = X

        if not Y.flags['C_CONTIGUOUS']:
            Y = np.ascontiguousarray(Y, dtype=np.int32) # Makes a contiguous copy of the numpy array.
        cdef int[::1] Ymv = Y

        if not W1.flags['C_CONTIGUOUS']:
            W1 = np.ascontiguousarray(W1, dtype=float) # Makes a contiguous copy of the numpy array.
        cdef double[::1] Wmv1 = W1

        if not W2.flags['C_CONTIGUOUS']:
            W2 = np.ascontiguousarray(W2, dtype=float) # Makes a contiguous copy of the numpy array.
        cdef double[::1] Wmv2 = W2

        return self.m.compareExact(n, &Xmv[0], &Ymv[0], &Wmv1[0], &Wmv2[0])


    def compareApprox(self, n, X, Y, W1, W2, L):
        if not X.flags['C_CONTIGUOUS']:
            X = np.ascontiguousarray(X, dtype=np.int32) # Makes a contiguous copy of the numpy array.
        cdef int[::1] Xmv = X

        if not Y.flags['C_CONTIGUOUS']:
            Y = np.ascontiguousarray(Y, dtype=np.int32) # Makes a contiguous copy of the numpy array.
        cdef int[::1] Ymv = Y

        if not W1.flags['C_CONTIGUOUS']:
            W1 = np.ascontiguousarray(W1, dtype=float) # Makes a contiguous copy of the numpy array.
        cdef double[::1] Wmv1 = W1

        if not W2.flags['C_CONTIGUOUS']:
            W2 = np.ascontiguousarray(W2, dtype=float) # Makes a contiguous copy of the numpy array.
        cdef double[::1] Wmv2 = W2

        return self.m.compareApprox(n, &Xmv[0], &Ymv[0], &Wmv1[0], &Wmv2[0], L)

    
    def compareApprox2(self, n, m, X, Y, W1, Ws, L):
        if not X.flags['C_CONTIGUOUS']:
            X = np.ascontiguousarray(X, dtype=np.int32) # Makes a contiguous copy of the numpy array.
        cdef int[::1] Xmv = X

        if not Y.flags['C_CONTIGUOUS']:
            Y = np.ascontiguousarray(Y, dtype=np.int32) # Makes a contiguous copy of the numpy array.
        cdef int[::1] Ymv = Y

        if not W1.flags['C_CONTIGUOUS']:
            W1 = np.ascontiguousarray(W1, dtype=float) # Makes a contiguous copy of the numpy array.
        cdef double[::1] Wmv1 = W1

        if not Ws.flags['C_CONTIGUOUS']:
            Ws = np.ascontiguousarray(Ws.flatten(), dtype=float) # Makes a contiguous copy of the numpy array.
        cdef double[::1] Wmvs = Ws

        return self.m.compareApprox(n, m, &Xmv[0], &Ymv[0], &Wmv1[0], &Wmvs[0], L)


    def compareApprox3(self, n, m, X, Y, Ws, L):
        if not X.flags['C_CONTIGUOUS']:
            X = np.ascontiguousarray(X, dtype=np.int32) # Makes a contiguous copy of the numpy array.
        cdef int[::1] Xmv = X

        if not Y.flags['C_CONTIGUOUS']:
            Y = np.ascontiguousarray(Y, dtype=np.int32) # Makes a contiguous copy of the numpy array.
        cdef int[::1] Ymv = Y

        if not Ws.flags['C_CONTIGUOUS']:
            Ws = np.ascontiguousarray(Ws.flatten(), dtype=float) # Makes a contiguous copy of the numpy array.
        cdef double[::1] Wmvs = Ws

        return self.m.compareApprox3(n, m, &Xmv[0], &Ymv[0], &Wmvs[0], L)
    
    def distance(self, Histogram2D A, Histogram2D B, L):
        return self.m.distance(A.mu, B.mu, L)

    def column_generation(self, Histogram2D A, Histogram2D B, L):
        return self.m.column_generation(A.mu, B.mu, L)

    def dense(self, Histogram2D A, Histogram2D B):
        return self.m.dense(A.mu, B.mu)

    def runtime(self):
        return self.m.runtime()
    
    def iterations(self):
        return self.m.iterations()

    def num_nodes(self):
        return self.m.num_nodes()

    def num_arcs(self):
        return self.m.num_arcs()

    def status(self):
        return self.m.status()
    
    def setStrParam(self, param, value):
        return self.m.setStrParam(param, value)

    def setDblParam(self, param, value):
        return self.m.setDblParam(param, value)

    def getStrParam(self, param):
        return self.m.getStrParam(param)

    def getDblParam(self, param):
        return self.m.getDblParam(param)

# Helper functions
import sys
def compareOneToOne(Coordinates, Weights, Options):
    L = Options.get('L', 3)  # L=3 default value
    method = Options.get('Method', 'Approx').encode('utf-8')
    model = Options.get('Model', 'mincostflow').encode('utf-8')
    algo = Options.get('Algorithm', 'colgen').encode('utf-8')
    verbosity = Options.get('Verbosity', 'info').encode('utf-8')
    time_limit = Options.get('TimeLimit', sys.maxsize)
    opt_tolerance = Options.get('OptTolerance', 1e-06)
    
    s = Solver()
    s.setStrParam('Method'.encode('utf-8'), method)
    s.setStrParam('Model'.encode('utf-8'), model)
    s.setStrParam('Algorithm'.encode('utf-8'), algo)
    s.setStrParam('Verbosity'.encode('utf-8'), verbosity)
    s.setDblParam('TimeLimit'.encode('utf-8'), time_limit)
    s.setDblParam('OptTolerance'.encode('utf-8'), opt_tolerance)
    
    n, m = Weights.shape
    X = Coordinates[:,0]
    Y = Coordinates[:,1]
    W1 = Weights[:,0]
    W2 = Weights[:,1]

    d = -1
    # if method == 'approx':
    d = s.compareApprox(n, X, Y, W1, W2, L)
    # else:
        # d = s.compareExact(n, X, Y, W1, W2)
    
    sol = {}
    sol['distance'] = d
    sol['runtime'] = s.runtime()
    sol['iterations'] = s.iterations()
    sol['nodes'] = s.num_nodes()
    sol['arcs'] = s.num_arcs()
    sol['status'] = s.status()
    
    return sol


def compareOneToMany(Coordinates, Weights, Options):
    L = Options.get('L', 3)  # L=3 default value
    method = Options.get('Method', 'Approx').encode('utf-8')
    model = Options.get('Model', 'mincostflow').encode('utf-8')
    algo = Options.get('Algorithm', 'colgen').encode('utf-8')
    verbosity = Options.get('Verbosity', 'info').encode('utf-8')
    time_limit = Options.get('TimeLimit', sys.maxsize)
    opt_tolerance = Options.get('OptTolerance', 1e-06)
    
    s = Solver()
    s.setStrParam('Method'.encode('utf-8'), method)
    s.setStrParam('Model'.encode('utf-8'), model)
    s.setStrParam('Algorithm'.encode('utf-8'), algo)
    s.setStrParam('Verbosity'.encode('utf-8'), verbosity)
    s.setDblParam('TimeLimit'.encode('utf-8'), time_limit)
    s.setDblParam('OptTolerance'.encode('utf-8'), opt_tolerance)
    
    n, m = Weights.shape
    m = m - 1
    X = Coordinates[:,0]
    Y = Coordinates[:,1]
    W1 = Weights[:,0]
    Ws = Weights[:,1:]

    d = -1
    # if method == 'approx':
    d = s.compareApprox2(n, m, X, Y, W1, Ws, L)
    # else:
    # d = s.compareApprox2(n, m, X, Y, W1, Ws, n-1)
    
    sol = {}
    sol['distance'] = d
    sol['runtime'] = s.runtime()
    sol['iterations'] = s.iterations()
    sol['nodes'] = s.num_nodes()
    sol['arcs'] = s.num_arcs()
    sol['status'] = s.status()
    
    return sol


def compareAll(Coordinates, Weights, Options):
    L = Options.get('L', 3)  # L=3 default value
    method = Options.get('Method', 'Approx').encode('utf-8')
    model = Options.get('Model', 'mincostflow').encode('utf-8')
    algo = Options.get('Algorithm', 'colgen').encode('utf-8')
    verbosity = Options.get('Verbosity', 'info').encode('utf-8')
    time_limit = Options.get('TimeLimit', sys.maxsize)
    opt_tolerance = Options.get('OptTolerance', 1e-06)
    
    s = Solver()
    s.setStrParam('Method'.encode('utf-8'), method)
    s.setStrParam('Model'.encode('utf-8'), model)
    s.setStrParam('Algorithm'.encode('utf-8'), algo)
    s.setStrParam('Verbosity'.encode('utf-8'), verbosity)
    s.setDblParam('TimeLimit'.encode('utf-8'), time_limit)
    s.setDblParam('OptTolerance'.encode('utf-8'), opt_tolerance)
    
    n, m = Weights.shape
    X = Coordinates[:,0]
    Y = Coordinates[:,1]
    Ws = Weights

    d = -1
    # if method == 'approx':
    d = s.compareApprox3(n, m, X, Y, Ws, L)
    # else:
    # d = s.compareApprox2(n, m, X, Y, W1, Ws, n-1)
    
    sol = {}
    sol['distance'] = d
    sol['runtime'] = s.runtime()
    sol['iterations'] = s.iterations()
    sol['nodes'] = s.num_nodes()
    sol['arcs'] = s.num_arcs()
    sol['status'] = s.status()
    
    return sol
