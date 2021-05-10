# @fileoverview Copyright (c) 2019-2021, Stefano Gualandi,
#               via Ferrata, 1, I-27100, Pavia, Italy
#
# @author stefano.gualandi@gmail.com (Stefano Gualandi)

# Library for reading input files
import numpy
import sys

from time import perf_counter

from KWD import *

X, Y, A, B = [], [], [], []

# open a file, where you stored the pickled data
with open(sys.argv[1], 'r') as fh:
    # Skip header line
    fh.readline()
    # Parse file
    for row in fh:
        l = row.replace('\n', '').split(',')
        X.append(np.int32(l[0]))
        Y.append(np.int32(l[1]))
        A.append(float(l[2]))
        B.append(float(l[3]))

Coordinates = np.transpose(np.array([X, Y]))
Weights = np.transpose(np.array([A, B]))

Options = {}
Options['Verbosity'] = 'info'
Options['Recode'] = 'True'
Options['L'] = 3

print("TESTING: ComparaOneToOne")
sol = compareOneToOne(Coordinates, Weights, Options)

for k in sol:
    print(k, sol[k])
print()

print("TESTING: focusArea")
x = len(X) // 2
y = len(Y) // 2
radius = 100
Options['Verbosity'] = 'debug'
sol = focusArea(Coordinates, Weights, X[x], Y[y], radius, Options)

for k in sol:
    print(k, sol[k])
print()
