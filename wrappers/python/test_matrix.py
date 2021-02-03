# @fileoverview Copyright (c) 2019-2020, Stefano Gualandi,
#               via Ferrata, 1, I-27100, Pavia, Italy
#
# @author stefano.gualandi@gmail.com (Stefano Gualandi)

# Library for reading input files
import numpy as np

from KWD import compareOneToOne, compareOneToMany, compareAll

np.random.seed(13)

N = 32 * 32
M = 3

# Random data
Coordinates = np.random.randint(0, 32, size=(N, 2), dtype=np.int32)
Weights = np.random.uniform(0, 100, size=(N, 2))

# Testing helper functions
print('-----------------------------\nTest one2one approx:')
Options = {}
Options['Verbosity'] = 'debug'
Options['Recode'] = 'True'
sol = compareOneToOne(Coordinates, Weights, Options)
for k in sol:
    print(k, sol[k])
print()

print('-----------------------------\nTest one2one exact:')
#Options = {'Method': 'exact'}
sol = compareOneToOne(Coordinates, Weights, Options)
for k in sol:
    print(k, sol[k])
print()

print('-----------------------------\nTest one2many approx:')
Weights = np.random.uniform(0, 100, size=(N, M))
#Options = {}
sol = compareOneToMany(Coordinates, Weights, Options)
for k in sol:
    print(k, sol[k])
print()

print('-----------------------------\nTest all approx:')
Weights = np.random.uniform(0, 100, size=(N, M))
#Options = {}
sol = compareAll(Coordinates, Weights, Options)
for k in sol:
    print(k, sol[k])
