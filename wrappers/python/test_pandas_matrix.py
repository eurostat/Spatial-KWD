# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 21:14:47 2021

@author: gualandi
"""

import pandas as pd
import numpy as np

from KWD import compareOneToOne, compareOneToMany, compareAll


N = 100

# Suppose you have your input data as a DataFrame:
Coordinates = pd.DataFrame(np.random.randint(0,32,size=(N, 2), dtype=np.int32), 
                           columns=list('XY'))

Weights = pd.DataFrame(np.random.uniform(0,100,size=(N, 2)), 
                           columns=['W1', 'W2'])

# Then you can "view it" (without copying), as "pointer" to numpy vectors/matrices
# Testing helper functions
print('-----------------------------\nTest one2one approx:')
Options = {}
sol = compareOneToOne(Coordinates.to_numpy(), Weights.to_numpy(), Options)
for k in sol:
    print(k, sol[k])
print()

print('-----------------------------\nTest one2many approx:')
M = 4
Weights = pd.DataFrame(np.random.uniform(0, 100, size=(N, M)))
Options = {}
sol = compareOneToMany(Coordinates.to_numpy(), Weights.to_numpy(), Options)
for k in sol:
    print(k, sol[k])
print()

print('-----------------------------\nTest all approx:')
Weights = pd.DataFrame(np.random.uniform(0, 100, size=(N, M)))
Options = {}
sol = compareAll(Coordinates.to_numpy(), Weights.to_numpy(), Options)
for k in sol:
    print(k, sol[k])
    
