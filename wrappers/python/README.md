# Spatial-KWD: Kantorovich-Wasserstein Distances for Large Spatial Maps

The Spatial-KWD package contains efficient implementations of Discrete Optimal Transport algorithms for the computation of Kantorovich-Wasserstein distances [1], which are customized for large spatial maps.
All the algorithms are based on an ad-hoc implementation of the Network Simplex algorithm [2].
Each implemented algorithm builds a different network, exploiting the special structure of spatial maps.

## Details
This library contains three helper functions and two main classes.

The three helper functions are `compareOneToOne`, `compareOneToMany`, and `compareAll`. All the functions take in input the data and an options list. Using the options is possible to configure the Kantorivich-Wasserstein solver, so that it uses different algorithms with different parameters.

The helper functions are built on top of two main classes: `Histogram2D` and `Solver`.

Note that in case of non-convex maps, the algorithm builds the convex-hull of the input bins and pads the weights with zeros.

## Prerequisities

You only need the following library:

* Cython (>= 0.23): tested on Windows 10 with v0.29.21

## Installation

You have to download this repository, and the from command line run the command:

```bash
python setup.py build_ext --inplace
```

This will compile the C++ code and build the python wrapper.

## Testing

For testing the library you can run the following command:

```bash
python test_matrix.py
```

The test program is the following

```python
import numpy as np

from KWD import compareOneToOne, compareOneToMany, compareAll


np.random.seed(13)

N = 32*32
M = 3

# Random data
Coordinates = np.random.randint(0, 32, size=(N, 2), dtype=np.int32)
Weights = np.random.uniform(0, 100, size=(N, 2))

# Testing helper functions
print('-----------------------------\nTest one2one approx:')
Options = {}
sol = compareOneToOne(Coordinates, Weights, Options)
for k in sol:
    print(k, sol[k])
print()

print('-----------------------------\nTest one2one exact:')
Options = {'Method': 'exact' }
sol = compareOneToOne(Coordinates, Weights, Options)
for k in sol:
    print(k, sol[k])
print()

print('-----------------------------\nTest one2many approx:')
Weights = np.random.uniform(0, 100, size=(N, M))
Options = {}
sol = compareOneToMany(Coordinates, Weights, Options)
for k in sol:
    print(k, sol[k])
print()

print('-----------------------------\nTest all approx:')
Weights = np.random.uniform(0, 100, size=(N, M))
Options = {}
sol = compareAll(Coordinates, Weights, Options)
for k in sol:
    print(k, sol[k])
```

### References
[1] Bassetti, F., Gualandi, S. and Veneroni, M., 2020. \emph{On the Computation of Kantorovich--Wasserstein Distances Between Two-Dimensional Histograms by Uncapacitated Minimum Cost Flows}. SIAM Journal on Optimization, 30(3), pp.2441-2469.

[2] Cunningham, W.H., 1976. \emph{A Network Simplex method}. Mathematical Programming, 11(1), pp.105-116.

[3] https://github.com/eurostat/Spatial-KWD

### Authors and maintainers
Stefano Gualandi, stefano.gualandi@gmail.com.

Maintainer: Stefano Gualandi <stefano.gualandi@gmail.com>
