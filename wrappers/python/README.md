# Basic istructions

### Prerequisities

You only need the following library:

* Cython (>= 0.23): tested on Windows 10 with v0.29.21

### Installation

You have to download this repository, and the from command line run the command:

```
python setup.py build_ext --inplace
```

This will compile the C++ code and build the python wrapper.

### Testing

For testing the library you can run the following command:

```
python test_kwd.py
```

And the output should be:

```
d(a,b) = 1.4142135623730951
d(a,c) = 4.035533905932738
d(b,c) = 2.8284271247461903
```

The test program is the following

```
from KWD import *

# Define first histogram
a = Histogram2D()
# Add at position (0,0) a unit of mass
a.add(0, 0, 1)
# Normalize the histogram
a.normalize()

# Define second histogram
b = Histogram2D()
# Add at position (1,1) two units of mass
b.add(1, 1, 2)
# Normalize the histogram
b.normalize()

# Define third histogram
c = Histogram2D()
# Add at position (1,0) and (0,1) an half unit of mass
c.add(1, 0, 0.5)
c.add(0, 1, 0.5)
# Add at position (5,5) a unit of mass
c.add(5, 5, 1)
# Normalize the histogram
c.normalize()

# Define a solver and compute the distance 
# Kantorovich-Wasserstein distance of order 1
# with L_2 as ground distance
s = Solver()

print("d(a,b) =", s.distance(a, b, 3))
print("d(a,c) =", s.distance(a, c, 3))
print("d(b,c) =", s.distance(b, c, 3))
```

