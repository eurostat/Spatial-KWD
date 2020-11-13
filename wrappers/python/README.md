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
d(a,c) = 3.0236892706218255
d(b,c) = 1.8856180831641267
```

The test program is the following

```
from KWD import *

# Define first histogram
a = Histogram2D()
a.add(0, 0, 1)
a.normalize()

# Define second histogram
b = Histogram2D()
b.add(1, 1, 1)
b.normalize()

# Define third histogram
c = Histogram2D()
c.add(1, 0, 1)
c.add(0, 1, 1)
c.add(5, 5, 1)
c.normalize()

# Define a solver and compute the distance
s = Solver()

print("d(a,b) =", s.distance(a, b, 3))
print("d(a,c) =", s.distance(a, c, 3))
print("d(b,c) =", s.distance(b, c, 3))
```

