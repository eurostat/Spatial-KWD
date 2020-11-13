# @fileoverview Copyright (c) 2019-2020, Stefano Gualandi,
#               via Ferrata, 1, I-27100, Pavia, Italy
#
# @author stefano.gualandi@gmail.com (Stefano Gualandi)


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
