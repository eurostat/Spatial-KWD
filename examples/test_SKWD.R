# @fileoverview Copyright (c) 2019-2020, Stefano Gualandi,
#               via Ferrata, 1, I-27100, Pavia, Italy
#
# @author stefano.gualandi@gmail.com (Stefano Gualandi)

library(SpatialKWD)

# Define first histogram
a <- new(SpatialKWD::Histogram2D)
# Add at position (0,0) a unit of mass
a$add(0, 0, 1)
# Normalize the histogram
a$normalize()

# Define second histogram
b <- new(SpatialKWD::Histogram2D)
# Add at position (1,1) two units of mass
b$add(1, 1, 2)
# Normalize the histogram
b$normalize()

# Define third histogram
c <- new(SpatialKWD::Histogram2D)
# Add at position (1,0) and (0,1) an half unit of mass
c$add(1, 0, 0.5)
c$add(0, 1, 0.5)
# Add at position (5,5) a unit of mass
c$add(5, 5, 1)
# Normalize the histogram
c$normalize()

# Define a solver and compute the distance 
# Kantorovich-Wasserstein distance of order 1
# with L_2 as ground distance
s <- new(SpatialKWD::Solver)

print(sprintf("d(a,b) = %f", s$distance(a, b, 3)))
print(sprintf("d(a,c) = %f", s$distance(a, c, 3)))
print(sprintf("d(b,c) = %f", s$distance(b, c, 3)))