library(SpatialKWD)

# @fileoverview Copyright (c) 2020, Stefano Gualandi,
#               via Ferrata, 1, I-27100, Pavia, Italy
#
# @author stefano.gualandi@gmail.com (Stefano Gualandi)

library(SpatialKWD)

# Random coordinates
N = 900
Xs <- as.integer(runif(N, 0, 31))
Ys <- as.integer(runif(N, 0, 31))
coordinates <- matrix(c(Xs, Ys), ncol=2)

# Random weights
test1 <- matrix(runif(2*N, 0, 1), ncol=2)
m <- 3
test2 <- matrix(runif((m+1)*N, 0, 1), ncol=(m+1))
test3 <- matrix(runif(m*N, 0, 1), ncol=m)

# Compute distance
print("Compare one-to-one with exact algorithm:")
opt = data.frame(Method="Exact")
d <- compareOneToOne(coordinates, Weights=test1, Options=opt)
cat("runtime:", d$runtime, " distance:", d$distance, "\n")

print("Compare one-to-one with approximate algorithm:")
opt = data.frame(Method="Approx", L=2)
d <- compareOneToOne(coordinates, Weights=test1, Options=opt)
cat("L: 2, runtime:", d$runtime, " distance:", d$distance, "\n")

opt = data.frame(Method="Approx", L=3)
d <- compareOneToOne(coordinates, Weights=test1, Options=opt)
cat("L: 3 runtime:", d$runtime, " distance:", d$distance, "\n")

opt = data.frame(Method="Approx", L=10)
d <- compareOneToOne(coordinates, Weights=test1, Options=opt)
cat("L: 10, runtime:", d$runtime, " distance:", d$distance, "\n")

print("Compare one-to-many with approximate algorithm:")
opt = data.frame(Method="Approx")
d <- compareOneToMany(coordinates, Weights=test2, Options=opt)
cat("L: 3, runtime:", d$runtime, " distances:", d$distance, "\n")

print("Compare all with approximate algorithm:")
opt = data.frame(Method="Approx")
d <- compareAll(coordinates, Weights=test3, Options=opt)
cat("L: 3, runtime:", d$runtime, " distances:", "\n")
m <- matrix(d$distance, ncol=3, nrow=3)
print(m)