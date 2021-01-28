library(SpatialKWD)

# @fileoverview Copyright (c) 2020, Stefano Gualandi,
#               via Ferrata, 1, I-27100, Pavia, Italy
#
# @author stefano.gualandi@gmail.com (Stefano Gualandi)

library(SpatialKWD)

# Random coordinates
N = 900
Xs <- 5*as.integer(runif(N, 0, 31))
Ys <- 5*as.integer(runif(N, 0, 31))
coordinates <- matrix(c(Xs, Ys), ncol=2)

# Random weights
test1 <- matrix(runif(2*N, 0, 1), ncol=2)
m <- 3
test2 <- matrix(runif((m+1)*N, 0, 1), ncol=(m+1))
test3 <- matrix(runif(m*N, 0, 1), ncol=m)

# Compute distance
print("Compare one-to-one with exact algorithm:")
d <- compareOneToOne(coordinates, Weights=test1, method="exact", recode=TRUE, verbosity = "debug")
cat("runtime:", d$runtime, " distance:", d$distance, " nodes:", d$nodes, " arcs:", d$arcs, "\n")

print("Compare one-to-one with approximate algorithm:")
d <- compareOneToOne(coordinates, Weights=test1, L=2, recode=TRUE)
cat("L: 2, runtime:", d$runtime, " distance:", d$distance, " nodes:", d$nodes, " arcs:", d$arcs, "\n")

d <- compareOneToOne(coordinates, Weights=test1, L=3)
cat("L: 3 runtime:", d$runtime, " distance:", d$distance, "\n")

d <- compareOneToOne(coordinates, Weights=test1, L=10)
cat("L: 10, runtime:", d$runtime, " distance:", d$distance, "\n")

print("Compare one-to-many with approximate algorithm:")
d <- compareOneToMany(coordinates, Weights=test2)
cat("L: 3, runtime:", d$runtime, " distances:", d$distance, "\n")

print("Compare all with approximate algorithm:")
d <- compareAll(coordinates, Weights=test3, verbosity="info")
cat("L: 3, runtime:", d$runtime, " distances:", "\n")
m <- matrix(d$distance, ncol=3, nrow=3)
print(m)