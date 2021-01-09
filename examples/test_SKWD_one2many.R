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

# Random weights
W1 <- runif(N, 0, 1)
W2 <- runif(N, 0, 1)

m <- 3
Ws <- matrix(runif(m*N, 0, 1), ncol=3)

# Data frame with four columns named: "x", "y", "h1", "h2"
one2one <- data.frame(Xs, Ys, W1, W2)
one2many <- data.frame(Xs=Xs, Ys=Ys, W1=W1)


# Compute distance
print("Compare one-to-one with exact algorithm:")
start_time <- Sys.time()
d <- compareExact(one2one)
print(c(d, Sys.time()-start_time))

print("Compare one-to-one with approximate algorithm:")
start_time <- Sys.time()
d <- compareApprox(one2one, 5)
print(c(d, Sys.time()-start_time))

print("Compare one-to-many with approximate algorithm:")
start_time <- Sys.time()
d <- compareOneToManyApprox(one2many, Ws, 3)
print(c(d, Sys.time()-start_time))
