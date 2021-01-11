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
opt = data.frame(Method="KWD_EXACT")
start_time <- Sys.time()
print(compareOneToOne(one2one, Options=opt))
print(c("R total time:", Sys.time()-start_time))

print("Compare one-to-one with approximate algorithm:")
opt = data.frame(Method="KWD_APPROX", L=2)
start_time <- Sys.time()
print(compareOneToOne(one2one, Options=opt))
print(c("R total time:", Sys.time()-start_time))

opt = data.frame(Method="KWD_APPROX", L=3)
print(compareOneToOne(one2one, Options=opt))
print(c("R total time:", Sys.time()-start_time))

opt = data.frame(Method="KWD_APPROX", L=10)
print(compareOneToOne(one2one, Options=opt))
print(c("R total time:", Sys.time()-start_time))

print("Compare one-to-many with approximate algorithm:")
opt = data.frame(Method="KWD_APPROX")
start_time <- Sys.time()
d <- compareOneToMany(one2many, Ws, Options=opt)
print(d)
print(c("R total time:", Sys.time()-start_time))
