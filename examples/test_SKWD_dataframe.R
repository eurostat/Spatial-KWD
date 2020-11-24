# @fileoverview Copyright (c) 2019-2020, Stefano Gualandi,
#               via Ferrata, 1, I-27100, Pavia, Italy
#
# @author stefano.gualandi@gmail.com (Stefano Gualandi)

library(SpatialKWD)

# Random coordinates
N = 900
x <- as.integer(runif(N, 0, 31))
y <- as.integer(runif(N, 0, 31))

# Random weights
h1 <- runif(N, 0, 1)
h2 <- runif(N, 0, 1)

# Data frame
df <- data.frame(x, y, h1, h2)

# Compute distance
d <- distanceDF(df, 3)

print(d)