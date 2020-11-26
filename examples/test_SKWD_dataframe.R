# @fileoverview Copyright (c) 2020, Stefano Gualandi,
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

# Data frame with four columns named: "x", "y", "h1", "h2"
df <- data.frame(x, y, h1, h2)

# Compute distance
print("Test with data frame:")
d <- distanceDF(df, 3)
print(d)


# Data table
require(data.table)
# The table needs to have at least the four columns: "x", "y", "h1", "h2"
dt <- data.table(x=x, y=y, h1=h1, h2=h2)

# Compute distance
print("Test with data table:")
d <- distanceDF(dt, 3)
print(d)


# Increasing the second value, it increase the precision of 
# the computed distance, but it can slow down the runtime speed
print("Benchmarking precision versus speed:")
start_time <- Sys.time()
d <- distanceDF(dt, 4)
print(c(d, Sys.time()-start_time))

start_time <- Sys.time()
d <- distanceDF(dt, 5)
print(c(d, Sys.time()-start_time))

start_time <- Sys.time()
d <- distanceDF(dt, 20)
print(c(d, Sys.time()-start_time))

start_time <- Sys.time()
d <- distanceDF(dt, 30)
print(c(d, Sys.time()-start_time))