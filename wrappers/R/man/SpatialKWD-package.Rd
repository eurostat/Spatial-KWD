\name{SpatialKWD-package}
\Rdversion{1.1}
\alias{SpatialKWD-package}
\alias{SpatialKWD}
\docType{package}
\title{
  Kantorovich-Wasserstein Distances for Large Spatial Maps
}
\description{
The Spatial-KWD package contains efficient implementations of Discrete Optimal Transport algorithms for the computation of Kantorovich-Wasserstein distances [1], which are customized for large spatial maps.
All the algorithms are based on a core implementation of the Network Simplex algorithm.
Each implemented algorithm builds a different network, exploiting the special structure of spatial maps.
}
\details{
This library contains three helper functions and two main classes.

The three helper functions are \code{\link{compareOneToOne}}, \code{\link{compareOneToMany}}, and \code{\link{compareAll}}. All the functions take in input the data and an options list. Using the options is possible to configure the Kantorivich-Wasserstein solver, so that it uses different algorithms with different algorithms.

The helper functions are built on top of two main classes: \code{\link{Histogram2D}} and \code{\link{Solver}}. If you need to specify your algorithms, you can also use directly those two classes.

% COMMENTS LEFT TO INTEGRATE:
%I would mention already in the description of helper-functions at page 2 that convex-hull thing
%(something like "in case of non-convex maps, the algorithm builds the convex-hull and pads the weights with zeros"
%or some similar sentence), and then skip referencing to convex hull in each bullet item in the "details" section.

%Say somewhere that the weights must be non negative, and clarify what happens if negative weights are given in input:
%error message? Warning?

%Adding  Table X  into the documentation, at least for the first few values of L, is a great a idea !

}
\seealso{
See also \code{\link{compareOneToOne}}, \code{\link{compareOneToMany}}, \code{\link{compareAll}}, \code{\link{Histogram2D}}, and \code{\link{Solver}}.
}
\examples{
  \dontrun{
library(SpatialKWD)

# Random coordinates
N = 900
Xs <- as.integer(runif(N, 0, 31))
Ys <- as.integer(runif(N, 0, 31))

# Random weights
W1 <- runif(N, 0, 1)
W2 <- runif(N, 0, 1)

# Data frame with four columns named: "Xs", "Ys", "W1", "W2"
one2one <- data.frame(Xs, Ys, W1, W2)

# Compute distance
print("Compare one-to-one with exact algorithm:")
opt = data.frame(Method="KWD_EXACT")
start_time <- Sys.time()
ls <- compareOneToOne(one2one, Options=opt))
print(c("R total time:", Sys.time()-start_time))
}
}
\author{
Stefano Gualandi, stefano.gualandi@gmail.com.

Maintainer: Stefano Gualandi <stefano.gualandi@gmail.com>
}
\references{
    \url{https://github.com/eurostat/Spatial-KWD}
}
\keyword{ package }
