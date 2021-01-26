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
All the algorithms are based on an ad-hoc implementation of the Network Simplex algorithm [2].
Each implemented algorithm builds a different network, exploiting the special structure of spatial maps.
}
\details{
This library contains three helper functions and two main classes.

The three helper functions are \code{\link{compareOneToOne}}, \code{\link{compareOneToMany}}, and \code{\link{compareAll}}. All the functions take in input the data and an options list. Using the options is possible to configure the Kantorivich-Wasserstein solver, so that it uses different algorithms with different parameters.

The helper functions are built on top of two main classes: \code{\link{Histogram2D}} and \code{\link{Solver}}.

Note that in case of non-convex maps, the algorithm builds the convex-hull of the input bins and pads the weights with zeros.

In case of spatial histograms with weights which do not sum up to 1, all the weights are automatically rescaled in such a way the overall sum of the weight of each single histogram is equal to 1.0: \eqn{w_i \gets \frac{w_i}{\sum_{i=1,\dots,n} w_i}}.
This way the spatial histograms have a natural interpretation as discrete probability measures.
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
coordinates <- matrix(c(Xs, Ys), ncol=2)

# Random weights
weights <- matrix(runif(2*N, 0, 1), ncol=2)

# Compute distance
print("Compare one-to-one with exact algorithm:")
d <- compareOneToOne(coordinates, weights, L=3)
cat("runtime:", d$runtime, " distance:", d$distance, "\n")
}
}
\author{
Stefano Gualandi, stefano.gualandi@gmail.com.

Maintainer: Stefano Gualandi <stefano.gualandi@gmail.com>
}
\references{
[1] Bassetti, F., Gualandi, S. and Veneroni, M., 2020. \emph{On the Computation of Kantorovich--Wasserstein Distances Between Two-Dimensional Histograms by Uncapacitated Minimum Cost Flows}. SIAM Journal on Optimization, 30(3), pp.2441-2469.

[2] Cunningham, W.H., 1976. \emph{A Network Simplex method}. Mathematical Programming, 11(1), pp.105-116.

[3] \url{https://github.com/eurostat/Spatial-KWD}
}
\keyword{ package }
