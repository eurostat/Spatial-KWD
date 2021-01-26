\name{CompareOneToMany-function}
\Rdversion{1.1}
\alias{compareOneToMany}
\docType{methods}
\title{
Compare a reference spatial histogram to other histograms
}
\description{
This function computes the Kantorovich-Wasserstein among a single reference histogram and a given list of other spatial histograms. All the histograms are defined over the same grid map.

The grid map is described by the two lists of \code{N} coordinates \code{Xs} and \code{Ys}, which specify the coordinates of the centroid of each tile of the map.
For each tile \code{i} with coordinates \code{Xs[i], Ys[i]}, we have a positive weight for each histogram.

The two lists of coordinates are passed to \code{compareOneToMany} as a matrix with \code{N} rows and two columns.
The weights of the histograms are passed as a single matrix with \code{N} rows and \code{M} columns, where the first column is the reference histogram.
}
\usage{
compareOneToMany(Coordinates, Weights, L = 3,
           method = "approx",    algorithm = "colgen",
           model="mincostflow",  verbosity = "silent",
           timelimit = 14400,    opt_tolerance = 1e-06)
}
\arguments{
  \item{Coordinates}{A \code{Matrix} with \code{N} rows and two columns:
    \itemize{
      \item{\code{Coordinates[,1]}: }{\emph{(First Column)} Vector of horizontal coordinates of the centroids of each tile of the map (\code{Xs}). Data type: vector of positive integers.}
      \item{\code{Coordinates[,2]}: }{\emph{(Second Column)} Vector of vertical coordinates of the centroids of each tile of the map (\code{Ys}). Data type: vector of positive integers.}
    }
  }

  \item{Weights}{A \code{Matrix} of positive weights of the tiles specified by \code{Coordinates}.
    \itemize{
      \item{\code{Weights[,1]}: }{\emph{(First Column)} Weights of the reference spatial histogram, a weight for each tile located at position \code{Xs[i], Ys[i]} for \code{i=1,...N}. Data type: vector of positive doubles.}
      \item{\code{Weights[,2:M]}: }{\emph{(Remaining Columns)} Weights of the spatial histograms to be compared with the reference histogram. Data type: vector of positive doubles.}
    }
  }

  \item{L}{Approximation parameter.
    Higher values of \emph{L} gives more accurate solution, but requires longer running time. Data type: positive integer.}

  \item{method}{Method for computing the KW distances: \code{exact} or \code{approx}.}

  \item{algorithm}{Algorithm for computing the KW distances: \code{fullmodel} or \code{colgen}.}

  \item{model}{Model for building the underlying network: \code{bipartite} or \code{mincostflow}.}

  \item{verbosity}{Level of verbosity of the log: \code{silent}, \code{info} or \code{debug}.}

  \item{timelimit}{Time limit in second for running the solver.}

  \item{opt_tolerance}{Numerical tolerance on the negative reduce cost for the optimal solution.}
}

\details{
The function \code{compareOneToMany(Coordinates, Weights, ...)} computes the distances among a reference spatial histogram and a given set of other histograms.
All the histograms are specified by the \code{M} columns of matrix \code{Weights}, and where the support points (i.e., centroids of each tile of the map)
are defined by the coordinates given in \code{Xs} and \code{Ys} in the two columns of matrix \code{Coordinates}.
The algorithm used to compute such distance depends on the parameters specified as optional arguments of the function.

The most important is the parameter \code{L}, which by default is equal to 3 (see \code{\link{compareOneToOne}}).
}
\value{
    Return an R List with the following named attributes:
  \itemize{
  \item{\code{distances}: }{An array of \code{M-1} KW-distances among the input histograms.}
  \item{\code{status}: }{Status of the solver used to compute the distances.}
  \item{\code{runtime}: }{Overall runtime in seconds to compute all the distances.}
  \item{\code{iterations}: }{Overall number of iterations of the Network Simplex algorithm.}
  \item{\code{nodes}: }{Number of nodes in the network model used to compute the distances.}
  \item{\code{arcs}: }{Number of arcs in the network model used to compute the distances.}
  }
}
\seealso{
See also \code{\link{compareOneToOne}}, \code{\link{compareAll}}, \code{\link{Histogram2D}}, and \code{\link{Solver}}.
}
\examples{
  \dontrun{
# Define a simple example
library(SpatialKWD)

# Random coordinates
N = 900
Xs <- as.integer(runif(N, 0, 31))
Ys <- as.integer(runif(N, 0, 31))
coordinates <- matrix(c(Xs, Ys), ncol=2, nrow=N)

# Random weights
m <- 3
test2 <- matrix(runif((m+1)*N, 0, 1), ncol=(m+1))

# Compute distance
print("Compare one-to-many with approximate algorithm:")
d <- compareOneToMany(coordinates, Weights=test2, L=3, method="approx")
cat("L: 3, runtime:", d$runtime, " distances:", d$distance, "\n")
}
}
