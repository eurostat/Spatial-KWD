\name{CompareOneToOne-function}
\Rdversion{1.1}
\alias{compareOneToOne}
\docType{methods}
\title{
Compare a pair of spatial histograms
}
\description{
This function computes the Kantorovich-Wasserstein between a pair of spatial histograms defined over the same grid map.

The grid map is described by the two lists of \code{N} coordinates \code{Xs} and \code{Ys}, which specify the coordinates of the centroid of each tile of the map.
For each tile \code{i} with coordinates \code{Xs[i], Ys[i]}, we have the two lists of weights, one for the first histograms, and the other for the second histogram.

The two lists of coordinates are passed to \code{compareOneToOne} as a matrix with \code{N} rows and two columns.
The two lists of weights are passed as a matrix with \code{N} rows and two columns, a column for each histogram.
}
\usage{
compareOneToOne(Coordinates, Weights, L = 3, recode = TRUE,
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
      \item{\code{Weights[,1]}: }{\emph{(First Column)} Weights of the first spatial histogram, a weight for each tile located at position \code{Xs[i], Ys[i]} for \code{i=1,...N}. Data type: vector of positive doubles.}
      \item{\code{Weights[,2]}: }{\emph{(Second Column)} Weights of the second spatial histogram, a weight for each tile located at position \code{Xs[i], Ys[i]} for \code{i=1,...N}. Data type: vector of positive doubles.}
    }
  }

  \item{L}{Approximation parameter.
    Higher values of \emph{L} gives more accurate solution, but requires longer running time. Data type: positive integer.}

  \item{recode}{If equal to \code{True}, recode the input coordinates as consecutive integers.}

  \item{method}{Method for computing the KW distances: \code{exact} or \code{approx}.}

  \item{algorithm}{Algorithm for computing the KW distances: \code{fullmodel} or \code{colgen}.}

  \item{model}{Model for building the underlying network: \code{bipartite} or \code{mincostflow}.}

  \item{verbosity}{Level of verbosity of the log: \code{silent}, \code{info} or \code{debug}.}

  \item{timelimit}{Time limit in second for running the solver.}

  \item{opt_tolerance}{Numerical tolerance on the negative reduce cost for the optimal solution.}
}

\details{
The function \code{compareOneToOne(Coordinates, Weights, ...)} computes the distance between the two histograms specified by the weights given in the two columns of matrix \code{Weights}, where the support points (i.e., centroids of each tile of the map) are defined by the coordinates given in \code{Xs} and \code{Ys} in the two columns of matrix \code{Coordinates}. The algorithm used to compute such distance depends on the parameters specified as optional arguments of the function.

The most important is the parameter \code{L}, which by default is equal to 3. The next table shows the worst-case approximation ratio as a function of the value assigned to \code{L}.
The table reports also the number of arcs in the network flow model as a function of the number of bins \emph{n} contained in the convex hull of the support points of the histograms given in input with matrix \code{Coordinates}.
    \tabular{lllllll}{
    \bold{L} \tab \bold{1} \tab \bold{2} \tab \bold{3} \tab \bold{5} \tab \bold{10}\tab \bold{15} \cr
    \code{Worst-case error}  \tab 7.61\% \tab  2.68\% \tab  1.29\% \tab 0.49\%  \tab 0.12\%  \tab   0.06\%  \cr
    \code{Number of arcs} \tab \emph{O(8n)} \tab \emph{O(16n)} \tab \emph{O(32n)}  \tab \emph{O(80n)} \tab \emph{O(256n)} \tab \emph{O(576n)} \cr
    }

  The next two figures show the network build on a grid with 8x8 nodes and using \emph{L=2} and \emph{L=3}.

  \if{html}{\figure{figL2.png}{options: width="35\%" alt="L=2"}}

  \if{html}{\figure{figL3.png}{options: width="35\%" alt="L=3"}}

  \if{latex}{\figure{figL2.png}{options: width=7cm}}

  \if{latex}{\figure{figL3.png}{options: width=7cm}}
}
\value{
    Return an R List with the following named attributes:
  \itemize{
  \item{\code{distance}: }{The value of the KW-distance between the two input histograms.}
  \item{\code{status}: }{Status of the solver used to compute the distances.}
  \item{\code{runtime}: }{Overall runtime in seconds to compute all the distances.}
  \item{\code{iterations}: }{Overall number of iterations of the Network Simplex algorithm.}
  \item{\code{nodes}: }{Number of nodes in the network model used to compute the distances.}
  \item{\code{arcs}: }{Number of arcs in the network model used to compute the distances.}
  }
}
\seealso{
See also \code{\link{compareOneToMany}}, \code{\link{compareAll}}, \code{\link{Histogram2D}}, and \code{\link{Solver}}.
}
\examples{
# Define a simple example
library(SpatialKWD)

# Random coordinates
N = 90
Xs <- as.integer(runif(N, 0, 31))
Ys <- as.integer(runif(N, 0, 31))
coordinates <- matrix(c(Xs, Ys), ncol=2, nrow=N)

# Random weights
test1 <- matrix(runif(2*N, 0, 1), ncol=2, nrow=N)

# Compute distance
print("Compare one-to-one with exact algorithm:")
d <- compareOneToOne(coordinates, Weights=test1, method="exact",
                     recode=TRUE, verbosity = "info")
cat("runtime:", d$runtime, " distance:", d$distance,
    " nodes:", d$nodes, " arcs:", d$arcs, "\n")

print("Compare one-to-one with approximate algorithm:")
d <- compareOneToOne(coordinates, Weights=test1, L=2, recode=TRUE)
cat("L: 2, runtime:", d$runtime, " distance:", d$distance,
    " nodes:", d$nodes, " arcs:", d$arcs, "\n")

d <- compareOneToOne(coordinates, Weights=test1, L=3)
cat("L: 3 runtime:", d$runtime, " distance:", d$distance, "\n")

d <- compareOneToOne(coordinates, Weights=test1, L=10)
cat("L: 10, runtime:", d$runtime, " distance:", d$distance, "\n")
}
