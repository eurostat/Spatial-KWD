\name{FocusArea-function}
\Rdversion{1.1}
\alias{focusArea}
\docType{methods}
\title{
Compute the KWD tranport distance within a given focus area
}
\description{
This function computes the Kantorovich-Wasserstein distance within a given focus area embedded into a large region described as a grid map.
Both the focus and the embedding areas are are described by spatial histograms, similarly to the input data of the other functions of this package.

The grid map is described by the two lists \code{Xs} and \code{Ys} of \code{N} coordinates, which specify the coordinates of the centroid of every single tile.
For each tile \code{i} with coordinates \code{Xs[i], Ys[i]}, we have an entry in the two lists of weights \code{W1} and \code{W2}, one for the first histograms, and the other for the second histogram.

The two lists of coordinates \code{Xs} and \code{Ys} are passed to the \code{focusArea} function as a matrix with \code{N} rows and two columns.
The two lists of weights \code{W1} and \code{W2} are passed as a matrix with \code{N} rows and two columns, a column for each histogram.

The focus area is specified by three parameters: the coordinates \code{x} and \code{y} of the center of the focus area, and the (circular) \code{radius} of the focus area.
The pair of coordinates (\code{x,y}) must correspond to a pair of coordinates contained in the vectors \code{Xs,Ys}.
Every tile whose distance is less or equal to the \code{radius} will be included in the focus area.

The focus area by default is circular, that is, the area is based on a $\ell_2$ norm. By setting the parameter \code{area} to the value \code{linf} it is possible to obtain
a squared focus area, induced by the norm $\ell_\infty$.
}
\usage{
focusArea(Coordinates, Weights, x, y, radius,
           L = 3, recode = TRUE,
           method = "approx",    algorithm = "colgen",
           model="mincostflow",  verbosity = "silent",
           timelimit = 14400,    opt_tolerance = 1e-06,
           area = "l2")
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
      \item{\code{Weights[,1]}: }{\emph{(First Column)} Weights of the embedding spatial histogram, a weight for each tile located at position \code{Xs[i], Ys[i]} for \code{i=1,...N}. Data type: vector of positive doubles.}
      \item{\code{Weights[,2]}: }{\emph{(Second Column)} Weights of the spatial histogram of the focus area, a weight for each tile located at position \code{Xs[i], Ys[i]} for \code{i=1,...N}. All the weights outside the focus area should be equal to zero. Data type: vector of positive doubles.}
    }
  }

  \item{x}{Horizontal coordinate of the centroid of the focus area.}

  \item{y}{Vertical coordinate of the centroid of the focus area.}

  \item{radius}{The radius of the focus area.}

  \item{L}{Approximation parameter.
    Higher values of \emph{L} give a more accurate solution, but they require a longer running time. Data type: positive integer.}

  \item{recode}{If equal to \code{True}, recode the input coordinates as consecutive integers.}

  \item{method}{Method for computing the KW distances: \code{exact} or \code{approx}.}

  \item{algorithm}{Algorithm for computing the KW distances: \code{fullmodel} or \code{colgen}.}

  \item{model}{Model for building the underlying network: \code{bipartite} or \code{mincostflow}.}

  \item{verbosity}{Level of verbosity of the log: \code{silent}, \code{info}, or \code{debug}.}

  \item{timelimit}{Time limit in second for running the solver.}

  \item{opt_tolerance}{Numerical tolerance on the negative reduced cost for the optimal solution.}

  \item{area}{Type of norm for delimiting the focus area: \code{l2} denotes a circular area of radius, \code{linf} denotes a squared area.}
}

\details{
The function \code{focusArea(Coordinates, Weights, x, y, radius, ...)} computes the KW distance within a focus area by implicitly considering the surrounding larger area.
The mass contained within the focus area is transported to a destination either within or outside the focus area.
All the mass contained outside the focus area could be used to balance the mass within the focus area.
}
\value{
    Return an R List with the following named attributes:
  \itemize{
  \item{\code{distance}: }{The value of the KW-distance between the two input areas.}
  \item{\code{status}: }{Status of the solver used to compute the distances.}
  \item{\code{runtime}: }{Overall runtime in seconds to compute all the distances.}
  \item{\code{iterations}: }{Overall number of iterations of the Capacitated Network Simplex algorithm.}
  \item{\code{nodes}: }{Number of nodes in the network model used to compute the distances.}
  \item{\code{arcs}: }{Number of arcs in the network model used to compute the distances.}
  }
}
\seealso{
See also \code{\link{compareOneToOne}}, \code{\link{compareOneToMany}}, \code{\link{compareAll}}, \code{\link{Histogram2D}}, and \code{\link{Solver}}.
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
d <- focusArea(coordinates, Weights=test1,
                x=15, y=15, radius=5,
                method="exact", recode=TRUE, verbosity = "info")
cat("runtime:", d$runtime, " distance:", d$distance,
    " nodes:", d$nodes, " arcs:", d$arcs, "\n")
}
