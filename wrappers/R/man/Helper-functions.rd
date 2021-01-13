\name{Helper-functions}
\Rdversion{1.1}
\alias{Helper-functions}
\alias{compareOneToOne}
\alias{compareOneToMany}
\alias{compareAll}
\alias{distanceDF}
\docType{methods}
\title{
  Spatial-KWD: Helper functions
}
\description{
The Spatial KWD package has a set of helper functions to ease the use of the library, without explicitly dealing with the internal data structures of the library specified later.
}
\usage{
compareOneToOne(Coordinates, Weights, Options)
compareOneToMany(Coordinates, Weights, Options)
compareAll(Coordinates, Weights, Options)
}
\arguments{
\item{Coordinates}{A \code{Matrix} object that contains two columns:
\itemize{
  \item{\code{Column 0 -> Xs}: }{Vector of horizontal coordinates the bins. Data type: vector of positive integers.}
  \item{\code{Column 1 -> Ys}: }{Vector of vertical coordinates the bins. Data type: vector of positive integers.}
  }
}

\item{Weights}{A \code{Matrix} of positive weights of the bins at the positions specified by \code{Xs} and \code{Ys}. Depending on the function, the \code{Weights} matrix is used as follows.
\itemize{
  \item{\code{Column 0 -> W1}: }{First spatial histogram, given as a vector of weights of the bins located at the positions specified by \code{Xs} and \code{Ys}. Data type: vector of positive doubles.}
  \item{\code{Column 1 -> W2}: }{First spatial histogram, given as a vector of weights of the bins located at the positions specified by \code{Xs} and \code{Ys}. Data type: vector of positive doubles.}
  \item{\code{Column 0 -> Ws}: }{Matrix of weights of the bins located at the positions specified by \code{Xs} and \code{Ys}. Data type: vector of positive doubles.}
}
}

\item{Options}{The main options are described below (see \code{\link{Solver}}).
  \itemize{
  \item{\code{L}: }{Approximation parameter.
    Higher values of \emph{L} gives more accurate solution, but requires longer running time. Data type: positive integer.}
  \item{\code{method}: }{Method for computing the KW distances: \code{exact} or \code{approx}.}
  \item{\code{model}: }{Model for building the underlying network: \code{bipartite} or \code{mincostflow}.}
  \item{\code{algorithm}: }{Algorithm for computing the KW distances: \code{fullmodel} or \code{colgen}.}
  \item{\code{verbosity}: }{Level of verbosity of the log: \code{silent}, \code{info} or \code{debug}.}
  \item{\code{timelimit}: }{Time limit in second for running the solver.}
  \item{\code{opt_tolerance}: }{Numerical tolerance on the negative reduce cost for the optimal solution.}
  }
}
}
\details{
The three functions behave as follows:
\itemize{
  \item{\code{compareOneToOne(Coordinates, Weights, Options)}: }{Compute the distance between the two histograms specified by the weights given in \code{W1} and \code{W2} contained in the matrix \code{Weights}, where the support points are defined by the coordinates given in \code{Xs} and \code{Ys} in the matrix \code{Coordinates}.}
  \item{\code{compareOneToMany(Coordinates, Weights, Options)}: }{Compute the distance among the reference histogram specified by the weights given in the first column of matrix \code{Weights} and all the histograms given as the remaining columns of the matrix \code{Weights} (from the second to the last column). Again the support points of each bin are defined by the two columns \code{Xs} and \code{Ys} in \code{Coordinates}.}
  \item{\code{compareAll(Coordinates, Weights, Options)}: }{Compute the distance among all the columns specified in matrix \code{Ws}, whose support points are in the columns \code{Xs} and \code{Ys} in \code{Coordinates}.}
  }

The next table shows the worst-case approximation ratio as a function of the parameter \code{L}. The table reports also the number of arcs in the network flow model as a function of the number of bins \emph{n} contained in the convex hull of the support points of the histograms given in input with the matrix.
    \tabular{lllllll}{
    \bold{L} \tab \bold{1} \tab \bold{2} \tab \bold{3} \tab \bold{5} \tab \bold{10}\tab \bold{15} \cr
    \code{Worst-case errror}  \tab 7.61\% \tab  2.68\% \tab  1.29\% \tab 0.49\%  \tab 0.12\%  \tab   0.06\%  \cr
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
  \item{\code{distances}: }{Vector of values, with the KW-distances betweeen the input histograms.}
  \item{\code{status}: }{Status of the solver used to compute the distances.}
  \item{\code{runtime}: }{Overall runtime in seconds to compute all the distances.}
  \item{\code{iterations}: }{Overall number of iterations of the Network Simplex algorithm.}
  \item{\code{nodes}: }{Number of nodes in the network model used to compute the distances.}
  \item{\code{arcs}:}{Number of arcs in the network model used to compute the distances.}
  }
}
\seealso{
See also \code{\link{Histogram2D}} and \code{\link{Solver}}.
}
\examples{
  \dontrun{
# Define a simple example
library(SpatialKWD)

# Random coordinates
N = 900
Xs <- as.integer(runif(N, 0, 31))
Ys <- as.integer(runif(N, 0, 31))
coordinates <- matrix(c(Xs, Ys), ncol=2)

# Random weights
test1 <- matrix(runif(2*N, 0, 1), ncol=2)
m <- 3
test2 <- matrix(runif((m+1)*N, 0, 1), ncol=(m+1))
test3 <- matrix(runif(m*N, 0, 1), ncol=m)

# Compute distance
print("Compare one-to-one with exact algorithm:")
opt = data.frame(Method="Exact")
d <- compareOneToOne(coordinates, Weights=test1, Options=opt)
cat("runtime:", d$runtime, " distance:", d$distance, "\n")

print("Compare one-to-one with approximate algorithm:")
opt = data.frame(Method="Approx", L=2)
d <- compareOneToOne(coordinates, Weights=test1, Options=opt)
cat("L: 2, runtime:", d$runtime, " distance:", d$distance, "\n")

opt = data.frame(Method="Approx", L=3)
d <- compareOneToOne(coordinates, Weights=test1, Options=opt)
cat("L: 3 runtime:", d$runtime, " distance:", d$distance, "\n")

opt = data.frame(Method="Approx", L=10)
d <- compareOneToOne(coordinates, Weights=test1, Options=opt)
cat("L: 10, runtime:", d$runtime, " distance:", d$distance, "\n")

print("Compare one-to-many with approximate algorithm:")
opt = data.frame(Method="Approx")
d <- compareOneToMany(coordinates, Weights=test2, Options=opt)
cat("L: 3, runtime:", d$runtime, " distances:", d$distance, "\n")

print("Compare all with approximate algorithm:")
opt = data.frame(Method="Approx")
d <- compareAll(coordinates, Weights=test3, Options=opt)
cat("L: 3, runtime:", d$runtime, " distances:", "\n")
m <- matrix(d$distance, ncol=3, nrow=3)
print(m)
}
}
