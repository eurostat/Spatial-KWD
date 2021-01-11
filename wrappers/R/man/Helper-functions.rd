\name{Helper-functions}
\Rdversion{1.1}
\alias{Helper-functions}
\alias{compareOneToOne}
\alias{compareOneToMany}
\alias{distanceDF}
\docType{methods}
\title{
  Spatial-KWD: Helper functions
}
\description{
The Spatial KWD package has a set of helper functions to ease the use of the library, without explicitly dealing with the internal data structures of the library specified later.
}
\usage{
compareOneToOne(Data, Options)
compareOneToMany(Data, Ws, Options)
}
\details{
The next table shows the worst-case approximation ratio as a function of the parameter \code{L}. The table reports also the number of arcs in the network flow model as a function of the number of bins \emph{n} contained in the convex hull of the support points of the histograms given in input with the dataframe.

    \tabular{lllllll}{
    \bold{L} \tab \bold{1} \tab \bold{2} \tab \bold{3} \tab \bold{5} \tab \bold{10}\tab \bold{15} \cr
    \code{Worst-case errror}  \tab 7.61\% \tab  2.68\% \tab  1.29\% \tab 0.49\%  \tab 0.12\%  \tab   0.06\%  \cr
    \code{Number of arcs} \tab \emph{O(8n)} \tab \emph{O(16n)} \tab \emph{O(32n)}  \tab \emph{O(80n)} \tab \emph{O(256n)} \tab \emph{O(576n)} \cr
    }

}
\arguments{
\item{Data}{It is a \code{Dataframe} object that contains the following fields:
\itemize{
  \item{\code{Xs}: }{Vector of horizontal coordinates the bins. Data type: vector of integers.}
  \item{\code{Ys}: }{Vector of vertical coordinates the bins. Data type: vector of integers.}
  \item{\code{W1}: }{Vector of positive weights of the bins located at the positions specified by \code{Xs} and \code{Ys}. Data type: vector of doubles.}
  }
}

\item{Ws}{Matrix of weights of the bin at the positions specified by \code{Xs} and \code{Ys}. Data type: matrix of doubles.}

\item{Options}{The main options are described below:
  \itemize{
  \item{\code{L}: }{Approximation parameter.
    Higher values of \emph{L} gives more accurate solution, but requires longer running time.
    The table below gives the percentage worst-case approximation error as a function of \emph{L}. Data type: positive integer}
  \item{\code{Method}: }{Method for KW distance algorithm: \code{KWD_EXACT} or \code{KWD_APPROX}.}
  }
}
}
\value{
    Return a dataframe with the following attributes:
  \itemize{
  \item{\code{distance}: }{Vector of values, with the KW-distances betweeen the input histograms.}
  \item{\code{status}: }{Status of the solver used to compute the distances.}
  \item{\code{runtime}: }{Overall runtime in seconds to compute all the distances.}
  \item{\code{iterations}: }{Overall number of iterations of the Network Simplex algorithm.}
  \item{\code{nodes}: }{Number of nodes in the network model used to compute the distances.}
  \item{\code{arcs}:}{Number of arcs in the network model used to compute the distances.}
  }
}
\examples{
  \dontrun{
# Define a simple example
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
}
}
