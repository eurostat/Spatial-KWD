\name{Helper-functions}
\alias{Helper-functions}
\alias{compareExact}
\alias{compareApprox}
\alias{distanceDF}
\docType{methods}
\title{
  Helper functions
}
\description{
The Spatial KWD package has a set of helper functions to ease the use of the library, without explicitly dealing with the internal data structures of the library specified later.
}
\usage{
# Exact methods
compareExact(Xs, Ys, W1, W2)
compareExact(Xs, Ys, W1, Ws)
compareExact(Xs, Ys, Ws)

# Approximation methods (with bound guarantee)
compareApprox(Xs, Ys, W1, W2, L)
compareApprox(Xs, Ys, W1, Ws, L)
compareApprox(Xs, Ys, Ws, L)
}

\arguments{
  \item{Xs}{Vector of horizontal coordinates the bins. Type: vector of integers.}
  \item{Ys}{Vector of vertical coordinates the bins. Type: vector of integers.}
  \item{W1}{Vector of weights of the bin at the positions specified by \code{Xs} and \code{Ys}. Type: vector of doubles.}
  \item{W2}{Vector of weights of the bin at the positions specified by \code{Xs} and \code{Ys}. Type: vector of doubles.}
  \item{Ws}{Matrix of weights of the bin at the positions specified by \code{Xs} and \code{Ys}. Type: matrix of doubles.}
  \item{L}{Approximation parameter. Higher values of \emph{L} gives more accurate solution, but require longer running time.
  Table X gives the guarantee approximation bound as a function of \emph{L}. Type: positive integer.}
}
\details{
    The detailed descriptions of the helper functions is as follows:
  \itemize{
    \item \code{compareExact(Xs, Ys, W1, W2)}: compute the exact distance between the two vector of weights \code{W1} and \code{W2}, on the convex hull of the points defined by the two vectors \code{Xs} and \code{Ys}.
    This method return a single value (double), which is the KW-distance between \code{W1} and \code{W2}.

    \item \code{compareExact(Xs, Ys, W1, Ws)}: compute the exact distances between the vector of weights \code{W1} and each of the vectors of weights in \code{Ws}, on the convex hull of the points defined by the two vectors \code{Xs} and \code{Ys}.
    This method returns a vector of double of the same size of \code{Ws}, representing the distance of \code{W1} to every element of \code{Ws}.

    \item \code{compareExact(Xs, Ys, Ws)}: compute a symmetric matrix of pairwise exact distances between all the possible pairs of the vector listed in \code{Ws}.
    The algorithm used by the solver is controlled by the parameter \code{ExactMethod} (see below).


    \item \code{compareApprox(Xs, Ys, W1, W2, L)}: compute the approximate distance between the two vector of weights \code{W1} and \code{W2}, on the convex hull of the points defined by the two vectors \code{Xs} and \code{Ys}.
    This method return a single value (double), which is the KW-distance between \code{W1} and \code{W2}.

    \item \code{compareApprox(Xs, Ys, W1, Ws, L)}: compute the approximate distances between the vector of weights \code{W1} and each of the vector of weights in \code{Ws}, on the convex hull of the points defined by the two vectors \code{Xs} and \code{Ys}.
    This method return a vector of double of the same size of \code{Ws}, representing the distance of \code{W1} to every element of \code{Ws}.

    \item \code{compareApprox(Xs, Ys, Ws, L)}: compute a symmetric matrix of pairwise approximate distances (which quality depends on the value of \emph{L}) between all the possible pairs of the vector listed in \code{Ws}.
    The algorithm used by the solver is controlled by the parameter \code{ApproxMethod} (see below).
  }
}
\examples{
  \dontrun{
# Define a simple example
}
}
