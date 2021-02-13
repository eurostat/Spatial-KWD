\name{Histogram2D-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_Histogram2D-class}
\alias{Histogram2D}
\alias{Rcpp_Histogram2D}
\title{
  Two Dimensional Histogram for Spatial Data
}
\description{
The \code{Histogram2D} class represents a single spatial 2-dimensional histograms. The class is mainly composed of three vectors of the same length \code{n}. The first two vectors of integers, called \code{Xs} and \code{Ys}, give the coordinates of each bin of the histogram, while the third
  vector of doubles, called \code{Ws}, gives the weight \code{Ws[i]} of the \emph{i}-th bin located at position \code{Xs[i]} and \code{Ys[i]}.

  A 2D histogram can be also defined by adding (or updating) a single element a the time (see second constructor).

  Note that the bins position are not required to lay on rectangular (or squared) grid, but they can lay every where in the plane. Before of computing the distance between a pair of algorithms, the solver will compute a convex hull of all non-empty bins.
}
\details{
  The public methods of the \code{Histogram2D} class are.
}
\seealso{
See also \code{\link{compareOneToOne}}, \code{\link{compareOneToMany}}, \code{\link{compareAll}}, and \code{\link{Solver}}.
}
\examples{
# Define a simple histogram
h <- new(SpatialKWD::Histogram2D)

# Add half unit of mass at positions (1,0) and (0,1)
h$add(1, 0, 0.5)
h$add(0, 1, 0.5)

# Add at position (5,5) a unit of mass
h$add(5, 5, 1)

# Normalize the histogram
h$normalize()

# Print the total weight (mass) of the histogram
print(sprintf("Histogram total weight = \%f", h$balance()))
}
\keyword{classes}
\section{Methods}{
  \describe{
    \item{\code{Histogram2D(n, Xs, Ys, Ws)}:}{ c'tor. }
    \item{\code{add(x, y, w)}:}{it adds a bin located at position \emph{x} and \emph{y} with weight \emph{w}. }
    \item{\code{update(x, y, u)}:}{ return the total mass balance of this histogram, that is, return the quantity \eqn{\sum_{i=1,\dots,n} w_i}. }
    \item{\code{size()}:}{return the number of non-empty bins \emph{n} of this histogram. }
    \item{\code{normalize()}:}{ normalize the weights of all non-empty bins, such that they all sum up to 1. Indeed, this method perform the operation: \eqn{w_i \gets \frac{w_i}{\sum_{i=1,\dots,n} w_i}}. }
    \item{\code{balance()}:}{ return the total mass balance of this histogram, that is, return the quantity \eqn{\sum_{i=1,\dots,n} w_i}.}
  }
}
\arguments{
  \item{n}{Number of non-empty bins. Type: positive integer.}
  \item{Xs}{Vector of horizontal coordinates the bins. Type: vector of integers.}
  \item{Ys}{Vector of vertical coordinates the bins. Type: vector of integers.}
  \item{Ws}{Vector of positive weights of the bin at position \emph{(x,y)}. Type: vector of positive doubles.}
  \item{x}{Horizontal coordinate of a bin. Type: integer.}
  \item{y}{Vertical coordinate of a bin. Type: integer.}
  \item{w}{Weight of the bin at position \emph{(x,y)}. Type: positive double.}
  \item{u}{Weight of the bin to be added to the weight at position \emph{(x,y)}. If a bin in position \emph{(x,y)} is absent, then it is added with weight equal to \code{u}. Type: positive double.}
}
\value{
  The \code{add}, \code{update}, and \code{normalize} does not return any value.

  The \code{size} method returns the number of non-empty bins in \code{h}.

  The \code{balance} method returns the sum of the weights in \code{h}.
}
