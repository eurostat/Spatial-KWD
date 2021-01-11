\name{compareOneToManyApprox}
\alias{compareOneToManyApprox}
\docType{methods}
\title{
  Compare an histogram with a matrix of histograms
}
\description{
This function permits to compare a single spatial histogram with a number of given spatial histograms. 
All histograms must have the same size of the reference spatial histogram.
}
\usage{
compareOneToManyApprox(DF, Mat, L)
}

\arguments{
\item{DF}{Dataframe.}
\item{Mat}{Matrix.}
\item{L}{Integer parameter, approximation factor.}
%  \item{n} Size of the reference spatial histogram.
%  \item{m} Number of histograms to compare with.
%  \item{Xs}{Vector of horizontal coordinates the bins. Type: vector of integers.}
%  \item{Ys}{Vector of vertical coordinates the bins. Type: vector of integers.}
%  \item{W1}{Vector of weights of the bin at the positions specified by \code{Xs} and \code{Ys}. Type: vector of doubles.}
%  \item{Ws}{Matrix of weights of the bin at the positions specified by \code{Xs} and \code{Ys}. Type: matrix of doubles.}
%  \item{Model}{Model to use in the bipartite netowrk (bipartite or uncapacitated network flow).}
%  \item{L}{Approximation parameter. Higher values of \emph{L} gives more accurate solution, but require longer running time.
%  Table X gives the guarantee approximation bound as a function of \emph{L}. Type: positive integer.}
%  \item{Options}{See definition at page X.}
}
\details{
  This function computes a vector \emph{Ds} of \emph{m} distances from the reference spatial histogram \emph{W1}.
  The element \code{Ds(i)} gives the KW-distance between \code{W1} and the vector \code{Ws(i)} corresponding to the \emph{i}-th column of \code{Ws}.
}
\examples{
  \dontrun{
# Define a simple example
}
}
