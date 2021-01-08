\name{SpatialKWD-package}
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
This library contains two main classes and a number of helper functions.

The helper functions are built on top of two main classes: \code{Histogram2D} and \code{Solver}. If you need to specify your algorithms, you can also use directly those two classes.

% COMMENTS LEFT TO INTEGRATE:
%I would mention already in the description of helper-functions at page 2 that convex-hull thing
%(something like "in case of non-convex maps, the algorithm builds the convex-hull and pads the weights with zeros"
%or some similar sentence), and then skip referencing to convex hull in each bullet item in the "details" section.

%Say somewhere that the weights must be non negative, and clarify what happens if negative weights are given in input:
%error message? Warning?

%Adding  Table X  into the documentation, at least for the first few values of L, is a great a idea !

}
\author{
Stefano Gualandi, stefano.gualandi@gmail.com.

Maintainer: Stefano Gualandi <stefano.gualandi@gmail.com>
}
\references{

}
\keyword{ package }
