  \name{Solver-class}
\Rdversion{1.1}
\docType{class}
\alias{Solver}
\alias{Rcpp_Solver-class}
\title{
  Spatial-KWD Solver
}
\description{
The \code{Solver} class is the main wrapper to the core algorithms implemented in the Spatial KWD package.
It has several methods that permit to compare two, or more, objects of type \code{Histogram2D}.
If you use the helper functions described at the begging of this document, you can avoid using this class directly
}
\details{
    The public methods of this class are:

  The \code{Solver} class can be controlled by the list of parameters given in the following table, which can be set with the \code{setParam(name, value)} method. A detailed description of each parameter is given below.

    \tabular{lll}{
    \bold{Parameter Name} \tab \bold{Possible Values} \tab \bold{Default Value} \cr
    \code{Method}  \tab \code{exact, approx} \tab \code{approx} \cr
    \code{Model}  \tab \code{bipartite, mincostflow} \tab \code{mincostflow} \cr
    \code{Algorithm} \tab \code{fullmodel, colgen} \tab \code{colgen}\cr
    \code{Verbosity}    \tab \code{silent, info, debug} \tab \code{info} \cr
    \code{TimeLimit}    \tab Any positive integer smaller than \code{INTMAX} \tab \code{INTMAX} \cr
    \code{OptTolerance} \tab Any value in \eqn{[10^{-9}, 10^{-1}]} \tab \eqn{10^{-6}}
    }

  \itemize{
    \item \code{Method}: set which method to use for computing the exact distance between a pair of histograms. The options for this parameter are:
      \itemize{
      \item \code{exact}: Compute the exact KW distance. This method is only helpful for small and sparse spatial maps.

      \item \code{approx}: Compute an approximation KW distance which depends on the parameter \emph{L}. This is the default value.
      }

\item \code{Model}: set which network model to use for computing the exact distance between a pair of histograms. The options for this parameter are:
      \itemize{
      \item \code{bipartite}: Build a complete bipartite graph. This method is only helpful for small and sparse spatial maps.

      \item \code{mincostflow}: Build an uncapacitated network flow. This is, in general, smaller than the \code{bipartite} model, except for very sparse histograms.
      }

    \item \code{Algorithm}: set which algorithm to use to compute an approximate distance between a pair of histograms, which depends on the parameter \emph{L}. The options for this parameter are:
      \itemize{
      \item \code{fullmodel}: Build a complete network model and solve the corresponding problem.

      \item \code{colgen}: Build the network model incrementally while computing the KW distance. It is the recommended method for very large dense spatial maps. On medium and small spatial maps, the fullmodel could be faster.
    }
      The default value is set to \code{colgen}.

    \item \code{Verbosity}: set the level of verbosity of the logs. Possible values are \code{silent}, \code{info}, \code{debug}. The last is more verbose than the other two.
      The default value is set to \code{info}.

    \item \code{TimeLimit}: set the time limit for computing the distance between a pair of spatial maps. Min values: \code{INTMAX}.
          The default value is set to \code{INTMAX}.

    \item \code{OptTolerance}: Optimality tolerance on negative reduced cost variables to enter the basis.
          Min value: \eqn{10^{-9}}, max value: \eqn{10^{-1}}.
          The default value is set to \eqn{10^{-6}}.
    }
}
\seealso{
See also \code{\link{compareOneToOne}}, \code{\link{compareOneToMany}}, \code{\link{compareAll}}, \code{\link{focusArea}}, and \code{\link{Histogram2D}}.
}
\keyword{classes}
\section{Methods}{
  \describe{
    \item{\code{compareExact(Xs, Ys, W1, W2)}:}{compute the exact distance between the two vector of weights \code{W1} and \code{W2}, on the convex hull of the points defined by the two vectors \code{Xs} and \code{Ys}.
    The algorithm used by the solver is controlled by the parameter \code{ExactMethod} (see below).
    This method returns a single value (double), which is the KW-distance between \code{W1} and \code{W2}.}

    \item{\code{compareExact(Xs, Ys, W1, Ws)}:}{compute the exact distances between the vector of weights \code{W1} and each of the vector of weights in \code{Ws}, on the convex hull of the points defined by the two vectors \code{Xs} and \code{Ys}.
    The algorithm used by the solver is controlled by the parameter \code{ExactMethod} (see below).
    This method returns a vector of double of the same size of \code{Ws}, representing the distance of \code{W1} to every element of \code{Ws}.}

    \item{\code{compareExact(Xs, Ys, Ws)}:}{compute a symmetric matrix of pairwise exact distances between all the possible pairs of the vector listed in \code{Ws}.
    The algorithm used by the solver is controlled by the parameter \code{ExactMethod} (see below).}

    \item{\code{compareApprox(Xs, Ys, W1, W2, L)}:}{compute the approximate distance between the two vector of weights \code{W1} and \code{W2}, on the convex hull of the points defined by the two vectors \code{Xs} and \code{Ys}.
    The parameter \code{ApproxMethod} (see below) controls the algorithm used by the solver.
    This method returns a single value (double), which is the KW-distance between \code{W1} and \code{W2}.}

    \item{\code{compareApprox(Xs, Ys, W1, Ws, L)}:}{compute the approximate distances between the vector of weights \code{W1} and each of the vector of weights in \code{Ws}, on the convex hull of the points defined by the two vectors \code{Xs} and \code{Ys}.
    The parameter \code{ApproxMethod} (see below) controls the algorithm used by the solver.
    This method returns a vector of double of the same size of \code{Ws}, representing the distance of \code{W1} to every element of \code{Ws}.}

    \item{\code{compareApprox(Xs, Ys, Ws, L)}:}{compute a symmetric matrix of pairwise approximate distances (which depends on the value of \emph{L}) between all the possible pairs of the vector listed in \code{Ws}.
    The parameter \code{ApproxMethod} (see below) controls the algorithm used by the solver.
    }
    
    \item{\code{runtime()}:}{return the runtime in seconds to the last call to one of the \emph{compare} methods. It reports the runtime of the execution of the Network Simplex algorithm.}

    \item{\code{preprocesstime()}:}{return the preprocessing time in seconds to the last call to one of the \emph{compare} methods. It reports the execution time to set up the main data structures and to compute the convex hull of all the input histograms.}

    \item{\code{setParam(name, value)}:}{set the parameter \code{name} to the new \code{value}. Every parameter has a default value. See below for the existing parameters.}

    \item{\code{getParam(name)}:}{return the current value of the parameter \code{name}.}
  }

}
\arguments{
  \item{n}{Number of bins in the histograms \code{Xs, Yw, W1, W2, Ws}.}
  \item{H1}{First object of type \code{Histogram2D}.}
  \item{H2}{Second object of type \code{Histogram2D}.}
  \item{L}{Approximation parameter. Higher values of \emph{L} give a more accurate solution, but they require a longer running time.
  Table X gives the guarantee approximation bound as a function of \emph{L}. Type: positive integer.}
  \item{Xs}{Vector of horizontal coordinates the bins. Type: vector of integers.}
  \item{Ys}{Vector of vertical coordinates the bins. Type: vector of integers.}
  \item{W1}{Vector of weights of the bin at the positions specified by \code{Xs} and \code{Ys}. Type: vector of doubles.}
  \item{W2}{Vector of weights of the bin at the positions specified by \code{Xs} and \code{Ys}. Type: vector of doubles.}
  \item{Ws}{Matrix of weights of the bin at the positions specified by \code{Xs} and \code{Ys}. Type: matrix of doubles.}
  \item{name}{Name of the parameter to set and/or get. Type: string.}
  \item{value}{Value to set the corresponding parameter specified by \code{name}. Type: double.}
}
