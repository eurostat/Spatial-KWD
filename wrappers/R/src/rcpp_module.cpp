// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil;
// -*-
//
// rcpp_module.cpp: Rcpp R/C++ interface class library -- Rcpp Module examples
//
// Copyright (C) 2010 - 2012  Dirk Eddelbuettel and Romain Francois
//
// This file is part of Rcpp.
//
// Rcpp is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Rcpp is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

#include <Rcpp.h>

#include "KWD_Histogram2D.h"

RCPP_EXPOSED_AS(KWD::Histogram2D)

double distanceDF(const Rcpp::DataFrame &DF, int L) {
  Rcpp::IntegerVector X = DF["x"];
  Rcpp::IntegerVector Y = DF["y"];
  Rcpp::NumericVector H1 = DF["h1"];
  Rcpp::NumericVector H2 = DF["h2"];

  KWD::Histogram2D a;
  KWD::Histogram2D b;

  for (int i = 0, i_max = X.size(); i < i_max; ++i) {
    a.add(X[i], Y[i], H1[i]);
    b.add(X[i], Y[i], H2[i]);
  }

  a.normalize();
  b.normalize();

  KWD::Solver s;

  try {
    return s.column_generation(a, b, L);
  } catch (...) {
    throw(Rcpp::exception("Error 13:", "KWD_NetSimplex", 13));
    return 0.0;
  }
}

Rcpp::NumericVector compareExact(const Rcpp::DataFrame &DF) {
  Rcpp::NumericVector ds;
  if (!(DF.containsElementNamed("Xs") && DF.containsElementNamed("Ys"))) {
    throw(Rcpp::exception("The Dataframe must contain at least vectors Xs and "
                          "Ys, and a subset of W1, W2, and/or Ws."));
  }

  Rcpp::IntegerVector Xs = DF["Xs"];
  Rcpp::IntegerVector Ys = DF["Ys"];

  int n = Xs.size();

  KWD::Solver s;
  s.setParam(KWD_METHOD, KWD_EXACT);
  s.setParam(KWD_ALGORITHM, KWD_COLGEN);

  try {
    if (DF.containsElementNamed("W1") && DF.containsElementNamed("W2")) {
      Rcpp::NumericVector W1 = DF["W1"];
      Rcpp::NumericVector W2 = DF["W2"];

      double d =
          s.compareExact(n, Xs.begin(), Ys.begin(), W1.begin(), W2.begin());
      ds.push_back(d);
      return ds;
    }
    if (DF.containsElementNamed("W1") && DF.containsElementNamed("Ws")) {
      Rcpp::NumericVector W1 = DF["W1"];
      Rcpp::NumericMatrix Ws = DF["Ws"];

      int m = Ws.ncol();
      vector<double> _ds = s.compareApprox(n, m, Xs.begin(), Ys.begin(),
                                           W1.begin(), Ws.begin(), n - 1);
      ds.import(_ds.begin(), _ds.end());
      return ds;
    }
    // TODO: add the case of a single matrix
    throw(Rcpp::exception(
        "The column format of this DataFrame is not supported: it must contain "
        "at least vectors Xs and Ys, and a subset of W1, W2, and/or Ws."));
  } catch (std::exception &e) {
    Rprintf("Error 13: Rcpp::NumericVector compareExact(const Rcpp::DataFrame "
            "&DF)\n");
    forward_exception_to_r(e);
  }
  return ds;
}

Rcpp::NumericVector compareApprox(const Rcpp::DataFrame &DF, int L) {
  Rcpp::NumericVector ds;
  if (!(DF.containsElementNamed("Xs") && DF.containsElementNamed("Ys"))) {
    throw(Rcpp::exception("The Dataframe must contain at least vectors Xs and "
                          "Ys, and a subset of W1, W2, and/or Ws."));
  }

  Rcpp::IntegerVector Xs = DF["Xs"];
  Rcpp::IntegerVector Ys = DF["Ys"];

  int n = Xs.size();

  KWD::Solver s;
  s.setParam(KWD_METHOD, KWD_APPROX);
  s.setParam(KWD_ALGORITHM, KWD_COLGEN);

  try {
    if (DF.containsElementNamed("W1") && DF.containsElementNamed("W2")) {
      Rcpp::NumericVector W1 = DF["W1"];
      Rcpp::NumericVector W2 = DF["W2"];

      double d =
          s.compareApprox(n, Xs.begin(), Ys.begin(), W1.begin(), W2.begin(), L);
      ds.push_back(d);
      return ds;
    }
    if (DF.containsElementNamed("W1") && DF.containsElementNamed("Ws")) {
      Rcpp::NumericVector W1 = DF["W1"];
      Rcpp::NumericMatrix Ws = DF["Ws"];

      int m = Ws.ncol();
      vector<double> _ds = s.compareApprox(n, m, Xs.begin(), Ys.begin(),
                                           W1.begin(), Ws.begin(), L);
      ds.import(_ds.begin(), _ds.end());
      return ds;
    }
    // TODO: add the case of a single matrix
    throw(Rcpp::exception(
        "The column format of this DataFrame is not supported: it must contain "
        "at least vectors Xs and Ys, and a subset of W1, W2, and/or Ws."));
  } catch (std::exception &e) {
    Rprintf("Error 13: Rcpp::NumericVector compareApprox(const Rcpp::DataFrame "
            "&DF, int L)\n");
    forward_exception_to_r(e);
  }
  return ds;
}

Rcpp::NumericVector compareOneToManyApprox(const Rcpp::DataFrame &DF,
                                           Rcpp::NumericMatrix &Mat, int L) {
  Rcpp::NumericVector ds;
  if (!(DF.containsElementNamed("Xs") && DF.containsElementNamed("Ys") &&
        DF.containsElementNamed("W1"))) {
    throw(Rcpp::exception(
        "The Dataframe must contain at vectors Xs, Ys and Ws.\n"));
  }

  Rcpp::IntegerVector Xs = DF["Xs"];
  Rcpp::IntegerVector Ys = DF["Ys"];
  Rcpp::NumericVector W1 = DF["W1"];
  Rcpp::NumericVector Ws = Mat.row(0);

  int n = Xs.size();
  int m = Mat.ncol();

  KWD::Solver s;
  s.setParam(KWD_METHOD, KWD_APPROX);
  s.setParam(KWD_ALGORITHM, KWD_COLGEN);

  try {
    vector<double> _ds = s.compareApprox(n, m, Xs.begin(), Ys.begin(),
                                         W1.begin(), Mat.begin(), L);
    for (auto v : _ds)
      ds.push_back(v);

    return ds;
  } catch (std::exception &e) {
    Rprintf("Error 13: Rcpp::NumericVector compareOneToMayApprox()\n");
    forward_exception_to_r(e);
  }
  return ds;
}

RCPP_MODULE(SKWD) {
  using namespace Rcpp;

  function("distanceDF", &distanceDF, "compare two Histogram2D");

  function("compareExact", &compareExact, "compare histograms with exact OT");

  function("compareApprox", &compareApprox,
           "compare histograms with compareApprox OT");

  function("compareOneToManyApprox", &compareOneToManyApprox,
           "compare one to many histograms with approx OT algorithm");

  class_<KWD::Histogram2D>("Histogram2D")
      // expose the default constructor
      .constructor()

      .method("add", &KWD::Histogram2D::add, "add an non empty support point")
      .method("update", &KWD::Histogram2D::update,
              "update an non empty support point")
      .method("size", &KWD::Histogram2D::size,
              "return the number of nonempty points")
      .method("balance", &KWD::Histogram2D::balance,
              "return the total sum of all the weights")
      .method("normalize", &KWD::Histogram2D::normalize,
              "normalize the weights to sum them up to one");

  class_<KWD::Solver>("Solver")
      // expose the default constructor
      .constructor()

      // Solution methods
      .method("distance", &KWD::Solver::distance,
              "compute the distance between a pair of histograms with given L")

      .method("column_generation", &KWD::Solver::column_generation,
              "compute the distance between a pair of histograms with given L "
              "using column generation")

      .method("dense", &KWD::Solver::dense,
              "compute the distance between a pair of histograms with given L "
              "using a bipartite graph (slow on large instances)")

      //.method("compareExact", &KWD::Solver::compareExact,
      //	"compare two histograms in vector format with an exact
      // algortihm")

      //.method("compareApprox", &KWD::Solver::compareApprox,
      //	"compare two histograms in vector format with an approximate
      // algortihm")

      // Paramaters and attributes
      .method("runtime", &KWD::Solver::runtime,
              "get the runtime in seconds of Network Simplex algorithm")

      .method("iterations", &KWD::Solver::iterations,
              "get the number of iterations of Network Simplex algorithm")

      .method("num_arcs", &KWD::Solver::num_arcs,
              "get the number of arcs in the Network model")

      .method("num_nodes", &KWD::Solver::num_nodes,
              "get the number of arcs in the Network model")

      .method("status", &KWD::Solver::status,
              "get the status of Network Simplex solver")

      .method("setParam", &KWD::Solver::setParam,
              "set a parameter of the Network Simplex solver")

      .method("getParam", &KWD::Solver::getParam,
              "set a parameter of the Network Simplex solver");
}
