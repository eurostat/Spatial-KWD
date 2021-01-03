// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
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


double distanceDF(const Rcpp::DataFrame& DF, int L) {
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
	}
	catch (...) {
		throw(Rcpp::exception("Error 13:", "KWD_NetSimplex", 13));
		return 0.0;
	}
}



RCPP_MODULE(SKWD) {
	using namespace Rcpp;

	function("distanceDF", &distanceDF, "documentation for distanceDF");


	class_<KWD::Histogram2D>("Histogram2D")
		// expose the default constructor
		.constructor()

		.method("add", &KWD::Histogram2D::add, "add an non empty support point")
		.method("update", &KWD::Histogram2D::update, "update an non empty support point")
		.method("size", &KWD::Histogram2D::size, "return the number of nonempty points")
		.method("balance", &KWD::Histogram2D::balance, "return the total sum of all the weights")
		.method("normalize", &KWD::Histogram2D::normalize, "normalize the weights to sum them up to one")
		;

	class_<KWD::Solver>("Solver")
		// expose the default constructor
		.constructor()

		.method("distance", &KWD::Solver::distance, 
			"compute the distance between a pair of histograms with given L")
		
		.method("column_generation", &KWD::Solver::column_generation, 
			"compute the distance between a pair of histograms with given L using column generation")

		.method("dense", &KWD::Solver::dense, 
			"compute the distance between a pair of histograms with given L using a bipartite graph (slow on large instances)")

		.method("runtime", &KWD::Solver::runtime, 
			"get the runtime in seconds of Network Simplex algorithm")
		; 
		
}
