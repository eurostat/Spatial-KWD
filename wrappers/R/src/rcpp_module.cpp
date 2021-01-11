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

Rcpp::DataFrame compareOneToOne(const Rcpp::DataFrame& Data, const Rcpp::DataFrame& Options) {
	Rcpp::DataFrame sol;
	if (!(Data.containsElementNamed("Xs") && Data.containsElementNamed("Ys") &&
		Data.containsElementNamed("W1") && Data.containsElementNamed("W2"))) {
		throw(Rcpp::exception("The Dataframe must contain the vectors Xs, "
			"Ys, W1,  and W2."));
	}

	Rcpp::IntegerVector Xs = Data["Xs"];
	Rcpp::IntegerVector Ys = Data["Ys"];
	Rcpp::NumericVector W1 = Data["W1"];
	Rcpp::NumericVector W2 = Data["W2"];

	int L = 3;
	if (Options.containsElementNamed("L")) {
		L = Rcpp::as<int>(Options["L"]);
		if (L < 1)
			Rprintf("WARNING: Paramater L can take only value greater than 1. Using default value L=3.");
	}

	std::string method = "KWD_APPROX";
	if (Options.containsElementNamed("Method"))
		method = Rcpp::as<std::string>(Options["Method"]);

	int n = Xs.size();

	KWD::Solver s;
	s.setParam(KWD_ALGORITHM, KWD_COLGEN);

	try {
		double d = -1;
		if (method == "KWD_APPROX") {
			Rprintf("Solution method: KWD_APPROX\n");
			s.setParam(KWD_METHOD, KWD_APPROX);

			d = s.compareApprox(n, Xs.begin(), Ys.begin(), W1.begin(), W2.begin(), L);
		}
		else {
			Rprintf("Solution method: KWD_EXACT\n");
			s.setParam(KWD_METHOD, KWD_EXACT);
			d = s.compareExact(n, Xs.begin(), Ys.begin(), W1.begin(), W2.begin());
		}
		sol = Rcpp::DataFrame::create(
			Rcpp::Named("distance") = d,
			Rcpp::Named("runtime") = s.runtime(),
			Rcpp::Named("iterations") = s.iterations(),
			Rcpp::Named("nodes") = s.num_nodes(),
			Rcpp::Named("arcs") = s.num_arcs(),
			Rcpp::Named("status") = s.status());
	}
	catch (std::exception& e) {
		Rprintf("Error 13: Rcpp::NumericVector compareOneToOne()\n");
		forward_exception_to_r(e);
	}
	return sol;
}

Rcpp::List compareOneToMany(const Rcpp::DataFrame& Data,
	Rcpp::NumericMatrix& Ws, const Rcpp::DataFrame& Options)
{
	Rcpp::List sol;
	if (!(Data.containsElementNamed("Xs") && Data.containsElementNamed("Ys") &&
		Data.containsElementNamed("W1"))) {
		throw(Rcpp::exception(
			"The Dataframe must contain at vectors Xs, Ys and Ws.\n"));
	}

	Rcpp::IntegerVector Xs = Data["Xs"];
	Rcpp::IntegerVector Ys = Data["Ys"];
	Rcpp::NumericVector W1 = Data["W1"];

	int L = 3;
	if (Options.containsElementNamed("L")) {
		L = Rcpp::as<int>(Options["L"]);
		if (L < 1)
			Rprintf("WARNING: Paramater L can take only value greater than 1. Using default value L=3.");
	}

	std::string method = "KWD_APPROX";
	if (Options.containsElementNamed("Method"))
		method = Rcpp::as<std::string>(Options["Method"]);

	int n = Xs.size();
	int m = Ws.ncol();

	KWD::Solver s;
	s.setParam(KWD_ALGORITHM, KWD_COLGEN);

	try {
		Rcpp::NumericVector ds;
		if (method == "KWD_APPROX") {
			Rprintf("Solution method: KWD_APPROX\n");
			s.setParam(KWD_METHOD, KWD_APPROX);
			vector<double> _ds = s.compareApprox(n, m, Xs.begin(), Ys.begin(),
				W1.begin(), Ws.begin(), L);
			for (auto v : _ds)
				ds.push_back(v);
		}
		else {
			Rprintf("Solution method: KWD_EXACT\n");
			s.setParam(KWD_METHOD, KWD_EXACT);
			vector<double> _ds = s.compareApprox(n, m, Xs.begin(), Ys.begin(),
				W1.begin(), Ws.begin(), n - 1); // Prepare the correct method
			for (auto v : _ds)
				ds.push_back(v);
		}
		sol = Rcpp::List::create(
			Rcpp::Named("distances") = ds,
			Rcpp::Named("runtime") = s.runtime(),
			Rcpp::Named("iterations") = s.iterations(),
			Rcpp::Named("nodes") = s.num_nodes(),
			Rcpp::Named("arcs") = s.num_arcs(),
			Rcpp::Named("status") = s.status());
	}
	catch (std::exception& e) {
		Rprintf("Error 13: Rcpp::NumericVector compareOneToMayApprox()\n");
		forward_exception_to_r(e);
	}
	return sol;
}

RCPP_MODULE(SKWD) {
	using namespace Rcpp;

	function("distanceDF", &distanceDF,
		List::create(_["DF"], _["L"]),
		"compare two stored in the Dataframe DF");

	/*function("compareExact", &compareExact, "compare histograms with exact OT");

	function("compareApprox", &compareApprox,
		List::create(_["Xs"], _["Ys"], _["W1"], _["W2"], _["L"]),
		"compare histograms with compareApprox OT");

	function("compareApprox", &compareApprox,
		List::create(_["Xs"], _["Ys"], _["W1"], _["Ws"], _["L"]),
		"compare histograms with compareApprox OT");

	function("compareApprox", &compareApprox,
		List::create(_["Xs"], _["Ys"], _["W2"], _["L"]),
		"compare histograms with compareApprox OT");*/

	function("compareOneToOne", &compareOneToOne,
		List::create(_["Data"], _["Options"]),
		"compare two histograms with given search options");

	function("compareOneToMany", &compareOneToMany,
		List::create(_["Data"], _["Ws"], _["Options"]),
		"compare one to many histograms with given search options");

	//function("compareManyToMany", &compareManyToMany,
	//	List::create(_["DF"], _["Mat"], _["L"]),
	//	"compare many to many histograms with given search options");

	class_<KWD::Histogram2D>("Histogram2D")
		// expose the default constructor
		.constructor()
		//.constructor<int, int*, int*, double*>()

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
