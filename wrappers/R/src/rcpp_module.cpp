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

//double distanceDF(const Rcpp::DataFrame& DF, int L) {
//	Rcpp::IntegerVector X = DF["x"];
//	Rcpp::IntegerVector Y = DF["y"];
//	Rcpp::NumericVector H1 = DF["h1"];
//	Rcpp::NumericVector H2 = DF["h2"];
//
//	KWD::Histogram2D a;
//	KWD::Histogram2D b;
//
//	for (int i = 0, i_max = X.size(); i < i_max; ++i) {
//		a.add(X[i], Y[i], H1[i]);
//		b.add(X[i], Y[i], H2[i]);
//	}
//
//	a.normalize();
//	b.normalize();
//
//	KWD::Solver s;
//
//	try {
//		return s.column_generation(a, b, L);
//	}
//	catch (...) {
//		throw(Rcpp::exception("Error 13:", "KWD_NetSimplex", 13));
//		return 0.0;
//	}
//}

Rcpp::List compareOneToOne(Rcpp::NumericMatrix& Coordinates, Rcpp::NumericMatrix& Weigths,
	int L = 3, bool recode = false,
	const std::string& method = "approx", const std::string& algorithm = "colgen", const std::string& model = "mincostflow", const std::string& verbosity = "silent",
	double timelimit = 14400, double opt_tolerance = 1e-06) {
	Rcpp::List sol;
	if (Coordinates.ncol() != 2)
		throw(Rcpp::exception("The Coordinates matrix must contain two columns for Xs and Ys."));

	if (Weigths.ncol() < 2)
		throw(Rcpp::exception("The Weigths matrix must contain two columns for W1 and W1."));

	if (Weigths.ncol() > 2)
		Rprintf("WARNING: only the first two columns of matrix Weights are used as histograms.");

	// Input data
	int n = Coordinates.nrow();

	vector<int> data1 = Rcpp::as<vector<int>>(Coordinates);
	int* Xs = &data1[0];
	int* Ys = &data1[n];
	vector<double> data2 = Rcpp::as<vector<double>>(Weigths);
	double* W1 = &data2[0];
	double* W2 = &data2[n];


	// Elaborate input parameters
	int _L = 3;
	if (L < 1)
		Rprintf("WARNING: Paramater L can take only value greater than 1. Using "
			"default value L=3.");
	else
		_L = L;

	KWD::Solver s;
	s.setStrParam(KWD_PAR_METHOD, method);
	s.setStrParam(KWD_PAR_MODEL, model);
	s.setStrParam(KWD_PAR_ALGORITHM, algorithm);
	s.setStrParam(KWD_PAR_VERBOSITY, verbosity);
	s.setDblParam(KWD_PAR_OPTTOLERANCE, opt_tolerance);
	s.setDblParam(KWD_PAR_TIMELIMIT, timelimit);
	if (recode)
		s.setStrParam(KWD_PAR_RECODE, "true");

	try {
		double d = -1;
		if (method == KWD_VAL_APPROX) {
			Rprintf("CompareOneToOne, Solution method: APPROX\n");
			d = s.compareApprox(n, Xs, Ys, W1, W2, _L);
		}
		else {
			Rprintf("CompareOneToOne, Solution method: EXACT\n");
			d = s.compareExact(n, Xs, Ys, W1, W2);
		}
		sol = Rcpp::List::create(
			Rcpp::Named("distance") = d, Rcpp::Named("runtime") = s.runtime(),
			Rcpp::Named("iterations") = s.iterations(),
			Rcpp::Named("nodes") = s.num_nodes(),
			Rcpp::Named("arcs") = s.num_arcs(), Rcpp::Named("status") = s.status());
	}
	catch (std::exception& e) {
		Rprintf("Error 13: Rcpp::NumericVector compareOneToOne()\n");
		forward_exception_to_r(e);
	}
	return sol;
}

Rcpp::List compareOneToMany(Rcpp::NumericMatrix& Coordinates, Rcpp::NumericMatrix& Weigths,
	int L = 3, bool recode = false,
	const std::string& method = "approx", const std::string& algorithm = "colgen",
	const std::string& model = "mincostflow", const std::string& verbosity = "silent",
	double timelimit = 14400, double opt_tolerance = 1e-06)
{
	Rcpp::List sol;
	if (Coordinates.ncol() != 2)
		throw(Rcpp::exception("The Coordinates matrix must contain two columns for Xs and Ys."));

	if (Weigths.ncol() < 2)
		throw(Rcpp::exception("The Weigths matrix must contain at least two columns."));

	// Input data
	int n = Coordinates.nrow();
	int m = Weigths.ncol() - 1;

	vector<int> data1 = Rcpp::as<vector<int>>(Coordinates);
	int* Xs = &data1[0];
	int* Ys = &data1[n];

	vector<double> data2 = Rcpp::as<vector<double>>(Weigths);
	double* W1 = &data2[0];
	double* Ws = &data2[n];

	// Elaborate input parameters
	int _L = 3;
	if (L < 1)
		Rprintf("WARNING: Paramater L can take only value greater than 1. Using "
			"default value L=3.");
	else
		_L = L;

	KWD::Solver s;
	s.setStrParam(KWD_PAR_METHOD, method);
	s.setStrParam(KWD_PAR_MODEL, model);
	s.setStrParam(KWD_PAR_ALGORITHM, algorithm);
	s.setStrParam(KWD_PAR_VERBOSITY, verbosity);
	s.setDblParam(KWD_PAR_OPTTOLERANCE, opt_tolerance);
	s.setDblParam(KWD_PAR_TIMELIMIT, timelimit);
	if (recode)
		s.setStrParam(KWD_PAR_RECODE, "true");

	try {
		Rcpp::NumericVector ds;
		if (method == KWD_VAL_APPROX) {
			Rprintf("CompareOneToMany, Solution method: APPROX\n");
			vector<double> _ds = s.compareApprox(n, m, Xs, Ys, W1, Ws, _L);
			for (auto v : _ds)
				ds.push_back(v);
		}
		else {
			Rprintf("CompareOneToMany, Solution method: EXACT\n");
			vector<double> _ds =
				s.compareApprox(n, m, Xs, Ys, W1, Ws, n - 1);
			for (auto v : _ds)
				ds.push_back(v);
		}
		sol = Rcpp::List::create(
			Rcpp::Named("distance") = ds, Rcpp::Named("runtime") = s.runtime(),
			Rcpp::Named("iterations") = s.iterations(),
			Rcpp::Named("nodes") = s.num_nodes(),
			Rcpp::Named("arcs") = s.num_arcs(), Rcpp::Named("status") = s.status());
	}
	catch (std::exception& e) {
		Rprintf("Error 13: Rcpp::NumericVector compareOneToMany()\n");
		forward_exception_to_r(e);
	}
	return sol;
}

Rcpp::List compareAll(Rcpp::NumericMatrix& Coordinates, Rcpp::NumericMatrix& Weigths,
	int L = 3, bool recode = false,
	const std::string& method = "approx", const std::string& algorithm = "colgen",
	const std::string& model = "mincostflow", const std::string& verbosity = "silent",
	double timelimit = 14400, double opt_tolerance = 1e-06)
{
	Rcpp::List sol;
	if (Coordinates.ncol() != 2)
		throw(Rcpp::exception("The Coordinates matrix must contain two columns for Xs and Ys."));

	if (Weigths.ncol() < 2)
		throw(Rcpp::exception("The Weigths matrix must contain at least two columns."));

	// Input data
	int n = Coordinates.nrow();
	int m = Weigths.ncol();

	vector<int> data1 = Rcpp::as<vector<int>>(Coordinates);
	int* Xs = &data1[0];
	int* Ys = &data1[n];

	vector<double> data2 = Rcpp::as<vector<double>>(Weigths);
	double* Ws = &data2[0];

	// Elaborate input parameters
	int _L = 3;
	if (L < 1)
		Rprintf("WARNING: Paramater L can take only value greater than 1. Using "
			"default value L=3.");
	else
		_L = L;

	KWD::Solver s;
	s.setStrParam(KWD_PAR_METHOD, method);
	s.setStrParam(KWD_PAR_MODEL, model);
	s.setStrParam(KWD_PAR_ALGORITHM, algorithm);
	s.setStrParam(KWD_PAR_VERBOSITY, verbosity);
	s.setDblParam(KWD_PAR_OPTTOLERANCE, opt_tolerance);
	s.setDblParam(KWD_PAR_TIMELIMIT, timelimit);
	if (recode)
		s.setStrParam(KWD_PAR_RECODE, "true");

	try {
		Rcpp::NumericMatrix ds;
		if (method == KWD_VAL_APPROX) {
			Rprintf("CompareAll, Solution method: APPROX\n");
			vector<double> _ds =
				s.compareApprox(n, m, Xs, Ys, Ws, _L);
			for (auto v : _ds)
				ds.push_back(v);
		}
		else {
			Rprintf("CompareAll, Solution method: EXACT\n");
			vector<double> _ds =
				s.compareApprox(n, m, Xs, Ys, Ws, n - 1);
			for (auto v : _ds)
				ds.push_back(v);
		}
		sol = Rcpp::List::create(
			Rcpp::Named("distance") = ds, Rcpp::Named("runtime") = s.runtime(),
			Rcpp::Named("iterations") = s.iterations(),
			Rcpp::Named("nodes") = s.num_nodes(),
			Rcpp::Named("arcs") = s.num_arcs(), Rcpp::Named("status") = s.status());
	}
	catch (std::exception& e) {
		Rprintf("Error 13: Rcpp::NumericVector compareAll()\n");
		forward_exception_to_r(e);
	}
	return sol;
}


RCPP_MODULE(SKWD) {
	using namespace Rcpp;

	function("compareOneToOne", &compareOneToOne,
		List::create(_["Coordinates"], _["Weights"], _["L"] = 3, _["recode"] = false,
			_["method"] = "approx", _["algorithm"] = "colgen", _["model"] = "mincostflow", _["verbosity"] = "silent",
			_["timelimit"] = 14400, _["opt_tolerance"] = 1e-06),
		"compare two histograms using the given search options");

	function("compareOneToMany", &compareOneToMany,
		List::create(_["Coordinates"], _["Weights"], _["L"] = 3, _["recode"] = false,
			_["method"] = "approx", _["algorithm"] = "colgen", _["model"] = "mincostflow", _["verbosity"] = "silent",
			_["timelimit"] = 14400, _["opt_tolerance"] = 1e-06),
		"compare one to many histograms using the given search options");

	function("compareAll", &compareAll,
		List::create(_["Coordinates"], _["Weights"], _["L"] = 3, _["recode"] = false,
			_["method"] = "approx", _["algorithm"] = "colgen", _["model"] = "mincostflow", _["verbosity"] = "silent",
			_["timelimit"] = 14400, _["opt_tolerance"] = 1e-06),
		"compare all histograms using the given search options");

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

		.method("setDblParam", &KWD::Solver::setDblParam,
			"set a double parameter of the Network Simplex solver")

		.method("getDblParam", &KWD::Solver::getDblParam,
			"get a double parameter of the Network Simplex solver")

		.method("setStrParam", &KWD::Solver::setStrParam,
			"set a string parameter of the Network Simplex solver")

		.method("getStrParam", &KWD::Solver::getStrParam,
			"get a string parameter of the Network Simplex solver");
}
