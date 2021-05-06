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

#ifdef __GLIBC__
#define _POSIX_C_SOURCE 200809L
#endif

#include <Rcpp.h>

#include "KWD_Histogram2D.h"

RCPP_EXPOSED_AS(KWD::Histogram2D)

RCPP_MODULE(SKWD) {
	Rcpp::class_<KWD::Histogram2D>("Histogram2D")
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

	Rcpp::class_<KWD::Solver>("Solver")
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
			"get the number of nodes in the Network model")

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

// [[Rcpp::export(compareOneToOne)]]
Rcpp::List compareOneToOne(
	Rcpp::NumericMatrix Coordinates, Rcpp::NumericMatrix Weights, int L = 3,
	bool recode = true, const std::string& method = "approx",
	const std::string& algorithm = "colgen",
	const std::string& model = "mincostflow",
	const std::string& verbosity = "silent", double timelimit = 14400,
	double opt_tolerance = 1e-06, bool unbalanced = false,
	double unbal_cost = 1e+09, bool convex = true) {
	if (Coordinates.ncol() != 2)
		throw(Rcpp::exception(
			"The Coordinates matrix must contain two columns for Xs and Ys."));

	if (Weights.ncol() < 2)
		throw(Rcpp::exception(
			"The Weigths matrix must contain two columns for W1 and W1."));

	if (Weights.ncol() > 2)
		Rprintf("WARNING: only the first two columns of matrix Weights are used as "
			"histograms.");

	// Input data
	int n = Coordinates.nrow();

	std::vector<int> data1 = Rcpp::as<std::vector<int>>(Coordinates);
	int* Xs = &data1[0];
	int* Ys = &data1[n];
	std::vector<double> data2 = Rcpp::as<std::vector<double>>(Weights);
	double* W1 = &data2[0];
	double* W2 = &data2[n];

	// Elaborate input parameters
	int LL = 3;
	if (L < 1)
		Rprintf("WARNING: Paramater L can take only value greater than 1. Using "
			"default value L=3.");
	else
		LL = L;

	KWD::Solver s;
	s.setStrParam(KWD_PAR_METHOD, method);
	s.setStrParam(KWD_PAR_MODEL, model);
	s.setStrParam(KWD_PAR_ALGORITHM, algorithm);
	s.setStrParam(KWD_PAR_VERBOSITY, verbosity);
	s.setDblParam(KWD_PAR_OPTTOLERANCE, opt_tolerance);
	s.setDblParam(KWD_PAR_TIMELIMIT, timelimit);
	if (recode)
		s.setStrParam(KWD_PAR_RECODE, KWD_VAL_TRUE);

	if (unbalanced) {
		s.setStrParam(KWD_PAR_UNBALANCED, KWD_VAL_TRUE);
		s.setDblParam(KWD_PAR_UNBALANCED_COST, unbal_cost);
	}

	if (convex)
		s.setStrParam(KWD_PAR_CONVEXHULL, KWD_VAL_TRUE);

	try {
		double d = -1;
		if (method == KWD_VAL_APPROX) {
			Rprintf("CompareOneToOne, Solution method: APPROX\n");
			d = s.compareApprox(n, Xs, Ys, W1, W2, LL);
		}
		else {
			Rprintf("CompareOneToOne, Solution method: EXACT\n");
			d = s.compareExact(n, Xs, Ys, W1, W2);
		}
		return Rcpp::List::create(
			Rcpp::Named("distance") = d, Rcpp::Named("runtime") = s.runtime(),
			Rcpp::Named("iterations") = s.iterations(),
			Rcpp::Named("nodes") = s.num_nodes(),
			Rcpp::Named("arcs") = s.num_arcs(), Rcpp::Named("status") = s.status());
	}
	catch (std::exception& e) {
		Rprintf("Error 13: Rcpp::NumericVector compareOneToOne()\n");
		forward_exception_to_r(e);
	}
	return Rcpp::List::create(Rcpp::Named("status") = "FAILED");
}

// [[Rcpp::export(compareOneToMany)]]
Rcpp::List compareOneToMany(
	Rcpp::NumericMatrix Coordinates, Rcpp::NumericMatrix Weights, int L = 3,
	bool recode = true, const std::string& method = "approx",
	const std::string& algorithm = "colgen",
	const std::string& model = "mincostflow",
	const std::string& verbosity = "silent", double timelimit = 14400,
	double opt_tolerance = 1e-06, bool unbalanced = false,
	double unbal_cost = 1e+09, bool convex = true) {
	if (Coordinates.ncol() != 2)
		throw(Rcpp::exception(
			"The Coordinates matrix must contain two columns for Xs and Ys."));

	if (Weights.ncol() < 2)
		throw(Rcpp::exception(
			"The Weights matrix must contain at least two columns."));

	// Input data
	int n = Coordinates.nrow();
	int m = Weights.ncol() - 1;

	std::vector<int> data1 = Rcpp::as<std::vector<int>>(Coordinates);
	int* Xs = &data1[0];
	int* Ys = &data1[n];

	std::vector<double> data2 = Rcpp::as<std::vector<double>>(Weights);
	double* W1 = &data2[0];
	double* Ws = &data2[n];

	// Elaborate input parameters
	int LL = 3;
	if (L < 1)
		Rprintf("WARNING: Paramater L can take only value greater than 1. Using "
			"default value L=3.");
	else
		LL = L;

	KWD::Solver s;
	s.setStrParam(KWD_PAR_METHOD, method);
	s.setStrParam(KWD_PAR_MODEL, model);
	s.setStrParam(KWD_PAR_ALGORITHM, algorithm);
	s.setStrParam(KWD_PAR_VERBOSITY, verbosity);
	s.setDblParam(KWD_PAR_OPTTOLERANCE, opt_tolerance);
	s.setDblParam(KWD_PAR_TIMELIMIT, timelimit);
	if (recode)
		s.setStrParam(KWD_PAR_RECODE, KWD_VAL_TRUE);

	if (unbalanced) {
		s.setStrParam(KWD_PAR_UNBALANCED, KWD_VAL_TRUE);
		s.setDblParam(KWD_PAR_UNBALANCED_COST, unbal_cost);
	}

	if (convex)
		s.setStrParam(KWD_PAR_CONVEXHULL, KWD_VAL_TRUE);

	try {
		Rcpp::NumericVector ds;
		if (method == KWD_VAL_APPROX) {
			Rprintf("CompareOneToMany, Solution method: APPROX\n");
			vector<double> _ds = s.compareApprox(n, m, Xs, Ys, W1, Ws, LL);
			for (auto v : _ds)
				ds.push_back(v);
		}
		else {
			Rprintf("CompareOneToMany, Solution method: EXACT\n");
			vector<double> _ds = s.compareApprox(n, m, Xs, Ys, W1, Ws, n - 1);
			for (auto v : _ds)
				ds.push_back(v);
		}
		return Rcpp::List::create(
			Rcpp::Named("distance") = ds, Rcpp::Named("runtime") = s.runtime(),
			Rcpp::Named("iterations") = s.iterations(),
			Rcpp::Named("nodes") = s.num_nodes(),
			Rcpp::Named("arcs") = s.num_arcs(), Rcpp::Named("status") = s.status());
	}
	catch (std::exception& e) {
		Rprintf("Error 13: Rcpp::NumericVector compareOneToMany()\n");
		forward_exception_to_r(e);
	}
	return Rcpp::List::create(Rcpp::Named("status") = "FAILED");
}

// [[Rcpp::export(compareAll)]]
Rcpp::List compareAll(Rcpp::NumericMatrix Coordinates,
	Rcpp::NumericMatrix Weights, int L = 3,
	bool recode = true, const std::string& method = "approx",
	const std::string& algorithm = "colgen",
	const std::string& model = "mincostflow",
	const std::string& verbosity = "silent",
	double timelimit = 14400, double opt_tolerance = 1e-06,
	bool unbalanced = false, double unbal_cost = 1e+09,
	bool convex = true) {
	if (Coordinates.ncol() != 2)
		throw(Rcpp::exception(
			"The Coordinates matrix must contain two columns for Xs and Ys."));

	if (Weights.ncol() < 2)
		throw(Rcpp::exception(
			"The Weights matrix must contain at least two columns."));

	// Input data
	int n = Coordinates.nrow();
	int m = Weights.ncol();

	std::vector<int> data1 = Rcpp::as<std::vector<int>>(Coordinates);
	int* Xs = &data1[0];
	int* Ys = &data1[n];

	std::vector<double> data2 = Rcpp::as<std::vector<double>>(Weights);
	double* Ws = &data2[0];

	// Elaborate input parameters
	int LL = 3;
	if (L < 1)
		Rprintf("WARNING: Paramater L can take only value greater than 1. Using "
			"default value L=3.");
	else
		LL = L;

	KWD::Solver s;
	s.setStrParam(KWD_PAR_METHOD, method);
	s.setStrParam(KWD_PAR_MODEL, model);
	s.setStrParam(KWD_PAR_ALGORITHM, algorithm);
	s.setStrParam(KWD_PAR_VERBOSITY, verbosity);
	s.setDblParam(KWD_PAR_OPTTOLERANCE, opt_tolerance);
	s.setDblParam(KWD_PAR_TIMELIMIT, timelimit);
	if (recode)
		s.setStrParam(KWD_PAR_RECODE, KWD_VAL_TRUE);

	if (unbalanced) {
		s.setStrParam(KWD_PAR_UNBALANCED, KWD_VAL_TRUE);
		s.setDblParam(KWD_PAR_UNBALANCED_COST, unbal_cost);
	}

	if (convex)
		s.setStrParam(KWD_PAR_CONVEXHULL, KWD_VAL_TRUE);

	try {
		Rcpp::NumericMatrix ds;
		if (method == KWD_VAL_APPROX) {
			Rprintf("CompareAll, Solution method: APPROX\n");
			vector<double> _ds = s.compareApprox(n, m, Xs, Ys, Ws, LL);
			for (auto v : _ds)
				ds.push_back(v);
		}
		else {
			Rprintf("CompareAll, Solution method: EXACT\n");
			vector<double> _ds = s.compareApprox(n, m, Xs, Ys, Ws, n - 1);
			for (auto v : _ds)
				ds.push_back(v);
		}
		return Rcpp::List::create(
			Rcpp::Named("distance") = ds, Rcpp::Named("runtime") = s.runtime(),
			Rcpp::Named("iterations") = s.iterations(),
			Rcpp::Named("nodes") = s.num_nodes(),
			Rcpp::Named("arcs") = s.num_arcs(), Rcpp::Named("status") = s.status());
	}
	catch (std::exception& e) {
		Rprintf("Error 13: Rcpp::NumericVector compareAll()\n");
		forward_exception_to_r(e);
	}
	return Rcpp::List::create(Rcpp::Named("status") = "FAILED");
}

// [[Rcpp::export(focusArea)]]
Rcpp::List focusArea(Rcpp::NumericMatrix Coordinates,
	Rcpp::NumericMatrix Weights, int x, int y, int radius,
	int L = 3, bool recode = true,
	const std::string& method = "approx",
	const std::string& algorithm = "colgen",
	const std::string& model = "mincostflow",
	const std::string& verbosity = "silent",
	double timelimit = 14400, double opt_tolerance = 1e-06,
	const std::string& area = "l2") {
	if (Coordinates.ncol() != 2)
		throw(Rcpp::exception(
			"The Coordinates matrix must contain two columns for Xs and Ys."));

	if (Weights.ncol() < 2)
		throw(Rcpp::exception(
			"The Weigths matrix must contain two columns for W1 and W1."));

	if (Weights.ncol() > 2)
		Rprintf("WARNING: only the first two columns of matrix Weights are used as "
			"histograms.");

	// Input data
	int n = Coordinates.nrow();

	std::vector<int> data1 = Rcpp::as<std::vector<int>>(Coordinates);
	int* Xs = &data1[0];
	int* Ys = &data1[n];
	std::vector<double> data2 = Rcpp::as<std::vector<double>>(Weights);
	double* W1 = &data2[0];
	double* W2 = &data2[n];

	// Elaborate input parameters
	int LL = 3;
	if (L < 1)
		Rprintf("WARNING: Paramater L can take only value greater than 1. Using "
			"default value L=3.");
	else
		LL = L;

	KWD::Solver s;
	s.setStrParam(KWD_PAR_METHOD, method);
	s.setStrParam(KWD_PAR_MODEL, model);
	s.setStrParam(KWD_PAR_ALGORITHM, algorithm);
	s.setStrParam(KWD_PAR_VERBOSITY, verbosity);
	s.setDblParam(KWD_PAR_OPTTOLERANCE, opt_tolerance);
	s.setDblParam(KWD_PAR_TIMELIMIT, timelimit);
	if (recode)
		s.setStrParam(KWD_PAR_RECODE, KWD_VAL_TRUE);
	s.setStrParam(KWD_PAR_AREA, area);

	try {
		Rprintf("FocusArea, Solution method: APPROX\n");

		double d = s.focusArea(n, Xs, Ys, W1, W2, x, y, radius, LL);

		return Rcpp::List::create(
			Rcpp::Named("distance") = d, Rcpp::Named("runtime") = s.runtime(),
			Rcpp::Named("iterations") = s.iterations(),
			Rcpp::Named("nodes") = s.num_nodes(),
			Rcpp::Named("arcs") = s.num_arcs(), Rcpp::Named("status") = s.status());
	}
	catch (std::exception& e) {
		Rprintf("Error 13: Rcpp::NumericVector FocusArea()\n");
		forward_exception_to_r(e);
	}

	return Rcpp::List::create(Rcpp::Named("status") = "FAILED");
}
