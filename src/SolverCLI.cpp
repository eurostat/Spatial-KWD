/*
 * @fileoverview Copyright (c) 2019-2020, Stefano Gualandi,
 *               via Ferrata, 5, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#include <random>

#include "KWD_Histogram2D.h"

// Basic test for focus area
bool TEST_1_FOCUS_AREA() {
  vector<vector<int>> A = {
      {1, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 1, 0}, {0, 1, 1, 0}};

  vector<vector<int>> B = {
      {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 1}};

  size_t n = 4;
  vector<int> Xs(n * n, 0);
  vector<int> Ys(n * n, 0);
  vector<double> W1(n * n, 0.0);
  vector<double> W2(n * n, 0.0);

  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j) {
      Xs[i * n + j] = i;
      Ys[i * n + j] = j;
      W1[i * n + j] = A[i][j];
      W2[i * n + j] = B[i][j];
    }

  PRINT("Start TEST_1_FOCUS AREA\n");
  KWD::Solver solver;
  solver.setStrParam(KWD_PAR_METHOD, KWD_VAL_APPROX);
  solver.setStrParam(KWD_PAR_ALGORITHM, KWD_VAL_COLGEN);
  // solver.setStrParam(KWD_PAR_ALGORITHM, KWD_VAL_FULLMODEL);
  solver.setStrParam(KWD_PAR_VERBOSITY, KWD_VAL_DEBUG);
  solver.setStrParam(KWD_PAR_RECODE, KWD_VAL_TRUE);

  solver.dumpParam();

  double d =
      solver.focusArea(Xs.size(), &Xs[0], &Ys[0], &W1[0], &W2[0], 1, 1, 1, 3);

  double r = solver.runtime();
  auto s = solver.status();
  int its = solver.iterations();

  PRINT("Dense => %d: fobj: %.6f, time: %.2f, status: %s, iter: %d, arcs: "
        "%d, nodes: %d\n",
        (int)n, d, r, s.c_str(), its, solver.num_arcs(), solver.num_nodes());

  return true;
}

bool TEST_2_FOCUS_AREA() {
  // Esempio di prima all'interno di un'area maggiore
  // vector<vector<int>> A = {
  //{1, 0,  0,  0,  0,  0, 0, 0},
  //{0, 0,  1,  0,  0,  0, 0, 0},
  //{0, 0,  1,  0,  0,  0, 0, 0},
  //{0, 1,  1,  0,  0,  0, 0, 0},
  //{0, 0,  0,  0,  0,  0, 0, 0},
  //{0, 0,  0,  0,  0,  0, 0, 0},
  //{0, 0,  0,  0,  0,  0, 0, 0},
  //{0, 0,  0,  0,  0,  0, 0, 0}
  //};

  // vector<vector<int>> B = {
  //{0, 0,  0,  0,  0,  0, 0, 0},
  //{0, 0,  0,  0,  0,  0, 0, 0},
  //{0, 1,  0,  0,  0,  0, 0, 0},
  //{0, 0,  1,  1,  0,  0, 0, 0},
  //{0, 0,  0,  0,  0,  0, 0, 0},
  //{0, 0,  0,  0,  0,  0, 0, 0},
  //{0, 0,  0,  0,  0,  0, 0, 0},
  //{0, 0,  0,  0,  0,  0, 0, 0}
  //};

  vector<vector<int>> A = {{0, 0, 0, 0, 0, 0, 0, 0},  {0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 0, +1, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 0, 0, 0, 0, +2, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 0, +1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}};

  vector<vector<int>> B = {
      {0, 0, 0, +1, 1, 0, 0, 0}, {0, 0, +1, +1, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},  {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 1, 1},  {0, 0, +1, 0, 0, +1, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0},  {0, 0, 0, 0, +1, +1, 0, 0}};

  size_t n = A[0].size();
  vector<int> Xs(n * n, 0);
  vector<int> Ys(n * n, 0);
  vector<double> W1(n * n, 0.0);
  vector<double> W2(n * n, 0.0);

  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j) {
      Xs[i * n + j] = i;
      Ys[i * n + j] = j;
      W1[i * n + j] = A[i][j];
      W2[i * n + j] = B[i][j];
    }

  PRINT("Start TEST_2_FOCUS AREA\n");
  KWD::Solver solver;
  solver.setStrParam(KWD_PAR_METHOD, KWD_VAL_APPROX);
  solver.setStrParam(KWD_PAR_ALGORITHM, KWD_VAL_COLGEN);
  solver.setStrParam(KWD_PAR_VERBOSITY, KWD_VAL_DEBUG);
  solver.setStrParam(KWD_PAR_RECODE, KWD_VAL_TRUE);

  solver.dumpParam();

  double d =
      solver.focusArea(Xs.size(), &Xs[0], &Ys[0], &W1[0], &W2[0], 3, 3, 2, 4);

  double r = solver.runtime();
  auto s = solver.status();
  int its = solver.iterations();

  PRINT("Dense => %d: fobj: %.6f, time: %.2f, status: %s, iter: %d, arcs: "
        "%d, nodes: %d\n",
        (int)n, d, r, s.c_str(), its, solver.num_arcs(), solver.num_nodes());

  return true;
}

bool TEST_3_FOCUS_AREA() {
  // Esempio di prima all'interno di un'area maggiore

  vector<vector<int>> A = {{0, 0, 0, 0, 0, +1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0},  {0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0},  {0, 0, 0, 0, 0, 0, 0, 0},
                           {0, 0, +1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}};

  vector<vector<int>> B = {
      {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, +1, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, +1, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, +1, +1, 0, 0}};

  size_t n = A[0].size();
  vector<int> Xs(n * n, 0);
  vector<int> Ys(n * n, 0);
  vector<double> W1(n * n, 0.0);
  vector<double> W2(n * n, 0.0);

  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j) {
      Xs[i * n + j] = i;
      Ys[i * n + j] = j;
      W1[i * n + j] = A[i][j];
      W2[i * n + j] = B[i][j];
    }

  PRINT("Start TEST_3_FOCUS AREA\n");
  KWD::Solver solver;
  solver.setStrParam(KWD_PAR_METHOD, KWD_VAL_APPROX);
  solver.setStrParam(KWD_PAR_ALGORITHM, KWD_VAL_COLGEN);
  solver.setStrParam(KWD_PAR_VERBOSITY, KWD_VAL_DEBUG);
  solver.setStrParam(KWD_PAR_RECODE, KWD_VAL_TRUE);

  solver.dumpParam();

  double d =
      solver.focusArea(Xs.size(), &Xs[0], &Ys[0], &W1[0], &W2[0], 3, 3, 2, 4);

  double r = solver.runtime();
  auto s = solver.status();
  int its = solver.iterations();

  PRINT("Dense => %d: fobj: %.6f, time: %.2f, status: %s, iter: %d, arcs: "
        "%d, nodes: %d\n",
        (int)n, d, r, s.c_str(), its, solver.num_arcs(), solver.num_nodes());

  return true;
}

// Main function
int main(int argc, char *argv[]) {
  int n = 32;
  int seed = 13;

  std::random_device
      rd; // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(seed); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> Uniform01(0, 1);
  std::uniform_int_distribution<> Uniform0N(0, n);

  // Test focus area
  if (false) {
    TEST_1_FOCUS_AREA();
    TEST_2_FOCUS_AREA();
    TEST_3_FOCUS_AREA();
  }

  // Test from command line
  if (true) {
    int L = 3;
    std::string filename = "";
    if (argc > 1)
      filename = std::string(argv[1]);
    else {
      fprintf(stdout, "Usage: %s <path_to_data_file>\n", argv[0]);
      return 1;
    }
    // fprintf(stdout, "%s\n", filename.c_str());

    if (argc > 2)
      L = std::stoi(argv[2]);

    std::ifstream in_file(filename);

    if (!in_file) {
      fprintf(stdout, "FATAL ERROR: Cannot open file\n");
      return 1;
    }

    vector<int> Xs;
    vector<int> Ys;
    vector<double> W1;
    vector<double> W2;
    std::string line;

    char sep = ',';
    std::getline(in_file, line);

    while (std::getline(in_file, line)) {
      std::stringstream lineStream(line);
      std::string cell;

      std::getline(lineStream, cell, sep);
      int x = std::stoi(cell);
      // fprintf(stdout, "%d\n", x);
      // fflush(stdout);
      std::getline(lineStream, cell, sep);
      int y = std::stoi(cell);
      // fprintf(stdout, "%d\n", y);
      // fflush(stdout);
      std::getline(lineStream, cell, sep);
      double a = 0.0;
      try {
        a = std::stod(cell);
      } catch (...) {
        a = 0.0;
      }
      // fprintf(stdout, "%g\n", a);
      // fflush(stdout);
      std::getline(lineStream, cell, sep);
      double b = 0.0;
      try {
        b = std::stod(cell);
      } catch (...) {
        b = 0.0;
      }
      // fprintf(stdout, "%g\n", b);
      // fflush(stdout);

      Xs.push_back(x);
      Ys.push_back(y);
      W1.push_back(a);
      W2.push_back(b);
    }

    int n = Xs.size();
    PRINT("start solver\n");
    KWD::Solver solver;

    // ----------------------------------------------------------------------------
    solver.setStrParam(KWD_PAR_METHOD, KWD_VAL_APPROX);
    solver.setStrParam(KWD_PAR_ALGORITHM, KWD_VAL_COLGEN);
    solver.setStrParam(KWD_PAR_MODEL, KWD_VAL_MINCOSTFLOW);
    solver.setStrParam(KWD_PAR_VERBOSITY, KWD_VAL_SILENT);
    solver.setStrParam(KWD_PAR_RECODE, KWD_VAL_TRUE);

    solver.setStrParam(KWD_PAR_UNBALANCED, KWD_VAL_FALSE);
    // solver.setDblParam(KWD_PAR_UNBALANCED_COST, 4); // TO BE SET

    solver.setStrParam(KWD_PAR_CONVEXHULL, KWD_VAL_TRUE);

    // solver.dumpParam();

    if (true) {
      auto dist =
          solver.compareApprox(Xs.size(), &Xs[0], &Ys[0], &W1[0], &W2[0], L);

      PRINT("Approx => %d: fobj: %.6f, time: %.4f, status: %s, iter: %d, "
            "arcs: "
            "%d, nodes: %d, L: %d\n",
            n, dist, solver.runtime(), solver.status().c_str(),
            solver.iterations(), solver.num_arcs(), solver.num_nodes(), L);
    }

    int radius = 50;
    if (argc > 3)
      radius = std::stoi(argv[3]);

    if (false) {
      solver.setStrParam(KWD_PAR_AREA, KWD_VAL_CIRCULAR);

      int xx = 0;
      double xv = W2[xx];
      for (int i = 0; i < (int)Xs.size(); i++) {
        if (W2[i] > xv) {
          xx = i;
          xv = W2[i];
        }
      }

      auto dist = solver.focusArea(Xs.size(), &Xs[0], &Ys[0], &W1[0], &W2[0],
                                   Xs[xx], Ys[xx], radius, L);

      PRINT("focusArea => %d: fobj: %.6f, time: %.4f, status: %s, iter: %d, "
            "arcs: "
            "%d, nodes: %d\n",
            n, dist, solver.runtime(), solver.status().c_str(),
            solver.iterations(), solver.num_arcs(), solver.num_nodes());
    }
  }

  if (false) {
    size_t samples = n * n / 2;
    vector<int> Xs(samples, 0);
    vector<int> Ys(samples, 0);
    vector<double> W1(samples, 0);

    size_t m = 3;

    vector<double> Ws(samples * m, 0);
    for (size_t i = 0; i < samples; i++) {
      Xs[i] = Uniform0N(gen);
      Ys[i] = Uniform0N(gen);
      W1[i] = Uniform01(gen);
      // Matrix as a flat array
      Ws[i] = Uniform01(gen);
      Ws[samples + i] = Uniform01(gen);
      Ws[2 * samples + i] = Uniform01(gen);
    }

    PRINT("start solver\n");
    KWD::Solver solver;
    // ----------------------------------------------------------------------------
    solver.setStrParam(KWD_PAR_METHOD, KWD_VAL_APPROX);

    for (auto algo : {KWD_VAL_COLGEN, KWD_VAL_MINCOSTFLOW}) {
      solver.setStrParam(KWD_PAR_ALGORITHM, algo);
      auto dist =
          solver.compareApprox(Xs.size(), m, &Xs[0], &Ys[0], &W1[0], &Ws[0], 3);

      for (double d : dist)
        PRINT("Approx => %d: fobj: %.6f, time: %.4f, status: %s, iter: %d, "
              "arcs: "
              "%d, nodes: %d\n",
              n, d, solver.runtime(), solver.status().c_str(),
              solver.iterations(), solver.num_arcs(), solver.num_nodes());

      dist = solver.compareApprox(Xs.size(), m, &Xs[0], &Ys[0], &Ws[0], 3);

      for (double d : dist)
        PRINT("AllCmp => %d: fobj: %.6f, time: %.4f, status: %s, iter: %d, "
              "arcs: "
              "%d, nodes: %d\n",
              n, d, solver.runtime(), solver.status().c_str(),
              solver.iterations(), solver.num_arcs(), solver.num_nodes());
    }

    KWD::Histogram2D A(Xs.size(), &Xs[0], &Ys[0], &W1[0]);

    KWD::Histogram2D B(Xs.size(), &Xs[0], &Ys[0], &Ws[0]);
    KWD::Histogram2D C(Xs.size(), &Xs[0], &Ys[0], &Ws[samples]);
    KWD::Histogram2D D(Xs.size(), &Xs[0], &Ys[0], &Ws[2 * samples]);

    double d = solver.distance(A, B, 3);
    double r = solver.runtime();
    auto s = solver.status();
    int its = solver.iterations();

    PRINT("Dista => %d: fobj: %.6f, time: %.4f, status: %s, iter: %d, arcs: "
          "%d, nodes: %d\n",
          n, d, r, s.c_str(), its, solver.num_arcs(), solver.num_nodes());

    d = solver.distance(A, C, 3);
    r += solver.runtime();
    s = solver.status();
    its += solver.iterations();

    PRINT("Dista => %d: fobj: %.6f, time: %.4f, status: %s, iter: %d, arcs: "
          "%d, nodes: %d\n",
          n, d, r, s.c_str(), its, solver.num_arcs(), solver.num_nodes());

    d = solver.distance(A, D, 3);
    r += solver.runtime();
    s = solver.status();
    its += solver.iterations();

    PRINT("Dista => %d: fobj: %.6f, time: %.4f, status: %s, iter: %d, arcs: "
          "%d, nodes: %d\n",
          n, d, r, s.c_str(), its, solver.num_arcs(), solver.num_nodes());
  }

  if (false) {
    size_t samples = n * n / 2;
    vector<int> Xs(samples, 0);
    vector<int> Ys(samples, 0);
    vector<double> W1(samples, 0);
    vector<double> W2(samples, 0);
    for (size_t i = 0; i < samples; i++) {
      Xs[i] = Uniform0N(gen);
      Ys[i] = Uniform0N(gen);
      W1[i] = Uniform01(gen);
      W2[i] = Uniform01(gen);

      // fprintf(stdout, "%d %d %.4f %.4f\n", Xs[i], Ys[i], W1[i], W2[i]);
    }

    // vector<int> Xs = {0, 5};
    // vector<int> Ys = {0, 7};
    // vector<double> W1 = {1, 0};
    // vector<double> W2 = {0, 1};
    //    size_t samples = Xs.size();

    PRINT("start solver\n");
    KWD::Solver solver;
    // solver.setParam(KWD_VERBOSITY, KWD_DEBUG);
    KWD::Histogram2D A(Xs.size(), &Xs[0], &Ys[0], &W1[0]);
    KWD::Histogram2D B(Xs.size(), &Xs[0], &Ys[0], &W2[0]);

    double d = solver.dense(A, B);
    double r = solver.runtime();
    auto s = solver.status();
    int its = solver.iterations();

    PRINT("Dense => %d: fobj: %.6f, time: %.2f, status: %s, iter: %d, arcs: "
          "%d, nodes: %d\n",
          n, d, r, s.c_str(), its, solver.num_arcs(), solver.num_nodes());

    d = solver.distance(A, B, 3);
    r = solver.runtime();
    s = solver.status();
    its = solver.iterations();

    PRINT("Dista => %d: fobj: %.6f, time: %.2f, status: %s, iter: %d, arcs: "
          "%d, nodes: %d\n",
          n, d, r, s.c_str(), its, solver.num_arcs(), solver.num_nodes());

    d = solver.column_generation(A, B, 3);
    r = solver.runtime();
    s = solver.status();
    its = solver.iterations();

    PRINT("ColGe => %d: fobj: %.6f, time: %.2f, status: %s, iter: %d, arcs: "
          "%d, nodes: %d\n",
          n, d, r, s.c_str(), its, solver.num_arcs(), solver.num_nodes());

    // ----------------------------------------------------------------------------
    solver.setStrParam(KWD_PAR_METHOD, KWD_VAL_EXACT);
    for (auto algo : {KWD_VAL_MINCOSTFLOW, KWD_VAL_COLGEN}) {
      solver.setStrParam(KWD_PAR_ALGORITHM, algo);
      double dist =
          solver.compareExact(Xs.size(), &Xs[0], &Ys[0], &W1[0], &W2[0]);

      PRINT("Exact => %d: fobj: %.6f, time: %.2f, status: %s, iter: %d, arcs: "
            "%d, nodes: %d\n",
            n, dist, solver.runtime(), solver.status().c_str(),
            solver.iterations(), solver.num_arcs(), solver.num_nodes());
    }

    // ----------------------------------------------------------------------------
    solver.setStrParam(KWD_PAR_METHOD, KWD_VAL_APPROX);
    for (auto algo : {KWD_VAL_MINCOSTFLOW, KWD_VAL_COLGEN}) {
      solver.setStrParam(KWD_PAR_ALGORITHM, algo);
      double dist =
          solver.compareApprox(Xs.size(), &Xs[0], &Ys[0], &W1[0], &W2[0], 3);

      PRINT("Approx => %d: fobj: %.6f, time: %.2f, status: %s, iter: %d, arcs: "
            "%d, nodes: %d\n",
            n, dist, solver.runtime(), solver.status().c_str(),
            solver.iterations(), solver.num_arcs(), solver.num_nodes());
    }
  }

  if (false) {
    size_t samples = n * n / 2;
    KWD::Histogram2D a;
    KWD::Histogram2D b;

    for (size_t i = 0; i < samples; i++) {
      // a.add(-100 + Uniform0N(gen), -10 + Uniform0N(gen), Uniform01(gen));
      // b.add(-20 + Uniform0N(gen), -30 + Uniform0N(gen), Uniform01(gen));
      a.add(Uniform0N(gen), Uniform0N(gen), Uniform01(gen));
      b.add(Uniform0N(gen), Uniform0N(gen), Uniform01(gen));
    }

    a.normalize();
    b.normalize();

    PRINT("start solver\n");
    KWD::Solver solver;

    for (int L = 2; L <= 3; ++L) {
      double dist = solver.distance(a, b, L);
      PRINT("Full => %d: %d %.6f %.3f sec\n", n, L, dist, solver.runtime());
    }

    for (int L = 2; L <= 3; ++L) {
      PRINT("CG %d\n", L);
      double dist = solver.column_generation(a, b, L);
      PRINT("ColG => %d: %d %.6f %.3f sec\n", n, L, dist, solver.runtime());
    }

    // double dist = solver.dense(a, b);
    // fprintf(stdout, "%d: %.6f %.3f ms\n", n, dist, solver.runtime());
  }

  return EXIT_SUCCESS;
}
