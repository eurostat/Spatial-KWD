/*
 * @fileoverview Copyright (c) 2019-2020, Stefano Gualandi,
 *               via Ferrata, 5, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#include <random>

#include "KWD_Histogram2D.h"

int main(int argc, char *argv[]) {
  int n = 128;

  if (argc > 1)
    n = atoi(argv[1]);

  int seed = 13;

  std::random_device
      rd; // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(seed); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> Uniform01(0, 1);
  std::uniform_int_distribution<> Uniform0N(0, n);

  if (true) {
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
    solver.setParam(KWD_METHOD, KWD_EXACT);
    // solver.setParam(KWD_VERBOSITY, KWD_DEBUG);
    KWD::Histogram2D A(Xs.size(), &Xs[0], &Ys[0], &W1[0]);
    KWD::Histogram2D B(Xs.size(), &Xs[0], &Ys[0], &W2[0]);

    double d = solver.dense(A, B);
    double r = solver.runtime();
    auto s = solver.status();
    uint64_t its = solver.iterations();

    PRINT("Dense => %d: fobj: %.6f, time: %.2f, status: %s, iter: %ld, arcs: "
          "%ld, nodes: %ld\n",
          n, d, r, s.c_str(), its, solver.num_arcs(), solver.num_nodes());

    d = solver.distance(A, B, 3);
    r = solver.runtime();
    s = solver.status();
    its = solver.iterations();

    PRINT("Dista => %d: fobj: %.6f, time: %.2f, status: %s, iter: %ld, arcs: "
          "%ld, nodes: %ld\n",
          n, d, r, s.c_str(), its, solver.num_arcs(), solver.num_nodes());

    d = solver.column_generation(A, B, 3);
    r = solver.runtime();
    s = solver.status();
    its = solver.iterations();

    PRINT("ColGe => %d: fobj: %.6f, time: %.2f, status: %s, iter: %ld, arcs: "
          "%ld, nodes: %ld\n",
          n, d, r, s.c_str(), its, solver.num_arcs(), solver.num_nodes());

    for (auto algo : {KWD_BIPARTITE, KWD_MINCOSTFLOW, KWD_COLGEN}) {
      solver.setParam(KWD_ALGORITHM, algo);
      double dist =
          solver.compareExact(Xs.size(), &Xs[0], &Ys[0], &W1[0], &W2[0]);

      PRINT("Exact => %d: fobj: %.6f, time: %.2f, status: %s, iter: %ld, arcs: "
            "%ld, nodes: %ld\n",
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