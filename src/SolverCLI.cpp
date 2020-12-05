/*
 * @fileoverview Copyright (c) 2019-2020, Stefano Gualandi,
 *               via Ferrata, 5, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#include "KWD_Histogram2D.h"

#include <random>


int main(int argc, char* argv[]) {
  int n = 32;

  if (argc > 1) 
     n = atoi(argv[1]);

  int algo = 0;
  if (argc > 2) 
     algo = atoi(argv[2]);

  int seed = 13;

  std::random_device rd;   // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(seed);  // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> Uniform01(0, 1);
  std::uniform_int_distribution<> Uniform0N(0, n);

  int samples = 1000;
  KWD::Histogram2D a;
  KWD::Histogram2D b;
  
  for (size_t i = 0; i < samples; i++) {
     a.add(-100 + Uniform0N(gen), -10 + Uniform0N(gen), Uniform01(gen));
     b.add(-20 + Uniform0N(gen), -30 + Uniform0N(gen), Uniform01(gen));
  }

  a.normalize();
  b.normalize();

  KWD::Solver solver(0);

  double dist = solver.dense(a, b);
  fprintf(stdout, "%d: %.6f %.3f ms\n", n, dist, solver.runtime());

  for (int L = 1; L <= 5; ++L) {
     dist = solver.distance(a, b, L);
     fprintf(stdout, "%d: %d %.6f %.3f ms\n", n, L, dist, solver.runtime());
  }

  return EXIT_SUCCESS;
}