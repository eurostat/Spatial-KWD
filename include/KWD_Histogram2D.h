/**
 * @fileoverview Copyright (c) 2019-2020, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#include <vector>
using std::vector;

#include <unordered_map>
using std::unordered_map;

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifndef __GCD_
#define __GCD_
int GCD(int _a, int _b) {
  int a = (_a >= 0 ? _a : -_a);
  int b = (_b >= 0 ? _b : -_b);
  while (b != 0) {
    int t = b;
    b = a % b;
    a = t;
  }
  return a;
}
#endif

// My Network Simplex
#include "KWD_NetSimplex.h"

struct coprimes_t {
public:
  coprimes_t(int _v, int _w, double _c) : v(_v), w(_w), c_vw(_c) {}
  int v;
  int w;
  double c_vw;
};

typedef std::pair<int, int> int_pair;

struct pair_hash {
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2> &pair) const {
    return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
  }
};

namespace std {
template <> struct hash<std::pair<int, int>> {
  inline size_t operator()(const std::pair<int, int> &v) const {
    std::hash<int> int_hasher;
    return int_hasher(v.first) ^ int_hasher(v.second);
  }
};

} // namespace std

typedef std::unordered_map<int_pair, double, pair_hash> int_pair_dict;

namespace KWD {

class Histogram2D {
public:
  // Standard c'tor
  Histogram2D() {}

  // Second c'tor
  Histogram2D(size_t n, int *X, int *Y, double *W) {
    for (size_t i = 0; i < n; i++)
      update(X[i], Y[i], W[i]);

    normalize();
  }

  // Add a new point
  void add(int _x, int _y, double _w) { Ws[std::make_pair(_x, _y)] = _w; }

  // Add or update a new point
  void update(int _x, int _y, double _w) {
    auto p = std::make_pair(_x, _y);
    auto it = Ws.find(p);
    if (it == Ws.end()) {
      Ws[std::make_pair(_x, _y)] = _w;
    } else {
      Ws[std::make_pair(_x, _y)] = _w + it->second;
    }
  }

  // Getters
  size_t size() const { return Ws.size(); }

  // Total Weigth
  double balance() {
    double t = 0;
    for (const auto &k : Ws)
      t += k.second;
    return t;
  }

  // Make all the weights sum up to 1
  void normalize() {
    double t = balance();
    for (auto &k : Ws)
      k.second = k.second / t;
  }

  // Support for loops
  std::unordered_map<int_pair, double, pair_hash>::iterator begin() {
    return Ws.begin();
  }
  std::unordered_map<int_pair, double, pair_hash>::const_iterator
  begin() const {
    return Ws.begin();
  }

  std::unordered_map<int_pair, double, pair_hash>::iterator end() {
    return Ws.end();
  }
  std::unordered_map<int_pair, double, pair_hash>::const_iterator end() const {
    return Ws.end();
  }

private:
  int_pair_dict Ws;
};

class PointCloud2D {
public:
  void remove(size_t i) {
    std::swap(X[i], X.back());
    std::swap(Y[i], Y.back());
    std::swap(B[i], B.back());
    X.resize(X.size() - 1);
    Y.resize(Y.size() - 1);
    B.resize(B.size() - 1);
  }

  void remove(int x, int y) {
    auto p = std::make_pair(x, y);
    if (M.find(p) != M.end()) {
      size_t i = M.at(p);
      std::swap(X[i], X.back());
      std::swap(Y[i], Y.back());
      std::swap(B[i], B.back());
      X.pop_back();
      Y.pop_back();
      B.pop_back();
      M.erase(p);
    }
  }

  void reserve(size_t t) {
    X.reserve(t);
    Y.reserve(t);
    B.reserve(t);
  }

  void pop_back() { remove(X.back(), Y.back()); }

  void shrink_to_fit() {
    X.shrink_to_fit();
    Y.shrink_to_fit();
    B.shrink_to_fit();
  }

  void add(int x, int y, double b = 0.0) {
    auto p = std::make_pair(x, y);
    if (M.find(p) == M.end()) {
      M[p] = X.size();
      X.push_back(x);
      Y.push_back(y);
      B.push_back(b);
    }
  }

  void update(int x, int y, double b = 0.0) {
    auto p = std::make_pair(x, y);
    auto k = M.find(p);
    if (k == M.end()) {
      M[p] = X.size();
      X.push_back(x);
      Y.push_back(y);
      B.push_back(-b);
    } else {
      B[k->second] = B[k->second] - b;
    }
  }

  void setX(size_t i, int x) { X[i] = x; }
  void setY(size_t i, int y) { Y[i] = y; }
  void setB(size_t i, double b) { B[i] = b; }

  // Merge all the points contained in "other" into this object.
  // The node balance are taken from the "other" object.
  void merge(const PointCloud2D &other) {
    std::unordered_map<std::pair<int, int>, size_t> O = other.getM();

    for (const auto &p : O) {
      size_t j = p.second;
      if (M.find(p.first) != M.end()) {
        size_t i = M.at(p.first);
        B[i] = other.getB(j);
      } else {
        throw std::runtime_error("ERROR 302: point missing");
      }
    }
  }

  void append(const PointCloud2D &other) {
    for (size_t i = 0, i_max = other.size(); i < i_max; ++i)
      add(other.getX(i), other.getY(i), other.getB(i));
  }

  double balance() {
    double t = 0;
    for (size_t i = 0, i_max = B.size(); i < i_max; ++i)
      t += B[i];
    return t;
  }

  bool empty() const { return X.empty(); }

  int getX(size_t i) const { return X[i]; }
  int getY(size_t i) const { return Y[i]; }
  double getB(size_t i) const { return B[i]; }

  size_t size(void) const { return X.size(); }

  // TODO: THIS IS NOT DEFINED RCPP !
  void dump(const std::string &msg = "") const {
    if (!msg.empty())
      fprintf(stdout, "%s\n", msg.c_str());
    for (size_t i = 0, i_max = X.size(); i < i_max; ++i)
      fprintf(stdout, "(%d, %d, %f)\n", X[i], Y[i], B[i]);
    fprintf(stdout, "\n");
    fflush(stdout);
  }

  const std::unordered_map<std::pair<int, int>, size_t> &getM() const {
    return M;
  }

private:
  // Point coordinates (integers)
  std::vector<int> X;
  std::vector<int> Y;
  // Pair to indices
  std::unordered_map<std::pair<int, int>, size_t> M;
  // Node balance
  std::vector<double> B;
};

// Class for computing the convex hull
class ConvexHull {
public:
  // Compute polar among between two points
  double PolarAngle(int ax, int ay, int bx = -1, int by = -1) const {
    int cx = bx, cy = by;
    if (bx == -1) {
      cx = anchor_x;
      cy = anchor_y;
    }
    int x_span = ax - cx;
    int y_span = ay - cy;
    return atan2(y_span, x_span);
  }

  // Square Euclidean distance
  int Distance(int ax, int ay, int bx = -1, int by = -1) const {
    int cx = bx, cy = by;
    if (bx == -1) {
      cx = anchor_x;
      cy = anchor_y;
    }
    int x_span = ax - cx;
    int y_span = ay - cy;
    return (y_span >> 2) + (x_span >> 2);
  }

  // Determinant to detect direction
  int Det(int ax, int ay, int bx, int by, int cx, int cy) const {
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
  }

  PointCloud2D PolarQuickSort(PointCloud2D &Ls) {
    if (Ls.size() <= 1)
      return Ls;

    PointCloud2D smaller, equal, larger;

    double pivot_ang = PolarAngle(Ls.getX(0), Ls.getY(0));

    for (size_t i = 0, i_max = Ls.size(); i < i_max; ++i) {
      double p_ang = PolarAngle(Ls.getX(i), Ls.getY(i));
      if (p_ang < pivot_ang) {
        smaller.add(Ls.getX(i), Ls.getY(i));
      } else {
        if (p_ang == pivot_ang)
          equal.add(Ls.getX(i), Ls.getY(i));
        else
          larger.add(Ls.getX(i), Ls.getY(i));
      }
    }

    auto l1 = PolarQuickSort(smaller);
    while (!equal.empty()) {
      size_t min_idx = 0;
      for (size_t i = 0, i_max = equal.size(); i < i_max; ++i)
        if (Distance(equal.getX(i), equal.getY(i)))
          min_idx = i;
      l1.add(equal.getX(min_idx), equal.getY(min_idx));
      equal.remove(min_idx);
    }
    auto l3 = PolarQuickSort(larger);
    for (size_t i = 0, i_max = l3.size(); i < i_max; ++i)
      l1.add(l3.getX(i), l3.getY(i));
    return l1;
  }

  // Filter the main point along a given direction
  PointCloud2D FilterAxis(const PointCloud2D &Ps) {
    // First filter along an axis
    std::unordered_map<int, std::vector<int>> Xs;
    for (size_t i = 0, i_max = Ps.size(); i < i_max; ++i) {
      int key = Ps.getY(i);
      if (Xs.find(key) == Xs.end()) {
        std::vector<int> tmp;
        Xs[key] = tmp;
      }
      Xs.at(key).push_back(Ps.getX(i));
    }

    PointCloud2D Bs;
    for (auto &k : Xs) {
      auto vet = k.second;
      std::sort(vet.begin(), vet.end());
      Bs.add(vet.front(), k.first);
      if (vet.size() > 1)
        Bs.add(vet.back(), k.first);
    }

    // Then, filter according to the second axis
    std::unordered_map<int, std::vector<int>> Ys;
    for (size_t i = 0, i_max = Bs.size(); i < i_max; ++i) {
      int key = Bs.getX(i);
      if (Ys.find(key) == Ys.end()) {
        std::vector<int> tmp;
        Ys[key] = tmp;
      }
      Ys.at(key).push_back(Bs.getY(i));
    }

    PointCloud2D Rs;
    for (auto &k : Ys) {
      auto vet = k.second;
      std::sort(vet.begin(), vet.end());
      Rs.add(k.first, vet.front());
      if (vet.size() > 1)
        Rs.add(k.first, vet.back());
    }

    return Rs;
  }

  // Find convex hull of given set of points
  PointCloud2D find(const PointCloud2D &Ps) {
    // Preprocessing
    auto Cs = FilterAxis(Ps);
    // Find anchor point
    constexpr size_t null_idx = std::numeric_limits<size_t>::max();
    size_t min_idx = null_idx;
    for (size_t i = 0, i_max = Cs.size(); i < i_max; ++i) {
      if (min_idx == null_idx || Cs.getY(i) < Cs.getY(min_idx))
        min_idx = i;
      if (Cs.getY(i) == Cs.getY(min_idx) && Cs.getX(i) < Cs.getX(min_idx))
        min_idx = i;
    }

    anchor_x = Cs.getX(min_idx);
    anchor_y = Cs.getY(min_idx);
    Cs.remove(anchor_x, anchor_y);

    PointCloud2D Ss = PolarQuickSort(Cs);

    PointCloud2D Hull;
    Hull.add(anchor_x, anchor_y);
    Hull.add(Ss.getX(0), Ss.getY(0));
    for (size_t i = 1, i_max = Ss.size(); i < i_max; ++i) {
      size_t cur = Hull.size();
      while (Det(Hull.getX(cur - 2), Hull.getY(cur - 2), Hull.getX(cur - 1),
                 Hull.getY(cur - 1), Ss.getX(i), Ss.getY(i)) <= 0) {
        Hull.pop_back();
        if (Hull.size() < 2)
          break;
        cur = Hull.size();
      }
      Hull.add(Ss.getX(i), Ss.getY(i));
    }

    // Add all points along the convex hull
    PointCloud2D Rs;
    size_t n = Hull.size();
    for (size_t i = 0, i_max = n - 1; i < i_max; i++) {
      auto tmp = WalkGrid(Hull.getX(i), Hull.getY(i), Hull.getX(i + 1),
                          Hull.getY(i + 1));
      Rs.append(tmp);
    }
    auto tmp = WalkGrid(Hull.getX(n - 1), Hull.getY(n - 1), Hull.getX(0),
                        Hull.getY(0));
    Rs.append(tmp);

    return Rs;
  }

  // Find all the points connecting two dots
  PointCloud2D WalkGrid(int ax, int ay, int bx, int by) {
    int dx = bx - ax;
    int dy = by - ay;

    if (dx == 0) {
      PointCloud2D ps;
      for (int i = std::min(ay, by), i_max = std::max(ay, by); i <= i_max; ++i)
        ps.add(ax, i);
      return ps;
    }

    if (dy == 0) {
      PointCloud2D ps;
      for (int i = std::min(ax, bx), i_max = std::max(ax, bx); i <= i_max; ++i)
        ps.add(i, ay);
      return ps;
    }

    int nx = abs(dx);
    int ny = abs(dy);

    int sign_x = (dx > 0 ? 1 : -1);
    int sign_y = (dy > 0 ? 1 : -1);

    int px = ax;
    int py = ay;
    PointCloud2D points;
    points.add(px, py);

    int ix = 0, iy = 0;
    while (ix < nx || iy < ny) {
      if ((0.5 + ix) / nx < (0.5 + iy) / ny) {
        px += sign_x;
        ix += 1;
      } else {
        py += sign_y;
        iy += 1;
      }
      points.add(px, py);
    }

    return points;
  }

  // Find all the point in the interior of the convex hull
  PointCloud2D FillHull(const PointCloud2D &Ps) const {
    int x_max = std::numeric_limits<int>::min();
    for (size_t i = 0, i_max = Ps.size(); i < i_max; ++i)
      x_max = std::max(x_max, Ps.getX(i));

    std::vector<int> Xmin(1 + x_max, std::numeric_limits<int>::max());
    std::vector<int> Xmax(1 + x_max, std::numeric_limits<int>::min());

    for (size_t i = 0, i_max = Ps.size(); i < i_max; ++i) {
      int x = Ps.getX(i);
      Xmin[x] = std::min(Xmin[x], Ps.getY(i));
      Xmax[x] = std::max(Xmax[x], Ps.getY(i));
    }

    for (int i = 0; i < x_max + 1; ++i)
      if (Xmax[i] == -1)
        throw std::runtime_error("ERROR 201: convex hull issue, Xmin[" +
                                 std::to_string(i) +
                                 "]=" + std::to_string(Xmin[i]));

    PointCloud2D Rs;
    for (int x = 0; x < x_max + 1; ++x)
      for (int y = Xmin[x], y_max = Xmax[x] + 1; y < y_max; ++y)
        Rs.add(x, y);

    return Rs;
  }

private:
  int anchor_x, anchor_y;
};

class Solver {
public:
  // Standard c'tor
  Solver(int n_log = 0) : _runtime(0.0), _n_log(n_log), L(-1) {}

  // Return runtime in milliseconds
  double runtime() const { return _runtime; }

  // Compute KWD distance between A and B with bipartite graph
  double dense(const Histogram2D &A, const Histogram2D &B) {
    // Compute xmin, xmax, ymin, ymax for each axis
    int xmax = std::numeric_limits<int>::min();
    int ymax = std::numeric_limits<int>::min();
    int xmin = std::numeric_limits<int>::max();
    int ymin = std::numeric_limits<int>::max();
    for (const auto &p : A) {
      xmax = std::max(xmax, p.first.first);
      ymax = std::max(ymax, p.first.second);
      xmin = std::min(xmin, p.first.first);
      ymin = std::min(ymin, p.first.second);
    }
    for (const auto &p : B) {
      xmax = std::max(xmax, p.first.first);
      ymax = std::max(ymax, p.first.second);
      xmin = std::min(xmin, p.first.first);
      ymin = std::min(ymin, p.first.second);
    }

    int l = xmax - xmin + 1;
    int w = ymax - ymin + 1;

    // Binary vector for positions
    auto ID = [&w](int x, int y) { return x * w + y; };

    std::vector<int> Ha(size_t(l * w), 0);
    std::vector<int> Hb(size_t(l * w), 0);

    {
      int a = 0;
      for (const auto &p : A) {
        Ha[ID(p.first.first - xmin, p.first.second - ymin)] = a++;
      }

      for (const auto &p : B) {
        Hb[ID(p.first.first - xmin, p.first.second - ymin)] = a++;
      }
    }

    // Network Simplex
    typedef double FlowType;
    typedef double CostType;

    // Build the graph for min cost flow
    NetSimplex<FlowType, CostType> simplex(
        'F', static_cast<int>(A.size() + B.size()),
        static_cast<int>(A.size() * B.size()), _n_log);

    // add first d source nodes
    {
      int i = 0;
      for (const auto &p : A)
        simplex.addNode(i++, p.second);
      for (const auto &p : B)
        simplex.addNode(i++, -p.second);
    }

    for (const auto &p : A) {
      for (const auto &q : B) {
        int v = p.first.first - q.first.first;
        int w = p.first.second - q.first.second;

        simplex.addArc(Ha[ID(p.first.first - xmin, p.first.second - ymin)],
                       Hb[ID(q.first.first - xmin, q.first.second - ymin)],
                       sqrt(pow(v, 2) + pow(w, 2)));
      }
    }

    // Solve the problem to compute the distance
    ProblemType status = simplex.run();

    _runtime = simplex.runtime();

    double distance = std::numeric_limits<CostType>::max();
    if (status != ProblemType::INFEASIBLE && status != ProblemType::UNBOUNDED)
      distance = simplex.totalCost();

    return distance;
  }

  // Compute KWD distance between A and B
  double distance(const Histogram2D &A, const Histogram2D &B, int _L = 3) {
    if (L != _L)
      init_coprimes(_L);

    PointCloud2D ps = mergeHistograms(A, B);

    // Compute convex hull
    ConvexHull ch;
    PointCloud2D As = ch.find(ps);
    PointCloud2D Rs = ch.FillHull(As);

    Rs.merge(ps);

    size_t n = Rs.size();

    // Compute xmax, ymax for each axis
    int xmax = std::numeric_limits<int>::min();
    int ymax = std::numeric_limits<int>::min();
    for (size_t i = 0; i < n; ++i) {
      xmax = std::max(xmax, Rs.getX(i));
      ymax = std::max(ymax, Rs.getY(i));
    }

    // Binary vector for positions
    auto ID = [&ymax](int x, int y) { return x * (ymax + 1) + y; };

    std::vector<bool> M((xmax + 1) * (ymax + 1), false);
    for (size_t i = 0; i < n; ++i)
      M[ID(Rs.getX(i), Rs.getY(i))] = true;

    std::vector<size_t> H((xmax + 1) * (ymax + 1), 0);
    for (size_t i = 0; i < n; ++i)
      H[ID(Rs.getX(i), Rs.getY(i))] = i;

    typedef double FlowType;
    typedef double CostType;

    // Build the graph for min cost flow
    NetSimplex<FlowType, CostType> simplex(
        'F', static_cast<int>(n), static_cast<int>(n * coprimes.size()),
        _n_log);

    // add first d source nodes
    for (size_t i = 0; i < n; ++i)
      simplex.addNode(int(i), Rs.getB(i));

    int m = 0;
    for (size_t h = 0; h < n; ++h) {
      int i = Rs.getX(h);
      int j = Rs.getY(h);
      for (const auto &p : coprimes) {
        int v = p.v;
        int w = p.w;
        if (i + v >= 0 && i + v <= xmax && j + w >= 0 && j + w <= ymax &&
            M[ID(i + v, j + w)]) {
          size_t ff = H[ID(i + v, j + w)];
          simplex.addArc(static_cast<int>(h), static_cast<int>(ff), p.c_vw);
          m++;
        }
      }
    }

    // Solve the problem to compute the distance
    ProblemType status = simplex.run();

    _runtime = simplex.runtime();

    double distance = std::numeric_limits<CostType>::max();
    if (status != ProblemType::INFEASIBLE && status != ProblemType::UNBOUNDED)
      distance = simplex.totalCost();

    return distance;
  }

  // Compute Kantorovich-Wasserstein distance between two measures
  double column_generation(const Histogram2D &A, const Histogram2D &B, int _L) {
    auto start_t = std::chrono::steady_clock::now();
    double _all_p;

    if (L != _L)
      init_coprimes(_L);

    PointCloud2D ps = mergeHistograms(A, B);

    // Compute convex hull
    ConvexHull ch;
    PointCloud2D As = ch.find(ps);
    fprintf(stdout, "find %d\n", As.size());

    PointCloud2D Rs = ch.FillHull(As);
    fprintf(stdout, "FillHull %d\n", As.size());

    Rs.merge(ps);

    size_t n = Rs.size();

    fprintf(stdout, "Convex hull ready %d\n", n);
    fflush(stdout);

    // Compute xmax, ymax for each axis
    int xmax = Rs.getX(0);
    int ymax = Rs.getY(0);
    for (size_t i = 0; i < n; ++i) {
      xmax = std::max(xmax, Rs.getX(i));
      ymax = std::max(ymax, Rs.getY(i));
    }

    fprintf(stdout, "xmax=%d, ymax=%d\n", xmax, ymax);
    fflush(stdout);

    // Binary vector for positions
    auto ID = [&ymax](int x, int y) { return x * (ymax + 1) + y; };

    std::vector<bool> M((xmax + 1) * (ymax + 1), false);
    for (size_t i = 0; i < n; ++i)
      M[ID(Rs.getX(i), Rs.getY(i))] = true;

    std::vector<size_t> H((xmax + 1) * (ymax + 1), 0);
    for (size_t i = 0; i < n; ++i)
      H[ID(Rs.getX(i), Rs.getY(i))] = i;

    typedef double FlowType;
    typedef double CostType;

    // Build the graph for min cost flow
    NetSimplex<FlowType, CostType> simplex('E', static_cast<int>(n), 0, 0);

    // add first d source nodes
    for (size_t i = 0; i < n; ++i)
      simplex.addNode(int(i), Rs.getB(i));

    int it = 0;
    int n_cuts = 0;
    CostType fobj = 0;

    int cmp_tot = 0;

    vector<double> PI(n, 0);

    Vars vars(n);
    for (int i = 0; i < n; ++i)
      vars[i].a = i;

    Vars vnew;
    vnew.reserve(n);

    // Init the simplex
    ProblemType status = simplex.run();

    // Start separation
    while (true) {
      ProblemType status = simplex.reRun();

      // Take the dual values
      double umin = std::numeric_limits<double>::infinity();
      for (int j = 0; j < n; ++j) {
        PI[j] = -simplex.potential(j);
        umin = std::min<double>(umin, PI[j]);
      }

      // Solve separation problem:
      //#pragma omp parallel
      auto start_tt = std::chrono::steady_clock::now();

      {
        //#pragma omp for schedule(dynamic, 8)
        for (size_t h = 0; h < n; ++h) {
          int a = Rs.getX(h);
          int b = Rs.getY(h);

          double best_v = -FEASIBILITY_TOL;
          double best_c = -1;
          int best_j = 0;

          for (const auto &p : coprimes) {
            int v = p.v;
            int w = p.w;
            double c_ij = p.c_vw; // TODO: Memomize!
            if (a + v >= 0 && a + v <= xmax && b + w >= 0 && b + w <= ymax &&
                M[ID(a + v, b + w)]) {
              size_t j = H[ID(a + v, b + w)];

              double violation = c_ij - PI[h] + PI[j];
              if (violation < best_v) {
                best_v = violation;
                best_c = c_ij;
                best_j = j;
              }
            }
          }
          // Store most violated cuts for element i
          vars[h].b = best_j;
          vars[h].c = best_c;
        }
      }
      auto end_tt = std::chrono::steady_clock::now();
      _all_p += double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_tt - start_tt)
                           .count()) /
                1000;

      // Take all negative reduced cost variables
      vnew.clear();
      for (auto &v : vars) {
        if (v.c > -1)
          vnew.push_back(v);
        v.c = -1;
      }

      if (vnew.empty())
        break;

      std::sort(vnew.begin(), vnew.end(),
                [](const Var &v, const Var &w) { return v.c > w.c; });

      // Replace old constraints with new ones
      int new_arcs = simplex.updateArcs(vnew);

      n_cuts += new_arcs;

      if (it % 100 == 0) {
        fobj = simplex.totalCost();
        fprintf(stdout, "it: %d, fobj: %f\n", it, fobj);
      }
      ++it;
    }

    _runtime = simplex.runtime();

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;
    fobj = simplex.totalCost();

    fprintf(stdout, "it: %d, all: %f, simplex: %f, all_p: %f\n", it, _all,
            _runtime, _all_p);

    return fobj;
  }

  void init_coprimes(int L) {
    coprimes.clear();
    for (int v = -L; v <= L; ++v)
      for (int w = -L; w <= L; ++w)
        if (!(v == 0 && w == 0) && GCD(v, w) == 1)
          coprimes.emplace_back(v, w, sqrt(pow(v, 2) + pow(w, 2)));
    coprimes.shrink_to_fit();
  }

  // Parse data from file, with format: i j b1 b1
  PointCloud2D parse(const std::string &filename, char sep = ' ', int off = 0) {
    std::ifstream in_file(filename);

    if (!in_file)
      throw std::runtime_error("FATAL ERROR: Cannot open file");

    PointCloud2D Rs;
    std::vector<double> Bs;
    std::string line;

    // Read first line
    double tot_a = 0;
    double tot_b = 0;
    while (std::getline(in_file, line)) {
      std::stringstream lineStream(line);
      std::string cell;

      std::getline(lineStream, cell, sep);
      int x = std::stoi(cell);
      std::getline(lineStream, cell, sep);
      int y = std::stoi(cell);
      std::getline(lineStream, cell, sep);
      double a = std::stof(cell);
      std::getline(lineStream, cell, sep);
      double b = std::stof(cell);

      //      if (fabs(a - b) > 1e-10) {
      tot_a += a;
      tot_b += b;
      // Check if grid start in position 1 or 0 with parameter "off"
      Rs.add(x - off, y - off, a);
      Bs.emplace_back(b);
      //    }
    };

    // Release resource as soon as possible
    in_file.close();
    // Use as few memory as possible
    Rs.shrink_to_fit();
    // normalize data (rescaling)
    if (Rs.size() != Bs.size())
      throw std::runtime_error(
          "ERROR 301: error in parsing an input file - PointCloud2D");

    double tot = 0;
    for (size_t i = 0, i_max = Bs.size(); i < i_max; ++i) {
      tot += Rs.getB(i) / tot_a - Bs[i] / tot_b;
      Rs.setB(i, Rs.getB(i) / tot_a - Bs[i] / tot_b);
    }

    return Rs;
  }

private:
  // Runtime in milliseconds
  double _runtime;
  // Interval for logging iterations in the simplex algorithm (if =0 no logs at
  // all)
  int _n_log;

  // Merge two historgram into a PointCloud
  PointCloud2D mergeHistograms(const Histogram2D &A, const Histogram2D &B) {
    int xmin = std::numeric_limits<int>::max();
    int ymin = std::numeric_limits<int>::max();
    for (const auto &p : A) {
      xmin = std::min(xmin, p.first.first);
      ymin = std::min(ymin, p.first.second);
    }
    for (const auto &p : B) {
      xmin = std::min(xmin, p.first.first);
      ymin = std::min(ymin, p.first.second);
    }

    PointCloud2D Rs;

    // Read first line
    for (const auto &k : A)
      Rs.add(k.first.first - xmin, k.first.second - ymin, k.second);

    for (const auto &k : B)
      Rs.update(k.first.first - xmin, k.first.second - ymin, k.second);

    // Use as few memory as possible
    Rs.shrink_to_fit();

    return Rs;
  }

  int L;
  // List of pair of coprimes number between (-L, L)
  std::vector<coprimes_t> coprimes;
};

} // end namespace KWD
