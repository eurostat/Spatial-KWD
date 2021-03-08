/**
 * @fileoverview Copyright (c) 2019-2021, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#include <omp.h>

#include <vector>
using std::vector;

#include <array>
using std::array;

#include <unordered_map>
using std::unordered_map;

#include <unordered_set>
using std::unordered_set;

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

void tolower(std::string &data) {
  std::transform(data.begin(), data.end(), data.begin(),
                 [](unsigned char c) { return std::tolower(c); });
}

// List of solver parameters:
//------------------------------

// (exact, approx)
constexpr auto KWD_PAR_METHOD = "Method";
constexpr auto KWD_VAL_EXACT = "exact";
constexpr auto KWD_VAL_APPROX = "approx";

// (bipartite, mincostflow)
constexpr auto KWD_PAR_MODEL = "Model";
constexpr auto KWD_VAL_BIPARTITE = "bipartite";
constexpr auto KWD_VAL_MINCOSTFLOW = "mincostflow";

// (fullmodel, colgen)
constexpr auto KWD_PAR_ALGORITHM = "Algorithm";
constexpr auto KWD_VAL_FULLMODEL = "fullmodel";
constexpr auto KWD_VAL_COLGEN = "colgen";

// ('SILENT', 'INFO', 'DEBUG')
constexpr auto KWD_PAR_VERBOSITY = "Verbosity";
constexpr auto KWD_VAL_SILENT = "silent";
constexpr auto KWD_VAL_INFO = "info";
constexpr auto KWD_VAL_DEBUG = "debug";

constexpr auto KWD_PAR_TIMELIMIT = "TimeLimit";
constexpr auto KWD_PAR_OPTTOLERANCE = "OptTolerance";

constexpr auto KWD_PAR_RECODE = "Recode";

constexpr auto KWD_PAR_UNBALANCED = "Unbalanced";
constexpr auto KWD_PAR_UNBALANCED_COST = "UnbalancedCost";
constexpr auto KWD_PAR_CONVEXHULL = "ConvexHull";

constexpr auto KWD_VAL_TRUE = "true";
constexpr auto KWD_VAL_FALSE = "false";

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
typedef std::unordered_map<int_pair, int, pair_hash> intpair2int;
typedef std::unordered_set<int_pair, pair_hash> int_pair_set;

namespace KWD {

class Histogram2D {
public:
  // Standard c'tor
  Histogram2D() {}

  // Second c'tor
  Histogram2D(int n, int *X, int *Y, double *W) {
    for (int i = 0; i < n; i++)
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
      B.push_back(b);
    } else {
      B[k->second] = B[k->second] + b;
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
      PRINT("%s\n", msg.c_str());
    for (size_t i = 0, i_max = X.size(); i < i_max; ++i)
      PRINT("(%d, %d, %f)\n", X[i], Y[i], B[i]);
    PRINT("\n");
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
  Solver()
      : _runtime(0.0), _n_log(0), L(-1), verbosity(KWD_VAL_INFO), recode(""),
        opt_tolerance(1e-06), timelimit(std::numeric_limits<double>::max()),
        unbalanced(false), unbal_cost(std::numeric_limits<double>::max()) {}

  // Setter/getter for parameters
  std::string getStrParam(const std::string &name) const {
    if (name == KWD_PAR_METHOD)
      return method;
    if (name == KWD_PAR_ALGORITHM)
      return algorithm;
    if (name == KWD_PAR_VERBOSITY)
      return verbosity;
    if (name == KWD_PAR_RECODE)
      return recode;
    if (name == KWD_PAR_UNBALANCED)
      return (unbalanced ? KWD_VAL_TRUE : KWD_VAL_FALSE);
    if (name == KWD_PAR_CONVEXHULL)
      return (convex_hull ? KWD_VAL_TRUE : KWD_VAL_FALSE);

    return "ERROR getStrParam: wrong parameter ->" + name;
  }

  double getDblParam(const std::string &name) const {
    if (name == KWD_PAR_TIMELIMIT)
      return timelimit;
    if (name == KWD_PAR_OPTTOLERANCE)
      return opt_tolerance;
    if (name == KWD_PAR_UNBALANCED_COST)
      return unbal_cost;
    return -1;
  }

  void setStrParam(const std::string &name, const std::string &_value) {
    std::string value(_value);
    tolower(value);

    if (name == KWD_PAR_METHOD)
      method = value;

    if (name == KWD_PAR_ALGORITHM)
      algorithm = value;

    if (name == KWD_PAR_VERBOSITY)
      verbosity = value;

    if (name == KWD_PAR_RECODE)
      recode = value;

    if (name == KWD_PAR_UNBALANCED)
      unbalanced = (value == KWD_VAL_TRUE ? true : false);

    if (name == KWD_PAR_CONVEXHULL)
      convex_hull = (value == KWD_VAL_TRUE ? true : false);
  }

  void setDblParam(const std::string &name, double value) {
    if (name == KWD_PAR_TIMELIMIT)
      timelimit = value;

    if (name == KWD_PAR_OPTTOLERANCE)
      opt_tolerance = value;

    if (name == KWD_PAR_UNBALANCED_COST)
      unbal_cost = value;
  }

  void dumpParam() const {
    PRINT("Internal parameters: %s %s %s %s %.3f %f %s %d %.1f, %d\n",
          method.c_str(), model.c_str(), algorithm.c_str(), verbosity.c_str(),
          timelimit, opt_tolerance, recode.c_str(), (int)unbalanced, unbal_cost,
          (int)convex_hull);
  }

  // Return status of the solver
  std::string status() const {
    if (_status == ProblemType::INFEASIBLE)
      return "Infeasible";
    if (_status == ProblemType::OPTIMAL)
      return "Optimal";
    if (_status == ProblemType::UNBOUNDED)
      return "Unbounded";
    if (_status == ProblemType::TIMELIMIT)
      return "TimeLimit";

    return "Undefined";
  }

  // Return runtime in milliseconds
  double runtime() const { return _runtime; }

  // Number of total iterations of simplex algorithms
  uint64_t iterations() const { return _iterations; }

  // Number of arcs in the model
  uint64_t num_arcs() const { return _num_arcs; }

  // Number of nodes in the model
  uint64_t num_nodes() const { return _num_nodes; }

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

    std::vector<int> Ha(size_t(l) * size_t(w), 0);
    std::vector<int> Hb(size_t(l) * size_t(w), 0);

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
        static_cast<int>(A.size() * B.size()));

    // Set the parameters
    simplex.setTimelimit(timelimit);
    simplex.setVerbosity(verbosity);
    simplex.setOptTolerance(opt_tolerance);

    // add first d source nodes
    {
      for (const auto &p : A)
        simplex.addNode(Ha[ID(p.first.first - xmin, p.first.second - ymin)],
                        p.second);

      for (const auto &q : B)
        simplex.addNode(Hb[ID(q.first.first - xmin, q.first.second - ymin)],
                        -q.second);
    }

    for (const auto &p : A) {
      for (const auto &q : B) {
        int v = p.first.first - q.first.first;
        int w = p.first.second - q.first.second;

        // fprintf(stdout, "%d %d %d %d %.4f\n",
        //        Ha[ID(p.first.first - xmin, p.first.second - ymin)],
        //        Hb[ID(q.first.first - xmin, q.first.second - ymin)], v, w,
        //        sqrt(pow(v, 2) + pow(w, 2)));

        simplex.addArc(Ha[ID(p.first.first - xmin, p.first.second - ymin)],
                       Hb[ID(q.first.first - xmin, q.first.second - ymin)],
                       sqrt(pow(v, 2) + pow(w, 2)));
      }
    }

    // Solve the problem to compute the distance
    if (verbosity == KWD_VAL_INFO)
      PRINT("INFO: running NetSimplex with V=%ld and E=%ld\n",
            simplex.num_nodes(), simplex.num_arcs());

    _status = simplex.run();

    _runtime = simplex.runtime();
    _iterations = simplex.iterations();
    _num_arcs = simplex.num_arcs();
    _num_nodes = simplex.num_nodes();

    double distance = std::numeric_limits<CostType>::max();
    if (_status != ProblemType::INFEASIBLE &&
        _status != ProblemType::UNBOUNDED && _status != ProblemType::TIMELIMIT)
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
    PointCloud2D As, Rs;
    if (convex_hull) {
      As = ch.find(ps);
      Rs = ch.FillHull(As);
    } else {
      Rs = ch.FillHull(ps);
    }

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
        'F', static_cast<int>(n), static_cast<int>(n * coprimes.size()));

    // Set the parameters
    simplex.setTimelimit(timelimit);
    simplex.setVerbosity(verbosity);
    simplex.setOptTolerance(opt_tolerance);

    // add first d source nodes
    for (size_t i = 0; i < n; ++i) {
      simplex.addNode(int(i), Rs.getB(i));
    }

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
        }
      }
    }

    // Solve the problem to compute the distance
    if (verbosity == KWD_VAL_INFO)
      PRINT("INFO: running NetSimplex with V=%ld and E=%ld\n",
            simplex.num_nodes(), simplex.num_arcs());

    _status = simplex.run();

    _runtime = simplex.runtime();
    _iterations = simplex.iterations();
    _num_arcs = simplex.num_arcs();
    _num_nodes = simplex.num_nodes();

    double distance = std::numeric_limits<CostType>::max();
    if (_status != ProblemType::INFEASIBLE &&
        _status != ProblemType::UNBOUNDED && _status != ProblemType::TIMELIMIT)
      distance = simplex.totalCost();

    return distance;
  }

  // Compute Kantorovich-Wasserstein distance between two measures
  double column_generation(const Histogram2D &A, const Histogram2D &B, int _L) {
    auto start_t = std::chrono::steady_clock::now();
    double _all_p = 0.0;

    if (L != _L)
      init_coprimes(_L);

    PointCloud2D ps = mergeHistograms(A, B);

    // Compute convex hull
    ConvexHull ch;
    PointCloud2D As, Rs;
    if (convex_hull) {
      As = ch.find(ps);
      Rs = ch.FillHull(As);
    } else {
      Rs = ch.FillHull(ps);
    }

    Rs.merge(ps);

    size_t n = Rs.size();

    // Compute xmax, ymax for each axis
    int xmax = Rs.getX(0);
    int ymax = Rs.getY(0);
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
    NetSimplex<FlowType, CostType> simplex('E', static_cast<int>(n), 0);

    // Set the parameters
    simplex.setTimelimit(timelimit);
    simplex.setVerbosity(verbosity);
    simplex.setOptTolerance(opt_tolerance);

    // add first d source nodes
    for (size_t i = 0; i < n; ++i)
      simplex.addNode(int(i), Rs.getB(i));

    int it = 0;
    int n_cuts = 0;
    CostType fobj = 0;

    vector<double> pi(n, 0);

    Vars vars(n);
    for (size_t i = 0; i < n; ++i)
      vars[i].a = int(i);

    Vars vnew;
    vnew.reserve(n);

    // Init the simplex
    simplex.run();
    _iterations = simplex.iterations();

    // Start separation
    while (true) {
      _status = simplex.reRun();

      if (_status == ProblemType::TIMELIMIT)
        break;

      // Take the dual values
      for (size_t j = 0; j < n; ++j)
        pi[j] = -simplex.potential(j);

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

              double violation = c_ij - pi[h] + pi[j];
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

      ++it;
    }

    _runtime = simplex.runtime();
    _iterations = simplex.iterations();
    _num_arcs = simplex.num_arcs();
    _num_nodes = simplex.num_nodes();

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    fobj = simplex.totalCost();

    if (_n_log > 0)
      PRINT("it: %d, fobj: %f, all: %f, simplex: %f, all_p: %f\n", it, fobj,
            _all, _runtime, _all_p);

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

  // In case of coordinates with regulare distance, recode the coordinate to
  // conescutive integers
  void recoding(int n, int *Xs) {
    int xmin = Xs[0];
    int xmax = Xs[0];

    std::unordered_set<int> unique;

    for (int i = 0; i < n; i++) {
      xmin = std::min(Xs[i], xmin);
      xmax = std::max(Xs[i], xmax);
      unique.insert(Xs[i]);
    }

    int nx = unique.size();
    int d = (xmax - xmin) / (nx - 1);

    for (int i = 0; i < n; i++) {
      Xs[i] = (Xs[i] - xmin) / d;
    }
  }

  // Check if input has all consecutive coordinates
  bool check_coding(int n, int *Xs) const {
    int xmin = Xs[0];
    int xmax = Xs[0];

    std::unordered_set<int> unique;

    for (int i = 0; i < n; i++) {
      xmin = std::min(Xs[i], xmin);
      xmax = std::max(Xs[i], xmax);
      unique.insert(Xs[i]);
    }

    int nx = unique.size();
    int d = (xmax - xmin) / (nx - 1);

    if (d != 1)
      return true;
    else
      return false;
  }

  // Reindex the input coordinates
  intpair2int reindex(int n, int *Xs, int *Ys) {
    intpair2int XY;

    int idx = 0;
    for (int i = 0; i < n; i++) {
      auto p = int_pair(Xs[i], Ys[i]);
      auto it = XY.find(p);
      if (it == XY.end()) {
        XY[p] = idx;
        idx++;
      }
    }

    return XY;
  }

  // New interface for the solver
  double compareExact(int _n, int *_Xs, int *_Ys, double *_W1, double *_W2) {
    // Check for correct input
    if (check_coding(_n, _Xs))
      PRINT(
          "WARNING: the Xs input coordinates are not consecutives integers.\n");
    if (check_coding(_n, _Ys))
      PRINT(
          "WARNING: the Ys input coordinates are not consecutives integers.\n");

    if (recode != "") {
      PRINT("INFO: Recoding the input coordinates to consecutive integers.\n");
      recoding(_n, _Xs);
      recoding(_n, _Ys);
    }

    if (verbosity == KWD_VAL_INFO)
      dumpParam();

    intpair2int XY = reindex(_n, _Xs, _Ys);
    int n = XY.size();

    vector<int> Xs, Ys;
    vector<double> W1, W2;
    Xs.resize(n, 0);
    Ys.resize(n, 0);
    W1.resize(n, 0);
    W2.resize(n, 0);
    for (int i = 0; i < _n; i++) {
      int idx = XY[int_pair(_Xs[i], _Ys[i])];
      Xs[idx] = _Xs[i];
      Ys[idx] = _Ys[i];
      W1[idx] += _W1[i];
      W2[idx] += _W2[i];

      if (_W1[i] < 0.0) {
        PRINT("WARNING: weight W1[%d]=%.4f is negative. Only positive weights "
              "are allowed.\n",
              i, _W1[i]);
        throw std::runtime_error(
            "FATAL ERROR: Input histogram with negative weigths");
      }
      if (_W2[i] < 0.0) {
        PRINT("WARNING: weight W2[%d]=%.4f is negative. Only positive weights "
              "are allowed.\n",
              i, _W2[i]);
        throw std::runtime_error(
            "FATAL ERROR: Input histogram with negative weigths");
      }
    }

    // Get the largest bounding box
    auto xy = getMinMax(n, &Xs[0], &Ys[0]);

    int xmin = xy[0];
    int ymin = xy[1];
    // int xmax = xy[2];
    // int ymax = xy[3];

    // Rescale all integers coordinates to (0,0)
    double tot_w1 = 0.0;
    double tot_w2 = 0.0;
    for (int i = 0; i < n; ++i) {
      Xs[i] = Xs[i] - xmin;
      Ys[i] = Ys[i] - ymin;

      tot_w1 += W1[i];
      tot_w2 += W2[i];
    }

    for (int i = 0; i < n; ++i) {
      W1[i] = W1[i] / tot_w1;
      W2[i] = W2[i] / tot_w2;
    }

    if (algorithm == KWD_VAL_BIPARTITE) {
      // Network Simplex: Build the bipartite graph
      NetSimplex<double, double> simplex('F', (n + n), (n * n));

      // Set the parameters
      simplex.setTimelimit(timelimit);
      simplex.setVerbosity(verbosity);
      simplex.setOptTolerance(opt_tolerance);

      for (int i = 0; i < n; ++i)
        simplex.addNode(i, W1[i]);

      for (int i = 0; i < n; ++i)
        simplex.addNode(n + i, -W2[i]);

      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
          int v = Xs[i] - Xs[j];
          int w = Ys[i] - Ys[j];

          // fprintf(stdout, "%d %d %d %d %.4f\n", i, j, v, w,
          //        sqrt(pow(v, 2) + pow(w, 2)));

          simplex.addArc(i, n + j, sqrt(pow(v, 2) + pow(w, 2)));
        }

      if (verbosity == KWD_VAL_INFO)
        PRINT("INFO: running NetSimplex with V=%ld and E=%ld\n",
              simplex.num_nodes(), simplex.num_arcs());

      // Solve the problem to compute the distance
      _status = simplex.run();

      _runtime = simplex.runtime();
      _iterations = simplex.iterations();
      _num_arcs = simplex.num_arcs();
      _num_nodes = simplex.num_nodes();

      double distance = std::numeric_limits<double>::max();
      if (_status != ProblemType::INFEASIBLE &&
          _status != ProblemType::UNBOUNDED &&
          _status != ProblemType::TIMELIMIT)

        distance = simplex.totalCost();

      return distance;
    }

    // Second option for algorithm
    if (algorithm == KWD_VAL_MINCOSTFLOW) {
      PointCloud2D ps = mergeHistograms(n, &Xs[0], &Ys[0], &W1[0], &W2[0]);

      // Compute convex hull
      ConvexHull ch;
      PointCloud2D As, Rs;
      if (convex_hull) {
        As = ch.find(ps);
        Rs = ch.FillHull(As);
      } else {
        Rs = ch.FillHull(ps);
      }

      Rs.merge(ps);

      int n = static_cast<int>(Rs.size());

      // Compute xmax, ymax for each axis
      int xmin = std::numeric_limits<int>::max();
      int ymin = std::numeric_limits<int>::max();
      int xmax = std::numeric_limits<int>::min();
      int ymax = std::numeric_limits<int>::min();
      for (int i = 0; i < n; ++i) {
        xmin = std::min(xmin, Rs.getX(i));
        ymin = std::min(ymin, Rs.getY(i));
        xmax = std::max(xmax, Rs.getX(i));
        ymax = std::max(ymax, Rs.getY(i));
      }
      xmax++;
      ymax++;

      init_coprimes(
          std::max(xmax - xmin, ymax - ymin)); // TODO: make it parallel

      // Binary vector for positions
      auto ID = [&ymax](int x, int y) { return x * ymax + y; };

      std::vector<bool> M(size_t(xmax) * size_t(ymax), false);
      for (int i = 0; i < n; ++i)
        M[ID(Rs.getX(i), Rs.getY(i))] = true;

      std::vector<int> H(size_t(xmax) * size_t(ymax), 0);
      for (int i = 0; i < n; ++i)
        H[ID(Rs.getX(i), Rs.getY(i))] = i;

      typedef double FlowType;
      typedef double CostType;

      // Build the graph for min cost flow
      NetSimplex<FlowType, CostType> simplex(
          'F', n, n * static_cast<int>(coprimes.size()));

      // Set the parameters
      simplex.setTimelimit(timelimit);
      simplex.setVerbosity(verbosity);
      simplex.setOptTolerance(opt_tolerance);

      // add first d source nodes
      for (int i = 0; i < n; ++i)
        simplex.addNode(i, Rs.getB(i));

      for (int h = 0; h < n; ++h) {
        int i = Rs.getX(h);
        int j = Rs.getY(h);
        for (const auto &p : coprimes) {
          int v = p.v;
          int w = p.w;
          if (i + v >= xmin && i + v < xmax && j + w >= ymin && j + w < ymax &&
              M[ID(i + v, j + w)]) {
            int ff = H[ID(i + v, j + w)];
            simplex.addArc(h, ff, p.c_vw);
          }
        }
      }

      // Solve the problem to compute the distance
      if (verbosity == KWD_VAL_INFO)
        PRINT("INFO: running NetSimplex with V=%ld and E=%ld\n",
              simplex.num_nodes(), simplex.num_arcs());

      _status = simplex.run();

      _runtime = simplex.runtime();
      _iterations = simplex.iterations();
      _num_arcs = simplex.num_arcs();
      _num_nodes = simplex.num_nodes();

      double distance = std::numeric_limits<CostType>::max();

      if (_status != ProblemType::INFEASIBLE &&
          _status != ProblemType::UNBOUNDED &&
          _status != ProblemType::TIMELIMIT)
        distance = simplex.totalCost();

      return distance;
    }

    // Second option for algorithm
    if (algorithm == KWD_VAL_COLGEN) {
      auto start_t = std::chrono::steady_clock::now();
      double _all_p = 0.0;

      PointCloud2D ps = mergeHistograms(n, &Xs[0], &Ys[0], &W1[0], &W2[0]);

      // Compute convex hull
      ConvexHull ch;
      PointCloud2D As, Rs;
      if (convex_hull) {
        As = ch.find(ps);
        Rs = ch.FillHull(As);
      } else {
        Rs = ch.FillHull(ps);
      }

      Rs.merge(ps);

      int n = Rs.size();

      // Compute xmax, ymax for each axis
      int xmax = std::numeric_limits<int>::min();
      int ymax = std::numeric_limits<int>::min();
      for (int i = 0; i < n; ++i) {
        xmax = std::max(xmax, Rs.getX(i));
        ymax = std::max(ymax, Rs.getY(i));
      }
      // The size is xmax+1, and ymax+1
      xmax++;
      ymax++;

      init_coprimes(std::max(xmax, ymax) - 1); // TODO: make it parallel

      // Binary vector for positions
      auto ID = [&ymax](int x, int y) { return x * ymax + y; };

      std::vector<bool> M(size_t(xmax) * size_t(ymax), false);
      for (int i = 0; i < n; ++i)
        M[ID(Rs.getX(i), Rs.getY(i))] = true;

      std::vector<int> H(size_t(xmax) * size_t(ymax), 0);
      for (int i = 0; i < n; ++i)
        H[ID(Rs.getX(i), Rs.getY(i))] = i;

      typedef double FlowType;
      typedef double CostType;

      // Build the graph for min cost flow
      NetSimplex<FlowType, CostType> simplex('E', n, 0);

      // Set the parameters
      simplex.setTimelimit(timelimit);
      simplex.setVerbosity(verbosity);
      simplex.setOptTolerance(opt_tolerance);

      // add first d source nodes
      for (int i = 0; i < n; ++i)
        simplex.addNode(i, Rs.getB(i));

      int it = 0;
      int n_cuts = 0;
      CostType fobj = 0;
      double negeps = std::nextafter(-opt_tolerance, -0.0);

      vector<double> pi(n, 0);

      Vars vars(n);
      for (int i = 0; i < n; ++i)
        vars[i].a = i;

      Vars vnew;
      vnew.reserve(n);

      // Init the simplex
      simplex.run();

      // Start separation
      while (true) {
        _status = simplex.reRun();
        if (_status == ProblemType::TIMELIMIT)
          break;

        // Take the dual values
        for (int j = 0; j < n; ++j)
          pi[j] = -simplex.potential(j);

        // Solve separation problem:
        auto start_tt = std::chrono::steady_clock::now();
#pragma omp parallel
        {
#pragma omp for schedule(dynamic, 1)
          for (int h = 0; h < n; ++h) {
            int a = Rs.getX(h);
            int b = Rs.getY(h);

            double best_v = negeps;
            double best_c = -1;
            int best_j = 0;

            for (const auto &p : coprimes) {
              int v = p.v;
              int w = p.w;
              double c_ij = p.c_vw;
              if (a + v >= 0 && a + v < xmax && b + w >= 0 && b + w < ymax &&
                  M[ID(a + v, b + w)]) {
                int j = H[ID(a + v, b + w)];

                double violation = c_ij - pi[h] + pi[j];
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

        // Take all negative reduced cost variables
        vnew.clear();
        for (auto &v : vars) {
          if (v.c > -1)
            vnew.push_back(v);
          v.c = -1;
        }

        auto end_tt = std::chrono::steady_clock::now();
        _all_p += double(std::chrono::duration_cast<std::chrono::milliseconds>(
                             end_tt - start_tt)
                             .count()) /
                  1000;

        if (vnew.empty())
          break;

        std::sort(vnew.begin(), vnew.end(),
                  [](const Var &v, const Var &w) { return v.c > w.c; });

        // Replace old constraints with new ones
        int new_arcs = simplex.updateArcs(vnew);

        n_cuts += new_arcs;

        ++it;
      }

      auto end_t = std::chrono::steady_clock::now();
      auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                             end_t - start_t)
                             .count()) /
                  1000;

      _runtime = _all;
      _iterations = simplex.iterations();
      _num_arcs = simplex.num_arcs();
      _num_nodes = simplex.num_nodes();

      fobj = simplex.totalCost();

      if (_n_log > 0)
        PRINT("it: %d, fobj: %f, all: %f, simplex: %f, all_p: %f\n", it, fobj,
              _all, _runtime, _all_p);

      return fobj;
    }

    return -1;
  }

  double compareApprox(int _n, int *_Xs, int *_Ys, double *_W1, double *_W2,
                       int _L) {
    // Check for correct input
    if (check_coding(_n, _Xs))
      PRINT(
          "WARNING: the Xs input coordinates are not consecutives integers.\n");
    if (check_coding(_n, _Ys))
      PRINT(
          "WARNING: the Ys input coordinates are not consecutives integers.\n");

    if (recode != "") {
      PRINT("INFO: Recoding the input coordinates to consecutive integers.\n");
      recoding(_n, _Xs);
      recoding(_n, _Ys);
    }

    if (verbosity == KWD_VAL_INFO)
      dumpParam();
    intpair2int XY = reindex(_n, _Xs, _Ys);
    int n = XY.size();

    vector<int> Xs, Ys;
    vector<double> W1, W2;
    Xs.resize(n, 0);
    Ys.resize(n, 0);
    W1.resize(n, 0);
    W2.resize(n, 0);
    for (int i = 0; i < _n; i++) {
      int idx = XY[int_pair(_Xs[i], _Ys[i])];
      Xs[idx] = _Xs[i];
      Ys[idx] = _Ys[i];
      W1[idx] += _W1[i];
      W2[idx] += _W2[i];
      if (_W1[i] < 0.0) {
        PRINT("WARNING: weight W1[%d]=%.4f is negative. Only positive weights "
              "are allowed.\n",
              i, _W1[i]);
        throw std::runtime_error(
            "FATAL ERROR: Input histogram with negative weigths");
      }
      if (_W2[i] < 0.0) {
        PRINT("WARNING: weight W2[%d]=%.4f is negative. Only positive weights "
              "are allowed.\n",
              i, _W2[i]);
        throw std::runtime_error(
            "FATAL ERROR: Input histogram with negative weigths");
      }
    }

    // Get the largest bounding box
    auto xy = getMinMax(n, &Xs[0], &Ys[0]);

    int xmin = xy[0];
    int ymin = xy[1];
    // int xmax = xy[2];
    // int ymax = xy[3];

    // Rescale all integers coordinates to (0,0)
    double tot_w1 = 0.0;
    double tot_w2 = 0.0;
    for (int i = 0; i < n; ++i) {
      Xs[i] = Xs[i] - xmin;
      Ys[i] = Ys[i] - ymin;

      tot_w1 += W1[i];
      tot_w2 += W2[i];
    }

    // Rebalance the total mass only if it is not an unbalanced probelm
    if (!unbalanced) {
      for (int i = 0; i < n; ++i) {
        W1[i] = W1[i] / tot_w1;
        W2[i] = W2[i] / tot_w2;
      }
    }

    // Second option for algorithm
    if (algorithm == KWD_VAL_FULLMODEL) {
      PointCloud2D ps = mergeHistograms(n, &Xs[0], &Ys[0], &W1[0], &W2[0]);

      // Compute convex hull
      ConvexHull ch;
      PointCloud2D As, Rs;
      if (convex_hull) {
        As = ch.find(ps);
        Rs = ch.FillHull(As);
      } else {
        Rs = ch.FillHull(ps);
      }

      Rs.merge(ps);
      Rs.dump();
      int n = static_cast<int>(Rs.size());

      // Compute xmax, ymax for each axis
      int xmin = std::numeric_limits<int>::max();
      int ymin = std::numeric_limits<int>::max();
      int xmax = std::numeric_limits<int>::min();
      int ymax = std::numeric_limits<int>::min();
      for (int i = 0; i < n; ++i) {
        xmin = std::min(xmin, Rs.getX(i));
        ymin = std::min(ymin, Rs.getY(i));
        xmax = std::max(xmax, Rs.getX(i));
        ymax = std::max(ymax, Rs.getY(i));
      }
      xmax++;
      ymax++;

      if (_L != L) {
        L = _L;
        init_coprimes(_L); // TODO: make it parallel
      }

      // Binary vector for positions
      auto ID = [&ymax](int x, int y) { return x * ymax + y; };

      std::vector<bool> M(size_t(xmax) * size_t(ymax), false);
      for (int i = 0; i < n; ++i)
        M[ID(Rs.getX(i), Rs.getY(i))] = true;

      std::vector<int> H(size_t(xmax) * size_t(ymax), 0);
      for (int i = 0; i < n; ++i)
        H[ID(Rs.getX(i), Rs.getY(i))] = i;

      typedef double FlowType;
      typedef double CostType;

      // Build the graph for min cost flow
      NetSimplex<FlowType, CostType> simplex(
          'F', n + int(unbalanced == true),
          n *static_cast<int>(coprimes.size()));

      // Set the parameters
      simplex.setTimelimit(timelimit);
      simplex.setVerbosity(verbosity);
      simplex.setOptTolerance(opt_tolerance);

      // add first d source nodes
      for (int i = 0; i < n; ++i)
        simplex.addNode(i, Rs.getB(i));

      for (int h = 0; h < n; ++h) {
        int i = Rs.getX(h);
        int j = Rs.getY(h);
        for (const auto &p : coprimes) {
          int v = p.v;
          int w = p.w;
          if (i + v >= xmin && i + v < xmax && j + w >= ymin && j + w < ymax &&
              M[ID(i + v, j + w)]) {
            int ff = H[ID(i + v, j + w)];
            simplex.addArc(h, ff, p.c_vw);
          }
        }
      }

      // Add noded for unbalanced transport, if parater is set
      if (unbalanced) {
        double bb = -Rs.balance();
        double c1 = (bb < 0 ? unbal_cost : 0);
        double c2 = (bb < 0 ? 0 : unbal_cost);
        simplex.addNode(n, bb);

        for (int i = 0; i < n; ++i)
          simplex.addArc(i, n, c1);

        for (int i = 0; i < n; ++i)
          simplex.addArc(n, i, c2);
      }

      // Solve the problem to compute the distance
      if (verbosity == KWD_VAL_INFO)
        PRINT("INFO: running NetSimplex with V=%ld and E=%ld\n",
              simplex.num_nodes(), simplex.num_arcs());
      _status = simplex.run();

      _runtime = simplex.runtime();
      _iterations = simplex.iterations();
      _num_arcs = simplex.num_arcs();
      _num_nodes = simplex.num_nodes();

      double distance = std::numeric_limits<CostType>::max();

      if (_status != ProblemType::INFEASIBLE &&
          _status != ProblemType::UNBOUNDED &&
          _status != ProblemType::TIMELIMIT)
        distance = simplex.totalCost();

      if (unbalanced)
        distance = distance / std::max(tot_w1, tot_w2);

      return distance;
    }

    // Second option for algorithm
    if (algorithm == KWD_VAL_COLGEN) {
      auto start_t = std::chrono::steady_clock::now();
      double _all_p = 0.0;

      PointCloud2D ps = mergeHistograms(n, &Xs[0], &Ys[0], &W1[0], &W2[0]);

      // Compute convex hull
      ConvexHull ch;
      PointCloud2D As, Rs;
      if (convex_hull) {
        As = ch.find(ps);
        Rs = ch.FillHull(As);
      } else {
        Rs = ch.FillHull(ps);
      }

      Rs.merge(ps);

      int n = Rs.size();

      // Compute xmax, ymax for each axis
      int xmax = std::numeric_limits<int>::min();
      int ymax = std::numeric_limits<int>::min();
      for (int i = 0; i < n; ++i) {
        xmax = std::max(xmax, Rs.getX(i));
        ymax = std::max(ymax, Rs.getY(i));
      }
      // The size is xmax+1, and ymax+1
      xmax++;
      ymax++;

      if (_L != L) {
        L = _L;
        init_coprimes(_L); // TODO: make it parallel
      }

      // Binary vector for positions
      auto ID = [&ymax](int x, int y) { return x * ymax + y; };

      std::vector<bool> M(size_t(xmax) * size_t(ymax), false);
      for (int i = 0; i < n; ++i)
        M[ID(Rs.getX(i), Rs.getY(i))] = true;

      std::vector<int> H(size_t(xmax) * size_t(ymax), 0);
      for (int i = 0; i < n; ++i)
        H[ID(Rs.getX(i), Rs.getY(i))] = i;

      typedef double FlowType;
      typedef double CostType;

      // Build the graph for min cost flow
      NetSimplex<FlowType, CostType> simplex('E', n + int(unbalanced == true),
                                             0);

      // Set the parameters
      simplex.setTimelimit(timelimit);
      simplex.setVerbosity(verbosity);
      simplex.setOptTolerance(opt_tolerance);

      // add first d source nodes
      for (int i = 0; i < n; ++i)
        simplex.addNode(i, Rs.getB(i));

      // Add noded for unbalanced transport, if parater is set
      if (unbalanced) {
        double bb = -Rs.balance();
        double c1 = (bb < 0 ? unbal_cost : 0);
        double c2 = (bb < 0 ? 0 : unbal_cost);
        simplex.addNode(n, bb);

        for (int i = 0; i < n; ++i)
          simplex.addArc(i, n, c1);

        for (int i = 0; i < n; ++i)
          simplex.addArc(n, i, c2);
      }

      int it = 0;
      int n_cuts = 0;
      CostType fobj = 0;
      double negeps = std::nextafter(-opt_tolerance, -0.0);

      vector<double> pi(n, 0);

      Vars vars(n);
      for (int i = 0; i < n; ++i)
        vars[i].a = i;

      Vars vnew;
      vnew.reserve(n);

      // Init the simplex
      simplex.run();

      // Start separation
      while (true) {
        _status = simplex.reRun();
        if (_status == ProblemType::TIMELIMIT)
          break;

        // Take the dual values
        for (int j = 0; j < n; ++j)
          pi[j] = -simplex.potential(j);

        // Solve separation problem:
        auto start_tt = std::chrono::steady_clock::now();
#pragma omp parallel
        {
#pragma omp for schedule(dynamic, 1)
          for (int h = 0; h < n; ++h) {
            int a = Rs.getX(h);
            int b = Rs.getY(h);

            double best_v = negeps;
            double best_c = -1;
            int best_j = 0;

            for (const auto &p : coprimes) {
              int v = p.v;
              int w = p.w;
              double c_ij = p.c_vw;
              if (a + v >= 0 && a + v < xmax && b + w >= 0 && b + w < ymax &&
                  M[ID(a + v, b + w)]) {
                int j = H[ID(a + v, b + w)];

                double violation = c_ij - pi[h] + pi[j];
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

        // Take all negative reduced cost variables
        vnew.clear();
        for (auto &v : vars) {
          if (v.c > -1)
            vnew.push_back(v);
          v.c = -1;
        }

        auto end_tt = std::chrono::steady_clock::now();
        _all_p += double(std::chrono::duration_cast<std::chrono::milliseconds>(
                             end_tt - start_tt)
                             .count()) /
                  1000;

        if (vnew.empty())
          break;

        std::sort(vnew.begin(), vnew.end(),
                  [](const Var &v, const Var &w) { return v.c > w.c; });

        // Replace old constraints with new ones
        int new_arcs = simplex.updateArcs(vnew);

        n_cuts += new_arcs;

        ++it;
      }

      auto end_t = std::chrono::steady_clock::now();
      auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                             end_t - start_t)
                             .count()) /
                  1000;

      _runtime = _all;
      _iterations = simplex.iterations();
      _num_arcs = simplex.num_arcs();
      _num_nodes = simplex.num_nodes();

      fobj = simplex.totalCost();
      if (unbalanced)
        fobj = fobj / std::max(tot_w1, tot_w2);

      if (_n_log > 0)
        PRINT("it: %d, fobj: %f, all: %f, simplex: %f, all_p: %f\n", it, fobj,
              _all, _runtime, _all_p);

      return fobj;
    }

    return -1;
  }

  vector<double> compareApprox(int _n, int _m, int *_Xs, int *_Ys, double *_W1,
                               double *_Ws, int _L) {
    // Check for correct input
    if (check_coding(_n, _Xs))
      PRINT(
          "WARNING: the Xs input coordinates are not consecutives integers.\n");
    if (check_coding(_n, _Ys))
      PRINT(
          "WARNING: the Ys input coordinates are not consecutives integers.\n");

    if (recode != "") {
      PRINT("INFO: Recoding the input coordinates to consecutive integers.\n");
      recoding(_n, _Xs);
      recoding(_n, _Ys);
    }

    if (verbosity == KWD_VAL_INFO)
      dumpParam();

    // Get the largest bounding box
    auto xy = getMinMax(_n, _Xs, _Ys);

    int xmin = xy[0];
    int ymin = xy[1];
    // int xmax = xy[2];
    // int ymax = xy[3];

    // Reindex and scale
    intpair2int XY = reindex(_n, _Xs, _Ys);
    int N = XY.size();

    // Binary vector for positions
    // auto WD = [&n](int x, int y) { return y * n + x; };

    vector<int> Xs, Ys;
    vector<double> W1, Ws;
    Xs.resize(N, 0);
    Ys.resize(N, 0);
    W1.resize(N, 0.0);
    Ws.resize(size_t(N) * size_t(_m), 0.0);

    double tot_w1 = 0.0;
    vector<double> tot_ws(_m, 0.0);

    intpair2int MXY;
    for (int i = 0; i < _n; i++) {
      int idx = XY.at(int_pair(_Xs[i], _Ys[i]));

      Xs[idx] = _Xs[i] - xmin;
      Ys[idx] = _Ys[i] - ymin;
      MXY[int_pair(Xs[idx], Ys[idx])] = idx;

      W1[idx] += _W1[i];
      tot_w1 += _W1[i];
      if (_W1[i] < 0.0) {
        PRINT("WARNING: weight W1[%d]=%.4f is negative. Only positive weights "
              "are allowed.\n",
              i, _W1[i]);
        throw std::runtime_error(
            "FATAL ERROR: Input histogram with negative weigths");
      }
      for (int j = 0; j < _m; ++j) {
        Ws[j * N + idx] += _Ws[j * _n + i];
        tot_ws[j] += _Ws[j * _n + i];
        if (_Ws[j * _n + i] < 0.0) {
          PRINT("WARNING: weight Ws[%d,%d]=%.4f is negative. Only positive "
                "weights are allowed.\n",
                i, j, _Ws[j * _n + i]);
          throw std::runtime_error(
              "FATAL ERROR: Input histogram with negative weigths");
        }
      }
    }

    // Rescale all integers coordinates to (0,0)
    if (!unbalanced) {
      for (int i = 0; i < N; ++i) {
        W1[i] = W1[i] / tot_w1;
        for (int j = 0; j < _m; ++j)
          Ws[j * N + i] = Ws[j * N + i] / tot_ws[j];
      }
    }

    // Set the coprimes set
    if (_L != L) {
      L = _L;
      init_coprimes(_L); // TODO: make it parallel
    }

    // Return value
    vector<double> Ds(_m, -1);

    PointCloud2D ps = mergeCoordinates(N, &Xs[0], &Ys[0]);

    // Compute convex hull
    ConvexHull ch;

    PointCloud2D As, Rs;
    if (convex_hull) {
      As = ch.find(ps);
      Rs = ch.FillHull(As);
    } else {
      Rs = ch.FillHull(ps);
    }

    int n = static_cast<int>(Rs.size());

    // Second option for algorithm
    if (algorithm == KWD_VAL_MINCOSTFLOW) {
      // Compute xmax, ymax for each axis
      int xmin = std::numeric_limits<int>::max();
      int ymin = std::numeric_limits<int>::max();
      int xmax = std::numeric_limits<int>::min();
      int ymax = std::numeric_limits<int>::min();
      for (int i = 0; i < n; ++i) {
        xmin = std::min(xmin, Rs.getX(i));
        ymin = std::min(ymin, Rs.getY(i));
        xmax = std::max(xmax, Rs.getX(i));
        ymax = std::max(ymax, Rs.getY(i));
      }
      xmax++;
      ymax++;

      // Binary vector for positions
      auto ID = [&ymax](int x, int y) { return x * ymax + y; };

      std::vector<bool> M(size_t(xmax) * size_t(ymax), false);
      for (int i = 0; i < n; ++i)
        M[ID(Rs.getX(i), Rs.getY(i))] = true;

      std::vector<int> H(size_t(xmax) * size_t(ymax), 0);
      for (int i = 0; i < n; ++i)
        H[ID(Rs.getX(i), Rs.getY(i))] = i;

      // Serve questo passaggio?
      // Rs.merge(ps); // CHECK: Posso mettere i pesi direttamente sul modello

      typedef double FlowType;
      typedef double CostType;
      // Build the graph for min cost flow
      NetSimplex<FlowType, CostType> simplex(
          'F', n + int(unbalanced == true),
          n *static_cast<int>(coprimes.size()));

      // Set the parameters
      simplex.setTimelimit(timelimit);
      simplex.setVerbosity(verbosity);
      simplex.setOptTolerance(opt_tolerance);

      // Initial set of arcs
      for (int h = 0; h < n; ++h) {
        int i = Rs.getX(h);
        int j = Rs.getY(h);
        for (const auto &p : coprimes) {
          int v = p.v;
          int w = p.w;
          if (i + v >= xmin && i + v < xmax && j + w >= ymin && j + w < ymax &&
              M[ID(i + v, j + w)]) {
            int ff = H[ID(i + v, j + w)];
            simplex.addArc(h, ff, p.c_vw);
          }
        }
      }

      // Add noded for unbalanced transport, if parater is set
      vector<size_t> lhs_arcs(n, 0);
      vector<size_t> rhs_arcs(n, 0);
      if (unbalanced) {
        for (int i = 0; i < n; ++i)
          lhs_arcs[i] = simplex.addArc(i, n, 0);

        for (int i = 0; i < n; ++i)
          rhs_arcs[i] = simplex.addArc(n, i, 0);
      }

      // Model attributes
      _num_arcs = simplex.num_arcs();

      if (verbosity == KWD_VAL_INFO)
        PRINT("INFO: running NetSimplex with V=%ld and E=%ld\n",
              simplex.num_nodes(), simplex.num_arcs());

      for (int jj = 0; jj < _m; ++jj) {

        //_runtime = 0.0;
        //_iterations = 0.0;

        // TODO: Devo ciclare sulla mappa iniziale (x,y)->idx
        // e usare quel bilancio ai nodi, se "i" è presente nella mappa
        // devo iterare su cosa?
        for (const auto &p : Rs.getM()) {
          auto q = MXY.find(p.first);
          if (q != MXY.end()) {
            simplex.addNode(p.second, W1[q->second] - Ws[jj * N + q->second]);
          } else
            simplex.addNode(p.second, 0.0);
        }

        if (unbalanced) {
          double bb = -tot_w1 + tot_ws[jj];
          simplex.addNode(n, bb); // Set the node value, it is not a true add

          double c1 = (bb < 0 ? unbal_cost : 0);
          double c2 = (bb < 0 ? 0 : unbal_cost);

          for (int i = 0; i < n; ++i)
            simplex.setArcCost(i, c1);

          for (int i = 0; i < n; ++i)
            simplex.setArcCost(i, c2);
        }

        _num_nodes = simplex.num_nodes();

        // Solve the problem to compute the distance
        _status = simplex.run();

        _runtime += simplex.runtime();
        _iterations += simplex.iterations();

        Ds[jj] = std::numeric_limits<CostType>::max();

        if (_status != ProblemType::INFEASIBLE &&
            _status != ProblemType::UNBOUNDED &&
            _status != ProblemType::TIMELIMIT) {
          Ds[jj] = simplex.totalCost();
          if (unbalanced)
            Ds[jj] = Ds[jj] / std::max(tot_w1, tot_ws[jj]);
        } else
          PRINT("ERROR 1001: Network Simplex wrong. Error code: %d\n",
                (int)_status);
      }

      return Ds;
    }

    // Second option for algorithm
    if (algorithm == KWD_VAL_COLGEN) {
      auto start_t = std::chrono::steady_clock::now();
      double _all_p = 0.0;

      // Compute xmax, ymax for each axis
      int xmax = std::numeric_limits<int>::min();
      int ymax = std::numeric_limits<int>::min();
      for (int i = 0; i < n; ++i) {
        xmax = std::max(xmax, Rs.getX(i));
        ymax = std::max(ymax, Rs.getY(i));
      }
      // The size is xmax+1, and ymax+1
      xmax++;
      ymax++;

      if (_L != L) {
        L = _L;
        init_coprimes(_L); // TODO: make it parallel
      }

      // Binary vector for positions
      auto ID = [&ymax](int x, int y) { return x * ymax + y; };

      std::vector<bool> M(size_t(xmax) * size_t(ymax), false);
      for (int i = 0; i < n; ++i)
        M[ID(Rs.getX(i), Rs.getY(i))] = true;

      std::vector<int> H(size_t(xmax) * size_t(ymax), 0);
      for (int i = 0; i < n; ++i)
        H[ID(Rs.getX(i), Rs.getY(i))] = i;

      typedef double FlowType;
      typedef double CostType;

      // Build the graph for min cost flow
      NetSimplex<FlowType, CostType> simplex('E', n + int(unbalanced == true),
                                             0);
      // Add noded for unbalanced transport, if parater is set
      vector<size_t> lhs_arcs(n, 0);
      vector<size_t> rhs_arcs(n, 0);
      if (unbalanced) {
        for (int i = 0; i < n; ++i)
          lhs_arcs[i] = simplex.addArc(i, n, 0);

        for (int i = 0; i < n; ++i)
          rhs_arcs[i] = simplex.addArc(n, i, 0);
      }

      // Set the parameters
      simplex.setTimelimit(timelimit);
      simplex.setVerbosity(verbosity);
      simplex.setOptTolerance(opt_tolerance);

      if (verbosity == KWD_VAL_INFO)
        PRINT("INFO: running NetSimplex with V=%ld and E=%ld\n",
              simplex.num_nodes(), simplex.num_arcs());

      for (int jj = 0; jj < _m; ++jj) {
        // TODO: Devo ciclare sulla mappa iniziale (x,y)->idx
        // e usare quel bilancio ai nodi, se "i" è presente nella mappa
        // devo iterare su cosa?
        for (const auto &p : Rs.getM()) {
          auto q = MXY.find(p.first);
          if (q != MXY.end()) {
            simplex.addNode(p.second, W1[q->second] - Ws[jj * N + q->second]);
          } else
            simplex.addNode(p.second, 0.0);
        }

        if (unbalanced) {
          double bb = -tot_w1 + tot_ws[jj];
          simplex.addNode(n, bb); // Set the node value, it is not a true add

          double c1 = (bb < 0 ? unbal_cost : 0);
          double c2 = (bb < 0 ? 0 : unbal_cost);

          for (int i = 0; i < n; ++i)
            simplex.setArcCost(i, c1);

          for (int i = 0; i < n; ++i)
            simplex.setArcCost(i, c2);
        }

        _num_nodes = simplex.num_nodes();

        int it = 0;
        int n_cuts = 0;
        // CostType fobj = 0;
        double negeps = std::nextafter(-opt_tolerance, -0.0);

        vector<double> pi(n, 0);

        Vars vars(n);
        for (int i = 0; i < n; ++i)
          vars[i].a = i;

        Vars vnew;
        vnew.reserve(n);

        // Init the simplex
        simplex.run();

        // Start separation
        while (true) {
          _status = simplex.reRun();
          if (_status == ProblemType::TIMELIMIT)
            break;

          // Take the dual values
          for (int j = 0; j < n; ++j)
            pi[j] = -simplex.potential(j);

          // Solve separation problem:
          auto start_tt = std::chrono::steady_clock::now();
#pragma omp parallel
          {
#pragma omp for schedule(dynamic, 1)
            for (int h = 0; h < n; ++h) {
              int a = Rs.getX(h);
              int b = Rs.getY(h);

              double best_v = negeps;
              double best_c = -1;
              int best_j = 0;

              for (const auto &p : coprimes) {
                int v = p.v;
                int w = p.w;
                double c_ij = p.c_vw;
                if (a + v >= 0 && a + v < xmax && b + w >= 0 && b + w < ymax &&
                    M[ID(a + v, b + w)]) {
                  int j = H[ID(a + v, b + w)];

                  double violation = c_ij - pi[h] + pi[j];
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

          // Take all negative reduced cost variables
          vnew.clear();
          for (auto &v : vars) {
            if (v.c > -1)
              vnew.push_back(v);
            v.c = -1;
          }

          auto end_tt = std::chrono::steady_clock::now();
          _all_p +=
              double(std::chrono::duration_cast<std::chrono::milliseconds>(
                         end_tt - start_tt)
                         .count()) /
              1000;

          if (vnew.empty())
            break;

          std::sort(vnew.begin(), vnew.end(),
                    [](const Var &v, const Var &w) { return v.c > w.c; });

          // Replace old constraints with new ones
          int new_arcs = simplex.updateArcs(vnew);

          n_cuts += new_arcs;

          ++it;
        }

        _iterations += simplex.iterations();
        _num_arcs = simplex.num_arcs();
        _num_nodes = simplex.num_nodes();
        auto end_t = std::chrono::steady_clock::now();
        auto _all = double(std::chrono::duration_cast<std::chrono::nanoseconds>(
                               end_t - start_t)
                               .count()) /
                    1000000000;

        Ds[jj] = simplex.totalCost();
        if (unbalanced)
          Ds[jj] = Ds[jj] / std::max(tot_w1, tot_ws[jj]);

        if (_n_log > 0)
          PRINT("it: %d, fobj: %f, all: %f, simplex: %f, all_p: %f\n", it,
                Ds[jj], _all, _runtime, _all_p);
      }
      auto end_t = std::chrono::steady_clock::now();
      auto _all = double(std::chrono::duration_cast<std::chrono::nanoseconds>(
                             end_t - start_t)
                             .count()) /
                  1000000000;
      _runtime = _all;

      return Ds;
    }

    return Ds;
  }

  // Alias for python module
  vector<double> compareApprox3(int _n, int _m, int *_Xs, int *_Ys, double *_Ws,
                                int _L) {
    return compareApprox(_n, _m, _Xs, _Ys, _Ws, _L);
  }

  vector<double> compareApprox(int _n, int _m, int *_Xs, int *_Ys, double *_Ws,
                               int _L) {
    // Check for correct input
    if (check_coding(_n, _Xs))
      PRINT(
          "WARNING: the Xs input coordinates are not consecutives integers.\n");
    if (check_coding(_n, _Ys))
      PRINT(
          "WARNING: the Ys input coordinates are not consecutives integers.\n");

    if (recode != "") {
      PRINT("INFO: Recoding the input coordinates to consecutive integers.\n");
      recoding(_n, _Xs);
      recoding(_n, _Ys);
    }

    if (verbosity == KWD_VAL_INFO)
      dumpParam();

    // Get the largest bounding box
    auto xy = getMinMax(_n, _Xs, _Ys);

    int xmin = xy[0];
    int ymin = xy[1];
    // int xmax = xy[2];
    // int ymax = xy[3];

    // Reindex and scale
    intpair2int XY = reindex(_n, _Xs, _Ys);
    int N = XY.size();

    // Binary vector for positions
    // auto WD = [&n](int x, int y) { return y * n + x; };

    vector<int> Xs, Ys;
    vector<double> W1, Ws;
    Xs.resize(N, 0);
    Ys.resize(N, 0);
    Ws.resize(size_t(N) * size_t(_m), 0.0);

    vector<double> tot_ws(_m, 0.0);

    intpair2int MXY;
    for (int i = 0; i < _n; i++) {
      int idx = XY.at(int_pair(_Xs[i], _Ys[i]));

      Xs[idx] = _Xs[i] - xmin;
      Ys[idx] = _Ys[i] - ymin;
      MXY[int_pair(Xs[idx], Ys[idx])] = idx;

      for (int j = 0; j < _m; ++j) {
        Ws[j * N + idx] += _Ws[j * _n + i];
        tot_ws[j] += _Ws[j * _n + i];
        if (_Ws[j * _n + i] < 0.0) {
          PRINT("WARNING: weight Ws[%d,%d]=%.4f is negative. Only positive "
                "weights are allowed.\n",
                i, j, _Ws[j * _n + i]);
          throw std::runtime_error(
              "FATAL ERROR: Input histogram with negative weigths");
        }
      }
    }

    // Rescale all integers coordinates to (0,0)
    if (!unbalanced) {
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < _m; ++j)
          Ws[j * N + i] = Ws[j * N + i] / tot_ws[j];
    }

    // Set the coprimes set
    if (_L != L) {
      L = _L;
      init_coprimes(_L); // TODO: make it parallel
    }

    // Return value
    vector<double> Ds(_m * _m, -1);

    for (int jj = 0; jj < _m; ++jj)
      Ds[jj * _m + jj] = 0.0;

    PointCloud2D ps = mergeCoordinates(N, &Xs[0], &Ys[0]);

    // Compute convex hull
    ConvexHull ch;
    PointCloud2D As, Rs;
    if (convex_hull) {
      As = ch.find(ps);
      Rs = ch.FillHull(As);
    } else {
      Rs = ch.FillHull(ps);
    }

    int n = static_cast<int>(Rs.size());

    // Second option for algorithm
    if (algorithm == KWD_VAL_MINCOSTFLOW) {
      // Compute xmax, ymax for each axis
      int xmin = std::numeric_limits<int>::max();
      int ymin = std::numeric_limits<int>::max();
      int xmax = std::numeric_limits<int>::min();
      int ymax = std::numeric_limits<int>::min();
      for (int i = 0; i < n; ++i) {
        xmin = std::min(xmin, Rs.getX(i));
        ymin = std::min(ymin, Rs.getY(i));
        xmax = std::max(xmax, Rs.getX(i));
        ymax = std::max(ymax, Rs.getY(i));
      }
      xmax++;
      ymax++;

      // Binary vector for positions
      auto ID = [&ymax](int x, int y) { return x * ymax + y; };

      std::vector<bool> M(size_t(xmax) * size_t(ymax), false);
      for (int i = 0; i < n; ++i)
        M[ID(Rs.getX(i), Rs.getY(i))] = true;

      std::vector<int> H(size_t(xmax) * size_t(ymax), 0);
      for (int i = 0; i < n; ++i)
        H[ID(Rs.getX(i), Rs.getY(i))] = i;

      // Serve questo passaggio?
      // Rs.merge(ps); // CHECK: Posso mettere i pesi direttamente
      // sul modello

      typedef double FlowType;
      typedef double CostType;
      // Build the graph for min cost flow
      NetSimplex<FlowType, CostType> simplex(
          'F', n + int(unbalanced == true),
          n *static_cast<int>(coprimes.size()));

      // Set the parameters
      simplex.setTimelimit(timelimit);
      simplex.setVerbosity(verbosity);
      simplex.setOptTolerance(opt_tolerance);

      // Initial set of arcs
      for (int h = 0; h < n; ++h) {
        int i = Rs.getX(h);
        int j = Rs.getY(h);
        for (const auto &p : coprimes) {
          int v = p.v;
          int w = p.w;
          if (i + v >= xmin && i + v < xmax && j + w >= ymin && j + w < ymax &&
              M[ID(i + v, j + w)]) {
            int ff = H[ID(i + v, j + w)];
            simplex.addArc(h, ff, p.c_vw);
          }
        }
      }

      // Add noded for unbalanced transport, if parater is set
      vector<size_t> lhs_arcs(n, 0);
      vector<size_t> rhs_arcs(n, 0);
      if (unbalanced) {
        for (int i = 0; i < n; ++i)
          lhs_arcs[i] = simplex.addArc(i, n, 0);

        for (int i = 0; i < n; ++i)
          rhs_arcs[i] = simplex.addArc(n, i, 0);
      }

      // Model attributes
      _num_arcs = simplex.num_arcs();

      if (verbosity == KWD_VAL_INFO)
        PRINT("INFO: running NetSimplex with V=%ld and E=%ld\n",
              simplex.num_nodes(), simplex.num_arcs());

      for (int ii = 0; ii < _m; ++ii) {
        for (int jj = ii + 1; jj < _m; ++jj) {
          // TODO: Devo ciclare sulla mappa iniziale (x,y)->idx
          // e usare quel bilancio ai nodi, se "i" è presente
          // nella mappa devo iterare su cosa?
          for (const auto &p : Rs.getM()) {
            auto q = MXY.find(p.first);
            if (q != MXY.end()) {
              simplex.addNode(p.second,
                              Ws[ii * N + q->second] - Ws[jj * N + q->second]);
            } else
              simplex.addNode(p.second, 0.0);
          }

          if (unbalanced) {
            double bb = -tot_ws[ii] + tot_ws[jj];
            simplex.addNode(n, bb); // Set the node value, it is not a true add

            double c1 = (bb < 0 ? unbal_cost : 0);
            double c2 = (bb < 0 ? 0 : unbal_cost);

            for (int i = 0; i < n; ++i)
              simplex.setArcCost(i, c1);

            for (int i = 0; i < n; ++i)
              simplex.setArcCost(i, c2);
          }

          _num_nodes = simplex.num_nodes();

          // Solve the problem to compute the distance
          _status = simplex.run();

          _runtime += simplex.runtime();
          _iterations += simplex.iterations();

          Ds[ii * _m + jj] = std::numeric_limits<CostType>::max();
          Ds[jj * _m + ii] = std::numeric_limits<CostType>::max();

          if (_status != ProblemType::INFEASIBLE &&
              _status != ProblemType::UNBOUNDED &&
              _status != ProblemType::TIMELIMIT) {
            Ds[ii * _m + jj] = simplex.totalCost();
            if (unbalanced)
              Ds[ii * _m + jj] =
                  Ds[ii * _m + jj] / std::max(tot_ws[ii], tot_ws[jj]);
            Ds[jj * _m + ii] = Ds[ii * _m + jj];
          } else
            PRINT("ERROR 1001: Network Simplex wrong. Error code: "
                  "%d\n",
                  (int)_status);
        }
      }

      return Ds;
    }

    // Second option for algorithm
    if (algorithm == KWD_VAL_COLGEN) {
      auto start_t = std::chrono::steady_clock::now();
      double _all_p = 0.0;

      // Compute xmax, ymax for each axis
      int xmax = std::numeric_limits<int>::min();
      int ymax = std::numeric_limits<int>::min();
      for (int i = 0; i < n; ++i) {
        xmax = std::max(xmax, Rs.getX(i));
        ymax = std::max(ymax, Rs.getY(i));
      }
      // The size is xmax+1, and ymax+1
      xmax++;
      ymax++;

      if (_L != L) {
        L = _L;
        init_coprimes(_L); // TODO: make it parallel
      }

      // Binary vector for positions
      auto ID = [&ymax](int x, int y) { return x * ymax + y; };

      std::vector<bool> M(size_t(xmax) * size_t(ymax), false);
      for (int i = 0; i < n; ++i)
        M[ID(Rs.getX(i), Rs.getY(i))] = true;

      std::vector<int> H(size_t(xmax) * size_t(ymax), 0);
      for (int i = 0; i < n; ++i)
        H[ID(Rs.getX(i), Rs.getY(i))] = i;

      typedef double FlowType;
      typedef double CostType;

      // Build the graph for min cost flow
      NetSimplex<FlowType, CostType> simplex('E', n + int(unbalanced == true),
                                             0);

      // Set the parameters
      simplex.setTimelimit(timelimit);
      simplex.setVerbosity(verbosity);
      simplex.setOptTolerance(opt_tolerance);

      // Add noded for unbalanced transport, if parater is set
      vector<size_t> lhs_arcs(n, 0);
      vector<size_t> rhs_arcs(n, 0);
      if (unbalanced) {
        for (int i = 0; i < n; ++i)
          lhs_arcs[i] = simplex.addArc(i, n, 0);

        for (int i = 0; i < n; ++i)
          rhs_arcs[i] = simplex.addArc(n, i, 0);
      }

      if (verbosity == KWD_VAL_INFO)
        PRINT("INFO: running NetSimplex with V=%ld and E=%ld\n",
              simplex.num_nodes(), simplex.num_arcs());

      for (int ii = 0; ii < _m; ++ii) {
        for (int jj = ii + 1; jj < _m; ++jj) {
          // TODO: Devo ciclare sulla mappa iniziale (x,y)->idx
          // e usare quel bilancio ai nodi, se "i" è presente
          // nella mappa devo iterare su cosa?
          for (const auto &p : Rs.getM()) {
            auto q = MXY.find(p.first);
            if (q != MXY.end()) {
              simplex.addNode(p.second,
                              Ws[ii * N + q->second] - Ws[jj * N + q->second]);
            } else
              simplex.addNode(p.second, 0.0);
          }

          if (unbalanced) {
            double bb = -tot_ws[ii] + tot_ws[jj];
            simplex.addNode(n, bb); // Set the node value, it is not a true add

            double c1 = (bb < 0 ? unbal_cost : 0);
            double c2 = (bb < 0 ? 0 : unbal_cost);

            for (int i = 0; i < n; ++i)
              simplex.setArcCost(i, c1);

            for (int i = 0; i < n; ++i)
              simplex.setArcCost(i, c2);
          }

          _num_nodes = simplex.num_nodes();

          int it = 0;
          int n_cuts = 0;
          // CostType fobj = 0;
          double negeps = std::nextafter(-opt_tolerance, -0.0);

          vector<double> pi(n, 0);

          Vars vars(n);
          for (int i = 0; i < n; ++i)
            vars[i].a = i;

          Vars vnew;
          vnew.reserve(n);

          // Init the simplex
          simplex.run();

          // Start separation
          while (true) {
            _status = simplex.reRun();
            if (_status == ProblemType::TIMELIMIT)
              break;

            // Take the dual values
            for (int j = 0; j < n; ++j)
              pi[j] = -simplex.potential(j);

            // Solve separation problem:
            auto start_tt = std::chrono::steady_clock::now();
#pragma omp parallel
            {
#pragma omp for schedule(dynamic, 1)
              for (int h = 0; h < n; ++h) {
                int a = Rs.getX(h);
                int b = Rs.getY(h);

                double best_v = negeps;
                double best_c = -1;
                int best_j = 0;

                for (const auto &p : coprimes) {
                  int v = p.v;
                  int w = p.w;
                  double c_ij = p.c_vw;
                  if (a + v >= 0 && a + v < xmax && b + w >= 0 &&
                      b + w < ymax && M[ID(a + v, b + w)]) {
                    int j = H[ID(a + v, b + w)];

                    double violation = c_ij - pi[h] + pi[j];
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

            // Take all negative reduced cost variables
            vnew.clear();
            for (auto &v : vars) {
              if (v.c > -1)
                vnew.push_back(v);
              v.c = -1;
            }

            auto end_tt = std::chrono::steady_clock::now();
            _all_p +=
                double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_tt - start_tt)
                           .count()) /
                1000;

            if (vnew.empty())
              break;

            std::sort(vnew.begin(), vnew.end(),
                      [](const Var &v, const Var &w) { return v.c > w.c; });

            // Replace old constraints with new ones
            int new_arcs = simplex.updateArcs(vnew);

            n_cuts += new_arcs;

            ++it;
          }

          _iterations += simplex.iterations();
          _num_arcs = simplex.num_arcs();
          _num_nodes = simplex.num_nodes();
          auto end_t = std::chrono::steady_clock::now();
          auto _all =
              double(std::chrono::duration_cast<std::chrono::nanoseconds>(
                         end_t - start_t)
                         .count()) /
              1000000000;

          Ds[jj * _m + ii] = simplex.totalCost();
          if (unbalanced)
            Ds[jj * _m + ii] =
                Ds[jj * _m + ii] / std::max(tot_ws[ii], tot_ws[jj]);

          Ds[ii * _m + jj] = Ds[jj * _m + ii];

          if (_n_log > 0)
            PRINT("it: %d, fobj: %f, all: %f, simplex: %f, all_p: "
                  "%f\n",
                  it, Ds[jj * _m + ii], _all, _runtime, _all_p);
        }
      }

      auto end_t = std::chrono::steady_clock::now();
      auto _all = double(std::chrono::duration_cast<std::chrono::nanoseconds>(
                             end_t - start_t)
                             .count()) /
                  1000000000;
      _runtime = _all;

      return Ds;
    }

    return Ds;
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
      Rs.update(k.first.first - xmin, k.first.second - ymin, -k.second);

    // Use as few memory as possible
    Rs.shrink_to_fit();

    return Rs;
  }

  PointCloud2D mergeCoordinates(size_t n, int *Xs, int *Ys) {
    // PREREQUISITIES: Ensures that Xs[i],Ys[i] have not duplicates

    PointCloud2D Rs;
    Rs.reserve(n);

    // Point cloud
    for (size_t i = 0; i < n; ++i)
      Rs.add(Xs[i], Ys[i], 0.0);

    // Use as few memory as possible
    Rs.shrink_to_fit();

    return Rs;
  }

  PointCloud2D mergeHistograms(size_t n, int *Xs, int *Ys, double *W1,
                               double *W2) {
    int xmin = std::numeric_limits<int>::max();
    int ymin = std::numeric_limits<int>::max();
    for (size_t i = 0; i < n; ++i) {
      xmin = std::min(xmin, Xs[i]);
      ymin = std::min(ymin, Ys[i]);
    }

    PointCloud2D Rs;
    Rs.reserve(n);

    // Read first line
    for (size_t i = 0; i < n; ++i)
      Rs.update(Xs[i] - xmin, Ys[i] - ymin, W1[i] - W2[i]);

    // Use as few memory as possible
    Rs.shrink_to_fit();

    return Rs;
  }

  // Get min and max of two coordinates
  std::array<int, 4> getMinMax(size_t n, int *Xs, int *Ys) {
    // Compute xmin, xmax, ymin, ymax for each axis
    int xmax = std::numeric_limits<int>::min();
    int ymax = std::numeric_limits<int>::min();
    int xmin = std::numeric_limits<int>::max();
    int ymin = std::numeric_limits<int>::max();
    for (size_t i = 0; i < n; ++i) {
      xmax = std::max(xmax, Xs[i]);
      ymax = std::max(ymax, Ys[i]);

      xmin = std::min(xmin, Xs[i]);
      ymin = std::min(ymin, Ys[i]);
    }
    return {xmin, ymin, xmax, ymax};
  }

  // Status of the solver
  ProblemType _status;

  // Runtime in milliseconds
  double _runtime;

  // Number of iterations
  uint64_t _iterations;
  uint64_t _num_nodes;
  uint64_t _num_arcs;

  // Interval for logging iterations in the simplex algorithm
  // (if _n_log=0 no logs at all)
  int _n_log;

  // Approximation parameter
  int L;

  // List of pair of coprimes number between (-L, L)
  std::vector<coprimes_t> coprimes;

  // Method to solve the problem
  std::string method;
  // Model to solve the problem
  std::string model;
  // Algorithm to solve the corresponding problem
  std::string algorithm;
  // Verbosity of the log
  std::string verbosity;
  // Recode the coordinates as consecutive integers
  std::string recode;
  // Tolerance for pricing
  double opt_tolerance;
  // Time limit for runtime of the algorithm
  double timelimit;
  // If the problem must be considered unbalanced
  bool unbalanced;
  // Cost for the unbalanced connection
  double unbal_cost;
  // Whether to compute the convex hull
  bool convex_hull;

}; // namespace KWD

} // end namespace KWD
