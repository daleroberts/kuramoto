/*
 * Kuramoto Model with Tempered Stable Noise
 *
 * Dale Roberts <dale.o.roberts@gmail.com>
 */

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <sstream>
#include <Eigen/Dense>

#include "graph.h"
#include "variates.h"
#include "statistics.h"

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0; }
inline omp_int_t omp_get_max_threads() { return 1; }
#endif

using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

template <typename Iterable> inline void elementwise_add(Iterable &in, Iterable &out) {
  auto in_iter = in.begin();
  for (auto &el : out)
    el += *(in_iter++);
}

inline double mod2pi(double theta) {
  while (theta > 2 * M_PI)
    theta -= 2 * M_PI;
  while (theta < 0.0)
    theta += 2 * M_PI;
  return theta;
}

template <typename Derived> inline double sum(const Eigen::MatrixBase<Derived> &v) {
  double total = 0.;
  for (int j = 0; j < v.size(); j++) {
    total += v[j];
  }
  return total;
}

template <typename Derived> inline double order_param(const Eigen::MatrixBase<Derived> &theta) {
  auto N             = theta.size();
  double sum_real    = 0.0;
  double sum_complex = 0.0;
  for (size_t j = 0; j < N; j++) {
    sum_real += cos(theta(j));
    sum_complex += sin(theta(j));
  }

  return sqrt(sum_real * sum_real + sum_complex * sum_complex) / N;
}

template <typename Distribution>
void one_path(Graph &G, Distribution &dist, const double alpha, const double a, const double b,
              const double K, const double max_t, const uint32_t nsteps, const uint32_t seed) {
  double dt = max_t / nsteps;
  auto N    = G.size();

  mt19937 rng(seed);
  uniform_real_distribution<> runif(-M_PI, M_PI);

  Vector theta(N);
  Vector theta_(N);
  Vector xi(N);

  // set initial condition
  for (size_t i = 0; i < N; ++i)
    theta[i] = runif(rng);

  theta[N - 1] = -(sum(theta) - theta[N - 1]); // adjust for zero expectation

  printf("{{%.6f", theta(0));
  for (size_t i = 1; i < N; ++i)
    printf(",%.6f", theta(i));
  printf("}");

  double drift;

  for (size_t k = 0; k <= nsteps; k++) {

    for (size_t i = 0; i < N; i++)
      xi[i] = dist(rng);

    xi[N - 1] = -(sum(xi) - xi[N - 1]); // adjust for zero expectation

    for (size_t i = 0; i < N; i++) {
      drift = 0;
      for (auto &j : G.neighbours(i))
        drift -= K / N * sin(theta[i] - theta[j]);

      theta_[i] += drift * dt + xi[i];
      //theta_[i] += xi[i];
    }

    theta = theta_;

    printf(",{%.6f", theta(0));
    for (size_t i = 1; i < N; ++i)
      printf(",%.6f", theta(i));
    printf("}");
  }

  printf("}\n");
}

int main(int argc, char const *argv[]) {
  const char *graphfile = argc > 1 ? argv[1] : "graphs.g6";
  int seed              = argc > 2 ? atol(argv[2]) : 1234;
  int nsteps            = argc > 3 ? atoi(argv[3]) : 5000;
  double alpha          = argc > 4 ? atof(argv[4]) : 1.7;  // stability
  double lambda         = argc > 5 ? atof(argv[5]) : 1.0;  // tempering
  double sigma          = argc > 6 ? atof(argv[6]) : 1.0;  // diffusivity
  double K              = argc > 7 ? atof(argv[7]) : 0.8;  // global coupling
  double max_t          = argc > 8 ? atof(argv[8]) : 30.0; // maximum time

  double dt = max_t / nsteps;
  double a  = dt * 0.5 * alpha * pow(sigma, 2.0) / (tgamma(1 - alpha) * cos(M_PI * alpha / 2));
  double b  = lambda;

  setbuf(stdout, NULL);
  ifstream infile(graphfile);
  string line;
  while (getline(infile, line)) {
    Graph G(line);
    if (alpha > 1.999) {
      NormalDistribution rnorm(0, dt * pow(sigma, 2) / 2); // scaled!
      one_path(G, rnorm, alpha, a, b, K, max_t, nsteps, seed);
    } else {
      if (lambda < 0.001) {
        StableDistribution rstable(alpha, dt * a);
        one_path(G, rstable, alpha, a, b, K, max_t, nsteps, seed);
      } else {
        // determine c
        double d = dt * tgamma(1 - a) * a * pow(b, alpha - 1);
        double c = -d - qtstable(0.9, alpha, a, b);
        TemperedStableDistribution rtstable(alpha, dt * a, b, c);
        one_path(G, rtstable, alpha, a, b, K, max_t, nsteps, seed);
      }
    }

    break; // one graph
  }

  return 0;
}
