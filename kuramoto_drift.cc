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

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

inline double mod2pi(double theta) {
  while (theta > 2 * M_PI)
    theta -= 2 * M_PI;
  while (theta < 0.0)
    theta += 2 * M_PI;
  return theta;
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

  Vector theta  = Vector::Zero(N);
  Vector theta_ = Vector::Zero(N);

  for (size_t i = 0; i < N; ++i)
    theta(i) = runif(rng);

  double xi, drift;

  for (size_t k = 0; k < nsteps; k++) {

    printf("{");

    for (size_t i = 0; i < N; i++) {
      drift = 0;
      for (auto &j : G.neighbours(i))
        drift -= K / N * sin(theta(i) - theta(j));

      if (i < N - 1)
        printf("{%.8f,%.8f,%.8f},", theta(i), drift * dt, xi);
      else
        printf("{%.8f,%.8f,%.8f}", theta(i), drift * dt, xi);

      xi = dist(rng);
      theta_(i) += drift * dt + xi;
    }

    printf("}\n");

    theta = theta_;
  }
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
  double d  = dt * tgamma(1 - a) * a * pow(b, alpha - 1);
  double c  = -d - qtstable(0.9, alpha, a, b);

  TemperedStableDistribution rtstable(alpha, dt * a, b, c);
  StableDistribution rstable(alpha, dt * a);
  NormalDistribution rnorm(0, dt * pow(sigma, 2) / 2);

  setbuf(stdout, NULL);
  ifstream infile(graphfile);
  string line;
  while (getline(infile, line)) {
    Graph G(line);
    if (alpha > 1.999) {
      one_path(G, rnorm, alpha, a, b, K, max_t, nsteps, seed);
    } else {
      if (lambda < 0.001) {
        one_path(G, rstable, alpha, a, b, K, max_t, nsteps, seed);
      } else {
        one_path(G, rtstable, alpha, a, b, K, max_t, nsteps, seed);
      }
    }

    break;
  }

  return 0;
}
