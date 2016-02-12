#include <complex>
#include <cmath>
#include <iostream>
#include <limits>
#include <boost/math/tools/roots.hpp>

#include "variates.h"

using namespace std;

inline Complex ctstable(double y, double alpha, double a, double b) {
  Complex I(0, 1);
  Complex A = a * tgamma(-alpha) *
              (pow(b - I * y, alpha) - pow(b, alpha) + I * y * alpha * pow(b, alpha - 1));
  return exp(A);
}

double ptstable(double x, double alpha, double a, double b) {
  // moments
  double m1, m2;
  if (alpha < 1) {
    m1 = a * pow(b, alpha - 1) * tgamma(1 - alpha);
    m2 =
        a * pow(b, alpha - 2) * (a * pow(b, alpha) * pow(tgamma(1 - alpha), 2) + tgamma(2 - alpha));
  } else {
    m1 = 0.0;
    m2 = a * pow(b, alpha - 2) * tgamma(2 - alpha);
  }

  // step size
  int q    = 5;
  double h = 2 * M_PI / (x + abs(m1) + q * abs(sqrt(m2 - m1)));

  double sum = h * x / M_PI;
  double term;
  Complex psi;
  int j = 1;

  do {
    psi  = ctstable(h * j, alpha, a, b);
    term = 2 / M_PI * sin(h * j * x) * psi.real() / j;
    sum += term;
    j++;
  } while ((abs(term) >= numeric_limits<double>::epsilon()) && (j < 50));

  return sum;
}

double qtstable(double p, double alpha, double a, double b) {
  auto f = [=](double x) { return ptstable(x, alpha, a, b) - p; };
  auto digits              = std::numeric_limits<double>::digits;
  auto tol                 = boost::math::tools::eps_tolerance<double>(digits);
  boost::uintmax_t maxiter = 100;
  auto r = boost::math::tools::bracket_and_solve_root(f, 0.05, 1.1, true, tol, maxiter);
  return (r.first + r.second) / 2.0;
}
