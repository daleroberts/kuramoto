#pragma once
#include <random>
#include <complex>
#include <cmath>

using namespace std;

typedef std::normal_distribution<> NormalDistribution;
typedef std::complex<double> Complex;

class StableDistribution {
 public:
  explicit StableDistribution(const double& alpha = double(0.5),
                              const double& a = double(1))
      : _alpha(alpha), _a(a), _runif(-M_PI_2, M_PI_2), _rexp(1) {
    _factor = pow(-a * tgamma(-alpha) * cos(M_PI_2 * alpha), 1. / alpha);
    _theta = atan(tan(M_PI_2 * alpha));
  }

  double alpha() const { return _alpha; }

  double a() const { return _a; }

  template <class _UniformRandomNumberGenerator>
  double operator()(_UniformRandomNumberGenerator& urng) {
    double U, E, X;
    U = _runif(urng);
    E = _rexp(urng);
    X = sin(_alpha * U + _theta) / pow(cos(U) * cos(_theta), 1. / _alpha) *
        pow(cos((1 - _alpha) * U - _theta) / E, (1 - _alpha) / _alpha);
    return _factor * X;
  }

 private:
  double _alpha;
  double _a;
  double _factor;
  double _theta;
  uniform_real_distribution<> _runif;
  exponential_distribution<> _rexp;
};

Complex ctstable(double y, double alpha, double a, double b);
double ptstable(double x, double alpha, double a, double b);
double qtstable(double p, double alpha, double a, double b);

class TemperedStableDistribution {
 public:
  explicit TemperedStableDistribution(const double& alpha = double(0.5),
                                      const double& a = double(1),
                                      const double& b = double(1),
                                      const double& c = double(1.1))
      : _alpha(alpha), _a(a), _b(b), _c(c), _rstable(alpha, a), _runif(0, 1) {}

  double alpha() const { return _alpha; }

  double a() const { return _a; }

  double b() const { return _b; }

  template <class _UniformRandomNumberGenerator>
  double operator()(_UniformRandomNumberGenerator& urng) {
    double U, V;
    if (_alpha < 1) {
      do {
        U = _runif(urng);
        V = _rstable(urng);
      } while (U > exp(-_b * V));
      return V;

    } else {
      do {
        U = _runif(urng);
        V = _rstable(urng);
      } while (U > exp(-_b * (V + _c)));
      return V - tgamma(1 - _alpha) * _a * pow(_b, _alpha - 1.);
    }
  }

 private:
  double _alpha;
  double _a;
  double _b;
  double _c;
  StableDistribution _rstable;
  uniform_real_distribution<> _runif;
};
