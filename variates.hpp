#pragma once
#include <random>

using namespace std;

class stable_distribution
{
public:
  explicit
  stable_distribution(const double& alpha = double(0.5),
                      const double& a = double(1))
  : _alpha(alpha), _a(a), _runif(-M_PI_2,M_PI_2), _rexp(1)
  { 
    _factor = pow(-a*tgamma(-alpha)*cos(M_PI_2*alpha), 1./alpha);
    _theta = atan(tan(M_PI_2*alpha));
  }

  double
  alpha() const
  { return _alpha; }

  double
  a() const
  { return _a; }

  template<class _UniformRandomNumberGenerator>
    double
    operator()(_UniformRandomNumberGenerator& urng) {
      double U, E, X;
      U = _runif(urng);
      E = _rexp(urng);
      X = sin(_alpha*U + _theta)/pow(cos(U)*cos(_theta), 1./_alpha) *
          pow(cos((1-_alpha)*U - _theta)/E,(1-_alpha)/_alpha);
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

class tempered_stable_distribution
{
public:
  explicit
  tempered_stable_distribution(const double& alpha = double(0.5),
                               const double& a = double(1),
                               const double& b = double(1))
  : _rstable(alpha, a), _runif(0,1), _b(b)
  { }

  double
  alpha() const
  { return _rstable.alpha(); }

  double
  a() const
  { return _rstable.a(); }

  double
  b() const
  { return _b; }

  template<class _UniformRandomNumberGenerator>
    double
    operator()(_UniformRandomNumberGenerator& urng) {
      double U, V, X;
      do {
        U = _runif(urng);
        V = _rstable(urng);
      } while (U > exp(-_b*V));
      return V;
    }

private:
  double _b;
  stable_distribution _rstable;
  uniform_real_distribution<> _runif;
};
