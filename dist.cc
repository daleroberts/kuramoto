#include <iostream>
#include <stdio.h>
#include "variates.h"

using namespace std;

int main(int argc, char const *argv[]) {
  int seed      = argc > 1 ? atol(argv[1]) : 1234;
  double alpha  = argc > 2 ? atof(argv[2]) : 1.5;
  double lambda = argc > 3 ? atof(argv[3]) : 1.0;
  double sigma  = argc > 4 ? atof(argv[4]) : 5.78115129;
  double dt     = argc > 5 ? atof(argv[5]) : 0.1;
  int N         = argc > 6 ? atof(argv[6]) : 1000;

  double a = alpha * pow(sigma, 2.0) / (2 * tgamma(1 - alpha) * cos(M_PI * alpha / 2));
  double b = lambda;

  double drift = dt * tgamma(1 - a) * a * pow(b, alpha - 1);
  double x     = qtstable(0.01, alpha, a, b);
  double c     = -drift - x;

  cerr << a << " " << b << " " << c << endl;

  mt19937 rng(seed);

  TemperedStableDistribution dist(alpha, dt * a, b, c);
  double X;

  setbuf(stdout, NULL);
  for (int i = 0; i < N; i++) {
    X = dist(rng);
    printf("%.4f\n", X);
  }

  return 0;
}
