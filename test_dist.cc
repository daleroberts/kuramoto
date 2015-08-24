#include "variates.h"

int main(int argc, char const* argv[]) {
  double alpha = argc > 1 ? atof(argv[1]) : 1.6;
  double beta = argc > 2 ? atof(argv[2]) : 0.0;
  double sigma = argc > 3 ? atof(argv[3]) : 2.0;
  int seed = argc > 4 ? atof(argv[4]) : 1231;

  double a = 0.5 * alpha * pow(sigma, 2.0) /
             (tgamma(1 - alpha) * cos(M_PI * alpha / 2));

  mt19937 rng(seed);
  StableDistribution dist(alpha, a);

  double xi;
  for (int i = 0; i < 10000; i++) {
    xi = ((1+beta)*dist(rng)-(1-beta)*dist(rng))/2;
    printf("%.6f\n", xi);
  }

  return 0;
}
