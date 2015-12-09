#include "variates.h"

int main(int argc, char const *argv[]) {
  double alpha = argc > 1 ? atof(argv[1]) : 1.6;
  double a     = argc > 2 ? atof(argv[2]) : 1.0;
  int J        = argc > 3 ? atoi(argv[3]) : 100;
  int N        = argc > 4 ? atoi(argv[4]) : 1000;
  int seed     = argc > 5 ? atof(argv[5]) : 1231;
  int sim      = argc > 6 ? atoi(argv[6]) : 0;

  mt19937 rng(seed);

  if (sim == 0) {
    double dt = 1.0 / J;
    StableDistribution dist(alpha, dt * a);
    double X;
    for (int i = 0; i < N; i++) {
      X = 0.0;
      for (int j = 0; j < J; j++)
        X += dist(rng);
      printf("%.4f\n", X);
    }
  } else {
    StableDistribution dist(alpha, a);
    double X;
    for (int i = 0; i < N; i++) {
      X = dist(rng);
      printf("%.4f\n", X);
    }
  }

  return 0;
}
