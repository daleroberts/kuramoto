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

template <typename Derived>
inline double order_param(const Eigen::MatrixBase<Derived>& theta) {
    auto N = theta.size();

    double sum_real = 0.0;
    double sum_complex = 0.0;
    for (size_t j = 0; j < N; j++) {
        sum_real += cos(theta(j));
        sum_complex += sin(theta(j));
    }

    return sqrt(sum_real * sum_real + sum_complex * sum_complex) / N;
}

template <typename Distribution>
void one_path(Graph& G,
              Distribution& dist,
              const double alpha,
              const double a,
              const double b,
              const double K,
              const double max_t,
              const uint32_t nsteps,
              const uint32_t seed) {
    double dt = max_t / nsteps;
    auto N = G.size();

    mt19937 rng(seed);
    uniform_real_distribution<> runif(-M_PI, M_PI);

    Vector theta = Vector::Zero(N);
    Vector theta_ = Vector::Zero(N);

    for (size_t i = 0; i < N; ++i)
        theta(i) = runif(rng);

    double xi, drift;
    double X = 0.0;

    printf("{");

    X = 0.0;
    //for (size_t i = 0; i < N; i++)
    //    X += theta(i);
    //X = X/N;
    printf("{0.0, %.4f},", X);

    for (size_t k = 1; k <= nsteps; k++) {

        for (size_t i = 0; i < N; i++) {
            drift = 0;
            for (auto& j : G.neighbours(i))
                drift -= K / N * sin(theta(i) - theta(j));

            xi = dist(rng);
            theta_(i) += drift * dt + xi;
            //theta_(i) += xi;
        }
    
        X = 0.0;
        for (size_t i = 0; i < N; ++i)
            X += theta_(i);
        X = X/N;

        if (k < nsteps)
            printf("{%.4f, %.4f},", k*dt, X);
        else
            printf("{%.4f, %.4f}}\n", k*dt, X);

        theta = theta_;
    }
}

int main(int argc, char const* argv[]) {
    //setbuf(stdout, NULL);

    const char* graphfile = argc > 1 ? argv[1] : "graphs.g6";
    int seed = argc > 2 ? atol(argv[2]) : 1234;
    int nsteps = argc > 3 ? atoi(argv[3]) : 5000;
    double alpha = argc > 4 ? atof(argv[4]) : 1.7;   // stability
    double lambda = argc > 5 ? atof(argv[5]) : 1.0;  // tempering
    double sigma = argc > 6 ? atof(argv[6]) : 1.0;   // diffusivity
    double K = argc > 7 ? atof(argv[7]) : 0.8;       // global coupling
    double max_t = argc > 8 ? atof(argv[8]) : 30.0;  // maximum time
    size_t M = argc > 9 ? atof(argv[9]) : 1;

    double dt = max_t / nsteps;
    double a = dt * 0.5 * alpha * pow(sigma, 2.0) /
        (tgamma(1 - alpha) * cos(M_PI * alpha / 2));
    double b = lambda;

    StableDistribution rstable(alpha, dt * a);
    NormalDistribution rnorm(0, dt * pow(sigma, 2));

    ifstream infile(graphfile);
    string line;
    while (getline(infile, line)) {
        Graph G(line);
        for (size_t m = 0; m < M; m++) {
            if (alpha > 1.999) {
                one_path(G, rnorm, alpha, a, b, K, max_t, nsteps, seed + m);
            } else {
                one_path(G, rstable, alpha, a, b, K, max_t, nsteps, seed + m);
            }
        }

        break;
    }

    return 0;
}
