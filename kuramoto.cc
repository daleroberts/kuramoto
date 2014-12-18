/*
 * Kuramoto Model with Tempered Stable Noise
 *
 * Dale Roberts <dale.o.roberts@gmail.com>
 */

#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <sstream>
#include <Eigen/Dense>

#include "graph.h"
#include "variates.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

inline double mod2pi(double theta) {
    while (theta > 2*M_PI)
        theta -= 2*M_PI;
    while (theta < 0.0)
        theta += 2*M_PI;
    return theta;
}

template <typename Derived>
inline double order_param(Eigen::MatrixBase<Derived>& theta) {
    auto N = theta.size();

    double sum_real = 0.0;
    double sum_complex = 0.0;
    for (size_t j = 0; j < N; j++) {
        sum_real += cos(theta(j));
        sum_complex += sin(theta(j));
    }

    return sqrt(sum_real*sum_real + sum_complex*sum_complex)/N;
}


inline double sinc(double x) {
    double ax = fabs(x);
    if (ax < 0.006) {
        if (x == 0.)
            return 1;
        double x2 = x*x;
        if(ax < 2e-4)
            return 1. - x2/6.;
        else 
            return 1. - x2/6.*(1 - x2/20.);
    }
    /* else */
    return sin(x)/x;
}

inline double A(double x, double rho) {
    double Irho = 1.-rho;
    return pow(Irho*sinc(Irho*x),Irho)*pow(rho*sinc(rho*x),rho)/sinc(x);
}

template <typename DerivedA, typename DerivedB>
void paths(const double globalCoupling,
           const double noiseScale,
           const Eigen::MatrixBase<DerivedA>& initialConditions,
           const Eigen::MatrixBase<DerivedB>& intrinsicFreqs,
           Graph& G,
           const double rho,
           const double c,
           const double alpha,
           const double p,
           const double t,
           const uint32_t timeSteps,
           const uint32_t seed,
           const uint32_t npaths) {

    mt19937 rng(seed);
    uniform_real_distribution<> runif(0,1);
    exponential_distribution<> rexp(1);

    double dt = t/timeSteps;

    double gamma = c/(1-rho) - p;
    double cf = pow(cos(M_PI_2*rho),-1./rho);
    double sigma = pow(-dt*c*cos(M_PI*rho/2.)*tgamma(-rho), 1./rho);
    double mu = -dt*p;

    auto N = initialConditions.size();

    Vector avg_order(timeSteps+1);

    for (size_t j=0; j < npaths; ++j) {
        Vector theta(N);
        for (size_t i = 0; i < N; ++i)
            theta(i) = initialConditions(i);

        avg_order(0) += order_param(theta);

        double U, V, W, xi, drift;
        double s = 0.0;
        long k = 0;

        do {
            s += dt;
            k++;

            for (size_t i = 0; i < N; i++) {
                // simulate a tempered stable random variable
                do {
                    U = runif(rng);
                    V = runif(rng);
                    do {
                        W = rexp(rng);
                    } while (W == 0.);
                    xi = mu + sigma*cf*pow(A(M_PI*U,rho)/pow(W,1.-rho),1./rho);
                } while (V > exp(-alpha*xi));

                // calculate the drift
                drift = intrinsicFreqs(i);
                for (auto &j: G.neighbours(i))
                    drift -= globalCoupling/N*sin(theta(i) - theta(j));

                // increment process
                theta(i) += drift*dt + noiseScale*xi;
                theta(i) = mod2pi(theta(i));
            }

            avg_order(k) += order_param(theta);

        } while (s < t);
    }

    avg_order /= npaths;

    for (size_t i = 0; i < avg_order.size(); ++i) {
        printf("%6.3f %6.3f\n", i*dt, avg_order(i));
    }
}


int main(int argc, char const *argv[]) {
    double rho = 0.99; // stability index
    double c = 0.01;
    double alpha = 1.0;
    double xi = 0.2;
    double p = (1+xi)*(-c*rho*tgamma(-rho)*pow(alpha,rho-1));

    double t = 30.0;
    int timeSteps = 3000;
    int seed = 1234;
    int npaths = 500;

    double globalCoupling = 0.8;
    double noiseScale = 0.05;

    string filename = "graphs.g6";
    ifstream graphfile(filename);
    string line;
    while (getline(graphfile, line)) {
        line.erase(line.find_last_not_of("\n\r")+1);
        Graph G(line);
        auto n = G.size();
        printf("rho: %f c: %f alpha: %f xi: %f p: %f graph size: %lu\n", rho, c, alpha, xi, p, n);
        
        Vector initialConditions = Vector::LinSpaced(n, 0, 3.14);
        Vector intrinsicFreqs = Vector::LinSpaced(n, -1, 1);

        paths(globalCoupling, noiseScale, initialConditions, intrinsicFreqs,
              G, rho, c, alpha, p, t, timeSteps, seed, npaths);

        cout << "----" << endl;
    }

    return 0;
}
