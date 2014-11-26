/* Kuramoto Model with Tempered Stable Noise
 * Dale Roberts <dale.o.roberts@gmail.com>
 */

#include <random>
#include <cmath>
#include <Eigen/Dense>
#include <stdio.h>
#include <iostream>
#include <unordered_map>
#include "rv.hpp"

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
    int N = theta.size();

    double sum_real = 0.0;
    double sum_complex = 0.0;
    for (int j = 0; j < N; j++) {
        sum_real += cos(theta(j));
        sum_complex += sin(theta(j));
    }

    return sqrt(sum_real*sum_real + sum_complex*sum_complex)/N;
}

template <typename DerivedA, typename DerivedB, typename DerivedC>
void path(const double globalCoupling,
        const double noiseScale,
        Eigen::MatrixBase<DerivedA>& initialConditions,
        Eigen::MatrixBase<DerivedB>& intrinsicFreqs,
        Eigen::MatrixBase<DerivedC>& adjacencyMatrix,
        const double rho,
        const double c,
        const double alpha,
        const double p,
        const double t,
        const int timeSteps,
        const int seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<> runif(0,1);
    std::exponential_distribution<> rexp(1);
    tempered_stable_distribution rtstable(alpha, 

    double dt = t/timeSteps;

    double gamma = c/(1-rho) - p;
    double cf = pow(cos(M_PI_2*rho),-1./rho);
    double sigma = pow(-dt*c*cos(M_PI*rho/2.)*tgamma(-rho), 1./rho);
    double mu = -dt*p;

    int N = initialConditions.size();

    // Initialize

    Vector theta(N);
    for (int i = 0; i < N; ++i)
        theta(i) = initialConditions(i);

    Vector r(timeSteps+1);
    r(0) = order_param(theta);

    // Time stepping

    double U, V, W, xi, drift;
    double s = 0.0;
    long k = 0;

    do 
    {
        s += dt;
        k++;

        for (int i = 0; i < N; i++)
        {
            // simulate a tempered stable random variable

            do 
            {
                U = runif(rng);
                V = runif(rng);
                do {
                    W = rexp(rng);
                } while (W == 0.);
                xi = mu + sigma*cf*pow(A(M_PI*U,rho)/pow(W,1.-rho),1./rho);
            } while (V > exp(-alpha*xi));

            // calculate the drift

            drift = intrinsicFreqs(i);
            for (int j = 0; j < N; j++)
                drift -= globalCoupling/N*adjacencyMatrix(i,j)*sin(theta(i) - theta(j));

            // increment process

            theta(i) += drift*dt + noiseScale*xi;

            theta(i) = mod2pi(theta(i));
        }

        r(k) = order_param(theta);

#ifdef DEBUG
        printf("%6.3f", s);
        printf(" %6.3f", order_param(theta));
        for (int i = 0; i < N; i++) {
            printf(" %6.3f", theta(i));
        }
        printf("\n");
#endif

    } while (s < t);

    //std::cout << r << std::endl;
}

template <typename DerivedA, typename DerivedB, typename DerivedC>
void paths(const double globalCoupling,
        const double noiseScale,
        Eigen::MatrixBase<DerivedA>& initialConditions,
        Eigen::MatrixBase<DerivedB>& intrinsicFreqs,
        Eigen::MatrixBase<DerivedC>& adjacencyMatrix,
        const double rho,
        const double c,
        const double alpha,
        const double p,
        const double t,
        const int timeSteps,
        const int seed,
        const int npaths) {

    std::mt19937 rng(seed);
    std::uniform_real_distribution<> runif(0,1);
    std::exponential_distribution<> rexp(1);

    double dt = t/timeSteps;

    double gamma = c/(1-rho) - p;
    double cf = pow(cos(M_PI_2*rho),-1./rho);
    double sigma = pow(-dt*c*cos(M_PI*rho/2.)*tgamma(-rho), 1./rho);
    double mu = -dt*p;

    int N = initialConditions.size();

    Vector avg_order(timeSteps+1);

    for (int j=0; j < npaths; ++j) {
        Vector theta(N);
        for (int i = 0; i < N; ++i)
            theta(i) = initialConditions(i);

        avg_order(0) += order_param(theta);

        double U, V, W, xi, drift;
        double s = 0.0;
        long k = 0;

        do {
            s += dt;
            k++;

            for (int i = 0; i < N; i++) {
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
                for (int j = 0; j < N; j++)
                    drift -= globalCoupling/N*adjacencyMatrix(i,j)*sin(theta(i) - theta(j));

                // increment process
                theta(i) += drift*dt + noiseScale*xi;
                theta(i) = mod2pi(theta(i));
            }

            avg_order(k) += order_param(theta);

        } while (s < t);
    }

    avg_order /= npaths;

    for (int i = 0; i < avg_order.size(); ++i) {
        printf("%6.3f %6.3f\n", i*dt, avg_order(i));
    }
}


int main(int argc, char const *argv[]) {
    double rho = 0.99; // stability index
    double c = 0.01;
    double alpha = 1.0;
    double xi = 0.2;
    double p = (1+xi)*(-c*rho*tgamma(-rho)*pow(alpha,rho-1));

    //printf("rho: %f c: %f alpha: %f xi: %f p: %f\n", rho, c, alpha, xi, p);

    double t = 30.0;
    int timeSteps = 5000;
    int seed = 1234;
    int npaths = 1000;

    double globalCoupling = 0.8;
    double noiseScale = 0.05;

    Vector initialConditions(3);
    initialConditions << 0.3, 6.0, 3.14;

    Vector intrinsicFreqs(3);
    intrinsicFreqs << 0.2, 0.2, 0.2;

    Matrix adjacencyMatrix(3,3);
    adjacencyMatrix <<
        0, 0, 1,
        0, 1, 0,
        0, 1, 1;

    path(globalCoupling, noiseScale, initialConditions, intrinsicFreqs,
            adjacencyMatrix, rho, c, alpha, p, t, timeSteps, seed);

    //paths(globalCoupling, noiseScale, initialConditions, intrinsicFreqs,
    //        adjacencyMatrix, rho, c, alpha, p, t, timeSteps, seed, npaths);

    return 0;
}
