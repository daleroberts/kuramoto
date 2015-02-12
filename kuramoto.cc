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
#include "process.h"

#if defined(_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_max_threads() { return 1;}
#endif

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

#pragma omp declare reduction (+ : Vector : omp_out += omp_in) \
                              initializer (omp_priv(omp_orig))

using namespace std;

inline double mod2pi(double theta) {
    while (theta > 2*M_PI)
        theta -= 2*M_PI;
    while (theta < 0.0)
        theta += 2*M_PI;
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
    
    return sqrt(sum_real*sum_real + sum_complex*sum_complex)/N;
}

template <typename DerivedA, typename DerivedB>
void paths(Graph& G,
           const Eigen::MatrixBase<DerivedA>& initial,
           const Eigen::MatrixBase<DerivedB>& freqs,
           const double alpha,
           const double a,
           const double b,
           const double K,
           const double max_t,
           const uint32_t nsteps,
           const uint32_t npaths,
           const uint32_t seed) {
    
    double dt = max_t/nsteps;
    auto N = initial.size();

    Vector avg_order = Vector::Zero(nsteps+1);
    
    #pragma omp parallel
    {
    int thread_seed = seed + omp_get_thread_num();
    mt19937 rng(thread_seed);
    TemperedStableDistribution rtstable(alpha, a, b, 5.0);
    
    #pragma omp for reduction(+:avg_order)
    for (size_t j=0; j < npaths; ++j) {
        Vector theta(N);
        for (size_t i = 0; i < N; ++i)
            theta(i) = initial(i);
        
        avg_order(0) += order_param(theta);
        
        double xi, drift;
        double s = 0.0;
        long k = 0;
        
        do {
            s += dt;
            k++;
            
            for (size_t i = 0; i < N; i++) {
                // simulate a tempered stable random variable
                xi = rtstable(rng);
                
                // calculate the drift
                drift = freqs(i);
                for (auto &j: G.neighbours(i))
                    drift -= K/N*sin(theta(i) - theta(j));
                
                // increment process
                theta(i) += drift*dt + xi;
                theta(i) = mod2pi(theta(i));
            }
            
            avg_order(k) += order_param(theta);
            
        } while (s < max_t);
    }

    }

    avg_order /= npaths;
    
    printf("{");
    for (size_t i = 0; i < avg_order.size()-1; ++i)
        printf("{%.6f, %.6f},", i*dt, avg_order(i));
    int i = avg_order.size() - 1;
    printf("{%.6f, %.6f}}\n", i*dt, avg_order(i));
}


int main(int argc, char const *argv[]) {
    int seed      = argc > 1 ? atol(argv[1]) : 1234;
    int npaths    = argc > 2 ? atoi(argv[2]) : 1000;
    int nsteps    = argc > 3 ? atoi(argv[3]) : 5000;
    double alpha  = argc > 4 ? atof(argv[4]) : 1.7; // stability
    double lambda = argc > 5 ? atof(argv[5]) : 1.0; // tempering
    double sigma  = argc > 6 ? atof(argv[6]) : 1.0; // diffusivity
    double K      = argc > 7 ? atof(argv[7]) : 0.8; // global coupling
    double max_t  = argc > 8 ? atof(argv[8]) : 30.0; // maximum time
    
    double a = 0.5*alpha*pow(sigma,2.0)/(tgamma(1-alpha)*cos(M_PI*alpha/2));
    double b = lambda;
    
    string filename = "graphs.g6";
    ifstream graphfile(filename);
    string line;
    while (getline(graphfile, line)) {
        Graph G(line);
        auto n = G.size();
        
        Vector initial = Vector::LinSpaced(n, -3.14, 3.14);
        Vector freqs = Vector::Zero(n);
        
        paths(G, initial, freqs, alpha, a, b, K, max_t, nsteps, npaths, seed);
    }
    
    return 0;
}
