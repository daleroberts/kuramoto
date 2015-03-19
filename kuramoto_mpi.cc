/*
 * Kuramoto Model with Tempered Stable Noise
 *
 * Dale Roberts <dale.o.roberts@gmail.com>
 */
#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdint>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <sstream>
#include <Eigen/Dense>
#include <boost/mpi.hpp>

#include "graph.h"
#include "variates.h"
#include "process.h"
#include "statistics.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

namespace mpi = boost::mpi;
using namespace std;

template <typename Iterable>
inline void elementwise_add(Iterable& in, Iterable& out) {
    auto in_iter = in.begin();
    for (auto& el : out)
        el += *(in_iter++);
}

inline double mod2pi(double theta) {
    while (theta > 2*M_PI)
        theta -= 2*M_PI;
    while (theta < 0.0)
        theta += 2*M_PI;
    return theta;
}

template <typename Derived>
inline double order_param(const Eigen::MatrixBase<Derived>& theta) {
    size_t N = theta.size();
    
    double sum_real = 0.0;
    double sum_complex = 0.0;
    for (size_t j = 0; j < N; j++) {
        sum_real += cos(theta(j));
        sum_complex += sin(theta(j));
    }
    
    return sqrt(sum_real*sum_real + sum_complex*sum_complex)/N;
}

template <typename Distribution>
vector<Statistics>
paths(Graph& G, Distribution& dist, const double alpha, const double a, const double b,
      const double K, const double max_t, const uint32_t nsteps,
      const uint32_t npaths, const uint32_t seed)
{
    double dt = max_t/nsteps;
    auto N = G.size();

    vector<Statistics> stats(nsteps+1);

    int thread_seed = seed + omp_get_thread_num();
    mt19937 rng(thread_seed);
    uniform_real_distribution<> runif(-M_PI, M_PI);
    
    for (size_t j=0; j < npaths; ++j) {
        Vector theta(N);

        // set initial condition
        for (size_t i = 0; i < N; ++i)
            theta(i) = mod2pi(runif(rng));
       
        stats[0].add(order_param(theta));
       
        double xi, drift;

        for (size_t k = 0; k <= nsteps; k++) {
            for (size_t i = 0; i < N; i++) {
                // simulate a tempered stable random variable
                xi = dist(rng);
               
                // calculate the drift
                drift = 0;
                for (auto &j: G.neighbours(i))
                    drift -= K/N*sin(theta(i) - theta(j));
               
                // increment process
                theta(i) += drift*dt + xi;
                theta(i) = mod2pi(theta(i));
            }
            stats[k].add(order_param(theta));
        }

        return stats;
    }
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
    string filename = "graphs.g6";
    
    double dt = max_t/nsteps;
    double a = 0.5*alpha*pow(sigma,2.0)/(tgamma(1-alpha)*cos(M_PI*alpha/2));
    double b = lambda;

    vector<Statistics> order_stats(nsteps+1);
    TemperedStableDistribution rtstable(alpha, dt*a, b, 1.1);
    NormalDistribution rnorm(0, dt*pow(sigma, 2));

    mpi::environment env;
    mpi::communicator world;

    seed += world.rank();
    npaths = npaths / world.size();
    
    ifstream graphfile(filename);
    string line;
    while (getline(graphfile, line)) {
        Graph G(line);
        vector<Statistics> stats;
        if (alpha > 1.999) {
            stats = paths(G, rnorm, alpha, a, b, K, max_t, nsteps, npaths, seed);
        } else {
            stats = paths(G, rtstable, alpha, a, b, K, max_t, nsteps, npaths, seed);
        }
        elementwise_add(stats, order_stats);
    }
    
    setbuf(stdout, NULL);
    printf("{");
    size_t i;
    for (i = 0; i < order_stats.size() - 1; ++i)
        printf("{%.6f, %.6f, %.6f, %.6f},", i*dt, order_stats[i].mean(),
               order_stats[i].stddev(), order_stats[i].variance());
    printf("{%.6f, %.6f, %.6f, %.6f}}\n", i*dt, order_stats[i].mean(),
           order_stats[i].stddev(), order_stats[i].variance());
    
    return 0;
}
