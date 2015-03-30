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
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

#include <boost/mpi.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/string.hpp>
#include <boost/lexical_cast.hpp>

#include "graph.h"
#include "variates.h"
#include "process.h"
#include "statistics.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

namespace mpi = boost::mpi;
using namespace std;

namespace boost { namespace mpi {
        template <>
        struct is_mpi_datatype<Statistics> : public mpl::true_ { };
    } }

std::ostream& operator<<(std::ostream& out, const Statistics& s) {
    return out << setprecision(6) << s.mean() << '\t' << s.stddev() << '\t' << s.variance();
}

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

    mt19937 rng(seed);
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
    }

    return stats;
}

int main(int argc, char const *argv[]) {
    const char *graphfile = argc > 1  ? argv[1]: "graphs.g6";
    uint32_t ngraphs      = argc > 2  ? atol(argv[2])  : 1;
    uint32_t seed         = argc > 3  ? atol(argv[3])  : 1234;
    uint32_t npaths       = argc > 4  ? atoi(argv[4])  : 1000;
    uint32_t nsteps       = argc > 5  ? atoi(argv[5])  : 5000;
    double alpha          = argc > 6  ? atof(argv[6])  : 1.7; // stability
    double lambda         = argc > 7  ? atof(argv[7])  : 1.0; // tempering
    double sigma          = argc > 8  ? atof(argv[8])  : 0.1; // diffusivity
    double K              = argc > 9  ? atof(argv[9])  : 0.8; // global coupling
    double max_t          = argc > 10 ? atof(argv[10]) : 30.0; // maximum time

    double dt = max_t/nsteps;
    double a = 0.5*alpha*pow(sigma,2.0)/(tgamma(1-alpha)*cos(M_PI*alpha/2));
    double b = lambda;

    vector<Statistics> order_stats(nsteps+1);
    TemperedStableDistribution rtstable(alpha, dt*a, b, 1.1);
    StableDistribution rstable(alpha, dt*a);
    NormalDistribution rnorm(0, dt*pow(sigma, 2));

    mpi::environment env;
    mpi::communicator world;

    seed += world.rank();
    int npaths_rank = npaths / world.size();

    ifstream infile(graphfile);
    string line;
    while (getline(infile, line))
    {
        Graph G(line);
        vector<Statistics> stats;
        K = (double) G.size();

        if (alpha > 1.999) {
            stats = paths(G, rnorm, alpha, a, b, K, max_t, nsteps, npaths_rank, seed);
        } else {
            if (lambda < 0.01) {
                stats = paths(G, rstable, alpha, a, b, K, max_t, nsteps, npaths_rank, seed);
            } else {
                stats = paths(G, rtstable, alpha, a, b, K, max_t, nsteps, npaths_rank, seed);
            }
        }

        if (world.rank() == 0) {
            mpi::reduce(world, stats, order_stats, std::plus<Statistics>(), 0);
        } else {
            mpi::reduce(world, stats, std::plus<Statistics>(), 0);
        }
    }

    if (world.rank() == 0) {
        cout << '"' << graphfile << " ngraphs:" << ngraphs << " npaths:"<< npaths_rank * world.size();
        cout << " nsteps:" << nsteps;
        cout << " alpha:" << alpha << " lambda:" << lambda << " sigma:" << sigma;
        cout << " K:" << K << " seed:" << seed << '"' << endl;
        cout << setiosflags(ios::fixed);
        for (size_t i = 0; i < order_stats.size(); ++i) {
            cout << setw(8) << setfill(' ') << setprecision(4) << i*dt << '\t';
            cout << order_stats[i] << endl;
        }
    }

    return 0;
}
