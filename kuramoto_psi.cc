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

#include <boost/mpi.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/string.hpp>
#include <boost/lexical_cast.hpp>

#include "graph.h"
#include "variates.h"
#include "statistics.h"

namespace boost {
namespace mpi {
template <> struct is_mpi_datatype<Statistics> : public mpl::true_ {};
}
}

namespace mpi = boost::mpi;
using namespace std;

std::ostream &operator<<(std::ostream &out, const Statistics &s) {
  return out << setprecision(6) << s.mean() << '\t' << s.stddev();
}

template <typename Iterable> inline void elementwise_add(Iterable &in, Iterable &out) {
  auto in_iter = in.begin();
  for (auto &el : out)
    el += *(in_iter++);
}

inline double sum(const vector<double> &v) {
  double total = 0.0;
  for (auto &el : v)
    total += el;
  return total;
}

inline double mod2pi(double theta) {
  while (theta > 2 * M_PI)
    theta -= 2 * M_PI;
  while (theta < 0.0)
    theta += 2 * M_PI;
  return theta;
}

inline double order_param(const vector<double> &theta) {
  size_t N           = theta.size();
  double sum_real    = 0.0;
  double sum_complex = 0.0;
  for (size_t j = 0; j < N; j++) {
    sum_real += cos(theta[j]);
    sum_complex += sin(theta[j]);
  }
  return sqrt(sum_real * sum_real + sum_complex * sum_complex) / N;
}

inline double order_psi(const vector<double> &theta) {
  size_t N           = theta.size();
  double sum_real    = 0.0;
  double sum_complex = 0.0;
  for (size_t j = 0; j < N; j++) {
    sum_real += cos(theta[j]);
    sum_complex += sin(theta[j]);
  }
  sum_real    = sum_real / N;
  sum_complex = sum_complex / N;

  return atan2(sum_complex, sum_real);
}

template <typename Distribution>
vector<Statistics> paths(Graph &G, Distribution &dist, const double alpha, const double a,
                         const double b, const double kappa, const double max_t,
                         const uint32_t nsteps, const uint32_t npaths, const uint32_t seed) {
  double dt = max_t / nsteps;
  auto N    = G.size();

  vector<Statistics> stats(nsteps + 1);

  mt19937 rng(seed);
  uniform_real_distribution<> runif(-M_PI, M_PI);

  vector<double> theta(N, 0.);
  vector<double> theta_(N, 0.);
  vector<double> xi(N, 0.);

  for (size_t j = 0; j < npaths; ++j) {

    // set initial condition
    for (size_t i = 0; i < N; ++i)
      theta[i] = runif(rng);

    theta[N - 1] = -(sum(theta) - theta[N - 1]); // adjust for zero expectation

    stats[0].add(order_psi(theta));

    double drift;

    for (size_t k = 0; k <= nsteps; k++) {

      for (size_t i = 0; i < N; i++)
        xi[i] = dist(rng);

      xi[N - 1] = -(sum(xi) - xi[N - 1]); // adjust for zero expectation

      for (size_t i = 0; i < N; i++) {
        drift = 0;
        for (auto &j : G.neighbours(i))
          drift -= kappa * sin(theta[i] - theta[j]);

        theta_[i] += drift * dt + xi[i];
      }
      theta = theta_;
      stats[k].add(order_psi(theta));
    }
  }

  return stats;
}

int main(int argc, char const *argv[]) {
  const char *graphfile = argc > 1 ? argv[1] : "graphs.g6";
  uint32_t ngraphs      = argc > 2 ? atol(argv[2]) : 1;
  uint32_t seed         = argc > 3 ? atol(argv[3]) : 0;
  uint32_t npaths       = argc > 4 ? atoi(argv[4]) : 10;
  uint32_t nsteps       = argc > 5 ? atoi(argv[5]) : 5;
  double alpha          = argc > 6 ? atof(argv[6]) : 1.5;     // stability
  double lambda         = argc > 7 ? atof(argv[7]) : 1.0;     // tempering
  double sigma          = argc > 8 ? atof(argv[8]) : 5.78115; // diffusivity
  double K              = argc > 9 ? atof(argv[9]) : 1.0;     // coupling
  double max_t          = argc > 10 ? atof(argv[10]) : 2.0;

  double dt = max_t / nsteps;
  double a  = 0.5 * alpha * pow(sigma, 2.0) / (tgamma(1 - alpha) * cos(M_PI * alpha / 2));
  double b  = lambda;
  double d  = dt * tgamma(1 - a) * a * pow(b, alpha - 1);
  double c  = -d - qtstable(0.9, alpha, a, b);

  vector<Statistics> global_stats(nsteps + 1);

  TemperedStableDistribution rtstable(alpha, dt * a, b, c);
  StableDistribution rstable(alpha, dt * a);
  NormalDistribution rnorm(0, dt * pow(sigma, 2) / 2); // scaled!

  mpi::environment env;
  mpi::communicator world;
  int npaths_rank;

  npaths_rank = npaths / world.size();
  seed        = seed + world.rank();
  ngraphs     = 0;

  bool first_graph = true;
  double kappa;

  ifstream infile(graphfile);
  string line;
  while (getline(infile, line)) {
    Graph G(line);
    vector<Statistics> stats;
    ngraphs++;

    if (first_graph) {
      K           = (double)G.max_degree();
      first_graph = false;
    }

    kappa = 1.0;

    if (alpha > 1.999) {
      // Gaussian noise
      stats = paths(G, rnorm, alpha, a, b, kappa, max_t, nsteps, npaths_rank, seed);
    } else {
      if (lambda < 0.001) {
        // Stable noise
        stats = paths(G, rstable, alpha, a, b, kappa, max_t, nsteps, npaths_rank, seed);
      } else {
        // Tempered stable noise
        stats = paths(G, rtstable, alpha, a, b, kappa, max_t, nsteps, npaths_rank, seed);
      }
    }

    if (world.rank() == 0) {
      mpi::reduce(world, stats, global_stats, std::plus<Statistics>(), 0);
    } else {
      mpi::reduce(world, global_stats, std::plus<Statistics>(), 0);
    }
  }

  if (world.rank() == 0) {
    // header
    cout << '"' << graphfile;
    cout << " ngraphs:" << ngraphs;
    cout << " npaths:" << npaths_rank * world.size();
    cout << " nsteps:" << nsteps;
    cout << " alpha:" << alpha;
    cout << " lambda:" << lambda;
    cout << " sigma:" << sigma;
    cout << " kappa:" << kappa;
    cout << " seed:" << seed;
    cout << '"' << endl;

    // rows
    cout << setiosflags(ios::fixed);
    for (size_t i = 0; i < global_stats.size(); ++i) {
      cout << setw(8) << setfill(' ') << setprecision(4);
      cout << i * dt << '\t' << global_stats[i] << endl;
    }
  }

  return 0;
}