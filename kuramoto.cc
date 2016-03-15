/*
 * Kuramoto Model with Levy Noise
 *
 * Dale Roberts <dale.o.roberts@gmail.com>
 */

#define _USE_MATH_DEFINES

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <iomanip>
#include <random>
#include <string>
#include <vector>
#include <cmath>

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

inline double mod2pi(double theta) {
  while (theta > 2 * M_PI)
    theta -= 2 * M_PI;
  while (theta < 0.0)
    theta += 2 * M_PI;
  return theta;
}

inline double sum(const vector<double> &v) {
  double total = 0.0;
  for (auto &el : v)
    total += el;
  return total;
}

inline double order_param(const vector<double> &theta) {
  size_t N = theta.size();
  double r = 0.0;
  double c = 0.0;
  for (size_t j = 0; j < N; j++) {
    r += cos(theta[j]);
    c += sin(theta[j]);
  }
  return sqrt(r * r + c * c) / N;
}

inline double order_psi(const vector<double> &theta) {
  size_t N = theta.size();
  double r = 0.0;
  double c = 0.0;
  for (size_t j = 0; j < N; j++) {
    r += cos(theta[j]);
    c += sin(theta[j]);
  }
  r = r / N;
  c = c / N;
  return atan(c / r);
}

template <typename Distribution>
vector<Statistics> paths(Graph &G, Distribution &dist, const double alpha, const double a,
                         const double b, const double kappa, const double max_t,
                         const uint32_t nsteps, const uint32_t npaths, const uint32_t seed) {

  double dt = max_t / nsteps;
  auto N    = G.size();

  vector<Statistics> stats(nsteps + 1);

  // Initialise the random number generator. Each rank's PRNG is
  // initialised with its own seed

  mt19937 rng(seed);

  // Initialise the distribution for the initial distribution of the
  // vertices

  uniform_real_distribution<> runif(-M_PI, M_PI);

  for (size_t j = 0; j < npaths; ++j) {
    // Start simulating a path

    vector<double> theta(N, 0.);
    vector<double> theta_(N, 0.);
    vector<double> xi(N, 0.);
    double drift;

    // Set initial condition

    for (size_t i = 0; i < N; ++i)
      theta[i] = runif(rng);

    // Adjust initial condition to ensure zero empirical expectation

    theta[N - 1] = -(sum(theta) - theta[N - 1]);

    stats[0].add(order_param(theta));

    for (size_t k = 0; k <= nsteps; k++) {
      for (size_t i = 0; i < N; i++)
        xi[i] = dist(rng);

      // Adjust for zero empirical expectation of noise

      // xi[N - 1] = -(sum(xi) - xi[N - 1]);

      // Perform an Euler time stepping

      for (size_t i = 0; i < N; i++) {
        drift = 0;
        for (auto &j : G.neighbours(i))
          drift -= kappa * sin(theta[i] - theta[j]);

        theta_[i] += drift * dt + xi[i];
      }

      // Update values in one go

      theta = theta_;

      // Save result for this time step

      stats[k].add(order_param(theta));
    }
  }

  return stats;
}

int main(int argc, char *argv[]) {
  // Setup the MPI environment

  mpi::environment env(argc, argv, true);
  mpi::communicator world;

  // Parse the command line parameters or set defaults

  const char *graphfn;
  uint32_t ngraphs, seed, npaths, nsteps, rpaths, rseed;
  double alpha, lambda, sigma, max_t, dt, a, b, c, d, kappa;

  graphfn = argc > 1 ? argv[1] : "graphs.g6";
  ngraphs = argc > 2 ? atol(argv[2]) : 1;
  seed    = argc > 3 ? atol(argv[3]) : 0;
  npaths  = argc > 4 ? atoi(argv[4]) : 1;
  nsteps  = argc > 5 ? atoi(argv[5]) : 100;
  alpha   = argc > 6 ? atof(argv[6]) : 1.5;     // stability
  lambda  = argc > 7 ? atof(argv[7]) : 1.0;     // tempering
  sigma   = argc > 8 ? atof(argv[8]) : 5.78115; // diffusivity
  kappa   = argc > 9 ? atof(argv[9]) : 1.0;     // coupling
  max_t   = argc > 10 ? atof(argv[10]) : 1.0;

  // Derive distributional and time step parameters

  dt = max_t / nsteps;
  a  = 0.5 * alpha * pow(sigma, 2.0) / (tgamma(1 - alpha) * cos(M_PI * alpha / 2));
  b  = lambda;
  d  = tgamma(1 - a) * a * pow(b, alpha - 1);

  // Determine number of paths and seed for each rank

  rpaths = (uint32_t)npaths / world.size();
  rseed  = seed + world.rank();

  // Print out header before we do anything so that output file is
  // generated even if the simulation fails

  if (world.rank() == 0) {
    cout << "\"graphfn:" << graphfn;
    cout << " ngraphs:" << ngraphs;
    cout << " npaths:" << rpaths * world.size();
    cout << " nsteps:" << nsteps;
    cout << " alpha:" << alpha;
    cout << " lambda:" << lambda;
    cout << " sigma:" << sigma;
    cout << " kappa:" << kappa;
    cout << " seed:" << seed;
    cout << "\"" << endl;
  }

  // Save statistics for each time step

  vector<Statistics> stats(nsteps + 1);

  try {
    // Attempt to determine optimal c parameter. This can fail as it
    // may be impossible to find the quantile for a certain choice of
    // parameters

    c = -d*dt - qtstable(0.9, alpha, a, b);

    // Initialise the distributions

    TemperedStableDistribution rtstable(alpha, dt * a, b, c);
    StableDistribution rstable(alpha, dt * a);
    NormalDistribution rnorm(0, dt * pow(sigma, 2) / 2);

    // Open file containing graphs to simulate on

    uint32_t k = 0;
    ifstream infile(graphfn);
    string line;

    while (getline(infile, line)) {

      // Read in a graph

      Graph G(line);

      // Each rank calculates a portion of the paths and saves the
      // statistics at each time step

      vector<Statistics> rstats;

      // Depending on the choice of parameters, perform the simulation
      // with a different distribution

      if (alpha > 1.999) {
        // Gaussian noise
        rstats = paths(G, rnorm, alpha, a, b, kappa, max_t, nsteps, rpaths, rseed);
      } else {
        if (lambda < 0.001) {
          // Stable noise
          rstats = paths(G, rstable, alpha, a, b, kappa, max_t, nsteps, rpaths, rseed);
        } else {
          // Tempered stable noise
          rstats = paths(G, rtstable, alpha, a, b, kappa, max_t, nsteps, rpaths, rseed);
        }
      }

      // Aggregate the statistics together

      if (world.rank() == 0) {
        mpi::reduce(world, rstats, stats, std::plus<Statistics>(), 0);
      } else {
        mpi::reduce(world, rstats, std::plus<Statistics>(), 0);
      }

      // Stop performing simulations if we exceed the desired number
      // of graphs

      if (k++ >= ngraphs)
        break;
    }

  } catch (const std::exception &e) {

    // Abort the simulation if an error is caught

    env.abort(1);
  }

  // Write out the results

  if (world.rank() == 0) {
    cout << setiosflags(ios::fixed);
    for (size_t i = 0; i < stats.size(); ++i) {
      cout << setw(8) << setfill(' ') << setprecision(4);
      cout << i * dt << '\t' << stats[i] << endl;
    }
  }

  return 0;
}
