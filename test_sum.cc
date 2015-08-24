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

int main(int argc, char const* argv[]) {
    setbuf(stdout, NULL);
    
    const char* graphfile = argc > 1 ? argv[1] : "graphs.g6";

    mt19937 rng(12312);
    uniform_real_distribution<> runif(-10*M_PI, 10*M_PI);
    StableDistribution rstable(1.6, 1.06);

    double total, drift, K;
    size_t N;
   
    ifstream infile(graphfile);
    string line;
    while (getline(infile, line)) {
        Graph G(line);

        N = G.size();
        K = G.size();
        
        Vector theta(N);
        for (size_t i = 0; i < N; ++i)
            theta(i) = rstable(rng);
         
        for (size_t k = 0; k < 10; k++) {
            total = 0.0;
            for (size_t i = 0; i < N; i++) {
                drift = 0.0;
                for (auto& j: G.neighbours(i))
                    drift -= K/N * sin(theta(i) - theta(j));
                total += drift;
                theta(i) += drift * 0.01;
            }
        }

        printf("%.8f\n", total);
    }

    return 0;
}
