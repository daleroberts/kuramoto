#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <Eigen/Dense>

using namespace std;

#include "graph.hpp"

using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;

int main(int argc, char *argv[]) {
  string filename = "test.g6";
  ifstream graphfile(filename);
  cout << "Parsing " << filename << endl;
  vector<Graph> graphs;
  string line;
  while (getline(graphfile, line)) {
    graphs.push_back(Graph(line));
  }

  for (auto &g: graphs) {
    auto n = g.size();
    Matrix A = Matrix::Zero(n,n);
    for (size_t i = 0; i < n; ++i) {
        for (auto &j : g.neighbours(i))
            A(i,j) = 1;
    }
    cout << A << endl;
    cout << "----" << endl;
   }
}
