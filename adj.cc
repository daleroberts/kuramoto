#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "graph.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  ifstream graphfile("graphs.g5");
  vector<Graph> graphs;
  string line;
  while (getline(graphfile, line)) {
    graphs.push_back(Graph(line));
  }

  for (auto &g: graphs) {
    cout << g.size() << endl;
  }
}
