#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <string>

#include "graph.h"

using namespace std;

template <typename Container> void sort(Container &c) { std::sort(c.begin(), c.end()); }

Graph::Graph(const string &str) {
  string s(str);

  if (s.find(">>graph6<<") != string::npos)
    s = s.substr(10);
  s.erase(s.find_last_not_of("\n\r") + 1);

  vector<int> data = graph6_to_data(s);

  size_t n;
  if (data[0] <= 62) {
    n = data[0];
    data.erase(data.begin());
  } else if (data[1] <= 62) {
    n = (data[1] << 12) + (data[2] << 6) + data[3];
    data.erase(data.begin(), data.begin() + 4);
  } else {
    n = (data[2] << 30) + (data[3] << 24) + (data[4] << 18) + (data[5] << 12) + (data[6] << 6) +
        data[7];
    data.erase(data.begin(), data.begin() + 8);
  }
  _n = n;

  size_t nd = (n * (n - 1) / 2 + 5) / 6;
  if (data.size() != nd)
    cerr << "error parsing graph." << endl;

  vector<size_t> nodes;
  for (size_t i = 0; i < n; ++i)
    nodes.push_back(i);

  size_t k = 5;
  auto d = data.begin();
  for (size_t j = 1; j < n; ++j) {
    for (size_t i = 0; i < j; ++i) {
      if ((*d >> k) & 1)
        add_edge(i, j);
      if (!k--) {
        k = 5;
        d++;
      }
    }
  }
}

size_t Graph::size() { return _n; }

vector<node_t> Graph::nodes(void) {
  vector<node_t> v;
  for (auto &node : _adj)
    v.push_back(node.first);
  sort(v);
  return v;
}

vector<node_t> Graph::neighbours(const node_t i) { return _adj[i]; }

void Graph::add_nodes_from(vector<node_t> nodes) {
  for (auto &n : nodes) {
    if (!has_node(n)) {
      _adj[n];
    }
  }
}

void Graph::add_edge(const node_t i, const node_t j) {
  if (!has_node(i)) {
    _adj[i];
  }
  if (!has_node(j)) {
    _adj[j];
  }
  _adj[i].push_back(j);
  _adj[j].push_back(i);
}

bool Graph::has_node(const node_t i) { return _adj.count(i) > 0; }

size_t Graph::max_degree(void) {
  size_t maximum = 1;
  for (auto &kv : _adj)
    maximum = std::max(kv.second.size(), maximum);
  return maximum;
}

vector<int> Graph::graph6_to_data(string s) {
  vector<int> v;
  for (auto &c : s)
    v.push_back(int(c) - 63);
  auto max = *max_element(v.begin(), v.end());
  auto min = *min_element(v.begin(), v.end());
  if (v.size() > 0 && (min < 0 || max > 63)) {
    v.clear();
    cerr << "error parsing graph6 file." << endl;
  }
  return v;
}
