#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>

typedef uint32_t node_t;

class Graph {
  public:
    Graph(const std::string &s);

    size_t size();
    std::vector<node_t> nodes();
    std::vector<node_t> neighbours(const node_t i);
    void add_nodes_from(std::vector<node_t> nodes);
    void add_edge(const node_t i, const node_t j);
    bool has_node(const node_t i);

  private:
    typedef std::vector<node_t> nodelist;
    std::unordered_map<node_t, nodelist> _adj;
    size_t _n;

    std::vector<int> graph6_to_data(std::string s);
};
