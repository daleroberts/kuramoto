#ifndef GRAPH_H
#define GRAPH_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

class Graph {
public:
    Graph(const std::string& str) {
        std::string s(str);
        if (s.find(">>graph6<<") != std::string::npos)
            s = s.substr(10);

        auto data = graph6_to_data(s);

        size_t n;
        if (data[0] <= 62) {
            n = data[0];
            data.erase(data.begin());
        } else if (data[1] <= 62) {
            n = (data[1]<<12) + (data[2]<<6) + data[3];
            data.erase(data.begin(), data.begin() + 4);
        } else {
            n = (data[2]<<30) + (data[3]<<24) + (data[4]<<18)
                + (data[5]<<12) + (data[6]<<6) + data[7];
            data.erase(data.begin(), data.begin() + 8);
        }
        n_ = n;

        size_t nd = (n*(n-1)/2 + 5) / 6;
        if (data.size() != nd)
            std::cout << "error parsing graph: " << s << std::endl;

        std::vector<size_t> nodes;
        for (size_t i = 0; i < n; ++i)
            nodes.push_back(i);

        size_t k = 5;
        auto d = data.begin();
        for (size_t j = 1; j < n; ++j) {
            for (size_t i = 0; i < j; ++i) {
                if ((*d>>k)&1)
                    add_edge(i, j);
                if (!k--) {
                    k = 5;
                    d++;
                }
            }
        }
    }

    size_t size() {
        return n_;
    }

    void add_nodes_from(std::vector<size_t> nodes) {
        for (auto &n : nodes) {
            if (!has_node(n)) {
                nodemap m;
                adj_[n] = m;
            }
        }
    }

    std::vector<size_t> nodes(void) {
        std::vector<size_t> v;
        for (auto &el : nodes_)
            v.push_back(el.first);
        std::sort(v.begin(), v.end());
        return v;
    }

    void add_edge(const size_t i, const size_t j) {
        if (!has_node(i)) {
            adj_[i];
            nodes_[i];
        }
        if (!has_node(j)) {
            adj_[j];
            nodes_[j];
        }
        adj_[i][j] = true; 
        adj_[j][i] = true; 
    }

    bool has_node(const size_t i) {
        return nodes_.count(i) > 0;
    }

    std::vector<size_t> neighbours(const size_t i) {
        std::vector<size_t> v;
        for (auto &a : adj_[i])
                v.push_back(a.first);
        std::sort(v.begin(), v.end());
        return v;
    }

private:
    typedef std::unordered_map<size_t, bool> nodemap;
    nodemap nodes_;
    std::unordered_map<size_t, nodemap> adj_;
    size_t n_;

    std::vector<int> graph6_to_data(std::string s) {
        std::vector<int> v;
        for (auto &c : s)
            v.push_back((int)c-63);
        int max = *std::max_element(v.begin(), v.end());
        int min = *std::min_element(v.begin(), v.end());
        if (v.size() > 0 && (min < 0 || max > 63))
            v.clear();
        return v;
    }
};

#endif
