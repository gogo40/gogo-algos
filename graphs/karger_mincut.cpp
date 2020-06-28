/*

(c) 2020, Pericles Lopes Machado  <pericles [at] newavesoft.com>

Implementation of Karger algorithm to mincut.
 
  This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 
 */

#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <map>
#include <random>
#include <sstream>
#include <functional>
#include <omp.h>
#include <mutex>
#include <cmath>

namespace graph {
struct Node {
    int64_t head = 0;
    std::unordered_map<int64_t, int64_t> conns;
    bool is_active = true;
    
    void merge(std::vector<int64_t>& parent, 
               std::vector<Node>& g, const Node& other) const {
        int64_t t = other.head;   
        int64_t s = head;
        parent[head] = s;
        for (int64_t k = 0; k < parent.size(); ++k) {
            if (parent[k] == t) {
                parent[k] = s;
            }
        }
        parent[t] = s;
        
        Node result;
        result.head = this->head;
        for (auto conn: this->conns) {
            int64_t v = std::get<0>(conn), w = std::get<1>(conn);
            int64_t tail = parent[v];
            
            if (tail != this->head && tail != parent[this->head] &&
                tail != other.head && tail != parent[other.head]) {
                result.conns[tail] += w;
            }
        }
        for (auto conn: other.conns) {
            int64_t v = std::get<0>(conn), w = std::get<1>(conn);
            int64_t tail = parent[v];
            
            if (tail != this->head && tail != parent[this->head] &&
                tail != other.head && tail != parent[other.head]) {
                result.conns[tail] += w;
            }
        }
        g[head] = result;
    }
};

using Graph = std::vector<Node>;

Graph contract(std::vector<int64_t>& parent,
                           int64_t g1, int64_t g2,
                           const Graph& nodes, 
                           int64_t seed) {
    Graph g = nodes;
    std::set<int64_t> active_set;
    int64_t u = 0;

    for (const auto& node: nodes) {
        if (node.is_active) {
            active_set.insert(u);
        }
        parent[u] = u;
        ++u;
    }
    std::mt19937_64 gen(seed);
    while (active_set.size() > 2) {
        std::uniform_int_distribution<int64_t> U(0,  active_set.size()-1);
        auto it_s = active_set.begin();
        int64_t s = U(gen);
        std::advance(it_s, s);
        s = parent[*it_s];
        
        if (g[s].conns.size() < 1) {
            continue;
        }
        std::uniform_int_distribution<int64_t> U2(0, g[s].conns.size() -1);
        int64_t t = U2(gen);
        auto it_t = g[s].conns.begin();
        std::advance(it_t, t);
        t = parent[std::get<0>(*it_t)];
        if ((s==g1 && t==g2 ) ||
             (s==g2 && t==g1)) {
            continue;
        }
        g[s].merge(parent, g, g[t]);
        parent[s] = s;
        
        g[t].is_active = false;
        active_set.erase(t);
    }
    for (auto& node: g) {
        node.head = parent[node.head];
    }
    return g;
}

int64_t mincut(const Graph& nodes,
               int64_t g1, int64_t g2,
               int64_t n_iterations, int64_t seed,
               int n_threads) {
    int64_t best = 1<<20;
    std::mutex m;
    std::mt19937_64 gen(seed);
#pragma omp parallel for num_threads(n_threads)
    for  (int64_t i =0; i < n_iterations; ++i) {
        int64_t s = 0;
        {
            std::lock_guard<std::mutex> _(m);
            s = gen() + 2 * i + 1;
        }
        std::vector<int64_t> parent(nodes.size());
        auto g = contract(parent, g1, g2, nodes, s);
        int64_t r = 0;
        for (const auto& node: g) {
            if (node.is_active && parent[node.head] == node.head) {
                int64_t conn_active = 0;
                for (auto conn: node.conns) {
                    auto  v = parent[std::get<0>(conn)];
                    int64_t w = std::get<1>(conn);
                    if (v != parent[node.head]) {
                        //++conn_active;
                        conn_active += w;
                    }
                }
                if (conn_active > 0){
                    r = std::max(conn_active,r);
                }
            }
        }
        {
            std::lock_guard<std::mutex> _(m);
            if (r < best) {
                best = r;
            }
        }
    }

    return best;
}

int64_t mincut(const Graph& nodes, int64_t ntries, int64_t n_threads) {
    std::mt19937_64 gen;
    int64_t N = nodes.size();
    int64_t mm = 1<<20;
    std::uniform_int_distribution<int64_t> U(1, N+ 1);
    for (int64_t i = 0; i <  ntries; ++i) {
        int64_t s  = U(gen);
        int64_t t = U(gen);
        mm = std::min(graph::mincut(nodes, s, t, 10, gen(), n_threads), mm);
        std::cerr<< s << " --> "<< t <<  "  | min_cut: " << mm << "\n";
    } 
    return mm;
}

} // namespace graph

int main() {
    int64_t N;
    {
        std::string line;
        std::getline(std::cin, line);
        std::istringstream in(line);
        in >> N;
    }
    
    graph::Graph nodes(N+1);
    nodes[0].is_active = false;
    int64_t M = 0;
    for (int64_t i = 0; i < N; ++i) {
        std::string line;
        std::getline(std::cin, line);
        std::istringstream in(line);
        int64_t v, u;
        in>>v;
        nodes[v].head = v;
        nodes[v].is_active = true;
        while (in>>u) {
            ++M;
            ++nodes[v].conns[u];
        }        
    }
    
   
    int ntries = N - 1;
    int n_threads = 8;
   
    std::cout << "min_cut:" << graph::mincut(nodes, ntries,  n_threads) << "\n";
    return 0;
}
