/*

(c) 2020, Pericles Lopes Machado  <pericles [at] newavesoft.com>

Implementation of strongly connect components.
 
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

#include <vector>
#include <set>
#include <iostream>
#include <unordered_map>
#include <map>
#include <algorithm>

void dfs_order(const std::map<int64_t, std::vector<int64_t>>& G,
int64_t u,
std::vector<int64_t>& order,
std::map<int64_t, bool>& vis) {
      vis[u] = true;
     if (G.count(u) > 0) {
          for (const auto& v: G.at(u)) {
            if (!vis[v]){
                vis[v] = true;
                dfs_order(G, v, order, vis);
            }
          }
      }
      order.push_back(u);
}

void dfs_scc(const std::map<int64_t, std::vector<int64_t>>& G, int64_t u,
std::map<int64_t, int64_t>& vis, int64_t id) {
      vis[u] = id;
      if (G.count(u) > 0) {
          for (const auto& v: G.at(u)) {
            if (vis[v] != id && vis[v] == 0){
                vis[v] = id;
                dfs_scc(G, v, vis, id);
            }
          }
      }
}

int main() {
    int64_t nlines = 0;
    int64_t s=0, t=0;
    std::map<int64_t, std::vector<int64_t>> G, Gr;
    std::vector<int64_t> order;
    std::map<int64_t, bool> vis;
    std::map<int64_t, int64_t> vis2;
    
    std::set<int64_t> nodes;
    
    while (std::cin >> s >> t) {
        G[s].push_back(t);
        Gr[t].push_back(s);
        nodes.insert(s);
        nodes.insert(t);
        ++nlines;
    }
    
    std::cout << "nlines:" << nlines << "\n";
    
    for (auto it: nodes) {
        if (vis[it]) {
            continue;
        }
        dfs_order(Gr, it, order, vis);
    }

    int64_t nc = 1;
    std::reverse(order.begin(), order.end());
    for (auto u: order) {
        if (vis2[u] != 0) {
            continue;
        }
        dfs_scc(G, u,vis2, nc);
        ++nc;
    }
    
    std::map<int64_t, int64_t> scc;
    for (auto it: vis2) {
        ++scc[it.second];
    }
    
    std::vector<int64_t> scc_sizes;
    for (auto it: scc) {
        scc_sizes.push_back(-it.second);
    }
    std::sort(scc_sizes.begin(), scc_sizes.end());
    std::cout << "scc size: " << scc.size() << "\n";
    
    for (int64_t i = 0; i < 5; ++i) {
        if (i > 0) {
            std::cout << ",";
        }
        if (i < scc_sizes.size()) {
            std::cout << (-scc_sizes[i]);
        } else {
            std::cout << "0";
        }
    }    
    std::cout << "\n";
    
    return 0;
}


