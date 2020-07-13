#include <iostream>
#include <map>
#include <queue>
#include <unordered_map>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <mutex>
#include <set>
#include <random>
#include <memory>
#include <algorithm>
#include <omp.h>

using namespace std;


struct State {
    State(const std::map<int, std::set<int>>& blocked_nodes = std::map<int, std::set<int>>(),
    const std::vector<int>& test_colors = std::vector<int>(),
    bool is_ok = true) : blocked_nodes(blocked_nodes), test_colors(test_colors), is_ok(is_ok) {}
    State(const State& other) : blocked_nodes(other.blocked_nodes), test_colors(other.test_colors), is_ok(other.is_ok) {}
    
    std::map<int, std::set<int>> blocked_nodes;
    std::vector<int> test_colors;
    bool is_ok = true;
  
    void copy(const std::map<int, std::set<int>>& blocked_nodes = std::map<int, std::set<int>>(),
    const std::vector<int>& test_colors = std::vector<int>()) {
    
        this->test_colors.reserve(test_colors.size());
        for (auto v: test_colors) {
            this->test_colors.push_back(v);
        }
        
        for (auto v: blocked_nodes) {
            auto& o = this->blocked_nodes[v.first];
            //o.reserve(v.second.size());
            for (auto p: v.second) {
                o.insert(p);
            }
        }
    }  
};
        

class Graph {
public:
    Graph(int64_t N, int64_t E) 
    : N(N), E(E), e(N), adj(N), conn(N), D(N) {
        for (auto& v: adj) {
            v = std::vector<bool>(N, false);
        }
    }   
    
    Graph(const Graph& other) 
    : r(other.r), N(other.N), E(other.E), 
    e(other.e), adj(other.adj), 
    conn(other.conn), D(other.D) {}
    
    void add_edge(int u, int v) {
        if (!adj[u][v]) {
            e[u].push_back(v);
            e[v].push_back(u);
            adj[u][v] = adj[v][u] = true;
            conn[u].insert(v);
            conn[v].insert(u);
            _edges.emplace_back(u, v);
        }
    }
    
    const std::vector<int>& conns(int u) const {
        return e[u];
    }
    
    int degree(int u) const {
        return e[u].size();
    }
    
    bool has_conn(int u, int v) const {
        return adj[u][v];
    }

    void update_distances() {
        r = -1;
        D.resize(N);
        for (int i = 0; i < N; ++i) {
            D[i].resize(N);
        }
        for (int i = 0; i < N; ++i) {
            for (int v = 0; v < N; ++v) {
                D[i][v] = -1;
            }
            
            auto& d = D[i];
            d[i] = 0;
            std::queue<int> Q;
            Q.push(i);
            
            int max_r = 0;
            while (!Q.empty()) {
                int u = Q.front();
                Q.pop();
                for (auto v: e[u]) {
                    if (d[v] < 0 || d[v] > d[u] + 1) {
                        d[v] = d[u] + 1;
                        Q.push(v);
                        if (d[v] > max_r) {
                            max_r = d[v];
                        }
                    }
                }
            }
            if (r < 0 || max_r < r) {
                max_r = r;
                ic = i;
            }
        }
    }
    
    int dist(int u, int v) const {
        return D[u][v];
    }
    
    int center() const {
        return ic;
    }
    
    int ratio() const {
        return r;
    }
    
    int n_nodes() const {
        return N;
    }
    
    int n_edges() const {
        return E;
    }
    
    const std::tuple<int, int> edge(int i) const {
        return _edges[i];
    }
    
    const std::vector<std::tuple<int, int>>& edges() const {
        return _edges;
    }
    
private:
    int64_t r;
    int64_t ic;
    int64_t N, E;
    std::vector<std::tuple<int, int>> _edges;
    std::vector<std::vector<int>> e;
    std::vector<std::vector<bool>> adj;
    std::vector<std::set<int>> conn;
    std::vector<std::vector<int>> D;
};

int main(int argc, char** argv) 
{
    //std::string _(argv[1]);
    auto& in = cin;
    //std::istringstream in(_);
    int64_t N, E;
    in >> N >> E;
    
    
    
    Graph G(N, E);
    for (int64_t i = 0; i < E; ++i) {
        int64_t u, v;
        in >> u >> v;
        G.add_edge(u, v);
    }
    G.update_distances();
    
    // identify 3-clique, 4-clique, 5-clique ...
    
    int n_cliques = 200;
    std::vector<std::vector<std::vector<int>>> clique(n_cliques);
    
    for (const auto& e: G.edges()) {
        int u1 = std::get<0>(e), u2 = std::get<1>(e);
        if (u1 < u2) {
            clique[0].push_back(std::vector<int>{u1, u2});
        }        
    }

    
    int max_clique = 2;
    int max_n_cliques_by_level = 100000;
    for (int ic = 1; ic < n_cliques; ++ic) {
        std::vector<int> in_clique(N+1,-1);
        int it = 0;
        for (const auto& e: clique[ic-1]) {
            std::vector<int> nei;
            for (auto u: e) {
                in_clique[u] = it;
            }
            
            for (auto u: e) {
                const auto& c = G.conns(u);
                for (auto v: c) {
                    if (in_clique[v] != ic) {
                        nei.push_back(v);
                        in_clique[v] = ic;
                    }
                }
            }
            
            for (auto p: nei) {
                bool ok = true;
                for (auto v: e) {
                    if (!G.has_conn(p, v) || p <= v) {
                        ok = false;
                    }
                }
                if (ok) {
                    std::vector<int> c;
                    for (auto v: e) {
                        c.push_back(v);
                    }
                    c.push_back(p);
                    clique[ic].push_back(c);
                }
            }
            if (clique[ic].size() > max_n_cliques_by_level) {
                break;
            }
            ++it;
        }
        if (clique[ic].size() > 0) {
            max_clique = ic + 2;
            std::cerr << "|" << (ic + 2) << "-clique| = " << clique[ic].size() << "\n";
        } else if (clique[ic-1].size() < 1) {
            break;
        } else if (ic > 2) {
            clique[ic-2].clear();
        }
    }
    
    bool ok_clique = true;
    
    std::vector<int> density(G.n_nodes(), 0), az(G.n_nodes(), 0);
    
    for (int64_t i = max_clique - 2; i > -1; --i) {
        for (const auto& c: clique[i]) {
            for (auto u: c) {
                if (density[u] <= 0) {
                    az[u] = density[u] = i + 2;
                    az[u] += G.conns(u).size();
                }
            }
        }
    }
    
    std::vector<std::tuple<int, int, int, int>> node_values;
    for (int i = 0; i < density.size(); ++i) {
        node_values.emplace_back(-density[i], -G.degree(i), -az[i],  i); 
    }
    std::sort(node_values.begin(), node_values.end());
    std::vector<int> nodes;
    for (int i = 0; i < node_values.size(); ++i) {
        nodes.push_back(std::get<3>(node_values[i]));
    }
    
    std::cerr << "max clique: " << max_clique << "\n";
    std::cerr << "# nodes: " << G.n_nodes() << "\n";
    std::cerr << "# edges: " << G.edges().size() << "\n";    
    std::cerr << "\n";
    
    
    std::vector<int> colors(G.n_nodes(), -1);
    
    int64_t color = 0;
    
    {
        std::map<int, std::set<int>> blocked_nodes;
        
        for (auto u: clique[max_clique - 2][0]) {
            colors[u] = color;
            for (auto v: G.conns(u)) {
                blocked_nodes[color].insert(v);
            }
            ++color;
        }
        
        // create base line using greedy
        for (int u: nodes) {
            if (colors[u] > -1) {
                continue;
            }
            
            bool ok = false;
            int selected_color = -1;
            for (int c = 0; c < color; ++c) {
                if (blocked_nodes[c].find(u) == blocked_nodes[c].end()) {
                    selected_color = c;
                    ok = true;
                    break;
                }
            }
            
            if (!ok) {
                selected_color = color;
                ++color;
            }
            
            colors[u] = selected_color;
            
            for (auto v: G.conns(u)) {
                blocked_nodes[selected_color].insert(v);
            }
        }
    }
    
    int n_iterations = 700;
    int n_threads = 8;
    std::mutex m;
    std::vector<std::default_random_engine> generators(n_threads);
    
    
    #pragma omp parallel for num_threads(n_threads)
    for (int i = 0; i < n_iterations; ++i) {
        int test_color = 0;
        int tid = omp_get_thread_num();
        auto& generator = generators[tid];

        {
            std::lock_guard<std::mutex> _(m);
            std::uniform_int_distribution<int> dec(color-n_threads,color-1);
            test_color = (tid < 0.5*n_threads)? (color - 1) : (dec(generator));
        }
        
        bool ok = true;
        int current_c = 0;
        std::priority_queue<std::tuple<int, int, int, int, int>> Q;
        
        std::uniform_int_distribution<int> cs(0, clique[max_clique - 2].size() - 1);
        int cid = cs(generator);
        
        int last_state = 1, max_n_states = 600;
        std::map<int, std::shared_ptr<State>> states;
        
        
        {
            std::vector<int> test_colors(G.n_nodes(), -1);
            std::map<int, std::set<int>> blocked_nodes;        
            std::set<int> nei;
            for (auto u: clique[max_clique - 2][cid]) {
                test_colors[u] = current_c;
                for (auto v: G.conns(u)) {
                    blocked_nodes[current_c].insert(v);
                    nei.insert(v);
                }
                ++current_c;
            }
                        
            for (auto u: nei) {
                if (test_colors[u] > -1) {
                    continue;
                }
                int n_opts = 0, last_opt=-1;
                for (int c = 0; c < test_color; ++c) {
                    if (blocked_nodes[c].find(u) == blocked_nodes[c].end()) {
                        ++n_opts;
                        last_opt = c;
                    }
                }
                if (n_opts == 0) {
                    ok = false;
                    break;
                }
                Q.push(std::make_tuple(-n_opts, density[u],  az[u],  0, -u));            
            }
            states[0] = std::make_shared< State>(
                blocked_nodes,
                test_colors,
                ok
            );
        }
        
        
        int sol_state = 0;
        bool is_solved = false;
        int max_iter = std::min(200000, G.n_nodes() * 1000);
        
        while (!Q.empty() && !is_solved && ok && max_iter > 0) {
            int is = -std::get<3>(Q.top());
            int u = -std::get<4>(Q.top());
            --max_iter;
            Q.pop();
            
            if (states.find(is) == states.end()) {
                continue;
            }
            auto& test_colors = states[is]->test_colors;
            auto& blocked_nodes = states[is]->blocked_nodes;
            
            is_solved = true;
            for (int u = 0; u < G.n_nodes(); ++u) {
                if (test_colors[u] < 0) {
                    is_solved = false;
                    break;
                }
            }
            
            if (is_solved) {
                ok = true;
                sol_state = is;
                break;
            }
            
            ok = false;
            for (auto its: states) {
                auto i = its.first;
                auto& s= its.second;
                if (s == nullptr) {
                    continue;
                }
                if (s->is_ok) {
                    ok = true;
                    break;
                }
            }
            
            if (!ok) {
                break;
            }

            
            if (test_colors[u] > -1) {
                continue;
            }
            
            std::vector<int> cands;
            int n_opts = 0;
            for (int c = 0; c < test_color; ++c) {
                if (blocked_nodes[c].find(u) == blocked_nodes[c].end()) {
                    ++n_opts;
                    cands.push_back(c);
                }
            }
            if (cands.size() < 1 && test_colors[u] < 0) {
                states[is]->is_ok = false;
                //states.erase(is);
                //states[is]=nullptr;
                continue;
            }
            
            std::uniform_int_distribution<int> U(0, cands.size() -1);
            int i = U(generator);
            int selected_color = cands[i];
            test_colors[u] = selected_color;
            for (auto v: G.conns(u)) {
                blocked_nodes[selected_color].insert(v);
            }
            
            for (auto v: G.conns(u)) {
                if (test_colors[v] > -1) {
                    continue;
                }
                n_opts = 0;
                int last_opt = -1;
                std::vector<int> opts;
                for (int c = 0; c < test_color; ++c) {
                    if (blocked_nodes[c].find(v) == blocked_nodes[c].end()) {
                        ++n_opts;
                        opts.push_back(c);
                        last_opt = c;
                    }
                }
                if (n_opts == 0) {
                    states[is]->is_ok = false;
                    //states.erase(is);
                    break;
                }
                
                Q.push(std::make_tuple(-n_opts, density[v], az[v], -is, -v));
                    
                if (n_opts < 3 && states.size() < max_n_states) {
                    int js = last_state;
                    for (auto c: opts) {
                        auto s = std::make_shared<State>();
                        s->copy(blocked_nodes, test_colors);
                        s->test_colors[v] = c;
                        s->is_ok = true;                          
                        for (auto p: G.conns(v)) {
                            if (s->test_colors[p] > -1) {
                                continue;
                            }
                            s->blocked_nodes[c].insert(p);
                            n_opts = 0;
                            for (int c2=0; c2 < test_color; ++c2) {
                                if (s->blocked_nodes[c2].find(p) == s->blocked_nodes[c2].end()) {
                                    ++n_opts;
                                }
                            }
                            if (n_opts > 0) {
                                Q.push(std::make_tuple(-n_opts, density[p],   az[p], -js, -p));
                                states[js] = s;
                                ++js;
                            }                                      
                         }
                    }
                    last_state = js;
                }   
            }        
        }
        
        if (i % 1500 == 0) {
            std::lock_guard<std::mutex> _(m);
            std::cerr << "iteration " << (i+1) << "/" << n_iterations << "\n";
        }
        
        if (ok && is_solved) {
            std::lock_guard<std::mutex> _(m);
            auto& test_colors = states[sol_state]->test_colors;
            if (test_color < color) {
                colors = test_colors;
                color = test_color;
                std::cerr << "sol state:" << sol_state << "\n";
                std::cerr << "new # color: " << color << " " << (i + 1) << "/" << n_iterations << "\n";
                for (int i = 0; i < colors.size(); ++i) {
                    if (i > 0) {
                        std::cerr << " ";
                    }
                    std::cerr << colors[i];
                }
                std::cerr << "\n";   
            }
        }   
    }
    
    std::cerr << "# colors: " << color << "\n";
    std::cout << color << " 0\n";
    
    
    std::set<int> cids;
    for (int i = 0; i < colors.size(); ++i) {
        cids.insert(colors[i]);
    } 
    std::cerr << "# colors ids: " << cids.size() << "\n";   
    for (int i = 0; i < colors.size(); ++i) {
        if (i > 0) {
            std::cout << " ";
        }
        std::cout << colors[i];
    }
    std::cout << "\n";
    
    
    return 0;
}


