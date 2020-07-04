#include <iostream>
#include <map>
#include <queue>
#include <unordered_map>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <set>
#include <random>
#include <algorithm>

using namespace std;

int main(int argc, char** argv) 
{
    std::string _(argv[1]);
    std::istringstream in(_);
    int64_t n, K;
    in >> n >> K;
        
    std::vector<int64_t> v_in(n), w_in(n);
    std::vector<int64_t> v(n), w(n);
    std::vector<std::tuple<double, int64_t>> v_w(n);
    for (int64_t i = 0; i < n; ++i) {
        in >> v_in[i] >> w_in[i];
        v_w[i] = std::make_tuple(-K * double(v_in[i])/double(w_in[i]), i);
    }
    std::sort(v_w.begin(), v_w.end());
    
    std::map<int64_t, int64_t> key_map;
    std::vector<double> acc_w(n, 0);
    std::vector<double> acc_v(n, 0);
    for (int64_t i = 0; i < n; ++i) {
        v[i] = v_in[std::get<1>(v_w[i])];
        w[i] = w_in[std::get<1>(v_w[i])];
        
        acc_w[i] = ((i > 0) ? acc_w[i-1]:0) + w[i];
        acc_v[i] = ((i > 0) ? acc_v[i-1]:0) + v[i];
        key_map[i] = std::get<1>(v_w[i]);
    }    
    
    int64_t n_iterations = 0;
    // (H, v, w, i)
    
    using heuristic = std::tuple<double, int64_t, int64_t, int64_t, int64_t>;
    std::set<heuristic> Q;
    
    // (w, i) --> parent
    std::map<std::tuple<int64_t, int64_t>, int64_t> parent;
    // (w, i) --> parent value
    std::map<std::tuple<int64_t, int64_t>, int64_t> value;
    int max_total = 0, max_total_w = 0, max_total_i = 0;
    int opt = 0;
    
    double max_HJ = 0;
    
    for (int64_t i = 0; i < n; ++i) {
       double Hi = v[i] + (K - w[i])*double(v[i]) / double(w[i]);
       auto key = std::make_tuple(w[i], i);
       Q.insert(std::make_tuple(-Hi, -v[i],  K-w[i],  w[i], i));
       if (v[i] > value[key] || (v[i] == value[key] && i < parent[key])) {
            parent[key] = i;
            value[key] = v[i]; 
       }
       
       if (Hi > max_HJ) {
            max_HJ = Hi;
       }
       if (v[i] > max_total) {
            max_total = v[i];
            max_total_w = w[i];
            max_total_i = i;
       }
    }
    
    
    std::default_random_engine generator;
    
    int64_t max_n_iterations = K + n + 54222111;
    int64_t max_q_size = 3000000;
    
    while (n_iterations < max_n_iterations && !Q.empty() && max_HJ >= max_total) {
        auto t = *Q.begin();
        int64_t i = std::get<4>(t);
        int64_t wi = std::get<3>(t);
        int64_t vi = -std::get<1>(t);
        int64_t hi = -std::get<0>(t);
        
        Q.erase(t);
        
        if (hi < max_total) {
            continue;
        }
       
        if (Q.size() > 1) {
            auto t = *Q.rbegin();
                    
            if (Q.size() > max_q_size || -std::get<0>(t) < max_total || -std::get<0>(t) < 0.6 * max_total) {
                bool ok = false;
                while (Q.size() > max_q_size / 2 || ok) {
                    auto t = *Q.rbegin();
                    Q.erase(t);
                    
                    if (Q.size() > 0) {
                        t = *Q.rbegin();
                        auto hj = -std::get<0>(t); 
                        ok = hj < max_total || hj < 0.6 * max_HJ;
                    }
                }
            }
        }
        
        int64_t greedy_hj = vi;
        int64_t greedy_w = wi;
        for (int64_t j = i + 1; j < n; ++j) {
            if (greedy_w + w[j] > K) {
                break;
            }
            greedy_hj += v[j];
            greedy_w += w[j];            
        }
        
        if (greedy_hj > max_total) {
            max_total = greedy_hj;
        }
                
        for (int64_t j = i + 1; j < n; ++j) {
            int64_t wj = wi + w[j];
            int64_t vj = vi + v[j];
            if (wj > K) {
                continue;
            }
            auto key = std::make_tuple(wj, j);
            double Hj = vi + (K - wi) * double(vj)/double(wj); 
            
            if (max_HJ < Hj) {
                max_HJ = Hj;
            }
            
            if (Hj < max_total || Hj < greedy_hj) {
                continue;
            }
            
            if (vi + acc_v[n-1] - acc_v[j] + v[j] < max_total) {
                continue;
            }
            
            ++n_iterations;
            
            if (n_iterations % 1000000 == 0) {
                std::cerr << "Hj: " << Hj << " max_total: " << max_total << " max_hj: " <<  max_HJ << " greedy w: " << greedy_hj << " queue size: " << Q.size() << "\n";
            }
            
            if (vj >= max_total) {
                max_total = vj;
                max_total_w = wj;
                max_total_i = j;
                //std::cerr << "max_total: " <<  vj << " max_hj: " << max_HJ << " greedy hj: " << greedy_hj << "\n";
            }
                      
            if (vj > value[key]) {
                value[key] = vj;
                parent[key] = i;
                Q.insert(std::make_tuple(-Hj,  -vj, -K+wj,  wj, j));
                //break;
            }
        }
    }
    
    std::cerr << "ok" << Q.size() << "\n";
    
    std::vector<bool> x(n, false);
    
    auto u = std::make_tuple(max_total_w, max_total_i);
    int64_t wu = max_total_w;
    while (wu > 0) {
        x[key_map[std::get<1>(u)]] = true;
        int64_t v = parent[u];
        x[key_map[v]] = true;
        wu = wu - w[std::get<1>(u)];
        u = std::make_tuple(wu, v);
    }
    
    int64_t acc = 0;
    for (int64_t i = 0; i < n; ++i) {
        acc += int64_t(x[i]) * v_in[i];
    }
    std::cerr << "is ok:" << (acc ==max_total)<< "\n";
    std::cerr << "acc:" << acc<< "\n";
    std::cerr << "max_total:" << max_total << "\n";
    
    max_total = acc;
    
    std::cerr << acc << " " << ((Q.empty() || max_HJ<=max_total)? 1: 0) << "\n"; 
    for (int64_t i = 0; i < n; ++i) {
        if (i > 0) {
            std::cerr << " ";
        }
        std::cerr << ((x[i])?1:0);
    }
    std::cerr << "\n";
    
    std::cout << max_total << " " << ((Q.empty())? 1: 0) << "\n"; 
    for (int64_t i = 0; i < n; ++i) {
        if (i > 0) {
            std::cout << " ";
        }
        std::cout << ((x[i])?1:0);
    }
    std::cout << "\n";
    std::cerr << "is ok:" << (acc ==max_total)<< "\n";
    std::cerr << "acc:" << acc<< "\n";
    std::cerr << "max_total:" << max_total << "\n";
    
    return 0;
}

