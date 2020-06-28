/*

(c) 2020, Pericles Lopes Machado  <pericles [at] newavesoft.com>

Implementation of merge sort count inversions.
 
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
#include <tuple>

std::tuple<std::vector<int>, int64_t> count_inversion(const std::vector<int>& input) {
    if (input.size() == 0) {
        return std::make_tuple(std::vector<int>(), 0);
    }
    if (input.size() == 1) {
        return std::make_tuple(input, 0);
    }
    
    std::vector<int> left(input.size()/2);
    std::vector<int> right(input.size() - left.size());
    
    for (int64_t i = 0; i < left.size(); ++i) {
        left[i]=input[i];
    }
    for (int64_t j = 0; j < right.size(); ++j) {
        right[j]=input[j+left.size()];
    }
    
    auto L = count_inversion(left);
    auto R = count_inversion(right);
    
    auto& l = std::get<0>(L);
    auto& r = std::get<0>(R);
    
    std::vector<int> D(l.size() + r.size());
    
    int64_t i = 0, j = 0;
    int64_t count = 0;
    int64_t k = 0;
    while (k < D.size()) {
        if (i < l.size() && j < r.size()) {
            if (l[i] < r[j]) {
                D[k] = l[i];
                ++i;
            } else {
                D[k] = r[j];
                ++j;
                count+=l.size()-i;
            }
        } else if (i < l.size()) {
            D[k] = l[i];
            ++i;
        } else if (j < r.size()) {
            D[k] = r[j];
            ++j;
        }
        
        ++k;
    }
    
    /*final exam
    std::cout << "count:" << count << "\n";
    for (int64_t i = 0; i < D.size(); ++i) {
        std::cout << D[i]<<" ";
    }
    std::cout << "\n";
    */
    return std::make_tuple(D, count + std::get<1>(L) + std::get<1>(R));
    
}


int main() {
    int N;
    std::cout << "# inputs:\n";
    std::cin >> N;
    std::vector<int> data(N);
    for(int i = 0; i < N; ++i) {
        std::cin >> data[i];
    }
    
    auto res = count_inversion(data);
    auto od = std::get<0>(res);
    for(int i = 0; i < N; ++i) {
        std::cout << od[i] << " ";
    }
    std::cout << "\n";
    std::cout << "# inversions:" << std::get<1>(res) << "\n";
    return 0;
}
