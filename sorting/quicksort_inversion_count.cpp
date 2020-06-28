/*

(c) 2020, Pericles Lopes Machado  <pericles [at] newavesoft.com>

Implementation of quick sort count inversions.
 
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

int64_t naive_choose_pivot(std::vector<int>& A, int64_t l, int64_t r) {
    return l;
}

int64_t naive_choose_pivot_2(std::vector<int>& A, int64_t l, int64_t r) {
    return r-1;
}

int64_t naive_choose_pivot_3(std::vector<int>& A, int64_t l, int64_t r) {

    int64_t N = r - l;
    int64_t middle = N / 2;
    if ((N % 2) == 0) {
        --middle;
    }
    int64_t min_v = (A[l] < A[l+N-1])?A[l]:A[l+N-1];
    int64_t max_v = (A[l] > A[l+N-1])?A[l]:A[l+N-1];
    int64_t median_v = A[l+middle];
    if (median_v < min_v) {
        std::swap(median_v, min_v);
    }
    if (median_v > max_v) {
        std::swap(median_v, max_v);
    }
    if (A[l] == median_v) {
        return l;
    } 
    if (A[l + N-1] == median_v) {
        return r-1;
    }
    return middle + l;
}

int64_t partition(std::vector<int>& A, int64_t l, int64_t r) {
    int64_t i = l + 1;
    for (int64_t j = l + 1; j < r; ++j) {
        if (A[j] < A[l]) {
            std::swap(A[i], A[j]);
            ++i;
        }
    }
    std::swap(A[i-1],A[l]);
    return i-1;
}

int64_t quick_sort(std::vector<int>& A, int64_t l, int64_t r) {
    if (l >= r) {
        return 0;
    }
    int64_t p = naive_choose_pivot_3(A, l, r);
    std::swap(A[p], A[l]);
    int64_t k = partition(A, l, r);
    int64_t cl = quick_sort(A, l, k);
    int64_t cr = quick_sort(A, k + 1, r);
    return cr + cl + r - l - 1;
}

int64_t quick_sort(std::vector<int>& A) {
    return quick_sort(A, 0, A.size());
}


int main() {
    int64_t N;
    std::cin >> N;
    std::vector<int> A(N);
    for (int64_t i = 0; i < N; ++i) {
        std::cin >> A[i];
    }
    int count = quick_sort(A);
    
    for (int64_t i = 0; i < N; ++i) {
        std::cout << A[i] << "  ";
    }
    std::cout << "\n";
    
    std::cout << "comparisons count: " << count << "\n";
    return 0;
}


