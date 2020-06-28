/*

(c) 2020, Pericles Lopes Machado  <pericles [at] newavesoft.com>

Implementation of Karatsuba fast multiplication.
 
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

#include <cmath>
#include <iostream>
#include <vector>

class Number {
private:
    std::vector<int64_t> n;
    bool reverse = false;
public:
    static const int64_t max_size;
    
    Number(const Number& other) : n(other.n), reverse(other.reverse) {
    
    }
    Number(const std::string& x) {
        if (x.size()%9) {
            n.resize(x.size() / 9+1);
        } else {
            n.resize(x.size() / 9);
        }
        for (int64_t i = 0; i < n.size(); ++i) {
            n[i] = 0;
        }
        for (int64_t i = 0; i <x.size() ; ++i) {
            int64_t p = (x.size()-i-1) / 9;
            n[n.size()-p-1] = static_cast<int64_t>(x[i] -'0') + 10 * n[n.size()-p-1];
        }
    }
    
    Number() {
    }
    
    Number(int64_t v) {
        n.push_back(v);
    }
    
    int64_t size() const {
        return n.size();
    }
    int64_t& operator[](int64_t i)  {
        if (reverse) {
            return n[n.size() - 1 - i];
        }
        return n[i];
    }
    int64_t operator[](int64_t i) const {
        if (reverse) {
            return n[n.size() - 1 - i];
        }
        return n[i];
    }
    
    Number operator+(const Number& other) const {
        Number sum;
        sum.reverse = true;
        int64_t k = other.size() - 1, j = this->size() - 1;
        int64_t remind = 0;
        while (true) {
            if (k > -1 && j > -1) {
                int64_t s = other[k] + (*this)[j] + remind;
                sum.n.push_back(s % max_size);
                remind = s / max_size;
            } else if (k > -1) {
                int64_t s = other[k] + remind;
                sum.n.push_back(s % max_size);
                remind = s / max_size;
            } else if (j > -1) {
                int64_t s = (*this)[j] + remind;
                sum.n.push_back(s % max_size);
                remind = s / max_size;
            } else {
                break;
            }
            --k; --j;
        }
        if (remind > 0){
            sum.n.push_back(remind);
        }
        
        return sum;
    }
    bool operator<(const Number& other) const {
        if (this->size() < other.size()) {
            return true;
        }
        if (this->size() > other.size()) {
            return false;
        }
        for (int64_t i = 0; i < this->size(); ++i) {
            if ((*this)[i] >= other[i]) {
                return false;
            }
        }
        return true;
    }
    bool operator==(const Number& other) const {
        if (this->size() != other.size()) {
            return false;
        }
        for (int64_t i = 0; i < this->size(); ++i) {
            if ((*this)[i] != other[i]) {
                return false;
            }
        }
        return true;
    }
    Number operator-(const Number& other) const {        
        if (*this < other) {
            auto v = other - *this;
            v[0] = -v[0];
            return v;
        }
        
        // this > other
        
        Number sub;
        
        if (*this == other) {
            sub.n.push_back(0);
            return sub;
        }
        
        sub.reverse = true;
        
        int64_t j = other.size() - 1, k = this->size() - 1;
        bool offset=false;
        
        while (true) {
            if (k > -1 && j > -1) {
                auto s =  (*this)[k] - other[j] - ((offset)?1:0); 
                if (s < 0) {
                    offset = true;
                    sub.n.push_back(max_size + s);
                } else {
                    offset = false;
                    sub.n.push_back(s);
                }
            } else if (k > -1) {
                auto s = (*this)[k] - ((offset)?1:0);
                if (s < 0) {
                    offset = true;
                    sub.n.push_back(max_size + s);
                } else {
                    offset = false;
                    sub.n.push_back(s);
                }
            }   else {
                break;
            }
            --k; --j;
        }
        
        return sub;
    }
    
    Number operator*(const Number& other) const {
        if (this->size() == 1) {
            Number res;
            res.reverse = true;
            int64_t v = (*this)[0];
            int64_t k = other.size() - 1;
            int64_t remind = 0;
            while (k>-1) {
                int64_t s = v * other[k] + remind;
                remind = s / max_size;
                res.n.push_back(s %max_size);
                --k;
            }
            if (remind > 0){
                res.n.push_back(remind);
            }
            return res;
        } else if (other.size() == 1) {
            return other * (*this);
        }
        
        
        int64_t m = this->size() / 2;
        int64_t p = this->size() - m;
        
        Number x1, x0;
        Number y1, y0;
        
        for (int64_t i = 0; i < p; ++i) {
            x1.n.push_back((*this)[i]);
        }
        for (int64_t i = p; i < this->size(); ++i) {
            x0.n.push_back((*this)[i]);
        }
        
        if (other.size() >= m) {
            int64_t  q = other.size() - m;
            for (int64_t i = 0; i < q; ++i) {
                y1.n.push_back(other[i]);
            } 
            if (q == 0) {
                y1.n.push_back(0);
            }
            for (int64_t i = q; i < other.size(); ++i){
                y0.n.push_back(other[i]);
            } 
        } else {
            y1.n.push_back(0);
            for (int64_t i = 0; i < other.size(); ++i){
                y0.n.push_back(other[i]);
            }
        }
        
        auto z2 = x1 * y1;
        auto z0 = x0 * y0;
        auto z1 = (x0+x1)*(y1+y0)-z2-z0;
        
        return z2.mul_pow10(2*m) + z1.mul_pow10(m) + z0;
    }
    
    // this * 10^m
    Number mul_pow10(int64_t m) const {
        Number res;
        for (int64_t i = 0; i < this->size(); ++i) {
            res.n.push_back((*this)[i]);
        }
        for (int64_t i = 0; i < m; ++i) {
            res.n.push_back(0);
        }
        return res;
    }
    
    
};

const int64_t Number::max_size = 1000000000;


std::ostream& operator<<(std::ostream& out, const Number& n) {
    for (int64_t i = 0; i < n.size() ; ++i) {
        if (i > 0) {
            int64_t nz = 8 - (n[i] <= 1?(0): static_cast<int64_t>(std::log10(n[i])));
            for (int64_t k = 0; k < nz; ++k) {
                out << "0";
            }
        } else if (n[i] == 0 && n.size() > 1) {
            continue;
        }
        out << n[i];
    }
    return out;
}

int main()
{
    std::string a, b;
    std::cout << "a>>  ";
    std::cin >> a;
    std::cout << "b>>  ";
    std::cin >> b;
    
    Number x(a), y(b);
    
    std::cout << "a:   " << x << "\n";
    std::cout << "b:   " << y << "\n";
    
    std::cout << "a+b: " << (x+y) << "\n";
    std::cout << "a-b: " << (x-y) << "\n";
    std::cout << "a*b: " << (x*y) << "\n";
    
	return 0;
}
