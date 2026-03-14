#pragma once
#include <cstdio>
#include <map>
#include <cmath>
#include <iostream>
#ifndef EPS
#define EPS 1e-9
#endif

namespace symbolic{
using PolyData = std::map<std::map<char, int>, double>;

class Poly{
public:
    PolyData data;
public:
    Poly(){}
    Poly(double r){ if(std::abs(r) > EPS) data ={{{}, r}}; }
    Poly(char a):data({{{{a, 1}}, 1}}){}
    ~Poly(){}
public:


public:
    template<class Archive>
    void serialize(Archive& ar, const uint version)
    { ar & data; }
public:
    const Poly& operator+=(const Poly& b){
        for(auto [t, c] : b.data){
            data[t] += c;
            if(std::abs(data[t]) < EPS) data.erase(t);
        }
        return *this;
    };
    const Poly& operator-=(const Poly& b){
        for(auto [t, c] : b.data){
            data[t] -= c;
            if(std::abs(data[t]) < EPS) data.erase(t);
        }
        return *this;
    };
    const Poly& operator*=(const Poly& b){
        Poly ret;
        auto prod = [](const std::map<char, int>& ta, const std::map<char, int>& tb){
            std::map<char, int> ret(ta);
            for(auto [eb, db] : tb){
                ret[eb] += db;
                if(ret[eb] == 0) ret.erase(eb);
            } 
            return ret;
        };
        for(auto [ta, ca] : data) for(auto [tb, cb] : b.data)
            ret.data[prod(ta, tb)] = ca * cb; 
        for(auto [t, c] : ret.data) if(std::abs(c) < EPS)
            ret.data.erase(t);
        data = ret.data;
        return *this;
    };
};
const Poly operator+(const Poly& a, const Poly& b){
    Poly ret(a);
    for(auto [t, c] : b.data){
        ret.data[t] += c;
        if(std::abs(ret.data[t]) < EPS) ret.data.erase(t);
    }
    return ret;
};
const Poly operator-(const Poly& a, const Poly& b){
    Poly ret(a);
    for(auto [t, c] : b.data){
        ret.data[t] -= c;
        if(std::abs(ret.data[t]) < EPS) ret.data.erase(t);
    }
    return ret;
};
const Poly operator*(const Poly& a, const Poly& b){
    Poly ret;
    auto prod = [](const std::map<char, int>& ta, const std::map<char, int>& tb){
        std::map<char, int> ret(ta);
        for(auto [eb, db] : tb){
            ret[eb] += db;
            if(ret[eb] == 0) ret.erase(eb);
        } 
        return ret;
    };
    for(auto [ta, ca] : a.data) for(auto [tb, cb] : b.data)
        ret.data[prod(ta, tb)] += ca * cb; 
    for(auto [t, c] : ret.data) if(std::abs(c) < EPS)
        ret.data.erase(t);
    return ret;
};
const double norm(const Poly& a){
    double ret = 0;
    for(auto [t, c] : a.data)
        ret += c * c;
    return std::sqrt(ret);
};
const double eval(const Poly& a, std::map<char, double> vs){
    double ret = 0;
    auto eval_terms = [&](const std::map<char, int>& t){
        double ret = 1;
        for(auto [e, d] : t) ret *= std::pow(vs[e], d);
        return ret;
    };
    for(auto [t, c] : a.data)
        ret += eval_terms(t) * c;
    return ret;
};

std::ostream& operator<<(std::ostream& os, const Poly& a){
    for(auto [t, c] : a.data){
        // os << "+(" << c << ")*";
        os << c << "*";
        if(t.size() == 0) os << "_1";
        else for(auto [e, d] : t)
            os << e << ":" << d << ",";  
            // os << e << "^(" << d << ")";  
    }
    return os;
}

void save_vec_coefi(std::vector<double>& os, const Poly& p){
    os.push_back(p.data.size());
    for(auto [term, coefi] : p.data){
        // save term
        os.push_back(term.size());
        for(auto [e, d] : term){
            os.push_back(e);
            os.push_back(d);
        }
        // save coefi
        os.push_back(coefi);
    }
}

void save_coefi(std::string file, const Poly& p){
    std::vector<double> data;
    save_vec_coefi(data, p);
    std::ofstream os(file);
    for(auto d : data) 
        os << std::scientific 
           << std::setprecision(std::numeric_limits<double>::max_digits10) << d << " ";
    os.close();
}


void save_vec_coefi(std::vector<double>& os, const double& p){
    os.push_back(p);
}


void read_vec_coefi(Poly& p, const std::vector<double>& is, int& idx){
    auto load = [&](){return is[idx++];};
    int data_size = load();
    p.data = {};
    for(int i = 0; i < data_size; i++){
        int term_size = load();
        std::map<char, int> term;
        for(int j = 0; j < term_size; j++){
            char e = load();
            int d = load();
            term[e] = d;
        }
        double coefi = load();
        p.data[term] = coefi;
    }
}

void read_vec_coefi(double& p, const std::vector<double>& is, int& idx){
    auto load = [&](){return is[idx++];};
    p = load(); 
}

void load_coefi(std::string file, Poly& p){
    std::vector<double> data;
    std::ifstream is(file);
    double d; int idx = 0;
    while(is>>d) data.push_back(d);
    read_vec_coefi(p, data, idx);
    is.close();
}

namespace all_polys{
    const Poly c_1 = Poly(1.); const Poly c_0 = Poly();
    const Poly c_a = Poly('a'), c_b = Poly('b'), c_c = Poly('c'), c_d = Poly('d');
    const Poly c_A = Poly('A'), c_B = Poly('B'), c_C = Poly('C'), c_D = Poly('D');
};
};