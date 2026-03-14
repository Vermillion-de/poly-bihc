#pragma once
#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cassert>
#include <functional>
#include <map>
#include <iomanip>
#include <tuple>
#include "poly.hpp"
#include "utils.hpp"

#ifndef EPS
#define EPS 1e-9
#endif

namespace symbolic{
enum Elem{ a__, x__, _y_, __z, ay_, _yz, ayz, _as, _ac, _ln};
typedef std::map<Elem, int> Term;

//// Functions
// From symb to Terms;
const Term symb() {return Term{};};
const Term symb(const Elem& e) {return Term{{e, 1}};}; 
const Term pow(const Elem& e, int n) {return Term{{e, n}};}; 
const Term pow(const Term& t, int n) {
    if(n == 0) return symb();
    Term ret;
    for(auto [e, d] : t) ret[e] = d * n;
    return ret;
};
const Term simplify_t(const Term& t){
    Term ret;
    for(auto [e, d] : t)
        if(d == 0 && e != _ln) continue;
        else ret[e] = d;
    return ret;   
};

const Term operator+(const Term& t){ return t; }
const Term operator*(const Term& ta, const Term& tb){
    Term ret(ta); 
    for(auto [e, d] : tb) 
        ret[e] += d;
    return simplify_t(ret);
};
const Term operator/(const Term& ta, const Term& tb){
    Term ret(ta); 
    for(auto [e, d] : tb)
        ret[e] -= d;
    return simplify_t(ret);
};
const Term& operator*=(Term& a, const Term& b){
    for(auto [e, d] : b){
        a[e] += d;
        if(a[e] == 0 && e != _ln) a.erase(e);
    }
    return a;
};
const Term operator/=(Term& a, const Term& b){
    for(auto [e, d] : b){
        a[e] -= d;
        if(a[e] == 0 && e != _ln) a.erase(e);
    }
    return a;
};

namespace all_elems{
    Elem elems[] = { a__, x__, _y_, __z, ay_, _yz, ayz, _as, _ac, _ln };
    std::string names[] = { "a", "x", "y", "z", "ay", "yz", "ayz", "asin", "acos", "ln" };
    Term s_a = symb(a__), s_x = symb(x__), s_y = symb(_y_), s_z = symb(__z);
    Term s_ay= symb(ay_), s_yz= symb(_yz),s_ayz= symb(ayz);
    Term s_as= symb(_as), s_ac= symb(_ac), s_ln= symb(_ln), s_ln_0 = {{_ln, 0}};
    bool is_sqrt[] = {false, false, false, false, true, true, true, false, false, false};
    bool has_a[]   = {true, false, false, false, true, false, true, true, true, true};
};

std::ostream& operator<<(std::ostream& os, const Term& t)
{ os << "["; for(auto [e, d] : t) os << all_elems::names[e] << ":" << d << ","; os << "]"; return os; };


template<typename T>
class _Terms{
public:
    std::map<Term, T> data;
public:
    _Terms(){}
    _Terms(const T& d):data({{{},d}}){};
    _Terms(const Term& t):data({{t, T(1.)}}){};
    _Terms(const Term& t, const T& d):data({{t, d}}){};
    ~_Terms(){}
public:

    const _Terms apply_term(std::function<const Term&(const Term&)> func){
        _Terms ret;
        for(auto [t, c] : data) 
            ret.data[simplify_t(func(t))] = c;
        return ret;
    }
    const _Terms simplify(){
        _Terms ret;
        for(auto [t, c] : data) if(norm(c) > EPS)
            ret.data[simplify_t(t)] = c;
        return ret;
    }
    void insert(const std::tuple<Term, T>& d)
    { data.insert({std::get<0>(d), std::get<1>(d)}); }
    uint size(){return data.size();}
    void clear(){data.clear();}

public:
    // "+=, -=, *=, /=", [_Terms] *** [_Terms, Term, T]
    const _Terms& operator+=(const _Terms& b){
        for(auto [t, c] : b.data){
            data[t] += c;
            if(norm(data[t]) < EPS) data.erase(t);
        }
        return *this;
    };
    const _Terms& operator-=(const _Terms& b){
        for(auto [t, c] : b.data){
            data[t] -= c;
            if(norm(data[t]) < EPS) data.erase(t);
        }
        return *this;
    };

    const _Terms& operator*=(const Term& tb){
        *this = apply_term([&](const Term& t){return t*tb;});
    };
    const _Terms& operator*=(const T& c){
        if(norm(c) < EPS) 
            data = {};
        else for(auto [t, _c] : data)
            data[t] = c * _c;
        return *this;
    };
    const _Terms& operator/=(const Term& tb){
        *this = apply_term([&](const Term& t){return t/tb;});
    };


friend const _Terms operator+(const _Terms& a, const _Terms& b){
    _Terms ret(a);
    for(auto [t, c] : b.data)
        ret.data[t] += c;
    return ret.simplify();
};
friend const _Terms operator-(const _Terms& a, const _Terms& b){
    _Terms ret(a);
    for(auto [t, c] : b.data)
        ret.data[t] -= c;
    return ret.simplify();
}
friend const _Terms operator*(const _Terms& a, const _Terms& b){
    _Terms ret; 
    for(auto [ta, ca] : a.data) for(auto [tb, cb] : b.data)
        ret.data[simplify_t(ta * tb)] += ca * cb;
    return ret.simplify();
}
friend const _Terms operator/(const _Terms& a, const Term& b){
    _Terms ret; 
    for(auto [t, c] : a.data) 
        ret.data[simplify_t(t/b)] = c;
    return ret;
};

friend std::ostream& operator<<(std::ostream& os, const _Terms& p){
    for(auto [t, coefi] : p.data)
        os << " + (" << coefi << ")*" << t; 
    return os;
};

public:

//"+": [_Terms] + [Term, T], [Term, T] + [Term, T]
friend const _Terms operator+(const _Terms& a, const Term& b) {return a+_Terms(b);};
friend const _Terms operator+(const Term& b, const _Terms& a) {return a+_Terms(b);};
friend const _Terms operator+(const _Terms& a, const T& b) {return a+_Terms(b);};
friend const _Terms operator+(const T& b, const _Terms& a) {return a+_Terms(b);};

//"-": [_Terms] + [Term, T], [Term, T] + [Term, T]  
friend const _Terms operator-(const _Terms& a, const Term& b) {return a-_Terms(b);};
friend const _Terms operator-(const Term& a, const _Terms& b) {return _Terms(a)-b;};
friend const _Terms operator-(const _Terms& a, const T& b) {return a-_Terms(b);};
friend const _Terms operator-(const T& a, const _Terms& b) {return _Terms(a)-b;};

//"*": [_Terms] + [Term, T], [Term, T] + [Term, T]   
friend const _Terms operator*(const _Terms& a, const Term& b) {return a*_Terms(b);};
friend const _Terms operator*(const Term& a, const _Terms& b) {return b*_Terms(a);};
friend const _Terms operator*(const _Terms& a, const T& b) {return a*_Terms(b);};
friend const _Terms operator*(const T& b, const _Terms& a) {return a*_Terms(b);};

public:
const _Terms diff(const _Terms& dpx, const _Terms& dpy, const _Terms& dpz){
    using namespace all_elems;
    auto diff_elem = [&](const Elem& e, int _d){
        _Terms t1, t2, t3; T d(_d*1.);
        _Terms ret;
        switch (e) {
        case a__: ret = d * dpx / s_a; break;
        case x__: ret = d * dpx / s_x; break;
        case _y_: ret = d * dpy / s_y; break;
        case __z: ret = d * dpz / s_z; break;
        case ay_: 
            ret = d * (dpx * s_a + dpy * s_y) / (s_ay * s_ay); break;
        case _yz: 
            ret = d * (dpy * s_y + dpz * s_z) / (s_yz * s_yz); break;
        case ayz:
            ret = d * (dpx * s_a + dpy * s_y + dpz * s_z) / (s_ayz * s_ayz); break;
        case _as:
            t1 = (s_z * dpx + s_a * dpz) / (s_y * s_ayz);
            t2 = -(s_y * dpy + s_z * dpz) * s_z * s_a / (s_y * s_ayz * s_yz * s_yz);
            t3 = -(s_a * dpx + s_y * dpy) * s_z * s_a / (s_y * s_ayz * s_ay * s_ay); 
            ret = d * (t1 + t2 + t3) / s_as; 
            break;
        case _ac:
            t1 = -(s_y * dpx - s_a * dpy) / pow(s_ay, 2);
            ret = d * t1 / s_ac;
            break;
        case _ln:
            if(_d == 0) d = T(1.);
            t1 = (s_a * dpx + s_y * dpy + s_z * dpz) / s_ayz;
            t2 = dpx;
            ret = d * (t1 + t2) / s_ln;
            break;
        }
        return ret.simplify();
    };

    auto diff_term = [&](const Term& t, const T& coefi){
        _Terms ret;
        for(auto [e, d] : t)
            ret += diff_elem(e, d) * coefi;
        return ret * t;
    };

    auto diff_terms = [&](const _Terms& p){
        _Terms ret;
        for(auto [t, coefi] : p.data){
            ret += diff_term(t, coefi);
        }
        return ret;
    };

    return diff_terms(*this);
}

public:
    void save(std::vector<double>& os){
        os.push_back(data.size());
        for(auto [term, coefi] : data){
            os.push_back(term.size());
            for(auto [e, d] : term){
                os.push_back(e);
                os.push_back(d);
            }
            save_vec_coefi(os, coefi);
        }
    }
    
    void save_vec(std::string filename){
        std::vector<double> ds;
        save(ds);
        std::ofstream os(filename);
        for(auto d : ds) 
            os << std::scientific 
               << std::setprecision(std::numeric_limits<double>::max_digits10) << d << " ";
        os.close();
        std::cout << "[INFO] save_vec into " << filename << ", size = " << ds.size() << ", map = " << data.size() << std::endl; 
    }

    void load_vec(std::string filename){
        std::cout << "[INFO] load_vec from " << filename << std::endl; 
        std::ifstream is(filename);
        double d; std::vector<double> ds;
        while (is >> d) ds.push_back(d);
        int idx = 0;
        read(ds, idx);
    }

    void read(std::vector<double>& is, int& idx){
        auto load = [&]() -> double {return is[idx++];};
        int data_size = load();
        data = {};
        for(int i = 0; i < data_size; i++){
            int term_size = load();
            Term term;
            for(int j = 0; j < term_size; j++){
                Elem e = all_elems::elems[int(load())];
                int d = load();
                term[e] = d;
            }
            T coefi;
            read_vec_coefi(coefi, is, idx);
            data[term] = coefi;
        }
    }

    void save(std::string filename){
        std::cout << "[save]: save to " << filename << " of size " << data.size() << std::endl;
        utils::save_binary(data, filename);
    }
    void load(std::string filename){
        data = utils::load_binary<std::map<Term, T>>(filename);
    }
    // 必须提供 begin() 和 end()
    auto begin() -> decltype(data.begin()) { return data.begin(); }
    auto end()   -> decltype(data.end())   { return data.end(); }
    // 如果是 const 对象，还需要 const 版本
    auto begin() const -> decltype(data.begin()) { return data.begin(); }
    auto end()   const -> decltype(data.end())   { return data.end(); }
};

template <typename T> const _Terms<T> operator+(const Term& b, const Term& a) {return _Terms<T>(a)+_Terms<T>(b);}; 
template <typename T> const _Terms<T> operator+(const Term& a, const T& b) {return _Terms<T>(a)+_Terms<T>(b);}; 
template <typename T> const _Terms<T> operator+(const T& b, const Term& a) {return _Terms<T>(a)+_Terms<T>(b);}; 
template <typename T> const _Terms<T> operator-(const Term& a, const Term& b) {return _Terms<T>(a)-_Terms<T>(b);}; 
template <typename T> const _Terms<T> operator-(const Term& a, const T& b) {return _Terms<T>(a)-_Terms<T>(b);}; 
template <typename T> const _Terms<T> operator-(const T& a, const Term& b) {return _Terms<T>(a)-_Terms<T>(b);}; 
template <typename T> const _Terms<T> operator*(const T& b, const Term& a) {return _Terms<T>(a, b);};
template <typename T> const _Terms<T> operator*(const Term& a, const T& b) {return _Terms<T>(a, b);};
template <typename T> const _Terms<T> operator/(const T& a, const Term& b) {return _Terms<T>(pow(b, -1), a);}; 
template <typename T> const _Terms<T> operator+(const _Terms<T>& p){return p;};
template <typename T> const _Terms<T> operator-(const _Terms<T>& p){return p*T(-1.);};
template <typename T> const _Terms<T> diff(_Terms<T>& p, const _Terms<T> dpx, const _Terms<T> dpy, const _Terms<T> dpz) { return p.diff(dpx, dpy, dpz); }
template <typename T> const _Terms<T> simplify(_Terms<T>& p) { return p.simplify(); }

// Terms of double as default
using Terms=_Terms<double>;
using Vec3d_Terms=std::tuple<Terms, Terms, Terms>;
double eval(const Terms& p, double x, double y, double z, double l){
    const double a = x, b = x - l;
    const double va__ = a, vb__ = b; // for alpha,
    const double vay_ = std::sqrt(a*a+y*y); // [alpha, y],
    const double vby_ = std::sqrt(b*b+y*y);
    const double vayz = std::sqrt(a*a+y*y+z*z); // [alpha, y, z]
    const double vbyz = std::sqrt(b*b+y*y+z*z);
    const double v_as_a = std::asin(a*z/std::sqrt((y*y+z*z)*(a*a+y*y))); // _as
    const double v_as_b = std::asin(b*z/std::sqrt((y*y+z*z)*(b*b+y*y)));
    const double v_ac_a = std::acos(a/std::sqrt(a*a+y*y)); // _ac
    const double v_ac_b = std::acos(b/std::sqrt(b*b+y*y));
    const double v_ln_a = std::sqrt(a*a+y*y+z*z)+a, v_log_a=std::log(v_ln_a); // _ln
    const double v_ln_b = std::sqrt(b*b+y*y+z*z)+b, v_log_b=std::log(v_ln_b);
    const double vx__ = x, v_y_ = y, v__z = z; // for x, y, z, yz
    const double v_yz = std::sqrt(y*y + z*z);

    auto eval_elem = [&](const Elem& e, int d, bool is_a){
        switch (e)
        {
        case a__: return std::pow(is_a? va__ : vb__, d); 
        case x__: return std::pow(vx__, d);
        case _y_: return std::pow(v_y_, d);
        case __z: return std::pow(v__z, d);
        case ay_: return std::pow(is_a? vay_ : vby_, d);
        case _yz: return std::pow(v_yz, d);
        case ayz: return std::pow(is_a? vayz : vbyz, d);
        case _as: return std::pow(is_a? v_as_a : v_as_b, d);
        case _ac: return std::pow(is_a? v_ac_a : v_ac_b, d);
        case _ln:
            if(d == 0) return (is_a? v_log_a : v_log_b);
            else return std::pow(is_a? v_ln_a : v_ln_b, d);
        default: return 1.;
        }
    };

    auto eval_term = [&](const Term& p, double coefi) {
        double ret = 1, temp_a=1, temp_b=1;
        for(auto [e, d] : p){
            if(!all_elems::has_a[e])
                ret *= eval_elem(e, d, false);
            else{
                temp_a *= eval_elem(e, d, true);
                temp_b *= eval_elem(e, d, false);
            }
        }
        if(std::abs(temp_a - temp_b) > EPS) 
           ret *= (temp_a - temp_b); 
        if(std::isnan(ret) || std::isinf(ret)){
            // std::cout << "Warning: bad value from " << p << " * " << coefi << " --> " << ret << std::endl;
            ret = 0;
        }
        return ret * coefi;
    };

    double ret = 0;
    for(auto [t, coefi] : p.data)
        ret += eval_term(t, coefi); 
    return ret;
};

const double norm(const double& r)
{ return std::abs(r); }

namespace all_terms{
    const Terms _1 = Terms({}, 1.);
    const Terms _0 = Terms();
}

// Terms of Poly
using TermsPoly = _Terms<Poly>;
using Vec3d_Poly = std::tuple<Poly, Poly, Poly>;
using Vec3d_TermsPoly = std::tuple<TermsPoly, TermsPoly, TermsPoly>;

namespace all_polys{
    const TermsPoly _1 = TermsPoly({}, 1.); 
    const TermsPoly _0 = TermsPoly();
}

Terms eval(const TermsPoly& p, std::map<char, double> vs){
    Terms ret;
    for(auto [t, c] : p.data){
        double _c = eval(c, vs);
        if(norm(_c) > EPS) ret.data[t] = _c;
    }
    return ret;
}

double eval(const TermsPoly& p, double x, double y, double z, double l, std::map<char, double> vs={}){
    // if(std::abs(x) < 1e-5 && std::abs(y) < 1e-5){ y += 1e-5, x-=1e-5; }
    const double a = x, b = x - l;
    const double va__ = a, vb__ = b; // for alpha,
    const double vay_ = std::sqrt(a*a+y*y); // [alpha, y],
    const double vby_ = std::sqrt(b*b+y*y);
    const double vayz = std::sqrt(a*a+y*y+z*z); // [alpha, y, z]
    const double vbyz = std::sqrt(b*b+y*y+z*z);
    const double v_as_a = std::asin(a*z/std::sqrt((y*y+z*z)*(a*a+y*y))); // _as
    const double v_as_b = std::asin(b*z/std::sqrt((y*y+z*z)*(b*b+y*y)));
    const double v_ac_a = std::acos(a/std::sqrt(a*a+y*y)); // _ac
    const double v_ac_b = std::acos(b/std::sqrt(b*b+y*y));
    const double v_ln_a = std::sqrt(a*a+y*y+z*z)+a, v_log_a=std::log(v_ln_a); // _ln
    const double v_ln_b = std::sqrt(b*b+y*y+z*z)+b, v_log_b=std::log(v_ln_b);
    const double vx__ = x, v_y_ = y, v__z = z; // for x, y, z, yz
    const double v_yz = std::sqrt(y*y + z*z);

    auto eval_elem = [&](const Elem& e, int d, bool is_a){
        switch (e)
        {
        case a__: return std::pow(is_a? va__ : vb__, d); 
        case x__: return std::pow(vx__, d);
        case _y_: return std::pow(v_y_, d);
        case __z: return std::pow(v__z, d);
        case ay_: return std::pow(is_a? vay_ : vby_, d);
        case _yz: return std::pow(v_yz, d);
        case ayz: return std::pow(is_a? vayz : vbyz, d);
        case _as: return std::pow(is_a? v_as_a : v_as_b, d);
        case _ac: return std::pow(is_a? v_ac_a : v_ac_b, d);
        case _ln:
            if(d == 0) return (is_a? v_log_a : v_log_b);
            else return std::pow(is_a? v_ln_a : v_ln_b, d);
        default: return 1.;
        }
    };

    auto eval_term = [&](const Term& p, double coefi) {
        double ret = 1, temp_a=1, temp_b=1;
        for(auto [e, d] : p){
            if(!all_elems::has_a[e])
                ret *= eval_elem(e, d, false);
            else{
                temp_a *= eval_elem(e, d, true);
                temp_b *= eval_elem(e, d, false);
            }
        }
        if(std::abs(temp_a - temp_b) > EPS) 
           ret *= (temp_a - temp_b); 
        if(std::isnan(ret) || std::isinf(ret) || std::abs(ret) > 1e20){
            // std::cout << "Nan from: " << p << " * " << coefi << " --> " << ret << std::endl;
            ret = 0;
        }
        return ret;
    };

    double ret = 0;
    for(auto [t, p] : p.data){
        double coefi = eval(p, vs);
        if(std::abs(coefi) > EPS)
            ret += eval_term(t, coefi) * coefi; 
    }
    return ret;
}

};