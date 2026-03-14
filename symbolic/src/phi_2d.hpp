#include <filesystem>
#include "terms.hpp"
#include "vec_3d.hpp"

namespace symbolic{
class phi_k_2d{
public:
    int _m, _n;
    int _k;
    TermsPoly phi_t;
public:
    phi_k_2d(){}
    phi_k_2d(int k) : _k(k), _m(0), _n(0) {initalize(k);}
    phi_k_2d(const phi_k_2d& p) : _k(p._k), _m(p._m), _n(p._n) { phi_t = p.phi_t;}
    ~phi_k_2d(){}
    void summary(){}
public:
    void initalize(int k){
        assert(std::abs(k) % 2 == 1 && "input k must be odd!");
        auto H = [](int k) {
            double ret = 1;
            for(int i=k; i > 0; i-=2)
                ret *= i/(i+1.);
            return ret;
        };
        using namespace symbolic::all_elems;
        using symbolic::all_polys::c_1;
        double c; Term t;
        if(k < 0 && (k=-k)){
            for(int i = 0; i <= k/2; i++){
                c = H(k) * H(k-2*i-1) / (2+k);
                t = s_a * pow(s_y, 1+2*i) * pow(s_ayz,k-2*i);
                phi_t.insert({t, c}); 
            }
            c = H(k) / (2+k);
            t = s_ln_0 * pow(s_y, k+2);
            phi_t.insert({t, c});
        }
        else if(k > 0 && k != 1) {
            for(int i = 1; i <= k/2; i++){
                c = H(k-3) * H(k-2*i-2) / (2-k);
                t = s_a * pow(s_y, 1-2*i) * pow(s_ayz,2*i-k);
                phi_t.insert({t, c}); 
            }
        }
        else phi_t.insert({s_y * s_ln_0, c_1});
        phi_t = simplify(phi_t);
    }

    double operator()(double x, double y, double l, std::map<char, double> vs={}) 
    { return eval(phi_t, x, y, 0, l, vs); };

    const phi_k_2d step(Poly dpx, Poly dpy, bool is_stepm){
        using all_polys::_1;
        using all_polys::_0;
        phi_t = diff(phi_t, dpx * _1, dpy * _1, _0)  * (1./_k);
        if(is_stepm) _m += 1; else _n += 1; _k += 2;
        phi_t = simplify(phi_t);
        return *this;
    };

    double operator()(Vec3d A, Vec3d B, Vec3d C, Vec3d X){
        assert( dot(cross(B-A, C-A), X-A) == 0 && "ABCX Must be coplanar!");
        if(_k > 2 && is_inside(A, B, C, X)) return INFINITY;
        Vec3d nz = normalize(cross(B-A, C-A));
        Vec3d P = A + proj(X-A, nz);
        auto cacu_edge = [&](Vec3d A_, Vec3d B_){
            auto nx = normalize(B_-A_);
            auto ny = cross(nz, nx); 
            double x = dot(P-A_, nx);
            double y = dot(P-A_, ny);
            double l = dist(B_, A_); 
            return this->operator()(x, y, l); 
        };
        return cacu_edge(A, B) + cacu_edge(B, C) + cacu_edge(C, A);
    };
    void save(std::string folder){
        std::filesystem::create_directories(folder);
        phi_t.save_vec(folder + "/phi_t.bin");
    }
    void info(){
        std::cout << "Size(" << phi_t.size() << ") = " << phi_t << std::endl;
    }
};

class phi_mnk_2d{
public: 
    int k, m, n; // target phi_mnk
    std::map<std::tuple<int, int, int>, phi_k_2d> phi_a, phi_b, phi_c;
public:
    int _k; 
    bool is_initialized = false;
public:
    phi_mnk_2d(int m_, int n_, int k_) : k(k_), m(m_), n(n_) { 
        for(int i = 0; i <= m; i++) for(int j = 0; j <= n; j++){
            const auto p = phi_k_2d(k-2*i-2*j); int pwd = k-2*i-2*j; 
            phi_a[{0, 0, pwd}] = p; phi_b[{0, 0, pwd}] = p; phi_c[{0, 0, pwd}] = p;
        }
        initialize();
    }
    phi_mnk_2d(const phi_mnk_2d& p) : _k(p._k), k(p.k), m(p.m), n(p.n), 
    phi_a(p.phi_a), phi_b(p.phi_b), phi_c(p.phi_c){}
    ~phi_mnk_2d(){}

    phi_mnk_2d initialize(){
        using namespace all_elems;
        using namespace all_polys;
        using std::get;
        using step_type = std::tuple<Vec3d_Poly, TermsPoly>;
        auto get_stepm = [&](Vec3d_Poly A, Vec3d_Poly B, Vec3d_Poly C){
            Poly vx = (get<1>(C)-get<1>(A)) * c_d;
            Poly vy = (get<0>(A)-get<0>(C)) * c_d;
            TermsPoly px = s_x - get<0>(A);
            TermsPoly py = s_y - get<1>(A);
            return step_type{{vx, vy, 0.}, vx * px + vy * py}; 
        };
        auto get_stepn = [&](Vec3d_Poly A, Vec3d_Poly B, Vec3d_Poly C){
            Poly vx = (get<1>(A)-get<1>(B)) * c_d;
            Poly vy = (get<0>(B)-get<0>(A)) * c_d;
            TermsPoly px = s_x - get<0>(A);
            TermsPoly py = s_y - get<1>(A);
            return step_type{{vx, vy, 0.}, vx * px + vy * py}; 
        };

        auto stepm_ABC = get_stepm({0.,0.,0.}, {c_a,0.,0.}, {c_b,c_c,0.});
        auto stepm_BCA = get_stepm({c_b,c_c,0.}, {0.,0.,0.}, {c_a,0.,0.});
        auto stepm_CAB = get_stepm({c_a,0.,0.}, {c_b,c_c,0.}, {0.,0.,0.});
        auto stepn_ABC = get_stepn({0.,0.,0.}, {c_a,0.,0.}, {c_b,c_c,0.});
        auto stepn_BCA = get_stepn({c_b,c_c,0.}, {0.,0.,0.}, {c_a,0.,0.});
        auto stepn_CAB = get_stepn({c_a,0.,0.}, {c_b,c_c,0.}, {0.,0.,0.});

        Vec3d_Poly stepAm = get<0>(stepm_ABC); TermsPoly pAm = get<1>(stepm_ABC);
        Vec3d_Poly stepBm = get<0>(stepm_BCA); TermsPoly pBm = get<1>(stepm_BCA);
        Vec3d_Poly stepCm = get<0>(stepm_CAB); TermsPoly pCm = get<1>(stepm_CAB);
        Vec3d_Poly stepAn = get<0>(stepn_ABC); TermsPoly pAn = get<1>(stepn_ABC);
        Vec3d_Poly stepBn = get<0>(stepn_BCA); TermsPoly pBn = get<1>(stepn_BCA);
        Vec3d_Poly stepCn = get<0>(stepn_CAB); TermsPoly pCn = get<1>(stepn_CAB);

        for(int i = 0; i < m; i++){
            for(int pwd = k-2*m-2*n+2*i+2; pwd <= k; pwd+=2){
                phi_a[{i+1, 0, pwd}] = phi_a[{i, 0, pwd-2}].step(get<0>(stepAm), get<1>(stepAm), true); 
                phi_a[{i+1, 0, pwd}].phi_t += phi_a[{i, 0, pwd}].phi_t * pAm;
                phi_b[{i+1, 0, pwd}] = phi_b[{i, 0, pwd-2}].step(get<0>(stepBm), get<1>(stepBm), true); 
                phi_b[{i+1, 0, pwd}].phi_t += phi_b[{i, 0, pwd}].phi_t * pBm;
                phi_c[{i+1, 0, pwd}] = phi_c[{i, 0, pwd-2}].step(get<0>(stepCm), get<1>(stepCm), true); 
                phi_c[{i+1, 0, pwd}].phi_t += phi_c[{i, 0, pwd}].phi_t * pCm;
            }
        }
        for(int i = 0; i < n; i++){
            for(int pwd = k-2*n+2*i+2; pwd <= k; pwd+=2){
                phi_a[{m, i+1, pwd}] = phi_a[{m, i, pwd-2}].step(get<0>(stepAn), get<1>(stepAn), false); 
                phi_a[{m, i+1, pwd}].phi_t += phi_a[{m, i, pwd}].phi_t * pAn;
                phi_b[{m, i+1, pwd}] = phi_b[{m, i, pwd-2}].step(get<0>(stepBn), get<1>(stepBn), false); 
                phi_b[{m, i+1, pwd}].phi_t += phi_b[{m, i, pwd}].phi_t * pBn;
                phi_c[{m, i+1, pwd}] = phi_c[{m, i, pwd-2}].step(get<0>(stepCn), get<1>(stepCn), false); 
                phi_c[{m, i+1, pwd}].phi_t += phi_c[{m, i, pwd}].phi_t * pCn;
            }
        }
        is_initialized = true;
        // std::cout << "========= Info Begin ==========" << std::endl;
        // phi_b[{m, n, k}].info(); // phi_b[{m, n, k}].info(); phi_c[{m, n, k}].info();
        // std::cout << "========= Info Ends ==========" << std::endl;
        return *this;
    }

    double operator()(Vec3d A, Vec3d B, Vec3d C, Vec3d X) {
        assert( std::abs(dot(X-A, cross(B-A, C-A))) < EPS && "Must be coplanar!");
        if(k > 2 && is_inside(A, B, C, X)) return INFINITY;
        Vec3d N = normalize(cross(B-A, C-A));
        double s = norm(cross(B-A, C-A));
        double a = norm(C-B), Aa = angle(B-A, C-A, N);
        double b = norm(A-C), Bb = angle(C-B, A-B, N); 
        double c = norm(B-A), Cc = angle(A-C, B-C, N);
        std::map<char, double> vs_ABC = {{'a', c}, {'b', b*cos(Aa)}, {'c', b*sin(Aa)}, {'d', 1/s}};
        std::map<char, double> vs_BCA = {{'a', a}, {'b', c*cos(Bb)}, {'c', c*sin(Bb)}, {'d', 1/s}};
        std::map<char, double> vs_CAB = {{'a', b}, {'b', a*cos(Cc)}, {'c', a*sin(Cc)}, {'d', 1/s}};

        auto xyl_edge = [&](Vec3d A_, Vec3d B_){
            auto nx = normalize(B_-A_);
            auto ny = cross(N, nx); 
            double x = dot(X-A_, nx);
            double y = dot(X-A_, ny);
            double l = dist(B_, A_); 
            return Vec3d{x, y, l};
        };

        Vec3d pA = xyl_edge(A, B);
        Vec3d pB = xyl_edge(B, C);
        Vec3d pC = xyl_edge(C, A);
        double retA = 0, retB = 0, retC = 0;
        retA = phi_a[{m, n, k}](pA[0], pA[1], pA[2], vs_ABC);
        retB = phi_b[{m, n, k}](pB[0], pB[1], pB[2], vs_BCA);
        retC = phi_c[{m, n, k}](pC[0], pC[1], pC[2], vs_CAB);
        return retA + retB + retC;
    }
    void save(std::string folder){
        using std::to_string;
        auto _m = to_string(m), _n = to_string(n), _k = to_string(k);
        auto path = folder+"/"+_m+"_"+_n+"_"+_k; 
        phi_a[{m, n, k}].save(path+"/phi_a/");
        phi_b[{m, n, k}].save(path+"/phi_b/");
        phi_c[{m, n, k}].save(path+"/phi_c/");
    }
    void info(){
        phi_a[{m, n, k}].info();
        phi_b[{m, n, k}].info();
        phi_c[{m, n, k}].info();
    }
};
}
