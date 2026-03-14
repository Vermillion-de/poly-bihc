#include <filesystem>
#include "terms.hpp"
#include "vec_3d.hpp"

namespace symbolic{
class phi_k
{
public:
    int _k;
    TermsPoly phi_t; // regular term
    TermsPoly phi_l; // ln term
    TermsPoly phi_c; // cos term
    TermsPoly phi_s; // sin term
    TermsPoly phi_z; // z^k / (k+2) term
public:
    phi_k(){}
    phi_k(int k) : _k(k){initalize(k);}
    phi_k(const phi_k& p) : _k(p._k)
    { phi_t = p.phi_t, phi_l = p.phi_l, phi_c = p.phi_c, phi_s = p.phi_s, phi_z = p.phi_z;}
    ~phi_k(){}
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
        double c; Term t; _k = k;
        if(k < 0 && (k=-k)){
            for(int i = 0; i <= k/2; i++) for(int j = 0; j <= k/2 - i; j++){
                c = H(k-2*i)*H(k-2*i-2*j-1)/(k+2); 
                t = s_a * s_y * pow(s_z,2*i) * pow(s_yz,2*j) * pow(s_ayz,k-2*i-2*j);
                phi_t.insert({t, c}); 
            }
            phi_t = simplify(phi_t);
            for(int i = 0; i <= k/2 + 1; i++){
                c = H(k-2*i) / (k+2);
                t = s_y * s_ln_0 * pow(s_z,2*i) * pow(s_yz,k+1-2*i);
                phi_l.insert({t, c}); 
            }
            phi_l = simplify(phi_l);
            phi_s.insert({s_as, c_1}); 
            phi_c.insert({s_ac, c_1});
            phi_z.insert({pow(s_z, k+2), 1./(k+2)});
        }
        else if(k > 0 && k!=1) {
            for(int i = 1; i <= k/2; i++) for(int j = 1; j <= k/2 - i; j++){
                c = -H(k-2*i-3)*H(k-2*i-2*j-2)/(2-k);
                t = s_a * s_y * pow(s_z,-2*i) * pow(s_yz,-2*j) * pow(s_ayz,2*i+2*j-k);
                phi_t.insert({t, c}); 
            }
            phi_t = simplify(phi_t);
            phi_s.insert({s_as, c_1}); 
            phi_c.insert({s_ac, c_1});
            phi_z.insert({pow(s_z, 2-k), 1./(2-k)});
        }
        else {
            phi_t.insert({s_y * s_ln_0, c_1});
            phi_s.insert({s_as, c_1}); 
            phi_c.insert({s_ac, c_1});
            phi_z.insert({s_z , c_1});
        }
    }
public:
    double operator()(double x, double y, double z, double l, std::map<char, double> vs={}) {
        auto sign = [](double x){return x>0? 1: (x<0? -1:0);};
        double v_tm = eval(phi_t, x, y, z, l, vs); 
        double v_ln = eval(phi_l, x, y, z, l, vs); 
        double v_as = eval(phi_s, x, y, z, l, vs); 
        double v_ac = eval(phi_c, x, y, z, l, vs); 
        double v_z  = eval(phi_z, x, y, z, l, vs); 
        double t1 = v_tm + v_ln; 
        double t2 = v_as * v_z * sign(y);
        double t3 = v_ac * v_z * sign(y*z); 
        std::cout << "t1 = "<< v_tm << " + " << v_ln << std::endl;
        std::cout << "t2 = "<< v_ac << " * " << v_z << "*" << sign(y) << std::endl;
        std::cout << "t3 = "<< v_as << " * " << v_z << "*" << sign(y*z) << std::endl;
        return t1 + t2 + t3;
    };

    double operator()(Vec3d A, Vec3d B, Vec3d C, Vec3d X){
        Vec3d nz = normalize(cross(B-A, C-A));
        Vec3d P = A + proj(X-A, nz);
        double z = dot(X-A, nz);
        auto cacu_edge = [&](Vec3d A_, Vec3d B_){
            auto nx = normalize(B_-A_);
            auto ny = cross(nz, nx); 
            double x = dot(P-A_, nx);
            double y = dot(P-A_, ny);
            double l = dist(B_, A_); 
            return this->operator()(x, y, z, l); 
        };
        return cacu_edge(A, B) + cacu_edge(B, C) + cacu_edge(C, A);
    };

    void step(const TermsPoly& dpx, const TermsPoly& dpy, const TermsPoly& dpz){
        phi_t = diff(phi_t, dpx, dpy, dpz);
        phi_l = diff(phi_l, dpx, dpy, dpz);
        phi_t += diff(phi_c, dpx, dpy, dpz) * phi_z;
        phi_t += diff(phi_s, dpx, dpy, dpz) * phi_z;
        phi_z = diff(phi_z, dpx, dpy, dpz);
    }

    void info(){
        std::cout << "==================" << std::endl;
        std::cout << phi_t.size() << " + " << phi_l.size() << " + " << phi_z.size() << " + " << std::endl;
        std::cout << "[" << phi_t << "] + [" << phi_l << "]" << " + " << std::endl;
        std::cout << "[" << phi_s * phi_z << "] + [" << phi_c * phi_z << "]" << " + " << std::endl;
    }
    void save(std::string folder){
        std::filesystem::create_directories(folder);
        phi_t.save_vec(folder + "/phi_t.bin");
        phi_l.save_vec(folder + "/phi_l.bin");
        // phi_c.save_vec(folder + "/phi_c.bin");
        // phi_s.save_vec(folder + "/phi_s.bin");
        phi_z.save_vec(folder + "/phi_z.bin");
    }

    void load(std::string folder){
        using namespace symbolic::all_elems;
        using symbolic::all_polys::c_1;
        phi_s.clear(); phi_s.insert({s_as, c_1}); 
        phi_c.clear(); phi_c.insert({s_ac, c_1});
        phi_t.load_vec(folder + "/phi_t.bin");
        phi_l.load_vec(folder + "/phi_l.bin");
        // phi_c.load_vec(folder + "/phi_c.bin");
        // phi_s.load_vec(folder + "/phi_s.bin");
        phi_z.load_vec(folder + "/phi_z.bin");
    }
};

class phi_mnk{ // this is specially designed for a specific triangle
public: 
    int k, m, n; // target phi_mnk
    phi_k phi_a, phi_b, phi_c;
    Poly r;
public:
    int _k;
    bool is_initialized = false;
public:
    phi_mnk(int m_, int n_, int k_) : _k(k_-2*m_-2*n_), k(k_), m(m_), n(n_),
    phi_a(k_-2*m_-2*n_), phi_b(k_-2*m_-2*n_), phi_c(k_-2*m_-2*n_), r(1.) { initialize(); }

    phi_mnk(int m_, int n_, int k_, std::string file) : _k(k_-2*m_-2*n_), k(k_), m(m_), n(n_),
    r(1.) { load(file); }

    phi_mnk(const phi_mnk& p) : _k(p._k), k(p.k), m(p.m), n(p.n),
    phi_a(p.phi_a), phi_b(p.phi_b), phi_c(p.phi_c) {}

    ~phi_mnk(){}

    phi_mnk initialize(){
        using namespace symbolic::all_polys;
        using namespace symbolic::all_elems;
        using std::get;
        // remark: c_a, c_b, c_c, specified, c_d = 1 / norm(cross(AB, AC)) 
        auto get_stepm = [&](Vec3d_Poly A, Vec3d_Poly B, Vec3d_Poly C){
            Vec3d_TermsPoly AX = {s_x - get<0>(A), s_y - get<1>(A), s_z - get<2>(A)};
            Poly ACx = get<0>(C)-get<0>(A), ACy = get<1>(C)-get<1>(A), ACz = get<2>(C)-get<2>(A);  
            Vec3d_TermsPoly AC = {ACx/s_z, ACy/s_z, ACz/s_z};
            Vec3d_TermsPoly _ret = cross(AC, AX); 
            return _ret;
        };
        auto get_stepn = [&](Vec3d_Poly A, Vec3d_Poly B, Vec3d_Poly C){
            Vec3d_TermsPoly AX = {s_x - get<0>(A), s_y - get<1>(A), s_z - get<2>(A)};
            Poly ABx = get<0>(B)-get<0>(A), ABy = get<1>(B)-get<1>(A), ABz = get<2>(B)-get<2>(A);  
            Vec3d_TermsPoly AB = {ABx/s_z, ABy/s_z, ABz/s_z};
            Vec3d_TermsPoly ret = cross(AX, AB); return ret;
        };

        auto stepm_ABC = get_stepm({0.,0.,0.}, {c_a, 0., 0.}, {c_b, c_c, 0.}); // for A B C, c_(a, b, c) = (c, b*cos(Aa), b*sin(Aa))  
        auto stepm_BCA = get_stepm({c_b, c_c, 0.}, {0.,0.,0.}, {c_a, 0., 0.}); // for B C A, c_(a, b, c) = (a, c*cos(Bb), c*sin(Bb)) 
        auto stepm_CAB = get_stepm({c_a, 0., 0.}, {c_b, c_c, 0.}, {0.,0.,0.}); // for C A B, c_(a, b, c) = (b, a*cos(Cc), a*sin(Cc))
        auto stepn_ABC = get_stepn({0.,0.,0.}, {c_a, 0., 0.}, {c_b, c_c, 0.}); // for A B C, c_(a, b, c) = (c, b*cos(Aa), b*sin(Aa))  
        auto stepn_BCA = get_stepn({c_b, c_c, 0.}, {0.,0.,0.}, {c_a, 0., 0.}); // for B C A, c_(a, b, c) = (a, c*cos(Bb), c*sin(Bb)) 
        auto stepn_CAB = get_stepn({c_a, 0., 0.}, {c_b, c_c, 0.}, {0.,0.,0.}); // for C A B, c_(a, b, c) = (b, a*cos(Cc), a*sin(Cc))
        for(int i = 0; i < m; i++) {
            r = r * (1./(_k+2*i)) * c_d;
            phi_a.step(get<0>(stepm_ABC), get<1>(stepm_ABC), get<2>(stepm_ABC));
            phi_b.step(get<0>(stepm_BCA), get<1>(stepm_BCA), get<2>(stepm_BCA));
            phi_c.step(get<0>(stepm_CAB), get<1>(stepm_CAB), get<2>(stepm_CAB));
        }
        for(int i = 0; i < n; i++) {
            r = r * (1./(_k+2*i+2*m)) * c_d;
            phi_a.step(get<0>(stepn_ABC), get<1>(stepn_ABC), get<2>(stepn_ABC));
            phi_b.step(get<0>(stepn_BCA), get<1>(stepn_BCA), get<2>(stepn_BCA));
            phi_c.step(get<0>(stepn_CAB), get<1>(stepn_CAB), get<2>(stepn_CAB));
        }
        // std::cout << "========= Info Begin ==========" << std::endl;
        // phi_a.info(); // phi_b.info(); phi_c.info();
        // std::cout << "========= Info Ends ==========" << std::endl;
        is_initialized = true;
        return *this;
    }

    double operator()(Vec3d A, Vec3d B, Vec3d C, Vec3d X) {
        Vec3d N = normalize(cross(B-A, C-A)), nz = N;
        double s = norm(cross(B-A, C-A));
        double a = norm(C-B), Aa = angle(B-A, C-A, N);
        double b = norm(A-C), Bb = angle(C-B, A-B, N); 
        double c = norm(B-A), Cc = angle(A-C, B-C, N);
        std::map<char, double> vs_ABC = {{'a', c}, {'b', b*cos(Aa)}, {'c', b*sin(Aa)}, {'d', eval(r, {{'d', 1/s}})}};
        std::map<char, double> vs_BCA = {{'a', a}, {'b', c*cos(Bb)}, {'c', c*sin(Bb)}, {'d', eval(r, {{'d', 1/s}})}};
        std::map<char, double> vs_CAB = {{'a', b}, {'b', a*cos(Cc)}, {'c', a*sin(Cc)}, {'d', eval(r, {{'d', 1/s}})}};
        // std::cout << "Fur phi_B, (a, b, c) = (" << a << ", " << c*cos(Bb) << ", " << c*sin(Bb) << ")" << std::endl; 
        Vec3d P = A + proj(X-A, nz);
        double z = dot(X-A, nz);
        auto xyl_edge = [&](Vec3d A_, Vec3d B_){
            auto nx = normalize(B_-A_);
            auto ny = cross(nz, nx); 
            double x = dot(P-A_, nx);
            double y = dot(P-A_, ny);
            double l = dist(B_, A_); 
            return Vec3d{x, y, l};
        };

        Vec3d pA = xyl_edge(A, B);
        Vec3d pB = xyl_edge(B, C);
        Vec3d pC = xyl_edge(C, A);
        double retA = 0, retB = 0, retC = 0, d = eval(r, {{'d', 1/s}});
        retA = phi_a(pA[0], pA[1], z, pA[2], vs_ABC) * d;
        retB = phi_b(pB[0], pB[1], z, pB[2], vs_BCA) * d;
        retC = phi_c(pC[0], pC[1], z, pC[2], vs_CAB) * d;
        std::cout << "(" << retA << " + " << retB << " + " << retC << ")*" << d << std::endl;
        return (retA + retB + retC);
    }
    void info(){
        phi_a.info();
        phi_b.info();
        phi_c.info();
    }

    void save(std::string folder){
        using std::to_string;
        using symbolic::all_polys::_1;
        auto _m = to_string(m), _n = to_string(n), _k = to_string(k);
        auto path = folder+"/"+_m+"_"+_n+"_"+_k; 
        save_coefi(path+"/r.bin", r);
        phi_a.save(path+"/phi_a/");
        phi_b.save(path+"/phi_b/");
        phi_c.save(path+"/phi_c/");
    }

    void load(std::string folder){
        using std::to_string;
        auto _m = to_string(m), _n = to_string(n), _k = to_string(k);
        auto path = folder+"/"+_m+"_"+_n+"_"+_k; 
        load_coefi(path+"/r.bin", r);
        phi_a.load(path+"/phi_a");
        phi_b.load(path+"/phi_b");
        phi_c.load(path+"/phi_c");
    }
};
}