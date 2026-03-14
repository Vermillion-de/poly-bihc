#include "types.hpp"
#include "./src/compute.h"
#include "surface.hpp"

void compute(const BezierLoader& _cage, BezierLoader& cage, const vector<real3>& vert, 
    MatrixXd& g1_, MatrixXd& g1n, MatrixXd& g2_, MatrixXd& g2n, int N=10, bool unify=true, bool verbose=true){
    // Cacu for general coordinates 
    if(verbose) std::cout << utils::prefix("d0: Preparing Coordaintes: ") << cage.get_num_p() * vert.size() << std::endl;
    vector<vector<real>> ABC(_cage.ps.size()), X(vert.size());
    int n = cage.get_num_p() * vert.size();
    vector<vector<int>> abc_x(n);
    vector<vector<int>> mnl(n);
    for(int i = 0; i < _cage.ps.size(); i++){
        real3 a = _cage.vs[_cage.ps[i][0]];
        real3 b = _cage.vs[_cage.ps[i][1]];
        real3 c = _cage.vs[_cage.ps[i][2]];
        ABC[i] = { a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]};
    }
    for(int i = 0; i < vert.size(); i++){
        double r = 1 + rand(0, 1) * 1e-4;
        // if(cuda::use_num == true) r = 1;
        r = 1;
        X[i] = {vert[i][0] * r, vert[i][1] * r, vert[i][2] * r};
    }


    for(int k = 0, idx=0; k < vert.size(); k++) for(int i = 0; i < cage.ps.size(); i++){
        for(int j = 0; j < cage.ps[i].size(); j++, idx++){
            auto [m, n, l] = to_mnl(j, cage.ds[i]);
            abc_x[idx] = {i, k};
            mnl[idx] = {m, n, l};
        }
    }


    if(verbose) std::cout << utils::prefix("d0: Begin Computing") << std::endl;
    cuda::quad::init(N);
    cuda::phi::init(cage.get_degree(), -1, 3, "/home/qiz/Desktop/project/bih_test/cuda/poly");
    vector<real> _1_ = cuda::interface(ABC, X, abc_x, mnl, 1);
    vector<real> _3_ = cuda::interface(ABC, X, abc_x, mnl, 3);
    vector<real> __1 = cuda::interface(ABC, X, abc_x, mnl, -1);
    cuda::quad::end();
    // cuda::phi::end(); // TODO

    if(verbose) std::cout << utils::prefix("d0: Coordinates Arrange") << std::endl;
    // std::cout << "Fine " << __LINE__ << std::endl;
    g1_ = MatrixXd::Constant(vert.size(), cage.get_num_p(), 0);
    g1n = MatrixXd::Constant(vert.size(), cage.get_num_p(), 0);
    g2_ = MatrixXd::Constant(vert.size(), cage.get_num_p(), 0);
    g2n = MatrixXd::Constant(vert.size(), cage.get_num_p(), 0);
    // std::cout << "Fine " << __LINE__ << std::endl;
    for(int k = 0, idx = 0; k < vert.size(); k++) for(int i = 0, idx_col=0; i < cage.ps.size(); i++){
        real3 A = _cage.vs[_cage.ps[i][0]];
        real3 B = _cage.vs[_cage.ps[i][1]];
        real3 C = _cage.vs[_cage.ps[i][2]];
        real3 N = (C - A).cross(B - A); N.normalize();
        real3 X = vert[k];
        real  z = N.dot(X - A);
        // std::cout << "Fine " << __LINE__ << std::endl;
        for(int j = 0; j < cage.ps[i].size(); j++, idx_col++, idx++){
            auto [m, n, l] = to_mnl(j, cage.ds[i]);
            real _d_  = {hat(cage.ds[i]) / (hat(m) * hat(n) * hat(l))};
            real _d1_ = _d_ / (4 * PI); 
            real _d2_ = _d_ / (8 * PI); 
            g1n(k, idx_col) += _d1_ *  z * filter(_3_[idx]);
            g1_(k, idx_col) += _d1_ *  1 * filter(_1_[idx]);
            g2n(k, idx_col) += _d2_ * -z * filter(_1_[idx]);
            g2_(k, idx_col) += _d2_ *  1 * filter(__1[idx]);
        }
    }

    // Important!
    if(unify)
    for(int idx_row = 0; idx_row < vert.size(); idx_row++){
        // std::cout << "Fine " << __LINE__ << std::endl;
        double r = filter(1 / g1n.row(idx_row).sum());
        std::cout << "Ratio " << r << std::endl;
        g1n.row(idx_row) *= r;
        g1_.row(idx_row) *= r;
        g2n.row(idx_row) *= r;
        g2_.row(idx_row) *= r;
    }

}

void green( const MatrixXd& G1n, const MatrixXd& G1_, 
            const MatrixXd& G2n, const MatrixXd& G2_, const std::string dir ){
    const MatrixXd& coefi_v_green = G1n;
    const MatrixXd& coefi_n_green = G1_;
    save(dir+"/coefi_v_green.csv", coefi_v_green);
    save(dir+"/coefi_n_green.csv", coefi_n_green);
}

void print(const MatrixXd& m, std::string name){
    std::cout << name << ": (" << m.rows() << ", " << m.cols() << ")" << std::endl;
}

void bih_method1( const MatrixXd& g1n, const MatrixXd& g1_, 
                  const MatrixXd& g2n, const MatrixXd& g2_, const MatrixXd& Ig,
                  const MatrixXd& G1n, const MatrixXd& G1_, 
                  const MatrixXd& G2n, const MatrixXd& G2_, const std::string dir ){
    print(Ig, "Ig");
    print(g1n, "g1n");
    MatrixXd A = inv(Ig - g1n) * g1_;
    MatrixXd B = g2n * A + g2_;
    MatrixXd C = G2n * A + G2_;
    MatrixXd D = C * inv(B);
    MatrixXd coefi_v_bih = G1n + D * (Ig - g1n);
    MatrixXd coefi_n_bih = G1_ - D * g1_;
    save(dir+"/coefi_v_bih.csv", coefi_v_bih);
    save(dir+"/coefi_n_bih.csv", coefi_n_bih);
}

void bih_method2( const MatrixXd& g1n, const MatrixXd& g1_, 
                  const MatrixXd& g2n, const MatrixXd& g2_, const MatrixXd& Ig,
                  const MatrixXd& G1n, const MatrixXd& G1_, 
                  const MatrixXd& G2n, const MatrixXd& G2_, const std::string dir ){
    MatrixXd A = inv(g1_) * (Ig - g1n); // d = Ac
    MatrixXd B = g2_ * A + g2n, B_inv = inv(B);
    MatrixXd C = G2n + G2_ * A, D = C * B_inv;
    MatrixXd coefi_v_bih = G1n + D * (Ig - g1n);
    MatrixXd coefi_n_bih = G1_ - D * g1_;
    save(dir+"/coefi_v_bih.csv", coefi_v_bih);
    save(dir+"/coefi_n_bih.csv", coefi_n_bih);
}

MatrixXd bih_mat_vn( const MatrixXd& g1n, const MatrixXd& g1_, 
                     const MatrixXd& g2n, const MatrixXd& g2_, 
                     const MatrixXd& Ig,  const vector<real> lamb){
    const int c_v = g1n.cols(), c_n = g1_.cols();
    const int c_l = g2n.cols(), c_m = g2_.cols();
    const MatrixXd mv_0 = MatrixXd::Zero(g1n.rows(), g1n.cols()); 
    const MatrixXd mn_0 = MatrixXd::Zero(g1_.rows(), g1_.cols()); 
    const MatrixXd m_vn = optimize({g1n-Ig, g2n}, {g1_, g2_}, {c_l, c_m}, {mv_0, Ig-g1n}, {mn_0, -g1_}, {c_v, c_n}, lamb);
    return m_vn;
}

void bih_reg_02(const MatrixXd& m_vn, const MatrixXd& G1n, const MatrixXd& G1_, 
                const MatrixXd& G2n, const MatrixXd& G2_, const std::string dir ){
    // Ig*L = g1n*L + g1_*M
    // Ig*V = g1n*V + g1_*N + g2n*L + g2_*M
    auto [coefi_v, coefi_n] = vnlm_to_vn({G1n, G1_, G2n, G2_}, m_vn);
    save(dir+"/coefi_v_bih_reg.csv", coefi_v);
    save(dir+"/coefi_n_bih_reg.csv", coefi_n);
}

void save_cage( const MatrixXd& g1n, const MatrixXd& g1_, const MatrixXd& g2n, const MatrixXd& g2_, const MatrixXd& Ig, const std::string dir){
    auto folder = dir+"/cage_cond";
    std::cout << "Saving (g1n, g1_, g2n, g2_) into " << folder << std::endl;
    std::filesystem::create_directories(folder);
    save(folder+"/Ig.csv", Ig);
    save(folder+"/g1n.csv", g1n); save(folder+"/g1_.csv", g1_); 
    save(folder+"/g2n.csv", g2n); save(folder+"/g2_.csv", g2_); 
}

void save_gG(const MatrixXd& g1n, const MatrixXd& g1_, const MatrixXd& g2n, const MatrixXd& g2_, const MatrixXd& Ig,
             const MatrixXd& G1n, const MatrixXd& G1_, const MatrixXd& G2n, const MatrixXd& G2_, const std::string dir ){
    std::filesystem::create_directories(dir+"/g");
    std::filesystem::create_directories(dir+"/G");
    std::cout << "Saving (g1n, g1_, g2n, g2_)  and (G1n, G1_, G2n, G2_) into [" << dir << "/g/, "<< dir << "/G/]" << std::endl;
    save(dir+"/g/Ig.csv", Ig);
    save(dir+"/g/g1n.csv", g1n); save(dir+"/g/g1_.csv", g1_); 
    save(dir+"/g/g2n.csv", g2n); save(dir+"/g/g2_.csv", g2_); 
    save(dir+"/G/G1n.csv", G1n); save(dir+"/G/G1_.csv", G1_); 
    save(dir+"/G/G2n.csv", G2n); save(dir+"/G/G2_.csv", G2_); 
}

int main(int argc, char** argv){
    int tgt_degree;
    std::string dir;
    bool compute_mesh;
    bool compute_p2p;
    bool as_linear = false;
    bool use_samples4 = false;

    { // argparser
        if(argc == 3){
            dir = std::string(argv[1]);
            tgt_degree = std::stoi(std::string(argv[2]));
            compute_mesh = true;
            compute_p2p = true;
        }
        else if(argc == 4){
            dir = std::string(argv[1]);
            tgt_degree = std::stoi(std::string(argv[2]));
            auto arg3 = std::string(argv[3]);
            if(arg3 == "compute_mesh") compute_mesh = true, compute_p2p = false;
            if(arg3 == "compute_p2p") compute_mesh = false, compute_p2p = true;
            if(arg3 == "as_linear") compute_mesh = true, compute_p2p = true, as_linear=true;
            if(arg3 == "use_samples4") compute_mesh = true, compute_p2p = true, use_samples4=true;
        } else {
            std::cout << "Usage: ./compute [PATH_TO_DATA] [Tgt_Degree] [optional: \"compute_mesh\", \"compute_p2p\" (only)]" << std::endl;
            return -1;
        }
    }
    BezierLoader _cage(dir+"/cage0.obj"); 
    BezierLoader cage(dir+"/cage0.obj");
    cage.change_degree(tgt_degree);

    if(as_linear){
        MeshLoader temp = cage.to_skeleton();
        cage = BezierLoader(temp);
        _cage = cage;
    }

    _cage.save(dir+"/cage0.obj");
    cage.save(dir+"/cage1.obj");
    MeshLoader mesh(dir+"/mesh.obj"); 

    MatrixXd g1n, g1_, g2n, g2_;
    int n = cage.get_num_p();
    vector<real3> pts;
    MatrixXd Ig;
    if(!use_samples4){
        Ig = MatrixXd::Identity(n, n);
        pts = cage.to_vertices(); //TODO Sampling on vertices
    }
    else {
        auto [_pts, _Ig] = cage.samples4();
        Ig = _Ig; pts = _pts;
    }

    compute(_cage, cage, pts, g1_, g1n, g2_, g2n);

    auto yield_coord = [&](std::string name, const vector<real3>& vs){
        auto folder = dir + "/" + name;
        std::filesystem::create_directories(folder);
        MatrixXd G1n, G1_, G2n, G2_;
        compute(_cage, cage, vs, G1_, G1n, G2_, G2n);
        save_gG(g1n, g1_, g2n, g2_, Ig,
                G1n, G1_, G2n, G2_, folder);
        bih_method1(g1n, g1_, g2n, g2_, Ig,
                  G1n, G1_, G2n, G2_, folder);
        green(G1n, G1_, G2n, G2_, folder);

    };

    if(compute_mesh) yield_coord("mesh", mesh.to_verts());
    return 0;
}
