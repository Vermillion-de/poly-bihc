#pragma once
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <tuple>
#include <map>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <filesystem>
#include <sstream>
#include <zlib.h>

#include <Eigen/Core>
#include <Eigen/Dense>

using std::cout;
using std::endl;
using std::vector;
using std::tuple;
using std::map;
using std::string;

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;

using real = double;
using real3 = Eigen::Vector3d;
using int3 = Eigen::Vector3i;

#ifndef PI 
#define PI 3.1415926535897932
#endif

Eigen::Matrix3d gen_rot(){
    using namespace Eigen;
    Matrix3d ret = Matrix3d::Random();
    JacobiSVD<Matrix3d> svd(ret, ComputeFullU | ComputeFullV);
    Matrix3d ret_ = svd.matrixV() * svd.matrixU().transpose(); 
    return ret_ / ret_.determinant();
}

Eigen::Matrix3d svd_mat(MatrixXd p){
    using namespace Eigen;
    JacobiSVD<Matrix3d> svd(p, ComputeFullU | ComputeFullV);
    Matrix3d ret_ = svd.matrixV() * svd.matrixU().transpose(); 
    return ret_ / ret_.determinant();
}

template<typename D>
vector<D> serialize(const vector<vector<D>>& vvd){
    vector<D> ret;
    for(auto vd : vvd) for(auto d : vd) ret.push_back(d);
    return ret;
}

MatrixXd to_mat(const vector<vector<real3>>& vvd){
    int n = 0; for(auto vd : vvd) n += vd.size();
    MatrixXd ret(n, 3); int idx = 0;
    for(auto vd : vvd) for(auto d : vd) 
        ret.row(idx++) << d.transpose();
    return ret;
}

MatrixXd to_mat(const vector<real3>& vd){
    int n = vd.size();
    MatrixXd ret(n, 3); int idx = 0;
    for(auto d : vd) 
        ret.row(idx++) << d.transpose();
    return ret;
}

MatrixXd to_mat(vector<vector<double>> vs, int rows, int cols){
    MatrixXd ret(rows, cols);
    for(int i = 0; i < rows; i++) for(int j = 0; j < cols; j++){
        ret(i, j) = vs[i][j];
    }
    return ret;
}


vector<real3> to_vec(const MatrixXd& vd){
    vector<real3> ret(vd.rows());
    for(int i = 0; i < vd.rows(); i++)
        ret[i] = vd.row(i).transpose();
    return ret;
}

MatrixXd to_mat(const vector<real3>& vp, const vector<vector<int>> vvid){
    int n = 0; for(auto vid : vvid) n += vid.size();
    MatrixXd ret(n, 3); int idx = 0;
    for(auto vid : vvid) for(auto id : vid) 
        ret.row(idx++) << vp[id].transpose();
    return ret;
}

MatrixXd to_mat(const MatrixXd& vp, const vector<int> vid){
    MatrixXd ret(vid.size(), 3); int idx = 0;
    for(auto id : vid) 
        ret.row(idx++) << vp.row(id);
    return ret;
}

vector<vector<double>> mat3(MatrixXd m, int r0, int c0){
    std::cout << "Matrix: (" << m.rows() << ", " << m.cols() << ")" << std::endl; 
    std::cout << "mat33 from (" << r0 << ", " << c0 << ")" << std::endl;
    vector<vector<double>> ret = {
        {m(r0+0, c0+0), m(r0+0, c0+1), m(r0+0, c0+2), 0},
        {m(r0+1, c0+0), m(r0+1, c0+1), m(r0+1, c0+2), 0}, 
        {m(r0+2, c0+0), m(r0+2, c0+1), m(r0+2, c0+2), 0}, 
        {0, 0, 0, 1}
    };
    std::cout << "Fine Here" << std::endl;
    return ret;
}

real rand(real l, real r) {
     return (r - l) / RAND_MAX * std::rand() + l;
}

real3 rand_uvw(){
    real u = rand(0, 1);
    real v = rand(0, 1);
    return {u, (1-u)*v, (1-u)*(1-v)};
}

void info_mat(std::string name, const MatrixXd& m){
    std::cout << name << ": " << m.rows() << ", " << m.cols() << std::endl;
}

double hat(int n){
    double ret = 1;
    for(int i = 1; i <= n; i++) ret *= i;
    return ret;
}
double filter(double d, double bound=1e3){
    if(std::isnan(d) || std::isinf(d) || d > bound || d < -bound)
        return 0;
    return d;
}

// xxx.eigen
void save(std::string file, const MatrixXd& d, string algorithm="txt", bool verbose=true){
    if(verbose) std::cout << "Saving into " << file << std::endl;
    if(algorithm == "txt"){
        std::ofstream os(file);
        os << d.rows() << " " << d.cols() << std::endl;
        MatrixXd d_ = MatrixXd(d);
        for(int i = 0; i < d.rows(); i++) for(int j = 0; j < d.cols(); j++){
            d_(i, j) = filter(d(i, j));
        }
        os << d_ << std::endl;
        os.close();
    } else if(algorithm == "zlib") {
            gzFile of = gzopen(file.c_str(), "wb9");
            if(!of)  throw std::runtime_error("Write "+ file + " error!");
            uint r = d.rows(), c = d.cols(); 
            gzwrite(of, &r, sizeof(uint));
            gzwrite(of, &c, sizeof(uint));
            gzwrite(of, d.data(), sizeof(double)*r*c);
            gzclose(of);
    } else if(algorithm == "binary") {
            std::ofstream os(file, std::ios::out | std::ios::binary | std::ios::trunc);
            uint r = d.rows(), c = d.cols();
            os.write((char*)(&r), sizeof(uint));
            os.write((char*)(&c), sizeof(uint));
            os.write((char*)d.data(), r * c * sizeof(double));
    }
}

// xxx.eigen
void load(std::string file, MatrixXd& mat, string algorithm="txt", bool verbose=true){
    uint r=0, c=0;
    if(verbose) std::cout << "Loading " << file << ", (" << r << ", " << c << ")." << std::endl;
    if(!std::filesystem::exists(file)) {
        std::cout << "File " << file << " not exist!" << std::endl;
        mat = MatrixXd::Constant(0, 0, 0);
        return ;
    }
    if(algorithm == "txt"){
        std::ifstream is(file);
        is >> r >> c; mat.resize(r, c);
        for(int i = 0; i < r; i++) for(int j = 0; j < c; j++){
            is >> mat(i, j);
        }
        // cout << file << endl << mat << endl;
        is.close(); 
    } else if(algorithm == "zlib"){
        gzFile in = gzopen(file.c_str(), "rb");
        if(!in) throw std::runtime_error("Can not open " + file + " !");
        gzread(in, &r, sizeof(uint));
        gzread(in, &c, sizeof(uint));
        mat.resize(r, c);
        gzread(in, mat.data(), sizeof(double) * r * c);
        gzclose(in);
    } else if(algorithm == "binary"){
        std::ifstream is(file, std::ios::in | std::ios::binary);
        is.read((char*)(&r), sizeof(uint));
        is.read((char*)(&c), sizeof(uint));
        mat.resize(r, c);
        is.read((char*)mat.data(), r * c * sizeof(double));
        is.close();
    }
    if(verbose) std::cout << "Loading " << file << ", (" << r << ", " << c << ")." << std::endl;
}

void save(std::string file, const vector<real3>& v){
    std::cout << "Saving into " << file << std::endl;
    std::ofstream os(file);
    for(auto d : v) os << d.transpose() << std::endl;
    os.close();
}

void load(std::string file, vector<real3>& vec){
    std::ifstream is(file);
    vec.clear();
    real x, y, z;
    while ((is >> x >> y >> z))
        vec.push_back({x, y, z});
    std::cout << "Load from " << file << " of " << vec.size() << std::endl;
    is.close(); 
}

string base_name(string path){
    return std::filesystem::path(path).filename().string();
}

vector<string> list_dir(string dir){
    vector<string> files;
    using std::filesystem::create_directories;
    using std::filesystem::directory_iterator;
    if(!std::filesystem::exists(dir)){
        std::cout << "Can not listdir " << dir << ", since non-existence" << std::endl;
        create_directories(dir);
    }
    for (const auto& entry : directory_iterator(dir)) {
        std::cout << "Found " << entry.path().filename().string() << " in " << dir << std::endl; 
        files.push_back(entry.path().filename().string());
    }
    return files;
}

namespace utils
{
    int idx = 0;
    void info(std::string s) { std::cout << "[INFO] " << s << std::endl; }
    std::string prefix(std::string name){
        auto now = std::chrono::system_clock::now();
        std::time_t time = std::chrono::system_clock::to_time_t(now);
        std::tm tm = *std::localtime(&time);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%H:%M:%S");
        return "[INFO " + oss.str() + " ] " + name;
    }
    std::string prefix(){
        auto now = std::chrono::system_clock::now();
        std::time_t time = std::chrono::system_clock::to_time_t(now);
        std::tm tm = *std::localtime(&time);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%H:%M:%S");
        return oss.str();
    }
} // namespace 

vector<vector<double>> to_double_(const vector<real3>& t){
    vector<vector<double>> ret;
    for(auto v : t) ret.push_back({v[0], v[1], v[2]});
    return ret;
}

vector<vector<double>> to_double(vector<real3> t){
    vector<vector<double>> ret;
    for(auto v : t) ret.push_back({v[0], v[1], v[2]});
    return ret;
}

vector<vector<double>> to_double(const MatrixXd& t){
    vector<vector<double>> ret;
    for(int i = 0; i < t.rows(); i++) 
        ret.push_back({t(i,0), t(i, 1), t(i, 2)});
    return ret;
}

tuple<int, int, int> to_mnl(int j, int d){
    int u = std::sqrt(2*j+0.25)-0.5;
    int v = j - ((u+1)*u)/2;
    int m = d - u;
    int n = u - v;
    int l = v;
    return {m, n, l};
}

int to_idx(int m, int n, int l){
    int u = n + l;
    int v = l; 
    return ((u+1)*u)/2+v;
}

MatrixXd inv(const MatrixXd& m){
    using namespace Eigen;
    double epsilon = std::numeric_limits<double>::epsilon(); 
    BDCSVD<MatrixXd> svd(m, ComputeFullU | ComputeFullV);
    const auto& S = svd.singularValues();
    MatrixXd sigma_pinv(m.cols(), m.rows());
    sigma_pinv.setZero();
    double tolerance = epsilon * std::max(m.rows(), m.cols()) * S(0); 
    for(int i = 0; i < S.size(); ++i){
        if(S(i) > tolerance)    
            sigma_pinv(i, i) = 1./S(i);
    }
    std::cout << "Inv: (" << m.rows() << ", " << m.cols() << ")" << std::endl;
    return svd.matrixV() * sigma_pinv * svd.matrixU().transpose();
}

MatrixXd rowwise_prod(const MatrixXd& a, const VectorXd& b){
    assert(a.rows() == b.rows());
    return a.array().colwise() * b.array();
}

MatrixXd rowwise_divide(const MatrixXd& a, const VectorXd& b){
    assert(a.rows() == b.rows());
    return a.array().colwise() / b.array();
}

MatrixXd solve(const MatrixXd& A, const MatrixXd& b, string solver="inv"){
    std::cout << "solving: (" << A.rows() << ", " << A.cols() << ")" << std::endl;
    if(solver == "LDLT"){
        Eigen::LDLT<MatrixXd> ldlt(A);
        return ldlt.solve(b);
    }
    else if(solver == "inv") {
        return inv(A) * b;
    }
    return {};
}

void filter(MatrixXd& A, MatrixXd& b, real threshold=5e3){
    vector<int> rows;
    for(int i = 0; i < A.rows(); i++){
        if(!(A.row(i).array().abs() > threshold).any() && 
           !(b.row(i).array().abs() > threshold).any()) {
            rows.push_back(i);
        }
    }
    MatrixXd A_(rows.size(), A.cols());
    MatrixXd b_(rows.size(), b.cols());
    for(int i = 0; i < rows.size(); i++){
        A_.row(i) = A.row(rows[i]);
        b_.row(i) = b.row(rows[i]);
    }
    A = A_; b = b_;
}

MatrixXd optimize(const vector<MatrixXd>& As, const vector<MatrixXd>& bs, int cA, int cb, vector<real> lambda, int filt=100){
    // optimze for  sum lambda||Ax - b||^2
    MatrixXd At_A = MatrixXd::Constant(cA, cA, 0);
    MatrixXd At_b = MatrixXd::Constant(cA, cb, 0);
    for(int i = 0; i < lambda.size(); i++){
        MatrixXd A = As[i];
        MatrixXd b = bs[i];
        assert(A.rows() == b.rows());
        assert(A.cols() == cA && b.cols() == cb);
        filter(A, b, filt);
        At_A += A.transpose() * A * lambda[i];
        At_b += A.transpose() * b * lambda[i];
    }
    return solve(At_A, At_b);
}

MatrixXd optimize(const vector<MatrixXd>& mat_l, const vector<MatrixXd>& mat_m, tuple<int, int> c_lm,
                  const vector<MatrixXd>& mat_v, const vector<MatrixXd>& mat_n, tuple<int, int> c_vn, 
                  vector<real> lambda){
    // optimize for mat_l * L + mat_m * M = mat_v * V + mat_n * N, weighted by lambda 
    auto [cl, cm] = c_lm; int clm = cl + cm;
    auto [cv, cn] = c_vn; int cvn = cv + cn;
    MatrixXd At_A = MatrixXd::Constant(clm, clm, 0);
    MatrixXd At_b = MatrixXd::Constant(clm, cvn, 0);
    for(int i = 0; i < lambda.size(); i++){
        const MatrixXd &m_l = mat_l[i], &m_m = mat_m[i];
        const MatrixXd &m_v = mat_v[i], &m_n = mat_n[i];
        int rl = m_l.rows(), rm = m_m.rows();
        int rv = m_v.rows(), rn = m_n.rows();
        assert(rv == rl && rl == rm && m_l.cols() == cl && m_m.cols() == cm);
        assert(rv == rl && rv == rn && m_v.cols() == cv && m_n.cols() == cn);
        MatrixXd A = MatrixXd::Constant(rl, clm, 0); A << m_l, m_m;
        MatrixXd b = MatrixXd::Constant(rl, cvn, 0); b << m_v, m_n; 
        At_A += A.transpose() * A * lambda[i];
        At_b += A.transpose() * b * lambda[i];
    };
    return solve(At_A, At_b); 
}

tuple<MatrixXd, MatrixXd> vnlm_to_vn(const vector<MatrixXd>& mat_vnlm, const MatrixXd& mat_vn){
    // mat_vn * (v; n) = (l; m)
    const MatrixXd &m_v = mat_vnlm[0], &m_n = mat_vnlm[1];
    const MatrixXd &m_l = mat_vnlm[2], &m_m = mat_vnlm[3];
    int rv = m_v.rows(), rn = m_n.rows(), cv = m_v.cols(), cn = m_n.cols();
    int rl = m_l.rows(), rm = m_m.rows(), cl = m_l.cols(), cm = m_m.cols();
    assert(rv == rn && rl == rm && rl == rv);
    MatrixXd m_vn(rv, cv+cn); m_vn << m_v, m_n;
    MatrixXd m_lm(rl, cl+cm); m_lm << m_l, m_m;
    m_vn = m_vn + m_lm * mat_vn;
    return {
        m_vn.block(0, 0,  rv, cv),
        m_vn.block(0, cv, rv, cn),
    };
}

tuple<MatrixXd, MatrixXd> lm_to_vn(const vector<MatrixXd>& coefi_vn, const MatrixXd& mat_vn){
    // mat_vn * (v; n) = (l; m)
    const MatrixXd &m_v = coefi_vn[0], &m_n = coefi_vn[1];
    int rn = m_n.rows(), rv = m_v.rows();
    int cn = m_n.cols(), cv = m_v.cols();
    assert(rv == rn);
    MatrixXd m_vn(rv, cv+cn); m_vn << m_v, m_n;
    m_vn = m_vn * mat_vn;
    return {
        m_vn.block(0, 0,  rv, cv),
        m_vn.block(0, cv, rv, cn),
    };
}