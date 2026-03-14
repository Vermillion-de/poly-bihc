#pragma once
#include "types.hpp"

class MeshLoader{
public:
MatrixXd vs;
vector<int3> fs; // patches 
public:
MeshLoader(){}
MeshLoader(std::string file){ load(file);}
MeshLoader(const vector<real3>& _vs, const vector<int3>& _fs)
  : fs(_fs) { vs.resize(_vs.size(), 3);
        for(int i = 0; i < _vs.size(); i++) 
            vs.row(i) = _vs[i].transpose(); }
MeshLoader(const MatrixXd& _vs, const vector<int3>& _fs)
  : fs(_fs) {vs = MatrixXd(_vs);}
~MeshLoader(){}
public:
void load(std::string file){
    std::ifstream is(file);
    std::string start, sline;
    std::istringstream ssline;
    vector<real3> _vs;
    int idx_min = 1;
    while (getline(is, sline)) {
        if(sline.size() == 0) continue;
        ssline.str(sline); ssline >> start;
        if(start == "v"){
            double x, y, z; ssline >> x >> y >> z;
            _vs.push_back({x, y, z});
        }  
        else if(start == "f") {
            int fx, fy, fz; ssline >> fx >> fy >> fz;
            fs.push_back({fx, fy, fz});
            idx_min = std::min({fx, fy, fz, idx_min});
        }
        else if(start[0] == '#') continue;
        ssline.clear();
    } is.close();
    for(int i=0;i<fs.size();i++) 
        for(int j=0;j<fs[i].size();j++) 
            fs[i][j] -= idx_min;
    vs = to_mat(_vs);
    std::cout << "Loaded Mesh: #v=" << vs.rows() << ", #f=" << fs.size() << std::endl;
}

void save(std::string file){
    std::ofstream os(file);
    for(int i = 0; i < vs.rows(); i++)
        os << "v " << vs.row(i) << endl;
    for(auto f : fs){
        os << "f " << f[0]+1 << " " << f[1]+1 << " " << f[2]+1 << endl;
    }
    os.close();
}
public:
std::vector<double> vert(int i){
    return {vs(i, 0), vs(i, 1), vs(i, 2)};
}

std::vector<int> face(int i){
    return {fs[i][0], fs[i][1], fs[i][2]};
}

public:

void apply(std::function<real3(real3, int)> func){
    for(int i = 0; i < vs.size(); i++){
        vs.row(i) = func(vs.row(i).transpose(), i).transpose();
    }
}

vector<real3> to_verts(){
    vector<real3> vert;
    for(int i = 0; i < vs.rows(); i++)
        vert.push_back(vs.row(i).transpose());
    return vert;
}

void info(){
    std::cout << "Loaded Mesh: #v=" << vs.rows() << ", #f=" << fs.size() << std::endl;
    for(int i = 0; i < vs.rows(); i++)
        std::cout << vs.row(i) << std::endl;
    for(int i = 0; i < fs.size(); i++)
        std::cout << fs[i][0] << " " << fs[i][1] << " " << fs[i][2] << std::endl;
}

};

class BezierLoader{
public:
vector<real3> vs;       // vertices
vector<vector<int>> ps; // patches 
vector<vector<real3>> uvw; // patch uvw
vector<vector<real3>> ns; // patch normal 
vector<int> ds;         // degrees

MatrixXd vs_mat;
MatrixXd ns_mat; 
int idx = 0;

std::string normal_type="uvw";

public:
BezierLoader() {}
BezierLoader(std::string file) { load(file); }
BezierLoader(MeshLoader& m){
    m.save("./temp.obj"); // 懒
    load("./temp.obj");
    std::filesystem::remove("./temp.obj");
}
void load(std::string file){
    // for obj mesh and obj-bzr standard
    std::ifstream is(file);
    std::string start, sline;
    std::istringstream ssline;
    int idx_min = 1;
    while (getline(is, sline)) {
        if(sline.size() == 0) continue;
        ssline.str(sline); ssline >> start;
        if(start == "v"){
            double x, y, z; ssline >> x >> y >> z;
            vs.push_back({x, y, z});
        }  
        else if(start == "bzp" || start == "f") {
            int t; ps.push_back({}); 
            while (ssline >> t) ps.back().push_back(t), idx_min = std::min(idx_min, t);
            ds.push_back(std::sqrt(2*ps.back().size()+0.25)-1.5);
        }
        else if(start[0] == '#') continue;
        ssline.clear();
    } is.close();
    for(int i=0;i<ps.size();i++) 
        for(int j=0;j<ps[i].size();j++) 
            ps[i][j] -= idx_min;
    normal_type = "uvw";
    uvw = get_uvw();
    ns = get_normal();
    vs_mat = to_mat(vs, ps);
    ns_mat = to_mat(ns);
    std::cout << "Loaded Bezier: #v=" << vs.size() << ", #p=" << ps.size() 
                           << ", d="  << ds[0]     << ", #np=" << get_num_p() << std::endl;
}

void save(std::string file){
    std::ofstream os(file);
    for(auto v : vs)
        os << "v " << v.transpose() << endl;
    for(auto p : ps){
        if(p.size() > 3) os << "bzp ";
        if(p.size() == 3) os << "f ";
        for(auto idx : p) os << idx+1 << " ";
        os << endl;
    }
    os.close();
}

public:
real3 handle(int i, double u, double v, double w){
    int d = ds[i];
    real3 ret{0, 0, 0}; 
    auto idx = [&](int i, int j) { return (i*i+i)/2+j;};
    for(int iu = 0; iu <= d; iu++) for(int iv = 0; iv <= d - iu; iv++){
        int iw = d - iu - iv;
        double temp = hat(d) / (hat(iu) * hat(iv) * hat(iw)); 
        temp *= std::pow(u, iu) * std::pow(v, iv) * std::pow(w, iw);
        ret += temp * vs[ps[i][idx(iu + iv, iv)]];
    }
    return ret;
}

real3 normal(int i, double u, double v, double w){
    auto d_u = du(i, u, v, w);
    auto d_v = dv(i, u, v, w);
    auto d_w = dw(i, u, v, w);
    auto ret_uv = d_u.cross(d_v); ret_uv.normalize(); 
    auto ret_vw = d_v.cross(d_w); ret_vw.normalize();
    auto ret_wu = d_w.cross(d_u); ret_wu.normalize();
    if(!std::isnan(ret_uv[0])) return ret_uv;
    if(!std::isnan(ret_vw[0])) return ret_vw;
    if(!std::isnan(ret_wu[0])) return ret_wu;
    return ret_uv;
}

real3 du(int i, double u, double v, double w){
    int d = ds[i];
    real3 ret{0, 0, 0}; 
    auto idx = [&](int i, int j) { return (i*i+i)/2+j;};
    for(int iu = 0; iu <= d; iu++) for(int iv = 0; iv <= d - iu; iv++){
        int iw = d - iu - iv;
        double temp = hat(d) / (hat(iu) * hat(iv) * hat(iw)); 
        if(iu - 1 < 0)
            if(iw - 1 < 0) temp = 0; 
            else temp *= -iw * std::pow(1-u-v, iw-1) * std::pow(v, iv);
        else if(iw - 1 < 0)
            temp *= iu * std::pow(u, iu-1) * std::pow(v, iv);
        else
            temp *= (iu*(1-u-v) - iw*u) * std::pow(u, iu-1) * std::pow(1-u-v, iw-1) * std::pow(v, iv);
        ret += temp * vs[ps[i][idx(iu + iv, iv)]];
    }
    return ret;
}

real3 dv(int i, double u, double v, double w){
    int d = ds[i];
    real3 ret{0, 0, 0}; 
    auto idx = [&](int i, int j) { return (i*i+i)/2+j;};
    for(int iu = 0; iu <= d; iu++) for(int iv = 0; iv <= d - iu; iv++){
        int iw = d - iu - iv;
        double temp = hat(d) / (hat(iu) * hat(iv) * hat(iw)); 
        if(iv - 1 < 0)
            if(iw - 1 < 0) temp = 0; 
            else temp *= -iw * std::pow(1-u-v, iw-1) * std::pow(u, iu);
        else if(iw - 1 < 0)
            temp *= iv * std::pow(v, iv-1) * std::pow(u, iu);
        else
            temp *= (iv*(1-u-v) - iw*v) * std::pow(1-u-v, iw-1) * std::pow(v, iv-1) * std::pow(u, iu);
        ret += temp * vs[ps[i][idx(iu + iv, iv)]];
    }
    return ret;
}

real3 dw(int i, double u, double v, double w){
    int d = ds[i];
    real3 ret{0, 0, 0}; 
    auto idx = [&](int i, int j) { return (i*i+i)/2+j;};
    for(int iu = 0; iu <= d; iu++) for(int iv = 0; iv <= d - iu; iv++){
        int iw = d - iu - iv;
        double temp = hat(d) / (hat(iu) * hat(iv) * hat(iw)); 
        if(iv - 1 < 0)
            if(iw - 1 < 0) temp = 0; 
            else temp *= iw * std::pow(u, iu) * std::pow(w, iw-1);
        else if(iw - 1 < 0)
            temp *= -iv * std::pow(u, iu) * std::pow(1-u-w, iv-1);
        else
            temp *= (iw*(1-u-w) - iv*w) * std::pow(u, iu) * std::pow(w, iw-1) * std::pow(1-u-w, iv-1);
        ret += temp * vs[ps[i][idx(iu + iv, iv)]];
    }
    return ret;
}

real3 vertex(int i, int iu, int iv, int iw){
    assert(iu + iv + iw == ds[i] && "Degree should be same!");
    auto idx = [&](int i, int j) { return (i*i+i)/2+j;};
    return vs[ps[i][idx(iu + iv, iv)]];
}

real3 arc_handle(int i, real3 p, double eps=1e-5, int step=100){
    double u = 1./3., v = 1./3., d_u = 1, d_v = 1;
    real3 vec = p - handle(i, u, v, 1-u-v);
    for(int s = 0; s<step && vec.norm() > eps; s++){
        auto du_ = du(i, u, v, 1-u-v);
        auto dv_ = dv(i, u, v, 1-u-v);
        vec = p - handle(i, u, v, 1-u-v);
        d_u = vec.dot(du_) / du_.squaredNorm();
        d_v = vec.dot(dv_) / dv_.squaredNorm();
        u += d_u; v+= d_v;
    }
    return {u, v, 1-u-v};
}

public:

int get_degree(){
    int d = ds[0];
    for(auto _d : ds) if(d != _d) return 0;
    return d;
}

int get_num_p(){
    int ret = 0;
    for(auto p : ps) ret += p.size();
    return ret;
}

vector<vector<real3>> get_uvw(){
    vector<vector<real3>> uvw;
    for(int i = 0; i < ps.size(); i++){
        uvw.push_back({});
        for(int k = 0; k < ps[i].size(); k++){
            int id = ps[i][k];
            auto _uvw = arc_handle(i, vs[id]);
            if(_uvw == real3{1, 0, 0}) _uvw = real3{0.9999, 5e-5, 5e-5}; 
            if(_uvw == real3{0, 1, 0}) _uvw = real3{5e-5, 0.9999, 5e-5}; 
            if(_uvw == real3{0, 0, 1}) _uvw = real3{5e-5, 5e-5, 0.9999}; 
            uvw.back().push_back(_uvw);
        }
    }
    return uvw;
}

vector<real3> get_normal(int i){
    vector<real3> ret;
    for(int k = 0; k < ps[i].size(); k++){
        if(normal_type == "uvw"){
            auto n = normal(i, uvw[i][k][0], uvw[i][k][1], uvw[i][k][2]);
            ret.push_back(n);
        }
        else if(normal_type == "vert"){
            auto [m, n, l] = to_mnl(k, ds[i]);
            vector<tuple<int, int>> ring_edge;
            auto add_edge = [&](vector<int> mnl1, vector<int> mnl2){
                int m1 = mnl1[0], n1 = mnl1[1], l1 = mnl1[2];
                int m2 = mnl2[0], n2 = mnl2[1], l2 = mnl2[2];
                if(m1 < 0 || n1 < 0 || l1 < 0) return;
                if(m2 < 0 || n2 < 0 || l2 < 0) return;
                int i1 = to_idx(m1, n1, l1);
                int i2 = to_idx(m2, n2, l2);
                ring_edge.push_back({i1, i2});
            };
            vector<int> p1{m-1, n+1, l}, p2{m-1, n, l+1}, p3{m, n-1, l+1};
            vector<int> p4{m+1, n-1, l}, p5{m+1, n, l-1}, p6{m, n+1, l-1};
            add_edge(p1, p2), add_edge(p2, p3), add_edge(p3, p4);
            add_edge(p4, p5), add_edge(p5, p6), add_edge(p6, p1);
            auto cacu_normal = [&](int iA, int iB, int iC){
                auto A = vs[ps[i][iA]];
                auto B = vs[ps[i][iB]];
                auto C = vs[ps[i][iC]];
                real3 N = (B - A).cross(C-A); 
                N.normalize();
                return N;
            };
            real3 N_ret{0, 0, 0};
            for(auto [jB, jC] : ring_edge){
                int jA = to_idx(m, n, l);
                N_ret += cacu_normal(jA, jB, jC);
            }
            N_ret /= ring_edge.size();
            ret.push_back(N_ret);
        } 
    }
    return ret;
}

void set_ns_mat(int i){
    for(int k = 0; k < ps[i].size(); k++){
        auto n = normal(i, uvw[i][k][0], uvw[i][k][1], uvw[i][k][2]);
        ns_mat.row(idx+k) = n.transpose(); 
    }
}

vector<vector<real3>> get_normal(){
    vector<vector<real3>> ret;
    for(int i = 0; i < ps.size(); i++){
        ret.push_back(get_normal(i));
    }
    return ret;
}

void change_normal_type(std::string tgt_type){
    if(tgt_type != normal_type){
        normal_type = tgt_type;
        ns = get_normal();
        ns_mat = to_mat(ns);
    }
}

void update_vert(int id, real3 pos){
    vs[id] = pos;
    idx = 0;
    for(int i = 0; i < ps.size(); i++){
        for(int j = 0; j < ps[i].size(); j++){
            if(ps[i][j] == id){
                vs_mat.row(idx+j) = pos.transpose();
                set_ns_mat(i);
                ns[i] = get_normal(i);
                break;
            }
        }
        idx += ps[i].size();
    }
}

void update_vert(const MatrixXd& _vs_mat){
    vs_mat = MatrixXd(_vs_mat);
    vector<real3> vs_new(vs.size(), {0, 0, 0});
    vector<int> vs_n(vs.size(), 0);
    int idx = 0;
    for(auto p : ps) for(auto v_id : p){
        vs_new[v_id] += vs_mat.row(idx++).transpose();
        vs_n[v_id] += 1;
    }
    for(int i = 0; i < vs.size(); i++){
        vs_new[i] /= vs_n[i];
    }
    vs = vs_new;
    ns = get_normal();
    ns_mat = to_mat(ns);
}

void update_vert(const vector<int> idx, const vector<vector<double>>& _vs){
    for(int i = 0; i < idx.size(); i++){
        int id = idx[i];
        vs[id] = {_vs[i][0], _vs[i][1], _vs[i][2]};
    }
    vs_mat = to_mat(vs, ps);
    ns = get_normal();
    ns_mat = to_mat(ns);
}

void change_degree(int d){
    double min_dist = get_min_dist();
    for(int i = 0; i < ps.size(); i++) 
        change_degree(i, d);
    std::cout << "Change Bezier: #v=" << vs.size() << ", #p=" << ps.size() << ", d=" << ds[0] << std::endl;
    confine(min_dist / (2 * d));
    uvw = get_uvw();
    ns = get_normal();
    vs_mat = to_mat(vs, ps);
    ns_mat = to_mat(ns);
}

void confine(double threshold){
    std::map<tuple<int, int, int>, int> idx;
    auto hash = [&](real3 v) -> tuple<int, int, int>
    { return {v[0]/threshold, v[1]/threshold, v[2]/threshold}; };

    vector<real3> vs_new; int _i = 0; 
    for(int i = 0; i < vs.size(); i++)
        if(idx.find(hash(vs[i])) == idx.end()){
            vs_new.push_back(vs[i]);
            idx[hash(vs[i])] = _i++;
        }
    for(int i = 0; i < ps.size(); i++)
        for(int j = 0; j < ps[i].size(); j++)
            ps[i][j] = idx[hash(vs[ps[i][j]])];
    vs = vs_new;
}

double get_min_dist(){
    double m = 1000;
    for(auto v1 : vs) for(auto v2 : vs) if(v1 != v2){
        m = std::sqrt(std::pow(v1[0] - v2[0], 2) +
                      std::pow(v1[1] - v2[1], 2) + 
                      std::pow(v1[2] - v2[2], 2)) / 3;
    }
    return m;
}

void change_degree(int i, int d){
    assert(i < ps.size() && "Can not change degree");
    auto p = ps[i]; auto _d = ds[i];
    if(_d == 1){
        assert(_d == 1 && "Currently not support high bases ");
        map<tuple<int, int>, int> index;
        index[{0, 0}] = p[0]; real3 A = vs[p[0]];
        index[{d, 0}] = p[1]; real3 B = vs[p[1]];
        index[{d, d}] = p[2]; real3 C = vs[p[2]]; 
        for(int u=0; u<=d; u++) for(int v=0; v<=u; v++){
            if(u == 0 && v == 0) continue;
            if(u == d && v == 0) continue;
            if(u == d && v == d) continue;
            index[{u, v}] = vs.size();
            vs.push_back(A + (B-A)/d*u + (C-B)/d*v);
        }
        ps[i].resize((d+2)*(d+1)/2); int _i=0;
        for(int u=0; u<=d; u++) for(int v=0; v<=u; v++)
            ps[i][_i++] = index[{u, v}];
        ds[i] = d;
        return ;
    } else if(d == 1){
        int d0 = ds[i];
        int ia = to_idx(d0, 0, 0);
        int ib = to_idx(0, d0, 0);
        int ic = to_idx(0, 0, d0);
        ps[i] = {ps[i][ia], ps[i][ib], ps[i][ic]};
        ds[i] = d;
        return ;
    } else {
        change_degree(i, 1);
        change_degree(i, d);
        return;
    }
}

public:
void to_mesh(int i, vector<real3>& vs, vector<int3>& fs, int res=20){
    // std::cout << "Processing patch " << i << endl;
    int n0 = vs.size();
    auto idx = [&](int i, int j) { return (i*i+i)/2+j+n0;};
    for(int u=0; u<=res; u++) for(int v=0; v<=u; v++)
        vs.push_back(handle(i, (u-v)*1./res, v*1./res, 1-u*1./res));
    for(int j = 0; j < res; j++){
        for(int k = 0; k < j; k++){
            fs.push_back({idx(j, k), idx(j+1, k), idx(j+1, k+1)});
            fs.push_back({idx(j, k), idx(j+1, k+1), idx(j, k+1)});
        }
        fs.push_back({idx(j, j), idx(j+1, j), idx(j+1, j+1)});
    }    
}

MeshLoader to_mesh(int res=2){
    vector<real3> vs; 
    vector<int3> fs;
    for(int i=0; i < ps.size(); i++)
        to_mesh(i, vs, fs, res);
    std::cout << "BezierLoader To MeshLoader Size: (v = " 
              << vs.size() << ", f = " << fs.size() << ")" 
              << std::endl;
    return {vs, fs};
}

MeshLoader to_skeleton(){
    int nf = 0; for(auto d : ds) nf += d*d;
    vector<int3> fs(nf); int idx = 0;
    for(int i = 0; i < ps.size(); i++){
        auto p = ps[i];
        auto index = [&](int i, int j) { return p[(i*i+i)/2+j];};
        for(int j = 0; j < ds[i]; j++){ for(int k = 0; k < j; k++){
            fs[idx++] = {index(j, k), index(j+1, k), index(j+1, k+1)};
            fs[idx++] = {index(j, k), index(j+1, k+1), index(j, k+1)};
        }
            fs[idx++] = {index(j, j), index(j+1, j), index(j+1, j+1)};
        }    
    }
    return MeshLoader(vs, fs);
}

// BezierLoader to_degree0(){
//     int nf = ds.size();
//     vector<real3> vs;
//     vector<int3> fs;
//     for(int i = 0; i < ps.size(); i++){
//         int d = ds[i];
//         vs.push_back();
//     }
//     return MeshLoader(vs, fs);
// }


void to_skeleton(int i, int res, vector<vector<double>>& points, vector<vector<int>>& edges){
    vector<tuple<double, double, double>> bA;
    vector<tuple<double, double, double>> bB;
    vector<tuple<double, double, double>> bC;
    for(int i = 0; i <= res; i++){
        double r0 = (i * 1.0 / res);
        double r1 = (1 - r0) ;
        bA.push_back({0, r0, r1});
        bB.push_back({r1, 0, r0});
        bC.push_back({r0, r1, 0});
    }
    auto add_line = [&](const vector<tuple<double, double, double>>& b){
        int idx = points.size();
        edges.push_back({});
        for(auto [u, v, w] : b){
            auto p = handle(i, u, v, w);
            points.push_back({p[0], p[1], p[2]});
            edges.back().push_back(idx++);
        }
    };
    add_line(bA);
    add_line(bB);
    add_line(bC);
}


tuple<vector<vector<double>>, vector<vector<int>>> to_skeleton(int res){
    vector<vector<double>> points;
    vector<vector<int>> edges;
    for(int i = 0; i < ps.size(); i++){
        to_skeleton(i, res, points, edges);
    }
    std::cout << "Boundary to (" << points.size() << ", "<< edges.size() << ", " << edges[0].size() << ")" << std::endl;
    return {points, edges};
}

vector<real3> to_vertices(){
    int n = get_num_p(), idx = 0;
    vector<real3> ret(n);
    for(int i = 0; i < ps.size(); i++) 
        for(int j = 0; j < ps[i].size(); j++, idx++){
            ret[idx] = vs[ps[i][j]]; 
    }
    return ret;
}


void print(vector<vector<int>> vvs){
    for(auto vs : vvs){
        for(auto s : vs)
            std::cout << s << " ";
        std::cout << std::endl;
    }
}

MatrixXd to_coefis(){
    int n = get_num_p();
    MatrixXd coefi = MatrixXd::Constant(n, n, 0);
    vector<vector<int>> idx = get_idx();
    for(int ir = 0; ir < ps.size(); ir++) for(int jr = 0; jr < ps[ir].size(); jr++){
        double u = uvw[ir][jr][0], v = uvw[ir][jr][1], w = uvw[ir][jr][2];
        auto idx_patch = [&](int i, int j) { return (i*i+i)/2+j;};
        for(int iu = 0; iu <= ds[ir]; iu++){
            for(int iv = 0; iv <= ds[ir] - iu; iv++){
                int iw = ds[ir] - iu - iv;
                double temp = hat(ds[ir]) / (hat(iu) * hat(iv) * hat(iw)); 
                temp *= std::pow(u, iu) * std::pow(v, iv) * std::pow(w, iw);
                coefi(idx[ir][jr], idx[ir][idx_patch(iu + iv, iv)]) += temp;
            }
        }
    }
    return coefi;
}

void give_coefi(int i, real3 uvw, const vector<vector<int>>& idx, std::function<void(int c, real)> eval){
    double u = uvw[0], v = uvw[1], w = uvw[2]; int d = ds[i];
    auto idx_patch = [&](int i, int j) { return (i*i+i)/2+j;};
    for(int iu = 0; iu <= d; iu++){
        for(int iv = 0; iv <= d - iu; iv++){
            int iw = d - iu - iv;
            double temp = hat(d) / (hat(iu) * hat(iv) * hat(iw)); 
            temp *= std::pow(u, iu) * std::pow(v, iv) * std::pow(w, iw);
            eval(idx[i][idx_patch(iu + iv, iv)], temp);
        }
    }
}

vector<vector<int>> get_idx(){
    vector<vector<int>> idx; int _i = 0;
    for(auto p : ps) { 
        idx.push_back({}); 
        for(auto _p : p) 
            idx.back().push_back(_i++); 
    } 
    return idx;
}

int get_n_face(){
    int nf = 0;
    for(auto d : ds) nf += d*d;
    return nf;
}

tuple<vector<real3>, MatrixXd> samples4(){ // TODO: May bugs be with you 
    MatrixXd coefi = MatrixXd::Constant(get_n_face()*4, get_num_p(), 0); 
    vector<real3> pts; int idx_row = 0;
    vector<vector<int>> idx = get_idx();
    for(int i = 0; i < ps.size(); i++){
        auto p = ps[i];
        auto index = [](int u, int v){ return (u*u+u)/2+v;};
        auto sample = [&](int a, int b, int c){
            real3 A = vs[p[a]], B = vs[p[b]], C = vs[p[c]];
            auto coord = [&](real u, real v, real w) {
                coefi(idx_row, idx[i][a]) = u; 
                coefi(idx_row, idx[i][b]) = v; 
                coefi(idx_row, idx[i][c]) = w;
                pts.push_back(A*u + B*v + C*w);
                idx_row++;
            };
            coord(2./3., 1./6., 1./6.); coord(1./6., 2./3., 1./6.);
            coord(1./6., 1./6., 2./3.); coord(1./3., 1./3., 1./3.);
        };
        for(int j = 0; j < ds[i]; j++) {
            for(int k = 0; k < j; k++){
                sample(index(j, k), index(j+1, k), index(j+1, k+1));
                sample(index(j, k), index(j+1, k+1), index(j, k+1));
            }
            sample(index(j, j), index(j+1, j), index(j+1, j+1));
        }
    }
    return {pts, coefi};
}

tuple<vector<vector<double>>, vector<vector<int>>> to_graph(){
    int nf = 0; for(auto d : ds) nf += d*d;
    vector<int3> fs(nf); int idx = 0;
    for(int i = 0; i < ps.size(); i++){
        auto p = ps[i];
        auto index = [&](int i, int j) { return p[(i*i+i)/2+j];};
        for(int j = 0; j < ds[i]; j++){ for(int k = 0; k < j; k++){
            fs[idx++] = {index(j, k), index(j+1, k), index(j+1, k+1)};
            fs[idx++] = {index(j, k), index(j+1, k+1), index(j, k+1)};
        }
            fs[idx++] = {index(j, j), index(j+1, j), index(j+1, j+1)};
        }    
    }
    vector<vector<double>> pts;
    vector<vector<int>> lines;

    map<tuple<int, int>, bool> edges;
    auto add_line = [&](int a, int b){
        if(edges[{a, b}] == false){
            edges[{a, b}] = true;
            lines.push_back({a, b});
        }
    };

    auto add_face = [&](int a, int b, int c){
        add_line(a, b); add_line(b, c); add_line(c, a);
    };

    for(auto v : vs) pts.push_back({v[0], v[1], v[2]});
    for(auto f : fs) add_face(f[0], f[1], f[2]); 
    return {pts, lines};
}

void to_graph(std::string file){
    auto [pts, ls] = to_graph();
    save_graph(pts, ls, file);
}


void to_skeleton(std::string file){
    auto [pts, ls] = to_skeleton(50);
    save_graph(pts, ls, file);
}

void save_graph(vector<vector<double>> pts, vector<vector<int>> lines, string file){
    std::ofstream os(file);
    for(auto pt : pts)
        os << "v " << pt[0] << " " << pt[1] << " " << pt[2] << endl;
    for(auto line : lines) {
        os << "l ";
        for(auto l : line) os << l+1 << " "; 
        os << endl;
    }
    os.close();
}

public:
void info(){
    for(auto v : vs) cout << "v " << v.transpose() << endl;
    for(auto p : ps){
        cout << "bzp ";
        for(auto idx : p) cout << idx << " ";
        cout << endl;
    }
    for(auto d : ds) cout << "d " << d << endl;
}
};