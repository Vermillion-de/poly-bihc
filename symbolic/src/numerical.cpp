#pragma once
#include <cstdio>
#include <vector>
#include <tuple>
#include <functional>

typedef std::tuple<double, double, double> real3;

std::vector<real3> sample_triangle_uv(int N){
    // TODO: use Gauss integration: https://zhilin.math.ncsu.edu/TEACHING/MA587/Gaussian_Quadrature2D.pdf
    std::vector<real3> ret; //u, v, weight
    double w0 = 2./(N*N)/2, w1 = 1./(N*N)/2;
    for(int i = 0; i < N; i++) for(int j = 0; j < N-i; j++){
        if (j != N-i-1)
            ret.push_back({(i+0.5)/N, (j+0.5)/N, w0});
        else ret.push_back({(i+0.5)/N, (j+0.5)/N, w1});
    }
    return ret;
}

double num_int_triangle_uv(std::function<double(double, double)> func, int N=100)
{ // Notes: func should self-complete, contain change of area.
    auto samples = sample_triangle_uv(N);
    double ret = 0;
    for(auto sample : samples){
        double u = std::get<0>(sample);
        double v = std::get<1>(sample);
        double w = std::get<2>(sample);
        ret += func(u, v) * w;
    }
    return ret;
} 

double num_int(std::function<double(double)> func, double l, double r, int N=100){
    double ret = 0, dx = (r - l) / N;
    for(int i = 0; i < N; i++)
        ret += dx * func(l + dx * i + dx / 2);
    return ret;
}

template <typename Vec3D> 
double num_diff_3d(std::function<double(Vec3D)> f, Vec3D p, Vec3D d, double h=1e-2) {
    auto eval = [&](double h){
        return f({p[0] + h*d[0], p[1] + h*d[1], p[2] + h*d[2]});
    };
    return (eval(-2*h)-eval(2*h)+8*eval(h)-8*eval(-h)) / (12*h);
}

double num_diff(std::function<double(double)> f, double x, double h=1e-2) { 
    return (f(x-2*h)-f(x+2*h)+8*f(x+h)-8*f(x-h)) / (12*h); 
}

template <typename Vec3D>
Vec3D num_nabla(std::function<double(Vec3D)> f, Vec3D X_, double h = 1e-2){
    double nabla_x = num_diff([&](double dx){ return f(X_+Vec3D{dx, 0, 0});}, 0, h);
    double nabla_y = num_diff([&](double dy){ return f(X_+Vec3D{0, dy, 0});}, 0, h);
    double nabla_z = num_diff([&](double dz){ return f(X_+Vec3D{0, 0, dz});}, 0, h);
    return Vec3D{nabla_x, nabla_y, nabla_z};
}