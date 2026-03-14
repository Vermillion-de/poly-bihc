#include <stdio.h>
#include <iostream>
#include <random>
#include "numerical.cpp"
#include "src/utils.hpp"
#include "src/phi_2d.hpp"
#include "src/phi_3d.hpp"

using namespace std;

double cacu_numerical(Vec3d A, Vec3d B, Vec3d C, Vec3d X, int m, int n, int k, int N=50)
{
    Vec3d AB = B - A, AC = C - A;
    double s_abc = norm(cross(AB, AC));
    auto func = [&](double u, double v){
        Vec3d Y = A + AB * u + AC * v;
        return std::pow(u, m) * std::pow(v, n) / pow(norm(Y-X), k) * s_abc;
    };
    return num_int_triangle_uv(func, N);
}

void save_main(){
    int N = 3;
    std::string dir="./poly";
    for(int m = 0; m <= N; m++) for(int n = 0; n <= N; n++){
        std::string folder = dir+"/poly2d";
        symbolic::phi_mnk_2d phi1(m, n, 1); phi1.save(folder); // phi1.info();
        symbolic::phi_mnk_2d phi3(m, n, 3); phi3.save(folder); // phi3.info();
        symbolic::phi_mnk_2d phi_1(m, n, -1); phi_1.save(folder); // phi_1.info();
        symbolic::phi_mnk_2d phi_3(m, n, -3); phi_3.save(folder);  // phi_3.info();
    }
    for(int m = 0; m <= N; m++) for(int n = 0; n <= N; n++){
        std::string folder = dir+"/poly3d";
        symbolic::phi_mnk phi1(m, n, 1); phi1.save(folder); // phi1.info();
        symbolic::phi_mnk phi3(m, n, 3); phi3.save(folder); // phi3.info();
        symbolic::phi_mnk phi_1(m, n, -1); phi_1.save(folder); // phi_1.info();
        symbolic::phi_mnk phi_3(m, n, -3); phi_3.save(folder);  // phi_3.info();
    }
}

int main(){
    save_main();
    return 0;
}
