#pragma once
#include <vector>
#include <math.h>
#include <cstdio>
#include <iostream>
#include <tuple>
#ifndef PI
#define PI 3.1415926535879932384626
#endif
#ifndef EPS 
#define EPS 1e-9
#endif

class Vec3d{
public:
    double x, y, z;
public:
    Vec3d(){}
    Vec3d(double _x, double _y, double _z) 
     : x(_x), y(_y), z(_z){}
public:
    const double dot(const Vec3d& b)
    { return x*b.x + y*b.y + z*b.z; }
    const Vec3d cross(const Vec3d& b)
    { return {y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x}; }
    const double norm() 
    {return std::sqrt(x*x + y*y + z*z);}
    const double dist(const Vec3d& b)
    {return std::sqrt((x-b.x)*(x-b.x) + (y-b.y)*(y-b.y) + (z-b.z)*(z-b.z));}
    const Vec3d& normalize()
    {double l = norm(); x/=l, y/=l, z/=l; return *this;}
    const double operator[](int n)
    {if(n==0) return x; if(n==1) return y; return z;}
};


std::ostream& operator<<(std::ostream& os, const Vec3d& a)
{ os << "(" << a.x << ", " << a.y << ", " << a.z << ")"; return os;}

const Vec3d operator+(const Vec3d& a, const Vec3d& b)
{ return {a.x + b.x, a.y + b.y, a.z + b.z}; }
const Vec3d operator-(const Vec3d& a, const Vec3d& b)
{ return {a.x - b.x, a.y - b.y, a.z - b.z}; }
const Vec3d operator-(const Vec3d& a)
{ return {-a.x, -a.y, -a.z}; }
const Vec3d operator*(const Vec3d& a, const double& r)
{ return {a.x*r, a.y*r, a.z*r}; }
const Vec3d operator*(const double& r, const Vec3d& a)
{ return {a.x*r, a.y*r, a.z*r}; }
const Vec3d operator/(const Vec3d& a, const double& r)
{ return {a.x/r, a.y/r, a.z/r}; }
const double dot(const Vec3d& a, const Vec3d& b)
{ return a.x*b.x + a.y*b.y + a.z*b.z; }
const Vec3d cross(const Vec3d& a, const Vec3d& b)
{ return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
const double norm(const Vec3d& a) 
{return std::sqrt(dot(a, a));}
const double dist(const Vec3d& a, const Vec3d& b)
{return std::sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));}
const Vec3d normalize(const Vec3d& a)
{double l = norm(a); return a/l; }
double angle(const Vec3d& a, const Vec3d& b)
{ return std::acos(dot(a, b) / (norm(a) * norm(b))); }
double angle(const Vec3d& a, const Vec3d& b, const Vec3d& n)
{
    if(dot(cross(a, b), n) > 0)
        return std::acos(dot(a, b) / (norm(a) * norm(b)));
    else return -std::acos(dot(a, b) / (norm(a) * norm(b)));
}

const Vec3d& operator+=(Vec3d& a, const Vec3d& b)
{ a.x += b.x, a.y += b.y, a.z += b.z; return a; }
const Vec3d& operator-=(Vec3d& a, const Vec3d& b)
{ a.x -= b.x, a.y -= b.y, a.z -= b.z; return a; }
const Vec3d& operator*=(Vec3d& a, const double& r)
{ a.x *= r, a.y *= r, a.z *= r; return a; }
const Vec3d& operator/=(Vec3d& a, const double& r)
{ a.x /= r, a.y /= r, a.z /= r; return a; }

const bool operator==(const Vec3d& a, const Vec3d& b)
{ return a.x==b.x && a.y==b.y && a.z==b.z;}

const Vec3d proj(const Vec3d& a, const Vec3d& b){
 // project a to a' _|_ b
 Vec3d n = normalize(b); return a - dot(a, n) * n;
}

const bool on_line(const Vec3d& x, const Vec3d& a, const Vec3d& b, double eps=EPS)
 { if(norm(cross(x-a, b-a)) < eps) return true; return false; }

const bool dist_line(const Vec3d& x, const Vec3d& a, const Vec3d& b){
    return norm(cross(x-a, b-a)) / norm(b-a); 
}

double sign(double z){
    if(z >= 0) return 1;
    else return -1;
}

bool is_inside(Vec3d A, Vec3d B, Vec3d C, Vec3d X){
    Vec3d N = normalize(cross(B-A, C-A));
    double s = norm(cross(B-A, C-A));
    double u = dot(cross(X-A, C-A), N) / s;
    double v = dot(cross(B-A, X-A), N) / s;
    double w = 1 - u - v;
    if( u>0 && v>0 && w>0) return true;
    return false;
}

bool on_bound(Vec3d A, Vec3d B, Vec3d C, Vec3d X){
    Vec3d N = normalize(cross(B-A, C-A));
    double s = norm(cross(B-A, C-A));
    double u = dot(cross(X-A, C-A), N) / s;
    double v = dot(cross(B-A, X-A), N) / s;
    double w = 1 - u - v;
    double eps = 1e-5;
    if( std::abs(u) < eps && v > 0 && w > 0) return true;
    if( std::abs(v) < eps && u > 0 && w > 0) return true;
    if( std::abs(w) < eps && u > 0 && v > 0) return true;
    return false;
}

const double rand(double l, double r)
{ return std::rand() * 1. / RAND_MAX * (r - l) + l;}

const Vec3d rand_vec(double l, double r)
{ return {rand(l, r), rand(l, r), rand(l, r)}; }

 template <typename T>
 const std::tuple<T, T, T> cross(const std::tuple<T,T,T>& a, const std::tuple<T,T,T>& b)
 {
    using std::get;
    return {
        get<1>(a) * get<2>(b) - get<1>(b) * get<2>(a),
        get<0>(b) * get<2>(a) - get<0>(a) * get<2>(b),
        get<0>(a) * get<1>(b) - get<0>(b) * get<1>(a),
    };
 }
