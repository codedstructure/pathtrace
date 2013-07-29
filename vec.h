#ifndef INCLUDE_PATHTRACE_VEC_H
#define INCLUDE_PATHTRACE_VEC_H

#include <cmath>

// TODO: make this a template, with 'double' and CARD being parameterized.

struct vec3 {
    enum {CARD=3};

    union { double x, s, r; };
    union { double y, t, g; };
    union { double z, u, b; };
    vec3(double x_=0, double y_=0, double z_=0) { x=x_; y=y_; z=z_; }
    double operator[] (const unsigned int idx) {
        if (idx < CARD) {
            return (&x)[idx];
        }
        return 0;
    }
    vec3 operator+(const vec3 &b) const { return vec3(
            x+b.x,
            y+b.y,
            z+b.z); }
    vec3 operator-(const vec3 &b) const { return vec3(
            x-b.x,
            y-b.y,
            z-b.z); }
    vec3 operator*(const vec3 &b) const { return vec3(
            x*b.x,
            y*b.y,
            z*b.z); }
    vec3 operator*(double b) const { return vec3(
            x*b,
            y*b,
            z*b); }
    vec3 operator/(const vec3 &b) const { return vec3(
            x/b.x,
            y/b.y,
            z/b.z); }
    vec3 operator/(double b) const { return vec3(
            x/b,
            y/b,
            z/b); }
    double dot(const vec3 &b) const { return x*b.x+y*b.y+z*b.z; } // cross:
    vec3 operator%(const vec3&b) const {return vec3(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
};

// for const multiplier on the LHS
vec3 operator*(const double t, const vec3 &v)
{
    return vec3(t*v.x, t*v.y, t*v.z);
}

double dot(const vec3 &p0, const vec3 &p1) {
    return p0.x*p1.x + p0.y*p1.y + p0.z*p1.z;
}

template <typename VT>
double length(const VT &x) {
    return sqrt(dot(x, x));
}

template <typename VT>
VT normalize(const VT x) {
    return x / length(x);
}

template <typename VT>
VT reflect(const VT I, const VT N) {
    return I - 2 * dot(N, I) * N;
}

template <typename VECTYPE> VECTYPE refract(VECTYPE I, VECTYPE N, double eta) {
    double dotNI = dot(N, I);
    double k = 1.0 - eta * eta * (1.0 - dotNI * dotNI);
    if (k < 0.0) {
        return VECTYPE();
    } else {
        return eta * I - (eta * dotNI + sqrt(k)) * N;
    }
}

#endif
