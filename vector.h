#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <iostream>

class Vector {
public:
    double data[3] = {0.0, 0.0, 0.0};

    Vector() = default;
    Vector(double x, double y, double z = 0.0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    double x() const { return data[0]; }
    double y() const { return data[1]; }
    double& x() { return data[0]; }
    double& y() { return data[1]; }

    double norm2() const {
        return data[0]*data[0] + data[1]*data[1] + data[2]*data[2];
    }

    double norm() const {
        return std::sqrt(norm2());
    }

    void normalize() {
        double len = norm();
        if (len > 0.0) {
            for (int i = 0; i < 3; ++i) {
                data[i] /= len;
            }
        }
    }

    Vector normalized() const {
        double len = norm();
        return (len > 0.0) ? (*this / len) : *this;
    }

    Vector operator/(double s) const {
        return Vector(data[0] / s, data[1] / s, data[2] / s);
    }

    double operator[](int idx) const { return data[idx]; }
    double& operator[](int idx) { return data[idx]; }

    void print() const {
        std::cout << "Vector("
                  << data[0] << ", "
                  << data[1] << ", "
                  << data[2] << ")";
    }

    int max_arg() const {
        int best = 0;
        for (int i = 1; i < 3; ++i) {
            if (data[i] > data[best]) best = i;
        }
        return best;
    }

    double dot(const Vector& v) const {
        return data[0]*v.data[0] + data[1]*v.data[1] + data[2]*v.data[2];
    }
};

inline Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

inline Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

inline Vector operator-(const Vector& v) {
    return Vector(-v[0], -v[1], -v[2]);
}

inline Vector operator*(double s, const Vector& v) {
    return Vector(v[0] * s, v[1] * s, v[2] * s);
}

inline Vector operator*(const Vector& v, double s) {
    return s * v;
}

inline Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

inline double dot(const Vector& a, const Vector& b) {
    return a.dot(b);
}

inline Vector cross(const Vector& a, const Vector& b) {
    return Vector(
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    );
}

#endif // VECTOR_H
