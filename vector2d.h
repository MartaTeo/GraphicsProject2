#pragma once
#include <cmath>

struct Vector2D {
    double x, y;
    Vector2D(): x(0), y(0) {}
    Vector2D(double _x, double _y): x(_x), y(_y) {}

    Vector2D operator+(const Vector2D& o) const { return {x+o.x, y+o.y}; }
    Vector2D operator-(const Vector2D& o) const { return {x-o.x, y-o.y}; }
    Vector2D operator*(double s)           const { return {x*s,   y*s  }; }

    double dot(const Vector2D& o) const { return x*o.x + y*o.y; }
    double norm2()            const { return x*x + y*y; }
    double norm()             const { return std::sqrt(norm2()); }
};
