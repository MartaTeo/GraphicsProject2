#pragma once
#include <vector>
#include "vector2d.h"

class Polygon {
public:
    std::vector<Vector2D> vertices;

    Polygon() {}
    explicit Polygon(const std::vector<Vector2D>& v): vertices(v) {}

    double area() const;

    Vector2D centroid() const;

    bool contains(const Vector2D& P) const;

    Polygon clipAgainst(const Polygon& clip) const;

    static Polygon box();
};
