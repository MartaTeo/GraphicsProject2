#include "polygon.h"
#include <cmath>

static Vector2D intersect(const Vector2D& A,
                          const Vector2D& B,
                          const Vector2D& P,
                          const Vector2D& N)
{
    Vector2D AB = B - A;
    double t = N.dot(P - A) / N.dot(AB);
    return A + AB * t;
}

double Polygon::area() const {
    double A = 0;
    int n = vertices.size();
    for (int i = 0; i < n; ++i) {
        auto& p = vertices[i];
        auto& q = vertices[(i+1)%n];
        A += (p.x * q.y - q.x * p.y);
    }
    return std::abs(A)*0.5;
}

Vector2D Polygon::centroid() const {
    double A = 0, cx = 0, cy = 0;
    int n = vertices.size();
    for (int i = 0; i < n; ++i) {
        auto& p = vertices[i];
        auto& q = vertices[(i+1)%n];
        double cr = p.x*q.y - q.x*p.y;
        A   += cr;
        cx  += (p.x + q.x)*cr;
        cy  += (p.y + q.y)*cr;
    }
    A *= 0.5;
    if (std::abs(A) < 1e-12) return {0,0};
    return { cx/(6*A), cy/(6*A) };
}

bool Polygon::contains(const Vector2D& P) const {
    bool inside = false;
    int n = vertices.size();
    for(int i=0,j=n-1;i<n;j=i++){
        auto& A = vertices[i];
        auto& B = vertices[j];
        if ( ((A.y> P.y) != (B.y> P.y)) &&
             (P.x < (B.x-A.x)*(P.y-A.y)/(B.y-A.y) + A.x) )
            inside = !inside;
    }
    return inside;
}

Polygon Polygon::clipAgainst(const Polygon& clip) const {
    Polygon out = *this;
    int m = clip.vertices.size();
    for(int i=0;i<m;++i){
        Vector2D A = clip.vertices[i];
        Vector2D B = clip.vertices[(i+1)%m];
        Vector2D edge = B - A;
        Vector2D N{-edge.y, edge.x};  

        std::vector<Vector2D> in = out.vertices;
        out.vertices.clear();
        int k = in.size();
        for(int e=0;e<k;++e){
            Vector2D P = in[e], Q = in[(e+1)%k];
            double dP = N.dot(P - A), dQ = N.dot(Q - A);
            if (dP >= 0) {
                out.vertices.push_back(P);
                if (dQ < 0) out.vertices.push_back(intersect(P,Q,A,N));
            } else if (dQ >= 0) {
                out.vertices.push_back(intersect(P,Q,A,N));
                out.vertices.push_back(Q);
            }
        }
    }
    return out;
}

Polygon Polygon::box() {
    return Polygon({{0,0},{1,0},{1,1},{0,1}});
}
