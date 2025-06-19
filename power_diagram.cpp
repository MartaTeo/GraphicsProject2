// power_diagram.cpp
#include "power_diagram.h"
#include <cmath>

static inline double dot(const Point& a, const Point& b){
    return a.x*b.x + a.y*b.y;
}

// Sutherland–Hodgman half‐plane clip
static void sutherlandHodgman(
    const std::vector<Point>& in,
    std::vector<Point>&       out,
    const Point&              n,
    double                    c)
{
    out.clear();
    size_t M = in.size();
    if (M==0) return;
    for (size_t i = 0; i < M; ++i) {
        const Point &P = in[i], &Q = in[(i+1)%M];
        double dP = dot(n,P) - c, dQ = dot(n,Q) - c;
        bool inP = (dP<=0), inQ = (dQ<=0);
        if (inP && inQ) {
            out.push_back(Q);
        } else if (inP && !inQ) {
            double t = dP/(dP-dQ);
            out.push_back({P.x + t*(Q.x-P.x), P.y + t*(Q.y-P.y)});
        } else if (!inP && inQ) {
            double t = dP/(dP-dQ);
            out.push_back({P.x + t*(Q.x-P.x), P.y + t*(Q.y-P.y)});
            out.push_back(Q);
        }
    }
}

std::vector<PowerCell> computePowerDiagram(
    const std::vector<Point>& sites,
    const std::vector<double>& weights,
    double x0, double y0,
    double x1, double y1)
{
    size_t N = sites.size();
    std::vector<PowerCell> cells(N);
    std::vector<Point> poly, tmp;
    poly.reserve(64); tmp.reserve(64);

    for (size_t i = 0; i < N; ++i) {
        poly = {{x0,y0},{x1,y0},{x1,y1},{x0,y1}};
        const auto& si = sites[i]; double wi = weights[i];
        for (size_t j = 0; j < N; ++j) {
            if (j==i) continue;
            const auto& sj = sites[j]; double wj = weights[j];
            Point n{sj.x-si.x, sj.y-si.y};
            double c = ((sj.x*sj.x+sj.y*sj.y)
                      - (si.x*si.x+si.y*si.y)
                      + wj - wi) * 0.5;
            sutherlandHodgman(poly, tmp, n, c);
            poly.swap(tmp);
            if (poly.empty()) break;
        }
        double A2=0,Cx=0,Cy=0,Jn=0;
        size_t M = poly.size();
        for (size_t k=0; k<M; ++k){
            const auto &P=poly[k], &Q=poly[(k+1)%M];
            double cr = P.x*Q.y - Q.x*P.y;
            A2  += cr;
            Cx  += (P.x+Q.x)*cr;
            Cy  += (P.y+Q.y)*cr;
            Jn  += (P.x*P.x + P.x*Q.x + Q.x*Q.x
                  + P.y*P.y + P.y*Q.y + Q.y*Q.y)*cr;
        }
        double A = 0.5*A2;
        Point cent{0,0};
        if (std::abs(A)>1e-16) {
            cent.x = Cx/(6*A); cent.y = Cy/(6*A);
        }
        cells[i] = { poly, A, cent, Jn/12.0 };
    }
    return cells;
}
