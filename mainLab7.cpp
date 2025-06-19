#include "liblbfgs/lib/lbfgs.h"
#include <vector>
#include <random>
#include <iostream>
#include <string>
#include "power_diagram.h"
#include "svg_writer.h"

int main(int argc, char* argv[]) {
    bool gaussian = (argc > 1 && std::string(argv[1]) == "--gaussian");
    if (gaussian) std::cout << "[Lab7] Gaussian mode\n";

    const int N = 1000;
    const int K = 7;
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> U(0.0, 1.0);

    std::vector<Point> sites(N);
    for (auto& p : sites) { p.x = U(rng); p.y = U(rng); }

    std::vector<double> w(N, 0.0);
    if (gaussian) {
        for (int i = 0; i < N; ++i) {
            double dx = sites[i].x - 0.5, dy = sites[i].y - 0.5;
            w[i] = -(dx * dx + dy * dy);
        }
    } else {
        struct TransportEvaluator {
            const std::vector<Point>& sites;
            TransportEvaluator(const std::vector<Point>& S) : sites(S) {}
            static lbfgsfloatval_t evaluate(void* instance, const lbfgsfloatval_t* w, lbfgsfloatval_t* grad, const int n, const lbfgsfloatval_t) {
                auto* self = reinterpret_cast<TransportEvaluator*>(instance);
                std::vector<double> weights(w, w + n);
                auto cells = computePowerDiagram(self->sites, weights);
                double fx = 0.0, uniform = 1.0 / n;
                for (int i = 0; i < n; ++i) {
                    auto& c = cells[i];
                    double A = c.area, cx = c.centroid.x, cy = c.centroid.y, J = c.J;
                    double sx = self->sites[i].x, sy = self->sites[i].y;
                    double cross = sx * (A * cx) + sy * (A * cy);
                    double Ii = J - 2 * cross + A * (sx * sx + sy * sy);
                    fx += w[i] * (A - uniform) - Ii;
                    grad[i] = A - uniform;
                }
                return fx;
            }
        };
        TransportEvaluator eval(sites);
        lbfgs_parameter_t param; lbfgs_parameter_init(&param);
        param.max_iterations = 200;
        double fx;
        int ret = lbfgs(N, w.data(), &fx, &TransportEvaluator::evaluate, nullptr, &eval, &param);
        std::cout << "[Lab7] L-BFGS code=" << ret << " fx=" << fx << "\n";
    }

    auto cells = computePowerDiagram(sites, w);
    std::vector<Point> samples;
    samples.reserve(N * K);
    for (const auto& cell : cells) {
        double minx = 1, miny = 1, maxx = 0, maxy = 0;
        for (auto& p : cell.poly) {
            minx = std::min(minx, p.x);
            miny = std::min(miny, p.y);
            maxx = std::max(maxx, p.x);
            maxy = std::max(maxy, p.y);
        }
        int cnt = 0;
        while (cnt < K) {
            double x = U(rng) * (maxx - minx) + minx;
            double y = U(rng) * (maxy - miny) + miny;
            bool inside = false;
            int N = cell.poly.size();
            for (int i = 0, j = N - 1; i < N; j = i++) {
                const auto& pi = cell.poly[i], &pj = cell.poly[j];
                if (((pi.y > y) != (pj.y > y)) &&
                    (x < (pj.x - pi.x) * (y - pi.y) / (pj.y - pi.y) + pi.x))
                    inside = !inside;
            }
            if (inside) {
                samples.push_back({x, y});
                cnt++;
            }
        }
    }

    savePowerSVG("lab7_with_centroids.svg", cells, sites, samples, true);
}
