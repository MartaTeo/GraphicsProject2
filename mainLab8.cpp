#include "power_diagram.h"
#include "save_frame.h"
#include "liblbfgs/lib/lbfgs.h"
#include <vector>
#include <random>
#include <iostream>
#include <sstream>
#include <algorithm>

struct Site {
    Point position;
    Vector2D velocity;
    double weight;
    double target_mass;
    bool is_fluid;
};

int main() {
    const int N = 300;
    const double FLUID_FRAC = 0.25;
    const double DT = 0.02;
    const int STEPS = 200;
    const double GRAVITY = -0.01;
    const double VELOCITY_GAIN = 2.0;
    const double DAMPING = 0.98;

    std::vector<Site> sites(N);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> Ux(0.3, 0.7);
    std::uniform_real_distribution<double> Uy(0.7, 0.95);

    for (int i = 0; i < N; ++i) {
        sites[i].position = {Ux(rng), Uy(rng)};
        sites[i].velocity = {0.0, 0.0};
        sites[i].weight = 0.0;
        sites[i].is_fluid = (i < N * FLUID_FRAC);
        sites[i].target_mass = sites[i].is_fluid ? 1.0 / (N * FLUID_FRAC) : 0.0;
    }

    struct FluidEvaluator {
        const std::vector<Point>& positions;
        const std::vector<double>& masses;
        FluidEvaluator(const std::vector<Point>& pos, const std::vector<double>& m)
            : positions(pos), masses(m) {}

        static lbfgsfloatval_t evaluate(void* instance, const lbfgsfloatval_t* w,
                                        lbfgsfloatval_t* g, const int n, const lbfgsfloatval_t) {
            auto* self = reinterpret_cast<FluidEvaluator*>(instance);
            std::vector<double> weights(w, w + n);
            auto cells = computePowerDiagram(self->positions, weights);
            double fx = 0;
            for (int i = 0; i < n; ++i) {
                double A = cells[i].area;
                fx += w[i] * (A - self->masses[i]) - cells[i].J;
                g[i] = A - self->masses[i];
            }
            return fx;
        }
    };

    for (int frame = 0; frame < STEPS; ++frame) {
        std::vector<Point> positions(N);
        std::vector<double> weights(N), masses(N);
        for (int i = 0; i < N; ++i) {
            positions[i] = sites[i].position;
            weights[i] = sites[i].weight;
            masses[i] = sites[i].target_mass;
        }

        FluidEvaluator eval(positions, masses);
        lbfgs_parameter_t param;
        lbfgs_parameter_init(&param);
        param.max_iterations = 100;
        double fx;
        lbfgs(N, weights.data(), &fx, FluidEvaluator::evaluate, nullptr, &eval, &param);

        auto cells = computePowerDiagram(positions, weights);

        for (int i = 0; i < N; ++i) {
            sites[i].weight = weights[i];
            if (sites[i].is_fluid) {
                Point c = cells[i].centroid;
                Vector2D disp = {c.x - sites[i].position.x, c.y - sites[i].position.y};
                sites[i].velocity = sites[i].velocity + disp * VELOCITY_GAIN;
                sites[i].velocity.y += GRAVITY;

                if (sites[i].position.y < 0.05)
                    sites[i].velocity = sites[i].velocity * 0.9;
                else
                    sites[i].velocity = sites[i].velocity * DAMPING;

                sites[i].position.x += DT * sites[i].velocity.x;
                sites[i].position.y += DT * sites[i].velocity.y;

                // Clamp Y and cancel bounce
                if (sites[i].position.y < 0.01) {
                    sites[i].position.y = 0.01;
                    if (sites[i].velocity.y < 0)
                        sites[i].velocity.y = 0;
                }
                if (sites[i].position.y > 0.99)
                    sites[i].position.y = 0.99;

                // Clamp X
                sites[i].position.x = std::clamp(sites[i].position.x, 0.01, 0.99);
            }
        }

        std::vector<Point> fluid_sites;
        for (int i = 0; i < N; ++i)
            if (sites[i].is_fluid)
                fluid_sites.push_back(sites[i].position);

        save_fluid_blobs(fluid_sites, "frame", frame);
    }
}
