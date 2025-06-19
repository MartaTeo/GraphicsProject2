
#include "liblbfgs/lib/lbfgs.h"
#include <vector>
#include <random>
#include <iostream>
#include <string>

#include "power_diagram.h"    
#include "svg_writer.h"       
#include "PastelColor.h"      


// L-BFGS evaluator for OT‐dual
struct TransportEvaluator {
    const std::vector<Point>& sites;
    TransportEvaluator(const std::vector<Point>& S): sites(S){}
    static lbfgsfloatval_t evaluate(
        void* instance,
        const lbfgsfloatval_t* w,
        lbfgsfloatval_t* grad,
        const int n,
        const lbfgsfloatval_t)
    {
        auto* self = reinterpret_cast<TransportEvaluator*>(instance);
        std::vector<double> weights(w, w+n);
        auto cells = computePowerDiagram(self->sites, weights);
        double uniform = 1.0/n, fx=0;
        for(int i=0;i<n;++i){
            auto &c=cells[i];
            double A=c.area, cx=c.centroid.x, cy=c.centroid.y, J=c.J;
            double sx=self->sites[i].x, sy=self->sites[i].y;
            double cross = sx*(A*cx)+sy*(A*cy);
            double Ii = J - 2*cross + A*(sx*sx+sy*sy);
            fx    += w[i]*(A-uniform) - Ii;
            grad[i] = A-uniform;
        }
        return fx;
    }
};


int main(int argc, char* argv[]){
    bool gaussian = (argc>1 && std::string(argv[1])=="--gaussian");
    if(gaussian) std::cout<<"[Pastel] Gaussian weights\n";

    constexpr int N=1000, W=700, H=700;
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> U(0,1);

    std::vector<Point> sites(N);
    for(auto &p:sites){ p.x=U(rng); p.y=U(rng); }

    std::vector<double> w(N,0.0);
    if(gaussian){
      for(int i=0;i<N;++i){
        double dx=sites[i].x-0.5, dy=sites[i].y-0.5;
        w[i]=-(dx*dx+dy*dy);
      }
    } else {
      TransportEvaluator eval(sites);
      lbfgs_parameter_t param; lbfgs_parameter_init(&param);
      param.max_iterations=200;
      double fx; int ret = lbfgs(
        N, w.data(), &fx,
        &TransportEvaluator::evaluate,
        nullptr, &eval, &param
      );
      std::cout<<"[Pastel] L-BFGS ret="<<ret<<" fx="<<fx<<"\n";
    }

    auto cells = computePowerDiagram(sites, w);

    save_svg_pastel(cells, "pastel_output.svg", W, H);
    std::cout<<"[Pastel] Wrote pastel_output.svg\n";
    return 0;
}
