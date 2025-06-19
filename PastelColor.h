// PastelColor.h
#pragma once
#include <tuple>
#include <random>


inline std::tuple<float,float,float> HSLtoRGB(float h, float s, float l) {
    auto hue2rgb = [&](float p, float q, float t){
        if(t<0) t+=1; if(t>1) t-=1;
        if(t<1/6.f) return p + (q-p)*6*t;
        if(t<1/2.f) return q;
        if(t<2/3.f) return p + (q-p)*(2/3.f-t)*6;
        return p;
    };
    float r,g,b;
    if(s==0) r=g=b=l;
    else {
        float q = l<.5f ? l*(1+s) : l+s-l*s;
        float p = 2*l-q;
        r = hue2rgb(p,q,h+1/3.f);
        g = hue2rgb(p,q,h);
        b = hue2rgb(p,q,h-1/3.f);
    }
    return {r,g,b};
}

inline std::tuple<float,float,float> randomPastelColor() {
    static std::mt19937 rng{std::random_device{}()};
    std::uniform_real_distribution<float> U(0,1);
    float h = U(rng);
    float s = 0.3f + 0.2f*U(rng);
    float l = 0.75f + 0.15f*U(rng);
    return HSLtoRGB(h,s,l);
}
