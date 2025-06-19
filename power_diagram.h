// power_diagram.h
#ifndef POWER_DIAGRAM_H
#define POWER_DIAGRAM_H

#include <vector>


struct Point { double x, y; };

struct PowerCell {
    std::vector<Point> poly;
    double             area;
    Point              centroid;
    double             J;
};


std::vector<PowerCell> computePowerDiagram(
    const std::vector<Point>& sites,
    const std::vector<double>& weights,
    double x0 = 0.0, double y0 = 0.0,
    double x1 = 1.0, double y1 = 1.0
);

#endif // POWER_DIAGRAM_H
