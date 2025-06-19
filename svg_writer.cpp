#include "svg_writer.h"
#include <cstdio>
#include <fstream>
#include <tuple>
#include <random>

// std::tuple<float,float,float> randomPastelColor() {
//     static std::mt19937 rng(42);
//     static std::uniform_real_distribution<float> hue(0.0f, 1.0f);
//     float h = hue(rng);
//     float s = 0.5f, v = 0.95f;
//     int i = int(h * 6.0f);
//     float f = h * 6.0f - i;
//     float p = v * (1.0f - s);
//     float q = v * (1.0f - f * s);
//     float t = v * (1.0f - (1.0f - f) * s);
//     float r, g, b;
//     switch (i % 6) {
//         case 0: r = v; g = t; b = p; break;
//         case 1: r = q; g = v; b = p; break;
//         case 2: r = p; g = v; b = t; break;
//         case 3: r = p; g = q; b = v; break;
//         case 4: r = t; g = p; b = v; break;
//         case 5: r = v; g = p; b = q; break;
//     }
//     return {r, g, b};
// }

void savePowerSVG(const std::string& filename,
                  const std::vector<PowerCell>& cells,
                  const std::vector<Point>& sites,
                  const std::vector<Point>& samples,
                  bool show_centroids)
{
    const int W = 512, H = 512;
    FILE* out = std::fopen(filename.c_str(), "w");
    if (!out) return;

    std::fprintf(out,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\">\n",
        W, H);

    for (auto& cell : cells) {
        auto [r, g, b] = randomPastelColor();
        std::fprintf(out, "<polygon points=\"");
        for (auto& p : cell.poly)
            std::fprintf(out, "%f,%f ", p.x * W, H - p.y * H);
        std::fprintf(out, "\" fill=\"rgb(%d,%d,%d)\" stroke=\"black\" stroke-width=\"0.5\"/>\n",
                     int(r * 255), int(g * 255), int(b * 255));
    }

    // for (auto& p : samples)
    //     std::fprintf(out, "<circle cx=\"%f\" cy=\"%f\" r=\"2\" fill=\"red\"/>\n",
    //                  p.x * W, H - p.y * H);

    // for (auto& s : sites)
    //     std::fprintf(out, "<circle cx=\"%f\" cy=\"%f\" r=\"2\" fill=\"black\"/>\n",
    //                  s.x * W, H - s.y * H);

    if (show_centroids) {
        for (const auto& cell : cells) {
            double cx = cell.centroid.x, cy = cell.centroid.y;
            std::fprintf(out, "<circle cx=\"%f\" cy=\"%f\" r=\"2\" fill=\"red\" />\n",
                         cx * W, H - cy * H);
        }
    }

    std::fprintf(out, "</svg>\n");
    std::fclose(out);
}
