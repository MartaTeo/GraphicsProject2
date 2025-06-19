#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include "vector2d.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "power_diagram.h"

inline void save_fluid_blobs(const std::vector<Point>& positions, std::string filename, int frameid = 0) {
    const int W = 1000, H = 1000;
    const int R = 30;
    std::vector<unsigned char> image(W * H * 3, 255);

    for (const auto& p : positions) {
        int cx = p.x * W;
        int cy = (1 - p.y) * H;  // flip vertically

        for (int dy = -R; dy <= R; ++dy) {
            for (int dx = -R; dx <= R; ++dx) {
                if (dx*dx + dy*dy > R*R) continue;
                int x = cx + dx;
                int y = cy + dy;
                if (x < 0 || x >= W || y < 0 || y >= H) continue;
                int idx = (y * W + x) * 3;
                image[idx + 0] = 0;
                image[idx + 1] = 0;
                image[idx + 2] = 200;
            }
        }
    }

    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, image.data(), 0);
}
