// svg_writer.h
#ifndef SVG_WRITER_H
#define SVG_WRITER_H

#include <vector>
#include <string>
#include "PastelColor.h"
#include "power_diagram.h"

/// Draw only the polygon outlines + red sample dots.
/// @param cells    power‐diagram cells (for outlines)
/// @param samples  sample points in [0,1]^2 (red dots)
/// @param filename output SVG file
/// @param width    canvas width (px)
/// @param height   canvas height (px)
/// @param fillcol  polygon fill color, e.g. "none"
void save_svg_with_samples(
    const std::vector<PowerCell>& cells,
    const std::vector<Point>&     samples,
    const std::string&            filename,
    int                            width,
    int                            height,
    const std::string&            fillcol = "none"
);

void save_svg_pastel(
    const std::vector<PowerCell>& cells,
    const std::string &filename,
    int W, int H
);

void savePowerSVG(const std::string& filename,
                  const std::vector<PowerCell>& cells,
                  const std::vector<Point>& sites,
                  const std::vector<Point>& samples,
                  bool show_centroids);

#endif // SVG_WRITER_H
