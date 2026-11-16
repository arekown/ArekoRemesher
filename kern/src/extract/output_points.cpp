#include "output_points.h"

namespace areko::remesher::extract {

void OutputPoints::resize(std::size_t count, PointSource newSource, float newWeight)
{
    onBorderEdge.resize(count, false);
    reachedBorder.resize(count, false);
    onCurve.resize(count, false);
    source.resize(count, newSource);
    weight.resize(count, newWeight);
    triangle.resize(count, -1);
    for (std::vector<int>& list : neighbours) {
        list.resize(count, -1);
    }
}

} // namespace areko::remesher::extract
