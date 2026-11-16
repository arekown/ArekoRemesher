#include "trace.h"

namespace areko::remesher::extract {

void LoopPolygons::close(FaceLoop& loop, Trace& trace)
{
    for (const int node : loop.nodes) {
        trace.polygons.entries.push_back(trace.nodes[node].id);
    }
    trace.polygons.starts.push_back(static_cast<int>(trace.polygons.entries.size()));
}

} // namespace areko::remesher::extract
