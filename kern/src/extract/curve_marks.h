#pragma once

#include "mesh/bit_array.h"
#include "mesh/csr.h"

#include <vector>

namespace areko::remesher::parts {
struct Part;
} // namespace areko::remesher::parts

namespace areko::remesher::extract {

// The crease curves of a part on its triangle mesh: their points and edges,
// the edges joined into chains, and per edge its chain and place in it. The
// steps after the trace read them where faces are merged and points weighed.
struct CurveMarks {
    mesh::BitArray points;
    mesh::BitArray edges;
    int edgeCount = 0;
    // Chain -> its edges; chain -> its points.
    mesh::Csr chainEdges;
    mesh::Csr chainPoints;
    // Per edge its chain (-1 for none) and its place in the chain, and
    // whether the chain passes it from its first end point to its second.
    std::vector<int> edgeChain;
    std::vector<int> edgePlace;
    mesh::BitArray alongFirstEnd;

    // Marks the active crease curves of `part`: every point of their cut
    // chains, and the edge between every two consecutive points, each once.
    [[nodiscard]] static CurveMarks mark(const parts::Part& part);
};

} // namespace areko::remesher::extract
