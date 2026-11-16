#include "connected_parts.h"

#include <vector>

namespace areko::remesher::mesh {

std::vector<int> ConnectedParts::faceComponents(const Topology& topology, const Adjacency& adjacency,
                                                const BitArray* blocked, int& componentCount)
{
    const int faces = topology.faceCount();
    std::vector<int> component(static_cast<std::size_t>(faces), -1);
    auto open = [blocked](int e) {
        return blocked == nullptr || static_cast<std::size_t>(e) >= blocked->size() || !blocked->test(static_cast<std::size_t>(e));
    };
    int part = 0;
    std::vector<int> queue; // FIFO: read position `head`
    for (int seed = 0; seed < faces; ++seed) {
        if (component[static_cast<std::size_t>(seed)] >= 0) {
            continue;
        }
        component[static_cast<std::size_t>(seed)] = part;
        queue.assign(1, seed);
        for (std::size_t head = 0; head < queue.size(); ++head) {
            const int f = queue[head];
            for (int c = topology.faceBegin(f); c < topology.faceEnd(f); ++c) {
                const int e = adjacency.cornerEdge(c);
                if (e < 0 || !open(e)) {
                    continue;
                }
                for (const int other : adjacency.edgeCorners(e)) {
                    const int g = adjacency.cornerFace(other);
                    if (g >= 0 && component[static_cast<std::size_t>(g)] < 0) {
                        component[static_cast<std::size_t>(g)] = part;
                        queue.push_back(g);
                    }
                }
            }
        }
        ++part;
    }
    componentCount = part;
    return component;
}

ConnectedParts::ConnectedParts(const Topology& topology, const Adjacency& adjacency)
{
    int part = 0;
    faceParts_ = faceComponents(topology, adjacency, nullptr, part);
    partFaces_ = Csr::fromValues(faceParts_, part);

    pointParts_.assign(static_cast<std::size_t>(adjacency.pointCount()), -1);
    for (int c = 0; c < topology.cornerCount(); ++c) {
        const int v = topology.vertex(c);
        const int f = adjacency.cornerFace(c);
        if (v < 0 || f < 0) {
            continue;
        }
        int& slot = pointParts_[static_cast<std::size_t>(v)];
        const int p = faceParts_[static_cast<std::size_t>(f)];
        if (slot == -1) {
            slot = p;
        } else if (slot != p) {
            slot = -2;
        }
    }
}

} // namespace areko::remesher::mesh
