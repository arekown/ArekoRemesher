#pragma once

#include "mesh/mesh.h"
#include "mesh/normals.h"
#include "util/reach_queue.h"

#include <cstddef>
#include <vector>

namespace areko::remesher::field {

// Shortest paths over the triangles of a mesh from a set of start triangles:
// a path steps from a triangle to its neighbours across inner edges, and a
// step is as long as the distance between the two triangle centres. First
// in, first out; a triangle that gets a strictly shorter path later is
// visited again.
class CentrePathSearch {
public:
    explicit CentrePathSearch(const mesh::Mesh& mesh)
        : adjacency_(mesh.adjacency()), centre_(mesh::triangleCentres(mesh)),
          search_(static_cast<std::size_t>(mesh.faceCount()))
    {
    }

    // A start triangle: path length 0.
    void start(int t) { search_.start(t); }
    // No triangle is waiting for its visit.
    [[nodiscard]] bool empty() const noexcept { return search_.empty(); }
    // Length of the shortest path to t found so far; FLT_MAX while t is not reached.
    [[nodiscard]] float length(int t) const { return search_.reach(t); }

    // Runs the search to its end. shorter(n, t, e) is called every time
    // triangle n takes a strictly shorter path, which comes from triangle t
    // across edge e.
    template <class Shorter>
    void run(Shorter shorter)
    {
        search_.run([&](int t) {
            const auto at = static_cast<std::size_t>(t);
            for (int c = 3 * t; c < 3 * t + 3; ++c) {
                const int e = adjacency_.cornerEdge(c);
                const int n = e >= 0 ? adjacency_.otherFace(e, t) : -1;
                if (n < 0) {
                    continue;
                }
                if (search_.offer(n, mesh::distance(centre_[at], centre_[static_cast<std::size_t>(n)]) +
                                         search_.reach(t))) {
                    shorter(n, t, e);
                }
            }
        });
    }

private:
    const mesh::Adjacency& adjacency_;
    std::vector<mesh::Vec3f> centre_;
    util::ReachQueue search_;
};

} // namespace areko::remesher::field
