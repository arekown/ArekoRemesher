#include "hole_filler.h"

#include "smoother.h"

#include "mesh/boundary_loops.h"
#include "mesh/modification.h"

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

int HoleFiller::fill(Mesh& mesh, const Options& options)
{
    const Adjacency& adjacency = mesh.adjacency();
    const auto& points = mesh.points();
    const BoundaryLoops loops(mesh);
    std::vector<int> faceStart{0};
    std::vector<int> cornerVertex;
    std::vector<int> origin;
    for (int r = 0; r < loops.count(); ++r) {
        const auto loop = loops.points().row(r);
        const auto n = static_cast<int>(loop.size());
        if (n > options.maxPoints || n < 3) {
            continue;
        }
        auto at = [&](int i) { return points[static_cast<std::size_t>(loop[static_cast<std::size_t>((i + n) % n)])]; };
        float perimeter = 0.0f;
        float area = 0.0f;
        for (int i = 0; i < n; ++i) {
            perimeter += distance(at(i), at(i + 1));
            area += length(cross(at(i) - at(i - 1), at(i + 1) - at(i)));
        }
        area = 0.25f * area;
        const bool small = (options.maxPerimeter > 0.0f && perimeter < options.maxPerimeter) ||
                           (options.maxArea > 0.0f && area < options.maxArea);
        if (!small) {
            continue;
        }
        cornerVertex.insert(cornerVertex.end(), loop.begin(), loop.end());
        faceStart.push_back(static_cast<int>(cornerVertex.size()));
        origin.push_back(adjacency.cornerFace(adjacency.pointCorners(loop[0])[0]));
    }
    const int closed = static_cast<int>(origin.size());
    if (closed == 0) {
        return 0;
    }

    // The new faces refer to the points of the mesh and bring none of their own.
    Mesh holes;
    std::vector<int> newCorners = cornerVertex;
    holes.editTopology().swapPolygons(faceStart, newCorners);
    Modification change;
    change.appendFaces = &holes;
    change.faces.newToOld = std::move(origin);
    change.points.state = 1;
    MeshModifier::apply(mesh, change);

    if (options.smoothIterations > 0) {
        Smoother::Options smoothing;
        smoothing.iterations = options.smoothIterations;
        smoothing.uniformWeight = options.smoothWeight;
        std::vector<Vec3f> smoothed = Smoother(mesh, smoothing).applyTo(cornerVertex, {});
        mesh.swapPoints(smoothed);
    }
    return closed;
}

} // namespace areko::remesher::prep
