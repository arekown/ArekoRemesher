#include "polyline_locator.h"

#include <cfloat>

namespace areko::remesher::curves {

using namespace areko::remesher::mesh;

LocatedPolyline PolylineLocator::locate(std::span<const Vec3f> points)
{
    LocatedPolyline out;
    const int n = static_cast<int>(points.size());
    if (n == 0) {
        return out;
    }
    const float diag = diagonal(mesh_.bounds());
    const float radius = diag / 40.0f;
    const float snap = diag / 3000.0f;
    const float snap2 = snap * snap;
    const bool closed = points.front() == points.back();
    const int last = closed ? n - 1 : n;
    if (!search_) {
        search_.emplace(mesh_);
    }
    const Topology& topology = mesh_.topology();
    for (int i = 0; i < last; ++i) {
        const ClosestFaceSearch::Hit hit = search_->nearest(points[static_cast<std::size_t>(i)], radius);
        if (hit.face == -1) {
            continue; // no face in reach: the point drops out
        }
        Location location{Location::Face, hit.face};
        Vec3f position = hit.point;
        if (snapToCorners_) {
            float best = FLT_MAX;
            for (int c = 3 * hit.face; c < 3 * hit.face + 3; ++c) {
                const int v = topology.vertex(c);
                const Vec3f d = mesh_.points()[static_cast<std::size_t>(v)] - hit.point;
                const float d2 = dot(d, d);
                if (d2 <= snap2 && d2 < best) {
                    best = d2;
                    location = {Location::Point, v};
                }
            }
            if (location.kind == Location::Point) {
                position = mesh_.points()[static_cast<std::size_t>(location.index)];
            }
        }
        out.points.push_back(position);
        out.locations.push_back(location);
        out.faces.push_back(hit.face);
    }
    if (closed && !out.points.empty()) {
        out.points.push_back(out.points.front());
        out.locations.push_back(out.locations.front());
        out.faces.push_back(out.faces.front());
    }
    return out;
}

} // namespace areko::remesher::curves
