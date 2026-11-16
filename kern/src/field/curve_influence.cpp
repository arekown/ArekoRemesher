#include "curve_influence.h"

#include "curves/polyline_locator.h"
#include "util/index_queue.h"

#include <cfloat>
#include <utility>

namespace areko::remesher::field {

using namespace areko::remesher::mesh;
using curves::GuideCurve;
using curves::isBand;

float CurveInfluence::length(const std::vector<Vec3f>& points)
{
    float sum = 0.0f;
    for (std::size_t i = 0; i + 1 < points.size(); ++i) {
        sum = mesh::distance(points[i], points[i + 1]) + sum;
    }
    return sum;
}

void CurveInfluence::assignWidths(std::vector<GuideCurve>& curves, float maxRadius)
{
    for (GuideCurve& curve : curves) {
        if (!curve.active || !(curve.kind < ARM_CURVE_OFF) || curve.width >= 0.0f) {
            continue;
        }
        const float quarter = 0.25f * length(isBand(curve.kind) ? curve.points : curve.cutPoints);
        const float b = quarter < maxRadius ? quarter : maxRadius;
        curve.width = isBand(curve.kind) ? b : b * 0.75f;
    }
}

void CurveInfluence::prepareBands(std::vector<GuideCurve>& curves, const Mesh& mesh)
{
    curves::PolylineLocator locator(mesh, false);
    for (GuideCurve& curve : curves) {
        if (!curve.active || !isBand(curve.kind) || curve.bandPrepared) {
            continue;
        }
        curve.bandPrepared = true;
        curves::LocatedPolyline laid = locator.locate(curve.source.points);
        curve.cutPoints = std::move(laid.points);
        curve.cutFaces = std::move(laid.faces);
    }
}

CurveDistance CurveInfluence::distance(const Mesh& mesh, const std::vector<GuideCurve>& curves)
{
    const Adjacency& adjacency = mesh.adjacency();
    const std::vector<Vec3f>& points = mesh.points();
    const std::size_t pointCount = points.size();
    CurveDistance out;
    out.pointDistance.assign(pointCount, FLT_MAX);
    out.pointSeed.assign(pointCount, -1);
    std::vector<float>& dist = out.pointDistance;
    std::vector<int>& seed = out.pointSeed;

    util::IndexQueue queue(pointCount);
    for (std::size_t i = 0; i < curves.size(); ++i) {
        const GuideCurve& curve = curves[i];
        if (!curve.active || !curves::counts(curve.kind)) {
            continue;
        }
        auto plant = [&](int v, int j) {
            const std::size_t vi = static_cast<std::size_t>(v);
            if (dist[vi] == 0.0f) {
                return; // the first seed wins
            }
            queue.pushBack(v);
            dist[vi] = 0.0f;
            seed[vi] = static_cast<int>(out.seedCurve.size());
            out.seedCurve.push_back(static_cast<int>(i));
            out.seedPoint.push_back(j);
        };
        if (isBand(curve.kind)) {
            for (std::size_t j = 0; j < curve.cutFaces.size(); ++j) {
                const int f = curve.cutFaces[j];
                if (f == -1) {
                    continue;
                }
                for (int c = 3 * f; c < 3 * f + 3; ++c) {
                    plant(adjacency.vertex(c), static_cast<int>(j));
                }
            }
        } else {
            for (std::size_t j = 0; j < curve.cutVertices.size(); ++j) {
                if (curve.cutVertices[j] != -1) {
                    plant(curve.cutVertices[j], static_cast<int>(j));
                }
            }
        }
    }
    while (!queue.empty()) {
        const int v = queue.pop();
        const std::size_t vi = static_cast<std::size_t>(v);
        const int s = seed[vi];
        const GuideCurve& curve = curves[static_cast<std::size_t>(out.seedCurve[static_cast<std::size_t>(s)])];
        const Vec3f from = curve.cutPoints[static_cast<std::size_t>(out.seedPoint[static_cast<std::size_t>(s)])];
        const float reach = curve.width + curve.width;
        for (const int w : adjacency.pointNeighbours(v)) {
            if (w == -1) {
                continue;
            }
            const std::size_t wi = static_cast<std::size_t>(w);
            const float a = mesh::distance(points[wi], from);
            const float grown = 1.005f * dist[vi];
            const float c = grown > a ? grown : a;
            if (c < dist[wi] && c < reach) {
                queue.pushBack(w);
                dist[wi] = c;
                seed[wi] = s;
            }
        }
    }

    const Topology& topology = mesh.topology();
    out.faceDistance.assign(static_cast<std::size_t>(mesh.faceCount()), FLT_MAX);
    for (int f = 0; f < mesh.faceCount(); ++f) {
        float sum = 0.0f;
        bool seeded = true;
        for (int c = 3 * f; c < 3 * f + 3 && seeded; ++c) {
            const std::size_t v = static_cast<std::size_t>(topology.vertex(c));
            seeded = seed[v] >= 0;
            if (seeded) {
                const float width = curves[static_cast<std::size_t>(out.seedCurve[static_cast<std::size_t>(seed[v])])].width;
                sum = c == 3 * f ? dist[v] / width : sum + dist[v] / width;
            }
        }
        if (seeded) {
            out.faceDistance[static_cast<std::size_t>(f)] = sum / 3.0f;
        }
    }
    return out;
}

std::vector<Vec3f> CurveInfluence::tangents(const std::vector<Vec3f>& points)
{
    const std::size_t n = points.size();
    std::vector<Vec3f> out(n);
    if (n < 2) {
        return out;
    }
    out[0] = points[1] - points[0];
    for (std::size_t i = 1; i + 1 < n; ++i) {
        out[i] = points[i + 1] - points[i - 1];
    }
    out[n - 1] = points[n - 1] - points[n - 2];
    for (Vec3f& t : out) {
        normalizeInPlace(t);
    }
    return out;
}

} // namespace areko::remesher::field
