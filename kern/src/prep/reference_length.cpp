#include "reference_length.h"

#include "freeze_marks.h"
#include "mesh/boundary_loops.h"

#include <cmath>
#include <limits>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

constexpr int largeMesh = 1000000;

// The mean over the divisor; a divisor of 0 leaves the mean as it is.
float divided(float mean, float divisor)
{
    return divisor == 0.0f ? mean : mean / divisor;
}

} // namespace

float ReferenceLength::normalRule(const Mesh& mesh, float targetEdgeLength, int quadCount, float* areaCache)
{
    if (targetEdgeLength > 0.0f) {
        return targetEdgeLength;
    }
    float area = areaCache != nullptr ? *areaCache : 0.0f;
    if (!(area > 0.0f)) {
        area = mesh.triangleGeometry().totalArea();
        if (areaCache != nullptr) {
            *areaCache = area;
        }
    }
    return std::sqrt(area / static_cast<float>(quadCount));
}

float ReferenceLength::compute(const Mesh& mesh, float targetEdgeLength, int quadCount, bool freezeBorders,
                               float* areaCache)
{
    if (!freezeBorders) {
        return normalRule(mesh, targetEdgeLength, quadCount, areaCache);
    }
    if (mesh.pointCount() > largeMesh) {
        return normalRule(mesh, targetEdgeLength, clampedQuadCount(quadCount), areaCache);
    }
    const auto* marks = FreezeMarks::find(mesh);
    return borderLength(mesh, marks != nullptr ? &marks->values() : nullptr);
}

int ReferenceLength::clampedQuadCount(int quadCount)
{
    if (quadCount < 3000) {
        return 3000;
    }
    if (quadCount < 8000) {
        return quadCount;
    }
    return quadCount < 15000
               ? truncToInt(((static_cast<float>(quadCount) - 8000.0f) * 4000.0f) / 7000.0f + 8000.0f)
               : 12000;
}

float ReferenceLength::decimationLength(const Mesh& mesh, float reference, bool freezeBorders, int workResolution)
{
    float length = 1.2f * reference;
    if (freezeBorders) {
        const float s = std::sqrt(mesh.triangleGeometry().totalArea() / static_cast<float>(2 * workResolution));
        length = s > length ? s : length;
    }
    return length;
}

float ReferenceLength::borderLength(const Mesh& mesh, const std::vector<FreezeMark>* marks, const BitArray* locked,
                                    float divisor)
{
    auto isLocked = [&](int e) { return locked != nullptr && locked->testBounded(static_cast<std::size_t>(e)); };
    if (marks == nullptr) {
        const BoundaryLoops loops(mesh);
        const auto& points = mesh.points();
        auto at = [&](int v) { return points[static_cast<std::size_t>(v)]; };
        if (mesh.connectedParts().partCount() == 0) {
            return 0.0f;
        }
        const Adjacency& adjacency = mesh.adjacency();
        float lo = 0.0f;
        float hi = std::numeric_limits<float>::max();
        float mean = 0.0f;
        for (int pass = 1; pass <= 2; ++pass) {
            float sum = 0.0f;
            int count = 0;
            for (int r = 0; r < loops.count(); ++r) {
                for (const int e : loops.edges().row(r)) {
                    if (isLocked(e)) {
                        continue;
                    }
                    const auto ends = adjacency.edgeVertices(e);
                    const float l = distance(at(ends[0]), at(ends[1]));
                    if (!(lo <= l && l <= hi)) {
                        continue;
                    }
                    sum += l;
                    ++count;
                }
            }
            if (count == 0) {
                return 0.0f;
            }
            mean = sum / static_cast<float>(count);
            lo = 0.5f * mean;
            hi = 3.0f * mean;
        }
        return divided(mean, divisor);
    }
    float sum = 0.0f;
    int count = 0;
    FreezeMarks::forEachKeptChord(mesh, *marks, [&](const FreezeMarks::KeptChord& chord) {
        // A chord counts unless one of the edges it spans is locked.
        bool spansLocked = false;
        chord.forEachEdge([&](int e) { spansLocked = spansLocked || isLocked(e); });
        if (!spansLocked) {
            sum += chord.length;
            ++count;
        }
    });
    if (count == 0) {
        return 0.0f;
    }
    return divided(sum / static_cast<float>(count), divisor);
}

} // namespace areko::remesher::prep
