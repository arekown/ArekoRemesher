#include "curve_distribution.h"

#include <optional>

namespace areko::remesher::curves {

using namespace areko::remesher::mesh;

LocatedPolyline CurveDistribution::locate(const GuideCurve& curve, PolylineLocator& locator)
{
    const int n = static_cast<int>(curve.points.size());
    if (!curve.locations.empty() && curve.pointLocationCount >= 1 && curve.pointLocationCount == n &&
        static_cast<int>(curve.locations.size()) == n) {
        LocatedPolyline out;
        out.points = curve.points;
        out.locations = curve.locations;
        return out;
    }
    return locator.locate(curve.points);
}

std::vector<GuideCurve> CurveDistribution::fromLines(std::span<const prep::NetLine> lines)
{
    std::vector<GuideCurve> out;
    out.reserve(lines.size());
    for (const prep::NetLine& line : lines) {
        GuideCurve curve;
        curve.index = -2;
        curve.points = line.points;
        curve.kind = line.type;
        curve.lineFlags = line.flags;
        for (const int v : line.vertices) {
            curve.locations.push_back({Location::Point, v});
        }
        if (line.flags == 1 || line.flags == 2) {
            curve.pointLocationCount = static_cast<int>(curve.locations.size());
        }
        out.push_back(std::move(curve));
    }
    return out;
}

void CurveDistribution::run(const Mesh& prepared, std::vector<GuideCurve> hostCurves,
                            std::span<const prep::NetLine> lines, parts::PartBuilder::Split& split)
{
    std::vector<GuideCurve> candidates = std::move(hostCurves);
    for (GuideCurve& curve : fromLines(lines)) {
        candidates.push_back(std::move(curve));
    }
    if (split.parts.empty()) {
        return;
    }
    PolylineLocator locator(prepared, true);
    const Adjacency& adjacency = prepared.adjacency();
    for (const GuideCurve& curve : candidates) {
        if (!curve.active) {
            continue;
        }
        const LocatedPolyline located = locate(curve, locator);
        int currentPart = -1;
        std::optional<GuideCurve> run;
        auto finish = [&]() {
            // Every run with at least two points stays with its part.
            if (run && currentPart >= 0 && run->source.points.size() >= 2) {
                split.parts[static_cast<std::size_t>(currentPart)].curves.push_back(std::move(*run));
            }
            run.reset();
        };
        for (std::size_t i = 0; i < located.locations.size(); ++i) {
            const Location& location = located.locations[i];
            int f = -1;
            if (location.kind == Location::Face) {
                f = location.index;
            } else if (location.kind == Location::Point) {
                const auto corners = adjacency.pointCorners(location.index);
                f = corners.empty() ? -1 : adjacency.cornerFace(corners.front());
            } else {
                continue;
            }
            if (f < 0 || f >= prepared.faceCount()) {
                continue;
            }
            const int part = split.faceParts[static_cast<std::size_t>(f)];
            const int place = split.facePlaces[static_cast<std::size_t>(f)];
            if (part < 0 || place < 0) {
                continue;
            }
            const Mesh& partMesh = split.parts[static_cast<std::size_t>(part)].mesh;
            Location local{Location::Face, place};
            if (location.kind == Location::Point) {
                const Vec3f p = prepared.points()[static_cast<std::size_t>(location.index)];
                for (int c = 3 * place; c < 3 * place + 3; ++c) {
                    const int v = partMesh.topology().vertex(c);
                    if (partMesh.points()[static_cast<std::size_t>(v)] == p) {
                        local = {Location::Point, v};
                        break;
                    }
                }
            }
            if (part != currentPart) {
                finish();
                currentPart = part;
                run.emplace();
                run->index = curve.index == -2 ? -1 : curve.index;
                run->points = curve.points;
                run->kind = curve.kind;
                run->active = curve.active;
                run->strength = curve.strength;
                run->width = curve.width;
                run->lineFlags = curve.lineFlags;
                run->locations = curve.locations;
                run->pointLocationCount = curve.pointLocationCount;
                run->sourcePrepared = true;
            }
            run->source.points.push_back(located.points[i]);
            run->source.locations.push_back(local);
            run->source.faces.push_back(place);
        }
        finish();
    }
}

} // namespace areko::remesher::curves
