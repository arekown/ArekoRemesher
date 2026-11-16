#include "degenerate_cleanup.h"

#include "bit_array.h"
#include "face_rebuild.h"
#include "modification.h"

#include <optional>
#include <vector>

namespace areko::remesher::mesh {

DegenerateCleanup::Result DegenerateCleanup::run(Mesh& mesh, const Switches& switches)
{
    Result result;
    const AttributeMaps& maps = mesh.maps();
    const bool faceMap = maps.hasKind(ElementKind::Face);
    // The origin lists of the plain way; FaceRebuild keeps its own.
    const bool plain = !switches.keepMaps;
    const bool faceOrigin = plain && switches.faceOrigin && faceMap;
    const bool halfEdgeOrigin =
        plain && switches.cornerOrigin && (maps.hasKind(ElementKind::Edge) || maps.hasKind(ElementKind::HalfEdge));
    const bool cornerOrigin = plain && (mesh.customNormals().has_value() || halfEdgeOrigin);

    const Topology& topology = mesh.topology();
    const int faceCount = topology.faceCount();
    const int oldCorners = topology.cornerCount();
    std::optional<FaceRebuild> rebuild;
    if (switches.keepMaps) {
        rebuild.emplace(mesh, topology.cornerVertex());
    }
    std::vector<int> faceStart{0};
    std::vector<int> corners;
    corners.reserve(static_cast<std::size_t>(oldCorners));
    std::vector<int> kept; // old corners of the face that stay
    BitArray dropped;
    if (cornerOrigin) {
        dropped.resize(static_cast<std::size_t>(oldCorners));
    }
    std::optional<std::vector<int>> faceSource; // laid out at the first removed face
    auto remove = [&](int f) {
        if (faceOrigin && !faceSource) {
            faceSource.emplace();
            for (int g = 0; g < f; ++g) {
                faceSource->push_back(g);
            }
        }
        if (cornerOrigin) {
            for (int c = topology.faceBegin(f); c < topology.faceEnd(f); ++c) {
                dropped.set(static_cast<std::size_t>(c));
            }
        }
        ++result.removedFaces;
    };

    for (int f = 0; f < faceCount; ++f) {
        const int begin = topology.faceBegin(f);
        const int n = topology.faceSize(f);
        auto point = [&](int k) { return topology.vertex(begin + k); };
        kept.clear();
        if (n < 3 || (n == 4 && (point(0) == point(2) || point(1) == point(3)))) {
            remove(f);
        } else {
            int last = -1;
            for (int k = n - 1; k >= 1 && last == -1; --k) {
                last = point(k);
            }
            if (point(0) == -1 || point(0) == last) {
                if (cornerOrigin) {
                    dropped.set(static_cast<std::size_t>(begin));
                }
            } else {
                kept.push_back(begin);
            }
            int previous = point(0);
            for (int k = 1; k < n; ++k) {
                const int v = point(k);
                if (v == previous || v == -1) {
                    if (cornerOrigin) {
                        dropped.set(static_cast<std::size_t>(begin + k));
                    }
                } else {
                    kept.push_back(begin + k);
                }
                previous = v;
            }
            if (kept.size() < 3) {
                kept.clear();
                remove(f);
            } else {
                for (const int c : kept) {
                    corners.push_back(topology.vertex(c));
                }
                faceStart.push_back(static_cast<int>(corners.size()));
                if (faceSource) {
                    faceSource->push_back(f);
                }
            }
        }
        if (rebuild) {
            rebuild->keep(f, kept);
        }
    }

    if (rebuild) {
        if (rebuild->reduced()) {
            // With a face map the face mode stays; otherwise the faces are set as polygons.
            rebuild->apply(mesh, faceMap ? topology.mode() : FaceMode::Polygons);
        }
        return result;
    }
    if (result.removedFaces == 0 && static_cast<int>(corners.size()) == oldCorners) {
        return result;
    }
    std::vector<int> cornerSource;
    if (cornerOrigin) {
        for (int c = 0; c < oldCorners; ++c) {
            if (!dropped.test(static_cast<std::size_t>(c))) {
                cornerSource.push_back(c);
            }
        }
        cornerSource.resize(corners.size(), 0);
    }
    Modification change;
    change.faceStart = std::move(faceStart);
    change.cornerVertex = std::move(corners);
    if (!faceOrigin) {
        // Set directly: polygons, nothing renumbered.
        change.faceMode = FaceMode::Polygons;
        change.points.state = 1;
    } else {
        if (faceSource) {
            change.faces.newToOld = std::move(faceSource);
        } else {
            change.faces.state = 1;
        }
        if (halfEdgeOrigin) {
            change.halfEdges.newToOld = cornerSource;
        }
        if (cornerOrigin) {
            change.corners.newToOld = std::move(cornerSource);
        }
    }
    MeshModifier::apply(mesh, change);
    return result;
}

} // namespace areko::remesher::mesh
