#include "modification.h"

#include "edge_chains.h"
#include "mesh.h"
#include "sharp_edges.h"

namespace areko::remesher::mesh {

namespace {

struct Merged {
    std::vector<int> faceStart{0};
    std::vector<int> cornerVertex;
    std::vector<std::vector<int>> uvCorners; // one list per UV set; empty = clear the set
    bool touchesUv = false;
    std::vector<int> faceOrigin;
    std::vector<int> cornerOrigin;
};

void appendCorners(const Topology& topology, int f, std::vector<int>& to)
{
    const auto& corners = topology.cornerVertex();
    to.insert(to.end(), corners.begin() + topology.faceBegin(f), corners.begin() + topology.faceEnd(f));
}

// Kept faces of `mesh` in their order, then the faces of `added` in theirs.
Merged merge(const Mesh& mesh, const Mesh* added, const BitArray* deleted,
             const std::optional<std::vector<int>>& addedOrigin, bool cornersNeeded)
{
    Merged out;
    const Topology& topology = mesh.topology();
    const bool polygons = topology.mode() == FaceMode::Polygons;
    const int faces = topology.faceCount();
    const int addedFaces = added != nullptr ? added->faceCount() : 0;
    const int addedCorners = added != nullptr ? added->cornerCount() : 0;
    const std::size_t deletedCount = deleted != nullptr ? deleted->count() : 0;
    auto kept = [deleted](int f) {
        return deleted == nullptr || static_cast<std::size_t>(f) >= deleted->size() || !deleted->test(static_cast<std::size_t>(f));
    };

    if (added != nullptr && faces == 0 && deletedCount == 0) {
        out.cornerVertex = added->topology().cornerVertex();
        if (polygons) {
            for (int f = 0; f < addedFaces; ++f) {
                out.faceStart.push_back(added->topology().faceEnd(f));
            }
        }
        return out;
    }

    const auto& uvSets = mesh.uvSets();
    out.touchesUv = true;
    out.uvCorners.assign(uvSets.size(), {});

    // Face origin (new -> old).
    if (addedOrigin && static_cast<int>(addedOrigin->size()) == addedFaces) {
        for (int f = 0; f < faces; ++f) {
            if (deletedCount == 0 || kept(f)) {
                out.faceOrigin.push_back(f);
            }
        }
        out.faceOrigin.insert(out.faceOrigin.end(), addedOrigin->begin(), addedOrigin->end());
    } else if (deletedCount > 0) {
        for (int f = 0; f < faces; ++f) {
            if (kept(f)) {
                out.faceOrigin.push_back(f);
            }
        }
        out.faceOrigin.insert(out.faceOrigin.end(), static_cast<std::size_t>(addedFaces), -1);
    }

    if (deletedCount == 0) {
        if (polygons) {
            out.faceStart = topology.faceStart();
        }
        out.cornerVertex = topology.cornerVertex();
        for (std::size_t k = 0; k < uvSets.size(); ++k) {
            out.uvCorners[k] = uvSets[k].cornerUv;
        }
    } else {
        for (int f = 0; f < faces; ++f) {
            if (!kept(f)) {
                continue;
            }
            appendCorners(topology, f, out.cornerVertex);
            if (polygons) {
                out.faceStart.push_back(static_cast<int>(out.cornerVertex.size()));
            }
            for (std::size_t k = 0; k < uvSets.size(); ++k) {
                const auto& uv = uvSets[k].cornerUv;
                if (!uv.empty()) {
                    out.uvCorners[k].insert(out.uvCorners[k].end(), uv.begin() + topology.faceBegin(f),
                                            uv.begin() + topology.faceEnd(f));
                }
            }
        }
    }
    if (cornersNeeded) {
        for (int f = 0; f < faces; ++f) {
            if (deletedCount == 0 || kept(f)) {
                for (int c = topology.faceBegin(f); c < topology.faceEnd(f); ++c) {
                    out.cornerOrigin.push_back(c);
                }
            }
        }
        // One entry per appended corner.
        out.cornerOrigin.insert(out.cornerOrigin.end(), static_cast<std::size_t>(addedCorners), -1);
    }

    if (added != nullptr) {
        const Topology& more = added->topology();
        for (int f = 0; f < addedFaces; ++f) {
            appendCorners(more, f, out.cornerVertex);
            if (polygons) {
                out.faceStart.push_back(static_cast<int>(out.cornerVertex.size()));
            }
        }
        for (std::size_t k = 0; k < uvSets.size(); ++k) {
            if (k >= added->uvSets().size()) {
                out.uvCorners[k].clear();
                continue;
            }
            const auto& uv = added->uvSets()[k].cornerUv;
            out.uvCorners[k].insert(out.uvCorners[k].end(), uv.begin(), uv.end());
            if (out.uvCorners[k].size() != out.cornerVertex.size()) {
                out.uvCorners[k].clear();
            }
        }
    }
    return out;
}

void followNormals(Mesh& mesh, const Modification& change)
{
    auto& normals = mesh.customNormals();
    if (!normals) {
        return;
    }
    if (change.faces.unchanged() && change.corners.unchanged()) {
        return;
    }
    const std::optional<std::vector<int>>& origin =
        change.corners.newToOld ? change.corners.newToOld : change.halfEdges.newToOld;
    if (!origin || origin->empty() || normals->cornerNormal.empty()) {
        normals.reset();
        return;
    }
    std::vector<int> renumbered(origin->size(), -1);
    bool missing = false;
    for (std::size_t c = 0; c < origin->size(); ++c) {
        const int o = (*origin)[c];
        if (o >= 0 && static_cast<std::size_t>(o) < normals->cornerNormal.size()) {
            renumbered[c] = normals->cornerNormal[static_cast<std::size_t>(o)];
        }
        missing = missing || renumbered[c] < 0;
    }
    normals->cornerNormal = std::move(renumbered);
    normals->state = missing ? 2 : 1;
}

} // namespace

bool MeshModifier::apply(Mesh& mesh, Modification& change)
{
    const Mesh* added = change.appendFaces;
    if (added != nullptr && added->pointCount() != 0) {
        return false;
    }
    const bool newPoints = change.replacePoints.has_value() || !change.appendPoints.empty();
    const bool deletion = change.deleteFaces.has_value() && change.deleteFaces->any();
    if (!deletion && added != nullptr && added->faceCount() == 0 && !newPoints) {
        return false;
    }
    const int oldPointCount = mesh.pointCount();
    const int oldCornerCount = mesh.cornerCount();
    if (!newPoints && !change.cornerVertex && !change.faceStart && added == nullptr && deletion) {
        change.points.state = 1; // only faces are removed: point maps stay
    }
    if (added != nullptr && added->topology().mode() != mesh.topology().mode()) {
        mesh.toPolygonMode();
    }

    // Sharp edges across the change.
    bool halfEdgeForm = false;
    std::optional<Csr> chains;
    if (change.faceStart || change.cornerVertex) {
        if (!change.halfEdges.oldToNew) {
            if (change.halfEdges.unchanged()) {
                halfEdgeForm = SharpEdges::toHalfEdgeForm(mesh);
            } else if (!change.halfEdges.newToOld) {
                if (const BoolMap* edges = SharpEdges::edgeForm(mesh)) {
                    chains = EdgeChains::build(mesh.adjacency(), edges->bits());
                }
            } else {
                // The half edges name their origin: the sharp edges follow it.
                halfEdgeForm = SharpEdges::toHalfEdgeForm(mesh);
            }
        } else if (change.halfEdges.unchanged()) {
            halfEdgeForm = SharpEdges::toHalfEdgeForm(mesh);
        }
    }

    // Points.
    if (change.replacePoints) {
        mesh.swapPoints(*change.replacePoints);
    } else if (!change.appendPoints.empty()) {
        mesh.appendPoints(change.appendPoints);
    }

    // Faces.
    Topology& topology = mesh.editTopology();
    const FaceMode previousMode = topology.mode();
    bool uvSetsKept = true;
    if (deletion || added != nullptr) {
        halfEdgeForm = SharpEdges::toHalfEdgeForm(mesh) || halfEdgeForm;
        const bool cornersNeeded =
            mesh.maps().hasKind(ElementKind::HalfEdge) || mesh.maps().hasKind(ElementKind::Edge);
        Merged merged = merge(mesh, added, deletion ? &*change.deleteFaces : nullptr, change.faces.newToOld, cornersNeeded);
        if (!merged.faceOrigin.empty()) {
            change.faces.newToOld = std::move(merged.faceOrigin);
        }
        if (!merged.cornerOrigin.empty() && !change.halfEdges.newToOld) {
            change.halfEdges.newToOld = std::move(merged.cornerOrigin);
        }
        if (topology.mode() == FaceMode::Polygons) {
            topology.swapPolygons(merged.faceStart, merged.cornerVertex);
        } else {
            topology.swapCorners(merged.cornerVertex, topology.mode());
        }
        if (merged.touchesUv) {
            auto& uvSets = mesh.uvSets();
            for (std::size_t k = 0; k < merged.uvCorners.size() && k < uvSets.size(); ++k) {
                if (merged.uvCorners[k].empty()) {
                    uvSets[k].clear();
                } else {
                    uvSets[k].cornerUv.swap(merged.uvCorners[k]);
                }
            }
        }
    } else if (change.faceStart && change.cornerVertex) {
        topology.swapPolygons(*change.faceStart, *change.cornerVertex);
        const FaceMode mode = change.faceMode.value_or(previousMode);
        if (mode != FaceMode::Polygons) {
            topology.trySetUniformMode(mode);
        }
    } else if (change.cornerVertex) {
        topology.swapCorners(*change.cornerVertex, change.faceMode.value_or(previousMode));
    } else if (!change.faceStart) {
        topology.clear();
        for (auto& uv : mesh.uvSets()) {
            uv.clear();
        }
        uvSetsKept = false;
    }
    mesh.topologyChanged();

    // A complete new corner list invalidates every UV set (not only the current one)
    // unless faces and corner count stay the same.
    if (change.cornerVertex && !deletion && added == nullptr &&
        !(change.faces.unchanged() && oldCornerCount == mesh.cornerCount())) {
        for (auto& uv : mesh.uvSets()) {
            uv.clear();
        }
        uvSetsKept = false;
    }

    followNormals(mesh, change);

    MapEdit edit;
    edit.points = &change.points;
    edit.faces = &change.faces;
    edit.corners = &change.corners;
    edit.halfEdges = &change.halfEdges;
    edit.edges = &change.edges;
    edit.uvEdges = &change.uvEdges;
    edit.pointSources = change.pointSources ? &*change.pointSources : nullptr;
    edit.pointWeights = change.pointWeights ? &*change.pointWeights : nullptr;
    edit.newPointCount = mesh.pointCount();
    edit.uvSetsKept = uvSetsKept;
    mesh.maps().applyEdit(edit);
    if (halfEdgeForm) {
        SharpEdges::toEdgeForm(mesh);
    }
    if (change.extraMaps != nullptr) {
        change.extraMaps->applyEdit(edit);
    }
    if (chains && chains->rows() > 0 && (change.faces.newToOld || change.faces.unchanged())) {
        SharpEdges::fromChains(mesh, *chains, change.points, oldPointCount);
    }
    return true;
}

} // namespace areko::remesher::mesh
