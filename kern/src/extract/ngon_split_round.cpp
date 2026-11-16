#include "ngon_split_round.h"

#include "ngon_polygon_splitter.h"

#include "mesh/modification.h"

#include <utility>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

SplitRound::SplitRound(Mesh& mesh)
    : mesh_(mesh), points_(mesh), done_(static_cast<std::size_t>(mesh.faceCount()))
{
}

void SplitRound::splitFace(int face)
{
    if (done(face)) {
        return;
    }
    const Topology& topology = mesh_.topology();
    const auto& corners = topology.cornerVertex();
    NgonPolygonList list(std::span<const int>(corners.data() + topology.faceBegin(face),
                                              static_cast<std::size_t>(topology.faceSize(face))));
    queueSplit(face, list);
}

void SplitRound::queueSplit(int face, NgonPolygonList& list)
{
    PolygonSplitter::split(list, points_);
    for (int p = 0; p < list.polygonCount(); ++p) {
        const auto polygon = list.polygon(p);
        pending_.entries.insert(pending_.entries.end(), polygon.begin(), polygon.end());
        pending_.starts.push_back(static_cast<int>(pending_.entries.size()));
    }
    markDone(face);
}

void SplitRound::queuePolygons(int face, const std::vector<std::vector<int>>& polygons)
{
    for (const auto& polygon : polygons) {
        pending_.entries.insert(pending_.entries.end(), polygon.begin(), polygon.end());
        pending_.starts.push_back(static_cast<int>(pending_.entries.size()));
    }
    markDone(face);
}

void SplitRound::markDone(int face)
{
    if (!done_.testBounded(static_cast<std::size_t>(face))) {
        done_.set(static_cast<std::size_t>(face));
        ++doneCount_;
    }
}

void SplitRound::clearDone()
{
    done_ = BitArray(static_cast<std::size_t>(mesh_.faceCount()));
}

void SplitRound::commit()
{
    Mesh added;
    added.editTopology().swapPolygons(pending_.starts, pending_.entries);
    added.topologyChanged();
    pending_.reset();

    Modification change;
    change.appendPoints = points_.positions();
    change.pointSources = points_.parents();
    if (!points_.weights().empty()) {
        change.pointWeights = points_.weights();
    }
    change.deleteFaces = done_;
    change.appendFaces = &added;
    MeshModifier::apply(mesh_, change);
}

void SplitRound::reset()
{
    pending_.reset();
    if (doneCount_ != 0) {
        done_ = BitArray(static_cast<std::size_t>(mesh_.faceCount()));
        doneCount_ = 0;
    }
    points_.reset();
}

} // namespace areko::remesher::extract
