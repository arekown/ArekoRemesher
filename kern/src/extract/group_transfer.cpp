#include "group_transfer.h"

#include "mesh/attribute_map.h"
#include "mesh/attribute_maps.h"

#include <memory>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

constexpr int maxRounds = 50;
constexpr float noValue = -1e20f;

void setFaceGroups(Mesh& mesh, std::vector<int> values)
{
    mesh.maps().set(std::make_unique<TypedMap<int>>(std::string(map_names::faceGroups), ElementKind::Face,
                                                    std::move(values)));
}

class GroupChoice {
public:
    GroupChoice(const Mesh& triangles, const std::vector<int>& groups)
        : adjacency_(triangles.adjacency()), points_(triangles.points()), groups_(groups)
    {
    }

    // Group of a corner whose point comes from `source`, looking along
    // `direction` where the triangles around it differ; false when the corner
    // gives none.
    bool cornerGroup(const PointSource& source, Vec3f direction, int& group) const
    {
        if (source.type == 4) {
            group = groups_[static_cast<std::size_t>(source.id)];
            return true;
        }
        if (source.type == 0 || source.id < 0) {
            return false;
        }
        const std::span<const int> corners =
            source.type == 1 ? adjacency_.pointCorners(source.id) : adjacency_.edgeCorners(source.id);
        if (shared(corners, group)) {
            return true;
        }
        const int t = source.type == 1 ? alongPoint(source.id, direction) : alongEdge(source.id, direction);
        if (t == -1) {
            return false;
        }
        group = groups_[static_cast<std::size_t>(t)];
        return true;
    }

private:
    // All triangles of the corners in one group (none: false).
    bool shared(std::span<const int> corners, int& group) const
    {
        if (corners.empty()) {
            return false;
        }
        group = groups_[static_cast<std::size_t>(corners.front() / 3)];
        for (const int c : corners.subspan(1)) {
            if (groups_[static_cast<std::size_t>(c / 3)] != group) {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] Vec3f at(int v) const { return points_[static_cast<std::size_t>(v)]; }

    // The triangle at point a whose two sides at a point most along `r`
    // (their sum, normalized in float); first one on ties.
    int alongPoint(int a, Vec3f r) const
    {
        const Vec3f pa = at(a);
        int best = -1;
        float bestValue = noValue;
        for (const int c : adjacency_.pointCorners(a)) {
            Vec3f s = (at(adjacency_.vertex(adjacency_.prev(c))) - pa) + (at(adjacency_.vertex(adjacency_.next(c))) - pa);
            normalizeInPlace(s);
            const float value = dot(s, r);
            if (bestValue < value) {
                bestValue = value;
                best = c / 3;
            }
        }
        return best;
    }

    // The triangle at edge e whose corner before the edge lies most along `r`,
    // seen from the middle of the edge across the edge direction.
    int alongEdge(int e, Vec3f r) const
    {
        const auto ends = adjacency_.edgeVertices(e);
        const Vec3f a = at(ends[0]);
        const Vec3f b = at(ends[1]);
        Vec3f axis = b - a;
        normalizeInPlace(axis);
        int best = -1;
        float bestValue = noValue;
        for (const int c : adjacency_.edgeCorners(e)) {
            const Vec3f p = at(adjacency_.vertex(adjacency_.prev(c)));
            const Vec3f h = projectOntoPlane({p.x - (a.x + b.x) * 0.5f, p.y - (a.y + b.y) * 0.5f, p.z - (a.z + b.z) * 0.5f},
                                             axis);
            const float value = dot(h, r);
            if (bestValue < value) {
                bestValue = value;
                best = c / 3;
            }
        }
        return best;
    }

    const Adjacency& adjacency_;
    const std::vector<Vec3f>& points_;
    const std::vector<int>& groups_;
};

// Open faces take the group most of their neighbours have that are not open
// (across edges with two corners that are no curve edges).
void resolveOpen(const Mesh& mesh, const BitArray& curveEdges, const std::vector<int>& open, std::vector<int>& groups)
{
    const Adjacency& adjacency = mesh.adjacency();
    const Topology& topology = mesh.topology();
    BitArray isOpen(groups.size());
    for (const int f : open) {
        isOpen.set(static_cast<std::size_t>(f));
    }
    auto remaining = open.size();
    int rounds = maxRounds;
    std::vector<int> seen;
    std::vector<int> votes;
    while (true) {
        bool changed = false;
        for (const int f : open) {
            if (!isOpen.test(static_cast<std::size_t>(f))) {
                continue;
            }
            seen.clear();
            votes.clear();
            for (int c = topology.faceBegin(f); c < topology.faceEnd(f); ++c) {
                const int e = adjacency.cornerEdge(c);
                if (e < 0 || curveEdges.testBounded(static_cast<std::size_t>(e))) {
                    continue;
                }
                const auto corners = adjacency.edgeCorners(e);
                if (corners.size() != 2) {
                    continue;
                }
                const int other = c == corners[0] ? corners[1] : (c == corners[1] ? corners[0] : -1);
                if (other == -1) {
                    continue;
                }
                const int g = adjacency.cornerFace(other);
                if (isOpen.testBounded(static_cast<std::size_t>(g))) {
                    continue;
                }
                const int value = groups[static_cast<std::size_t>(g)];
                std::size_t k = 0;
                while (k < seen.size() && seen[k] != value) {
                    ++k;
                }
                if (k == seen.size()) {
                    seen.push_back(value);
                    votes.push_back(1);
                } else {
                    ++votes[k];
                }
            }
            int most = 0;
            int chosen = -1;
            for (std::size_t k = 0; k < votes.size(); ++k) {
                if (most < votes[k]) {
                    most = votes[k];
                    chosen = static_cast<int>(k);
                }
            }
            if (chosen != -1) {
                groups[static_cast<std::size_t>(f)] = seen[static_cast<std::size_t>(chosen)];
                isOpen.set(static_cast<std::size_t>(f), false);
                --remaining;
                changed = true;
            }
        }
        if (remaining == 0 || !changed || --rounds == 0) {
            break;
        }
    }
}

} // namespace

void GroupTransfer::uniform(const Mesh& triangles, RawMesh& raw)
{
    const TypedMap<int>* groups = triangles.maps().faceGroups();
    if (groups == nullptr || groups->values().empty()) {
        return;
    }
    setFaceGroups(raw.mesh, std::vector<int>(static_cast<std::size_t>(raw.mesh.faceCount()), groups->values().front()));
}

void GroupTransfer::run(const Mesh& triangles, RawMesh& raw)
{
    const TypedMap<int>* map = triangles.maps().faceGroups();
    if (map == nullptr) {
        return;
    }
    const GroupChoice choice(triangles, map->values());
    const Mesh& mesh = raw.mesh;
    const Topology& topology = mesh.topology();
    const auto& points = mesh.points();
    const int faceCount = topology.faceCount();
    std::vector<int> groups(static_cast<std::size_t>(faceCount), 0);
    std::vector<int> open;
    for (int f = 0; f < faceCount; ++f) {
        const int begin = topology.faceBegin(f);
        const int end = topology.faceEnd(f);
        int found = -1;
        bool differs = false;
        for (int c = begin; c < end && !differs; ++c) {
            const int v = topology.vertex(c);
            const Vec3f p = points[static_cast<std::size_t>(v)];
            const Vec3f next = points[static_cast<std::size_t>(topology.vertex(c + 1 == end ? begin : c + 1))];
            const Vec3f previous = points[static_cast<std::size_t>(topology.vertex(c == begin ? end - 1 : c - 1))];
            const Vec3f direction = ((next - p) + previous) - p;
            int group = 0;
            if (!choice.cornerGroup(raw.points.source[static_cast<std::size_t>(v)], direction, group)) {
                continue;
            }
            if (found == -1) {
                found = group;
            } else if (group != found) {
                differs = true;
            }
        }
        if (differs || found < 0) {
            open.push_back(f);
        } else {
            groups[static_cast<std::size_t>(f)] = found;
        }
    }
    if (!open.empty()) {
        resolveOpen(mesh, raw.curveEdges, open, groups);
    }
    if (!groups.empty()) {
        setFaceGroups(raw.mesh, std::move(groups));
    }
}

} // namespace areko::remesher::extract
