#include "small_region_merge.h"

#include "face_borders.h"
#include "mesh/connected_parts.h"
#include "mesh/measures.h"
#include "mesh/sharp_edges.h"

#include <algorithm>
#include <bit>
#include <cmath>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

// Neighbour groups and shape values of one region.
struct RegionShape {
    std::vector<int> groups;
    std::vector<int> regions;
    std::vector<float> shared;
    float angleSum = 0.0f;   // sum of the corner angles above 20 degrees
    int above20 = 0;
    int above110 = 0;
    float maxAngle = 0.0f;
};

class RegionAnalysis {
public:
    RegionAnalysis(const Mesh& mesh, const BitArray& borders, const std::vector<int>& faceRegion,
                   const std::vector<int>& groups)
        : mesh_(mesh), adjacency_(mesh.adjacency()), borders_(borders), faceRegion_(faceRegion), groups_(groups),
          lengths_(edgeLengths(mesh))
    {
    }

    // Shape of the region with the given faces; the border measure runs on.
    RegionShape analyse(std::span<const int> faces)
    {
        RegionShape s;
        const auto& points = mesh_.points();
        std::vector<int> at;
        for (const int f : faces) {
            for (int c = 3 * f; c < 3 * f + 3; ++c) {
                const int v = adjacency_.vertex(c);
                at.clear();
                for (const int e : adjacency_.pointEdges(v)) {
                    if (!borders_.testBounded(static_cast<std::size_t>(e))) {
                        continue;
                    }
                    at.push_back(e);
                    const float l = lengths_[static_cast<std::size_t>(e)];
                    border_ += l;
                    neighbour(s, e, f, l);
                }
                const Vec3f p = points[static_cast<std::size_t>(v)];
                for (std::size_t i = 0; i < at.size(); ++i) {
                    for (std::size_t j = i + 1; j < at.size(); ++j) {
                        angle(s, p, at[i], at[j], v);
                    }
                }
            }
        }
        border_ = border_ * 0.5f;
        return s;
    }

    [[nodiscard]] float border() const noexcept { return border_; }

private:
    void neighbour(RegionShape& s, int e, int f, float l) const
    {
        const int other = adjacency_.otherFace(e, f);
        if (other < 0) {
            return;
        }
        const int group = groups_[static_cast<std::size_t>(other)];
        for (std::size_t i = 0; i < s.groups.size(); ++i) {
            if (s.groups[i] == group) {
                s.shared[i] += l;
                return;
            }
        }
        s.groups.push_back(group);
        s.regions.push_back(faceRegion_[static_cast<std::size_t>(other)]);
        s.shared.push_back(l);
    }

    void angle(RegionShape& s, Vec3f v, int ei, int ej, int vertex) const
    {
        const auto& points = mesh_.points();
        const Vec3f ci = points[static_cast<std::size_t>(adjacency_.otherEnd(ei, vertex))];
        const Vec3f cj = points[static_cast<std::size_t>(adjacency_.otherEnd(ej, vertex))];
        const float li = lengths_[static_cast<std::size_t>(ei)];
        const float lj = lengths_[static_cast<std::size_t>(ej)];
        const float x = ((v.x - cj.x) / lj) * ((ci.x - v.x) / li);
        const float y = ((v.y - cj.y) / lj) * ((ci.y - v.y) / li);
        const float z = ((v.z - cj.z) / lj) * ((ci.z - v.z) / li);
        const double cosine = static_cast<double>((x + y) + z);
        const float phi = static_cast<float>(angleOfCosine(cosine));
        if (static_cast<double>(phi) > 0.3490658503988659) {
            ++s.above20;
            s.angleSum += phi;
        }
        if (static_cast<double>(phi) > 1.9198621771937625) {
            ++s.above110;
        }
        if (phi > s.maxAngle) {
            s.maxAngle = phi;
        }
    }

    const Mesh& mesh_;
    const Adjacency& adjacency_;
    const BitArray& borders_;
    const std::vector<int>& faceRegion_;
    const std::vector<int>& groups_;
    std::vector<float> lengths_;
    float border_ = 0.0f;
};

int clampedQuadCount(int n)
{
    if (n == 0) {
        return 5000;
    }
    if (n < 1000) {
        return 1000;
    }
    return n > 9999 ? 10000 : n;
}

} // namespace

void SmallRegionMerge::run(Mesh& mesh, const PrepState& state)
{
    if (mesh.topology().mode() != FaceMode::Triangles) {
        return;
    }
    BoolMap* sharp = SharpEdges::edgeForm(mesh);
    auto* groupMap = mesh.maps().findAs<TypedMap<int>>(map_names::faceGroups, ElementKind::Face);
    if (groupMap == nullptr) {
        return;
    }
    std::vector<int>& groups = groupMap->values();
    const Adjacency& adjacency = mesh.adjacency();

    BitArray borders = FaceBorders::valueBorders(adjacency, groups, false);
    if (sharp != nullptr) {
        borders.orWith(sharp->bits());
    }
    int regionCount = 0;
    const std::vector<int> faceRegion =
        ConnectedParts::faceComponents(mesh.topology(), adjacency, &borders, regionCount);
    const Csr regionFaces = Csr::fromValues(faceRegion, regionCount);
    const auto& doubleAreas = mesh.triangleGeometry().doubleAreas();
    std::vector<float> area(static_cast<std::size_t>(regionCount), 0.0f);
    for (int r = 0; r < regionCount; ++r) {
        for (const int f : regionFaces.row(r)) {
            area[static_cast<std::size_t>(r)] += doubleAreas[static_cast<std::size_t>(f)];
        }
    }

    const float total = mesh.triangleGeometry().totalArea();
    const float quads = static_cast<float>(clampedQuadCount(state.quadCount));
    const float threshold =
        static_cast<float>((2.0 * static_cast<double>(total)) / static_cast<double>(quads * 15.0f));

    RegionAnalysis analysis(mesh, borders, faceRegion, groups);
    int targetGroup = 0;
    for (int r = 0; r < regionCount; ++r) {
        const float a = area[static_cast<std::size_t>(r)];
        if (a > threshold) {
            continue;
        }
        const auto faces = regionFaces.row(r);
        const RegionShape s = analysis.analyse(faces);
        const std::size_t k = s.groups.size();
        if (k <= 1) {
            if (k == 1) {
                targetGroup = s.groups[0];
            }
            continue;
        }
        auto large = [&](std::size_t i) {
            const float neighbour = area[static_cast<std::size_t>(s.regions[i])];
            return neighbour > threshold && neighbour > a + a;
        };
        const std::size_t faceLimit = static_cast<double>(threshold) * 0.01 > static_cast<double>(a) ? 20 : 5;
        int target = -1;
        if (k <= 3 && faces.size() <= faceLimit) {
            float best = -1.0f;
            for (std::size_t i = 0; i < k; ++i) {
                if (large(i) && best < s.shared[i]) {
                    target = s.regions[i];
                    best = s.shared[i];
                    targetGroup = s.groups[i];
                }
            }
        } else {
            bool ambiguous = false;
            for (std::size_t i = 0; i < k && !ambiguous; ++i) {
                if (!large(i)) {
                    continue;
                }
                if (target >= 0) {
                    ambiguous = true;
                } else {
                    target = s.regions[i];
                    targetGroup = s.groups[i];
                }
            }
            if (ambiguous) {
                continue;
            }
        }
        if (target == -1) {
            continue;
        }
        const float p = analysis.border();
        const float p2 = p * p;
        const bool thin = p2 / 200.0f > a;
        const bool jagged = p2 / 40.0f > a &&
                            (s.above110 > 2 || (static_cast<double>(s.maxAngle) > 2.6179938779914944 &&
                                                static_cast<double>(s.angleSum) > 8.726646259971648));
        if (!thin && !jagged) {
            continue;
        }
        for (const int f : faces) {
            groups[static_cast<std::size_t>(f)] = targetGroup;
        }
        if (sharp == nullptr) {
            continue;
        }
        BitArray& bits = sharp->bits();
        for (const int f : faces) {
            for (int c = 3 * f; c < 3 * f + 3; ++c) {
                const int e = adjacency.cornerEdge(c);
                if (e < 0 || static_cast<std::size_t>(e) >= bits.size()) {
                    continue;
                }
                const int other = adjacency.otherFace(e, f);
                if (other < 0) {
                    continue;
                }
                if (groups[static_cast<std::size_t>(other)] == targetGroup) {
                    bits.set(static_cast<std::size_t>(e), false);
                }
            }
        }
    }
}

} // namespace areko::remesher::prep
