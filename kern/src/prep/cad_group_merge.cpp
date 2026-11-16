#include "cad_group_merge.h"

#include "face_borders.h"
#include "mesh/connected_parts.h"
#include "mesh/edge_chains.h"

#include <array>
#include <cmath>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

constexpr int maxPasses = 5;
constexpr float halfPi = 1.5707964f;
constexpr float rightAngleTolerance = 0.34906584f; // 20 degrees
constexpr float bothTolerance = 0.5235988f;        // 30 degrees
constexpr float pi = 3.1415927f;
constexpr float straightTolerance = 0.38397244f;   // 22 degrees
constexpr float shortShare = 0.25f;

// Chain ends: nodes numbered by first appearance (chain 0 start, chain 0 end,
// chain 1 start ...); every node lists its chains in chain order, a closed
// chain twice.
struct ChainNodes {
    std::vector<std::vector<int>> chainsAt;
    std::vector<std::array<int, 2>> endsOf;
    std::vector<int> pointOf;

    explicit ChainNodes(const Csr& chains)
    {
        std::vector<int> nodeOfPoint;
        auto nodeFor = [&](int p) {
            if (static_cast<std::size_t>(p) >= nodeOfPoint.size()) {
                nodeOfPoint.resize(static_cast<std::size_t>(p) + 1, -1);
            }
            int& n = nodeOfPoint[static_cast<std::size_t>(p)];
            if (n < 0) {
                n = static_cast<int>(chainsAt.size());
                chainsAt.emplace_back();
                pointOf.push_back(p);
            }
            return n;
        };
        for (int c = 0; c < chains.rows(); ++c) {
            const auto row = chains.row(c);
            const int a = nodeFor(row.front());
            const int b = nodeFor(row.back());
            endsOf.push_back({a, b});
            chainsAt[static_cast<std::size_t>(a)].push_back(c);
            chainsAt[static_cast<std::size_t>(b)].push_back(c);
        }
    }

    [[nodiscard]] int degree(int node) const { return static_cast<int>(chainsAt[static_cast<std::size_t>(node)].size()); }
};

double cosine(Vec3f a, Vec3f b, double lengths)
{
    return static_cast<double>(dot(a, b)) / lengths;
}

// Direction of the chain away from the node point.
Vec3f directionAt(const std::vector<Vec3f>& points, std::span<const int> chain, int nodePoint)
{
    const std::size_t n = chain.size();
    Vec3f v = chain[0] == nodePoint ? points[static_cast<std::size_t>(chain[1])] - points[static_cast<std::size_t>(chain[0])]
                                    : points[static_cast<std::size_t>(chain[n - 2])] - points[static_cast<std::size_t>(chain[n - 1])];
    normalizeInPlace(v);
    return v;
}

// The two other chains leave the node at about right angles to the chain and
// run on in a straight line. NaN values pass every limit test.
bool crossesStraight(Vec3f c, Vec3f a, Vec3f b)
{
    const double lc = lengthAsDouble(c);
    const double la = lengthAsDouble(a);
    const double lb = lengthAsDouble(b);
    auto angleToChain = [&](Vec3f d, double ld, float& deviation) {
        const double lengths = lc * ld;
        if (lengths == 0.0) {
            return false;
        }
        const double q = cosine(c, d, lengths);
        if (q <= -1.0 || q >= 1.0) {
            return false;
        }
        deviation = std::fabs(static_cast<float>(std::acos(q)) - halfPi);
        return !(deviation > rightAngleTolerance);
    };
    float da = 0.0f;
    float db = 0.0f;
    if (!angleToChain(a, la, da) || !angleToChain(b, lb, db)) {
        return false;
    }
    if (da + db > bothTolerance) {
        return false;
    }
    const double lengths = la * lb;
    if (lengths == 0.0) {
        return false;
    }
    const double q = cosine(a, b, lengths);
    if (q <= -1.0) {
        return true;
    }
    if (q >= 1.0) {
        return false;
    }
    const float gamma = static_cast<float>(std::acos(q));
    return !(pi - gamma > straightTolerance);
}

// Restore table entries of the new groups: the old group of their smallest face.
void followRenumbering(GroupRestoreTable& table, const std::vector<int>& oldGroups, const std::vector<int>& newGroups,
                       int newCount)
{
    if (!table.restore) {
        return;
    }
    std::vector<int> oldOfNew(static_cast<std::size_t>(newCount), -1);
    for (std::size_t f = 0; f < newGroups.size() && f < oldGroups.size(); ++f) {
        int& slot = oldOfNew[static_cast<std::size_t>(newGroups[f])];
        if (slot == -1) {
            slot = oldGroups[f];
        }
    }
    auto follow = [&](auto& values) {
        if (values.empty()) {
            return;
        }
        std::decay_t<decltype(values)> out(oldOfNew.size());
        for (std::size_t g = 0; g < oldOfNew.size(); ++g) {
            const int o = oldOfNew[g];
            if (o >= 0 && static_cast<std::size_t>(o) < values.size()) {
                out[g] = values[static_cast<std::size_t>(o)];
            }
        }
        values.swap(out);
    };
    follow(table.groupFaceGroup);
    follow(table.groupMaterial);
    follow(table.groupSmoothing);
    if (table.materialsAreGroups) {
        // The groups were the materials themselves: keep them as explicit entries.
        table.groupMaterial = oldOfNew;
        table.materialsAreGroups = false;
    }
}

} // namespace

void CadGroupMerge::apply(Mesh& mesh, GroupRestoreTable& table)
{
    for (int pass = 0; pass < maxPasses; ++pass) {
        const auto* groupMap = mesh.maps().findAs<TypedMap<int>>(map_names::faceGroups, ElementKind::Face);
        if (groupMap == nullptr) {
            return;
        }
        const Adjacency& adjacency = mesh.adjacency();
        const auto& points = mesh.points();
        BitArray borders = FaceBorders::valueBorders(adjacency, groupMap->values(), true);
        Csr chainEdges;
        const Csr chains = EdgeChains::build(adjacency, borders, &chainEdges);
        const ChainNodes nodes(chains);

        std::vector<float> chainLength(static_cast<std::size_t>(chains.rows()), 0.0f);
        for (int c = 0; c < chains.rows(); ++c) {
            const auto row = chains.row(c);
            float sum = 0.0f;
            for (std::size_t i = 0; i + 1 < row.size(); ++i) {
                sum += distance(points[static_cast<std::size_t>(row[i])], points[static_cast<std::size_t>(row[i + 1])]);
            }
            chainLength[static_cast<std::size_t>(c)] = sum;
        }

        std::vector<int> longChains;
        std::vector<int> strips;
        for (int c = 0; c < chains.rows(); ++c) {
            const auto row = chains.row(c);
            const auto ends = nodes.endsOf[static_cast<std::size_t>(c)];
            if (row.size() < 2 || nodes.degree(ends[0]) >= 4 || nodes.degree(ends[1]) >= 4) {
                continue;
            }
            float sum = 0.0f;
            int fulfilled = 0;
            for (const int node : ends) {
                std::vector<int> others;
                for (const int x : nodes.chainsAt[static_cast<std::size_t>(node)]) {
                    if (x != c) {
                        others.push_back(x);
                    }
                }
                if (others.size() != 2) {
                    break;
                }
                const int point = nodes.pointOf[static_cast<std::size_t>(node)];
                const Vec3f a = directionAt(points, chains.row(others[0]), point);
                const Vec3f b = directionAt(points, chains.row(others[1]), point);
                const Vec3f own = directionAt(points, row, point);
                if (!crossesStraight(own, a, b)) {
                    break;
                }
                ++fulfilled;
                sum = (chainLength[static_cast<std::size_t>(others[0])] + chainLength[static_cast<std::size_t>(others[1])]) + sum;
            }
            const float len = chainLength[static_cast<std::size_t>(c)];
            if (shortShare * sum < len) {
                longChains.push_back(c);
            } else if (fulfilled == 2) {
                strips.push_back(c);
            }
        }
        if (strips.empty()) {
            return;
        }

        for (const int c : strips) {
            for (const int e : chainEdges.row(c)) {
                borders.set(static_cast<std::size_t>(e), false);
            }
        }
        int groupCount = 0;
        std::vector<int> regrouped = ConnectedParts::faceComponents(mesh.topology(), adjacency, &borders, groupCount);
        followRenumbering(table, groupMap->values(), regrouped, groupCount);
        mesh.maps().set(std::make_unique<TypedMap<int>>(std::string(map_names::faceGroups), ElementKind::Face, std::move(regrouped)));

        std::vector<bool> marked(nodes.chainsAt.size(), false);
        for (const int c : strips) {
            for (const int node : nodes.endsOf[static_cast<std::size_t>(c)]) {
                marked[static_cast<std::size_t>(node)] = true;
            }
        }
        bool again = false;
        for (const int l : longChains) {
            for (const int node : nodes.endsOf[static_cast<std::size_t>(l)]) {
                for (const int x : nodes.chainsAt[static_cast<std::size_t>(node)]) {
                    const auto ends = nodes.endsOf[static_cast<std::size_t>(x)];
                    if (x != l && (marked[static_cast<std::size_t>(ends[0])] || marked[static_cast<std::size_t>(ends[1])])) {
                        again = true;
                    }
                }
            }
        }
        if (!again) {
            return;
        }
    }
}

} // namespace areko::remesher::prep
