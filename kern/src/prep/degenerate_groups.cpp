#include "degenerate_groups.h"

#include "mesh/attribute_maps.h"
#include "mesh/measures.h"
#include "mesh/sharp_edges.h"

#include <bit>
#include <vector>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

constexpr float needleAngle = std::bit_cast<float>(0x3bab92a6u); // 0.3 degrees
constexpr float capAngle = std::bit_cast<float>(0x3d0efa35u);    // 2 degrees
constexpr int maxRounds = 10;

// One run over a triangle mesh with face groups.
class Reassignment {
public:
    Reassignment(Mesh& mesh, std::vector<int>& groups, BoolMap* sharp, float length)
        : adjacency_(mesh.adjacency()), groups_(groups), sharp_(sharp), lengths_(edgeLengths(mesh)),
          bad_(static_cast<std::size_t>(mesh.faceCount())), done_(static_cast<std::size_t>(mesh.faceCount()))
    {
        const std::vector<float> angles = cornerAngles(mesh);
        const std::vector<float>& doubleAreas = mesh.triangleGeometry().doubleAreas();
        const auto limit = static_cast<float>(static_cast<double>(length * length) * 0.2);
        for (int t = 0; t < mesh.faceCount(); ++t) {
            bool needle = false;
            int caps = 0;
            for (int i = 0; i < 3; ++i) {
                const float a = angles[static_cast<std::size_t>(3 * t + i)];
                needle = needle || a < needleAngle;
                caps += a < capAngle ? 1 : 0;
            }
            if (doubleAreas[static_cast<std::size_t>(t)] < limit && (needle || caps >= 2)) {
                bad_.set(static_cast<std::size_t>(t));
            }
        }
    }

    void run()
    {
        for (int round = 0; round < maxRounds; ++round) {
            changed_ = false;
            for (std::ptrdiff_t t = bad_.findNext(0); t >= 0; t = bad_.findNext(static_cast<std::size_t>(t) + 1)) {
                visit(static_cast<int>(t));
            }
            if (!changed_) {
                return;
            }
        }
    }

private:
    // Index -1 reads the flags of triangle 7 (signed bit access).
    [[nodiscard]] bool isBad(int t) const { return bad_.testSigned(t); }
    [[nodiscard]] bool isDone(int t) const { return done_.testSigned(t); }
    [[nodiscard]] int group(int t) const { return groups_[static_cast<std::size_t>(t)]; }

    // The other triangle at an edge with exactly two corners, or -1.
    [[nodiscard]] int neighbour(int e, int t) const { return adjacency_.otherFace(e, t); }

    void clearSharp(int e)
    {
        if (sharp_ != nullptr && e >= 0 && static_cast<std::size_t>(e) < sharp_->bits().size()) {
            sharp_->bits().set(static_cast<std::size_t>(e), false);
        }
    }

    // Every assignment marks the triangle as done; it counts as a change only
    // when the group really changes.
    void assign(int t, int value)
    {
        changed_ = changed_ || group(t) != value;
        groups_[static_cast<std::size_t>(t)] = value;
        done_.set(static_cast<std::size_t>(t));
    }

    // Takes the group of the neighbour x and clears the edge shared with it.
    void assignTo(int t, int x)
    {
        assign(t, group(x));
        clearSharp(adjacency_.commonEdge(x, t));
    }

    // Fallback: the group of the last valid neighbour, unless n is bad and open.
    void fallback(int t, int n, int longest)
    {
        if (isBad(n) && !isDone(n)) {
            return;
        }
        assign(t, last_);
        clearSharp(longest);
    }

    void visit(int t)
    {
        const int e[3] = {adjacency_.cornerEdge(3 * t), adjacency_.cornerEdge(3 * t + 1),
                          adjacency_.cornerEdge(3 * t + 2)};
        // Longest edge; ties keep the smaller corner.
        int k = -1;
        float best = -1.0f;
        if (e[0] != -1 && lengths_[static_cast<std::size_t>(e[0])] > -1.0f) {
            k = 0;
            best = lengths_[static_cast<std::size_t>(e[0])];
        }
        if (e[1] != -1 && lengths_[static_cast<std::size_t>(e[1])] > best) {
            k = 1;
            best = lengths_[static_cast<std::size_t>(e[1])];
        }
        if (e[2] != -1 && lengths_[static_cast<std::size_t>(e[2])] > best) {
            k = 2;
        }
        if (k == -1) {
            return;
        }
        const int longest = e[k];
        const int opposite = (k + 2) % 3;
        const int shortA = e[(k + 1) % 3];
        const int shortB = e[opposite];

        const int n = neighbour(longest, t);
        if (n != -1) {
            last_ = group(n);
            if (group(t) == last_) {
                if (!(isBad(n) && !isDone(n))) {
                    done_.set(static_cast<std::size_t>(t));
                }
                return;
            }
        }

        // Fan around the opposite corner: several groups, or only bad triangles besides t.
        const auto fan = adjacency_.pointCorners(adjacency_.vertex(3 * t + opposite));
        const int first = fan[0] / 3;
        bool othersBad = first == t || isBad(first);
        bool severalGroups = false;
        for (std::size_t j = 1; j < fan.size(); ++j) {
            const int f = fan[j] / 3;
            severalGroups = severalGroups || group(f) != group(first);
            othersBad = othersBad && (f == t || isBad(f));
        }
        if (!severalGroups && !othersBad) {
            return;
        }

        if (shortA < 0 || shortB < 0) {
            return;
        }
        const int a = neighbour(shortA, t);
        if (adjacency_.edgeCorners(shortB).size() != 2) {
            fallback(t, n, longest);
            return;
        }
        const int b = neighbour(shortB, t);
        if (b == -1 || a == -1) {
            fallback(t, n, longest);
            return;
        }
        const int ga = group(a);
        const int gb = group(b);
        if (ga == gb) {
            if (group(t) == ga) {
                fallback(t, n, longest);
                return;
            }
            assign(t, ga);
            clearSharp(shortA);
            clearSharp(shortB);
            return;
        }
        if (ga == group(t) || gb == group(t) || last_ == group(t)) {
            fallback(t, n, longest);
            return;
        }
        if (last_ != ga && last_ != gb && isBad(n) && !isDone(n)) {
            // n is bad and open: a good or finished short-edge neighbour, else wait.
            if (!isBad(a) || isDone(a)) {
                assignTo(t, a);
            } else if (!isBad(b) || isDone(b)) {
                assignTo(t, b);
            }
            return;
        }
        if (n == -1) {
            fallback(t, n, longest);
            return;
        }
        assignTo(t, n);
    }

    const Adjacency& adjacency_;
    std::vector<int>& groups_;
    BoolMap* sharp_;
    std::vector<float> lengths_;
    BitArray bad_;
    BitArray done_;
    int last_ = 0; // group of the last valid neighbour, kept across triangles and rounds
    bool changed_ = false;
};

} // namespace

void DegenerateGroups::run(Mesh& mesh, float length)
{
    if (mesh.topology().mode() != FaceMode::Triangles) {
        return;
    }
    // Brings the sharp edges into edge form even when the groups are missing.
    BoolMap* sharp = SharpEdges::edgeForm(mesh);
    auto* groups = mesh.maps().findAs<TypedMap<int>>(map_names::faceGroups, ElementKind::Face);
    if (groups == nullptr || groups->values().size() != static_cast<std::size_t>(mesh.faceCount())) {
        return;
    }
    Reassignment(mesh, groups->values(), sharp, length).run();
}

} // namespace areko::remesher::prep
