#include "sliver_fixer.h"

#include "mesh/measures.h"
#include "smoother.h"

#include <bit>
#include <vector>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

constexpr int outerRounds = 4;
constexpr int innerRounds = 4;
constexpr float baseStep = std::bit_cast<float>(0x3e99999au); // 0.3
constexpr double smallAngle = std::bit_cast<double>(0x3f9acee9f37bebd6ull); // 1.5 degrees
constexpr double largeAngle = std::bit_cast<double>(0x40086f4291ecf2f4ull); // 175 degrees
constexpr float sameSide = 0.75f;
constexpr double keptNormal = 0.9;
constexpr int replacementMode = 2;

// Triangle geometry as the fixer sees it: values follow moved points only
// when refreshed.
class TrackedGeometry {
public:
    explicit TrackedGeometry(const Mesh& mesh) : geometry_(mesh) {}

    [[nodiscard]] const TriangleGeometry& values() const noexcept { return geometry_; }

    void moved(const std::vector<int>& points) { pending_.insert(pending_.end(), points.begin(), points.end()); }

    void refresh(const Mesh& mesh)
    {
        if (!pending_.empty()) {
            geometry_.updateAt(mesh, pending_);
            pending_.clear();
        }
    }

    void replaceDegenerateNormals(const Mesh& mesh)
    {
        geometry_.replaceDegenerateNormals(mesh, mesh.adjacency(), replacementMode);
    }

private:
    TriangleGeometry geometry_;
    std::vector<int> pending_;
};

// Points chosen in one round with their triangle and its normal at choice time.
struct Choice {
    std::vector<int> points;
    std::vector<int> faces;
    std::vector<Vec3f> normals;
};

class Fixer {
public:
    Fixer(Mesh& mesh, const BitArray& lockedPoints)
        : mesh_(mesh), adjacency_(mesh.adjacency()), locked_(lockedPoints), geometry_(mesh)
    {
    }

    void run()
    {
        for (int outer = 0; outer < outerRounds; ++outer) {
            geometry_.refresh(mesh_);
            geometry_.replaceDegenerateNormals(mesh_);
            const float limit = geometry_.values().smallAreaLimit();
            for (int inner = 1; inner <= innerRounds; ++inner) {
                // Choose with the values of the current positions.
                geometry_.refresh(mesh_);
                const Choice choice = choose(limit);
                if (choice.points.empty()) {
                    return;
                }
                const float step = inner == 1 ? baseStep : (inner == 2 ? 0.5f * baseStep : baseStep / 3.0f);
                move(choice.points, step);
                if (!turned(choice)) {
                    break;
                }
            }
        }
    }

private:
    [[nodiscard]] Choice choose(float limit) const
    {
        Choice out;
        BitArray chosen(static_cast<std::size_t>(mesh_.pointCount()));
        const auto& doubleAreas = geometry_.values().doubleAreas();
        const auto& normals = geometry_.values().normals();
        for (int f = 0; f < mesh_.faceCount(); ++f) {
            if (doubleAreas[static_cast<std::size_t>(f)] > limit) {
                continue;
            }
            float smallest = 10.0f;
            float widest = -10.0f;
            int widestCorner = 3 * f;
            for (int c = 3 * f; c < 3 * f + 3; ++c) {
                const float a = cornerAngle(mesh_, c);
                smallest = a < smallest ? a : smallest;
                if (a > widest) {
                    widest = a;
                    widestCorner = c;
                }
            }
            if (!(static_cast<double>(smallest) <= smallAngle || static_cast<double>(widest) >= largeAngle)) {
                continue;
            }
            const int v = adjacency_.vertex(widestCorner);
            if (chosen.test(static_cast<std::size_t>(v)) || !movable(v, normals[static_cast<std::size_t>(f)], normals)) {
                continue;
            }
            chosen.set(static_cast<std::size_t>(v));
            out.points.push_back(v);
            out.faces.push_back(f);
            out.normals.push_back(normals[static_cast<std::size_t>(f)]);
        }
        return out;
    }

    // Manifold neighbourhood, not locked, all triangles on the side of `normal`.
    [[nodiscard]] bool movable(int v, Vec3f normal, const std::vector<Vec3f>& normals) const
    {
        for (const int e : adjacency_.pointEdges(v)) {
            if (adjacency_.edgeCorners(e).size() != 2) {
                return false;
            }
        }
        if (locked_.testBounded(static_cast<std::size_t>(v))) {
            return false;
        }
        for (const int c : adjacency_.pointCorners(v)) {
            if (!(dot(normals[static_cast<std::size_t>(c / 3)], normal) > sameSide)) {
                return false;
            }
        }
        return true;
    }

    // Smooths the listed points together; the fold guard sees the triangle
    // normals brought up to date at the start of the step.
    void move(const std::vector<int>& points, float step)
    {
        geometry_.refresh(mesh_);
        Smoother::Options options;
        options.iterations = 1;
        options.uniformWeight = step;
        options.foldGuard = true;
        const std::vector<Vec3f> smoothed =
            Smoother(mesh_, options).applyTo(points, geometry_.values().normals());
        std::vector<Vec3f> positions;
        positions.reserve(points.size());
        for (const int v : points) {
            positions.push_back(smoothed[static_cast<std::size_t>(v)]);
        }
        mesh_.movePoints(points, positions);
        geometry_.moved(points);
    }

    // True when a chosen triangle's normal turned too far since its choice.
    [[nodiscard]] bool turned(const Choice& choice) const
    {
        const Topology& topology = mesh_.topology();
        const auto& points = mesh_.points();
        for (std::size_t i = 0; i < choice.faces.size(); ++i) {
            const int c = 3 * choice.faces[i];
            Vec3f normal;
            float doubleArea = 0.0f;
            floatTriangleNormal(points[static_cast<std::size_t>(topology.vertex(c))],
                                points[static_cast<std::size_t>(topology.vertex(c + 1))],
                                points[static_cast<std::size_t>(topology.vertex(c + 2))], normal, doubleArea);
            if (static_cast<double>(dot(normal, choice.normals[i])) < keptNormal) {
                return true;
            }
        }
        return false;
    }

    Mesh& mesh_;
    const Adjacency& adjacency_;
    const BitArray& locked_;
    TrackedGeometry geometry_;
};

} // namespace

void SliverFixer::run(Mesh& mesh, const BitArray& lockedPoints)
{
    if (mesh.topology().mode() != FaceMode::Triangles) {
        return;
    }
    Fixer(mesh, lockedPoints).run();
}

} // namespace areko::remesher::prep
