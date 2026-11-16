#include "smoother.h"

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

// One application: weights, normals and the list of points the guard reset.
struct Smoother::Sweep {
    const Adjacency& adjacency;
    std::span<const float> weights;
    float uniformWeight;
    bool guard;
    const BitArray* guardPoints;
    std::span<const Vec3f> normals;
    std::vector<int> repair;

    [[nodiscard]] float weight(int v) const
    {
        return weights.size() == static_cast<std::size_t>(adjacency.pointCount()) ? weights[static_cast<std::size_t>(v)]
                                                                                  : uniformWeight;
    }

    // A point the sweep leaves where it is (weight 0 from the per point weights,
    // or no edges) is copied without the fold guard; the uniform weight always
    // counts as movement.
    [[nodiscard]] bool moves(int v) const
    {
        if (adjacency.pointNeighbours(v).empty()) {
            return false;
        }
        return weights.size() != static_cast<std::size_t>(adjacency.pointCount()) ||
               weights[static_cast<std::size_t>(v)] != 0.0f;
    }

    [[nodiscard]] Vec3f laplace(int v, const std::vector<Vec3f>& source) const
    {
        const Vec3f p = source[static_cast<std::size_t>(v)];
        const float w = weight(v);
        const auto neighbours = adjacency.pointNeighbours(v);
        if (w == 0.0f || neighbours.empty()) {
            return p;
        }
        Vec3f sum;
        for (const int u : neighbours) {
            const Vec3f q = source[static_cast<std::size_t>(u)];
            sum.x += q.x;
            sum.y += q.y;
            sum.z += q.z;
        }
        const float keep = 1.0f - w;
        const float c = w / static_cast<float>(neighbours.size());
        return {p.x * keep + sum.x * c, p.y * keep + sum.y * c, p.z * keep + sum.z * c};
    }

    // `p` may be an element of `positions` itself.
    void foldGuard(Vec3f& p, int v, const std::vector<Vec3f>& positions, const std::vector<Vec3f>& before)
    {
        if (guardPoints != nullptr && !guardPoints->testBounded(static_cast<std::size_t>(v))) {
            return;
        }
        int tries = 2;
        bool halved = false;
        while (true) {
            bool corrected = false;
            for (const int c : adjacency.pointCorners(v)) {
                const Vec3f a = positions[static_cast<std::size_t>(adjacency.vertex(adjacency.next(c)))];
                const Vec3f b = positions[static_cast<std::size_t>(adjacency.vertex(adjacency.prev(c)))];
                Vec3f e = b - a;
                const float l = normalizeWithLength(e);
                if (static_cast<double>(l) < 1e-12) {
                    continue;
                }
                Vec3f inward = cross(normals[static_cast<std::size_t>(c / 3)], e);
                normalizeInPlace(inward);
                const float s = dot(p - a, inward);
                if (!(s <= 0.0f)) {
                    continue;
                }
                corrected = true;
                const float t = s - l * 0.01f;
                p = {p.x - inward.x * t, p.y - inward.y * t, p.z - inward.z * t};
            }
            if (!corrected) {
                return;
            }
            if (--tries > 0) {
                continue;
            }
            const Vec3f q = before[static_cast<std::size_t>(v)];
            if (!halved) {
                p = {(p.x + q.x) * 0.5f, (p.y + q.y) * 0.5f, (p.z + q.z) * 0.5f};
                halved = true;
                tries = 2;
            } else {
                p = q;
                repair.push_back(v);
                return;
            }
        }
    }
};

Smoother::Smoother(const Mesh& mesh, const Options& options)
    : mesh_(mesh), options_(options),
      guard_(options.foldGuard && mesh.topology().mode() == FaceMode::Triangles)
{
}

std::vector<Vec3f> Smoother::applyAll(std::span<const float> weights, std::span<const Vec3f> triangleNormals) const
{
    return run(weights, nullptr, triangleNormals);
}

std::vector<Vec3f> Smoother::applyTo(std::span<const int> points, std::span<const Vec3f> triangleNormals) const
{
    return run({}, &points, triangleNormals);
}

std::vector<Vec3f> Smoother::run(std::span<const float> weights, const std::span<const int>* list,
                                 std::span<const Vec3f> triangleNormals) const
{
    const std::vector<Vec3f>& meshPoints = mesh_.points();
    std::vector<Vec3f> out = meshPoints;
    Sweep sweep{mesh_.adjacency(), weights, options_.uniformWeight, guard_, options_.guardPoints, triangleNormals, {}};
    std::vector<Vec3f> copy;
    for (int it = 0; it < options_.iterations; ++it) {
        if (it > 0) {
            copy = out;
        }
        const std::vector<Vec3f>& source = it == 0 ? meshPoints : copy;
        if (list == nullptr) {
            const int n = mesh_.pointCount();
            for (int v = 0; v < n; ++v) {
                Vec3f& p = out[static_cast<std::size_t>(v)];
                p = sweep.laplace(v, source);
                if (guard_ && sweep.moves(v)) {
                    sweep.foldGuard(p, v, out, source);
                }
            }
        } else {
            for (const int v : *list) {
                Vec3f p = sweep.laplace(v, source);
                if (guard_ && sweep.moves(v)) {
                    sweep.foldGuard(p, v, source, source);
                }
                out[static_cast<std::size_t>(v)] = p;
            }
        }
        if (sweep.repair.empty() || options_.maxRepairRuns < 0) {
            continue;
        }
        for (int k = 0;; ++k) {
            const std::vector<int> pending = std::move(sweep.repair);
            sweep.repair.clear();
            std::vector<Vec3f> moved;
            moved.reserve(pending.size());
            for (const int v : pending) {
                Vec3f p = sweep.laplace(v, out);
                if (sweep.moves(v)) {
                    sweep.foldGuard(p, v, out, out);
                }
                moved.push_back(p);
            }
            for (std::size_t i = 0; i < pending.size(); ++i) {
                out[static_cast<std::size_t>(pending[i])] = moved[i];
            }
            if (sweep.repair.size() == pending.size() || options_.maxRepairRuns <= k || sweep.repair.empty()) {
                break;
            }
        }
    }
    return out;
}

} // namespace areko::remesher::prep
