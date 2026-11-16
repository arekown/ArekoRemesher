#include "source_weights.h"

#include "trace_node.h"

namespace areko::remesher::extract {

using mesh::Vec3f;

SourceWeights SourceWeights::of(const mesh::Mesh& triangles, std::span<const Vec3f> points,
                                std::span<const PointSource> sources)
{
    SourceWeights out;
    const std::vector<Vec3f>& at = triangles.points();
    const mesh::Topology& topology = triangles.topology();
    const mesh::Adjacency& adjacency = triangles.adjacency();
    auto add = [&](int point, float weight) {
        out.rows.entries.push_back(point);
        out.weights.push_back(weight);
    };
    for (std::size_t i = 0; i < points.size(); ++i) {
        const PointSource s = i < sources.size() ? sources[i] : PointSource{};
        const Vec3f p = points[i];
        if (s.type == node_kind::point) {
            add(s.id, 1.0f);
        } else if (s.type == node_kind::edge) {
            const auto ends = adjacency.edgeVertices(s.id);
            const Vec3f a = at[static_cast<std::size_t>(ends[0])];
            const Vec3f d = at[static_cast<std::size_t>(ends[1])] - a;
            const float t = mesh::dot(p - a, d) / mesh::dot(d, d);
            add(ends[0], 1.0f - t);
            add(ends[1], t);
        } else if (s.type == 4) {
            const int c = 3 * s.id;
            const int q0 = topology.vertex(c);
            const int q1 = topology.vertex(c + 1);
            const int q2 = topology.vertex(c + 2);
            const Vec3f o = at[static_cast<std::size_t>(q0)];
            const Vec3f e1 = at[static_cast<std::size_t>(q1)] - o;
            const Vec3f e2 = at[static_cast<std::size_t>(q2)] - o;
            const Vec3f w = p - o;
            const float d11 = mesh::dot(e1, e1);
            const float b1 = mesh::dot(e1, w);
            const float d12 = mesh::dot(e1, e2);
            const float d22 = mesh::dot(e2, e2);
            const float b2 = mesh::dot(e2, w);
            const float inverse = 1.0f / (d11 * d22 - d12 * d12);
            const float u = (d22 * b1 - d12 * b2) * inverse;
            const float v = (b2 * d11 - d12 * b1) * inverse;
            add(q0, (1.0f - u) - v);
            add(q1, u);
            add(q2, v);
        }
        out.rows.starts.push_back(static_cast<int>(out.rows.entries.size()));
    }
    return out;
}

} // namespace areko::remesher::extract
