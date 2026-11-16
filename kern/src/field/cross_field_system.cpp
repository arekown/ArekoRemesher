#include "cross_field_system.h"

#include "mesh/measures.h"

#include <bit>
#include <cstdint>

namespace areko::remesher::field {

using namespace areko::remesher::mesh;
using solver::ListMatrix;
using solver::SparseEntry;

namespace {

static_assert(wide(piF) == 3.1415927410125732 && wide(fullTurnF) == 6.2831854820251465);
static_assert(CrossFieldSystem::quarterTurn == 1.5707963705062866);

} // namespace

static_assert(std::bit_cast<std::uint32_t>(2e-6f) == 0x360637bdu);

std::vector<double> CrossFieldSystem::targetAngles(const Mesh& mesh, const Curvature& curvature,
                                                   const std::vector<Vec3f>& frameE1, const std::vector<int>* faceMirror)
{
    const std::vector<Vec3f>& normals = mesh.triangleGeometry().normals();
    const std::size_t faceCount = static_cast<std::size_t>(mesh.faceCount());
    std::vector<double> out(faceCount, 0.0);
    for (std::size_t t = 0; t < faceCount; ++t) {
        const Vec3f d = projectOntoPlane(curvature.direction[t], normals[t]);
        // The angle as float, wrapped with pi and 2 pi as float, all widened.
        double angle = wide(signedVectorAngle(frameE1[t], d, cross(frameE1[t], d), normals[t]));
        if (angle > wide(piF)) {
            angle = angle - wide(fullTurnF);
        }
        out[t] = angle;
    }
    if (faceMirror != nullptr && curvature.mirrorFaceCount != 0) {
        for (std::size_t t = 0; t < faceCount && t < faceMirror->size(); ++t) {
            const int m = (*faceMirror)[t];
            if (m != -1 && curvature.mirrorFaces.testBounded(static_cast<std::size_t>(m))) {
                out[static_cast<std::size_t>(m)] = piD - out[t];
            }
        }
    }
    return out;
}

std::vector<int> CrossFieldSystem::anchorJumps(const Input& input)
{
    const Adjacency& adjacency = input.mesh.adjacency();
    const BitArray& notComputed = input.curvature.mirrorFaces;
    std::vector<int> out(static_cast<std::size_t>(adjacency.edgeCount()), 0);
    for (int e = 0; e < adjacency.edgeCount(); ++e) {
        const auto faces = adjacency.innerEdgeFaces(e);
        if (!faces) {
            continue;
        }
        const std::size_t f = static_cast<std::size_t>((*faces)[0]);
        const std::size_t g = static_cast<std::size_t>((*faces)[1]);
        if (notComputed.testBounded(f) && notComputed.testBounded(g)) {
            continue;
        }
        if (input.anchors.testBounded(f) && input.anchors.testBounded(g)) {
            out[static_cast<std::size_t>(e)] = nearestQuarterTurns(input.targetAngle[f], input.targetAngle[g],
                                                                   input.connection[static_cast<std::size_t>(e)]);
        }
    }
    return out;
}

CrossFieldSystem::CrossFieldSystem(const Input& input, const std::vector<int>& periodJump) : input_(input)
{
    const Adjacency& adjacency = input_.mesh.adjacency();
    const BitArray& notComputed = input_.curvature.mirrorFaces;
    const int faceCount = input_.mesh.faceCount();
    unknownOfFace_.assign(static_cast<std::size_t>(faceCount), -1);
    for (int t = 0; t < faceCount; ++t) {
        if (!notComputed.testBounded(static_cast<std::size_t>(t))) {
            unknownOfFace_[static_cast<std::size_t>(t)] = static_cast<int>(unknownFaces_.size());
            unknownFaces_.push_back(t);
        }
    }
    for (std::ptrdiff_t e = input_.freeEdges.findNext(0); e != -1;
         e = input_.freeEdges.findNext(static_cast<std::size_t>(e) + 1)) {
        const auto faces = adjacency.innerEdgeFaces(static_cast<int>(e));
        if (faces && unknownOf((*faces)[0]) >= 0 && unknownOf((*faces)[1]) >= 0) {
            candidates_.push_back(static_cast<int>(e));
        }
    }
    const int n = static_cast<int>(unknownFaces_.size());
    matrix_ = ListMatrix(n, n);
    rhs_.assign(static_cast<std::size_t>(n), 0.0);
    buildRows(periodJump);
}

std::array<int, 2> CrossFieldSystem::edgeUnknowns(int edge) const
{
    // A candidate edge is an inner edge.
    const auto [f, g] = input_.mesh.adjacency().innerEdgeFaces(edge).value();
    return {unknownOf(f), unknownOf(g)};
}

void CrossFieldSystem::buildRows(const std::vector<int>& periodJump)
{
    const Adjacency& adjacency = input_.mesh.adjacency();
    const double w = smoothWeight;
    for (std::size_t u = 0; u < unknownFaces_.size(); ++u) {
        const int t = unknownFaces_[u];
        ListMatrix::List& row = matrix_.list(static_cast<int>(u));
        double b = 0.0;
        double diagonal = 0.0;
        // Couples to the triangle across an edge; one that is not computed counts with the angle 0.
        auto couple = [&](int other) {
            const int v = unknownOf(other);
            if (v >= 0) {
                row.push_back({v, -w});
            } else {
                b = b - (-w) * 0.0;
            }
        };
        for (int i = 0; i < 3; ++i) {
            const int e = adjacency.cornerEdge(3 * t + i);
            if (e < 0 || input_.freeEdges.test(static_cast<std::size_t>(e))) {
                continue;
            }
            const auto faces = adjacency.innerEdgeFaces(e);
            if (!faces) {
                continue;
            }
            const auto [first, second] = *faces;
            const double jump = static_cast<double>(periodJump[static_cast<std::size_t>(e)]);
            const double kappa = wide(input_.connection[static_cast<std::size_t>(e)]);
            if (first == t) {
                couple(second);
                b = (b - (quarterTurn * w) * jump) + kappa * w;
            }
            if (second == t) {
                couple(first);
                b = ((quarterTurn * w) * jump + b) - kappa * w;
            }
            diagonal = diagonal + w;
        }
        const float weight = input_.curvature.weight[static_cast<std::size_t>(t)];
        if (input_.anchors.testBounded(static_cast<std::size_t>(t)) && carriesWeight(weight)) {
            diagonal = diagonal + wide(weight);
            b = wide(weight) * input_.targetAngle[static_cast<std::size_t>(t)] + b;
        }
        row.push_back({static_cast<int>(u), diagonal});
        // Sorted by column, entries of the same column added up in that order.
        ListMatrix::order(row, ListMatrix::SameIndex::Summed);
        rhs_[u] = b;
    }
}

void CrossFieldSystem::add(int row, int column, double value)
{
    ListMatrix::List& list = matrix_.list(row);
    auto at = list.begin();
    while (at != list.end() && at->index < column) {
        ++at;
    }
    if (at != list.end() && at->index == column) {
        at->value = at->value + value;
    } else {
        list.insert(at, SparseEntry{column, value});
    }
}

} // namespace areko::remesher::field
