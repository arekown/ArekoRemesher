#include "face_closing.h"

#include <algorithm>
#include <limits>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

// Whether a node lies on an edge of a crease curve.
bool onCurveEdge(const TraceNode& node, const CurveMarks& curves)
{
    return node.kind == node_kind::edge && curves.edges.testBounded(static_cast<std::size_t>(node.reference));
}

// Removes the entries at `first` and at `second` (two different places).
template <class T>
void eraseTwo(std::vector<T>& list, int first, int second)
{
    if (first < second) {
        list.erase(list.begin() + first, list.begin() + second + 1);
    } else {
        list.erase(list.begin() + first);
        list.erase(list.begin() + second);
    }
}

} // namespace

int FaceClosing::number(const TraceNode& node) noexcept
{
    return node.mergeTarget >= 0 ? node.mergeTarget : node.id;
}

void FaceClosing::close(FaceLoop& loop, Trace& trace)
{
    if (!merge(loop, trace)) {
        return;
    }
    std::vector<int>& k = loop.nodes;
    const int taken = static_cast<int>(std::count_if(k.begin(), k.end(), [&](int node) {
        return (trace.nodes[node].flags & node_flag::frozenBorder) != 0;
    }));
    int n = static_cast<int>(k.size());
    if (n < 3) {
        return;
    }
    if (n - taken < 5) {
        emitSmall(k, trace);
        return;
    }
    while (n > 5) {
        int spike = -1;
        for (int i = 0; i < n && spike < 0; ++i) {
            if (k[static_cast<std::size_t>(i)] == k[static_cast<std::size_t>((i + 2) % n)]) {
                spike = i;
            }
        }
        if (spike < 0) {
            break;
        }
        if (n == 6) {
            eraseTwo(k, (spike + 1) % 6, (spike + 2) % 6);
            emitRange(k, 0, static_cast<int>(k.size()) - 1, nullptr, trace);
            return;
        }
        std::vector<bool> inner(static_cast<std::size_t>(n), false);
        emitRange(k, spike + 1, spike + 4, &inner, trace);
        emitRange(k, spike - 2, spike + 1, &inner, trace);
        std::vector<int> kept;
        for (int i = 0; i < n; ++i) {
            if (!inner[static_cast<std::size_t>(i)]) {
                kept.push_back(k[static_cast<std::size_t>(i)]);
            }
        }
        k.swap(kept);
        n = static_cast<int>(k.size());
        if (n - taken < 5) {
            emitSmall(k, trace);
            return;
        }
    }
    emitSplit(k, trace);
}

bool FaceClosing::merge(FaceLoop& loop, Trace& trace) const
{
    std::vector<int>& k = loop.nodes;
    std::vector<std::uint8_t>& degenerate = loop.degenerate;
    const CurveMarks& curves = trace.curves;
    int n = static_cast<int>(k.size());
    int i = 0;
    while (i < n) {
        if (degenerate[static_cast<std::size_t>(i)] == 0) {
            ++i;
            continue;
        }
        const int ia = i;
        const int ib = (i + 1) % n;
        TraceNode& a = trace.nodes[k[static_cast<std::size_t>(ia)]];
        TraceNode& b = trace.nodes[k[static_cast<std::size_t>(ib)]];
        const bool curveA = onCurveEdge(a, curves);
        const bool curveB = onCurveEdge(b, curves);
        if (curveA && curveB &&
            curves.edgeChain[static_cast<std::size_t>(a.reference)] !=
                curves.edgeChain[static_cast<std::size_t>(b.reference)]) {
            ++i;
            continue;
        }
        bool smaller = a.id < b.id;
        bool keepA = false;
        if (a.seed()) {
            // A walk from a seed that comes back over its second node keeps that one.
            if (n > 2 && i == 0 && k[1] == k[static_cast<std::size_t>(n - 1)]) {
                smaller = false;
            }
            keepA = curveB ? curveA : (smaller || curveA);
        } else {
            keepA = curveB ? (!smaller && curveA) : (smaller || curveA);
        }
        if (frozenMarks_ != nullptr) {
            auto kept = [&](const TraceNode& node) {
                return node.kind == node_kind::point && (*frozenMarks_)[static_cast<std::size_t>(node.reference)].value == 1;
            };
            const bool keptA = kept(a);
            const bool keptB = kept(b);
            if (keptA && keptB) {
                ++i;
                continue;
            }
            if (keptA != keptB) {
                if (keptA ? curveB : curveA) {
                    ++i;
                    continue;
                }
                keepA = keptA;
            }
        }
        if (keepA) {
            b.mergeTarget = number(a);
            k.erase(k.begin() + ib);
        } else {
            a.mergeTarget = number(b);
            k.erase(k.begin() + ia);
        }
        degenerate.erase(degenerate.begin() + i);
        --n;
    }
    if (n == 3) {
        int onPoints = 0;
        int degenerateSides = 0;
        for (int j = 0; j < n; ++j) {
            degenerateSides += degenerate[static_cast<std::size_t>(j)] != 0 ? 1 : 0;
            onPoints += trace.nodes[k[static_cast<std::size_t>(j)]].kind == node_kind::point ? 1 : 0;
        }
        if (onPoints == 0 && degenerateSides > 0) {
            return false;
        }
    }
    return true;
}

void FaceClosing::emitSmall(const std::vector<int>& nodes, Trace& trace)
{
    const std::size_t n = nodes.size();
    std::vector<int> numbers;
    numbers.reserve(n);
    for (const int node : nodes) {
        numbers.push_back(number(trace.nodes[node]));
    }
    // A face that doubles back on itself is left out.
    for (std::size_t i = 0; i < n; ++i) {
        if (numbers[i] == numbers[(i + 2) % n]) {
            return;
        }
    }
    trace.polygons.entries.insert(trace.polygons.entries.end(), numbers.begin(), numbers.end());
    endPolygon(trace);
}

void FaceClosing::emitRange(const std::vector<int>& nodes, int from, int to, std::vector<bool>* inner, Trace& trace)
{
    const int n = static_cast<int>(nodes.size());
    const int first = (from + n) % n;
    const int last = (to + n) % n;
    for (int i = first; i != last; i = (i + 1) % n) {
        trace.polygons.entries.push_back(number(trace.nodes[nodes[static_cast<std::size_t>(i)]]));
        if (inner != nullptr && i != first) {
            (*inner)[static_cast<std::size_t>(i)] = true;
        }
    }
    trace.polygons.entries.push_back(number(trace.nodes[nodes[static_cast<std::size_t>(last)]]));
    endPolygon(trace);
}

void FaceClosing::emitSplit(const std::vector<int>& nodes, Trace& trace)
{
    std::vector<int> numbers;
    std::vector<Vec3f> positions;
    numbers.reserve(nodes.size());
    positions.reserve(nodes.size());
    for (const int node : nodes) {
        numbers.push_back(number(trace.nodes[node]));
        positions.push_back(trace.nodes[node].position);
    }
    for (const std::array<int, 4>& quad : splitIntoQuads(numbers, std::move(positions))) {
        trace.polygons.entries.insert(trace.polygons.entries.end(), quad.begin(), quad.end());
        endPolygon(trace);
    }
    trace.polygons.entries.insert(trace.polygons.entries.end(), numbers.begin(), numbers.end());
    endPolygon(trace);
}

void FaceClosing::endPolygon(Trace& trace)
{
    trace.polygons.starts.push_back(static_cast<int>(trace.polygons.entries.size()));
}

std::vector<std::array<int, 4>> FaceClosing::splitIntoQuads(std::vector<int>& numbers, std::vector<Vec3f> positions)
{
    std::vector<std::array<int, 4>> quads;
    for (;;) {
        const int n = static_cast<int>(numbers.size());
        if (n < 1) {
            return quads;
        }
        int best = -1;
        float shortest = std::numeric_limits<float>::max();
        for (int i = 0; i < n; ++i) {
            const Vec3f d = positions[static_cast<std::size_t>(i)] - positions[static_cast<std::size_t>((i + 3) % n)];
            const float d2 = dot(d, d);
            if (d2 < shortest) {
                shortest = d2;
                best = i;
            }
        }
        if (best < 0) {
            return quads;
        }
        auto at = [&](int k) { return numbers[static_cast<std::size_t>((best + k) % n)]; };
        quads.push_back({at(0), at(1), at(2), at(3)});
        eraseTwo(numbers, (best + 1) % n, (best + 2) % n);
        // The positions are cut together with the numbers.
        eraseTwo(positions, (best + 1) % n, (best + 2) % n);
        if (numbers.size() < 5) {
            return quads;
        }
    }
}

} // namespace areko::remesher::extract
