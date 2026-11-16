#include "border_sizes.h"

#include "centre_path_search.h"

#include "prep/freeze_marks.h"

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

namespace {

// The seam ends at this many target edge lengths from the border.
constexpr float seamWidth = 6.0f;

} // namespace

BorderSizes::BorderSizes(const parts::Part& part, const std::vector<FreezeMark>& marks, float borderDensityInverse)
    : part_(part), adjacency_(part.mesh.adjacency()), cuts_(*part.cutGraph), edgeLength_(part.edgeLength)
{
    collectGiven(marks, borderDensityInverse);
    measureDistances();
    fix();
    markSwapped();
}

void BorderSizes::collectGiven(const std::vector<FreezeMark>& marks, float borderDensityInverse)
{
    const float factor = borderDensityInverse != 0.0f ? 1.0f / borderDensityInverse : 1.0f;
    prep::FreezeMarks::forEachKeptChord(part_.mesh, marks, [&](const prep::FreezeMarks::KeptChord& chord) {
        chord.forEachEdge([&](int e) {
            if (part_.lockedBorderEdges.testBounded(static_cast<std::size_t>(e))) {
                return;
            }
            const int c = adjacency_.edgeCorners(e)[0];
            given_.push_back({c / 3, cuts_.curveFollowsU.testBounded(static_cast<std::size_t>(c)), chord.length * factor});
        });
    });
}

void BorderSizes::measureDistances()
{
    const std::size_t faceCount = static_cast<std::size_t>(part_.mesh.faceCount());
    data_.previousFace.assign(faceCount, -1);
    data_.borderCorner.assign(faceCount, -1);
    CentrePathSearch search(part_.mesh);
    // A triangle with two border edges keeps the later one.
    for (const int e : adjacency_.boundaryEdges()) {
        const int c = adjacency_.edgeCorners(e)[0];
        const int t = c / 3;
        search.start(t);
        data_.previousFace[static_cast<std::size_t>(t)] = t;
        data_.borderCorner[static_cast<std::size_t>(t)] = c;
    }
    search.run([&](int n, int t, int) {
        data_.previousFace[static_cast<std::size_t>(n)] = t;
        data_.borderCorner[static_cast<std::size_t>(n)] = data_.borderCorner[static_cast<std::size_t>(t)];
    });
    data_.distance.resize(faceCount);
    for (std::size_t t = 0; t < faceCount; ++t) {
        data_.distance[t] = search.length(static_cast<int>(t));
    }
}

void BorderSizes::fix()
{
    // A size stays exactly where the border gave it.
    const std::size_t faceCount = static_cast<std::size_t>(part_.mesh.faceCount());
    data_.fixedU = BitArray(faceCount);
    data_.fixedV = BitArray(faceCount);
    for (const Given& g : given_) {
        const auto t = static_cast<std::size_t>(g.face);
        if (!cuts_.mirrorHalf.testBounded(t)) {
            (g.alongU ? data_.fixedU : data_.fixedV).set(t);
        }
    }
}

bool BorderSizes::inSeam(std::size_t t) const
{
    return !cuts_.mirrorHalf.testBounded(t) && data_.distance[t] < edgeLength_ * seamWidth && data_.borderCorner[t] != -1;
}

void BorderSizes::markSwapped()
{
    // The paths are final once the distances are measured.
    const std::size_t faceCount = static_cast<std::size_t>(part_.mesh.faceCount());
    swapped_ = BitArray(faceCount);
    for (std::size_t t = 0; t < faceCount; ++t) {
        if (inSeam(t) && pathSwaps(static_cast<int>(t), data_.borderCorner[t] / 3)) {
            swapped_.set(t);
        }
    }
}

bool BorderSizes::pathSwaps(int t, int borderFace) const
{
    bool odd = false;
    int from = t;
    int to = data_.previousFace[static_cast<std::size_t>(t)];
    if (t == borderFace || to == -1) {
        return odd;
    }
    while (true) {
        const int e = adjacency_.commonEdge(from, to);
        if (e == -1) {
            break;
        }
        if ((cuts_.residualJump[static_cast<std::size_t>(e)] & 1) != 0) {
            odd = !odd;
        }
        // The path ends at its border triangle, which is its own predecessor.
        const int next = data_.previousFace[static_cast<std::size_t>(to)];
        if (to == borderFace || next == -1 || next == to) {
            break;
        }
        from = to;
        to = next;
    }
    return odd;
}

void BorderSizes::apply(std::vector<float>& sizeU, std::vector<float>& sizeV) const
{
    for (const Given& g : given_) {
        (g.alongU ? sizeU : sizeV)[static_cast<std::size_t>(g.face)] = g.size;
    }
    const float h = edgeLength_;
    const float two = h + h;
    const float six = h * seamWidth;
    for (std::size_t t = 0; t < sizeU.size(); ++t) {
        if (!inSeam(t)) {
            continue;
        }
        const float d = data_.distance[t];
        const auto f0 = static_cast<std::size_t>(data_.borderCorner[t] / 3);
        const bool swap = swapped(t);
        // The border size as it is close to the border, fading into h farther away.
        auto seam = [&](float border) { return d < two ? border : ((d - two) * (h - border)) / (six - two) + border; };
        if (data_.fixedU.test(f0)) {
            (swap ? sizeV : sizeU)[t] = seam(sizeU[f0]);
        }
        if (data_.fixedV.test(f0)) {
            (swap ? sizeU : sizeV)[t] = seam(sizeV[f0]);
        }
    }
}

} // namespace areko::remesher::field
