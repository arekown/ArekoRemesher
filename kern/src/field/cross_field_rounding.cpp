#include "cross_field_rounding.h"

#include "util/rounding.h"

#include <cfloat>
#include <cmath>
#include <numeric>

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

CrossFieldRounding::OpenIntegers::OpenIntegers(std::size_t slots, std::vector<int> open)
    : open_(std::move(open)), closed_(slots), error_(slots, 0.0), rounded_(slots, 0.0), total_(open_.size())
{
}

std::size_t CrossFieldRounding::OpenIntegers::binOf(double error) noexcept
{
    // Upper limits of the first seven bins, as fractions of the budget.
    static constexpr std::array<double, binCount - 1> limits = {
        budget * 1e-4, budget * 1e-3, budget * 0.01, budget * 0.03, budget * 0.07, budget * 0.1, budget * 0.25};
    // An error that is not a number falls through to the last bin.
    std::size_t bin = 0;
    while (bin < binCount - 1 && !(error < limits[bin])) {
        ++bin;
    }
    return bin;
}

void CrossFieldRounding::OpenIntegers::beginRound()
{
    for (std::vector<int>& bin : bins_) {
        bin.clear();
    }
    binSum_.fill(0.0);
}

void CrossFieldRounding::OpenIntegers::rate(int index, double value)
{
    const std::size_t i = static_cast<std::size_t>(index);
    const double rounded = static_cast<double>(util::roundHalfAway(value));
    const double error = std::fabs(value - rounded);
    error_[i] = error;
    rounded_[i] = rounded;
    const std::size_t bin = binOf(error);
    bins_[bin].push_back(index);
    binSum_[bin] = error + binSum_[bin];
}

void CrossFieldRounding::OpenIntegers::close(int index)
{
    closed_.set(static_cast<std::size_t>(index));
    ++closedCount_;
}

void CrossFieldRounding::OpenIntegers::endRound()
{
    std::erase_if(open_, [this](int index) { return closed(index); });
}

CrossFieldRounding::CrossFieldRounding(CrossFieldSystem& system, const run::Progress& progress,
                                       const solver::DriverOptions& options)
    : system_(system), progress_(progress), driver_(system.matrix(), system.rhs(), *this, progress, options)
{
    const CrossFieldSystem::Input& in = system_.input();
    const std::size_t faceCount = static_cast<std::size_t>(in.mesh.faceCount());
    const std::size_t candidateCount = system_.candidateEdges().size();
    std::vector<int> candidates(candidateCount);
    std::iota(candidates.begin(), candidates.end(), 0);
    edges_ = OpenIntegers(candidateCount, std::move(candidates));
    // Weak triangles: computed, no anchor, but with a direction worth following.
    std::vector<int> weak;
    for (std::size_t t = 0; t < faceCount; ++t) {
        if (system_.unknownOf(static_cast<int>(t)) >= 0 && !in.anchors.testBounded(t) &&
            CrossFieldSystem::carriesWeight(in.curvature.weight[t])) {
            weak.push_back(static_cast<int>(t));
        }
    }
    faces_ = OpenIntegers(faceCount, std::move(weak));
}

void CrossFieldRounding::onRequest(solver::DriverRequest& request)
{
    request.noMark = true;
    const std::vector<double>& x = *request.x;
    edges_.beginRound();
    faces_.beginRound();
    touched_.clear();

    rate(x);
    if (request.round == 0) {
        freezeAnchors(x);
    }
    const int edges = roundEdges();
    const int faces = roundFaces();
    edges_.endRound();
    faces_.endRound();

    // A cancelled run ends here as cancelled.
    progress_.check();
    if (edges + faces == 0) {
        request.finished = true;
    } else {
        request.extra = &touched_;
        request.systemChanged = true;
        if (edges_.allClosed() && faces_.allClosed()) {
            request.finished = true;
        }
    }
    request.answer = solver::DriverRequest::Answer::Delivered;
}

void CrossFieldRounding::rate(const std::vector<double>& x)
{
    const CrossFieldSystem::Input& in = system_.input();
    // Across a candidate edge: the turn from the cross of its first triangle to that of its second.
    for (const int position : edges_.open()) {
        const int e = system_.candidateEdges()[static_cast<std::size_t>(position)];
        const auto [u0, u1] = system_.edgeUnknowns(e);
        const double angle0 = x[static_cast<std::size_t>(u0)];
        const double angle1 = x[static_cast<std::size_t>(u1)];
        edges_.rate(position, ((wide(in.connection[static_cast<std::size_t>(e)]) - angle0) + angle1) /
                                  CrossFieldSystem::quarterTurn);
    }
    // In a weak triangle: the turn of its cross against its curvature direction.
    for (const int face : faces_.open()) {
        faces_.rate(face, (x[static_cast<std::size_t>(system_.unknownOf(face))] -
                           in.targetAngle[static_cast<std::size_t>(face)]) /
                              CrossFieldSystem::quarterTurn);
    }
}

void CrossFieldRounding::freezeAnchors(const std::vector<double>& x)
{
    const std::vector<float>& weight = system_.input().curvature.weight;
    for (std::size_t f = 0; f < weight.size(); ++f) {
        const int u = system_.unknownOf(static_cast<int>(f));
        if (u >= 0 && weight[f] >= CrossFieldSystem::anchorWeight) {
            driver_.fixVariable(u, x[static_cast<std::size_t>(u)]);
        }
    }
}

int CrossFieldRounding::roundEdges()
{
    int count = 0;
    double used = 0.0;
    for (std::size_t b = 0; b < OpenIntegers::binCount; ++b) {
        const double sum = used + edges_.binSum(b);
        if (sum <= budget) {
            // The whole bin fits.
            for (const int position : edges_.bin(b)) {
                roundEdge(position);
                ++count;
            }
            used = sum;
            continue;
        }
        // The bin that does not fit: its edges one by one, smallest error first,
        // until the budget is used up; at least one edge per round.
        while (true) {
            int pick = -1;
            double best = DBL_MAX;
            for (const int position : edges_.bin(b)) {
                if (!edges_.closed(position) && edges_.error(position) < best) {
                    pick = position;
                    best = edges_.error(position);
                }
            }
            if (pick == -1) {
                break;
            }
            used = used + best;
            if (used > budget) {
                if (count == 0) {
                    roundEdge(pick);
                    count = 1;
                }
                break;
            }
            roundEdge(pick);
            ++count;
        }
        break;
    }
    return count;
}

int CrossFieldRounding::roundFaces()
{
    int count = 0;
    double used = 0.0;
    for (std::size_t b = 0; b < OpenIntegers::binCount; ++b) {
        for (const int face : faces_.bin(b)) {
            roundFace(face);
            ++count;
        }
        // The budget is looked at only after a whole bin.
        used = used + faces_.binSum(b);
        if (used > budget + budget) {
            break;
        }
    }
    return count;
}

void CrossFieldRounding::roundEdge(int position)
{
    const CrossFieldSystem::Input& in = system_.input();
    const int e = system_.candidateEdges()[static_cast<std::size_t>(position)];
    const auto [u0, u1] = system_.edgeUnknowns(e);
    const double w = CrossFieldSystem::smoothWeight;
    const double c =
        wide(in.connection[static_cast<std::size_t>(e)]) - edges_.rounded(position) * CrossFieldSystem::quarterTurn;
    std::vector<double>& rhs = system_.rhs();
    // Both rows take the term, also the row of an anchor frozen in round 0: the
    // anchor then follows its neighbours a little. Keeping anchors strictly at
    // their angle was measured to give clearly more triangles and irregular
    // points (crease and seam spheres), so they are left free to follow.
    const auto r0 = static_cast<std::size_t>(u0);
    const auto r1 = static_cast<std::size_t>(u1);
    system_.add(u0, u0, w);
    system_.add(u0, u1, -w);
    rhs[r0] = c * w + rhs[r0];
    touched_.push_back(u0);
    system_.add(u1, u0, -w);
    system_.add(u1, u1, w);
    rhs[r1] = rhs[r1] - c * w;
    touched_.push_back(u1);
    edges_.close(position);
}

void CrossFieldRounding::roundFace(int face)
{
    const std::size_t f = static_cast<std::size_t>(face);
    const int u = system_.unknownOf(face);
    std::vector<double>& rhs = system_.rhs();
    system_.add(u, u, 1.0);
    rhs[static_cast<std::size_t>(u)] = (faces_.rounded(face) * CrossFieldSystem::quarterTurn + system_.input().targetAngle[f]) +
                                       rhs[static_cast<std::size_t>(u)];
    touched_.push_back(u);
    faces_.close(face);
}

} // namespace areko::remesher::field
