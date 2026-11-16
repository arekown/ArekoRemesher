#pragma once

#include "cross_field_system.h"

#include "mesh/bit_array.h"
#include "run/progress.h"
#include "solver/constrained_driver.h"

#include <array>
#include <vector>

namespace areko::remesher::field {

// Rounds the integers of the cross field, a batch per round of its driver:
// the quarter turn across each free edge and, for each weak triangle, the
// quarter turn against its curvature direction. Every round rates all open
// integers by their distance to the nearest integer, sorts them into bins and
// rounds the nearest ones until an error budget is used up; each rounded
// integer becomes a fixed term of the system. The first round also freezes
// the anchors at their solved angles; edges rounded later still add their
// term to an anchor's row, so the anchor follows its neighbours a little.
class CrossFieldRounding : private solver::DriverCallback {
public:
    // The driver that solves `system` is set up with `options`.
    CrossFieldRounding(CrossFieldSystem& system, const run::Progress& progress, const solver::DriverOptions& options);

    // Solves and rounds until no integer is left open; `x` gets the angle of
    // every unknown. False when the system cannot be solved.
    [[nodiscard]] bool run(std::vector<double>& x) { return driver_.run(x); }

private:
    // Error budget of the edges of one round; the weak triangles get twice as much.
    static constexpr double budget = 0.5235987901687622;

    // The integers of one kind that are still to be rounded, with what a
    // round finds out about them: the nearest whole number of each, its
    // distance to it, and the integers sorted into bins by that distance.
    class OpenIntegers {
    public:
        static constexpr std::size_t binCount = 8;

        OpenIntegers() = default;
        // The integers are numbered 0 .. slots - 1; `open` lists the ones to round, ascending.
        OpenIntegers(std::size_t slots, std::vector<int> open);

        // Empties the bins for a new round.
        void beginRound();
        // Takes the value of an open integer in this round: keeps the nearest
        // whole number and the distance to it, and puts the integer into the
        // bin of that distance.
        void rate(int index, double value);
        // The integer has been rounded.
        void close(int index);
        // Drops the rounded integers from the open ones.
        void endRound();

        [[nodiscard]] const std::vector<int>& open() const noexcept { return open_; }
        [[nodiscard]] const std::vector<int>& bin(std::size_t b) const { return bins_[b]; }
        // Sum of the distances in a bin.
        [[nodiscard]] double binSum(std::size_t b) const { return binSum_[b]; }
        [[nodiscard]] bool closed(int index) const { return closed_.test(static_cast<std::size_t>(index)); }
        [[nodiscard]] double error(int index) const { return error_[static_cast<std::size_t>(index)]; }
        [[nodiscard]] double rounded(int index) const { return rounded_[static_cast<std::size_t>(index)]; }
        [[nodiscard]] bool allClosed() const noexcept { return closedCount_ == total_; }

    private:
        // First bin whose upper limit lies above the distance; the last one otherwise.
        [[nodiscard]] static std::size_t binOf(double error) noexcept;

        std::vector<int> open_;
        mesh::BitArray closed_;
        std::vector<double> error_;
        std::vector<double> rounded_;
        std::array<std::vector<int>, binCount> bins_;
        std::array<double, binCount> binSum_{};
        std::size_t total_ = 0;
        std::size_t closedCount_ = 0;
    };

    void onRequest(solver::DriverRequest& request) override;

    void rate(const std::vector<double>& x);
    void freezeAnchors(const std::vector<double>& x);
    [[nodiscard]] int roundEdges();
    [[nodiscard]] int roundFaces();
    void roundEdge(int position);
    void roundFace(int face);

    CrossFieldSystem& system_;
    const run::Progress& progress_;
    solver::ConstrainedDriver driver_;
    OpenIntegers edges_;                 // candidate edges, by position in the candidate list
    OpenIntegers faces_;                 // weak triangles, by triangle
    std::vector<int> touched_;           // unknowns changed in this round, in order
};

} // namespace areko::remesher::field
