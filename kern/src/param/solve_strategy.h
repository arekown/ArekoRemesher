#pragma once

#include "param/final_rounding.h"
#include "param/rounding_step.h"

#include <functional>

namespace areko::remesher::param {

// Decides round by round how the unknowns of a part become whole numbers
// while the driver solves its system:
//  1. reweighting: the first solution shows where the grid fits badly; the
//     weights are raised there and the system is solved again,
//  2. regions: asked every round until it delivers; a round without its
//     values rounds a single unknown, the one nearest to a whole number.
//     Values only for positions fixed already count as delivered, and the
//     next step goes on in the same round,
//  3. borders: once, for a part with frozen borders,
//  4. everything left is rounded at once.
// A stop request ends the run as cancelled.
class SolveStrategy final : public solver::DriverCallback {
public:
    // Rounds of reweighting before the first unknown is rounded.
    static constexpr int reweightRounds = 1;

    // `report` is called once per round with the share of the integer list
    // that is fixed (may be empty). All arguments outlive the strategy.
    SolveStrategy(const SystemContext& context, System& system, RoundingStep& regions, RoundingStep& borders,
                  const run::Progress& progress, const std::function<void(float)>& report);

    // The driver the strategy answers to; set before its run.
    void attach(const solver::ConstrainedDriver& driver) noexcept { driver_ = &driver; }

    void onRequest(solver::DriverRequest& request) override;

    // Rounds asked so far.
    [[nodiscard]] int rounds() const noexcept { return rounds_; }

private:
    enum class State { Reweight, Regions, Borders, Final };

    [[nodiscard]] SolveState state() const { return {context_, system_, *driver_, progress_}; }

    const SystemContext& context_;
    System& system_;
    RoundingStep& regions_;
    RoundingStep& borders_;
    FinalRounding final_;
    const run::Progress& progress_;
    const std::function<void(float)>& report_;
    const solver::ConstrainedDriver* driver_ = nullptr;
    State state_ = State::Reweight;
    int rounds_ = 0;
};

} // namespace areko::remesher::param
