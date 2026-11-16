#include "solve_strategy.h"

#include "param/reweighting.h"

#include <cassert>

namespace areko::remesher::param {

using solver::DriverRequest;

SolveStrategy::SolveStrategy(const SystemContext& context, System& system, RoundingStep& regions,
                             RoundingStep& borders, const run::Progress& progress,
                             const std::function<void(float)>& report)
    : context_(context), system_(system), regions_(regions), borders_(borders), progress_(progress), report_(report)
{
}

void SolveStrategy::onRequest(DriverRequest& request)
{
    // A stop request ends the run as cancelled, in whatever step it arrives.
    progress_.check();
    // attach() comes before the run.
    assert(driver_ != nullptr);
    rounds_ = request.round + 1;
    if (report_) {
        const std::size_t positions = driver_->integerUnknowns().size();
        report_(positions == 0 ? 1.0f : static_cast<float>(request.fixedCount) / static_cast<float>(positions));
    }

    if (state_ == State::Reweight) {
        if (request.changedRounds < reweightRounds) {
            Reweighting::apply(context_, system_, *request.x, progress_);
            request.systemChanged = true;
            request.allowLocal = false;
            if (reweightRounds <= request.changedRounds + 1) {
                state_ = State::Regions;
            }
            request.answer = DriverRequest::Answer::Delivered;
            return;
        }
        state_ = State::Regions;
    }
    if (state_ == State::Regions) {
        const RoundingStep::Outcome outcome = regions_.select(state(), request);
        if (outcome == RoundingStep::Outcome::Failed) {
            // Without values of the regions the driver rounds one unknown
            // itself; the regions are asked again in the next round.
            request.answer = DriverRequest::Answer::None;
            return;
        }
        state_ = context_.part().frozenBorderMarks() != nullptr ? State::Borders : State::Final;
        if (outcome == RoundingStep::Outcome::Delivered) {
            request.answer = DriverRequest::Answer::Delivered;
            return;
        }
        // Values of the regions only for positions fixed already: the next
        // step goes on in the same round.
    }
    if (state_ == State::Borders) {
        // Only a part with frozen borders gets here. What counts is the
        // selection, not the answer of the step.
        (void)borders_.select(state(), request);
        state_ = State::Final;
        if (!request.selection.empty()) {
            request.answer = DriverRequest::Answer::Delivered;
            return;
        }
    }
    (void)final_.select(state(), request);
    request.answer = DriverRequest::Answer::Delivered;
}

} // namespace areko::remesher::param
