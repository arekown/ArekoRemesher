#pragma once

#include "areko_remesher.h"

#include <atomic>
#include <span>
#include <vector>

namespace areko::remesher::run {

// Progress of one run: phases with weights, a callback and the cancel flag
// that another thread may set. Once the flag is set no further progress is
// reported and the next check stops the run.
class Progress {
public:
    Progress(arm_progress_callback callback, void* user, const std::atomic<bool>& cancel) noexcept
        : callback_(callback), user_(user), cancel_(cancel)
    {
    }

    // Relative weights of the phases of this run.
    void plan(std::span<const float> weights);
    void beginPhase(std::size_t phase, arm_stage stage);
    // Fraction 0..1 of the current phase.
    void report(float fraction);
    void finish();

    [[nodiscard]] bool cancelled() const noexcept { return cancel_.load(std::memory_order_relaxed); }
    // Throws a StatusError (cancelled) when the flag is set.
    void check() const;

private:
    arm_progress_callback callback_;
    void* user_;
    const std::atomic<bool>& cancel_;
    std::vector<float> weights_;
    float total_ = 0.0f;
    float done_ = 0.0f;
    std::size_t phase_ = 0;
    arm_stage stage_ = ARM_STAGE_INPUT;
};

} // namespace areko::remesher::run
