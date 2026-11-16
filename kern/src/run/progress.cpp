#include "progress.h"

#include "api/status.h"

#include <algorithm>

namespace areko::remesher::run {

void Progress::plan(std::span<const float> weights)
{
    weights_.assign(weights.begin(), weights.end());
    total_ = 0.0f;
    for (const float w : weights_) {
        total_ += w;
    }
    done_ = 0.0f;
    phase_ = 0;
}

void Progress::beginPhase(std::size_t phase, arm_stage stage)
{
    done_ = 0.0f;
    for (std::size_t i = 0; i < phase && i < weights_.size(); ++i) {
        done_ += weights_[i];
    }
    phase_ = phase;
    stage_ = stage;
    report(0.0f);
}

void Progress::report(float fraction)
{
    if (callback_ == nullptr || cancelled() || !(total_ > 0.0f)) {
        return;
    }
    const float weight = phase_ < weights_.size() ? weights_[phase_] : 0.0f;
    const float value = (done_ + weight * std::clamp(fraction, 0.0f, 1.0f)) / total_;
    callback_(user_, value, static_cast<int32_t>(stage_));
}

void Progress::finish()
{
    if (callback_ != nullptr && !cancelled()) {
        callback_(user_, 1.0f, static_cast<int32_t>(stage_));
    }
}

void Progress::check() const
{
    if (cancelled()) {
        throw api::StatusError(ARM_CANCELLED, "cancelled");
    }
}

} // namespace areko::remesher::run
