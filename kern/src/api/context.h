#pragma once

#include "areko_remesher.h"
#include "input/mesh_input.h"
#include "settings/settings.h"
#include "stage_view.h"

#include <atomic>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

namespace areko::remesher::run {
class PrepCache;
}

namespace areko::remesher::api {

// State of one session behind the C handle: settings, input, the cancel flag,
// the remembered preparation it shares with the host and the data of the last
// computed stage.
class Context {
public:
    Context() = default;
    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;

    void setError(std::string_view message);
    void clearError() noexcept;
    [[nodiscard]] const std::string& lastError() const noexcept { return lastError_; }

    [[nodiscard]] settings::Settings& settings() noexcept { return settings_; }
    // A new input drops the last stage; further input parts keep it.
    [[nodiscard]] input::MeshInput& startInput() noexcept
    {
        stage_.reset();
        return input_;
    }
    [[nodiscard]] input::MeshInput& editInput() noexcept { return input_; }

    // Runs from now on use this cache (null: none); the context holds it.
    void usePrepCache(std::shared_ptr<run::PrepCache> cache) noexcept { prepCache_ = std::move(cache); }

    void run(arm_stage until, arm_progress_callback callback, void* user);
    void requestCancel() noexcept { cancel_.store(true, std::memory_order_relaxed); }

    // Throws a StatusError (stage not available) when no stage was computed.
    [[nodiscard]] const StageView& stage() const;
    // The last computed stage took the stages "input" and "prepared" from the cache.
    [[nodiscard]] bool preparationReused() const;

private:
    std::string lastError_;
    settings::Settings settings_;
    input::MeshInput input_;
    std::atomic<bool> cancel_{false};
    std::shared_ptr<run::PrepCache> prepCache_;
    std::optional<StageView> stage_;
    bool preparationReused_ = false;
};

} // namespace areko::remesher::api
