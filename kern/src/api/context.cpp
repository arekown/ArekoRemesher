#include "context.h"

#include "run/pipeline.h"
#include "run/prep_cache.h"
#include "run/progress.h"
#include "status.h"

namespace areko::remesher::api {

void Context::setError(std::string_view message)
{
    lastError_.assign(message);
}

void Context::clearError() noexcept
{
    lastError_.clear();
}

void Context::run(arm_stage until, arm_progress_callback callback, void* user)
{
    stage_.reset();
    preparationReused_ = false;
    cancel_.store(false, std::memory_order_relaxed);
    run::Progress progress(callback, user, cancel_);
    run::Pipeline pipeline(input_, settings_.block(), progress, prepCache_.get());
    run::StageResult result = pipeline.run(until);
    preparationReused_ = result.preparationReused;
    stage_.emplace(result.stage, result.mesh, result.lines, result.mirrorSides);
    stage_->setParts(result.faceParts, result.partValues);
    stage_->setLineKinds(result.lineKinds);
    stage_->setFaceDirections(result.faceDirections);
    stage_->setCrossField(result.crossDirections, result.singularities);
    stage_->setCutEdges(result.cutEdges);
    if (result.stage == ARM_STAGE_SIZE_FIELD) {
        stage_->setSizeField(result.faceSizes, result.partSizes);
    }
    if (result.stage == ARM_STAGE_PARAMETERIZATION) {
        stage_->setGrid(result.gridCoordinates, result.gridParts);
    }
}

const StageView& Context::stage() const
{
    if (!stage_) {
        throw StatusError(ARM_ERROR_STAGE_NOT_AVAILABLE, "no stage computed");
    }
    return *stage_;
}

bool Context::preparationReused() const
{
    static_cast<void>(stage());
    return preparationReused_;
}

} // namespace areko::remesher::api
