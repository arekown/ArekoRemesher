#include "status.h"

namespace areko::remesher::api {

const char* statusText(arm_status status) noexcept
{
    switch (status) {
    case ARM_OK:
        return "ok";
    case ARM_ERROR_INVALID_ARGUMENT:
        return "invalid argument";
    case ARM_ERROR_OUT_OF_MEMORY:
        return "out of memory";
    case ARM_ERROR_INTERNAL:
        return "internal error";
    case ARM_ERROR_NO_MESH:
        return "no mesh";
    case ARM_ERROR_NO_RESULT:
        return "no result";
    case ARM_CANCELLED:
        return "cancelled";
    case ARM_ERROR_INVALID_SETTINGS:
        return "invalid settings";
    case ARM_ERROR_STAGE_NOT_AVAILABLE:
        return "stage not available";
    }
    return "unknown status";
}

} // namespace areko::remesher::api
