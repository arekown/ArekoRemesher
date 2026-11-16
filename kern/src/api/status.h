#pragma once

#include "areko_remesher.h"
#include "context.h"

#include <new>
#include <stdexcept>
#include <string>
#include <utility>

namespace areko::remesher::api {

// Exception carrying an interface status; thrown inside the core,
// converted to a return code at the C boundary by guarded().
class StatusError : public std::runtime_error {
public:
    StatusError(arm_status status, const std::string& message)
        : std::runtime_error(message), status_(status) {}

    [[nodiscard]] arm_status status() const noexcept { return status_; }

private:
    arm_status status_;
};

[[nodiscard]] const char* statusText(arm_status status) noexcept;

// The single place where exceptions are translated into status codes.
// Every C entry point that can fail runs its body through this function.
// The context may be null (e.g. while creating it); the message is then dropped.
template <class Body>
[[nodiscard]] arm_status guarded(Context* context, Body&& body) noexcept
{
    auto fail = [context](arm_status status, const char* message) noexcept {
        if (context != nullptr) {
            try {
                context->setError(message);
            } catch (...) {
                context->clearError();
            }
        }
        return status;
    };

    try {
        if (context != nullptr) {
            context->clearError();
        }
        std::forward<Body>(body)();
        return ARM_OK;
    } catch (const StatusError& error) {
        return fail(error.status(), error.what());
    } catch (const std::bad_alloc&) {
        return fail(ARM_ERROR_OUT_OF_MEMORY, statusText(ARM_ERROR_OUT_OF_MEMORY));
    } catch (const std::exception& error) {
        return fail(ARM_ERROR_INTERNAL, error.what());
    } catch (...) {
        return fail(ARM_ERROR_INTERNAL, statusText(ARM_ERROR_INTERNAL));
    }
}

} // namespace areko::remesher::api
