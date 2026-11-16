#pragma once

#include "index_queue.h"

#include <cfloat>
#include <cstddef>
#include <vector>

namespace areko::remesher::util {

// Search for the smallest reach of the indices 0 … n-1 from a set of starts:
// first in, first out, and an index whose reach gets strictly smaller later
// is queued again. Which neighbours an index has and what the step to them
// costs is up to the visit of run().
class ReachQueue {
public:
    explicit ReachQueue(std::size_t count) : reach_(count, FLT_MAX), queue_(count) {}

    // Reach of i; FLT_MAX while i is not reached.
    [[nodiscard]] float reach(int i) const { return reach_[static_cast<std::size_t>(i)]; }
    [[nodiscard]] bool reached(int i) const { return reach(i) != FLT_MAX; }
    [[nodiscard]] bool empty() const noexcept { return queue_.empty(); }

    // A start: reach 0, queued at the back.
    void start(int i)
    {
        reach_[static_cast<std::size_t>(i)] = 0.0f;
        queue_.pushBack(i);
    }

    // Offers i the reach `offered`. It is taken when it is strictly smaller
    // than the reach of i (a reach that is not a number never is); i is then
    // queued at the back, or at the front. Returns whether it was taken.
    bool offer(int i, float offered, bool atFront = false)
    {
        float& reach = reach_[static_cast<std::size_t>(i)];
        if (!(offered < reach)) {
            return false;
        }
        reach = offered;
        if (atFront) {
            queue_.pushFront(i);
        } else {
            queue_.pushBack(i);
        }
        return true;
    }

    // Calls visit(i) for the queued indices in order until none is left; the
    // visit offers reaches to the neighbours of i.
    template <class Visit>
    void run(Visit visit)
    {
        while (!queue_.empty()) {
            visit(queue_.pop());
        }
    }

private:
    std::vector<float> reach_;
    IndexQueue queue_;
};

} // namespace areko::remesher::util
