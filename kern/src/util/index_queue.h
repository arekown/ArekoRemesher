#pragma once

#include <cstddef>
#include <vector>

namespace areko::remesher::util {

// First-in-first-out queue of distinct indices 0 … n-1 kept as successor
// links: next(i) is -1 while i is not queued, -2 for the last element,
// otherwise its successor. Indices already queued are not inserted again.
class IndexQueue {
public:
    IndexQueue() = default;
    explicit IndexQueue(std::size_t n) { reset(n); }

    // Empties the queue for indices 0 … n-1.
    void reset(std::size_t n);

    [[nodiscard]] bool empty() const noexcept { return head_ < 0; }
    [[nodiscard]] bool contains(int i) const noexcept { return next_[static_cast<std::size_t>(i)] != -1; }

    void pushBack(int i);
    void pushFront(int i);
    // Removes and returns the front index; the queue must not be empty.
    int pop();
    // Replaces `out` with the queued indices from the front; the queue is empty afterwards.
    void drainTo(std::vector<int>& out);
    void clear() noexcept;

private:
    // Visits the queued indices from the front and unlinks them.
    template <class Visit>
    void drain(Visit visit)
    {
        for (int i = head_; i >= 0;) {
            visit(i);
            int& link = next_[static_cast<std::size_t>(i)];
            i = link;
            link = -1;
        }
        head_ = -1;
        tail_ = -1;
    }

    std::vector<int> next_;
    int head_ = -1;
    int tail_ = -1;
};

} // namespace areko::remesher::util
