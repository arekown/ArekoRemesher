#include "index_queue.h"

namespace areko::remesher::util {

void IndexQueue::reset(std::size_t n)
{
    next_.assign(n, -1);
    head_ = -1;
    tail_ = -1;
}

void IndexQueue::pushBack(int i)
{
    int& link = next_[static_cast<std::size_t>(i)];
    if (link != -1) {
        return;
    }
    link = -2;
    if (head_ >= 0) {
        next_[static_cast<std::size_t>(tail_)] = i;
    } else {
        head_ = i;
    }
    tail_ = i;
}

void IndexQueue::pushFront(int i)
{
    int& link = next_[static_cast<std::size_t>(i)];
    if (link != -1) {
        return;
    }
    if (head_ < 0) {
        link = -2;
        tail_ = i;
    } else {
        link = head_;
    }
    head_ = i;
}

int IndexQueue::pop()
{
    const int i = head_;
    int& link = next_[static_cast<std::size_t>(i)];
    head_ = link >= 0 ? link : -1;
    if (head_ < 0) {
        tail_ = -1;
    }
    link = -1;
    return i;
}

void IndexQueue::drainTo(std::vector<int>& out)
{
    out.clear();
    drain([&out](int i) { out.push_back(i); });
}

void IndexQueue::clear() noexcept
{
    drain([](int) {});
}

} // namespace areko::remesher::util
