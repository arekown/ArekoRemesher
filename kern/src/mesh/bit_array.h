#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace areko::remesher::mesh {

// Fixed-size bit set; bit i lives in byte i >> 3 at position i & 7.
class BitArray {
public:
    BitArray() = default;
    explicit BitArray(std::size_t size, bool value = false) { resize(size, value); }

    [[nodiscard]] std::size_t size() const noexcept { return size_; }
    [[nodiscard]] bool empty() const noexcept { return size_ == 0; }

    // New bits get `value`; existing bits keep theirs.
    void resize(std::size_t size, bool value = false);
    void clear() noexcept;

    [[nodiscard]] bool test(std::size_t i) const noexcept { return (bytes_[i >> 3] >> (i & 7)) & 1u; }
    void set(std::size_t i, bool value = true) noexcept;

    // Bounded access: indices at or above size() read false and are not written.
    [[nodiscard]] bool testBounded(std::size_t i) const noexcept { return i < size_ && test(i); }
    void setBounded(std::size_t i, bool value = true) noexcept
    {
        if (i < size_) {
            set(i, value);
        }
    }

    // Sets every bit that is set in `other`, up to the smaller of both sizes.
    void orWith(const BitArray& other) noexcept;
    // Inverts every bit below size().
    void flipAll() noexcept;

    // Signed access: an index below size() addresses its bit; a negative index
    // addresses bit (i & 7) of the first byte (-1 therefore means bit 7), which
    // counts only while that bit lies below size(). Other indices read false and
    // are not written.
    [[nodiscard]] bool testSigned(std::ptrdiff_t i) const noexcept;
    void setSigned(std::ptrdiff_t i, bool value = true) noexcept;

    [[nodiscard]] std::size_t count() const noexcept;
    [[nodiscard]] bool any() const noexcept;

    // Index of the first set bit at or after `from`, or -1.
    [[nodiscard]] std::ptrdiff_t findNext(std::size_t from) const noexcept;
    // Calls visit(i) for every set bit i, ascending. The search goes on
    // behind i after each call: a bit the visit changes further up counts as
    // it is when the search gets there.
    template <class Visit>
    void forEachSet(Visit visit) const
    {
        for (std::ptrdiff_t i = findNext(0); i != -1; i = findNext(static_cast<std::size_t>(i) + 1)) {
            visit(static_cast<std::size_t>(i));
        }
    }

    friend bool operator==(const BitArray& a, const BitArray& b);

private:
    std::vector<std::uint8_t> bytes_;
    std::size_t size_ = 0;
};

} // namespace areko::remesher::mesh
