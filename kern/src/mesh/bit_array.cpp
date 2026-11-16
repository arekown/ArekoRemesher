#include "bit_array.h"

#include <bit>

namespace areko::remesher::mesh {

void BitArray::resize(std::size_t size, bool value)
{
    const std::size_t old = size_;
    bytes_.resize((size + 7) >> 3, 0);
    size_ = size;
    if (size < old) {
        // Keep the unused tail bits of the last byte at zero so that count() stays exact.
        if ((size & 7) != 0) {
            bytes_[size >> 3] &= static_cast<std::uint8_t>((1u << (size & 7)) - 1u);
        }
        return;
    }
    for (std::size_t i = old; i < size; ++i) {
        set(i, value);
    }
}

void BitArray::clear() noexcept
{
    bytes_.clear();
    size_ = 0;
}

void BitArray::set(std::size_t i, bool value) noexcept
{
    const auto mask = static_cast<std::uint8_t>(1u << (i & 7));
    if (value) {
        bytes_[i >> 3] |= mask;
    } else {
        bytes_[i >> 3] &= static_cast<std::uint8_t>(~mask);
    }
}

void BitArray::orWith(const BitArray& other) noexcept
{
    const std::size_t n = size_ < other.size_ ? size_ : other.size_;
    for (std::size_t i = 0; i < n; ++i) {
        if (other.test(i)) {
            set(i);
        }
    }
}

namespace {

// Bit addressed by a signed index, or size when it addresses none.
std::size_t signedIndex(std::ptrdiff_t i, std::size_t size) noexcept
{
    const auto index = static_cast<std::size_t>(i < 0 && i >= -7 ? (i & 7) : i);
    return i < -7 || index >= size ? size : index;
}

} // namespace

bool BitArray::testSigned(std::ptrdiff_t i) const noexcept
{
    const std::size_t index = signedIndex(i, size_);
    return index < size_ && test(index);
}

void BitArray::setSigned(std::ptrdiff_t i, bool value) noexcept
{
    const std::size_t index = signedIndex(i, size_);
    if (index < size_) {
        set(index, value);
    }
}

void BitArray::flipAll() noexcept
{
    for (std::uint8_t& b : bytes_) {
        b = static_cast<std::uint8_t>(~b);
    }
    if ((size_ & 7) != 0) {
        bytes_[size_ >> 3] &= static_cast<std::uint8_t>((1u << (size_ & 7)) - 1u);
    }
}

std::size_t BitArray::count() const noexcept
{
    std::size_t n = 0;
    for (const std::uint8_t b : bytes_) {
        n += static_cast<std::size_t>(std::popcount(b));
    }
    return n;
}

bool BitArray::any() const noexcept
{
    for (const std::uint8_t b : bytes_) {
        if (b != 0) {
            return true;
        }
    }
    return false;
}

std::ptrdiff_t BitArray::findNext(std::size_t from) const noexcept
{
    for (std::size_t i = from; i < size_; ++i) {
        if (test(i)) {
            return static_cast<std::ptrdiff_t>(i);
        }
    }
    return -1;
}

bool operator==(const BitArray& a, const BitArray& b)
{
    return a.size_ == b.size_ && a.bytes_ == b.bytes_;
}

} // namespace areko::remesher::mesh
