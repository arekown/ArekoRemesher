#pragma once

#include <algorithm>
#include <bit>
#include <cstddef>
#include <iterator>
#include <utility>

namespace areko::remesher::util {

// Introsort with a fixed course of steps, so that elements with equal keys end
// up in the same order with every compiler and standard library:
//  - depth limit 2 * floor(log2 n),
//  - while a range has more than 16 elements: median of (first+1, middle,
//    last-1) is swapped to first, unguarded partition around *first, the right
//    part recurses, the loop continues with the left part; at depth 0 the
//    range is heap sorted,
//  - finally insertion sort over the first 16 elements and unguarded insertion
//    for the rest (n <= 16: insertion sort over everything, which is stable).
class Introsort {
public:
    template <class It, class Less>
    static void sort(It first, It last, Less less)
    {
        const auto n = last - first;
        if (n <= 0) {
            return;
        }
        const int depth = 2 * (static_cast<int>(std::bit_width(static_cast<std::size_t>(n))) - 1);
        loop(first, last, depth, less);
        if (n > threshold) {
            insertion(first, first + threshold, less);
            for (It i = first + threshold; i != last; ++i) {
                unguardedInsert(i, less);
            }
        } else {
            insertion(first, last, less);
        }
    }

private:
    static constexpr std::ptrdiff_t threshold = 16;

    template <class It, class Less>
    static void loop(It first, It last, int depth, Less less)
    {
        while (last - first > threshold) {
            if (depth == 0) {
                heapSort(first, last, less);
                return;
            }
            --depth;
            const It middle = first + (last - first) / 2;
            medianToFirst(first, first + 1, middle, last - 1, less);
            const It cut = partition(first + 1, last, first, less);
            loop(cut, last, depth, less);
            last = cut;
        }
    }

    template <class It, class Less>
    static void medianToFirst(It result, It a, It b, It c, Less less)
    {
        using std::iter_swap;
        if (less(*a, *b)) {
            if (less(*b, *c)) {
                iter_swap(result, b);
            } else if (less(*a, *c)) {
                iter_swap(result, c);
            } else {
                iter_swap(result, a);
            }
        } else if (less(*a, *c)) {
            iter_swap(result, a);
        } else if (less(*b, *c)) {
            iter_swap(result, c);
        } else {
            iter_swap(result, b);
        }
    }

    template <class It, class Less>
    static It partition(It lo, It hi, It pivot, Less less)
    {
        while (true) {
            while (less(*lo, *pivot)) {
                ++lo;
            }
            --hi;
            while (less(*pivot, *hi)) {
                --hi;
            }
            if (!(lo < hi)) {
                return lo;
            }
            std::iter_swap(lo, hi);
            ++lo;
        }
    }

    template <class It, class Less>
    static void adjustHeap(It first, std::ptrdiff_t hole, std::ptrdiff_t len,
                           typename std::iterator_traits<It>::value_type value, Less less)
    {
        const std::ptrdiff_t top = hole;
        std::ptrdiff_t child = hole;
        while (child < (len - 1) / 2) {
            child = 2 * (child + 1);
            if (less(first[child], first[child - 1])) {
                --child;
            }
            first[hole] = std::move(first[child]);
            hole = child;
        }
        if ((len & 1) == 0 && child == (len - 2) / 2) {
            child = 2 * (child + 1);
            first[hole] = std::move(first[child - 1]);
            hole = child - 1;
        }
        std::ptrdiff_t parent = (hole - 1) / 2;
        while (hole > top && less(first[parent], value)) {
            first[hole] = std::move(first[parent]);
            hole = parent;
            parent = (hole - 1) / 2;
        }
        first[hole] = std::move(value);
    }

    template <class It, class Less>
    static void heapSort(It first, It last, Less less)
    {
        const std::ptrdiff_t len = last - first;
        if (len >= 2) {
            for (std::ptrdiff_t parent = (len - 2) / 2;; --parent) {
                adjustHeap(first, parent, len, std::move(first[parent]), less);
                if (parent == 0) {
                    break;
                }
            }
        }
        while (last - first > 1) {
            --last;
            auto value = std::move(*last);
            *last = std::move(*first);
            adjustHeap(first, 0, last - first, std::move(value), less);
        }
    }

    template <class It, class Less>
    static void unguardedInsert(It i, Less less)
    {
        auto value = std::move(*i);
        It next = i;
        --next;
        while (less(value, *next)) {
            *i = std::move(*next);
            i = next;
            --next;
        }
        *i = std::move(value);
    }

    template <class It, class Less>
    static void insertion(It first, It last, Less less)
    {
        if (first == last) {
            return;
        }
        for (It i = first + 1; i != last; ++i) {
            if (less(*i, *first)) {
                auto value = std::move(*i);
                std::move_backward(first, i, i + 1);
                *first = std::move(value);
            } else {
                unguardedInsert(i, less);
            }
        }
    }
};

} // namespace areko::remesher::util
