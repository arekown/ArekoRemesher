#include "sparse_matrix.h"

#include "util/introsort.h"

#include <algorithm>

namespace areko::remesher::solver {

namespace {

// Counting sort of the entries of `outerCount` outer lists into compressed
// columns of `columns` columns: column c receives (outer, value) for every entry
// with index c, outer lists in ascending order, each in its stored order.
template <class Entries>
CompactMatrix transposeLists(int outerCount, int columns, Entries entries)
{
    std::vector<int> start(static_cast<std::size_t>(columns) + 1, 0);
    for (int o = 0; o < outerCount; ++o) {
        entries(o, [&](int index, double) { ++start[static_cast<std::size_t>(index) + 1]; });
    }
    for (std::size_t c = 1; c < start.size(); ++c) {
        start[c] += start[c - 1];
    }
    std::vector<int> rows(static_cast<std::size_t>(start.back()));
    std::vector<double> values(rows.size());
    std::vector<int> next(start.begin(), start.end() - 1);
    for (int o = 0; o < outerCount; ++o) {
        entries(o, [&](int index, double value) {
            const auto k = static_cast<std::size_t>(next[static_cast<std::size_t>(index)]++);
            rows[k] = o;
            values[k] = value;
        });
    }
    CompactMatrix out(outerCount);
    for (int c = 0; c < columns; ++c) {
        for (int k = start[static_cast<std::size_t>(c)]; k < start[static_cast<std::size_t>(c) + 1]; ++k) {
            out.append(rows[static_cast<std::size_t>(k)], values[static_cast<std::size_t>(k)]);
        }
        out.closeColumn();
    }
    return out;
}

} // namespace

CompactMatrix CompactMatrix::transposed() const
{
    return transposeLists(columns(), rows_, [this](int c, auto&& visit) {
        for (int k = columnStart(c); k < columnEnd(c); ++k) {
            visit(row(k), value(k));
        }
    });
}

CompactMatrix CompactMatrix::fromListsTransposed(const ListMatrix& lists)
{
    return transposeLists(lists.listCount(), lists.otherDimension(), [&lists](int o, auto&& visit) {
        for (const SparseEntry& e : lists.list(o)) {
            visit(e.index, e.value);
        }
    });
}

ListMatrix::ListMatrix(int listCount, int otherDimension)
    : lists_(static_cast<std::size_t>(listCount)), other_(otherDimension)
{
}

int ListMatrix::find(const List& list, int index) noexcept
{
    const auto it = std::lower_bound(list.begin(), list.end(), index,
                                     [](const SparseEntry& e, int i) { return e.index < i; });
    return it != list.end() && it->index == index ? static_cast<int>(it - list.begin()) : -1;
}

void ListMatrix::order(List& list, SameIndex sameIndex)
{
    util::Introsort::sort(list.begin(), list.end(),
                          [](const SparseEntry& a, const SparseEntry& b) { return a.index < b.index; });
    if (sameIndex == SameIndex::Separate) {
        return;
    }
    const bool withoutZeros = sameIndex == SameIndex::SummedNonZero;
    std::size_t kept = 0;
    for (std::size_t i = 0; i < list.size(); ++i) {
        const SparseEntry entry = list[i];
        if (withoutZeros && entry.value == 0.0) {
            continue;
        }
        if (kept > 0 && list[kept - 1].index == entry.index) {
            const double sum = entry.value + list[kept - 1].value;
            list[kept - 1].value = sum;
            if (withoutZeros && sum == 0.0) {
                --kept;
            }
        } else {
            list[kept++] = entry;
        }
    }
    list.resize(kept);
}

} // namespace areko::remesher::solver
