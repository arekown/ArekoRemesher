#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace areko::remesher::solver {

// One stored entry: index along the list or column, and its value.
struct SparseEntry {
    int index = 0;
    double value = 0.0;
};

class ListMatrix;

// Sparse matrix in compressed columns: column c holds the entries
// [columnStart(c), columnEnd(c)) in storage order. Columns are filled in order
// with append() and closeColumn().
class CompactMatrix {
public:
    CompactMatrix() = default;
    explicit CompactMatrix(int rows) : rows_(rows) {}

    [[nodiscard]] int rows() const noexcept { return rows_; }
    [[nodiscard]] int columns() const noexcept { return static_cast<int>(start_.size()) - 1; }
    [[nodiscard]] std::size_t nonZeros() const noexcept { return rowIndex_.size(); }
    [[nodiscard]] int columnStart(int c) const { return start_[static_cast<std::size_t>(c)]; }
    [[nodiscard]] int columnEnd(int c) const { return start_[static_cast<std::size_t>(c) + 1]; }
    [[nodiscard]] int row(int k) const { return rowIndex_[static_cast<std::size_t>(k)]; }
    [[nodiscard]] double value(int k) const { return values_[static_cast<std::size_t>(k)]; }
    double& value(int k) { return values_[static_cast<std::size_t>(k)]; }

    [[nodiscard]] std::span<const int> starts() const noexcept { return start_; }
    [[nodiscard]] std::span<const int> rowIndices() const noexcept { return rowIndex_; }
    [[nodiscard]] std::span<const double> values() const noexcept { return values_; }
    // Row indices of column c in storage order.
    [[nodiscard]] std::span<const int> columnRows(int c) const
    {
        return std::span<const int>(rowIndex_).subspan(static_cast<std::size_t>(columnStart(c)),
                                                        static_cast<std::size_t>(columnEnd(c) - columnStart(c)));
    }

    // Appends an entry to the open column; closeColumn() ends it.
    void append(int row, double value)
    {
        rowIndex_.push_back(row);
        values_.push_back(value);
    }
    void closeColumn() { start_.push_back(static_cast<int>(rowIndex_.size())); }

    // Transposed copy by counting sort: column r of the result holds the
    // columns c of this matrix that have row r, ascending.
    [[nodiscard]] CompactMatrix transposed() const;
    // Lists (e.g. rows) turned into columns the same way: column c holds
    // (list number, value) of every entry with index c, lists ascending.
    [[nodiscard]] static CompactMatrix fromListsTransposed(const ListMatrix& lists);

private:
    int rows_ = 0;
    std::vector<int> start_{0};
    std::vector<int> rowIndex_;
    std::vector<double> values_;
};

// Sparse matrix as one entry list per row or column; `otherDimension` is the
// size along the lists. Symmetric systems keep every list sorted by index.
class ListMatrix {
public:
    using List = std::vector<SparseEntry>;

    ListMatrix() = default;
    ListMatrix(int listCount, int otherDimension);

    [[nodiscard]] int listCount() const noexcept { return static_cast<int>(lists_.size()); }
    [[nodiscard]] int otherDimension() const noexcept { return other_; }
    [[nodiscard]] bool square() const noexcept { return listCount() == other_; }
    [[nodiscard]] const List& list(int i) const { return lists_[static_cast<std::size_t>(i)]; }
    List& list(int i) { return lists_[static_cast<std::size_t>(i)]; }

    // Position of the entry with `index` in a sorted list (binary search), or -1.
    [[nodiscard]] static int find(const List& list, int index) noexcept;

    // What becomes of entries with the same index when a list is ordered.
    enum class SameIndex {
        Separate,     // they stay apart, in the order the sorting leaves them
        Summed,       // they are added up in that order, each to the sum so far
        SummedNonZero // the same, but an entry equal to 0 is skipped and a sum equal to 0 removed
    };
    // Sorts a list by index (with the fixed course of steps of
    // util::Introsort, so entries with the same index always end up in the
    // same order) and treats entries with the same index as asked. Zeros take
    // part in the sorting. A value that is not a number is not 0: it is never
    // skipped and stays as a sum.
    static void order(List& list, SameIndex sameIndex);

private:
    std::vector<List> lists_;
    int other_ = 0;
};

} // namespace areko::remesher::solver
