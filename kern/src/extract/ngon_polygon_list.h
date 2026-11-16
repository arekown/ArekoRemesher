#pragma once

#include "mesh/csr.h"

#include <span>
#include <vector>

namespace areko::remesher::extract {

// Working list of polygons (point ids per polygon). Positions are indices
// into the joint corner list; neighbours wrap around inside their polygon.
// Removing a polygon moves the later ones forward by one; new polygons are
// appended behind all others.
class NgonPolygonList {
public:
    NgonPolygonList() = default;
    explicit NgonPolygonList(std::span<const int> polygon) { append(polygon); }

    [[nodiscard]] int polygonCount() const noexcept { return rows_.rows(); }
    [[nodiscard]] int cornerCount() const noexcept { return static_cast<int>(rows_.entries.size()); }
    [[nodiscard]] int begin(int p) const { return rows_.starts[static_cast<std::size_t>(p)]; }
    [[nodiscard]] int end(int p) const { return rows_.starts[static_cast<std::size_t>(p) + 1]; }
    [[nodiscard]] int size(int p) const { return end(p) - begin(p); }
    [[nodiscard]] int vertex(int position) const { return rows_.entries[static_cast<std::size_t>(position)]; }
    [[nodiscard]] std::span<const int> polygon(int p) const { return rows_.row(p); }
    [[nodiscard]] const mesh::Csr& rows() const noexcept { return rows_; }

    // Cyclic neighbours of a position of polygon p.
    [[nodiscard]] int next(int p, int position) const { return position + 1 >= end(p) ? begin(p) : position + 1; }
    [[nodiscard]] int prev(int p, int position) const { return position - 1 < begin(p) ? end(p) - 1 : position - 1; }
    // A position past the end of polygon p brought back into it (repeatedly minus its size).
    [[nodiscard]] int wrap(int p, int position) const;

    void remove(int p);
    void append(std::span<const int> polygon);

private:
    mesh::Csr rows_;
};

// Angle and flag per corner of a polygon list (same positions).
// Flags: -1 concave, 0 straight, 1 corner, 2 spike, 3 unclear.
struct NgonCornerData {
    std::vector<float> angle;
    std::vector<int> flag;

    // Drops the entries of the positions [first, last).
    void erase(int first, int last);
};

} // namespace areko::remesher::extract
