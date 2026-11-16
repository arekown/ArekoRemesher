#include "final_cleanup.h"

#include "mesh/degenerate_cleanup.h"
#include "prep/mesh_repair.h"

namespace areko::remesher::extract {

namespace {

// How many of the first `length` old points stay (the renumbering keeps the
// order, so they are the new points 0 .. count - 1).
std::size_t keptAmong(std::size_t length, const std::vector<int>& renumber)
{
    std::size_t count = 0;
    for (std::size_t p = 0; p < length && p < renumber.size(); ++p) {
        count += renumber[p] >= 0 ? 1 : 0;
    }
    return count;
}

// Keeps the entries of the points that stay, in their new order; a list
// shorter than the points stays as short over the points it covered.
template <class T>
void keepKept(std::vector<T>& values, const std::vector<int>& renumber)
{
    std::vector<T> out;
    out.reserve(keptAmong(values.size(), renumber));
    for (std::size_t p = 0; p < values.size() && p < renumber.size(); ++p) {
        if (renumber[p] >= 0) {
            out.push_back(values[p]);
        }
    }
    values = std::move(out);
}

void keepKept(mesh::BitArray& bits, const std::vector<int>& renumber)
{
    mesh::BitArray out(keptAmong(bits.size(), renumber));
    for (std::size_t p = 0; p < bits.size() && p < renumber.size(); ++p) {
        if (renumber[p] >= 0 && bits.test(p)) {
            out.set(static_cast<std::size_t>(renumber[p]));
        }
    }
    bits = std::move(out);
}

// A point number in its new numbering (-1 stays -1, a removed point gives -1).
int renumbered(int point, const std::vector<int>& renumber)
{
    return point >= 0 && static_cast<std::size_t>(point) < renumber.size() ? renumber[static_cast<std::size_t>(point)]
                                                                           : -1;
}

} // namespace

void FinalCleanup::run(RawMesh& raw)
{
    // Every map is kept (point colours, hard edges along the half edges).
    mesh::DegenerateCleanup::run(raw.mesh, {.keepMaps = true});
    std::vector<int> renumber;
    if (prep::MeshRepair::removeUnusedPoints(raw.mesh, &renumber) > 0) {
        OutputPoints& points = raw.points;
        keepKept(points.onBorderEdge, renumber);
        keepKept(points.reachedBorder, renumber);
        keepKept(points.onCurve, renumber);
        keepKept(points.source, renumber);
        keepKept(points.weight, renumber);
        keepKept(points.triangle, renumber);
        for (std::vector<int>& list : points.neighbours) {
            keepKept(list, renumber);
            for (int& n : list) {
                n = renumbered(n, renumber);
            }
        }
        for (int& out : raw.originalToOutput) {
            out = renumbered(out, renumber);
        }
    }
    raw.curveEdges.clear();
}

} // namespace areko::remesher::extract
