#include "ngon_polygon_list.h"

namespace areko::remesher::extract {

int NgonPolygonList::wrap(int p, int position) const
{
    const int n = size(p);
    while (position >= end(p)) {
        position -= n;
    }
    return position;
}

void NgonPolygonList::remove(int p)
{
    const int first = begin(p);
    const int n = size(p);
    rows_.entries.erase(rows_.entries.begin() + first, rows_.entries.begin() + first + n);
    auto& starts = rows_.starts;
    for (std::size_t q = static_cast<std::size_t>(p) + 1; q + 1 < starts.size(); ++q) {
        starts[q] = starts[q + 1] - n;
    }
    starts.pop_back();
}

void NgonPolygonList::append(std::span<const int> polygon)
{
    rows_.entries.insert(rows_.entries.end(), polygon.begin(), polygon.end());
    rows_.starts.push_back(static_cast<int>(rows_.entries.size()));
}

void NgonCornerData::erase(int first, int last)
{
    angle.erase(angle.begin() + first, angle.begin() + last);
    flag.erase(flag.begin() + first, flag.begin() + last);
}

} // namespace areko::remesher::extract
