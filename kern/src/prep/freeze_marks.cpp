#include "freeze_marks.h"

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

const TypedMap<FreezeMark>* FreezeMarks::find(const Mesh& mesh)
{
    const AttributeMap* map = mesh.maps().find(map_names::freezeMarks);
    return dynamic_cast<const TypedMap<FreezeMark>*>(map);
}

std::vector<int> FreezeMarks::keptPositions(std::span<const int> loopPoints, const std::vector<FreezeMark>& marks)
{
    std::vector<int> kept;
    const int n = static_cast<int>(loopPoints.size());
    for (int step = 1; step <= n; ++step) {
        const int i = step % n;
        if (marks[static_cast<std::size_t>(loopPoints[static_cast<std::size_t>(i)])].value == 1) {
            kept.push_back(i);
        }
    }
    if (kept.size() < 2) {
        kept.clear();
    }
    return kept;
}

void FreezeMarks::mark(Mesh& mesh, bool onlyIfMissing)
{
    std::vector<FreezeMark> values;
    const auto* existing = find(mesh);
    if (existing != nullptr && existing->size() == static_cast<std::size_t>(mesh.pointCount())) {
        if (onlyIfMissing) {
            return;
        }
        values = existing->values();
    } else {
        values.assign(static_cast<std::size_t>(mesh.pointCount()), FreezeMark{});
    }
    for (const int v : mesh.adjacency().boundaryPoints()) {
        FreezeMark& mark = values[static_cast<std::size_t>(v)];
        if (mark.value == 0) {
            mark.value = 1;
        }
    }
    mesh.maps().set(std::make_unique<TypedMap<FreezeMark>>(std::string(map_names::freezeMarks), ElementKind::Point,
                                                           std::move(values)));
}

} // namespace areko::remesher::prep
