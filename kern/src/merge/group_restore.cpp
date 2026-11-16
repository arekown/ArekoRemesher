#include "group_restore.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace areko::remesher::merge {

using namespace areko::remesher::mesh;

namespace {

// Value of region g in `table`, 0 for a region outside it.
template <class T>
T valueOf(const std::vector<T>& table, int g)
{
    return g >= 0 && static_cast<std::size_t>(g) < table.size() ? table[static_cast<std::size_t>(g)] : T{};
}

// One value per face from the table, looked up by the region of the face.
template <class T>
std::vector<T> lookUp(const std::vector<T>& table, const std::vector<int>& regions)
{
    std::vector<T> out;
    out.reserve(regions.size());
    for (const int g : regions) {
        out.push_back(valueOf(table, g));
    }
    return out;
}

} // namespace

void GroupRestore::run(Mesh& mesh, const prep::GroupRestoreTable& table)
{
    AttributeMaps& maps = mesh.maps();
    const TypedMap<int>* groups = maps.faceGroups();
    const auto faces = static_cast<std::size_t>(mesh.faceCount());
    if (groups == nullptr || groups->size() < faces) {
        return;
    }
    // Every value below is looked up by the region the face had before.
    const std::vector<int> regions(groups->values().begin(), groups->values().begin() + static_cast<std::ptrdiff_t>(faces));
    if (!table.groupMaterial.empty()) {
        maps.set(std::make_unique<TypedMap<int>>(std::string(map_names::materials), ElementKind::Face,
                                                 lookUp(table.groupMaterial, regions)));
    } else if (table.materialsAreGroups) {
        maps.set(std::make_unique<TypedMap<int>>(std::string(map_names::materials), ElementKind::Face, regions));
    }
    if (!table.groupFaceGroup.empty()) {
        maps.set(std::make_unique<TypedMap<int>>(std::string(map_names::faceGroups), ElementKind::Face,
                                                 lookUp(table.groupFaceGroup, regions)));
    }
    if (!table.groupSmoothing.empty()) {
        maps.set(std::make_unique<TypedMap<std::uint32_t>>(std::string(map_names::smoothingGroups), ElementKind::Face,
                                                           lookUp(table.groupSmoothing, regions)));
    }
}

} // namespace areko::remesher::merge
