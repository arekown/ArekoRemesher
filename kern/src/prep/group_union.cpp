#include "group_union.h"

#include "face_borders.h"
#include "mesh/connected_parts.h"
#include "mesh/sharp_edges.h"

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

void unite(BitArray& into, const BitArray& from)
{
    for (std::size_t e = 0; e < from.size() && e < into.size(); ++e) {
        if (from.test(e)) {
            into.set(e);
        }
    }
}

} // namespace

GroupRestoreTable GroupUnion::apply(Mesh& mesh)
{
    GroupRestoreTable table;
    AttributeMaps& maps = mesh.maps();
    const auto* materials = maps.findAs<TypedMap<int>>(map_names::materials, ElementKind::Face);
    const auto* smoothing = maps.findAs<TypedMap<std::uint32_t>>(map_names::smoothingGroups, ElementKind::Face);
    const auto* groups = maps.findAs<TypedMap<int>>(map_names::faceGroups, ElementKind::Face);
    if (materials == nullptr && smoothing == nullptr) {
        return table;
    }
    if (groups == nullptr && smoothing == nullptr) {
        maps.set(std::make_unique<TypedMap<int>>(std::string(map_names::faceGroups), ElementKind::Face, materials->values()));
        table.materialsAreGroups = true;
        table.restore = true;
        return table;
    }

    const Adjacency& adjacency = mesh.adjacency();
    BitArray borders(static_cast<std::size_t>(adjacency.edgeCount()));
    if (groups != nullptr) {
        unite(borders, FaceBorders::valueBorders(adjacency, groups->values(), false));
    }
    if (materials != nullptr) {
        unite(borders, FaceBorders::valueBorders(adjacency, materials->values(), false));
    }
    const BitArray groupBorders = borders;
    if (smoothing != nullptr) {
        unite(borders, FaceBorders::smoothingBorders(adjacency, smoothing->values()));
    }
    int regionCount = 0;
    std::vector<int> regions = ConnectedParts::faceComponents(mesh.topology(), adjacency, &borders, regionCount);

    // Host values of the smallest face of every region.
    std::vector<int> smallest(static_cast<std::size_t>(regionCount), -1);
    for (std::size_t f = 0; f < regions.size(); ++f) {
        int& slot = smallest[static_cast<std::size_t>(regions[f])];
        if (slot < 0) {
            slot = static_cast<int>(f);
        }
    }
    auto fill = [&](auto& out, const auto* map) {
        if (map == nullptr) {
            return;
        }
        out.resize(smallest.size());
        for (std::size_t k = 0; k < smallest.size(); ++k) {
            const auto f = static_cast<std::size_t>(smallest[k]);
            out[k] = f < map->values().size() ? map->values()[f] : typename std::decay_t<decltype(out)>::value_type{};
        }
    };
    fill(table.groupFaceGroup, groups);
    fill(table.groupMaterial, materials);
    fill(table.groupSmoothing, smoothing);

    BitArray sharp;
    int added = 0;
    if (smoothing != nullptr) {
        sharp = SharpEdges::extendableEdgeBits(mesh);
        for (std::size_t e = 0; e < borders.size(); ++e) {
            if (!borders.test(e) || groupBorders.test(e)) {
                continue;
            }
            const auto corners = adjacency.edgeCorners(static_cast<int>(e));
            if (corners.size() != 2) {
                continue;
            }
            if (regions[static_cast<std::size_t>(adjacency.cornerFace(corners[0]))] !=
                regions[static_cast<std::size_t>(adjacency.cornerFace(corners[1]))]) {
                continue;
            }
            if (!sharp.test(e)) {
                sharp.set(e);
                ++added;
            }
        }
    }
    maps.set(std::make_unique<TypedMap<int>>(std::string(map_names::faceGroups), ElementKind::Face, std::move(regions)));
    if (added > 0) {
        SharpEdges::setEdgeBits(mesh, std::move(sharp));
    }
    table.restore = true;
    return table;
}

} // namespace areko::remesher::prep
