#include "import_attributes.h"

#include "cad_group_merge.h"
#include "normal_split_edges.h"
#include "seam_welder.h"

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

// Removes the first map of this name, whatever its kind.
void removeMap(AttributeMaps& maps, std::string_view name)
{
    if (const AttributeMap* map = maps.find(name)) {
        maps.remove(name, map->kind());
    }
}

} // namespace

GroupRestoreTable ImportAttributes::apply(Mesh& mesh, const arm_settings& settings)
{
    AttributeMaps& maps = mesh.maps();
    if (!settings.follow_materials) {
        removeMap(maps, map_names::materials);
    }
    if (!settings.follow_face_groups) {
        removeMap(maps, map_names::faceGroups);
    }
    if (!settings.follow_smoothing_groups) {
        removeMap(maps, map_names::smoothingGroups);
    }
    if (!settings.follow_sharp_edges) {
        removeMap(maps, map_names::sharpEdges);
    }
    GroupRestoreTable table = GroupUnion::apply(mesh);
    if (settings.merge_cad_groups) {
        CadGroupMerge::apply(mesh, table);
    }
    if (settings.follow_normal_splits) {
        NormalSplitEdges::apply(mesh);
    }
    if (settings.weld_seams) {
        SeamWelder::apply(mesh, settings.weld_as_crease != 0);
    }
    return table;
}

} // namespace areko::remesher::prep
