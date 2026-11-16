#pragma once

#include "attribute_map.h"

#include <memory>
#include <string_view>
#include <vector>

namespace areko::remesher::mesh {

// Names of the maps the core knows.
namespace map_names {
inline constexpr std::string_view faceGroups = "face_groups";
inline constexpr std::string_view materials = "materials";
inline constexpr std::string_view smoothingGroups = "smoothing_groups";
inline constexpr std::string_view sharpEdges = "sharp_edges";
inline constexpr std::string_view colors = "colors";
inline constexpr std::string_view density = "density";
inline constexpr std::string_view freezeMarks = "freeze_marks";
} // namespace map_names

// What a mesh change tells the maps (see MeshModifier).
struct MapEdit {
    const IndexMap* points = nullptr;
    const IndexMap* faces = nullptr;
    const IndexMap* corners = nullptr;
    const IndexMap* halfEdges = nullptr;
    const IndexMap* edges = nullptr;
    const IndexMap* uvEdges = nullptr;
    const Csr* pointSources = nullptr;
    const std::vector<float>* pointWeights = nullptr;
    int newPointCount = 0;
    bool uvSetsKept = true;
};

// Element counts of a mesh that is appended to another one.
struct ElementCounts {
    int points = 0;
    int faces = 0;
    int corners = 0;
    int edges = 0;
};

// All maps of a mesh in insertion order; looked up linearly by name.
class AttributeMaps {
public:
    AttributeMaps() = default;
    AttributeMaps(const AttributeMaps& other);
    AttributeMaps& operator=(const AttributeMaps& other);
    AttributeMaps(AttributeMaps&&) noexcept = default;
    AttributeMaps& operator=(AttributeMaps&&) noexcept = default;

    [[nodiscard]] std::size_t count() const noexcept { return maps_.size(); }
    [[nodiscard]] AttributeMap& at(std::size_t i) { return *maps_[i]; }
    [[nodiscard]] const AttributeMap& at(std::size_t i) const { return *maps_[i]; }

    [[nodiscard]] AttributeMap* find(std::string_view name);
    [[nodiscard]] const AttributeMap* find(std::string_view name) const;
    [[nodiscard]] AttributeMap* find(std::string_view name, ElementKind kind);
    [[nodiscard]] const AttributeMap* find(std::string_view name, ElementKind kind) const;
    template <class M>
    [[nodiscard]] M* findAs(std::string_view name, ElementKind kind)
    {
        return dynamic_cast<M*>(find(name, kind));
    }
    template <class M>
    [[nodiscard]] const M* findAs(std::string_view name, ElementKind kind) const
    {
        return dynamic_cast<const M*>(find(name, kind));
    }
    [[nodiscard]] bool hasKind(ElementKind kind) const;
    // The face group map (int values per face), or null.
    [[nodiscard]] TypedMap<int>* faceGroups() { return findAs<TypedMap<int>>(map_names::faceGroups, ElementKind::Face); }
    [[nodiscard]] const TypedMap<int>* faceGroups() const
    {
        return findAs<TypedMap<int>>(map_names::faceGroups, ElementKind::Face);
    }

    // Replaces a map of the same name in place, otherwise appends.
    AttributeMap& set(std::unique_ptr<AttributeMap> map);
    void remove(std::string_view name, ElementKind kind);
    void removeKind(ElementKind kind);
    void clear() noexcept { maps_.clear(); }

    // Follows a mesh change, from the last map to the first:
    //  points: unchanged state keeps; interpolation rows that match the number
    //    of new points append; otherwise int maps grow with -1, others with
    //    their default; without sources a renumbering applies, else removed,
    //  faces / half edges: unchanged state keeps, renumbering applies, else removed,
    //  edges: kept only with unchanged state,
    //  corners: renumbered by their origin (corner or half edge renumbering),
    //  UV points: kept while the UV sets are kept,
    //  UV edges: while the UV sets are kept, like faces (unchanged keeps,
    //    renumbering applies, else removed); otherwise removed.
    void applyEdit(const MapEdit& edit);

    // Appends the maps of another mesh: each own map takes the values of the
    // same-named map of the same kind; missing or rejected ones grow by the
    // element count for points, faces, corners and edges. Maps only present in
    // `other` are not taken over, unless the own element counts are given:
    // then such a map of points, faces, corners, edges or half edges joins
    // first with default values for the own elements, and an own half edge
    // map grows by the corner count of `other` as well. UV edge maps are
    // removed afterwards.
    void appendFrom(const AttributeMaps& other, const ElementCounts& otherCounts,
                    const ElementCounts* ownCounts = nullptr);

private:
    std::vector<std::unique_ptr<AttributeMap>> maps_;
};

} // namespace areko::remesher::mesh
