#include "attribute_maps.h"

namespace areko::remesher::mesh {

AttributeMaps::AttributeMaps(const AttributeMaps& other)
{
    *this = other;
}

AttributeMaps& AttributeMaps::operator=(const AttributeMaps& other)
{
    if (this != &other) {
        std::vector<std::unique_ptr<AttributeMap>> copy;
        copy.reserve(other.maps_.size());
        for (const auto& map : other.maps_) {
            copy.push_back(map->clone());
        }
        maps_ = std::move(copy);
    }
    return *this;
}

AttributeMap* AttributeMaps::find(std::string_view name)
{
    for (auto& map : maps_) {
        if (map->name() == name) {
            return map.get();
        }
    }
    return nullptr;
}

const AttributeMap* AttributeMaps::find(std::string_view name) const
{
    return const_cast<AttributeMaps*>(this)->find(name);
}

AttributeMap* AttributeMaps::find(std::string_view name, ElementKind kind)
{
    for (auto& map : maps_) {
        if (map->name() == name && map->kind() == kind) {
            return map.get();
        }
    }
    return nullptr;
}

const AttributeMap* AttributeMaps::find(std::string_view name, ElementKind kind) const
{
    return const_cast<AttributeMaps*>(this)->find(name, kind);
}

bool AttributeMaps::hasKind(ElementKind kind) const
{
    for (const auto& map : maps_) {
        if (map->kind() == kind) {
            return true;
        }
    }
    return false;
}

AttributeMap& AttributeMaps::set(std::unique_ptr<AttributeMap> map)
{
    for (auto& slot : maps_) {
        if (slot->name() == map->name()) {
            slot = std::move(map);
            return *slot;
        }
    }
    maps_.push_back(std::move(map));
    return *maps_.back();
}

void AttributeMaps::remove(std::string_view name, ElementKind kind)
{
    std::erase_if(maps_, [&](const auto& map) { return map->name() == name && map->kind() == kind; });
}

void AttributeMaps::removeKind(ElementKind kind)
{
    std::erase_if(maps_, [kind](const auto& map) { return map->kind() == kind; });
}

namespace {

// true = keep the map, false = remove it.
bool followIndexMap(AttributeMap& map, const IndexMap* indexMap)
{
    if (indexMap == nullptr) {
        return false;
    }
    if (indexMap->unchanged()) {
        return true;
    }
    return indexMap->present() && map.reorder(*indexMap);
}

bool followPoints(AttributeMap& map, const MapEdit& edit)
{
    if (edit.points != nullptr && edit.points->unchanged()) {
        return true;
    }
    if (edit.pointSources != nullptr) {
        const Csr& sources = *edit.pointSources;
        const auto added = static_cast<long long>(edit.newPointCount) - static_cast<long long>(map.size());
        if (added == sources.rows()) {
            const bool useWeights =
                edit.pointWeights != nullptr && edit.pointWeights->size() == sources.entries.size();
            map.interpolateAppend(sources, useWeights ? edit.pointWeights : nullptr);
        } else if (auto* ints = dynamic_cast<TypedMap<int>*>(&map)) {
            ints->resizeFilled(static_cast<std::size_t>(edit.newPointCount), -1);
        } else {
            map.resize(static_cast<std::size_t>(edit.newPointCount));
        }
        return true;
    }
    if (edit.points != nullptr && edit.points->present()) {
        map.reorder(*edit.points);
        return true;
    }
    return false;
}

bool followCorners(AttributeMap& map, const MapEdit& edit)
{
    if (edit.corners != nullptr && (edit.corners->unchanged() || edit.corners->present())) {
        return followIndexMap(map, edit.corners);
    }
    return followIndexMap(map, edit.halfEdges);
}

} // namespace

void AttributeMaps::applyEdit(const MapEdit& edit)
{
    for (std::size_t i = maps_.size(); i-- > 0;) {
        AttributeMap& map = *maps_[i];
        bool keep = false;
        switch (map.kind()) {
        case ElementKind::Point:
            keep = followPoints(map, edit);
            break;
        case ElementKind::Face:
            keep = followIndexMap(map, edit.faces);
            break;
        case ElementKind::HalfEdge:
            keep = followIndexMap(map, edit.halfEdges);
            break;
        case ElementKind::Edge:
            keep = edit.edges != nullptr && edit.edges->unchanged();
            break;
        case ElementKind::Corner:
            keep = followCorners(map, edit);
            break;
        case ElementKind::UvPoint:
            keep = edit.uvSetsKept;
            break;
        case ElementKind::UvEdge:
            keep = edit.uvSetsKept && followIndexMap(map, edit.uvEdges);
            break;
        }
        if (!keep) {
            maps_.erase(maps_.begin() + static_cast<std::ptrdiff_t>(i));
        }
    }
}

namespace {

// Element count of a map kind that grows when a mesh is appended, or -1;
// half edge maps grow (by the corner count) only when the maps unite.
int growingCount(ElementKind kind, const ElementCounts& counts, bool unite)
{
    switch (kind) {
    case ElementKind::Point:
        return counts.points;
    case ElementKind::Face:
        return counts.faces;
    case ElementKind::Corner:
        return counts.corners;
    case ElementKind::Edge:
        return counts.edges;
    case ElementKind::HalfEdge:
        return unite ? counts.corners : -1;
    default:
        return -1;
    }
}

} // namespace

void AttributeMaps::appendFrom(const AttributeMaps& other, const ElementCounts& otherCounts,
                               const ElementCounts* ownCounts)
{
    const bool unite = ownCounts != nullptr;
    if (unite) {
        // Maps only `other` has join first, with default values for the own elements.
        for (const auto& map : other.maps_) {
            const int own = growingCount(map->kind(), *ownCounts, unite);
            if (own < 0 || find(map->name(), map->kind()) != nullptr) {
                continue;
            }
            std::unique_ptr<AttributeMap> added = map->clone();
            added->resize(0);
            added->resize(static_cast<std::size_t>(own));
            maps_.push_back(std::move(added));
        }
    }
    for (auto& map : maps_) {
        const AttributeMap* partner = other.find(map->name(), map->kind());
        if (partner != nullptr && map->appendFrom(*partner)) {
            continue;
        }
        const int grow = growingCount(map->kind(), otherCounts, unite);
        if (grow >= 0) {
            map->resize(map->size() + static_cast<std::size_t>(grow));
        }
    }
    removeKind(ElementKind::UvEdge);
}

} // namespace areko::remesher::mesh
