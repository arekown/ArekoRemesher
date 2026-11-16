#include "part.h"

#include "prep/freeze_marks.h"
#include "prep/symmetry_detector.h"

namespace areko::remesher::parts {

using namespace areko::remesher::mesh;
using prep::SymmetryCache;

SymmetryCache* PartMirror::symmetric(const Mesh& mesh, bool refresh)
{
    if (!current_ && refresh) {
        current_ = true;
        std::optional<SymmetryCache> found = prep::SymmetryDetector::run(mesh, stored());
        if (found) {
            cache_ = std::move(found);
            symmetric_ = true;
        } else {
            symmetric_ = false;
            if (cache_) {
                // Without a find only the point pairs of the same point count remain.
                const bool keep = static_cast<int>(cache_->pointMirror().size()) == mesh.pointCount();
                cache_ = SymmetryCache::fromPoints(keep ? cache_->pointMirror() : std::vector<int>{},
                                                   keep ? cache_->pointSide() : BitArray{}, cache_->lockedPlane());
            }
        }
    }
    return symmetric_ && cache_ ? &*cache_ : nullptr;
}

SymmetryCache* PartMirror::evaluate(const Mesh& mesh, bool use)
{
    SymmetryCache* pairs = use ? symmetric(mesh, true) : nullptr;
    classes = pairs != nullptr ? prep::PartSymmetry::classify(mesh, *pairs) : prep::PartSymmetry::Classes{};
    return pairs;
}

const std::vector<FreezeMark>* Part::frozenBorderMarks() const
{
    const auto* marks = bordersFrozen ? prep::FreezeMarks::find(mesh) : nullptr;
    return marks != nullptr ? &marks->values() : nullptr;
}

} // namespace areko::remesher::parts
