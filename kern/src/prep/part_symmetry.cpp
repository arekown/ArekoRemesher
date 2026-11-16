#include "part_symmetry.h"

#include <cstddef>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

PartSymmetry::Classes PartSymmetry::classify(const Mesh& mesh, SymmetryCache& cache)
{
    Classes classes;
    const ConnectedParts& parts = mesh.connectedParts();
    if (parts.partCount() == 0 || parts.partFaces().rowSize(0) == 0) {
        return classes;
    }
    const std::vector<int>& faceMirror = cache.faceMirror(mesh);
    int unpaired = 0;
    int group = -1; // -2: several groups
    for (const int f : parts.partFaces().row(0)) {
        const int m = static_cast<std::size_t>(f) < faceMirror.size() ? faceMirror[static_cast<std::size_t>(f)] : -1;
        if (m == -1) {
            ++unpaired;
        } else if (group != -2) {
            const int g = parts.partOfFace(m);
            if (g != group) {
                group = group == -1 ? g : -2;
            }
        }
    }
    if (group >= 0 && unpaired == 0) {
        classes.valid = true;
        classes.mirrorGroup = group;
    }
    return classes;
}

BitArray PartSymmetry::half(const Mesh& mesh, const SymmetryCache& cache, const Classes& classes)
{
    BitArray half(static_cast<std::size_t>(mesh.faceCount()));
    if (!classes.valid || classes.mirrorGroup != 0) {
        return half;
    }
    const ConnectedParts& parts = mesh.connectedParts();
    int zero = 0;
    int one = 0;
    const auto faces = parts.partFaces().row(0);
    for (const int f : faces) {
        const int side = cache.sideOfFace(mesh.topology(), f);
        if (side == 0) {
            half.set(static_cast<std::size_t>(f));
            ++zero;
        } else if (side == 1) {
            ++one;
        }
    }
    if (!(one == zero && static_cast<int>(faces.size()) == 2 * one)) {
        for (const int f : faces) {
            half.set(static_cast<std::size_t>(f), false);
        }
    }
    return half;
}

} // namespace areko::remesher::prep
