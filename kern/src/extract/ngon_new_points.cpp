#include "ngon_new_points.h"

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

NgonNewPoints::NgonNewPoints(const Mesh& mesh) : mesh_(mesh), base_(mesh.pointCount()) {}

Vec3f NgonNewPoints::position(int id) const
{
    if (id < base_) {
        return mesh_.points()[static_cast<std::size_t>(id)];
    }
    return positions_[static_cast<std::size_t>(id - base_)];
}

int NgonNewPoints::add(Vec3f position)
{
    positions_.push_back(position);
    return base_ + count() - 1;
}

void NgonNewPoints::addParent(int id, float weight)
{
    if (id < base_) {
        parents_.entries.push_back(id);
        weights_.push_back(weight);
        return;
    }
    const int row = id - base_;
    for (int r = parents_.starts[static_cast<std::size_t>(row)]; r < parents_.starts[static_cast<std::size_t>(row) + 1];
         ++r) {
        const int parent = parents_.entries[static_cast<std::size_t>(r)];
        const float product = weights_[static_cast<std::size_t>(r)] * weight;
        parents_.entries.push_back(parent);
        weights_.push_back(product);
    }
}

void NgonNewPoints::endParents()
{
    parents_.starts.push_back(static_cast<int>(parents_.entries.size()));
}

void NgonNewPoints::move(int id, Vec3f position)
{
    positions_[static_cast<std::size_t>(id - base_)] = position;
}

void NgonNewPoints::reset()
{
    base_ = mesh_.pointCount();
    positions_.clear();
    parents_.reset();
    // The weights go with their parents, so rows of later rounds read their own weights.
    weights_.clear();
}

} // namespace areko::remesher::extract
