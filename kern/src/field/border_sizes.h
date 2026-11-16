#pragma once

#include "size_field.h"

#include "parts/part.h"

#include <vector>

namespace areko::remesher::field {

// Sizes at frozen borders of one part. Between two consecutive kept border
// points the quads are to be as long as the chord between them, so the
// triangles along that stretch of border get the chord as their size along
// the border. From there the size passes into the part: unchanged up to two
// target edge lengths from the border, then fading linearly into the target
// edge length at six.
class BorderSizes {
public:
    // The part has a cut graph and its target edge length; `marks` are its
    // freeze marks. Sizes given by the border are the chords times the border
    // density factor, i.e. over `borderDensityInverse` unless that is 0.
    BorderSizes(const parts::Part& part, const std::vector<mesh::FreezeMark>& marks, float borderDensityInverse);

    // Writes the sizes given by the border, then the sizes of the seam beside
    // it. The seam leaves out the triangles of the mirrored half.
    void apply(std::vector<float>& sizeU, std::vector<float>& sizeV) const;

    // Triangle t lies in the seam: it is not in the mirrored half, has a path
    // to the border and that path is shorter than six target edge lengths.
    [[nodiscard]] bool inSeam(std::size_t t) const;
    // Whether the two directions change places between a triangle of the seam
    // and the border triangle its path leads to: the quarter turns left
    // across the edges of the path add up to an odd number.
    [[nodiscard]] bool swapped(std::size_t t) const { return swapped_.test(t); }

    [[nodiscard]] const SizeField::BorderData& data() const noexcept { return data_; }
    [[nodiscard]] SizeField::BorderData release() noexcept { return std::move(data_); }

private:
    // One size given by the border: a triangle, the direction and the value.
    struct Given {
        int face = 0;
        bool alongU = false;
        float size = 0.0f;
    };

    void collectGiven(const std::vector<mesh::FreezeMark>& marks, float borderDensityInverse);
    void measureDistances();
    void fix();
    void markSwapped();
    [[nodiscard]] bool pathSwaps(int t, int borderFace) const;

    const parts::Part& part_;
    const mesh::Adjacency& adjacency_;
    const CutGraph& cuts_;
    float edgeLength_;
    std::vector<Given> given_;
    SizeField::BorderData data_;
    mesh::BitArray swapped_;
};

} // namespace areko::remesher::field
