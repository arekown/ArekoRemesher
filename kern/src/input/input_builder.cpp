#include "input_builder.h"

#include "field/density_factors.h"
#include "mesh/bit_array.h"
#include "mesh/rank_index.h"

#include <bit>
#include <cstdint>
#include <numeric>
#include <utility>

namespace areko::remesher::input {

using namespace areko::remesher::mesh;

namespace {

float positiveZero(float v) { return v == 0.0f ? 0.0f : v; }

void buildTopology(const MeshInput& input, Mesh& out)
{
    std::vector<int> faceStart;
    faceStart.reserve(static_cast<std::size_t>(input.faceCount()) + 1);
    faceStart.push_back(0);
    for (const std::int32_t size : input.faceSizes()) {
        faceStart.push_back(faceStart.back() + size);
    }
    std::vector<int> corners(input.cornerPoints().begin(), input.cornerPoints().end());
    out.editTopology().swapPolygons(faceStart, corners);
    out.topologyChanged();
}

void buildPoints(const MeshInput& input, Mesh& out)
{
    std::vector<Vec3f> points(static_cast<std::size_t>(input.pointCount()));
    const auto& raw = input.points();
    for (std::size_t i = 0; i < points.size(); ++i) {
        points[i] = {raw[3 * i], raw[3 * i + 1], raw[3 * i + 2]};
    }
    out.swapPoints(points);
}

void buildUvSets(const MeshInput& input, Mesh& out)
{
    const auto& corners = input.cornerPoints();
    for (const auto& uv : input.uvSets()) {
        using Key = std::pair<std::uint64_t, int>;
        std::vector<Key> keys(corners.size());
        for (std::size_t c = 0; c < corners.size(); ++c) {
            const auto u = std::bit_cast<std::uint32_t>(positiveZero(uv.values[2 * c]));
            const auto v = std::bit_cast<std::uint32_t>(positiveZero(uv.values[2 * c + 1]));
            keys[c] = {(static_cast<std::uint64_t>(v) << 32) | u, corners[c]};
        }
        std::vector<int> representative;
        UvSet set;
        set.name = uv.name;
        set.cornerUv = rankByKey(keys, representative);
        set.coords.reserve(representative.size());
        for (const int c : representative) {
            set.coords.push_back({positiveZero(uv.values[2 * static_cast<std::size_t>(c)]),
                                  positiveZero(uv.values[2 * static_cast<std::size_t>(c) + 1])});
        }
        out.uvSets().push_back(std::move(set));
    }
}

void buildNormals(const MeshInput& input, Mesh& out)
{
    if (!input.cornerNormals()) {
        return;
    }
    const auto& raw = *input.cornerNormals();
    const std::size_t corners = static_cast<std::size_t>(input.cornerCount());
    std::vector<int> representative;
    CustomNormals normals;
    if (input.cornerNormalsShared()) {
        using Key = std::array<std::uint32_t, 6>;
        std::vector<Key> keys(corners);
        for (std::size_t c = 0; c < corners; ++c) {
            for (int k = 0; k < 3; ++k) {
                const double value = static_cast<double>(positiveZero(raw[3 * c + static_cast<std::size_t>(k)]));
                const auto bits = std::bit_cast<std::uint64_t>(value);
                keys[c][static_cast<std::size_t>(2 * k)] = static_cast<std::uint32_t>(bits);
                keys[c][static_cast<std::size_t>(2 * k + 1)] = static_cast<std::uint32_t>(bits >> 32);
            }
        }
        normals.cornerNormal = rankByKey(keys, representative);
    } else {
        representative.resize(corners);
        std::iota(representative.begin(), representative.end(), 0);
        normals.cornerNormal = representative;
    }
    normals.normals.reserve(representative.size());
    for (const int c : representative) {
        const std::size_t i = 3 * static_cast<std::size_t>(c);
        normals.normals.push_back({positiveZero(raw[i]), positiveZero(raw[i + 1]), positiveZero(raw[i + 2])});
    }
    normals.state = 1;
    out.customNormals() = std::move(normals);
}

void addFaceMap(const MeshInput& input, arm_face_values kind, std::string_view name, Mesh& out)
{
    const auto& values = input.faceValues(kind);
    if (!values) {
        return;
    }
    if (kind == ARM_FACE_SMOOTHING_GROUPS) {
        std::vector<std::uint32_t> bits(values->size());
        for (std::size_t i = 0; i < bits.size(); ++i) {
            bits[i] = static_cast<std::uint32_t>((*values)[i]);
        }
        out.maps().set(std::make_unique<TypedMap<std::uint32_t>>(std::string(name), ElementKind::Face, std::move(bits)));
    } else {
        out.maps().set(std::make_unique<TypedMap<int>>(std::string(name), ElementKind::Face,
                                                       std::vector<int>(values->begin(), values->end())));
    }
}

void addSharpEdges(const MeshInput& input, Mesh& out)
{
    const auto& pairs = input.sharpEdges();
    if (!pairs || pairs->size() < 2) {
        return;
    }
    const Adjacency& adjacency = out.adjacency();
    BitArray bits(static_cast<std::size_t>(adjacency.edgeCount()));
    for (std::size_t i = 0; i + 1 < pairs->size(); i += 2) {
        const int e = adjacency.edgeBetween((*pairs)[i], (*pairs)[i + 1]);
        if (e >= 0) {
            bits.set(static_cast<std::size_t>(e));
        }
    }
    out.maps().set(std::make_unique<BoolMap>(std::string(map_names::sharpEdges), ElementKind::Edge, std::move(bits)));
}

std::uint8_t colorChannel(float value)
{
    double v = value;
    if (!(v > 0.0)) {
        v = 0.0;
    } else if (v > 1.0) {
        v = 1.0;
    }
    return static_cast<std::uint8_t>(static_cast<int>(v * 255.0));
}

void addColors(const MeshInput& input, const arm_settings& settings, Mesh& out)
{
    if (!input.cornerColors() || !(settings.color_density || settings.keep_colors)) {
        return;
    }
    const bool opaque = settings.color_density && !settings.keep_colors;
    const auto& raw = *input.cornerColors();
    std::vector<Color4ub> colors(static_cast<std::size_t>(input.pointCount()), Color4ub{0, 0, 0, 255});
    const auto& corners = input.cornerPoints();
    for (std::size_t c = 0; c < corners.size(); ++c) {
        const float* rgba = raw.data() + 4 * c;
        colors[static_cast<std::size_t>(corners[c])] = {colorChannel(rgba[0]), colorChannel(rgba[1]),
                                                        colorChannel(rgba[2]),
                                                        opaque ? std::uint8_t{255} : colorChannel(rgba[3])};
    }
    out.maps().set(std::make_unique<TypedMap<Color4ub>>(std::string(map_names::colors), ElementKind::Point, std::move(colors)));
}

void addDensity(const MeshInput& input, Mesh& out)
{
    if (!input.pointDensity()) {
        return;
    }
    // The map holds the factors the given values stand for, so that mixing
    // them later mixes what scales the sizes.
    std::vector<float> factors = *input.pointDensity();
    for (float& value : factors) {
        value = field::DensityFactors::fromDensity(value);
    }
    out.maps().set(std::make_unique<TypedMap<float>>(std::string(map_names::density), ElementKind::Point,
                                                     std::move(factors)));
}

} // namespace

Mesh InputBuilder::build(const MeshInput& input, const arm_settings& settings)
{
    Mesh out;
    buildPoints(input, out);
    buildTopology(input, out);
    buildUvSets(input, out);
    buildNormals(input, out);
    addFaceMap(input, ARM_FACE_GROUP, map_names::faceGroups, out);
    addFaceMap(input, ARM_FACE_MATERIAL, map_names::materials, out);
    addFaceMap(input, ARM_FACE_SMOOTHING_GROUPS, map_names::smoothingGroups, out);
    addSharpEdges(input, out);
    addColors(input, settings, out);
    addDensity(input, out);
    return out;
}

} // namespace areko::remesher::input
