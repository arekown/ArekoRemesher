#include "vertex_color_transfer.h"

#include "source_weights.h"

#include "mesh/attribute_maps.h"

#include <cstdint>
#include <memory>
#include <string>

namespace areko::remesher::extract {

namespace {

// One channel: a sum above 255 is cut to 255 (a sum without number stays),
// one below 0 gives 0, otherwise its truncated value (without number 0).
std::uint8_t channelByte(float sum)
{
    const float capped = 255.0f < sum ? 255.0f : sum;
    if (0.0f > capped) {
        return 0;
    }
    return static_cast<std::uint8_t>(static_cast<std::uint32_t>(mesh::truncToInt(capped)) & 0xffu);
}

template <class Color, int Channels>
std::vector<Color> mixed(const std::vector<Color>& source, const SourceWeights& weights)
{
    std::vector<Color> out(static_cast<std::size_t>(weights.rows.rows()), Color{});
    for (int i = 0; i < weights.rows.rows(); ++i) {
        float sum[Channels] = {};
        for (int k = weights.rows.starts[static_cast<std::size_t>(i)];
             k < weights.rows.starts[static_cast<std::size_t>(i) + 1]; ++k) {
            const Color& c = source[static_cast<std::size_t>(weights.rows.entries[static_cast<std::size_t>(k)])];
            const float w = weights.weights[static_cast<std::size_t>(k)];
            for (int ch = 0; ch < Channels; ++ch) {
                sum[ch] = sum[ch] + static_cast<float>((&c.r)[ch]) * w;
            }
        }
        Color& target = out[static_cast<std::size_t>(i)];
        for (int ch = 0; ch < Channels; ++ch) {
            (&target.r)[ch] = channelByte(sum[ch]);
        }
    }
    return out;
}

template <class Color, int Channels>
void transfer(mesh::Mesh& output, const std::vector<Color>& source, std::span<const PointSource> sources,
              const mesh::Mesh& triangles)
{
    const SourceWeights weights = SourceWeights::of(triangles, output.points(), sources);
    output.maps().set(std::make_unique<mesh::TypedMap<Color>>(std::string(mesh::map_names::colors),
                                                              mesh::ElementKind::Point,
                                                              mixed<Color, Channels>(source, weights)));
}

} // namespace

void VertexColorTransfer::run(mesh::Mesh& output, std::span<const PointSource> sources, const mesh::Mesh& triangles)
{
    using mesh::ElementKind;
    using mesh::TypedMap;
    const mesh::AttributeMaps& maps = triangles.maps();
    if (const auto* rgb = maps.findAs<TypedMap<mesh::Color3ub>>(mesh::map_names::colors, ElementKind::Point)) {
        transfer<mesh::Color3ub, 3>(output, rgb->values(), sources, triangles);
    } else if (const auto* rgba = maps.findAs<TypedMap<mesh::Color4ub>>(mesh::map_names::colors, ElementKind::Point)) {
        transfer<mesh::Color4ub, 4>(output, rgba->values(), sources, triangles);
    }
}

} // namespace areko::remesher::extract
