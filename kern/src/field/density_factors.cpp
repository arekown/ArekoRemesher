#include "density_factors.h"

namespace areko::remesher::field {

using namespace areko::remesher::mesh;

namespace {

// Range of the factor: 4 at full red, 1 / 4 at full cyan.
constexpr float finest = 4.0f;
constexpr float coarsest = 0.25f;

float colorFactor(std::uint8_t r, std::uint8_t g, std::uint8_t b)
{
    const float fine = finest - 1.0f;
    const float coarse = 1.0f / coarsest - 1.0f;
    return ((1.0f - static_cast<float>(static_cast<int>(g) + static_cast<int>(b)) / 510.0f) * fine + 1.0f) *
           (1.0f / ((1.0f - static_cast<float>(r) / 255.0f) * coarse + 1.0f));
}

template <class Color>
std::vector<float> colorFactors(const std::vector<Color>& colors)
{
    std::vector<float> out(colors.size());
    for (std::size_t i = 0; i < colors.size(); ++i) {
        out[i] = colorFactor(colors[i].r, colors[i].g, colors[i].b);
    }
    return out;
}

} // namespace

float DensityFactors::fromDensity(float density)
{
    if (density == 1.0f) {
        return 1.0f;
    }
    // How far the value lies on the way from 1 to its limit; not above 0
    // (also NaN) gives 1, not below 1 the limit.
    const bool finer = 1.0f < density;
    const float share = finer ? (density - 1.0f) / (finest - 1.0f) : (1.0f / density - 1.0f) / (1.0f / coarsest - 1.0f);
    if (!(0.0f < share)) {
        return 1.0f;
    }
    if (!(share < 1.0f)) {
        return finer ? finest : coarsest;
    }
    return density;
}

std::optional<std::vector<float>> DensityFactors::fromColors(const Mesh& mesh)
{
    const AttributeMaps& maps = mesh.maps();
    const auto points = static_cast<std::size_t>(mesh.pointCount());
    if (const auto* rgb = maps.findAs<TypedMap<Color3ub>>(map_names::colors, ElementKind::Point)) {
        return rgb->size() == points ? std::optional(colorFactors(rgb->values())) : std::nullopt;
    }
    if (const auto* rgba = maps.findAs<TypedMap<Color4ub>>(map_names::colors, ElementKind::Point)) {
        return rgba->size() == points ? std::optional(colorFactors(rgba->values())) : std::nullopt;
    }
    return std::nullopt;
}

std::optional<std::vector<float>> DensityFactors::fromValues(const Mesh& mesh)
{
    const auto* values = mesh.maps().findAs<TypedMap<float>>(map_names::density, ElementKind::Point);
    if (values == nullptr || values->size() != static_cast<std::size_t>(mesh.pointCount())) {
        return std::nullopt;
    }
    std::vector<float> out(values->size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = fromDensity(values->values()[i]);
    }
    return out;
}

void DensityFactors::apply(const Mesh& mesh, const std::vector<float>& factors, std::vector<float>& sizeU,
                           std::vector<float>& sizeV)
{
    const Topology& topology = mesh.topology();
    for (std::size_t t = 0; t < sizeU.size(); ++t) {
        const int c = 3 * static_cast<int>(t);
        const float mean = ((factors[static_cast<std::size_t>(topology.vertex(c))] + 0.0f) +
                            factors[static_cast<std::size_t>(topology.vertex(c + 1))] +
                            factors[static_cast<std::size_t>(topology.vertex(c + 2))]) /
                           3.0f;
        sizeU[t] = sizeU[t] / mean;
        sizeV[t] = sizeV[t] / mean;
    }
}

} // namespace areko::remesher::field
