#include "size_smoothing.h"

namespace areko::remesher::field {

void SizeSmoothing::mean(std::vector<float>& sizeU, std::vector<float>& sizeV) const
{
    std::vector<float> newU(sizeU.size());
    std::vector<float> newV(sizeV.size());
    for (std::size_t t = 0; t < sizeU.size(); ++t) {
        float sumU = sizeU[t];
        float sumV = sizeV[t];
        float countU = 1.0f;
        float countV = 1.0f;
        for (int k = 0; k < 3; ++k) {
            const auto across = neighbours_.across(3 * static_cast<int>(t) + k);
            if (!across) {
                continue;
            }
            const auto n = static_cast<std::size_t>(across->triangle);
            if (across->take.u) {
                countU += 1.0f;
                sumU += across->odd ? sizeV[n] : sizeU[n];
            }
            if (across->take.v) {
                countV += 1.0f;
                sumV += across->odd ? sizeU[n] : sizeV[n];
            }
        }
        newU[t] = sumU / countU;
        newV[t] = sumV / countV;
    }
    sizeU.swap(newU);
    sizeV.swap(newV);
}

void SizeSmoothing::minimum(int rounds, std::vector<float>& sizeU, std::vector<float>& sizeV) const
{
    std::vector<float> newU(sizeU.size());
    std::vector<float> newV(sizeV.size());
    for (int round = 0; round < rounds; ++round) {
        for (std::size_t t = 0; t < sizeU.size(); ++t) {
            float minU = sizeU[t];
            float minV = sizeV[t];
            for (int k = 0; k < 3; ++k) {
                const auto across = neighbours_.across(3 * static_cast<int>(t) + k);
                if (!across) {
                    continue;
                }
                const auto n = static_cast<std::size_t>(across->triangle);
                const float u = across->odd ? sizeV[n] : sizeU[n];
                const float v = across->odd ? sizeU[n] : sizeV[n];
                if (across->take.u && u < minU) {
                    minU = u;
                }
                if (across->take.v && v < minV) {
                    minV = v;
                }
            }
            newU[t] = minU;
            newV[t] = minV;
        }
        sizeU.swap(newU);
        sizeV.swap(newV);
    }
}

} // namespace areko::remesher::field
