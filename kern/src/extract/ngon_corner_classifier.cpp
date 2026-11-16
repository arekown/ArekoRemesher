#include "ngon_corner_classifier.h"

#include "mesh/measures.h"
#include "prep/polygon_triangulator.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

namespace areko::remesher::extract {

using namespace areko::remesher::mesh;

namespace {

constexpr float bits(std::uint32_t pattern) { return std::bit_cast<float>(pattern); }

// Bands of the folded angle min(theta, 2 pi - theta).
constexpr float bandLow = bits(0x3f71463bu);       // 54 degrees
constexpr float bandHigh = bits(0x400cbe4cu);      // 126 degrees
constexpr float straightLimit = bits(0x3ef1463bu); // 27 degrees
constexpr float spikeLimit = bits(0x402ae714u);    // 153 degrees
// Hexagon with one reflex corner.
constexpr float hexReflex = bits(0x4071463bu);     // 216 degrees
constexpr float hexConvexMax = bits(0x4024deecu);  // 147.6 degrees
constexpr float hexMin = bits(0x3f10c3bdu);        // 32.4 degrees
constexpr float rightAngle = bits(0x3fc90fdbu);    // 90 degrees
constexpr float sumLow = bits(0x3f490fdbu);        // 45 degrees
constexpr float sumHigh = bits(0x4016cbe4u);       // 135 degrees
constexpr float wideStraight = bits(0x407b53d2u);  // 225 degrees
// 20 degrees in double, limit of |angle sum - 2 pi|.
constexpr double sumTolerance = 0.3490658503988659;
static_assert(std::bit_cast<std::uint64_t>(sumTolerance) == 0x3fd657184ae74487ull);

constexpr float initialAngle = -100.0f;
constexpr int unclear = 3;

float folded(float theta) { return theta > piF ? fullTurnF - theta : theta; }

} // namespace

void CornerClassifier::classify(const NgonPolygonList& list, const NgonNewPoints& points, NgonCornerData& data)
{
    const int old = static_cast<int>(data.angle.size());
    const auto corners = static_cast<std::size_t>(list.cornerCount());
    data.angle.resize(corners, initialAngle);
    data.flag.resize(corners, unclear);
    for (int p = 0; p < list.polygonCount(); ++p) {
        if (list.end(p) > old) {
            classifyPolygon(list, points, p, data);
        }
    }
}

void CornerClassifier::classifyPolygon(const NgonPolygonList& list, const NgonNewPoints& points, int p, NgonCornerData& data)
{
    const int s = list.begin(p);
    const int e = list.end(p);
    const int n = e - s;
    std::vector<Vec3f> q(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        q[static_cast<std::size_t>(i)] = points.position(list.vertex(s + i));
    }
    const Vec3f normal = prep::PolygonTriangulator::planeNormal(q);
    auto& angle = data.angle;
    auto& flag = data.flag;
    auto at = [](int i) { return static_cast<std::size_t>(i); };
    int turn = 0;     // +1 per corner, -1 per concave corner, +2 per spike
    int unclearCount = 0;
    int spikes = 0;
    int concave = 0;
    for (int i = s; i < e; ++i) {
        const Vec3f here = q[static_cast<std::size_t>(i - s)];
        const Vec3f a = here - points.position(list.vertex(list.prev(p, i)));
        const Vec3f b = points.position(list.vertex(list.next(p, i))) - here;
        const float theta = static_cast<float>(signedVectorAngleAsDouble(a, b, cross(a, b), normal));
        const float phi = folded(theta);
        if (bandLow <= phi && phi <= bandHigh) {
            if (theta <= piF) {
                flag[at(i)] = 1;
                ++turn;
            } else {
                flag[at(i)] = -1;
                --turn;
                ++concave;
            }
        } else if (!(phi <= straightLimit)) {
            if (!(phi >= spikeLimit)) {
                ++unclearCount;
            } else {
                flag[at(i)] = 2;
                ++spikes;
                turn += 2;
            }
        } else {
            flag[at(i)] = 0;
        }
        angle[at(i)] = theta;
    }

    // Hexagon with one clearly reflex corner: the others are corners.
    if (n == 6 && unclearCount != 0) {
        float maxAngle = 0.0f;
        int maxAt = -1;
        float minAngle = std::numeric_limits<float>::max();
        float maxConvex = 0.0f;
        for (int i = s; i < e; ++i) {
            const float theta = angle[at(i)];
            if (theta > maxAngle) {
                maxAngle = theta;
                maxAt = i;
            }
            if (theta < minAngle) {
                minAngle = theta;
            }
            if (piF > theta && theta > maxConvex) {
                maxConvex = theta;
            }
        }
        if (maxAngle > hexReflex && maxConvex < hexConvexMax && minAngle > hexMin) {
            for (int i = s; i < e; ++i) {
                flag[at(i)] = 1;
            }
            flag[at(maxAt)] = -1;
        }
    }

    // Three corners: the unclear corner nearest to a right angle becomes the fourth.
    if (turn == 3) {
        float best = std::numeric_limits<float>::max();
        int bestAt = -1;
        for (int i = s; i < e; ++i) {
            const float theta = angle[at(i)];
            const float phi = folded(theta);
            if (flag[at(i)] == unclear && straightLimit <= phi && phi <= spikeLimit && theta >= 0.0f) {
                const float off = std::fabs(phi - rightAngle);
                if (off < best) {
                    best = off;
                    bestAt = i;
                }
            }
        }
        if (bestAt >= 0) {
            flag[at(bestAt)] = 1;
            --unclearCount;
            turn = 4;
        }
    }

    auto unclearToStraight = [&] {
        for (int i = s; i < e; ++i) {
            if (flag[at(i)] == unclear && folded(angle[at(i)]) <= bandLow) {
                flag[at(i)] = 0;
            }
        }
    };
    if (concave != 0 || spikes != 0) {
        if ((turn & 3) == 0 && unclearCount > 0) {
            unclearToStraight();
        }
    } else if (unclearCount < 5) {
        int counted = 0;
        float sum = 0.0f;
        for (int i = s; i < e; ++i) {
            const float theta = angle[at(i)];
            if (flag[at(i)] == unclear && !(theta < sumLow) && !(theta > sumHigh)) {
                sum += theta;
                ++counted;
            }
        }
        if (counted == unclearCount && wide(std::fabs(sum - fullTurnF)) < sumTolerance) {
            for (int i = s; i < e; ++i) {
                if (flag[at(i)] == unclear) {
                    flag[at(i)] = 1;
                }
            }
        } else if ((turn & 3) == 0 && unclearCount > 0) {
            unclearToStraight();
        }
    } else if ((turn & 3) == 0) {
        unclearToStraight();
    }
}

FlagCounts CornerClassifier::count(const NgonCornerData& data, int first, int last)
{
    FlagCounts k;
    for (int i = first; i < last; ++i) {
        switch (data.flag[static_cast<std::size_t>(i)]) {
        case 1:
            ++k.corners;
            break;
        case 0:
            ++k.straights;
            k.lastStraight = i;
            break;
        case -1:
            ++k.concave;
            break;
        case 2:
            ++k.spikes;
            break;
        case unclear:
            ++k.unclear;
            break;
        default:
            break;
        }
    }
    return k;
}

void CornerClassifier::reclassify(NgonCornerData& data, FlagCounts& k, int first, int last)
{
    auto& flag = data.flag;
    const auto& angle = data.angle;
    if (k.corners + k.unclear == 4 && k.concave == 0) {
        for (int i = first; i < last; ++i) {
            if (flag[static_cast<std::size_t>(i)] == unclear && angle[static_cast<std::size_t>(i)] < spikeLimit) {
                flag[static_cast<std::size_t>(i)] = 1;
                --k.unclear;
                ++k.corners;
            }
        }
    }
    if (k.corners == 4 && k.concave == 0) {
        for (int i = first; i < last; ++i) {
            const float theta = angle[static_cast<std::size_t>(i)];
            if (flag[static_cast<std::size_t>(i)] == unclear && sumHigh < theta && theta < wideStraight) {
                flag[static_cast<std::size_t>(i)] = 0;
                --k.unclear;
                ++k.straights;
                k.lastStraight = i;
            }
        }
    }
}

} // namespace areko::remesher::extract
