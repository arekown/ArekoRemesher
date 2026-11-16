#include "feature_chains.h"

#include "mesh/edge_chains.h"
#include "mesh/measures.h"
#include "mesh/vertex_fan.h"

#include <algorithm>
#include <bit>
#include <cfloat>
#include <cmath>
#include <numbers>
#include <utility>

namespace areko::remesher::prep {

using namespace areko::remesher::mesh;

namespace {

const float rampDivisor = std::bit_cast<float>(0x3eb2b8c3u); // 20 degrees, one step above degF(20)
const float jumpRampWidth = std::bit_cast<float>(0x3ee85697u); // 26 degrees, one step above degF(26)
const float rivalFactor = std::bit_cast<float>(0x3f866666u);   // 1.05

constexpr int maxRounds = 10;
constexpr int lateRounds = 3;

// Measures of one chain (edges e0 .. e(n-1)); float unless noted.
struct Measures {
    float length = 0.0f;        // L
    float bendMin = 0.0f;       // smallest bend angle at an inner point
    float bendMax = 0.0f;
    float jumpMin = 0.0f;       // smallest |edge angle difference| of neighbouring edges
    float jumpMax = 0.0f;
    float rivalMax = 0.0f;      // largest edge angle of a free edge at an inner point
    float ratioMax = 0.0f;      // rival / mean edge angle of the two chain edges
    float ratioMean = 0.0f;     // sum of the ratios divided by the edge count
    float smoothness = 0.0f;    // length weighted point weight
    float weightMin = 0.0f;
    float weightMax = 0.0f;
    float angleMax = 0.0f;
    float angleMin = 0.0f;
    float angleMean = 0.0f;
    float angleWeighted = 0.0f; // length weighted edge angle
    int preClass = 0;           // 0, 1 or 3
    int decision = 0;           // 0, 1 reject, 2 soft line, 3 crease line
    bool classOne = false;      // crease line of class 1 instead of 2
    bool rescue = false;        // keep the chain aside when it is rejected
};

enum class Next { Reject, Late, Assign, Soft };

struct FreeAngles {
    float max = 0.0f;
    float sum = 0.0f;
};

// Measures of edge chains on the current geometry. Free edges at a point are
// those that are neither candidates nor protected.
class ChainMeasure {
public:
    ChainMeasure(const Mesh& mesh, const BitArray& candidates, const BitArray* protectedEdges)
        : adjacency_(mesh.adjacency()), points_(mesh.points()), angles_(mesh.dihedralAngles()),
          lengths_(edgeLengths(mesh)), candidates_(candidates), protected_(protectedEdges)
    {
    }

    [[nodiscard]] const Adjacency& adjacency() const noexcept { return adjacency_; }
    [[nodiscard]] float angle(int e) const { return angles_[static_cast<std::size_t>(e)]; }
    [[nodiscard]] Vec3f point(int v) const { return points_[static_cast<std::size_t>(v)]; }
    [[nodiscard]] bool isProtected(int e) const
    {
        return protected_ != nullptr && protected_->testBounded(static_cast<std::size_t>(e));
    }
    [[nodiscard]] bool isCandidate(int e) const { return candidates_.test(static_cast<std::size_t>(e)); }

    // Edge angles of the free edges at v.
    [[nodiscard]] FreeAngles freeAngles(int v) const
    {
        FreeAngles out;
        for (const int f : adjacency_.pointEdges(v)) {
            if (isCandidate(f) || isProtected(f)) {
                continue;
            }
            const float t = angle(f);
            out.sum += t;
            if (t > out.max) {
                out.max = t;
            }
        }
        return out;
    }

    // The measures of a chain with at least one edge; point i + 1 of `chain`
    // joins edges i and i + 1.
    [[nodiscard]] Measures measure(std::span<const int> edges, std::span<const int> chain) const
    {
        const std::size_t n = edges.size();
        const float count = static_cast<float>(n);
        const float t0 = angle(edges[0]);
        Measures m;
        m.length = edgeLength(edges[0]);
        m.bendMin = FLT_MAX;
        m.jumpMin = FLT_MAX;
        m.weightMin = FLT_MAX;
        m.angleMax = t0;
        m.angleMin = t0;
        m.angleMean = t0;
        m.angleWeighted = m.length * t0;
        if (n < 2) {
            m.angleMean = m.angleMean / count;
            m.angleWeighted = m.angleWeighted / m.length;
            m.ratioMean = m.ratioMean / count;
            m.bendMin = 0.0f;
            m.weightMin = 0.0f;
            m.jumpMin = 0.0f;
            return m;
        }
        float lengthSum = 0.0f;
        for (std::size_t i = 0; i + 1 < n; ++i) {
            const int ea = edges[i];
            const int eb = edges[i + 1];
            // The common end point of two consecutive chain edges.
            const int v = chain[i + 1];
            const Vec3f p = point(v);
            const float bend = vectorAngle(point(adjacency_.otherEnd(ea, v)) - p, p - point(adjacency_.otherEnd(eb, v)));
            const float ta = angle(ea);
            const float tb = angle(eb);
            const float jump = std::fabs(ta - tb);
            if (tb > m.angleMax) {
                m.angleMax = tb;
            }
            if (tb < m.angleMin) {
                m.angleMin = tb;
            }
            const float sumBefore = m.angleMean;
            m.angleMean = m.angleMean + tb;
            m.angleWeighted += edgeLength(eb) * tb;
            m.length += edgeLength(eb);
            if (bend > m.bendMax) {
                m.bendMax = bend;
            }
            if (bend < m.bendMin) {
                m.bendMin = bend;
            }
            if (jump > m.jumpMax) {
                m.jumpMax = jump;
            }
            if (jump < m.jumpMin) {
                m.jumpMin = jump;
            }
            const float rival = freeAngles(v).max;
            const float mid = (tb + ta) * 0.5f;
            const float ratio = rival / mid;
            if (rival > m.rivalMax) {
                m.rivalMax = rival;
            }
            if (ratio > m.ratioMax) {
                m.ratioMax = ratio;
            }
            m.ratioMean += ratio;
            const float w = pointWeight(bend, jump, ratio, mid, rival, sumBefore + tb);
            const float s = edgeLength(ea) + edgeLength(eb);
            lengthSum += s;
            m.smoothness += s * w;
            if (w > m.weightMax) {
                m.weightMax = w;
            }
            if (w < m.weightMin) {
                m.weightMin = w;
            }
        }
        if (lengthSum != 0.0f) {
            m.smoothness = m.smoothness / lengthSum;
        }
        m.angleMean = m.angleMean / count;
        m.angleWeighted = m.angleWeighted / m.length;
        m.ratioMean = m.ratioMean / count;
        return m;
    }

private:
    [[nodiscard]] float edgeLength(int e) const { return lengths_[static_cast<std::size_t>(e)]; }

    // Weight of an inner chain point: 1 smooth continuation, 0 break. A break
    // only lowers the smoothness of the chain; the chain is not split there,
    // since splitting at kinks was measured to give more irregular points and
    // triangles.
    [[nodiscard]] static float pointWeight(float bend, float jump, float ratio, float mid, float rival, float angleSum)
    {
        const double k = wide(bend);
        const double d = wide(jump);
        auto byMix = [&] { return wide((2.0f * jump + bend) / 3.0f) < degD(10) ? 1.0f : 0.5f; };
        if ((k < degD(8) && d < degD(5)) || (ratio < 2.5f && k < 0.7 * wide(mid))) {
            return 1.0f;
        }
        if (!(d < degD(5))) {
            if (k < degD(25) && d < degD(20)) {
                return byMix();
            }
            if (k >= degD(36) || d >= degD(12) || wide(angleSum) <= degD(60)) {
                return 0.0f;
            }
            if (wide(rival) < degD(24)) {
                return byMix();
            }
            return 0.0f;
        }
        if (bend < degF(10)) {
            return 1.0f;
        }
        if (k < degD(32)) {
            return 1.0f + ((bend - degF(10)) * -0.8f) / rampDivisor;
        }
        return 0.0f;
    }

    const Adjacency& adjacency_;
    const std::vector<Vec3f>& points_;
    const std::vector<float>& angles_;
    std::vector<float> lengths_;
    const BitArray& candidates_;
    const BitArray* protected_;
};

class Classifier {
public:
    Classifier(const Mesh& mesh, const FeatureChains::Options& options)
        : candidates_(static_cast<std::size_t>(mesh.adjacency().edgeCount())),
          soft_(static_cast<std::size_t>(mesh.adjacency().edgeCount())),
          measure_(mesh, candidates_, options.protectedEdges), adjacency_(measure_.adjacency()),
          creaseAngle_(options.creaseAngle), h_(options.lengths.reference), hQuad_(options.lengths.quad)
    {
        for (int e = 0; e < adjacency_.edgeCount(); ++e) {
            if (!isProtected(e) && !(creaseAngle_ > measure_.angle(e))) {
                candidates_.set(static_cast<std::size_t>(e));
            }
        }
    }

    FeatureChains::Result run()
    {
        FeatureChains::Result result;
        for (int round = 0;; ++round) {
            chainEdges_.reset();
            chains_ = EdgeChains::build(adjacency_, candidates_, &chainEdges_);
            classes_.assign(static_cast<std::size_t>(chains_.rows()), 0);
            if (round == maxRounds) {
                break;
            }
            measures_.clear();
            for (int j = 0; j < chains_.rows(); ++j) {
                measures_.push_back(measure_.measure(chainEdges_.row(j), chains_.row(j)));
            }
            buildNodes();
            rejected_ = false;
            for (int j = 0; j < chains_.rows(); ++j) {
                firstPass(j);
            }
            for (int j = 0; j < chains_.rows(); ++j) {
                secondPass(j, round);
            }
            if (!rejected_) {
                break;
            }
        }
        result.points = std::move(chains_);
        result.edges = std::move(chainEdges_);
        result.classes = std::move(classes_);
        for (int r = 0; r < rescuedPoints_.rows(); ++r) {
            append(result.points, rescuedPoints_.row(r));
            append(result.edges, rescuedEdges_.row(r));
            result.classes.push_back(3);
        }
        return result;
    }

private:
    static void append(Csr& table, std::span<const int> row)
    {
        table.entries.insert(table.entries.end(), row.begin(), row.end());
        table.starts.push_back(static_cast<int>(table.entries.size()));
    }

    [[nodiscard]] Vec3f point(int v) const { return measure_.point(v); }
    [[nodiscard]] bool isProtected(int e) const { return measure_.isProtected(e); }
    [[nodiscard]] bool isCandidate(int e) const { return measure_.isCandidate(e); }
    [[nodiscard]] FreeAngles freeAngles(int v) const { return measure_.freeAngles(v); }

    // Nodes are the chain end points in order of first appearance; a closed
    // chain counts twice at its node.
    void buildNodes()
    {
        nodeOfPoint_.assign(static_cast<std::size_t>(adjacency_.pointCount()), -1);
        nodeChains_.clear();
        chainNodes_.clear();
        auto nodeOf = [&](int v) {
            int& k = nodeOfPoint_[static_cast<std::size_t>(v)];
            if (k < 0) {
                k = static_cast<int>(nodeChains_.size());
                nodeChains_.push_back(0);
            }
            ++nodeChains_[static_cast<std::size_t>(k)];
            return k;
        };
        for (int j = 0; j < chains_.rows(); ++j) {
            const auto chain = chains_.row(j);
            const int a = nodeOf(chain.front());
            const int b = nodeOf(chain.back());
            chainNodes_.emplace_back(a, b);
        }
    }

    [[nodiscard]] int chainsAt(int node) const
    {
        return node < 0 ? 0 : nodeChains_[static_cast<std::size_t>(node)];
    }

    [[nodiscard]] bool isClosed(int j) const
    {
        const auto chain = chains_.row(j);
        return chain.front() == chain.back();
    }

    void reject(int j, bool rescue)
    {
        const auto edges = chainEdges_.row(j);
        for (const int e : edges) {
            candidates_.set(static_cast<std::size_t>(e), false);
        }
        if (rescue && soft_.test(static_cast<std::size_t>(edges[0]))) {
            append(rescuedPoints_, chains_.row(j));
            append(rescuedEdges_, edges);
        }
        rejected_ = true;
    }

    // Preliminary class (or a first decision).
    void firstPass(int j)
    {
        Measures& m = measures_[static_cast<std::size_t>(j)];
        m.preClass = 0;
        m.decision = 0;
        m.classOne = false;
        m.rescue = false;
        const bool closed = isClosed(j);
        const float L = m.length;
        const double angleMax = wide(m.angleMax);
        const double angleMin = wide(m.angleMin);
        const double jumpMax = wide(m.jumpMax);
        const double bendMax = wide(m.bendMax);
        auto shortCheck = [&] {
            if (L < 2.0f * h_ && m.angleMax * 0.5f < m.rivalMax) {
                reject(j, false);
            }
        };
        if (!(L < h_)) {
            if (angleMin > degD(120) && L < 2.0f * h_) {
                m.decision = 1;
                return;
            }
            if (angleMax < degD(44.5) && L < 5.0f * h_ && jumpMax >= degD(2) && bendMax >= degD(2)) {
                m.decision = 1;
                return;
            }
            if ((!closed && jumpMax < degD(25) && bendMax < degD(5) && angleMin > degD(44)) ||
                (jumpMax < degD(4) && bendMax < degD(30) && angleMin > degD(60)) ||
                (jumpMax < degD(3) && bendMax < degD(2) && m.angleMin > creaseAngle_) ||
                (jumpMax < degD(25) && bendMax < degD(30) && angleMin > degD(60))) {
                m.preClass = 3;
            } else if (m.weightMin != 1.0f) {
                m.preClass = 1;
            } else if (angleMax <= degD(60) && !closed) {
                shortCheck();
            }
        } else if (closed) {
            m.decision = 1;
            if (angleMin > degD(50) && wide(m.angleMax - m.angleMin) < degD(10) && jumpMax < degD(5) &&
                bendMax < degD(30)) {
                m.decision = 0;
            }
        } else if (L < hQuad_ / 10.0f || angleMax > degD(100) || angleMin < degD(60) || jumpMax > degD(3)) {
            m.decision = 1;
        } else if (angleMax <= degD(60)) {
            shortCheck();
        }
    }

    // Counts at one chain end.
    struct End {
        int point = 0;
        int edge = 0;
        int chains = 0;     // chains at its node, the own one included
        int protectedEdges = 0;
        int softEdges = 0;
    };

    [[nodiscard]] End endAt(int point, int edge, int node) const
    {
        End out{point, edge, chainsAt(node), 0, 0};
        for (const int f : adjacency_.pointEdges(point)) {
            out.protectedEdges += isProtected(f) ? 1 : 0;
            out.softEdges += soft_.test(static_cast<std::size_t>(f)) ? 1 : 0;
        }
        return out;
    }

    // Angles from the chain end edge to the next candidate (or protected) edge
    // around the end point, one in each fan direction.
    [[nodiscard]] std::vector<float> continuationAngles(int edge, int v) const
    {
        std::vector<float> out;
        const std::vector<int> fan = vertexFanEdges(adjacency_, v, edge);
        const std::size_t m = fan.size();
        if (m < 2) {
            return out;
        }
        auto chosen = [&](int f) { return isCandidate(f) || isProtected(f); };
        const Vec3f p = point(v);
        const Vec3f along = point(adjacency_.otherEnd(edge, v)) - p;
        auto angleTo = [&](int f) { return vectorAngle(along, point(adjacency_.otherEnd(f, v)) - p); };
        int first = -1;
        for (std::size_t i = 1; i + 1 < m; ++i) {
            if (chosen(fan[i])) {
                first = fan[i];
                out.push_back(angleTo(first));
                break;
            }
        }
        if (first != -1 || fan.front() != fan.back()) {
            for (std::size_t i = m - 2; i >= 1; --i) {
                if (chosen(fan[i]) && fan[i] != first) {
                    out.push_back(angleTo(fan[i]));
                    break;
                }
            }
        }
        return out;
    }

    void secondPass(int j, int round)
    {
        Measures& m = measures_[static_cast<std::size_t>(j)];
        Next next = decide(j, m);
        if (next == Next::Late) {
            late(j, m, round);
            next = Next::Assign;
        }
        if (next == Next::Assign) {
            if (m.decision == 1) {
                next = Next::Reject;
            } else if (m.decision == 2) {
                next = Next::Soft;
            } else {
                classes_[static_cast<std::size_t>(j)] = m.classOne ? 1 : 2;
                return;
            }
        }
        if (next == Next::Reject) {
            classes_[static_cast<std::size_t>(j)] = -1;
            reject(j, m.rescue);
            return;
        }
        m.decision = 2;
        classes_[static_cast<std::size_t>(j)] = 3;
        for (const int e : chainEdges_.row(j)) {
            soft_.set(static_cast<std::size_t>(e));
        }
    }

    [[nodiscard]] static Next rejectNow(Measures& m)
    {
        m.decision = 1;
        m.rescue = false;
        return Next::Reject;
    }

    Next decide(int j, Measures& m)
    {
        if (m.decision == 1) {
            return Next::Reject;
        }
        if (isClosed(j)) {
            return Next::Late;
        }
        const auto chain = chains_.row(j);
        const auto edges = chainEdges_.row(j);
        const auto [nodeA, nodeB] = chainNodes_[static_cast<std::size_t>(j)];
        const End ends[2] = {endAt(chain.front(), edges.front(), nodeA), endAt(chain.back(), edges.back(), nodeB)};
        const int n0 = ends[0].protectedEdges + ends[0].chains;
        const int n1 = ends[1].protectedEdges + ends[1].chains;
        const int c = n0 + n1;
        const float L = m.length;
        const double angleMax = wide(m.angleMax);
        const double angleMin = wide(m.angleMin);
        const double jumpMax = wide(m.jumpMax);
        const double bendMax = wide(m.bendMax);

        const bool twoEdges = edges.size() >= 2;
        if (m.preClass == 3 && wide(m.angleMean) >= degD(44) && wide(m.smoothness) >= 0.9 && L >= h_ &&
            ((bendMax < degD(2) && jumpMax < degD(4) && twoEdges && L >= 2.0f * h_) ||
             (bendMax < degD(10) && jumpMax < degD(2.5) && twoEdges && L >= 2.0f * h_) ||
             (angleMax < degD(91) && angleMin > degD(89)))) {
            m.decision = 3;
            return Next::Assign;
        }

        float s = 0.0f;
        if (n0 < 2 && n1 < 2) {
            if (L < h_) {
                reject(j, false);
            } else if (!(L < 2.5f * h_)) {
                s = 0.0f;
            } else {
                s = 0.5f;
            }
        }
        int straight = 0;
        int open = 0;
        bool quietEnd = false;
        for (const End& e : ends) {
            const FreeAngles free = freeAngles(e.point);
            if (wide(free.max) < degD(5) && e.chains >= 2 && e.chains <= 4) {
                quietEnd = true;
            } else if (e.chains != 1) {
                quietEnd = false;
            } else {
                ++open;
                continue;
            }
            const std::vector<float> turns = continuationAngles(e.edge, e.point);
            if (e.chains == 2 && L < h_) {
                if (wide(free.max) <= 0.7 * angleMax) {
                    if (wide(free.sum) > degD(30) && angleMax < degD(60)) {
                        s += 0.5f;
                    }
                } else {
                    return rejectNow(m);
                }
            }
            for (const float w : turns) {
                const float dev = degD(135) > wide(w) ? static_cast<float>(std::fabs(wide(w) - degD(90))) : piF - w;
                const double d = wide(dev);
                if (L < h_) {
                    if (d >= degD(20)) {
                        if (d < degD(40)) {
                            if (!quietEnd) {
                                s += 0.3f;
                            }
                        } else {
                            s += quietEnd ? 0.5f : 1.0f;
                        }
                    }
                } else if (d < degD(3)) {
                    ++straight;
                    s -= 0.2f;
                } else if (d >= degD(20)) {
                    if (d < degD(40)) {
                        if (!quietEnd) {
                            s += 0.3f;
                        }
                    } else if (d < degD(60)) {
                        s += quietEnd ? 0.5f : 1.0f;
                    } else {
                        s += 1.5f;
                    }
                }
            }
        }
        if (L < h_ && open > 0) {
            return rejectNow(m);
        }

        bool border;
        if (!(L < 3.0f * h_)) {
            border = c > 3;
        } else if (c >= 4) {
            border = true;
        } else {
            // Short chains with few neighbours: the smallest edge angle in
            // degrees must reach a length dependent limit (80, from h on
            // falling towards 50 at 3h; a NaN passes), and the chain must run
            // almost straight without a strong rival.
            const float limit = L < h_ ? 80.0f : ((L - h_) * -30.0f) / (3.0f * h_ - h_) + 80.0f;
            const float smallest = degreesOf(m.angleMin);
            const bool rivalBelowBend = m.rivalMax < rivalFactor * m.bendMax;
            const bool straight = (jumpMax < degD(15) && bendMax < degD(4) && rivalBelowBend) ||
                                  (jumpMax < degD(4) && bendMax < degD(15) && angleMin > degD(60) && rivalBelowBend);
            if (limit > smallest || !straight) {
                return rejectNow(m);
            }
            border = false;
        }
        if (!border) {
            const float limit = ((m.angleMax - degF(30)) * (3.0f * h_ - 7.0f * h_)) / degF(30) + 7.0f * h_;
            if (angleMax < degD(60) && L < limit) {
                return rejectNow(m);
            }
            if (!(wide(m.angleMean) >= degD(90) && angleMin >= degD(60))) {
                if (jumpMax > degD(4)) {
                    s += ((m.jumpMax - degF(4)) * 0.5f) / jumpRampWidth + 0.0f;
                }
                int smooth = (m.angleMean / 15.0f <= m.jumpMax || angleMin <= 0.9 * angleMax) ? 0 : 1;
                if (angleMin < 0.8 * angleMax) {
                    s += (m.angleMax / m.angleMin - 1.0f) * 0.8f;
                }
                const float r = (m.ratioMax + m.ratioMean) * 0.5f;
                if (r > 0.15f) {
                    s += ((r - 0.15f) * 0.5f) / 0.18f + 0.0f;
                }
                if (bendMax < degD(2)) {
                    smooth += 1;
                    s -= 0.5f;
                }
                if (m.preClass == 3) {
                    s -= 0.5f;
                }
                if (smooth < 1) {
                    if (wide(s) > 0.4) {
                        return Next::Soft;
                    }
                    if (c != 3) {
                        m.classOne = true;
                        return Next::Assign;
                    }
                } else if (c != 3 && m.preClass != 3 && wide(m.angleMean) < degD(60)) {
                    return Next::Soft;
                }
            }
        }

        // Border-like chains.
        if (n0 - ends[0].softEdges <= 1 && n1 - ends[1].softEdges <= 1 && adjacency_.firstBoundaryEdgeAt(chain.front()) < 0 &&
            adjacency_.firstBoundaryEdgeAt(chain.back()) < 0 && angleMin < degD(50) &&
            wide(m.angleMax - m.angleMin) > degD(15) && m.ratioMax > 0.25f) {
            m.classOne = true;
        }
        if (m.preClass == 3) {
            if (s >= 1.5f) {
                m.decision = 1;
                m.rescue = true;
            } else if (s >= 1.0f) {
                m.decision = 2;
            } else {
                m.decision = 3;
            }
            return Next::Late;
        }
        if (s >= 1.0f) {
            m.decision = 1;
            m.rescue = false;
            return Next::Late;
        }
        if (m.preClass == 0 && (s <= 0.0f || straight > 2) && m.weightMin == 1.0f) {
            m.decision = 3;
            return Next::Assign;
        }
        if (!(s <= 0.4f) || m.weightMin < 0.7f || m.smoothness < 0.9f) {
            if (!(wide(s) <= 0.2) || wide(m.jumpMin) >= degD(5) || std::fabs(wide(m.angleMean) - degD(90)) >= degD(5)) {
                m.decision = 2;
                return Next::Late;
            }
            m.decision = 3;
            return Next::Assign;
        }
        return Next::Soft;
    }

    // Late rejection of uneven chains in the first rounds (the class is assigned afterwards).
    void late(int j, const Measures& m, int round)
    {
        if (round < lateRounds && m.smoothness < 0.8f) {
            const bool uneven = m.jumpMax < degF(1) ? m.bendMax >= degF(30) : true;
            if (uneven && (m.weightMax - m.weightMin >= 0.15f || m.weightMin <= 0.4f || m.jumpMax >= degF(4) ||
                           m.bendMax >= degF(30))) {
                reject(j, true);
            }
        }
    }

    BitArray candidates_;
    BitArray soft_; // edges of chains that got class 3, kept over all rounds
    ChainMeasure measure_;
    const Adjacency& adjacency_;
    float creaseAngle_;
    float h_;
    float hQuad_;
    Csr chains_;
    Csr chainEdges_;
    std::vector<int> classes_;
    std::vector<Measures> measures_;
    std::vector<int> nodeOfPoint_;
    std::vector<int> nodeChains_;
    std::vector<std::pair<int, int>> chainNodes_;
    Csr rescuedPoints_;
    Csr rescuedEdges_;
    bool rejected_ = false;
};

} // namespace

FeatureChains::Lengths FeatureChains::referenceLengths(float area, int quadCount)
{
    int divisor = quadCount < 1000 ? 1000 : quadCount;
    divisor = divisor < 10001 ? divisor : 10000;
    Lengths out;
    out.reference = std::sqrt(area / static_cast<float>(divisor));
    out.quad = std::sqrt(area / std::max(static_cast<float>(quadCount), 1.0f));
    return out;
}

FeatureChains::Result FeatureChains::classify(const Mesh& mesh, const Options& options)
{
    return Classifier(mesh, options).run();
}

// The chain measures of the check: no candidates, no protected edges.
class FeatureChains::CreaseLineCheck::Measure {
public:
    explicit Measure(const Mesh& mesh)
        : none_(static_cast<std::size_t>(mesh.adjacency().edgeCount())), chains_(mesh, none_, nullptr)
    {
    }

    [[nodiscard]] const ChainMeasure& chains() const noexcept { return chains_; }

private:
    BitArray none_; // declared before chains_, which keeps a reference to it
    ChainMeasure chains_;
};

FeatureChains::CreaseLineCheck::CreaseLineCheck(const Mesh& mesh) : measure_(std::make_unique<const Measure>(mesh)) {}

FeatureChains::CreaseLineCheck::~CreaseLineCheck() = default;

bool FeatureChains::CreaseLineCheck::operator()(std::span<const int> points) const
{
    if (points.size() < 2) {
        return false;
    }
    const ChainMeasure& chains = measure_->chains();
    const Adjacency& adjacency = chains.adjacency();
    std::vector<int> edges;
    edges.reserve(points.size() - 1);
    for (std::size_t i = 1; i < points.size(); ++i) {
        const int e = adjacency.edgeBetween(points[i - 1], points[i]);
        if (e < 0) {
            return false;
        }
        edges.push_back(e);
    }
    const Measures m = chains.measure(edges, points);
    return wide(m.angleMin) > degD(45) ||
           (!(degD(30) > wide(m.angleMax)) && !(degD(30) > wide(m.angleWeighted)) && wide(m.bendMax) < degD(8) &&
            wide(m.jumpMax) < degD(10));
}

} // namespace areko::remesher::prep
