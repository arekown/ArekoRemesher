#include "edge_chains.h"

#include <utility>
#include <vector>

namespace areko::remesher::mesh {

Csr EdgeChains::build(const Adjacency& adjacency, const BitArray& chosen, Csr* chainEdgeRows,
                      std::span<const int> separators)
{
    Options options;
    options.separators = separators;
    Result result = run(adjacency, chosen, {}, options);
    if (chainEdgeRows != nullptr) {
        *chainEdgeRows = std::move(result.edges);
    }
    return std::move(result.points);
}

EdgeChains::Result EdgeChains::build(const Adjacency& adjacency, const BitArray& chosen, const Options& options)
{
    return run(adjacency, chosen, {}, options);
}

EdgeChains::Result EdgeChains::build(const Adjacency& adjacency, std::span<const int> edgeList,
                                     const Options& options)
{
    BitArray chosen(static_cast<std::size_t>(adjacency.edgeCount()));
    std::vector<int> order;
    order.reserve(edgeList.size());
    for (const int e : edgeList) {
        if (!chosen.test(static_cast<std::size_t>(e))) {
            chosen.set(static_cast<std::size_t>(e));
            order.push_back(e);
        }
    }
    return run(adjacency, chosen, order, options);
}

EdgeChains::Result EdgeChains::run(const Adjacency& adjacency, const BitArray& chosen, std::span<const int> order,
                                   const Options& options)
{
    Result result;
    const int edges = adjacency.edgeCount();
    const int points = adjacency.pointCount();
    auto isChosen = [&](int e) { return static_cast<std::size_t>(e) < chosen.size() && chosen.test(static_cast<std::size_t>(e)); };

    std::vector<int> degree(static_cast<std::size_t>(points), 0);
    std::vector<int> touched;
    BitArray seen(static_cast<std::size_t>(points));
    auto take = [&](int e) {
        for (const int v : adjacency.edgeVertices(e)) {
            ++degree[static_cast<std::size_t>(v)];
            if (!seen.test(static_cast<std::size_t>(v))) {
                seen.set(static_cast<std::size_t>(v));
                touched.push_back(v);
            }
        }
    };
    if (order.empty()) {
        for (int e = 0; e < edges; ++e) {
            if (isChosen(e)) {
                take(e);
            }
        }
    } else {
        for (const int e : order) {
            take(e);
        }
    }
    if (touched.empty()) {
        return result;
    }
    std::vector<int> ends;
    std::vector<int> branches;
    for (const int v : touched) {
        const int d = degree[static_cast<std::size_t>(v)];
        if (d == 1) {
            ends.push_back(v);
        } else if (d >= 3) {
            branches.push_back(v);
        }
    }
    const bool pureLoops = ends.empty() && branches.empty();
    if (pureLoops) {
        branches.push_back(touched.front());
    }
    auto separate = [&](int t) {
        int& d = degree[static_cast<std::size_t>(t)];
        if (d == 2) {
            branches.push_back(t);
            d = 3;
        }
    };
    for (const int t : options.separators) {
        separate(t);
    }
    if (options.borderPointsSeparate) {
        for (const int v : adjacency.boundaryPoints()) {
            separate(v);
        }
    }
    result.nodes = ends;
    for (std::size_t i = 0; i < branches.size(); ++i) {
        // The start of pure loops stands once, also when it separates now.
        if (!(pureLoops && i > 0 && branches[i] == branches.front())) {
            result.nodes.push_back(branches[i]);
        }
    }
    auto copiesOf = [&](int v) {
        const int d = degree[static_cast<std::size_t>(v)];
        return d == 1 ? 1 : d - 1;
    };
    std::vector<int> stack;
    auto addStarts = [&](int v, int count) { stack.insert(stack.end(), static_cast<std::size_t>(count), v); };
    if (!options.startsFollowChains) {
        for (const int v : branches) {
            addStarts(v, copiesOf(v));
        }
    }

    BitArray used(static_cast<std::size_t>(edges));
    std::size_t freeEdges = 0;
    for (int e = 0; e < edges; ++e) {
        freeEdges += isChosen(e) ? 1 : 0;
    }
    std::size_t endsLeft = ends.size();
    std::size_t kept = 0; // entries at the front of the stack that are taken and not taken again
    bool refilled = false;

    auto firstFreeEdgeAt = [&](int v) {
        for (const int e : adjacency.pointEdges(v)) {
            if (isChosen(e) && !used.test(static_cast<std::size_t>(e))) {
                return e;
            }
        }
        return -1;
    };

    std::vector<int> chain;
    std::vector<int> chainEdges;
    while (freeEdges > 0) {
        int start;
        if (endsLeft > 0) {
            start = ends[--endsLeft];
        } else if (stack.size() > kept) {
            if (options.startsFollowChains) {
                start = stack[kept++];
            } else {
                start = stack.back();
                stack.pop_back();
            }
        } else if (!refilled && !branches.empty()) {
            refilled = true;
            for (const int v : branches) {
                addStarts(v, copiesOf(v));
            }
            start = stack[kept++];
        } else {
            int smallest = -1;
            for (int e = 0; e < edges && smallest < 0; ++e) {
                if (isChosen(e) && !used.test(static_cast<std::size_t>(e))) {
                    smallest = e;
                }
            }
            start = adjacency.edgeVertices(smallest)[0];
        }

        chain.clear();
        chainEdges.clear();
        int v = start;
        bool finished = false;
        while (!finished) {
            const int e = firstFreeEdgeAt(v);
            if (e < 0) {
                // Closed walk: the start lies on the last edge.
                if (chain.size() >= 3 && v == start) {
                    chain.push_back(start);
                }
                break;
            }
            used.set(static_cast<std::size_t>(e));
            --freeEdges;
            chain.push_back(v);
            chainEdges.push_back(e);
            v = adjacency.otherEnd(e, v);
            if (const int d = degree[static_cast<std::size_t>(v)]; d != 2) {
                chain.push_back(v);
                finished = true;
                if (options.startsFollowChains && d >= 3) {
                    addStarts(v, d - 1);
                }
            }
        }
        if (!chain.empty()) {
            result.points.entries.insert(result.points.entries.end(), chain.begin(), chain.end());
            result.points.starts.push_back(static_cast<int>(result.points.entries.size()));
            result.edges.entries.insert(result.edges.entries.end(), chainEdges.begin(), chainEdges.end());
            result.edges.starts.push_back(static_cast<int>(result.edges.entries.size()));
        }
    }
    return result;
}

} // namespace areko::remesher::mesh
