#pragma once

#include "ngon_split_round.h"

#include "mesh/bit_array.h"

#include <vector>

namespace areko::remesher::extract {

// Makes faces with an odd corner count even in pairs: a breadth-first search
// from every odd face of the set grows regions over inner edges; for every
// pair of touching regions the cheapest edge between them closes a path of
// faces from one odd face to the other. Along the path, each crossed edge gets
// its middle point: the two end faces gain one corner (and are split), the
// faces between are cut along the path. Paths through a face already handled
// and paths with an end face of eight or more corners are skipped.
class ParityBalance {
public:
    static void run(SplitRound& round, const mesh::BitArray& faces);

private:
    struct Search {
        std::vector<float> distance;
        std::vector<int> source;
        std::vector<int> entry; // edge the face was reached over, -1 for sources
    };

    static Search search(const mesh::Mesh& mesh, const mesh::BitArray& faces);
    // Faces from `start` back to its source and the edges between them.
    static void path(const mesh::Mesh& mesh, const Search& found, int start, std::vector<int>& faces,
                     std::vector<int>& edges);
    static void splitPath(SplitRound& round, const std::vector<int>& faces, const std::vector<int>& edges);
    static void endFace(SplitRound& round, int edge, int face, int& middle);
    static void innerFace(SplitRound& round, int entry, int exit, int face, int entryMiddle, int& exitMiddle);
};

} // namespace areko::remesher::extract
