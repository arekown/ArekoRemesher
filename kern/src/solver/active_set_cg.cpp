#include "active_set_cg.h"

namespace areko::remesher::solver {

namespace {

constexpr double nothingToDo = 1e-5;
constexpr int defaultMaxIterations = 1000;

// q = A p over the columns of the active set; returns the grown set: every row
// reached, in the order found, then the members of the set not reached
// (their q is 0).
class FrontProduct {
public:
    explicit FrontProduct(std::size_t n) : marks_(n, 0u) {}

    void apply(const ListMatrix& matrix, const std::vector<double>& p, const std::vector<int>& set,
               std::vector<double>& q, std::vector<int>& grown)
    {
        grown.clear();
        if (++stamp_ == 0u) {
            marks_.assign(marks_.size(), 0u);
            stamp_ = 1u;
        }
        for (const int j : set) {
            const double pj = p[static_cast<std::size_t>(j)];
            for (const SparseEntry& e : matrix.list(j)) {
                const auto i = static_cast<std::size_t>(e.index);
                if (marks_[i] == stamp_) {
                    q[i] = e.value * pj + q[i];
                } else {
                    q[i] = e.value * pj;
                    marks_[i] = stamp_;
                    grown.push_back(e.index);
                }
            }
        }
        for (const int j : set) {
            const auto u = static_cast<std::size_t>(j);
            if (marks_[u] != stamp_) {
                q[u] = 0.0;
                marks_[u] = stamp_;
                grown.push_back(j);
            }
        }
    }

private:
    std::vector<unsigned> marks_;
    unsigned stamp_ = 0u;
};

} // namespace

bool ActiveSetCg::solve(const Options& options, std::vector<double>& x, std::vector<int> set,
                        const ListMatrix& matrix, std::span<const double> rhs)
{
    const auto n = static_cast<std::size_t>(matrix.listCount());
    if (!matrix.square() || rhs.size() != n || x.size() != n) {
        return false;
    }
    std::vector<double> r(n, 0.0);
    std::vector<double> q(n, 0.0);
    for (const int i : set) {
        double s = 0.0;
        for (const SparseEntry& e : matrix.list(i)) {
            s = s + e.value * x[static_cast<std::size_t>(e.index)];
        }
        r[static_cast<std::size_t>(i)] = rhs[static_cast<std::size_t>(i)] - s;
    }
    const int maxIterations = options.maxIterations == -1 ? defaultMaxIterations : options.maxIterations;
    std::vector<double> p = r;
    double rho = 0.0;
    for (const int i : set) {
        rho = r[static_cast<std::size_t>(i)] * r[static_cast<std::size_t>(i)] + rho;
    }
    if (set.empty() || nothingToDo > rho) {
        return true;
    }
    const double goal = (options.tolerance * options.tolerance) * rho;
    int k = 0;
    if (maxIterations >= 1 && !(goal > rho)) {
        FrontProduct product(n);
        std::vector<int> grown;
        do {
            product.apply(matrix, p, set, q, grown);
            double sigma = 0.0;
            for (const int i : set) {
                sigma = sigma + p[static_cast<std::size_t>(i)] * q[static_cast<std::size_t>(i)];
            }
            if (sigma == 0.0) {
                return false;
            }
            const double alpha = rho / sigma;
            for (const int i : set) {
                const auto u = static_cast<std::size_t>(i);
                x[u] = p[u] * alpha + x[u];
            }
            set.swap(grown);
            double next = 0.0;
            for (const int i : set) {
                const auto u = static_cast<std::size_t>(i);
                r[u] = r[u] - q[u] * alpha;
                next = next + r[u] * r[u];
            }
            const double beta = next / rho;
            for (const int i : set) {
                const auto u = static_cast<std::size_t>(i);
                p[u] = p[u] * beta + r[u];
            }
            ++k;
            rho = next;
        } while (k != maxIterations && !(goal > rho));
    }
    // Convergence in the last allowed iteration counts as success.
    return k < maxIterations || goal > rho;
}

} // namespace areko::remesher::solver
