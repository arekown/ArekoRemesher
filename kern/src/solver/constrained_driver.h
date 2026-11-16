#pragma once

#include "mesh/bit_array.h"
#include "solver/active_set_cg.h"
#include "solver/direct_solver.h"
#include "solver/sparse_matrix.h"
#include "util/index_queue.h"

#include <variant>
#include <vector>

namespace areko::remesher::run {
class Progress;
}

namespace areko::remesher::solver {

// One round of the driver as seen by its callback. The callback fills the
// answer part; `finished` is never reset by the driver.
struct DriverRequest {
    enum class Answer { None, Delivered, Stop };
    // One selected position of the integer list and its value.
    struct Selection {
        int position = 0;
        double value = 0.0;
    };

    int round = 0;
    const std::vector<double>* x = nullptr;  // latest solution
    int fixedCount = 0;                      // selections fixed so far
    int changedRounds = 0;                   // rounds in which the callback changed the system

    // Answer.
    std::vector<Selection> selection;
    const std::vector<int>* extra = nullptr; // unknowns touched by the callback (start of the local update)
    Answer answer = Answer::None;
    bool systemChanged = false;
    bool allowLocal = true;
    bool noMark = false;                     // do not mark the delivered positions as fixed
    bool finished = false;
};

class DriverCallback {
public:
    virtual ~DriverCallback() = default;
    // Asked once per round before the driver picks positions itself.
    virtual void onRequest(DriverRequest& request) = 0;
    // A selection is about to be fixed.
    virtual void onSelection(const DriverRequest& /*request*/) {}
    // Every position of the integer list has been fixed; the system is solved once more.
    virtual void onAllFixed() {}
};

struct DriverOptions {
    bool keepGoing = false;        // go on after the last fixing until the callback reports finished
    bool localUpdate = false;      // CG on a growing active set instead of a full solve (list systems)
    ActiveSetCg::Options local;    // CG of the local update
    int fullSolveEvery = 80;       // with local updates: every n-th round solves fully (<= 0: never)
    bool batchMode = false;        // batch rounding instead of greedy rounding
    double batchBudget = 0.0;      // error budget of one batch
};

// Solves a symmetric system and rounds the unknowns of an integer list one
// selection at a time: each round the callback may deliver a selection or
// change the system, otherwise the driver picks by greedy or batch rounding;
// selected unknowns are fixed in matrix and right-hand side, then the system
// is solved again (fully or by a local update).
class ConstrainedDriver {
public:
    enum class Status { Ok, SolveFailed };

    ConstrainedDriver(CompactMatrix& matrix, std::vector<double>& rhs, DriverCallback& callback,
                      const run::Progress& progress, const DriverOptions& options = {});
    ConstrainedDriver(ListMatrix& matrix, std::vector<double>& rhs, DriverCallback& callback,
                      const run::Progress& progress, const DriverOptions& options = {});

    // Unknowns to round (their index is the position) and optional periods per position.
    void setIntegerUnknowns(std::vector<int> unknowns, std::vector<double> periods = {});

    // Runs until every position is fixed (true) or the run ends early (false).
    // Cancelling throws from the progress check.
    bool run(std::vector<double>& x);

    // Turns unknown `var` into x[var] = value: its row and column become the
    // unit row, the old couplings move to the right-hand side.
    void fixVariable(int var, double value);

    [[nodiscard]] Status status() const noexcept { return status_; }
    [[nodiscard]] const std::vector<int>& integerUnknowns() const noexcept { return unknowns_; }
    [[nodiscard]] const std::vector<double>& periods() const noexcept { return periods_; }
    [[nodiscard]] const mesh::BitArray& fixedPositions() const noexcept { return fixed_; }
    // Whether a position of the integer list is still to be fixed.
    [[nodiscard]] bool isFree(int position) const { return !fixed_.testBounded(static_cast<std::size_t>(position)); }
    // The whole number nearest to the value of a position in `x`, half away
    // from zero; periods do not count. The whole number goes straight to double.
    [[nodiscard]] double nearestWhole(const std::vector<double>& x, int position) const;
    // The whole multiple of a period nearest to a value (the number of
    // periods rounded half away from zero) and its distance from the value,
    // counted in periods.
    struct PeriodRounding {
        double value = 0.0;
        double error = 0.0;
    };
    [[nodiscard]] static PeriodRounding toPeriod(double value, double period) noexcept;
    // Matrix of the path chosen at construction (CompactMatrix or ListMatrix);
    // asking for the other type throws std::bad_variant_access.
    template <class Matrix>
    [[nodiscard]] const Matrix& matrix() const
    {
        return *std::get<Matrix*>(matrix_);
    }

private:
    bool solveSystem(std::vector<double>& x);
    bool solveIn(const CompactMatrix& matrix, std::vector<double>& x);
    bool solveIn(const ListMatrix& matrix, std::vector<double>& x);
    void fixIn(CompactMatrix& matrix, int var, double value);
    void fixIn(ListMatrix& matrix, int var, double value);
    bool pickGreedy(const std::vector<double>& x, DriverRequest& request);
    void pickBatch(const std::vector<double>& x, DriverRequest& request);
    void batchStage(const std::vector<double>& x, double tolerance, double& batchError, DriverRequest& request);
    // Rounding error of a position (in units of its period when `withPeriods`
    // and periods are set) and its rounded value.
    void roundingOf(const std::vector<double>& x, int position, bool withPeriods, double& error,
                    double& value) const;
    // Free position with the strictly smallest rounding error (the first one
    // on ties) with its error and rounded value, or -1.
    int bestFree(const std::vector<double>& x, bool withPeriods, double& error, double& value) const;
    // Adds the position with its value to the selection and marks it fixed.
    void take(DriverRequest& request, int position, double value);
    bool localUpdate(const DriverRequest& request, std::vector<double>& x);

    std::variant<CompactMatrix*, ListMatrix*> matrix_;
    std::vector<double>& rhs_;
    DriverCallback& callback_;
    const run::Progress& progress_;
    DriverOptions options_;

    std::vector<int> unknowns_;
    std::vector<double> periods_;
    mesh::BitArray fixed_;
    Status status_ = Status::Ok;
    int zeroCounter_ = 0;          // couplings set to zero by fixing; >= 1 drops zeros before solving
    double batchMemo_ = 0.0;       // tolerance memo of the batch rounding (0 = run the stages)
    util::IndexQueue queue_;
    DirectSolver direct_;
};

} // namespace areko::remesher::solver
