#include "elimination_columns.h"

#include <utility>

namespace areko::remesher::param {

EliminationColumns::EliminationColumns(std::vector<TermList>& rows, int variableCount, const EliminationSetup& setup)
{
    // Column of every variable, -1 for variables without one.
    std::vector<int> columnOf(static_cast<std::size_t>(variableCount), -1);
    int count = 0;
    const auto number = [&](int variable) {
        int& column = columnOf[static_cast<std::size_t>(variable)];
        if (column == -1) {
            column = count++;
        }
        return column;
    };
    const auto flag = [&](int variable) { return setup.flags[static_cast<std::size_t>(variable)]; };
    // Every term of every row, rows and terms in order.
    const auto forEachTerm = [&](auto visit) {
        for (TermList& row : rows) {
            for (Term& term : row) {
                visit(term);
            }
        }
    };

    if (setup.mode == EliminationMode::IsoLines) {
        forEachTerm([&](const Term& term) {
            if (flag(term.index) == continuousFlag) {
                number(term.index);
            }
        });
    }
    for (const int variable : setup.preferred) {
        number(variable);
    }
    if (setup.mode == EliminationMode::SoftSystem) {
        // Each pass rewrites in place and tests the entry as it stands. The
        // new numbers stay below the constant column as long as fewer
        // variables are numbered than it is large.
        const int constant = setup.constantColumn;
        const auto rewrite = [&](auto wanted) {
            forEachTerm([&](Term& term) {
                if (wanted(term.index)) {
                    term.index = number(term.index);
                }
            });
        };
        rewrite([&](int index) { return index < constant; });
        rewrite([&](int index) { return constant < index; });
        split_ = count;
        rewrite([&](int index) { return index == constant; });
    } else {
        forEachTerm([&](const Term& term) {
            if (flag(term.index) != keptFreeFlag) {
                number(term.index);
            }
        });
        split_ = count;
        forEachTerm([&](Term& term) { term.index = number(term.index); });
    }

    variableOf_.assign(static_cast<std::size_t>(count), 0);
    for (int variable = 0; variable < variableCount; ++variable) {
        const int column = columnOf[static_cast<std::size_t>(variable)];
        if (column != -1) {
            variableOf_[static_cast<std::size_t>(column)] = variable;
        }
    }
}

void EliminationColumns::swap(std::vector<TermList>& rows, int a, int b)
{
    for (TermList& row : rows) {
        for (Term& term : row) {
            if (term.index == a) {
                term.index = b;
            } else if (term.index == b) {
                term.index = a;
            }
        }
        solver::ListMatrix::order(row, solver::ListMatrix::SameIndex::Separate);
    }
    std::swap(variableOf_[static_cast<std::size_t>(a)], variableOf_[static_cast<std::size_t>(b)]);
}

} // namespace areko::remesher::param
