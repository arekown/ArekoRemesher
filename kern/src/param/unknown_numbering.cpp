#include "unknown_numbering.h"

namespace areko::remesher::param {

void UnknownNumbering::assign(const SystemContext& context, System& system)
{
    const int coordinateEnd = 2 * system.copyPointCount;
    system.variableToUnknown.assign(static_cast<std::size_t>(system.variableCount), -1);
    system.unknownToVariable.clear();
    for (int variable = 0; variable < system.variableCount; ++variable) {
        const bool mirrored = variable < coordinateEnd && context.mirrored(variable / 2);
        if (mirrored || system.constraints.ofVariable(variable) != -1) {
            continue;
        }
        system.variableToUnknown[static_cast<std::size_t>(variable)] = system.unknownCount();
        system.unknownToVariable.push_back(variable);
    }
}

int UnknownNumbering::renumber(const System& system, std::span<TermList> rows)
{
    int withoutUnknown = 0;
    for (TermList& row : rows) {
        for (Term& term : row) {
            term.index = system.variableToUnknown[static_cast<std::size_t>(term.index)];
            if (term.index == -1) {
                ++withoutUnknown;
            }
        }
    }
    return withoutUnknown;
}

void UnknownNumbering::listIntegers(System& system)
{
    system.integerUnknowns.clear();
    system.integerPeriods.clear();
    for (int unknown = 0; unknown < system.unknownCount(); ++unknown) {
        const int variable = system.unknownToVariable[static_cast<std::size_t>(unknown)];
        if (system.flags[static_cast<std::size_t>(variable)] != variable_flag::continuous) {
            system.integerUnknowns.push_back(unknown);
        }
    }
    if (system.periodChanges > 0) {
        for (const int unknown : system.integerUnknowns) {
            const int variable = system.unknownToVariable[static_cast<std::size_t>(unknown)];
            system.integerPeriods.push_back(static_cast<double>(system.periods[static_cast<std::size_t>(variable)]));
        }
    }
}

} // namespace areko::remesher::param
