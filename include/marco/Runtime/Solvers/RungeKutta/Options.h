#ifndef MARCO_RUNTIME_SOLVERS_RUNGEKUTTA_OPTIONS_H
#define MARCO_RUNTIME_SOLVERS_RUNGEKUTTA_OPTIONS_H

#include "marco/Runtime/Support/Mangling.h"

namespace marco::runtime::rungekutta
{
struct Options
{
  double timeStep = 0.1;
  double tolerance = 1e-6;
  uint64_t maxIterations = 10;
};

Options& getOptions();
}

#endif // MARCO_RUNTIME_SOLVERS_RUNGEKUTTA_OPTIONS_H
