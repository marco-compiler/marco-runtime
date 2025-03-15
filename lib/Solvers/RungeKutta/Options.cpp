#include "marco/Runtime/Solvers/RungeKutta/Options.h"

using namespace ::marco::runtime::rungekutta;

namespace marco::runtime::rungekutta {
Options &getOptions() {
  static Options options;
  return options;
}
} // namespace marco::runtime::rungekutta
