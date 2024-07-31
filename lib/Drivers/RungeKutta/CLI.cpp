#ifdef CLI_ENABLE

#include "marco/Runtime/Drivers/RungeKutta/CLI.h"
#include "marco/Runtime/Solvers/RungeKutta/Options.h"

namespace marco::runtime::rungekutta {
std::string CommandLineOptions::getTitle() const { return "Runge-Kutta"; }

void CommandLineOptions::printCommandLineOptions(std::ostream &os) const {
  // clang-format off
  os << "  --time-step=<value>        Set the time step (in seconds). Defaults to " << getOptions().timeStep << "." << std::endl;
  os << "  --tolerance=<value>        Set the tolerance for adaptive Runge-Kutta methods. If set, the provided time step is used only for the first step. Defaults to " << getOptions().tolerance << "." << std::endl;
  os << "  --max-iterations=<value>   Set the maximum number of iterations for adaptive Runge-Kutta methods. Defaults to " << getOptions().maxIterations << "." << std::endl;
  // clang-format on
}

void CommandLineOptions::parseCommandLineOptions(
    const argh::parser &options) const {
  // clang-format off
  options("time-step") >> getOptions().timeStep;
  options("tolerance") >> getOptions().tolerance;
  options("max-iterations") >> getOptions().maxIterations;
  // clang-format on
}
} // namespace marco::runtime::rungekutta

#endif // CLI_ENABLE
