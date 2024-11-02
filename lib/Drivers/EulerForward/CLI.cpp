#ifdef CLI_ENABLE

#include "marco/Runtime/Drivers/EulerForward/CLI.h"
#include "marco/Runtime/Solvers/EulerForward/Options.h"

namespace marco::runtime::eulerforward {
std::string CommandLineOptions::getTitle() const { return "Euler forward"; }

void CommandLineOptions::printCommandLineOptions(std::ostream &os) const {
  // clang-format off
  os << "  --time-step=<value>         Set the time step (in seconds). Defaults to " << getOptions().timeStep << "." << std::endl;
  os << "  --snapshot-steps=<value>    Print simulation state every <value> time steps. Defaults to " << getOptions().snapshotSteps << "." << std::endl;
  // clang-format on
}

void CommandLineOptions::parseCommandLineOptions(
    const argh::parser &options) const {
  // clang-format off
  options("time-step") >> getOptions().timeStep;
  options("snapshot-steps") >> getOptions().snapshotSteps;
  // clang-format on
}
} // namespace marco::runtime::eulerforward

#endif // CLI_ENABLE
