#ifndef MARCO_RUNTIME_DRIVERS_RUNGEKUTTA_DRIVER_H
#define MARCO_RUNTIME_DRIVERS_RUNGEKUTTA_DRIVER_H

#include "marco/Runtime/Drivers/Driver.h"

namespace marco::runtime
{
class RungeKutta : public Driver
{
public:
  RungeKutta(Simulation* simulation);

#ifdef CLI_ENABLE
  std::unique_ptr<cli::Category> getCLIOptions() override;
#endif // CLI_ENABLE

  int run() override;
};
} // namespace marco::runtime

//===---------------------------------------------------------------------===//
// Functions defined inside the module of the compiled model
//===---------------------------------------------------------------------===//

extern "C"
{
void tryStep(double timeStep);
double estimateError();
void acceptStep();
void updateNonStateVariables();
}

#endif // MARCO_RUNTIME_DRIVERS_RUNGEKUTTA_DRIVER_H
