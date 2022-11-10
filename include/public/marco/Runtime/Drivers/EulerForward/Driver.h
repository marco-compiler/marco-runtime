#ifndef MARCO_RUNTIME_DRIVERS_EULERFORWARD_DRIVER_H
#define MARCO_RUNTIME_DRIVERS_EULERFORWARD_DRIVER_H

#include "marco/Runtime/Drivers/Driver.h"

namespace marco::runtime
{
  class EulerForward : public Driver
  {
    public:
      EulerForward(Simulation* simulation);

      std::unique_ptr<cli::Category> getCLIOptions() override;

      int run() override;
  };
}

//===---------------------------------------------------------------------===//
// Functions defined inside the module of the compiled model
//===---------------------------------------------------------------------===//

extern "C"
{
  void calcIC(void* data);
  void updateNonStateVariables(void* data);
  void updateStateVariables(void* data, double timeStep);
}

#endif // MARCO_RUNTIME_DRIVERS_EULERFORWARD_DRIVER_H