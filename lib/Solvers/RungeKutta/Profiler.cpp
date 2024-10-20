#include "marco/Runtime/Solvers/RungeKutta/Profiler.h"
#include <iostream>

#ifdef MARCO_PROFILING

namespace marco::runtime::profiling
{
RungeKuttaProfiler::RungeKuttaProfiler()
    : Profiler("Runge-Kutta")
{
  registerProfiler(*this);
}

void RungeKuttaProfiler::reset()
{
  std::lock_guard<std::mutex> lockGuard(mutex);

  variablesCopy.reset();
  stateVariables.reset();
  nonStateVariables.reset();
}

void RungeKuttaProfiler::print() const
{
  std::lock_guard<std::mutex> lockGuard(mutex);

  std::cerr << "Time spent on copying variables: "
            << variablesCopy.totalElapsedTime<std::milli>() << " ms\n";

  std::cerr << "Time spent on state variables computation: "
            << stateVariables.totalElapsedTime<std::milli>() << " ms\n";

  std::cerr << "Time spent on non-state variables computation: "
            << nonStateVariables.totalElapsedTime<std::milli>() << " ms\n";
}

RungeKuttaProfiler& rungeKuttaProfiler()
{
  static RungeKuttaProfiler obj;
  return obj;
}
} // namespace marco::runtime::profiling

#endif // MARCO_PROFILING