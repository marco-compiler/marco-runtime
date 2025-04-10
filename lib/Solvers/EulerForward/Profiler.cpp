#include "marco/Runtime/Solvers/EulerForward/Profiler.h"
#include <iostream>

#ifdef MARCO_PROFILING

namespace marco::runtime::profiling {
EulerForwardProfiler::EulerForwardProfiler() : Profiler("Euler forward") {
  registerProfiler(*this);
}

void EulerForwardProfiler::reset() {
  std::lock_guard<std::mutex> lockGuard(mutex);

  nonStateVariables.reset();
  stateVariables.reset();
}

void EulerForwardProfiler::print() const {
  std::lock_guard<std::mutex> lockGuard(mutex);

  std::cerr << "Time spent on state variables computation: "
            << stateVariables.totalElapsedTime<std::milli>() << " ms\n";

  std::cerr << "Time spent on non-state variables computation: "
            << nonStateVariables.totalElapsedTime<std::milli>() << " ms\n";
}

EulerForwardProfiler &eulerForwardProfiler() {
  static EulerForwardProfiler obj;
  return obj;
}
} // namespace marco::runtime::profiling

#endif // MARCO_PROFILING
