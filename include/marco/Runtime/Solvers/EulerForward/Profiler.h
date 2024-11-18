#ifndef MARCO_RUNTIME_SOLVERS_EULERFORWARD_PROFILER_H
#define MARCO_RUNTIME_SOLVERS_EULERFORWARD_PROFILER_H

#ifdef MARCO_PROFILING

#include "marco/Runtime/Profiling/Profiling.h"
#include "marco/Runtime/Profiling/Timer.h"
#include "marco/Runtime/Simulation/Options.h"
#include <mutex>

namespace marco::runtime::profiling {
class EulerForwardProfiler : public Profiler {
public:
  EulerForwardProfiler();

  void reset() override;

  void print() const override;

public:
  Timer stateVariables;
  Timer nonStateVariables;

  mutable std::mutex mutex;
};

EulerForwardProfiler &eulerForwardProfiler();
} // namespace marco::runtime::profiling

// clang-format off
#define EULER_FORWARD_PROFILER_STATEVAR_START                                   \
  if (::marco::runtime::simulation::getOptions().profiling) {                   \
    ::marco::runtime::profiling::eulerForwardProfiler().stateVariables.start(); \
  }

#define EULER_FORWARD_PROFILER_STATEVAR_STOP                                   \
  if (::marco::runtime::simulation::getOptions().profiling) {                  \
    ::marco::runtime::profiling::eulerForwardProfiler().stateVariables.stop(); \
  }

#define EULER_FORWARD_PROFILER_NONSTATEVAR_START                                   \
  if (::marco::runtime::simulation::getOptions().profiling) {                      \
    ::marco::runtime::profiling::eulerForwardProfiler().nonStateVariables.start(); \
  }

#define EULER_FORWARD_PROFILER_NONSTATEVAR_STOP                                   \
  if (::marco::runtime::simulation::getOptions().profiling) {                     \
    ::marco::runtime::profiling::eulerForwardProfiler().nonStateVariables.stop(); \
  }
// clang-format on

#else

#define EULER_FORWARD_PROFILER_DO_NOTHING static_assert(true);

#define EULER_FORWARD_PROFILER_STATEVAR_START EULER_FORWARD_PROFILER_DO_NOTHING
#define EULER_FORWARD_PROFILER_STATEVAR_STOP EULER_FORWARD_PROFILER_DO_NOTHING

#define EULER_FORWARD_PROFILER_NONSTATEVAR_START                               \
  EULER_FORWARD_PROFILER_DO_NOTHING
#define EULER_FORWARD_PROFILER_NONSTATEVAR_STOP                                \
  EULER_FORWARD_PROFILER_DO_NOTHING

#endif // MARCO_PROFILING

#endif // MARCO_RUNTIME_SOLVERS_EULERFORWARD_PROFILER_H
