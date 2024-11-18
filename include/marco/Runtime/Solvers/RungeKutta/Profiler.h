#ifndef MARCO_RUNTIME_SOLVERS_RUNGEKUTTA_PROFILER_H
#define MARCO_RUNTIME_SOLVERS_RUNGEKUTTA_PROFILER_H

#ifdef MARCO_PROFILING

#include "marco/Runtime/Profiling/Profiling.h"
#include "marco/Runtime/Profiling/Timer.h"
#include "marco/Runtime/Simulation/Options.h"
#include <mutex>

namespace marco::runtime::profiling {
class RungeKuttaProfiler : public Profiler {
public:
  RungeKuttaProfiler();

  void reset() override;

  void print() const override;

public:
  Timer stateVariables;
  Timer error;
  Timer variablesCopy;
  Timer nonStateVariables;

  mutable std::mutex mutex;
};

RungeKuttaProfiler &rungeKuttaProfiler();
} // namespace marco::runtime::profiling

// clang-format off
#define RUNGEKUTTA_PROFILER_STATEVAR_START                                    \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::rungeKuttaProfiler().stateVariables.start(); \
  }

#define RUNGEKUTTA_PROFILER_STATEVAR_STOP                                     \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::rungeKuttaProfiler().stateVariables.stop();  \
  }

#define RUNGEKUTTA_PROFILER_ERROR_START                                       \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::rungeKuttaProfiler().error.start();          \
  }

#define RUNGEKUTTA_PROFILER_ERROR_STOP                                        \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::rungeKuttaProfiler().error.stop();           \
  }

#define RUNGEKUTTA_PROFILER_VARCOPY_START                                     \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::rungeKuttaProfiler().variablesCopy.start();  \
  }

#define RUNGEKUTTA_PROFILER_VARCOPY_STOP                                      \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::rungeKuttaProfiler().variablesCopy.stop();   \
  }

#define RUNGEKUTTA_PROFILER_NONSTATEVAR_START                                    \
  if (::marco::runtime::simulation::getOptions().profiling) {                    \
    ::marco::runtime::profiling::rungeKuttaProfiler().nonStateVariables.start(); \
  }

#define RUNGEKUTTA_PROFILER_NONSTATEVAR_STOP                                    \
  if (::marco::runtime::simulation::getOptions().profiling) {                   \
    ::marco::runtime::profiling::rungeKuttaProfiler().nonStateVariables.stop(); \
  }
// clang-format on

#else

#define RUNGEKUTTA_PROFILER_DO_NOTHING static_assert(true);

#define RUNGEKUTTA_PROFILER_STATEVAR_START RUNGEKUTTA_PROFILER_DO_NOTHING
#define RUNGEKUTTA_PROFILER_STATEVAR_STOP RUNGEKUTTA_PROFILER_DO_NOTHING

#define RUNGEKUTTA_PROFILER_ERROR_START RUNGEKUTTA_PROFILER_DO_NOTHING
#define RUNGEKUTTA_PROFILER_ERROR_STOP RUNGEKUTTA_PROFILER_DO_NOTHING

#define RUNGEKUTTA_PROFILER_VARCOPY_START RUNGEKUTTA_PROFILER_DO_NOTHING
#define RUNGEKUTTA_PROFILER_VARCOPY_STOP RUNGEKUTTA_PROFILER_DO_NOTHING

#define RUNGEKUTTA_PROFILER_NONSTATEVAR_START RUNGEKUTTA_PROFILER_DO_NOTHING
#define RUNGEKUTTA_PROFILER_NONSTATEVAR_STOP RUNGEKUTTA_PROFILER_DO_NOTHING

#endif // MARCO_PROFILING

#endif // MARCO_RUNTIME_SOLVERS_RUNGEKUTTA_PROFILER_H
