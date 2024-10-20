#ifndef MARCO_RUNTIME_SOLVERS_RUNGEKUTTA_PROFILER_H
#define MARCO_RUNTIME_SOLVERS_RUNGEKUTTA_PROFILER_H

#ifdef MARCO_PROFILING

#include "marco/Runtime/Profiling/Profiling.h"
#include "marco/Runtime/Profiling/Timer.h"
#include <mutex>

namespace marco::runtime::profiling
{
class RungeKuttaProfiler : public Profiler
{
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

RungeKuttaProfiler& rungeKuttaProfiler();
} // namespace marco::runtime::profiling

#define RUNGEKUTTA_PROFILER_STATEVAR_START ::marco::runtime::profiling::rungeKuttaProfiler().stateVariables.start();
#define RUNGEKUTTA_PROFILER_STATEVAR_STOP ::marco::runtime::profiling::rungeKuttaProfiler().stateVariables.stop();

#define RUNGEKUTTA_PROFILER_ERROR_START ::marco::runtime::profiling::rungeKuttaProfiler().error.start();
#define RUNGEKUTTA_PROFILER_ERROR_STOP ::marco::runtime::profiling::rungeKuttaProfiler().error.stop();

#define RUNGEKUTTA_PROFILER_VARCOPY_START ::marco::runtime::profiling::rungeKuttaProfiler().variablesCopy.start();
#define RUNGEKUTTA_PROFILER_VARCOPY_STOP ::marco::runtime::profiling::rungeKuttaProfiler().variablesCopy.stop();

#define RUNGEKUTTA_PROFILER_NONSTATEVAR_START ::marco::runtime::profiling::rungeKuttaProfiler().nonStateVariables.start();
#define RUNGEKUTTA_PROFILER_NONSTATEVAR_STOP ::marco::runtime::profiling::rungeKuttaProfiler().nonStateVariables.stop();

#else

#define RUNGEKUTTA_PROFILER_STATEVAR_START static_assert(true)
#define RUNGEKUTTA_PROFILER_STATEVAR_STOP static_assert(true)

#define RUNGEKUTTA_PROFILER_ERROR_START static_assert(true)
#define RUNGEKUTTA_PROFILER_ERROR_STOP static_assert(true)

#define RUNGEKUTTA_PROFILER_VARCOPY_START static_assert(true)
#define RUNGEKUTTA_PROFILER_VARCOPY_STOP static_assert(true)

#define RUNGEKUTTA_PROFILER_NONSTATEVAR_START static_assert(true)
#define RUNGEKUTTA_PROFILER_NONSTATEVAR_STOP static_assert(true)

#endif // MARCO_PROFILING

#endif // MARCO_RUNTIME_SOLVERS_RUNGEKUTTA_PROFILER_H
