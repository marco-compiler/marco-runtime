#ifndef MARCO_RUNTIME_SOLVERS_KINSOL_PROFILER_H
#define MARCO_RUNTIME_SOLVERS_KINSOL_PROFILER_H

#ifdef MARCO_PROFILING

#include "marco/Runtime/Profiling/Profiler.h"
#include "marco/Runtime/Profiling/Timer.h"
#include "marco/Runtime/Simulation/Options.h"
#include <mutex>

namespace marco::runtime::profiling {
class KINSOLProfiler : public Profiler {
public:
  KINSOLProfiler();

  void reset() override;

  void print() const override;

  void incrementResidualsCallCounter();

  void incrementPartialDerivativesCallCounter();

public:
  int64_t residualsCallCounter{0};
  Timer residualsTimer;
  int64_t partialDerivativesCallCounter{0};
  Timer partialDerivativesTimer;
  Timer copyVarsFromMARCOTimer;
  Timer copyVarsIntoMARCOTimer;

  mutable std::mutex mutex;
};

KINSOLProfiler &kinsolProfiler();
} // namespace marco::runtime::profiling

// clang-format off
#define KINSOL_PROFILER_IC_START                                                  \
  if (::marco::runtime::simulation::getOptions().profiling) {                     \
    ::marco::runtime::profiling::kinsolProfiler().initialConditionsTimer.start(); \
  }

#define KINSOL_PROFILER_IC_STOP                                                  \
  if (::marco::runtime::simulation::getOptions().profiling) {                    \
    ::marco::runtime::profiling::kinsolProfiler().initialConditionsTimer.stop(); \
  }

#define KINSOL_PROFILER_RESIDUALS_CALL_COUNTER_INCREMENT                           \
  if (::marco::runtime::simulation::getOptions().profiling) {                      \
    ::marco::runtime::profiling::kinsolProfiler().incrementResidualsCallCounter(); \
  }

#define KINSOL_PROFILER_RESIDUALS_START                                       \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::kinsolProfiler().residualsTimer.start();     \
  }

#define KINSOL_PROFILER_RESIDUALS_STOP                                        \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::kinsolProfiler().residualsTimer.stop();      \
  }

#define KINSOL_PROFILER_PARTIAL_DERIVATIVES_CALL_COUNTER_INCREMENT                          \
  if (::marco::runtime::simulation::getOptions().profiling) {                               \
    ::marco::runtime::profiling::kinsolProfiler().incrementPartialDerivativesCallCounter(); \
  }

#define KINSOL_PROFILER_PARTIAL_DERIVATIVES_START                                  \
  if (::marco::runtime::simulation::getOptions().profiling) {                      \
    ::marco::runtime::profiling::kinsolProfiler().partialDerivativesTimer.start(); \
  }

#define KINSOL_PROFILER_PARTIAL_DERIVATIVES_STOP                                  \
  if (::marco::runtime::simulation::getOptions().profiling) {                     \
    ::marco::runtime::profiling::kinsolProfiler().partialDerivativesTimer.stop(); \
  }

#define KINSOL_PROFILER_COPY_VARS_FROM_MARCO_START                                \
  if (::marco::runtime::simulation::getOptions().profiling) {                     \
    ::marco::runtime::profiling::kinsolProfiler().copyVarsFromMARCOTimer.start(); \
  }

#define KINSOL_PROFILER_COPY_VARS_FROM_MARCO_STOP                                \
  if (::marco::runtime::simulation::getOptions().profiling) {                    \
    ::marco::runtime::profiling::kinsolProfiler().copyVarsFromMARCOTimer.stop(); \
  }

#define KINSOL_PROFILER_COPY_VARS_INTO_MARCO_START                                \
  if (::marco::runtime::simulation::getOptions().profiling) {                     \
    ::marco::runtime::profiling::kinsolProfiler().copyVarsIntoMARCOTimer.start(); \
  }

#define KINSOL_PROFILER_COPY_VARS_INTO_MARCO_STOP                                \
  if (::marco::runtime::simulation::getOptions().profiling) {                    \
    ::marco::runtime::profiling::kinsolProfiler().copyVarsIntoMARCOTimer.stop(); \
  }
// clang-format on

#else

#define KINSOL_PROFILER_DO_NOTHING static_assert(true);

#define KINSOL_PROFILER_IC_START KINSOL_PROFILER_DO_NOTHING
#define KINSOL_PROFILER_IC_STOP KINSOL_PROFILER_DO_NOTHING

#define KINSOL_PROFILER_RESIDUALS_CALL_COUNTER_INCREMENT                       \
  KINSOL_PROFILER_DO_NOTHING

#define KINSOL_PROFILER_RESIDUALS_START KINSOL_PROFILER_DO_NOTHING
#define KINSOL_PROFILER_RESIDUALS_STOP KINSOL_PROFILER_DO_NOTHING

#define KINSOL_PROFILER_PARTIAL_DERIVATIVES_CALL_COUNTER_INCREMENT             \
  KINSOL_PROFILER_DO_NOTHING

#define KINSOL_PROFILER_PARTIAL_DERIVATIVES_START KINSOL_PROFILER_DO_NOTHING
#define KINSOL_PROFILER_PARTIAL_DERIVATIVES_STOP KINSOL_PROFILER_DO_NOTHING

#define KINSOL_PROFILER_COPY_VARS_FROM_MARCO_START KINSOL_PROFILER_DO_NOTHING
#define KINSOL_PROFILER_COPY_VARS_FROM_MARCO_STOP KINSOL_PROFILER_DO_NOTHING

#define KINSOL_PROFILER_COPY_VARS_INTO_MARCO_START KINSOL_PROFILER_DO_NOTHING
#define KINSOL_PROFILER_COPY_VARS_INTO_MARCO_STOP KINSOL_PROFILER_DO_NOTHING

#endif // MARCO_PROFILING

#endif // MARCO_RUNTIME_SOLVERS_KINSOL_PROFILER_H
