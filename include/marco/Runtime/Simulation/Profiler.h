#ifndef MARCO_RUNTIME_SIMULATION_PROFILER_H
#define MARCO_RUNTIME_SIMULATION_PROFILER_H

#ifdef MARCO_PROFILING

#include "marco/Runtime/Profiling/Profiling.h"
#include "marco/Runtime/Profiling/Timer.h"
#include "marco/Runtime/Simulation/Options.h"
#include <mutex>

namespace marco::runtime::profiling {
class SimulationProfiler : public Profiler {
public:
  SimulationProfiler();

  void reset() override;

  void print() const override;

public:
  Timer initialization;
  Timer initialModel;
  Timer dynamicModel;
  Timer printing;

  mutable std::mutex mutex;
};

SimulationProfiler &simulationProfiler();
} // namespace marco::runtime::profiling

// clang-format off
#define SIMULATION_PROFILER_INIT_START                                        \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::simulationProfiler().initialization.start(); \
  }

#define SIMULATION_PROFILER_INIT_STOP                                         \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::simulationProfiler().initialization.stop();  \
  }

#define SIMULATION_PROFILER_INITIAL_MODEL_START                               \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::simulationProfiler().initialModel.start();   \
  }

#define SIMULATION_PROFILER_INITIAL_MODEL_STOP                                \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::simulationProfiler().initialModel.stop();    \
  }

#define SIMULATION_PROFILER_DYNAMIC_MODEL_START                               \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::simulationProfiler().dynamicModel.start();   \
  }

#define SIMULATION_PROFILER_DYNAMIC_MODEL_STOP                                \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::simulationProfiler().dynamicModel.stop();    \
  }

#define SIMULATION_PROFILER_PRINTING_START                                    \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::simulationProfiler().printing.start();       \
  }

#define SIMULATION_PROFILER_PRINTING_STOP                                     \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ::marco::runtime::profiling::simulationProfiler().printing.stop();        \
  }
// clang-format on

#else

#define SIMULATION_PROFILER_DO_NOTHING static_assert(true);

#define SIMULATION_PROFILER_INIT_START SIMULATION_PROFILER_DO_NOTHING
#define SIMULATION_PROFILER_INIT_STOP SIMULATION_PROFILER_DO_NOTHING

#define SIMULATION_PROFILER_INITIAL_MODEL_START SIMULATION_PROFILER_DO_NOTHING
#define SIMULATION_PROFILER_INITIAL_MODEL_STOP SIMULATION_PROFILER_DO_NOTHING

#define SIMULATION_PROFILER_DYNAMIC_MODEL_START SIMULATION_PROFILER_DO_NOTHING
#define SIMULATION_PROFILER_DYNAMIC_MODEL_STOP SIMULATION_PROFILER_DO_NOTHING

#define SIMULATION_PROFILER_PRINTING_START SIMULATION_PROFILER_DO_NOTHING
#define SIMULATION_PROFILER_PRINTING_STOP SIMULATION_PROFILER_DO_NOTHING

#endif // MARCO_PROFILING

#endif // MARCO_RUNTIME_SIMULATION_PROFILER_H
