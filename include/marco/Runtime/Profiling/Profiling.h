#ifndef MARCO_RUNTIME_PROFILING_PROFILING_H
#define MARCO_RUNTIME_PROFILING_PROFILING_H

#include "marco/Runtime/Profiling/Profiler.h"
#include <memory>

namespace marco::runtime {
void registerProfiler(profiling::Profiler &profiler);
void registerProfiler(std::shared_ptr<profiling::Profiler> profiler);

void printProfilingStats();
} // namespace marco::runtime

#endif // MARCO_RUNTIME_PROFILING_PROFILING_H
