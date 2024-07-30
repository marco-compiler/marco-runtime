marco_add_runtime_static_library(SolverRungeKutta
    STATIC
    Options.cpp
    Profiler.cpp)

target_link_libraries(SolverRungeKutta
    PUBLIC
    Profiling)
