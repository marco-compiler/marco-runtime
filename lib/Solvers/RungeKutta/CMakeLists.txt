# Options libraries.
set(OPTIONS_SRC
    Options.cpp)

marco_add_runtime_shared_library(OptionsRungeKutta
    ${OPTIONS_SRC})

marco_add_runtime_static_library(OptionsRungeKutta_static
    ${OPTIONS_SRC})

# Solver libraries.
marco_add_runtime_static_library(SolverRungeKutta
    Profiler.cpp)

target_link_libraries(SolverRungeKutta
    PUBLIC
    OptionsRungeKutta)

if (MARCO_PROFILING)
    target_link_libraries(SolverRungeKutta
        PUBLIC
        Profiling)
endif()
