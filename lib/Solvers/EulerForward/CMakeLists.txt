# Options libraries.
set(OPTIONS_SRC
    Options.cpp)

marco_add_runtime_shared_library(OptionsEulerForward
    ${OPTIONS_SRC})

marco_add_runtime_static_library(OptionsEulerForward_static
    ${OPTIONS_SRC})

# Solver library.
marco_add_runtime_static_library(SolverEulerForward
    Profiler.cpp)

target_link_libraries(SolverEulerForward
    PUBLIC
    OptionsEulerForward)

if (MARCO_PROFILING)
    target_link_libraries(SolverEulerForward
        PUBLIC
        Profiling)
endif()
