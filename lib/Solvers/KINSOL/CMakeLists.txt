# Options libraries.
set(OPTIONS_SRC
    Options.cpp)

marco_add_runtime_shared_library(OptionsKINSOL
    ${OPTIONS_SRC})

marco_add_runtime_static_library(OptionsKINSOL_static
    ${OPTIONS_SRC})

set(OPTIONS_LINKED_LIBS
    SUNDIALS::kinsol)

target_link_libraries(OptionsKINSOL
    PUBLIC
    ${OPTIONS_LINKED_LIBS})

target_link_libraries(OptionsKINSOL_static
    PUBLIC
    ${OPTIONS_LINKED_LIBS})

# Solver libraries.
set(SOLVER_SRC
    Instance.cpp
    Profiler.cpp)

marco_add_runtime_shared_library(SolverKINSOL
    ${SOLVER_SRC})

marco_add_runtime_static_library(SolverKINSOL_static
    ${SOLVER_SRC})

if (MARCO_ENABLE_SUNDIALS)
    set(SOLVER_LINKED_LIBS
        SUNDIALS::KLU
        SUNDIALS::kinsol
        SUNDIALS::nvecserial
        SUNDIALS::sunlinsolklu)

    target_link_libraries(SolverKINSOL
        PUBLIC
        ${SOLVER_LINKED_LIBS})

    target_link_libraries(SolverKINSOL_static
        PUBLIC
        ${SOLVER_LINKED_LIBS})
endif()

set(SOLVER_LINKED_LIBS
    OptionsKINSOL
    SolverSUNDIALS
    Modeling
    Multithreading
    Support)

target_link_libraries(SolverKINSOL
    PUBLIC
    ${SOLVER_LINKED_LIBS}
    OptionsSimulation)

target_link_libraries(SolverKINSOL_static
    PUBLIC
    ${SOLVER_LINKED_LIBS}
    OptionsSimulation_static)

if (MARCO_PROFILING)
    target_link_libraries(SolverKINSOL
        PUBLIC
        Profiling)

    target_link_libraries(SolverKINSOL_static
        PUBLIC
        Profiling)
endif()
