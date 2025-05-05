# Options libraries.
set(OPTIONS_SRC
    Options.cpp)

marco_add_runtime_shared_library(OptionsIDA
    ${OPTIONS_SRC})

marco_add_runtime_static_library(OptionsIDA_static
    ${OPTIONS_SRC})

set(OPTIONS_LINKED_LIBS
    SUNDIALS::ida)

target_link_libraries(OptionsIDA
    PUBLIC
    ${OPTIONS_LINKED_LIBS})

target_link_libraries(OptionsIDA_static
    PUBLIC
    ${OPTIONS_LINKED_LIBS})

# Solver libraries.
set(SOLVER_SRC
    Instance.cpp
    Profiler.cpp)

marco_add_runtime_shared_library(SolverIDA
    ${SOLVER_SRC})

marco_add_runtime_static_library(SolverIDA_static
    ${SOLVER_SRC})

if (MARCO_ENABLE_SUNDIALS)
    set(SOLVER_LINKED_LIBS
        SUNDIALS::KLU
        SUNDIALS::ida
        SUNDIALS::kinsol
        SUNDIALS::nvecserial
        SUNDIALS::sunlinsolklu)

    target_link_libraries(SolverIDA
        PUBLIC
        ${SOLVER_LINKED_LIBS})

    target_link_libraries(SolverIDA_static
        PUBLIC
        ${SOLVER_LINKED_LIBS})
endif()

set(SOLVER_LINKED_LIBS
    OptionsIDA
    SolverSUNDIALS
    Modeling
    Multithreading
    Support)

target_link_libraries(SolverIDA
    PUBLIC
    ${SOLVER_LINKED_LIBS}
    OptionsSimulation)

target_link_libraries(SolverIDA_static
    PUBLIC
    ${SOLVER_LINKED_LIBS}
    OptionsSimulation_static)

if (MARCO_PROFILING)
    target_link_libraries(SolverIDA
        PUBLIC
        Profiling)

    target_link_libraries(SolverIDA_static
        PUBLIC
        Profiling)
endif()
