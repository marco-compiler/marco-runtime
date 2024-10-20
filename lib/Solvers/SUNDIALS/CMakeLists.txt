marco_add_runtime_shared_library(SolverSUNDIALS
    Instance.cpp
    Profiler.cpp)

marco_add_runtime_static_library(SolverSUNDIALS_static
    Instance.cpp
    Profiler.cpp)

set(LINKED_LIBS
    Modeling)

target_link_libraries(SolverSUNDIALS
    PUBLIC
    ${LINKED_LIBS})

target_link_libraries(SolverSUNDIALS_static
    PUBLIC
    ${LINKED_LIBS})

if (MARCO_ENABLE_SUNDIALS)
    set(LINKED_LIBS
        SUNDIALS::KLU
        SUNDIALS::nvecserial
        SUNDIALS::sunmatrixsparse)

    target_link_libraries(SolverSUNDIALS
        PUBLIC
        ${LINKED_LIBS})

    target_link_libraries(SolverSUNDIALS_static
        PUBLIC
        ${LINKED_LIBS})
endif()

if (MARCO_PROFILING)
    target_link_libraries(SolverSUNDIALS
        PUBLIC
        Profiling)

    target_link_libraries(SolverSUNDIALS_static
            PUBLIC
            Profiling)
endif()
