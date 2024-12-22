# Declare the libraries holding the options.
set(OPTIONS_SRC
    Options.cpp)

marco_add_runtime_shared_library(OptionsSimulation
    ${OPTIONS_SRC})

marco_add_runtime_static_library(OptionsSimulation_static
    ${OPTIONS_SRC})

# Declare the library driving the simulation.
# Both the versions lead to static objects, but differ in the linked libraries.
# The '_static' version depends only on static libraries.
marco_add_runtime_static_library(Simulation
    CLI.cpp
    Profiler.cpp
    Runtime.cpp
    Scheduler.cpp)

target_compile_definitions(Simulation PUBLIC SHARED_DEPS)

marco_add_runtime_static_library(Simulation_static
    CLI.cpp
    Profiler.cpp
    Runtime.cpp
    Scheduler.cpp)

set(SIMULATION_LINKED_LIBS
    CLI
    Modeling
    Multithreading
    OptionsSimulation)

target_link_libraries(Simulation
    PUBLIC
    ${SIMULATION_LINKED_LIBS}
    DriverKINSOL)

target_link_libraries(Simulation_static
    PUBLIC
    ${SIMULATION_LINKED_LIBS})

if (MARCO_PROFILING)
    target_link_libraries(Simulation
        PUBLIC
        Profiling)

    target_link_libraries(Simulation_static
        PUBLIC
        Profiling)
endif()

# Declare the library used to make the simulation executable.
marco_add_runtime_static_library(Starter
    Starter.cpp)

target_link_libraries(Starter PUBLIC Simulation)
