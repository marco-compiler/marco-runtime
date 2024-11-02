# Declare the libraries holding the options.
set(OPTIONS_SRC
    Options.cpp)

marco_add_runtime_shared_library(OptionsSimulation
    ${OPTIONS_SRC})

marco_add_runtime_static_library(OptionsSimulation_static
    ${OPTIONS_SRC})

# Declare the library driving the simulation.
marco_add_runtime_static_library(Simulation
    CLI.cpp
    Profiler.cpp
    Runtime.cpp
    Scheduler.cpp)

target_link_libraries(Simulation
    PUBLIC
    CLI
    Modeling
    Multithreading
    OptionsSimulation)

if (MARCO_PROFILING)
    target_link_libraries(Simulation
        PUBLIC
        Profiling)
endif()

# Declare the library used to make the simulation executable.
marco_add_runtime_static_library(Starter
    Starter.cpp)

target_link_libraries(Starter PUBLIC Simulation)
