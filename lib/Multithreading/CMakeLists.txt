# Options libraries.
set(OPTIONS_SRC
    Options.cpp)

marco_add_runtime_shared_library(OptionsMultithreading
    ${OPTIONS_SRC})

marco_add_runtime_static_library(OptionsMultithreading_static
    ${OPTIONS_SRC})

# Main library.
marco_add_runtime_static_library(Multithreading
    Barrier.cpp
    CLI.cpp
    ThreadPool.cpp)

target_link_libraries(Multithreading
    PUBLIC
    CLI
    OptionsMultithreading)

if (Threads_FOUND)
    target_link_libraries(Multithreading
        PUBLIC
        Threads::Threads)
endif()
