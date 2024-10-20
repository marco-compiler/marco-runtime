set(SRC
    Profiler.cpp
    Profiling.cpp
    Statistics.cpp
    Timer.cpp)

marco_add_runtime_shared_library(Profiling
    ${SRC})

marco_add_runtime_static_library(Profiling_static
    ${SRC})
