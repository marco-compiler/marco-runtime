set(SRC
    BuiltInFunctions.cpp
    Math.cpp
    MemoryManagement.cpp
    Print.cpp
    UtilityFunctions.cpp)

marco_add_runtime_shared_library(Support
    ${SRC})

marco_add_runtime_static_library(Support_static
    ${SRC})

target_link_libraries(Support
    PUBLIC
    OptionsSimulation)

target_link_libraries(Support_static
    PUBLIC
    OptionsSimulation_static)

if (MARCO_PROFILING)
    target_link_libraries(Support
        PUBLIC
        Profiling)

    target_link_libraries(Support_static
        PUBLIC
        Profiling_static)
endif()

marco_add_runtime_shared_library(JITSupport
        JIT_LIB
        BuiltInFunctions.cpp
        Math.cpp
        Print.cpp
        UtilityFunctions.cpp)
