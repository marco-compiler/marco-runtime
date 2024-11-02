marco_add_runtime_static_library(Support
    BuiltInFunctions.cpp
    Math.cpp
    MemoryManagement.cpp
    Print.cpp
    UtilityFunctions.cpp)

if (MARCO_PROFILING)
    target_link_libraries(Support
        PUBLIC
        Profiling)
endif()
