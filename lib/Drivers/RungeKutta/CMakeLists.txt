marco_add_runtime_static_library(DriverRungeKutta
    STATIC
    ../Driver.cpp
    CLI.cpp
    Driver.cpp)

target_link_libraries(DriverRungeKutta
    PUBLIC
    CLI)
