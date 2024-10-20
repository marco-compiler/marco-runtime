marco_add_runtime_static_library(DriverIDA
    ../Driver.cpp
    CLI.cpp
    Driver.cpp)

target_link_libraries(DriverIDA
    PUBLIC
    CLI
    OptionsIDA)
