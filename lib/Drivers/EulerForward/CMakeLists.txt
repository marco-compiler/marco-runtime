marco_add_runtime_static_library(DriverEulerForward
    ../Driver.cpp
    CLI.cpp
    Driver.cpp)

target_link_libraries(DriverEulerForward
    PUBLIC
    CLI
    OptionsEulerForward)
