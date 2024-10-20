marco_add_runtime_static_library(DriverKINSOL
    CLI.cpp)

target_link_libraries(DriverKINSOL
    PUBLIC
    CLI
    OptionsKINSOL)
