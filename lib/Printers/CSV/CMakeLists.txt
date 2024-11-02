# Options libraries.
set(OPTIONS_SRC
    Options.cpp)

marco_add_runtime_shared_library(OptionsCSVPrinter
    ${OPTIONS_SRC})

marco_add_runtime_static_library(OptionsCSVPrinter_static
    ${OPTIONS_SRC})

# Main library.
marco_add_runtime_static_library(PrinterCSV
    ../Printer.cpp
    CLI.cpp
    Printer.cpp
    Profiler.cpp)

target_link_libraries(PrinterCSV
    PUBLIC
    CLI
    OptionsCSVPrinter)

if (MARCO_PROFILING)
    target_link_libraries(PrinterCSV
        PUBLIC
        Profiling)
endif()
