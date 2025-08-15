# Fetch the Argh! project, which is used to handle command-line arguments.
include(FetchContent)

FetchContent_Declare(
    argh
    GIT_REPOSITORY https://github.com/adishavit/argh.git
    GIT_TAG c3f0d8c8a6dacb00df626b409248a34e3bcd15f5)

FetchContent_MakeAvailable(argh)

install(TARGETS argh
        EXPORT MARCORuntimeTargets)

# Declare the library.
marco_add_runtime_static_library(CLI
    Category.cpp
    CLI.cpp)

target_link_libraries(CLI PUBLIC argh)
