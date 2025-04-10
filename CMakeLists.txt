#------------------------------------------------------------------------------
# Project setup
#------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.9)

project(marco-runtime LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)

# Set the default build type if none was specified
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

include(GNUInstallDirs)
include(FetchContent)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# CMake configuration files.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")
include(AddMARCORuntime)

#------------------------------------------------------------------------------
# CLI
#------------------------------------------------------------------------------

option(MARCO_CLI "Enable the command-line interface" ON)

if (MARCO_CLI)
    add_compile_definitions(CLI_ENABLE)
endif()

#------------------------------------------------------------------------------
# Multithreading
#------------------------------------------------------------------------------

# Find the thread library.
find_package(Threads)

if (Threads_FOUND)
    add_compile_definitions(THREADS_ENABLE)
endif()

#------------------------------------------------------------------------------
# LLVM / MLIR configuration
#------------------------------------------------------------------------------

if (DEFINED LLVM_PATH)
    get_filename_component(LLVM_ABS_PATH ${LLVM_PATH} ABSOLUTE)
    set(LLVM_DIR "${LLVM_ABS_PATH}/lib/cmake/llvm")
    set(MLIR_DIR "${LLVM_ABS_PATH}/lib/cmake/mlir")
endif()

find_package(LLVM REQUIRED CONFIG)
find_package(MLIR REQUIRED CONFIG)

message(STATUS "Using LLVMConfig.cmake in: ${LLVM_DIR}")
message(STATUS "Using MLIRConfig.cmake in: ${MLIR_DIR}")

# These are a bit of a hack, because we're hijacking a lot of LLVM machinery.
set(LLVM_LIBRARY_OUTPUT_INTDIR ${CMAKE_BINARY_DIR}/lib)
set(LLVM_RUNTIME_OUTPUT_INTDIR ${CMAKE_BINARY_DIR}/bin)
set(MLIR_BINARY_DIR ${CMAKE_BINARY_DIR})

list(APPEND CMAKE_MODULE_PATH "${LLVM_CMAKE_DIR}")
list(APPEND CMAKE_MODULE_PATH "${MLIR_CMAKE_DIR}")

include(HandleLLVMOptions)
include(TableGen)
include(AddLLVM)
include(AddMLIR)

link_directories(${LLVM_BUILD_LIBRARY_DIR})
add_definitions(${LLVM_DEFINITIONS})

#------------------------------------------------------------------------------
# SUNDIALS libraries
#------------------------------------------------------------------------------

# Find the SUNDIALS libraries.
option(MARCO_ENABLE_SUNDIALS "Enable the usage of the SUNDIALS suite" ON)
option(MARCO_USE_BUILTIN_SUNDIALS "Use built-in SUNDIALS libraries (to be built separately, see README)" OFF)

if (MARCO_ENABLE_SUNDIALS)
    if (MARCO_USE_BUILTIN_SUNDIALS)
        set(SUNDIALS_DIR "${SUNDIALS_PATH}/${CMAKE_INSTALL_LIBDIR}/cmake/sundials")
    endif()

    find_package(SUNDIALS)

    if (SUNDIALS_FOUND)
        set(IDA_ENABLE ON)
    else()
        set(IDA_ENABLE OFF)
    endif()

    if (MARCO_USE_BUILTIN_SUNDIALS)
        set(SUNDIALS_LIBRARY_DIR ${SUNDIALS_PATH}/${CMAKE_INSTALL_LIBDIR})
    else()
        # Compute SUNDIALS installation prefix.
        # The CMake configuration file is inside lib/cmake/sundials, so we need
        # to go two folders up in the folders tree to get the path of the
        # libraries.
        get_filename_component(SUNDIALS_PATH "${SUNDIALS_DIR}" PATH)
        get_filename_component(SUNDIALS_PATH "${SUNDIALS_PATH}" PATH)

        set(SUNDIALS_LIBRARY_DIR ${SUNDIALS_PATH})
    endif()

    message(STATUS "SUNDIALS libraries path: ${SUNDIALS_LIBRARY_DIR}")
    link_directories("${SUNDIALS_LIBRARY_DIR}")
endif()

#------------------------------------------------------------------------------
# Runtime libraries
#------------------------------------------------------------------------------

# Profiling.
option(MARCO_PROFILING "Generate code for runtime profiling." ON)

# Add sanitizer
if (DEFINED MARCO_SANITIZER)
    message(STATUS "Using sanitizer: ${MARCO_SANITIZER}")
    add_compile_options("-fsanitize=${MARCO_SANITIZER}")
    add_link_options("-fsanitize=${MARCO_SANITIZER}")
endif ()

# Include directories.
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${MLIR_INCLUDE_DIRS})

# Set the RPATH to include the paths to external libraries
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Add a target that depends on all the libraries, in order to ease the
# declaration of tests.
add_custom_target(marco-runtime)

# Discover the libraries.
add_subdirectory(lib)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/marco
        DESTINATION include
        COMPONENT marco-headers
        FILES_MATCHING
        PATTERN "*.h")

#------------------------------------------------------------------------------
# Configuration export
#------------------------------------------------------------------------------

install(EXPORT MARCORuntimeTargets
        FILE MARCORuntimeTargets.cmake
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake"
        NAMESPACE MARCORuntime::)

add_subdirectory(cmake/modules)

#------------------------------------------------------------------------------
# Tests
#------------------------------------------------------------------------------

enable_testing()

if (NOT CMAKE_CROSSCOMPILING)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/f7902802f1a61140e188223fb6d1c95925cbec4a.zip
    )

    set(INSTALL_GTEST OFF CACHE INTERNAL "")

    # For Windows: prevent overriding the parent project's compiler / linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    add_subdirectory(unittest)
endif()
