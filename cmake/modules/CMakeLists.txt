include(CMakePackageConfigHelpers)

set(INSTALL_PATH "")
set(TARGETS_PATH "${CMAKE_INSTALL_LIBDIR}/cmake")

set(INCLUDE_PATH "include")
set(LIB_PATH ${CMAKE_INSTALL_LIBDIR})
set(SUNDIALS_LIB_PATH "dependencies/lib/sundials")

list(APPEND EXTRA_LINK_ARGS
        ${CMAKE_THREAD_LIBS_INIT})

configure_package_config_file(
    MARCORuntimeConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/MARCORuntimeConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/MARCORuntime"
    PATH_VARS
        TARGETS_PATH
        INCLUDE_PATH
        LIB_PATH
        SUNDIALS_LIBRARY_DIR)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/MARCORuntimeConfig.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/MARCORuntime")
