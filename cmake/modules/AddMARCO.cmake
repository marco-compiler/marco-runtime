include_guard()

# Convert the "friendly" name of a variable into its target name.
# For example: utils -> MARCOUtils
function(marco_canonize_library_name canonical_name name)
  # Get first letter and capitalize.
  string(SUBSTRING ${name} 0 1 first-letter)
  string(TOUPPER ${first-letter} first-letter)

  # Get the rest of the name.
  string(LENGTH ${name} length)
  execute_process(COMMAND expr ${length} - 1 OUTPUT_VARIABLE length)
  string(SUBSTRING ${name} 1 ${length} rest)

  # Compose macro-name with first letter uppercase.
  set(${canonical_name} "MARCO${first-letter}${rest}" PARENT_SCOPE)
endfunction()

# Declare a MARCO library
macro(marco_add_library name)
  marco_canonize_library_name(canonized_name ${name})
  set_property(GLOBAL APPEND PROPERTY MARCO_LIBS ${canonized_name})

  cmake_parse_arguments(ARG
      "SHARED"
      ""
      "ADDITIONAL_HEADERS"
      ${ARGN})
  set(srcs)

  if (srcs OR ARG_ADDITIONAL_HEADERS)
    set(srcs
        ADDITIONAL_HEADERS
        ${srcs}
        ${ARG_ADDITIONAL_HEADERS}) # It may contain unparsed unknown args.
  endif()

  if (ARG_SHARED)
    set(LIBTYPE SHARED)
  else()
    # llvm_add_library ignores BUILD_SHARED_LIBS if STATIC is explicitly set,
    # so we need to handle it here.
    if (BUILD_SHARED_LIBS)
      set(LIBTYPE SHARED OBJECT)
    else()
      set(LIBTYPE STATIC OBJECT)
    endif()

    set_property(GLOBAL APPEND PROPERTY MARCO_STATIC_LIBS ${name})
  endif()

  llvm_add_library(${name} ${LIBTYPE} OUTPUT_NAME ${canonized_name} ${ARG_UNPARSED_ARGUMENTS} ${srcs})
  add_library(marco::${name} ALIAS ${name})

  include(GNUInstallDirs)

  install(TARGETS ${name}
      COMPONENT ${name}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  install(TARGETS ${name} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endmacro()

# Declare a MARCO tool
macro(marco_add_tool name)
  add_llvm_executable(${name} ${ARGN})
  add_executable(marco::${name} ALIAS ${name})

  llvm_update_compile_flags(${name})

  include(GNUInstallDirs)
  install(TARGETS ${name} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endmacro()

# Convert the "friendly" names of MARCO libraries into the ones to be used for linking
function(marco_map_components_to_libnames out_libs)
  set(link_components ${ARGN})

  if(NOT MARCO_AVAILABLE_LIBS)
    # Inside MARCO itself available libs are in a global property.
    get_property(MARCO_AVAILABLE_LIBS GLOBAL PROPERTY MARCO_LIBS)
  endif()

  string(TOUPPER "${MARCO_AVAILABLE_LIBS}" capitalized_libs)

  foreach(c ${link_components})
    get_property(c_rename GLOBAL PROPERTY LLVM_COMPONENT_NAME_${c})

    if(c_rename)
      set(c ${c_rename})
    endif()

    # Canonize the component name
    marco_canonize_library_name(canonized_name ${c})
    string(TOUPPER "${canonized_name}" capitalized)
    list(FIND capitalized_libs ${capitalized} lib_idx)

    if(lib_idx LESS 0)
      marco_canonize_library_name(canonical_name ${c})
      list(APPEND expanded_components ${canonical_name})
    else()
      list(GET MARCO_AVAILABLE_LIBS ${lib_idx} canonical_lib)
      list(APPEND expanded_components ${canonical_lib})
    endif()
  endforeach(c)

  set(${out_libs} ${expanded_components} PARENT_SCOPE)
endfunction()

# Declare a MARCO unit test leveraging Google Test
macro(marco_add_unittest test_suite test_name)
  add_executable(${test_name} ${ARGN})
  target_compile_features(${test_name} PUBLIC cxx_std_17)

  include(GoogleTest)
  target_link_libraries(${test_name} PRIVATE gmock gtest_main)
  gtest_discover_tests(${test_name})
endmacro()

###
###  marco_link_llvm_libs
###
###  Same as llvm_config, but automatically turns on shared libs if the
###  LLVM we are linking with has them available. Everything is linked as
###  PUBLIC so that everything works with the way we use CMake.
###    "Shouldn't have this been done by LLVM's Cmake stuff already?" You bet!
###  "Does it do that then?" NO!!
###
###             ----->>>> DO NOT REMOVE THIS FUNCTION! <<<<----
###
###   EVEN IF IT SEEMS NOT NECESSARY ON ***YOUR*** MACHINE IT DOES NOT MEAN
###                     IT IS USELESS FOR EVERYBODY
###
###    ----->>>> NEVER USE llvm_config, USE THIS FUNCTION INSTEAD <<<<-----
###
###     IF YOU DON'T, THE BUILD WILL BREAK FOR EVERYBODY WHO HAS BOTH THE
###                STATIC AND DYNAMIC LLVM LIBRARIES AVAILABLE
###
function(marco_link_llvm_libs target)
  set(link_components ${ARGN})

  if(LLVM IN_LIST LLVM_AVAILABLE_LIBS)
    if (DEFINED link_components AND DEFINED LLVM_DYLIB_COMPONENTS)
      if("${LLVM_DYLIB_COMPONENTS}" STREQUAL "all")
        set(link_components "")
      else()
        list(REMOVE_ITEM link_components ${LLVM_DYLIB_COMPONENTS})
      endif()
    endif()

    target_link_libraries(${target} PUBLIC LLVM)
  endif()

  llvm_map_components_to_libnames(libs ${link_components})
  target_link_libraries(${target} PUBLIC ${libs})
endfunction(marco_link_llvm_libs)