# - Example project for accessing SuperNEMO RTD data (Raw Trigger Data)

cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(snfee_examples_read_rtd LANGUAGES CXX)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(STATUS "${PROJECT_NAME} requires an out-of-source build.")
  message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
  message(STATUS "  CMakeCache.txt")
  message(STATUS "  CMakeFiles")
  message(STATUS "Once these files are removed, create a separate directory")
  message(STATUS "and run CMake from there, pointing it to:")
  message(STATUS "  ${CMAKE_SOURCE_DIR}")
  message(FATAL_ERROR "in-source build detected")
endif()

# # - Load Custom Modules
# list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)
# message(STATUS "[info] CMAKE_MODULE_PATH = '${CMAKE_MODULE_PATH}'")


find_package(SNFrontEndElectronics)
message(STATUS  "SNFrontEndElectronics_INCLUDE_DIR=${SNFrontEndElectronics_INCLUDE_DIR}")
include_directories(${SNFrontEndElectronics_INCLUDE_DIRS})

add_compile_options ( -Wall) 


# LI system program
add_executable(Li_system Li_system.cxx)
target_link_libraries(Li_system SNFrontEndElectronics::snfee)


# LI system program
add_executable(charge_amplitude_energie charge_amplitude_energie.cxx)
target_link_libraries(charge_amplitude_energie SNFrontEndElectronics::snfee)


# end
