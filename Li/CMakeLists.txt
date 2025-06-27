
# - Example project for accessing SuperNEMO RTD data (Raw Trigger Data)

cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(Li_system LANGUAGES CXX)

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

find_package(SNFrontEndElectronics 1.0.0 REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Minuit Spectrum RooFit RooFitCore)
message(STATUS  "SNFrontEndElectronics_INCLUDE_DIR=${SNFrontEndElectronics_INCLUDE_DIR}")
include_directories(${SNFrontEndElectronics_INCLUDE_DIRS})

add_compile_options ( -Wall) 

#executable will be created 1 folder before build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../Li)


add_executable(Li_system src/Li_system.cxx)
add_executable(Li src/new_histo_Li.cc)

target_link_libraries(Li_system SNFrontEndElectronics::snfee)
target_link_libraries(Li ${ROOT_LIBRARIES})

include_directories(${ROOT_INCLUDE_DIRS})
# end
