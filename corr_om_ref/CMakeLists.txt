# - Example project for accessing SuperNEMO RTD data (Raw Trigger Data)

cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(corr_om_ref_exe LANGUAGES CXX)

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
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../corr_om_ref)

add_executable(charge_amplitude_energie_RHD src/charge_amplitude_energie_RHD.cxx)
target_link_libraries(charge_amplitude_energie_RHD SNFrontEndElectronics::snfee)

add_executable(charge_amplitude_energie src/charge_amplitude_energie.cxx)
target_link_libraries(charge_amplitude_energie SNFrontEndElectronics::snfee)

add_executable(corr_om_ref_exe src/new_Loop_fit_remake.cc)
target_link_libraries(corr_om_ref_exe ${ROOT_LIBRARIES})

add_executable(Modele_OM src/Modele_OM.cc)
target_link_libraries(Modele_OM ${ROOT_LIBRARIES})

include_directories(${ROOT_INCLUDE_DIRS})
# end
