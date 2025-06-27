
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(bi LANGUAGES CXX)

find_package(ROOT REQUIRED COMPONENTS Minuit Spectrum RooFit RooFitCore Gpad)
add_compile_options ( -Wall)

#executable will be created 1 folder before build                                                       
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../comparison_Bi_Li)

add_executable(Bi src/Bi.cxx)

target_link_libraries(Bi ${ROOT_LIBRARIES})

include_directories(${ROOT_INCLUDE_DIRS})
# end                   
