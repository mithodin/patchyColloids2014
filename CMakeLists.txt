cmake_minimum_required(VERSION 2.6)
project(patchy)

#set(CMAKE_BUILD_TYPE Debug)

include_directories("/Users/ltref/Documents/bachelor_project/libraries/include")
link_directories("/Users/ltref/Documents/bachelor_project/libraries/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/")
set(CMAKE_C_FLAGS "-msse2 -O4")
add_definitions(-DDSFMT_MEXP=19937 -DHAVE_SSE2)
add_executable(patchy_main dSFMT/dSFMT.c patchy_main.c load_config.c initialize.c distance.c colloid.c monte_carlo.c statistics.c threaded.c)
target_link_libraries(patchy_main config m pthread)

set_property( TARGET patchy_main APPEND_STRING PROPERTY COMPILE_FLAGS -Wall )
