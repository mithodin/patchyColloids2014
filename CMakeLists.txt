cmake_minimum_required(VERSION 2.6)
project(patchy)

#set(CMAKE_BUILD_TYPE Debug)

include_directories("/Users/ltref/Documents/bachelor_project/libraries/include")
link_directories("/Users/ltref/Documents/bachelor_project/libraries/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/")
add_executable(patchy_main patchy_main.c load_config.c initialize.c distance.c mt19937ar.c colloid.c monte_carlo.c statistics.c threaded.c)
target_link_libraries(patchy_main config m pthread)

set_property( TARGET patchy_main APPEND_STRING PROPERTY COMPILE_FLAGS -Wall )
