cmake_minimum_required(VERSION 2.6)
project(patchy)

include_directories("/Users/ltref/Documents/bachelor_project/libraries/include")
link_directories("/Users/ltref/Documents/bachelor_project/libraries/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin/")
add_executable(patchy_main patchy_main.c load_config.c mt19937ar.c)
target_link_libraries(patchy_main config m)
