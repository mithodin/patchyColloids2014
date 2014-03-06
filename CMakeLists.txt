cmake_minimum_required(VERSION 2.6)
project(patchy)

link_libraries(m)
add_executable(patchy_main patchy_main.c mt19937ar.c)
