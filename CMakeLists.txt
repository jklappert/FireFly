cmake_minimum_required(VERSION 2.6)
project(firefly)

add_executable(firefly main.cpp)

install(TARGETS firefly RUNTIME DESTINATION bin)
