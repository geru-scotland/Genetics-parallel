cmake_minimum_required(VERSION 3.21)
project(genetics C)

set(CMAKE_C_STANDARD 11)

include_directories(include)
add_subdirectory(src)

add_executable(genetics
        src/fun_s.c
        src/gengrupos_s.c)

target_link_libraries(genetics m)