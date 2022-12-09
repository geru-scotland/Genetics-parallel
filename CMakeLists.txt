cmake_minimum_required(VERSION 3.23)
project(genetics C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(genetics
        defineg.h
        fun.h
        fun_s.c
        gengrupos_s.c)

target_link_libraries(genetics m)