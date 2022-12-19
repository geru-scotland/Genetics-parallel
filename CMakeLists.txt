cmake_minimum_required(VERSION 3.23)
project(genetics C)

set(CMAKE_C_STANDARD 11)

include_directories(include)

add_subdirectory(src/serie/
        src/parallel/)

add_executable(genetics_s
        src/serie/fun_s.c
        src/serie/gengrupos_s.c)

add_executable(genetics_p
        src/parallel/fun_s.c
        src/parallel/gengrupos_s.c)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

target_link_libraries(genetics_s m)
target_link_libraries(genetics_p m)
