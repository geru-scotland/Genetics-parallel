cmake_minimum_required(VERSION 3.23)
project(genetics C)

set(CMAKE_C_STANDARD 11)

include_directories(serie/include
                    parallel/include)
add_subdirectory(serie/src
                 parallel/src)

add_executable(genetics_s
        serie/src/fun_s.c
        serie/src/gengrupos_s.c)

add_executable(genetics_p
        parallel/src/fun_s.c
        parallel/src/gengrupos_s.c)

find_package(OpenMP)
if (OPENMP_FOUND)
    # La flag OpenMP_C_FLAGS contiene -Xclang en MacOS, fuerzo su eliminación
    # a la hora de pasarle las flags de OpenMP al compilador
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

target_link_libraries(genetics_s m)
target_link_libraries(genetics_p m)
