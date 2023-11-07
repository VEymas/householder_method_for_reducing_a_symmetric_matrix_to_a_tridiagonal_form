cmake_minimum_required(VERSION 3.23)
project(vtm_prac)

set(CMAKE_CXX_STANDARD 17)

add_executable(vtm_prac 
        src/main.cpp
        src/householder_for_sym_matrix.cpp
        src/householder_for_sym_matrix.h
)
