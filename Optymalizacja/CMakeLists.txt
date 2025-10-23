cmake_minimum_required(VERSION 4.0)
project(Optymalizacja)

set(CMAKE_CXX_STANDARD 20)

add_executable(Optymalizacja main.cpp
        matrix.cpp
        opt_alg.cpp
        ode_solver.cpp
        solution.cpp
        user_funs.cpp)
