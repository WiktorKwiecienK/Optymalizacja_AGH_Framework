cmake_minimum_required(VERSION 4.0)
project(Optymalizacja_AGH_Framework)

set(CMAKE_CXX_STANDARD 14)

include_directories(Optymalizacja)

add_executable(Optymalizacja_AGH_Framework
        Optymalizacja/main.cpp
        Optymalizacja/matrix.cpp
        Optymalizacja/matrix.h
        Optymalizacja/ode_solver.cpp
        Optymalizacja/ode_solver.h
        Optymalizacja/opt_alg.cpp
        Optymalizacja/opt_alg.h
        Optymalizacja/solution.cpp
        Optymalizacja/solution.h
        Optymalizacja/user_funs.cpp
        Optymalizacja/user_funs.h)
