cmake_minimum_required(VERSION 3.10)

# The EvoLearn program runs evolutionary simulations
# Copyright (C) 2025  Olof Leimar
# See Readme.md in top repo directory for copyright notice

# set the project name
project(EvoLearn VERSION 1.0 LANGUAGES CXX C)

# export compile commands to build directory
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set c++ standard
set(CMAKE_CXX_STANDARD 20)

# we need HighFive to handle hdf5 files
find_package(HighFive REQUIRED)

find_package(OpenMP)

# add the executable
add_executable(EvoLearn Evo.cpp EvoSim.cpp hdf5code.cpp Utils.cpp)

# link target
target_link_libraries(EvoLearn PUBLIC HighFive OpenMP::OpenMP_CXX)
