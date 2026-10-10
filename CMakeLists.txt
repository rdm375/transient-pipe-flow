cmake_minimum_required(VERSION 3.20)

project(
    IsothermalPipeSim
    VERSION 0.12.0
    LANGUAGES CXX
)

include(CTest)

# Native C++20 simulator.
include(cmake/CppSimulator.cmake)
