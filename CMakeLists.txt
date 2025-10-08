cmake_minimum_required(VERSION 3.22.1)
project(linalg)
set(CMAKE_CXX_STANDARD 17)

add_library(${PROJECT_NAME} source/matrix.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC include/)