cmake_minimum_required(VERSION 3.22.1)
project(linalg)
set(CMAKE_CXX_STANDARD 17)

add_library(${PROJECT_NAME} source/matrix.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC include/)

find_package(Doxygen)

if(DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(docs
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM
    )

    message(STATUS "Doxygen found. Use 'make docs' to generate documentation.")
else()
    message(WARNING "Doxygen not found. Documentation will not be generated.")
endif()