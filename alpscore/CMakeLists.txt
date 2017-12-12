# This is a minimal CMake file that illustrates how to link your
# project to ALPSCore
# Several additional options are possible (or even recommended), 
# this tutorial only illustrates the minimal options

# we need at least cmake version 3.0
cmake_minimum_required(VERSION 3.0)

set (CMAKE_CXX_STANDARD 11)

# the project is called 'linking', and it is a C++ project
project(alpscoretest CXX)

# the project relies on the ALPSCore package. If ALPSCore is not found
# automatically, specify its location using:
#export ALPSCore_DIR=/usr/local/share/ALPSCore/
find_package(ALPSCore REQUIRED)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

find_package(Curses REQUIRED)
find_package(GSL REQUIRED)

# we only have one cpp file: it's called linking.cpp
add_executable(${PROJECT_NAME} alpscoretest.cpp)

# Use ALPSCore_LIBRARIES variable to link to ALPSCore 
target_link_libraries(${PROJECT_NAME} ${ALPSCore_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${CURSES_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${CMAKE_SOURCE_DIR}/../libangulon.a)
