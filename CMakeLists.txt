cmake_minimum_required(VERSION 3.0)
project(rma)

include_directories(include)
file(GLOB SOURCES "src/*.cpp")

set( CMAKE_CXX_STANDARD 11)
set( CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall" )

FIND_PACKAGE( PEBBL REQUIRED )
FIND_PACKAGE( MPI REQUIRED )

add_executable(boosting ${SOURCES})

install(TARGET boosting DESTINATION lib)
