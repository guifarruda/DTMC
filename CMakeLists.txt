cmake_minimum_required(VERSION 2.6)
project(dtmc)

# Uncomment this to compile using the intel compiler
# set(CMAKE_C_COMPILER "icc")
# set(CMAKE_CXX_COMPILER "icpc")

set(CMAKE_CXX_FLAGS "-std=gnu++98 -O3 -pipe -Wall -funroll-loops -lgsl -lgslcblas -fopenmp -lpthread -ligraph -static")

add_executable(dtmc main.cpp Markov.cpp MonteCarlo.cpp)

target_link_libraries(dtmc gsl pthread igraph)

# install(TARGETS dtmc RUNTIME DESTINATION bin)
