cmake_minimum_required(VERSION 3.14)
project(fasta_parser)

set(CMAKE_CXX_STANDARD 14)

set(COMPILE_FLAGS "-Wall -fPIC -O0 -std=c++11 -lpthread ")

add_executable(fasta_parser fasta_parser.cpp)