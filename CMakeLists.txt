cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(example-multithreading)

add_subdirectory(example)

#############################################
# Unit tests
set(BUILD_TESTING ON)

if(BUILD_TESTING)
    add_subdirectory(third_party/googletest)
    enable_testing()
    add_subdirectory(tests)
endif()