cmake_minimum_required(VERSION 3.12)
project(la_vector)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)

add_subdirectory(lib/sdsl-lite)

add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(${PROJECT_NAME} INTERFACE sdsl)
target_include_directories(${PROJECT_NAME} INTERFACE include)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} INTERFACE OpenMP::OpenMP_CXX)
endif ()

add_executable(example example.cpp)
target_link_libraries(example la_vector)

enable_testing()
add_subdirectory(test)