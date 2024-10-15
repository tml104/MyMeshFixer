cmake_minimum_required(VERSION 3.10)
 
project(OpenMesh-Example)
 
 # ���� C++ ��׼
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find OpenMesh
find_package(OpenMesh REQUIRED)
if(OpenMesh_FOUND)
    message(STATUS "OpenMesh found: ${OpenMesh_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR "OpenMesh not found")
endif()

# Find spdlog
find_package(spdlog REQUIRED)
if(spdlog_FOUND)
    message(STATUS "spdlog found: ${spdlog_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR "spdlog not found")
endif()

 
set (targetName MyMeshFixer)

# Include directories
include_directories(${OpenMesh_INCLUDE_DIRS} ${spdlog_INCLUDE_DIRS})
include_directories(./include)
 
add_executable (${targetName} main.cpp)
target_compile_definitions(${targetName} PRIVATE _USE_MATH_DEFINES)
target_link_libraries(${targetName} OpenMeshCore OpenMeshTools spdlog::spdlog)