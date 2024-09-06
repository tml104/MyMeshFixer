cmake_minimum_required(VERSION 3.25)
 
project(OpenMesh-Example)
 
 # ���� C++ ��׼
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# ���� spdlog �İ�װ·��
set(SPDLOG_DIR "E:/spdlog")

find_package(OpenMesh HINTS "E:/OpenMesh")
find_package(spdlog REQUIRED PATHS ${SPDLOG_DIR}/lib/cmake/spdlog)
 
set (targetName MyOwnProject)
 
add_executable (${targetName} main.cpp)
target_compile_definitions(${targetName} PRIVATE _USE_MATH_DEFINES)
target_link_libraries(${targetName} PRIVATE OpenMeshCore OpenMeshTools)

# ���� spdlog ��
target_link_libraries(${targetName} PRIVATE spdlog::spdlog)