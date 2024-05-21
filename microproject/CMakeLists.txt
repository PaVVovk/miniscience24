cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

project(gmshdemo)
add_executable(cylinder cylinder.cpp)
add_executable(rect rect.cpp)

# find_package() не находит gmsh, поэтому приходится искать через find_library() и find_path()
find_library(GMSH_LIB gmsh PATH_SUFFIXES lib)
find_path(GMSH_INC "gmsh.h" PATH_SUFFIXES include gmsh include/gmsh)

if(GMSH_LIB AND GMSH_INC)
  target_link_libraries(cylinder ${GMSH_LIB})
  target_link_libraries(rect ${GMSH_LIB})
endif()


