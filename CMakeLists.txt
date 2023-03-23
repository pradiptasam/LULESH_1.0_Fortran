cmake_minimum_required (VERSION 2.8.11)

#add_subdirectory (src)

#include_directories(include)

#find_package(BLAS)
#if (NOT BLAS_LIBRARIES)
#  message(FATAL_ERROR "A required library with BLAS API not found.")
#else()
#  message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
#endif()


#set( WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src)

SET(MY_COMP_OPT "Intel")
IF (MY_COMP_OPT STREQUAL "Intel")
  MESSAGE( STATUS "** Compiling with Intel settings **" )
  SET( CMAKE_FC_COMPILER "ifort" )
# SET( CMAKE_CXX_COMPILER "icpc" )
ENDIF ()

project (lulesh Fortran)
set(${PROJECT_NAME}_CONFIGURATIONS lulesh)

#IF (MY_COMP_OPT STREQUAL "Intel")
#  MESSAGE( STATUS "** Compiling with Intel settings **" )
#  SET( CMAKE_CXX_FLAGS "-O3 -std=c++11 -fopenmp" )
##  SET( CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++11 -fopenmp" )
#  SET( CMAKE_CXX_FLAGS_DEBUG "-g" )
#ENDIF ()
#include(CMakeForceCompiler)
#CMAKE_FORCE_C_COMPILER(icc "Intel C Compiler")
#CMAKE_FORCE_CXX_COMPILER(icpc "Intel C++ Compiler")

#set(CMAKE_C_FLAGS "-Wall ${CMAKE_C_FLAGS}")
#set(CMAKE_Fortran_FLAGS "--warn -fpp -O5")
set(CMAKE_Fortran_FLAGS "-cpp -O5")
#set(CMAKE_Fortran_FLAGS "-O5")
#set(CMAKE_Fortran_FLAGS_DEBUG "-fbounds-check -check bounds -O5")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -check all -fpe0 -warn -traceback -debug extended -O5")
#set(CMAKE_EXE_LINKER_FLAGS "-qopenmp -lblas -llapack")
set(CMAKE_EXE_LINKER_FLAGS "-qopenmp")
#set(CMAKE_EXE_LINKER_FLAGS "-qopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -O5")

file(GLOB SOURCES "src/*.f90")

# The file(GLOB...) allows for wildcard additions:
#file(GLOB SOURCES "src/*/*.*pp")
#include_directories(${PROJECT_SOURCE_DIR}/src)

add_executable(test ${SOURCES})
target_link_libraries(test ${PROJECT_LINK_LIBS} )
