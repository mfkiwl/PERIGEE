# =========================================================
# 1. VTK VARIABLES
# =========================================================
set(VTK_DIR /home/groups/amarsden/lib-perigee/VTK-7.1.1-shared/lib/cmake/vtk-7.1)

set(PETSC_DIR /home/groups/amarsden/lib-perigee/petsc-3.11.3-opt)
set(PETSC_ARCH .)

set(METIS_DIR /home/groups/amarsden/lib-perigee/metis-5.0.3)

set(HDF5_ROOT /home/groups/amarsden/lib-perigee/hdf5-1.8.16)

# ========================================================
# setup the libraries
# ========================================================

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

find_package(VTK REQUIRED)
find_package(PETSc REQUIRED)
find_package(HDF5 REQUIRED)

include_directories(${VTK_INCLUDE_DIRS})
include_directories(${PETSC_INC})

set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_LIB})

if(PETSC_METIS)
  set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_METIS_LIB})
  message(STATUS "Use METIS in PETSc: " ${PETSC_METIS_LIB})
else(PETSC_METIS)
  find_package(METIS)
  INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIRS})
  set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${METIS_LIBRARIES})
endif(PETSC_METIS)

include_directories(${HDF5_INCLUDE_DIRS})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LIBRARIES})

message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# ========================================================
# Compiler options 
# ========================================================
SET(CMAKE_C_COMPILER  /home/groups/amarsden/lib-perigee/petsc-3.11.3-opt/bin/mpicc)
SET(CMAKE_CXX_COMPILER /home/groups/amarsden/lib-perigee/petsc-3.11.3-opt/bin/mpicxx)
SET(CMAKE_CXX_STANDARD 11)

IF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_CXX_FLAGS "-O3 -Wall")
ELSE( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_CXX_FLAGS "-O0 -Wall")
ENDIF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )



# EOF
