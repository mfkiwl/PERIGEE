# Configuration setup for Qi-Ming of J. Liu
set(HOME_DIR /home/maad-liuj)

# ========================================================
# Specify the library locations
# ========================================================
set(VTK_DIR ${HOME_DIR}/lib/VTK-7.1.1-SHARED/lib/cmake/vtk-7.1)

set(PETSC_DIR ${HOME_DIR}/lib/petsc-3.15.3-opt)
set(PETSC_ARCH .)

set(HDF5_ROOT ${HOME_DIR}/lib/hdf5-1.12.2)

# ========================================================
# Setup the libraries
# ========================================================
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

find_package(VTK REQUIRED)
find_package(PETSc REQUIRED)
find_package(HDF5 REQUIRED)

include_directories(${VTK_INCLUDE_DIRS})
include_directories(${PETSC_INC})
include_directories(${HDF5_INCLUDE_DIRS})

set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${VTK_LIBRARIES})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_LIB})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${PETSC_METIS_LIB})
message(STATUS "Use METIS in PETSc: " ${PETSC_METIS_LIB})
set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${HDF5_LIBRARIES})

message(STATUS "External Libraries: " ${EXTRA_LINK_LIBS})

# ===================================================================
# Compiler options
# ===================================================================
set(CMAKE_C_COMPILER /opt/software/intel/2017/compilers_and_libraries_2017.8.262/linux/mpi/intel64/bin/mpicc)
set(CMAKE_CXX_COMPILER /opt/software/intel/2017/compilers_and_libraries_2017.8.262/linux/mpi/intel64/bin/mpiicpc)
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_FLAGS "-xHOST -Wall")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_VERBOSE_MAKEFILE OFF)

# EOF
