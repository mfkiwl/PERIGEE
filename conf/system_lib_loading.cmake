# Load appropriate system files into the configuration file

IF( $ENV{MACHINE_NAME} MATCHES "poincare")
  MESSAGE(STATUS "JL's Mac laptop banach")
  INCLUDE(${CMAKE_CURRENT_LIST_DIR}/poincare_PETSc_VTK.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "lax")
  MESSAGE(STATUS "Linux desktop Lax")
  INCLUDE(${CMAKE_CURRENT_LIST_DIR}/lax.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "sp2")
  MESSAGE(STATUS "TACC Stampede2")
  INCLUDE(${CMAKE_CURRENT_LIST_DIR}/stampede2_juliu.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "LS5")
  MESSAGE(STATUS "TACC Lonestar5")
  INCLUDE(${CMAKE_CURRENT_LIST_DIR}/ls5.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "chern")
  MESSAGE(STATUS "Stanford desktop chern")
  INCLUDE(${CMAKE_CURRENT_LIST_DIR}/stanford_bacon_v2018.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "comet")
  MESSAGE(STATUS "SDSC Comet")
  INCLUDE(${CMAKE_CURRENT_LIST_DIR}/comet_PETSc_v375_shared_VTK.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "sherlock")
  MESSAGE(STATUS "Stanford Sherlock")
  INCLUDE(${CMAKE_CURRENT_LIST_DIR}/stanford_sherlock.cmake)
elseif( $ENV{MACHINE_NAME} MATCHES "ingridxlan")
  MESSAGE(STATUS "Stanford ingridxlan")
  INCLUDE(${CMAKE_CURRENT_LIST_DIR}/marsden_linux_ingridxlan.cmake)
else($ENV{MACHINE_NAME} MATCHES "poincare")
  MESSAGE(STATUS "The system cannot be identified.")
endif( $ENV{MACHINE_NAME} MATCHES "poincare")

# End of the file
