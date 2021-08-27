PROJECT(perigee)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  cmake_policy(SET CMP0074 NEW)
ENDIF(COMMAND cmake_policy)

SET(perigee_VERSION_MAJOR 0)
SET(perigee_VERSION_MINOR 0)

IF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "set to release" FORCE)
ELSE( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  SET(CMAKE_BUILD_TYPE Debug CACHE STRING "default to debug" FORCE)
ENDIF( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

# Libraries location
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# binary location
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# print makefile details or not
SET(CMAKE_VERBOSE_MAKEFILE FALSE)

INCLUDE(../../conf/basic_variable_log.cmake)
INCLUDE(../../conf/system_lib_loading.cmake)
#INCLUDE(../../conf/setup_external_libs.cmake)

# ===============================================
# Set my own source code
# ===============================================
SET( perigee_source ${perigee_SOURCE_DIR}/../../src )
# 1. Preprocessor source cpp
SET( perigee_preprocess_lib_src 
  ${perigee_source}/System/Sys_Tools.cpp
  ${perigee_source}/System/Vec_Tools.cpp
  ${perigee_source}/Mesh/Gmsh_FileIO.cpp
  ${perigee_source}/Mesh/Tet_Tools.cpp
  ${perigee_source}/System/HDF5_Writer.cpp
  ${perigee_source}/System/HDF5_Reader.cpp
  ${perigee_source}/System/Math_Tools.cpp
  ${perigee_source}/System/Matrix_double_3by3_Array.cpp
  ${perigee_source}/Mesh/IEN_Line_P1.cpp
  ${perigee_source}/Mesh/IEN_Tetra_P1.cpp
  ${perigee_source}/Mesh/Mesh_Line_3D.cpp
  ${perigee_source}/Mesh/Mesh_Tet4.cpp
  ${perigee_source}/Mesh/ElemBC_3D_tet4.cpp
  ${perigee_source}/Mesh/ElemBC_3D_Line.cpp
  ${perigee_source}/Mesh/NBC_Partition_3D.cpp
  ${perigee_source}/Mesh/NBC_Partition_3D_inflow.cpp
  ${perigee_source}/Mesh/EBC_Partition_vtp.cpp
  ${perigee_source}/Mesh/NodalBC_3D_vtp.cpp
  ${perigee_source}/Mesh/NodalBC_Line_3D_vtp.cpp
  ${perigee_source}/Mesh/NodalBC_3D_inflow.cpp
  ${perigee_source}/Mesh/NodalBC_Line_3D_stimulus.cpp
  ${perigee_source}/Mesh/INodalBC.cpp
  ${perigee_source}/Mesh/Map_Node_Index.cpp
  ${perigee_source}/Mesh/Global_Part_METIS.cpp
  ${perigee_source}/Mesh/Global_Part_Serial.cpp
  ${perigee_source}/Mesh/Part_Tet.cpp
  ${perigee_source}/Mesh/Part_Line.cpp
  ${perigee_source}/Mesh/IPart.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Triangle.cpp
  ${perigee_source}/Element/FEAElement_Triangle3_3D_der0.cpp
  ${perigee_source}/Element/FEAElement_Triangle6_3D_der0.cpp
  ${perigee_SOURCE_DIR}/src/IEN_Mixed.cpp
  ${perigee_SOURCE_DIR}/src/Mesh_Mixed.cpp
  ${perigee_SOURCE_DIR}/src/Global_Part_METIS_Mixed.cpp
  ${perigee_SOURCE_DIR}/src/Part_Mixed_Mesh.cpp
  )

# 2. Anlysis source cpp
SET( perigee_analysis_lib_src
  ${perigee_source}/System/Sys_Tools.cpp
  ${perigee_source}/System/Vec_Tools.cpp
  ${perigee_source}/System/Math_Tools.cpp
  ${perigee_source}/System/Matrix_double_3by3_Array.cpp
  ${perigee_source}/System/Matrix_double_6by6_Array.cpp
  ${perigee_source}/System/HDF5_PartReader.cpp
  ${perigee_source}/System/HDF5_Writer.cpp
  ${perigee_source}/System/HDF5_Reader.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Triangle.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss_Tet.cpp
  ${perigee_source}/Mesh/QuadPts_Gauss.cpp
  ${perigee_source}/Analysis_Tool/AGlobal_Mesh_Info_FEM_3D.cpp  
  ${perigee_source}/Analysis_Tool/APart_Basic_Info.cpp
  ${perigee_source}/Analysis_Tool/ALocal_Elem.cpp
  ${perigee_source}/Analysis_Tool/ALocal_IEN.cpp
  ${perigee_source}/Analysis_Tool/ALocal_BC_3D.cpp
  ${perigee_source}/Analysis_Tool/AInt_Weight.cpp
  ${perigee_source}/Analysis_Tool/APart_Node.cpp
  ${perigee_source}/Analysis_Tool/ALocal_NodalBC.cpp
  ${perigee_source}/Analysis_Tool/ALocal_Inflow_NodalBC.cpp
  ${perigee_source}/Analysis_Tool/ALocal_EBC.cpp
  ${perigee_source}/Element/FEANode.cpp
  ${perigee_source}/Element/FEAElement_Triangle3_3D_der0.cpp
  #${perigee_source}/Element/FEAElement_Triangle6_3D_der0.cpp
  ${perigee_source}/Element/FEAElement_Tet4.cpp
  #${perigee_source}/Element/FEAElement_Tet10_v2.cpp
  ${perigee_source}/Element/FEAElement_Line2_3D_der1.cpp
  ${perigee_source}/Solver/PDNSolution.cpp
  ${perigee_source}/Solver/PDNTimeStep.cpp
  ${perigee_source}/Solver/TimeMethod_GenAlpha.cpp
  ${perigee_source}/Solver/PLinear_Solver_PETSc.cpp
  ${perigee_SOURCE_DIR}/src/ALocal_IEN_Mixed.cpp
  ${perigee_SOURCE_DIR}/src/PLocAssem_EP_3D.cpp
  ${perigee_SOURCE_DIR}/src/PDNSolution_EP.cpp
  ${perigee_SOURCE_DIR}/src/PTime_Solver_EP_OperatorSplit.cpp
  ${perigee_SOURCE_DIR}/src/PNonlinear_Solver_EP.cpp
  ${perigee_SOURCE_DIR}/src/PGAssem_EP.cpp
  ${perigee_SOURCE_DIR}/src/IonicModel_AP.cpp
  ${perigee_SOURCE_DIR}/src/IonicModel_TTP.cpp
  ${perigee_SOURCE_DIR}/src/IonicModel_FHN.cpp
  ${perigee_SOURCE_DIR}/src/IonicModel_Test.cpp
  ${perigee_SOURCE_DIR}/src/IonicModel_Purkinje.cpp
  ${perigee_SOURCE_DIR}/src/IonicModel.cpp
  ${perigee_SOURCE_DIR}/src/AGlobal_Mesh_Info_Mixed.cpp
  ${perigee_SOURCE_DIR}/src/ALocal_Elem_Fiber.cpp
  )

# 3. Postprocess cpp
SET( perigee_postprocess_lib_src
  ${perigee_source}/Analysis_Tool/AInt_Weight.cpp
  ${perigee_source}/Postproc_Tool/PostVectSolution.cpp
  ${perigee_source}/Postproc_Tool/Vis_Tools.cpp
  #${perigee_SOURCE_DIR}/src/VTK_Writer_EP_3D.cpp
  ${perigee_SOURCE_DIR}/src/VTK_Writer_EP_Mixed.cpp
  #${perigee_SOURCE_DIR}/src/Post_error_tet4_ns.cpp
  #${perigee_source}/Postproc_Tool/Interpolater.cpp
  ${perigee_source}/Element/FEAElement_Tet4.cpp
  ${perigee_source}/Mesh/QuadPts_vis_tet4.cpp
  ${perigee_source}/Mesh/QuadPts_vis.cpp
  ##
  #${perigee_source}/Postproc_Tool/PostVectSolution.cpp
  #${perigee_source}/Postproc_Tool/VTK_Writer.cpp
  ${perigee_SOURCE_DIR}/src/VisDataPrep_3DEP.cpp
  ${perigee_SOURCE_DIR}/src/Post_tools.cpp
  ${perigee_source}/Postproc_Tool/Interpolater.cpp
  #${perigee_SOURCE_DIR}/src/ALocal_IEN_Mixed.cpp
  )
# ===============================================
# Finish setting my own source code
# ===============================================

# INCLUDE MY OWN HEADER FILES
INCLUDE_DIRECTORIES( ../../include )
INCLUDE_DIRECTORIES( ./include )

MESSAGE(STATUS "=================== Misc setup ==========================")
MESSAGE( STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )
MESSAGE( STATUS "CMAKE_VERBOSE_MAKEFILE: " ${CMAKE_VERBOSE_MAKEFILE} )
MESSAGE( STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )
# -------------------------------------------------------------------
# MAKE MY OWN LIBRARIES
# 1. Preprocess libs
ADD_LIBRARY( perigee_preprocess ${perigee_preprocess_lib_src} )
TARGET_LINK_LIBRARIES( perigee_preprocess ${EXTRA_LINK_LIBS} )

# 2. Analysis libs
ADD_LIBRARY( perigee_analysis ${perigee_analysis_lib_src} )
TARGET_LINK_LIBRARIES( perigee_analysis ${EXTRA_LINK_LIBS} )

# 3. Postprocess lib
ADD_LIBRARY( perigee_postprocess ${perigee_postprocess_lib_src} )
TARGET_LINK_LIBRARIES( perigee_postprocess ${EXTRA_LINK_LIBS} perigee_analysis )
# -------------------------------------------------------------------

# Print compiler setup
INCLUDE(../../conf/message_compiler_setup.cmake)

# ===================================================================
# BUILD THE MAIN DRIVERS
# ===================================================================
# 1. Make executable
# 1.1 Preprocess
ADD_EXECUTABLE( gmshIO3d gmsh_process_3d.cpp)
ADD_EXECUTABLE( preprocess3d preprocess3d_main.cpp)

# 1.2 Analysis and tests
ADD_EXECUTABLE( analysis_3d_EP EP_3d_driver.cpp)

## 1.3 Post-process
#
ADD_EXECUTABLE( prepostproc3d Pre_postprocess3d.cpp)
#ADD_EXECUTABLE( post_compare_manu_3d post_compare_manu_3d.cpp)
ADD_EXECUTABLE( vis_3d_EP vis_3DEP.cpp)
ADD_EXECUTABLE( compute_ECG compute_ECG.cpp)
ADD_EXECUTABLE( activation_times activation_times.cpp)
#
## -------------------------------------------------------------------
# 2. Link libraries
# 2.1 Preprocess
TARGET_LINK_LIBRARIES( gmshIO3d perigee_preprocess )
TARGET_LINK_LIBRARIES( preprocess3d ${EXTRA_LINK_LIBS} perigee_preprocess )

# 2.2 Analysis and tests
TARGET_LINK_LIBRARIES( analysis_3d_EP ${EXTRA_LINK_LIBS} perigee_analysis )

## 2.3 Postprocess
TARGET_LINK_LIBRARIES( prepostproc3d ${EXTRA_LINK_LIBS} perigee_preprocess
  perigee_analysis )
#
#TARGET_LINK_LIBRARIES( post_compare_manu_3d ${EXTRA_LINK_LIBS} 
#  perigee_preprocess perigee_analysis perigee_postprocess )
#
TARGET_LINK_LIBRARIES( vis_3d_EP ${EXTRA_LINK_LIBS}
  perigee_preprocess perigee_analysis perigee_postprocess )

TARGET_LINK_LIBRARIES( compute_ECG ${EXTRA_LINK_LIBS}
  perigee_preprocess perigee_analysis perigee_postprocess )

TARGET_LINK_LIBRARIES( activation_times ${EXTRA_LINK_LIBS}
  perigee_preprocess perigee_analysis perigee_postprocess )  


#EOF
