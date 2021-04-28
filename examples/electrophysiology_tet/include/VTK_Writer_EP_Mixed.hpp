#ifndef VTK_WRITER_EP_MIXED_HPP
#define VTK_WRITER_EP_MIXED_HPP
// ==================================================================
// VTK_Writer_EP_Mixed.hpp
//
// This is the class that visualizes 3D Electophysiology results using
// tetrahedral elements.
//
// Date Created: June 27 2017
// Author: Ju Liu
// ==================================================================
#include "ALocal_Elem.hpp"
#include "ALocal_IEN_Mixed.hpp"
#include "IVisDataPrep.hpp"
#include "Interpolater.hpp"
#include "Vis_Tools.hpp"

#include "vtkIntArray.h"
#include "vtkCellData.h"

class VTK_Writer_EP_Mixed
{
public:
  VTK_Writer_EP_Mixed( const IAGlobal_Mesh_Info * const &GMIptr,
		       const std::string &epart_file);

  ~VTK_Writer_EP_Mixed();

  //    // --------------------------------------------------------------
  //    // This is the most original VTK writer that follows MBorden
  //    // style, that saves each sampling small tetrahedral as a 
  //    // separate element. This will lead to very huge number of pts,
  //    // very costly in memeory. This function will be used to test
  //    // and compare with future VTK writer implementations.
  //    // --------------------------------------------------------------
  //    void writeOutput(
  //        const FEANode * const &fnode_ptr,
  //        const ALocal_IEN * const &lien_ptr,
  //        const ALocal_Elem * const &lelem_ptr,
  //        const IVisDataPrep * const &vdata_ptr,
  //        FEAElement * const &elemptr,
  //        const IQuadPts * const &quad,
  //        const double * const * const &pointArrays,
  //        const int &rank, const int &size,
  //        const double &sol_time,
  //        const std::string &basename,
  //        const std::string &outputBName,
  //        const std::string &outputName,
  //        const bool &isXML );
  //
  //    // --------------------------------------------------------------
  //    // writeOutput_compact : this is a vtk writer that writes the grid
  //    //                       as a whole unstructured body, instead of
  //    //                       individual tet blocks. This function
  //    //                       significantly saves the memory size.
  //    // Note: On March 2, 2019, I add the Q-criterion visualization
  //    // into this function. This requires the VisDataPrep has the 
  //    // allocation of a scalar value named as Q-criterion/vortex 
  //    // identification.
  //    // --------------------------------------------------------------

  void writeOutput_compact(const IAGlobal_Mesh_Info * const &GMIptr,
			   const FEANode * const &fnode_ptr,
			   const ALocal_IEN_Mixed * const &lien_ptr,
			   const ALocal_Elem * const &lelem_ptr,
			   const IVisDataPrep * const &vdata_ptr,
			   std::vector<FEAElement *>  &elemArray,
			   const std::vector< IQuadPts * >  &quadArray,
			   const double * const * const &pointArrays,
			   const int &rank, const int &size,
			   const int &num_of_nodes,
			   const double &sol_time,
			   const std::string &basename,
			   const std::string &outputBName,
			   const std::string &outputName,
			   const bool &isXML );

private:
  int nLocBas;
  const int nElem;

  //Interpolater intep;

  std::vector<int> epart_map;

  int IEN_e[4];
  double ectrl_x[4];
  double ectrl_y[4];
  double ectrl_z[4];
};

#endif
