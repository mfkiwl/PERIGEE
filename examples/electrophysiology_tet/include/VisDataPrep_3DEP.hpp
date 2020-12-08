#ifndef VISDATAPREP_3DEP_HPP
#define VISDATAPREP_3DEP_HPP
// ==================================================================
// VisDataPrep_3DEP.hpp
// This is the data preparation for 3D Nonlinear Heat problem.
//
// Dec. 17 2013
// Author: Ju Liu, liujuy@gmail.com
// Modified: Oguz Ziya Tikenogullari, o.z.tikenogullari@gmail.com
// ==================================================================

#include "IVisDataPrep.hpp"

class VisDataPrep_3DEP : public IVisDataPrep
{
  public:
    VisDataPrep_3DEP();
    virtual ~VisDataPrep_3DEP();

    virtual void get_pointArray(
        const std::string solution_file_name,
        const std::string analysis_node_mapping_file,
        const std::string post_node_mapping_file,
        const APart_Node * const &nNode_ptr,
        const IAGlobal_Mesh_Info * const &gInfo_ptr,
        const int &input_dof,
        double ** &pointArrays ) const;
};
#endif
