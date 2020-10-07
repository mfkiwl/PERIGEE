#ifndef VISDATAPREP_2DEP_HPP
#define VISDATAPREP_2DEP_HPP
// ==================================================================
// VisDataPrep_2DEP.hpp
// This is the data preparation for 2D Nonlinear Ep problem.
//
// Date: oct 5 2020
// ==================================================================

#include "IVisDataPrep.hpp"

class VisDataPrep_2DEP : public IVisDataPrep
{
  public:
    VisDataPrep_2DEP();
    virtual ~VisDataPrep_2DEP();

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
