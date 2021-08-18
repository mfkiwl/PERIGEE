#include "VisDataPrep_3DEP.hpp"

VisDataPrep_3DEP::VisDataPrep_3DEP()
{
  arrayCompSize = 1;
  arrayNames.push_back("Transmembrane V.");
  arraySizes.push_back(1);
}

VisDataPrep_3DEP::VisDataPrep_3DEP( const std::string array_name)
{
  arrayCompSize = 1;
  arrayNames.push_back(array_name);
  arraySizes.push_back(1);
}


VisDataPrep_3DEP::~VisDataPrep_3DEP()
{}


void VisDataPrep_3DEP::get_pointArray(
    const std::string solution_file_name,
    const std::string analysis_node_mapping_file,
    const std::string post_node_mapping_file,
    const APart_Node * const &nNode_ptr,
    const IAGlobal_Mesh_Info * const &gInfo_ptr,
    const int &input_dof,
    double ** &pointArrays ) const
{
  int ntotal = nNode_ptr->get_nlocghonode();
  
  // Read in the postprocessor's local solution from petsc binary
  // this constructor automatically handle the different index method
  PostVectSolution pvsolu(solution_file_name, analysis_node_mapping_file,
      post_node_mapping_file, nNode_ptr, gInfo_ptr, input_dof);

  for(int ii=0; ii<ntotal; ++ii)
    pointArrays[0][ii] = pvsolu.get_locsol(ii*input_dof);

}

//EOF
