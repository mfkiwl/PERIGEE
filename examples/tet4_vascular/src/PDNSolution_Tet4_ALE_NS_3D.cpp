#include "PDNSolution_Tet4_ALE_NS_3D.hpp"

PDNSolution_Tet4_ALE_NS_3D::PDNSolution_Tet4_ALE_NS_3D( 
    const APart_Node * const &pNode,
    const FEANode * const &fNode_ptr,
    const ALocal_Inflow_NodalBC * const &infbc,
    const int &type, const bool &isprint ) 
: PDNSolution( pNode ), is_print(isprint)
{
  if( pNode->get_dof() != 7 ) SYS_T::print_fatal("Error: PDNSolution_Tet4_ALE_NS_3D : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    case 1:
      Init_flow_parabolic( pNode, fNode_ptr, infbc );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Tet4_ALE_NS_3D: No such type of initional condition. \n");
      break;
  }
}


PDNSolution_Tet4_ALE_NS_3D::PDNSolution_Tet4_ALE_NS_3D( 
    const APart_Node * const &pNode,
    const int &type, const bool &isprint ) 
: PDNSolution( pNode ), is_print( isprint )
{
  if( pNode->get_dof() != 7 ) SYS_T::print_fatal("Error: PDNSolution_Tet4_ALE_NS_3D : the APart_Node gives wrong dof number. \n");

  switch(type)
  {
    case 0:
      Init_zero( pNode );
      break;
    default:
      SYS_T::print_fatal("Error: PDNSolution_Tet4_ALE_NS_3D: No such type of initional condition. \n");
      break;
  }
}


PDNSolution_Tet4_ALE_NS_3D::~PDNSolution_Tet4_ALE_NS_3D()
{}


void PDNSolution_Tet4_ALE_NS_3D::Init_zero(const APart_Node * const &pNode_ptr)
{
  int location[7];
  const double value[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const int nlocalnode = pNode_ptr->get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 7;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;
    location[6] = location[0] + 6;

    VecSetValues(solution, 7, location, value, INSERT_VALUES);
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();

  if( is_print )
  {
    SYS_T::commPrint("===> Initial solution: mesh_x = 0.0 \n");
    SYS_T::commPrint("                       mesh_y = 0.0 \n");
    SYS_T::commPrint("                       mesh_z = 0.0 \n");
    SYS_T::commPrint("                       pres   = 0.0 \n");
    SYS_T::commPrint("                       velo_x = 0.0 \n");
    SYS_T::commPrint("                       velo_y = 0.0 \n");
    SYS_T::commPrint("                       velo_z = 0.0 \n");
  }
}


void PDNSolution_Tet4_ALE_NS_3D::Init_flow_parabolic(
    const APart_Node * const &pNode_ptr,
    const FEANode * const &fNode_ptr,
    const ALocal_Inflow_NodalBC * const &infbc )
{
  int location[7];
  double value[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int nlocalnode = pNode_ptr->get_nlocalnode();

  // First enforce everything to be zero
  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location[0] = pNode_ptr->get_node_loc(ii) * 7;
    location[1] = location[0] + 1;
    location[2] = location[0] + 2;
    location[3] = location[0] + 3;
    location[4] = location[0] + 4;
    location[5] = location[0] + 5;
    location[6] = location[0] + 6;

    VecSetValues(solution, 7, location, value, INSERT_VALUES);
  }

  // Maximum speed formula is 
  //             2.0 x flow rate (1.0) / surface area
  // Here I use the unit flow rate, and the actual flow rate is adjusted
  // based on the CVFlowRate class.
  const double vmax = 2.0 / infbc->get_fularea();

  const double out_nx = infbc->get_outvec(0);
  const double out_ny = infbc->get_outvec(1);
  const double out_nz = infbc->get_outvec(2);

  if( infbc->get_Num_LD() > 0)
  {
    for(int ii=0; ii<nlocalnode; ++ii)
    {
      if( infbc->is_inLDN(pNode_ptr->get_node_loc(ii)) )
      {
        location[0] = pNode_ptr->get_node_loc(ii) * 7;
        location[1] = location[0] + 1;
        location[2] = location[0] + 2;
        location[3] = location[0] + 3;
        location[4] = location[0] + 4;
        location[5] = location[0] + 5;
        location[6] = location[0] + 6;

        const double x = fNode_ptr->get_ctrlPts_x(ii);
        const double y = fNode_ptr->get_ctrlPts_y(ii);
        const double z = fNode_ptr->get_ctrlPts_z(ii);

        const double r = infbc->get_radius(x,y,z);
        const double vel = vmax * (1.0 - r*r);

        // -1.0 is multiplied to make the flow direction inward
        value[4] = vel * (-1.0) * out_nx;
        value[5] = vel * (-1.0) * out_ny;
        value[6] = vel * (-1.0) * out_nz;

        VecSetValues(solution, 7, location, value, INSERT_VALUES);
      }
    }
  }

  VecAssemblyBegin(solution); VecAssemblyEnd(solution);
  GhostUpdate();

  if(is_print)
  {
    SYS_T::commPrint("===> Initial solution: pres   = 0.0 \n");
    SYS_T::commPrint("                       velo_x = parabolic \n");
    SYS_T::commPrint("                       velo_y = parabolic \n");
    SYS_T::commPrint("                       velo_z = parabolic \n");
    SYS_T::commPrint("                       flow rate 1.0 .\n");
    SYS_T::commPrint("                       max speed %e.\n", vmax);
    SYS_T::commPrint("                       active area is %e.\n", infbc->get_actarea() );
    SYS_T::commPrint("                       full area is %e.\n", infbc->get_fularea() );
    SYS_T::commPrint("                       direction [%e %e %e].\n", out_nx, out_ny, out_nz);
  }
}

// EOF
