#include "PDNSolution_EP.hpp"

PDNSolution_EP::PDNSolution_EP(const class APart_Node * const &pNode,
					 const FEANode * const &fNode,
					 const class ALocal_NodalBC * const &locbc,
					 int type )
  : PDNSolution(pNode)
{
  switch (type)
  {
    case 0:
      Init_ZeroTemp(locbc);
      SYS_T::commPrint("===> Initial solution: Zero temperature for heat equation. \n");
      break;
    case 1:
      Init_OneTemp(locbc);
      SYS_T::commPrint("===> Initial solution: Interior 1.0, boundary 0.0, for heat equation. \n");
      break;
    case 2:
      Init_Partial(pNode, fNode, locbc);
      SYS_T::commPrint("===> Initial solution: -80 and 0 volts partially. \n");
      break;
    case 3:
      Init_Rest(locbc);
      SYS_T::commPrint("===> Initial solution: -80 milivolt overall. \n");
      break;            
    default:
      SYS_T::commPrint("ERROR: PDNSolution_EP: No such type of initial solution. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1); 
  }
}

PDNSolution_EP::~PDNSolution_EP()
{}

void PDNSolution_EP::Init_ZeroTemp( const class ALocal_NodalBC * const &LBC )
{
  VecSet(solution, 0.0);
  GhostUpdate();
}

void PDNSolution_EP::Init_OneTemp( const class ALocal_NodalBC * const &LBC )
{
  VecSet(solution, 1.0);
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  int num = LBC->get_Num_LD(0);
  int * index = new int [num];
  double * value = new double [num];

  for(int ii=0; ii<num; ++ii)
  {
    index[ii] = LBC->get_LDN(0, ii);
    value[ii] = 0.0;
  }
  VecSetValues(solution, num, index, value, INSERT_VALUES);
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  delete [] index; delete [] value;
}

//case 2
void PDNSolution_EP::Init_Partial( const class APart_Node * const &pNode,
					const FEANode * const &fNode,
					const class ALocal_NodalBC * const &LBC)
{
  //VecSet(solution, 1.0);
  int location;
  double value, x_coor, y_coor, z_coor; 
  const int nlocalnode = pNode -> get_nlocalnode();

  for(int ii=0; ii<nlocalnode; ++ii)
  {
    location = pNode -> get_node_loc(ii);
    x_coor = fNode ->  get_ctrlPts_x(ii);
    y_coor = fNode ->  get_ctrlPts_y(ii);
    z_coor = fNode ->  get_ctrlPts_z(ii);

    //set some nodes to 0 and some to -80
    if ( std::sqrt(std::pow(x_coor-104.9, 2.0) + std::pow(y_coor-301.5, 2.0) + std::pow(z_coor-248.24, 2.0)) <= 2.0 ) { //
      value = 0.0;
    }
    else {
      value = -80.0 ;
    }
    VecSetValue(solution, location, value, INSERT_VALUES);
  }
  
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();

  int num = LBC->get_Num_LD(0);
  int * index = new int [num];
  double * value_bc = new double [num];
  
  for(int ii=0; ii<num; ++ii)
  {
    index[ii] = LBC->get_LDN(0, ii);
    value_bc[ii] = 0.0;
  }
  VecSetValues(solution, num, index, value_bc, INSERT_VALUES);
  VecAssemblyBegin(solution);
  VecAssemblyEnd(solution);
  GhostUpdate();
  delete [] index; delete [] value_bc;
}

void PDNSolution_EP::Init_Rest( const class ALocal_NodalBC * const &LBC )
{
  VecSet(solution, -80.0);
  GhostUpdate();
}

int PDNSolution_EP::GetSize() const
{
  SYS_T::commPrint("GetSize implemented. \n");
  int size;
  VecGetSize(solution, &size);
  
  return size;
}

// EOF
