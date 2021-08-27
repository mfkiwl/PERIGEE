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
      SYS_T::commPrint("===> Initial solution: 0mV overall. \n");
      break;
    case 1:
      Init_OneTemp(locbc);
      SYS_T::commPrint("===> Initial solution: Interior 1.0, boundary 0.0. \n");
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

PDNSolution_EP::PDNSolution_EP(const class APart_Node * const &pNode,
			       const int &input_dof_num,
			       const FEANode * const &fNode,
			       const class ALocal_NodalBC * const &locbc,
			       const int type,
			       const ALocal_IEN_Mixed * const &lien_ptr,
			       const std::vector< IonicModel * > &ionicmodel_array)
  : PDNSolution(pNode, input_dof_num)
{
  switch (type){
  case 0:{
    
    VecSet(solution, 0.0);
    //int nlocghonode = pNode->get_nlocghonode();
    int nlocghonode = pNode->get_nlocalnode();
    std::vector<int> node_to_elem;
    PetscScalar * val = new PetscScalar[input_dof_num];
    PetscInt * idx = new PetscInt[input_dof_num];

    for (int count{ 0 }; count < nlocghonode; ++count){
      
      lien_ptr->get_node_to_elem(count, node_to_elem);

      for (int ii=0; ii<input_dof_num ; ++ii) val[ii] = 0.0;
      ionicmodel_array[node_to_elem[0]]->get_int_vars(val);

      int glo_node_idx= pNode->get_local_to_global(count);
      
      for (int ii=0; ii<input_dof_num ; ++ii){
	//idx[ii] = input_dof_num*count + ii;
	idx[ii] = input_dof_num* glo_node_idx + ii;
      }
      VecSetValues(solution,   // insert values into solution vec. 
      		   input_dof_num, // insert this many values
      		   idx,  // location of solution vector to inset values
      		   val,     // values to isnert
      		   INSERT_VALUES);// insert; not add. 
    }

    delete [] val; delete [] idx;
    
    VecAssemblyBegin(solution);
    VecAssemblyEnd(solution);
    GhostUpdate();
    SYS_T::commPrint
      ("===> Initial solution: internal variables (%d per node) are set to zero . \n",
       input_dof_num);
    break;
  }
    //
  default:{
    SYS_T::commPrint("ERROR: PDNSolution_EP: No such type of initial solution. \n");
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
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
    if ( std::sqrt(  std::pow(x_coor-(-107.110), 2.0)
		   + std::pow(y_coor-(-301.313), 2.0)
    		   + std::pow(z_coor-( 248.233), 2.0)  ) <= 10 ) {
    //if(x_coor < 20.0) {
    //if ((x_coor < 2.1) && (x_coor > 1.9)) {
    //if ( std::sqrt(  std::pow(x_coor-(3.0), 2.0)
    //		   + std::pow(y_coor-(2.0), 2.0)
    //		   + std::pow(z_coor-(0.0), 2.0)  ) <= 0.2 ) { 
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
