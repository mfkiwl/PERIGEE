#include "PGAssem_EP.hpp"

PGAssem_EP::PGAssem_EP( const std::vector< IPLocAssem * > &locassem_array,
			const IAGlobal_Mesh_Info * const &agmi_ptr,
			const ALocal_Elem * const &alelem_ptr,
			const ALocal_IEN_Mixed * const &aien_ptr,
			const APart_Node * const &pnode_ptr,
			const ALocal_NodalBC * const &part_bc )
{
  //member variables  that are not assigned here and should be assigned in
  // element for loops :
  // nLocBas , row_index, col_index, local_a b c d, IEN_e, ectrl_xyz
  
  dof = locassem_array[0]->get_dof();

  const int nlocalnode = pnode_ptr->get_nlocalnode();
  const int nlocrow = dof * nlocalnode;
  const int nElem = alelem_ptr->get_nlocalele();

  int * dnnz = new int [nlocrow];
  int * onnz = new int [nlocrow];

  SYS_T::commPrint("===> Estimate sparse nonzero structure. \n");
  Get_dnz_onz(nElem, aien_ptr, pnode_ptr, part_bc, dnnz, onnz);
  
  Init_petsc_35(dnnz, onnz, nlocrow);

  delete [] dnnz; delete [] onnz;

  // PETSc 3.5.3 and 3.6.0 use the same function call for creating Mat object  
  SYS_T::commPrint("===> PETSc-3.6.0: MatCreateAIJ called. \n");

  Release_nonzero_err_str();

  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE. \n");

  // allocate the frequently used arrays in global assembly
  int nlgn = pnode_ptr->get_nlocghonode();
  array_a = new double [nlgn * dof];
  array_b = new double [nlgn * dof];
  array_c = new double [nlgn * dof];
  array_d = new double [nlgn * dof];

//  nLocBas = agmi_ptr->get_nLocBas();
//
//  row_index = new PetscInt [dof * nLocBas];
//  col_index = new PetscInt [dof * nLocBas];
//  
//  local_a = new double [dof * nLocBas];
//  local_b = new double [dof * nLocBas];
//  local_c = new double [dof * nLocBas];  
//  local_d = new double [dof * nLocBas];
//
//  IEN_e = new int [nLocBas];
//
//  ectrl_x = new double [nLocBas];
//  ectrl_y = new double [nLocBas];
//  ectrl_z = new double [nLocBas];
//  
}

PGAssem_EP::PGAssem_EP( const IPLocAssem * const &locassem_ptr,
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const APart_Node * const &pnode_ptr,
    const int &petsc_version_type )
{
  dof = locassem_ptr->get_dof();
  nLocBas = agmi_ptr->get_nLocBas();

  row_index = new PetscInt [dof * nLocBas];
  col_index = new PetscInt [dof * nLocBas];

  int sdegree = agmi_ptr->get_xdegree();
  int tdegree = agmi_ptr->get_ydegree();
  int udegree = agmi_ptr->get_zdegree();
  int nlocalnode = pnode_ptr->get_nlocalnode();

  //std::cout << "PGAssem constructor \n"
  //	    << "sdegree" << sdegree <<"\n"
  //	    << "tdegree" << tdegree <<"\n"
  //	    << "udegree" << udegree <<"\n"
  //	    << std::endl;

  int nz_prow;
  if((tdegree==0)&&(udegree==0)) {
    nz_prow = 2 * dof * (2*sdegree+1);
  }else if ((tdegree==0)&&(udegree==0)&&(sdegree==0)) {
    SYS_T::commPrint("element degrees in space are all zero");
  }else{
    nz_prow = dof * (2*sdegree+1) * (2*tdegree+1) * (2*udegree+1);
  }
  int nlocrow = dof * nlocalnode;

  switch(petsc_version_type)
  {
    case 0:
      Init_petsc_32(nz_prow, nlocrow);
      SYS_T::commPrint("===> PETSc-3.5.3: MatCreateAIJ called. \n");
      break;
    default:
      SYS_T::commPrint("Error: given petsc_version_type not implemented. \n");
      MPI_Abort(PETSC_COMM_WORLD, 1);
      break;
  }

  // allocate the frequently used arrays in global assembly
  int nlgn = pnode_ptr->get_nlocghonode(); //number of local ghost nodes
  array_a = new double [nlgn * dof];
  array_b = new double [nlgn * dof];
  array_c = new double [nlgn * dof];
  array_d = new double [nlgn * dof];  

  local_a = new double [dof * nLocBas];
  local_b = new double [dof * nLocBas];
  local_c = new double [dof * nLocBas];
  local_d = new double [dof * nLocBas];  

  IEN_e = new int [nLocBas];

  ectrl_x = new double [nLocBas];
  ectrl_y = new double [nLocBas];
  ectrl_z = new double [nLocBas];
}


PGAssem_EP::PGAssem_EP( const IPLocAssem * const &locassem_ptr,
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const APart_Node * const &pnode_ptr )
{
  dof = locassem_ptr->get_dof();
  nLocBas = agmi_ptr->get_nLocBas();

  row_index = new PetscInt [dof * nLocBas];
  col_index = new PetscInt [dof * nLocBas];

  const int sdegree = agmi_ptr->get_xdegree();
  const int tdegree = agmi_ptr->get_ydegree();
  const int udegree = agmi_ptr->get_zdegree();
  const int nlocalnode = pnode_ptr->get_nlocalnode();

  const int dnz = int ( 1.2 * dof * (2*sdegree+1) * (2*tdegree+1) * (2*udegree+1) );
  const int onz = int ( 1.2 * dof * (2*sdegree+1) * (2*tdegree+1) * (2*udegree+1) );
  const int nlocrow = dof * nlocalnode;

  Init_petsc_35(dnz, onz, nlocrow);
  SYS_T::commPrint("===> PETSc-3.5.3: MatCreateAIJ called. \n");

  Release_nonzero_err_str();

  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE. \n");
  
  // allocate the frequently used arrays in global assembly
  int nlgn = pnode_ptr->get_nlocghonode(); //number of local ghost nodes
  array_a = new double [nlgn * dof];
  array_b = new double [nlgn * dof];

  local_a = new double [dof * nLocBas];
  local_b = new double [dof * nLocBas];

  IEN_e = new int [nLocBas];

  ectrl_x = new double [nLocBas];
  ectrl_y = new double [nLocBas];
  ectrl_z = new double [nLocBas];
}


PGAssem_EP::PGAssem_EP( const IPLocAssem * const &locassem_ptr,
    const IAGlobal_Mesh_Info * const &agmi_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &aien_ptr,
    const APart_Node * const &pnode_ptr,
    const ALocal_NodalBC * const &part_bc )
{
  dof = locassem_ptr->get_dof();
  nLocBas = agmi_ptr->get_nLocBas();

  const int nlocalnode = pnode_ptr->get_nlocalnode();
  const int nlocrow = dof * nlocalnode;
  const int nElem = alelem_ptr->get_nlocalele();

  int * dnnz = new int [nlocrow];
  int * onnz = new int [nlocrow];

  SYS_T::commPrint("===> Estimate sparse nonzero structure. \n");
  Get_dnz_onz(nElem, aien_ptr, pnode_ptr, part_bc, dnnz, onnz);
  
  Init_petsc_35(dnnz, onnz, nlocrow);

  delete [] dnnz; delete [] onnz;

  // PETSc 3.5.3 and 3.6.0 use the same function call for creating Mat object  
  SYS_T::commPrint("===> PETSc-3.6.0: MatCreateAIJ called. \n");

  Release_nonzero_err_str();

  SYS_T::commPrint("===> MAT_NEW_NONZERO_ALLOCATION_ERR = FALSE. \n");

  // allocate the frequently used arrays in global assembly
  int nlgn = pnode_ptr->get_nlocghonode();
  array_a = new double [nlgn * dof];
  array_b = new double [nlgn * dof];

  row_index = new PetscInt [dof * nLocBas];
  col_index = new PetscInt [dof * nLocBas];
  
  local_a = new double [dof * nLocBas];
  local_b = new double [dof * nLocBas];

  IEN_e = new int [nLocBas];

  ectrl_x = new double [nLocBas];
  ectrl_y = new double [nLocBas];
  ectrl_z = new double [nLocBas];
}


PGAssem_EP::~PGAssem_EP()
{
  VecDestroy(&G);
  MatDestroy(&K);
  delete [] row_index;
  delete [] col_index;
  delete [] array_a;
  delete [] array_b;
  delete [] array_c;
  delete [] array_d;  
  delete [] local_a;
  delete [] local_b;
  delete [] local_c;
  delete [] local_d;  
  delete [] IEN_e;
  delete [] ectrl_x;
  delete [] ectrl_y;
  delete [] ectrl_z;
}

void PGAssem_EP::Init_petsc_32(const int &nonzero_per_row, const int &num_loc_row)
{
  MatCreateAIJ(PETSC_COMM_WORLD, num_loc_row, num_loc_row, PETSC_DECIDE,
      PETSC_DECIDE, nonzero_per_row, PETSC_NULL, nonzero_per_row, PETSC_NULL, &K);

  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, num_loc_row, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}


void PGAssem_EP::Init_petsc_35(const int &nonzero_per_row, const int &num_loc_row)
{
  MatCreateAIJ(PETSC_COMM_WORLD, num_loc_row, num_loc_row, PETSC_DECIDE,
      PETSC_DECIDE, nonzero_per_row, PETSC_NULL, nonzero_per_row, PETSC_NULL, &K);

  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, num_loc_row, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}


void PGAssem_EP::Init_petsc_35(const int &dnz, const int &onz, const int &num_loc_row)
{
  MatCreateAIJ(PETSC_COMM_WORLD, num_loc_row, num_loc_row, PETSC_DECIDE,
      PETSC_DECIDE, dnz, PETSC_NULL, onz, PETSC_NULL, &K);

  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, num_loc_row, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}


void PGAssem_EP::Init_petsc_35(const PetscInt * const &dnz,
    const PetscInt * const &onz, const int &num_loc_row)
{
  MatCreateAIJ(PETSC_COMM_WORLD, num_loc_row, num_loc_row, PETSC_DECIDE,
      PETSC_DECIDE, 0, dnz, 0, onz, &K);

  VecCreate(PETSC_COMM_WORLD, &G);
  VecSetSizes(G, num_loc_row, PETSC_DECIDE);

  VecSetFromOptions(G);
  VecSet(G, 0.0);
  VecSetOption(G, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}


void PGAssem_EP::EssBC_KG(const ALocal_NodalBC * const &bc_part, const int &field)
{
  // Check if dirichlet nodes exists within this partition
  // NOTE: we use ADD_VALUES here for matrix assembly, where we assumes that
  // the matrix is assemblyed with LID, which does nothing to the essential
  // boundary nodes, i.e., the boundary nodes' rows are zero rows.
  int local_dir = bc_part->get_Num_LD(field);
  if(local_dir > 0)
  {
    int row, col;
    for(int i=0; i<local_dir; ++i)
    {
      row = bc_part->get_LDN(field, i) * dof + field;
      col = row;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
    }
  }

  // Check if periodic slave nodes exists in this partition
  int local_sla = bc_part->get_Num_LPS(field);
  if(local_sla > 0)
  {
    int row, col;
    for(int i=0; i<local_sla; ++i)
    {
      row = bc_part->get_LPSN(field, i) * dof + field;
      col = bc_part->get_LPMN(field, i) * dof + field;
      MatSetValue(K, row, col, 1.0, ADD_VALUES);
      MatSetValue(K, row, row, -1.0, ADD_VALUES);
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}

void PGAssem_EP::EssBC_G(const ALocal_NodalBC * const &bc_part, const int &field)
{
  int local_dir = bc_part->get_Num_LD(field);
  int local_sla = bc_part->get_Num_LPS(field);
  int row;
  if( local_dir > 0 )
  {
    for(int ii=0; ii<local_dir; ++ii)
    {
      row = bc_part->get_LDN(field, ii) * dof + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }

  if( local_sla > 0 )
  {
    for(int ii=0; ii<local_sla; ++ii)
    {
      row = bc_part->get_LPSN(field, ii) * dof + field;
      VecSetValue(G, row, 0.0, INSERT_VALUES);
    }
  }
}


// // ! GetLocal is moved to hpp to make it inline.
// void PGAssem_EP::GetLocal(const double * const &array, const int * const &IEN,
//     double * const &local_array) const
// {
//   for(int ii=0; ii<nLocBas; ++ii)
//   {
//     for(int jj=0; jj<dof; ++jj)
//       local_array[ii*dof+jj] = array[IEN[ii]*dof + jj];
//   }
// }


void PGAssem_EP::Assem_nonzero_estimate(
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr, 
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const ALocal_NodalBC * const &bc_part )
{
  int nElem = alelem_ptr->get_nlocalele();
  int loc_dof = dof * nLocBas;
  int loc_index, lrow_index; // lcol_index;

  // loop over elements and insert 1.0 to every possible slots
  lassem_ptr->Assem_Estimate();

  for(int e=0; e<nElem; ++e)
  {
    for(int i=0; i<nLocBas; ++i)
    {
      loc_index  = lien_ptr->get_LIEN(e, i);
      for(int m=0; m<dof; ++m)
      {
        lrow_index = bc_part->get_LID( m, loc_index );
        row_index[dof * i + m] = dof * lrow_index + m;
        col_index[dof * i + m] = dof * lrow_index + m;
      }
    }
    MatSetValues(K, loc_dof, row_index, loc_dof, col_index,
        lassem_ptr->Tangent, ADD_VALUES);
  }

  for( int fie=0; fie<dof; ++fie)
    EssBC_KG( bc_part, fie );

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}

void PGAssem_EP::Assem_nonzero_estimate(
     const ALocal_Elem * const &alelem_ptr,
     std::vector< IPLocAssem * > &locassem_array,
     const ALocal_IEN_Mixed * const &lien_ptr,
     const APart_Node * const &node_ptr,
     const ALocal_NodalBC * const &bc_part )
{
  int nElem = alelem_ptr->get_nlocalele();
  int loc_dof;
  int loc_index, lrow_index; // lcol_index;
  
  // loop over elements and insert 1.0 to every possible slots
  if (locassem_array.size() != nElem) {
    SYS_T::print_fatal("Error: PGAssem_EP::Assem_nonzero_estimate locassem size != nElem.\n");
  }
  
  for(int e=0; e<nElem; ++e) {
  
    (locassem_array.at(e))->Assem_Estimate();
    loc_dof = dof * (lien_ptr->get_nLocBas_loc(e));

    row_index = new PetscInt [dof * (lien_ptr->get_nLocBas_loc(e))];
    col_index = new PetscInt [dof * (lien_ptr->get_nLocBas_loc(e))];
    
    for(int i=0; i<(lien_ptr->get_nLocBas_loc(e)); ++i) {
      
      loc_index  = lien_ptr->get_LIEN(e, i);

      for(int m=0; m<dof; ++m) {
	
	lrow_index = bc_part->get_LID( m, loc_index );
	row_index[dof * i + m] = dof * lrow_index + m;
	col_index[dof * i + m] = dof * lrow_index + m;
      }
    }

    MatSetValues(K, loc_dof, row_index, loc_dof, col_index,
    		 (locassem_array.at(e))->Tangent, ADD_VALUES);
    
    delete row_index; delete col_index;
    row_index = nullptr; col_index = nullptr; 
  }
  
  for( int fie=0; fie<dof; ++fie)
    EssBC_KG( bc_part, fie );
  
  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_EP::Get_dnz_onz( const int &nElem,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const ALocal_NodalBC * const &bc_part,
    PetscInt * const &dnz, PetscInt * const &onz ) const
{
  const int nlocalnode = node_ptr->get_nlocalnode();
  const int nnode = node_ptr->get_nlocghonode();

  std::vector<int> numLocNode;
  VEC_T::read_int_h5("NumLocalNode", "/", "nln", numLocNode);

  std::vector<unsigned int> nlist;
  nlist.clear();
  nlist.resize(numLocNode.size()+1);
  nlist[0] = 0;
  for(unsigned int ii=1; ii<=numLocNode.size(); ++ii)
    nlist[ii] = nlist[ii-1] + numLocNode[ii-1];

  // numlocNode is not needed anymore. nlist will provide the nodal 
  // partition info.
  VEC_T::clean(numLocNode);

  //for(int ii=0; ii<dof*nlocalnode; ++ii)
  //{
  //  dnz[ii] = 0; // Initialization: the num of nz are zero in each row.
  //  onz[ii] = 0;
  //}

  // This vector stores each row's diagonal col index and off-diagonal col index
  std::vector<int> interfun_d, interfun_o;

  // This MPI vector stores globally collected diagonal and off-diagonal col
  // number
  Vec vdnz, vonz;

  VecCreateMPI(PETSC_COMM_WORLD, dof * nlocalnode, PETSC_DETERMINE, &vdnz);
  VecCreateMPI(PETSC_COMM_WORLD, dof * nlocalnode, PETSC_DETERMINE, &vonz);

  VecSet(vdnz, 0.0);
  VecSet(vonz, 0.0);

  int row, col, ien_index, part_id_row, part_id_col;

  std::vector<int> elem4node;

  // loop for each row 
  for( int ii=0; ii<nnode; ++ii )
  {
    elem4node.clear();
    for(int ee=0; ee<nElem; ++ee)
    {
      if( lien_ptr->isNode_in_Elem(ee, ii) )
        elem4node.push_back(ee);
    }

    for( int mm=0; mm<dof; ++mm )
    {
      row = bc_part->get_LID( mm, ii );
      interfun_d.clear();
      interfun_o.clear();

      // only calculate for nondirichlet node
      if(row >= 0)
      {
        // based on nlist, find the row node's corresponding processor id.
        part_id_row = Get_part_id(nlist, row);
        for(unsigned int ei=0; ei<elem4node.size(); ++ei)
        {
          int ee = elem4node[ei];
          for(int kk=0; kk<nLocBas; ++kk)
          {
            ien_index = lien_ptr->get_LIEN(ee,kk);
            for( int nn=0; nn<dof; ++nn )
            {
              col = bc_part->get_LID( nn, ien_index );

              if( col>=0 )
              {
                // based on nlist, find the processor that the col node 
                // belong to. If they belong to the same processor, add 
                // to diagonal, otherwise, add to off-diagonal.
                part_id_col = Get_part_id(nlist, col);

                if( part_id_row == part_id_col )
                  interfun_d.push_back(dof*col + nn);
                else
                  interfun_o.push_back(dof*col + nn);
              }
            }
          }
        }
        // now the interfun_d and interfun_o have cached all the col that is
        // nonzero for this row, we sort them and remove repeated col number.
        VEC_T::sort_unique_resize(interfun_d);
        VEC_T::sort_unique_resize(interfun_o);

        // Finish calculating for each row, the real row index is
        // dof * row + mm
        VecSetValue(vdnz, dof*row+mm, double(interfun_d.size()), ADD_VALUES);
        VecSetValue(vonz, dof*row+mm, double(interfun_o.size()), ADD_VALUES);
      }
    }
  }

  VecAssemblyBegin(vdnz);
  VecAssemblyEnd(vdnz);

  VecAssemblyBegin(vonz);
  VecAssemblyEnd(vonz);

  // We need to handle the Dirichlet and Periodic Slave nodes
  for(int mm=0; mm<dof; ++mm)
  {
    int local_dir = bc_part->get_Num_LD(mm);
    if(local_dir > 0)
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        int row = bc_part->get_LDN(mm, ii) * dof + mm;
        VecSetValue(vdnz, row, 1.0, INSERT_VALUES);
        VecSetValue(vonz, row, 0.0, INSERT_VALUES);
      }
    }

    int local_sla = bc_part->get_Num_LPS(mm);

    if(local_sla > 0)
    {
      for(int ii=0; ii<local_sla; ++ii)
      {
        int row = bc_part->get_LPSN(mm,ii) * dof + mm;
        VecSetValue(vdnz, row, 2.0, INSERT_VALUES);
        VecSetValue(vonz, row, 2.0, INSERT_VALUES);
      }
    }
  }

  VecAssemblyBegin(vdnz);
  VecAssemblyEnd(vdnz);

  VecAssemblyBegin(vonz);
  VecAssemblyEnd(vonz);

  // Get the global vector size
  PetscInt vec_size;
  VecGetSize(vdnz, &vec_size);

  // Now we get the globally collected dnz and onz number
  PetscScalar * array_d;
  PetscScalar * array_o;

  VecGetArray(vdnz, &array_d);
  for(int ii=0; ii<dof*nlocalnode; ++ii)
  {
    dnz[ii] = int(array_d[ii]);
    // the estimator above may overestimate for periodic master nodes.
    // if the number of nonzeros is greater than the dof * nlocalnode
    // reduce it to full diagonal rows. Otherwise PETSc will throw an
    // error message.
    if(dnz[ii] > dof*nlocalnode)
      dnz[ii] = dof * nlocalnode;                                   
  }
  VecRestoreArray(vdnz, &array_d);

  const int max_onz = vec_size - dof * nlocalnode;

  VecGetArray(vonz, &array_o);
  for(int ii=0; ii<dof*nlocalnode; ++ii)
  {
    onz[ii] = int(array_o[ii]);

    if(onz[ii] > max_onz)
      onz[ii] = max_onz;
  }
  VecRestoreArray(vonz, &array_o);

  VecDestroy(&vdnz);
  VecDestroy(&vonz);
}


void PGAssem_EP::Get_dnz_onz( const int &nElem,
			      const ALocal_IEN_Mixed * const &lien_ptr,
			      const APart_Node * const &node_ptr,
			      const ALocal_NodalBC * const &bc_part,
			      PetscInt * const &dnz, PetscInt * const &onz ) const
{
  const int nlocalnode = node_ptr->get_nlocalnode();
  const int nnode = node_ptr->get_nlocghonode();

  std::vector<int> numLocNode;
  VEC_T::read_int_h5("NumLocalNode", "/", "nln", numLocNode);

  std::vector<unsigned int> nlist;
  nlist.clear();
  nlist.resize(numLocNode.size()+1);
  nlist[0] = 0;
  for(unsigned int ii=1; ii<=numLocNode.size(); ++ii)
    nlist[ii] = nlist[ii-1] + numLocNode[ii-1];

  // numlocNode is not needed anymore. nlist will provide the nodal 
  // partition info.
  VEC_T::clean(numLocNode);

  //for(int ii=0; ii<dof*nlocalnode; ++ii)
  //{
  //  dnz[ii] = 0; // Initialization: the num of nz are zero in each row.
  //  onz[ii] = 0;
  //}

  // This vector stores each row's diagonal col index and off-diagonal col index
  std::vector<int> interfun_d, interfun_o;

  // This MPI vector stores globally collected diagonal and off-diagonal col
  // number
  Vec vdnz, vonz;

  VecCreateMPI(PETSC_COMM_WORLD, dof * nlocalnode, PETSC_DETERMINE, &vdnz);
  VecCreateMPI(PETSC_COMM_WORLD, dof * nlocalnode, PETSC_DETERMINE, &vonz);

  VecSet(vdnz, 0.0);
  VecSet(vonz, 0.0);

  int row, col, ien_index, part_id_row, part_id_col;

  std::vector<int> elem4node;

  // loop for each row 
  for( int ii=0; ii<nnode; ++ii )
  {
    elem4node.clear();
    for(int ee=0; ee<nElem; ++ee)
    {
      if( lien_ptr->isNode_in_Elem(ee, ii) )
        elem4node.push_back(ee);
    }

    for( int mm=0; mm<dof; ++mm )
    {
      row = bc_part->get_LID( mm, ii );
      interfun_d.clear();
      interfun_o.clear();

      // only calculate for nondirichlet node
      if(row >= 0)
      {
        // based on nlist, find the row node's corresponding processor id.
        part_id_row = Get_part_id(nlist, row);
        for(unsigned int ei=0; ei<elem4node.size(); ++ei)
        { 
          int ee = elem4node[ei];
          for(int kk=0; kk<(lien_ptr->get_nLocBas_loc(ee)); ++kk)
          {
            ien_index = lien_ptr->get_LIEN(ee,kk);
            for( int nn=0; nn<dof; ++nn )
            {
              col = bc_part->get_LID( nn, ien_index );

              if( col>=0 )
              {
                // based on nlist, find the processor that the col node 
                // belong to. If they belong to the same processor, add 
                // to diagonal, otherwise, add to off-diagonal.
                part_id_col = Get_part_id(nlist, col);

                if( part_id_row == part_id_col )
                  interfun_d.push_back(dof*col + nn);
                else
                  interfun_o.push_back(dof*col + nn);
              }
            }
          }
        }
        // now the interfun_d and interfun_o have cached all the col that is
        // nonzero for this row, we sort them and remove repeated col number.
        VEC_T::sort_unique_resize(interfun_d);
        VEC_T::sort_unique_resize(interfun_o);

        // Finish calculating for each row, the real row index is
        // dof * row + mm
        VecSetValue(vdnz, dof*row+mm, double(interfun_d.size()), ADD_VALUES);
        VecSetValue(vonz, dof*row+mm, double(interfun_o.size()), ADD_VALUES);
      }
    }
  }

  VecAssemblyBegin(vdnz);
  VecAssemblyEnd(vdnz);

  VecAssemblyBegin(vonz);
  VecAssemblyEnd(vonz);

  // We need to handle the Dirichlet and Periodic Slave nodes
  for(int mm=0; mm<dof; ++mm)
  {
    int local_dir = bc_part->get_Num_LD(mm);
    if(local_dir > 0)
    {
      for(int ii=0; ii<local_dir; ++ii)
      {
        int row = bc_part->get_LDN(mm, ii) * dof + mm;
        VecSetValue(vdnz, row, 1.0, INSERT_VALUES);
        VecSetValue(vonz, row, 0.0, INSERT_VALUES);
      }
    }

    int local_sla = bc_part->get_Num_LPS(mm);

    if(local_sla > 0)
    {
      for(int ii=0; ii<local_sla; ++ii)
      {
        int row = bc_part->get_LPSN(mm,ii) * dof + mm;
        VecSetValue(vdnz, row, 2.0, INSERT_VALUES);
        VecSetValue(vonz, row, 2.0, INSERT_VALUES);
      }
    }
  }

  VecAssemblyBegin(vdnz);
  VecAssemblyEnd(vdnz);

  VecAssemblyBegin(vonz);
  VecAssemblyEnd(vonz);

  // Get the global vector size
  PetscInt vec_size;
  VecGetSize(vdnz, &vec_size);

  // Now we get the globally collected dnz and onz number
  PetscScalar * array_dnz;
  PetscScalar * array_onz;

  VecGetArray(vdnz, &array_dnz);
  for(int ii=0; ii<dof*nlocalnode; ++ii)
  {
    dnz[ii] = int(array_dnz[ii]);
    // the estimator above may overestimate for periodic master nodes.
    // if the number of nonzeros is greater than the dof * nlocalnode
    // reduce it to full diagonal rows. Otherwise PETSc will throw an
    // error message.
    if(dnz[ii] > dof*nlocalnode)
      dnz[ii] = dof * nlocalnode;                                   
  }
  VecRestoreArray(vdnz, &array_dnz);

  const int max_onz = vec_size - dof * nlocalnode;

  VecGetArray(vonz, &array_onz);
  for(int ii=0; ii<dof*nlocalnode; ++ii)
  {
    onz[ii] = int(array_onz[ii]);

    if(onz[ii] > max_onz)
      onz[ii] = max_onz;
  }
  VecRestoreArray(vonz, &array_onz);

  VecDestroy(&vdnz);
  VecDestroy(&vonz);
}


//assemble when history variables exist
void PGAssem_EP::Assem_tangent_residual(
    const PDNSolution * const &sol_a, //velo
    const PDNSolution * const &sol_b, //disp
    //const PDNSolution * const &sol_c, //pre_hist
    //PDNSolution * const &sol_d, //new hist
    const double &t_n,
    const double &dt,
    //const double &dt_ion,
    //const IonicModel * const &ionicmodel_ptr,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr, 
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const IQuadPts * const &quad,
    //const AInt_Weight * const &wei_ptr,
    std::vector<FEAElement*> &eptr_array,
    const ALocal_NodalBC * const &bc_part )
{
  int nElem = alelem_ptr->get_nlocalele();
  int loc_dof = dof * nLocBas;
  int loc_index, lrow_index; // lcol_index;

  //int node_num = sol_a->get_nlocal();
  //int node_num = node_ptr->get_nlocghonode();
    
  sol_a->GetLocalArray( array_a, node_ptr );//velo
  sol_b->GetLocalArray( array_b, node_ptr );//disp
  //sol_c->GetLocalArray( array_c, node_ptr );//pre_hist
  //sol_d->GetLocalArray( array_d, node_ptr );//new hist

  //int node_locgho =node_ptr->get_nlocghonode();
  //int dof=  lassem_ptr->get_dof();
  //double *array_iion = new double [ dof * node_locgho ];
  //double *array_dphi = new double [ dof * node_locgho ];
  
  //double r_new_tmp, r_old, new_soln, dPhi_Iion_tmp, Iion_tmp;
  //for (int count{ 0 }; count < node_num; ++count)
  //  {
  //    new_soln = array_b[count];      
  //    r_old    = array_c[count];
  //    
  //    ionicmodel_ptr->
  //	material_routine(r_old, dt_ion,new_soln,
  //			 Iion_tmp, dPhi_Iion_tmp,
  //			 r_new_tmp);
  //
  //    //use negative below, to be consistent with krishnamoorthi
  //    //2013 quadrature paper and goktepe 2009 paper.
  //    array_d   [count] = r_new_tmp;
  //    array_iion[count] = -Iion_tmp;
  //    array_dphi[count] = -dPhi_Iion_tmp;
  //  }

  for( int ee=0; ee<nElem; ++ee )
  {
    if( 1 )  //eptr_array[ee]->is_sizeNonzero()
    {
      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a); //velo
      GetLocal(array_b, IEN_e, local_b); //disp
      //GetLocal(array_c, IEN_e, local_c); //old_hist
      //GetLocal(array_d, IEN_e, local_d); //new hist

      //double *local_iion = new double [dof * nLocBas];
      //double *local_dphi = new double [dof * nLocBas];
      //GetLocal(array_iion, IEN_e, local_iion); 
      //GetLocal(array_dphi, IEN_e, local_dphi); 
  
      fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      lassem_ptr->
	Assem_Tangent_Residual(t_n, dt, local_a, local_b,
			       eptr_array[ee],//local_iion, local_dphi,
			       ectrl_x, ectrl_y, ectrl_z, quad);
  
      for(int i=0; i<nLocBas; ++i)
      {
        loc_index = IEN_e[i];
        //lcol_index = node_ptr->get_local_to_global( loc_index );
  
        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);
  
          row_index[dof * i + m] = dof * lrow_index + m;
          col_index[dof * i + m] = dof * lrow_index + m;
        }
      }
      
      MatSetValues(K, loc_dof, row_index, loc_dof, col_index,
          lassem_ptr->Tangent, ADD_VALUES);
  
      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);

      //VecSetValues(sol_d->solution, loc_dof,
      //		   row_index, local_d, INSERT_VALUES);

      //std::cout << "local assem residual" << std::endl ;
      //PetscScalarView(nLocBas, lassem_ptr->Residual ,PETSC_VIEWER_STDOUT_WORLD);
	
    }
  }
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
  //VecAssemblyBegin(sol_d->solution);
  //VecAssemblyEnd(sol_d->solution);

  for(int fie = 0; fie<dof; ++fie)
    EssBC_KG( bc_part, fie);

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  //std::cout<< "residual and tangent the end of PGassem tangent calc" << std::endl;
  //Print_G();
  //MatView(K, PETSC_VIEWER_STDOUT_WORLD);
}


//this version is for mixed mesh 
void PGAssem_EP::Assem_tangent_residual(const PDNSolution * const &sol_a,
					const PDNSolution * const &sol_b,
					//const PDNSolution * const &sol_c,//pre_hist
					//PDNSolution * const &sol_d,//new hist
					const double &t_n,
					const double &dt,
					//const double &dt_ion,
					//const IonicModel * const &ionicmodel_ptr,
					const ALocal_Elem * const &alelem_ptr,
					std::vector< IPLocAssem * > &lassem_array,
					const ALocal_IEN_Mixed * const &lien_ptr,
					const APart_Node * const &node_ptr,
					const FEANode * const &fnode_ptr,
					const std::vector< IQuadPts * > &quad_array,
					//const AInt_Weight * const &wei_ptr,
					std::vector<FEAElement*> &eptr_array,
					const ALocal_NodalBC * const &bc_part)
{
  std::cout << "start assem tangent res." << std::endl ;  
  int nElem = alelem_ptr->get_nlocalele();
  int loc_dof, nlocbas_ee;
  int loc_index, lrow_index; // lcol_index;
    
  sol_a->GetLocalArray( array_a, node_ptr );//velo
  sol_b->GetLocalArray( array_b, node_ptr );//disp

  for( int ee=0; ee<nElem; ++ee )
  {
    if( 1 )  //eptr_array[ee]->is_sizeNonzero()
    {
      nlocbas_ee = lien_ptr->get_nLocBas_loc(ee);
      loc_dof= dof * nlocbas_ee;
      IEN_e = new int [nlocbas_ee];
      local_a = new double [dof * nlocbas_ee];
      local_b = new double [dof * nlocbas_ee];
      ectrl_x = new double [nlocbas_ee];
      ectrl_y = new double [nlocbas_ee];
      ectrl_z = new double [nlocbas_ee];
      row_index = new PetscInt [dof * nlocbas_ee];
      col_index = new PetscInt [dof * nlocbas_ee];
    
      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a, nlocbas_ee); //velo
      GetLocal(array_b, IEN_e, local_b, nlocbas_ee); //disp

  
      fnode_ptr->get_ctrlPts_xyz(nlocbas_ee, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      std::cout << "local_a : " <<"\n";
      for (int i = 0; i < nlocbas_ee; ++i) 
	std::cout << local_a[i] << "\t";
      std::cout << std::endl;

      fnode_ptr->get_ctrlPts_xyz(nlocbas_ee, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      std::cout << "element " << ee<< " coords : " <<"\n"
		<< "x \t y \t z \n";
      for (int i=0; i<nlocbas_ee; ++i){
	std::cout << ectrl_x[i] << "\t" << ectrl_y[i] <<"\t"<< ectrl_z[i] << "\n";
      }

      (lassem_array.at(ee))->
	Assem_Tangent_Residual(t_n, dt, local_a, local_b,
			       eptr_array[ee],//local_iion, local_dphi,
			       ectrl_x, ectrl_y, ectrl_z, quad_array.at(ee));
  
      for(int i=0; i<nlocbas_ee; ++i)
      {
        loc_index = IEN_e[i];
  
        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);
  
          row_index[dof * i + m] = dof * lrow_index + m;
          col_index[dof * i + m] = dof * lrow_index + m;
        }
      }

      //set residual values from element to global
      MatSetValues(K, loc_dof, row_index, loc_dof, col_index,
		   (lassem_array.at(ee))->Tangent, ADD_VALUES);
  
      VecSetValues(G, loc_dof, row_index, (lassem_array.at(ee))->Residual, ADD_VALUES);

      std::cout << "local assem" << std::endl ;
      PetscScalarView(nlocbas_ee, (lassem_array.at(ee))->Residual ,PETSC_VIEWER_STDOUT_WORLD);
      std::cout << "local tangent " << std::endl ;
      PetscScalarView(nlocbas_ee*nlocbas_ee, (lassem_array.at(ee))->Tangent,PETSC_VIEWER_STDOUT_WORLD);
    }
  }
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  for(int fie = 0; fie<dof; ++fie)
    EssBC_KG( bc_part, fie);

  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

  std::cout<< "residual and tangent the end of PGassem tangent calc" << std::endl;
  Print_G();
  MatView(K, PETSC_VIEWER_STDOUT_WORLD);
  std::cout << "end assem tangent-residual ." << std::endl ;  
}


//assemble when history variables exist
void PGAssem_EP::Assem_residual(
    const PDNSolution * const &sol_a, //velo
    const PDNSolution * const &sol_b, //disp
    //const PDNSolution * const &sol_c, //pre_hist
    //PDNSolution * const &sol_d, //new hist
    const double &t_n,
    const double &dt,
    //const double &dt_ion,
    //const IonicModel * const &ionicmodel_ptr,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr, 
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const IQuadPts * const &quad,
    //const AInt_Weight * const &wei_ptr,
    std::vector<FEAElement*> &eptr_array,
    const ALocal_NodalBC * const &bc_part )
{
  int nElem = alelem_ptr->get_nlocalele();
  int loc_dof = dof * nLocBas;
  int loc_index, lrow_index;

  //int node_num = sol_a->get_nlocal();
  //int node_num = node_ptr->get_nlocghonode();
  
  sol_a->GetLocalArray( array_a, node_ptr );//velo
  sol_b->GetLocalArray( array_b, node_ptr );//disp
  //sol_c->GetLocalArray( array_c, node_ptr );//pre_hist
  //sol_d->GetLocalArray( array_d, node_ptr );//new hist

  //int node_locgho =node_ptr->get_nlocghonode();
  //int dof=  lassem_ptr->get_dof();
  //double *array_iion = new double [ dof * node_locgho ];
  //double *array_dphi = new double [ dof * node_locgho ];
  //
  //double r_new_tmp, r_old, new_soln, dPhi_Iion_tmp, Iion_tmp;
  //for (int count{ 0 }; count < node_num; ++count)
  //  {
  //    new_soln = array_b[count];      
  //    r_old    = array_c[count];
  //    
  //    ionicmodel_ptr->
  //	material_routine(r_old, dt_ion,new_soln,
  //			 Iion_tmp, dPhi_Iion_tmp,
  //			 r_new_tmp);
  //
  //    //use negative below, to be consistent with krishnamoorthi
  //    //2013 quadrature paper and goktepe 2009 paper.
  //    array_d   [count] = r_new_tmp;
  //    array_iion[count] = -Iion_tmp;
  //    array_dphi[count] = -dPhi_Iion_tmp;
  //  }

  for( int ee=0; ee<nElem; ++ee )
  {
    if( 1 ) //eptr_array[ee]->is_sizeNonzero()
    {
      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a); //velo
      GetLocal(array_b, IEN_e, local_b); //disp
      //GetLocal(array_c, IEN_e, local_c); //old_hist
      //GetLocal(array_d, IEN_e, local_d); //new hist

      //double *local_iion = new double [dof * nLocBas];
      //double *local_dphi = new double [dof * nLocBas];
      //GetLocal(array_iion, IEN_e, local_iion); 
      //GetLocal(array_dphi, IEN_e, local_dphi);
      
      fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      lassem_ptr->Assem_Residual(t_n, dt, local_a, local_b,
				 eptr_array[ee], //local_iion, local_dphi,
				 ectrl_x, ectrl_y, ectrl_z, quad);

      for(int i=0; i<nLocBas; ++i)
      {
        loc_index = IEN_e[i];

        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);
          row_index[dof * i + m] = dof * lrow_index + m;
        }
      }
      //set residual values from element to global
      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
      ////set hist values from element to global
      //VecSetValues(sol_d->solution, loc_dof,
      //              row_index, local_d, INSERT_VALUES);
    }
  }

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
  //VecAssemblyBegin(sol_d->solution);
  //VecAssemblyEnd(sol_d->solution);

  //std::cout<< "residual at the end of PGassem residual calc" << std::endl;
  //Print_G();
}


//this version is for mixed mesh 
void PGAssem_EP::Assem_residual(
    const PDNSolution * const &sol_a, //velo
    const PDNSolution * const &sol_b, //disp
    const double &t_n,
    const double &dt,
    const ALocal_Elem * const &alelem_ptr,
    std::vector< IPLocAssem * > &lassem_array,
    const ALocal_IEN_Mixed * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const std::vector< IQuadPts * > &quad_array,
    std::vector<FEAElement*> &eptr_array,
    const ALocal_NodalBC * const &bc_part )
{
  std::cout << "start assem  residual" << std::endl ;  
  int nElem = alelem_ptr->get_nlocalele();
  int loc_dof, nlocbas_ee;
  int loc_index, lrow_index;

  sol_a->GetLocalArray( array_a, node_ptr );//velo
  sol_b->GetLocalArray( array_b, node_ptr );//disp

  for( int ee=0; ee<nElem; ++ee )
  {
    if( 1 ) //eptr_array[ee]->is_sizeNonzero()
    {
      nlocbas_ee = lien_ptr->get_nLocBas_loc(ee);
      loc_dof= dof * nlocbas_ee;
      IEN_e = new int [nlocbas_ee];
      local_a = new double [dof * nlocbas_ee];
      local_b = new double [dof * nlocbas_ee];
      ectrl_x = new double [nlocbas_ee];
      ectrl_y = new double [nlocbas_ee];
      ectrl_z = new double [nlocbas_ee];
      row_index = new PetscInt [dof * nlocbas_ee];
      col_index = new PetscInt [dof * nlocbas_ee];

      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a, nlocbas_ee); //velo
      GetLocal(array_b, IEN_e, local_b, nlocbas_ee); //disp
      
      fnode_ptr->get_ctrlPts_xyz(nlocbas_ee, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      (lassem_array.at(ee))->
	Assem_Residual(t_n, dt, local_a, local_b,
		       eptr_array[ee], //local_iion, local_dphi,
		       ectrl_x, ectrl_y, ectrl_z, quad_array.at(ee));

      for(int i=0; i<nlocbas_ee; ++i)
      {
        loc_index = IEN_e[i];

        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);
          row_index[dof * i + m] = dof * lrow_index + m;
        }
      }

      //set residual values from element to global
      VecSetValues(G, loc_dof, row_index, lassem_array.at(ee)->Residual, ADD_VALUES);
    }
  }

  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
  //VecAssemblyBegin(sol_d->solution);
  //VecAssemblyEnd(sol_d->solution);

  std::cout<< "residual at the end of PGassem residual calc" << std::endl;
  Print_G();
  std::cout << "end assem residual." << std::endl ;  
}


void PGAssem_EP::Assem_mass_residual(
    const PDNSolution * const &sol_a,
    const PDNTimeStep * const &time_info,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const IQuadPts * const &quad, //changed weight with quad
    std::vector<FEAElement*> &eptr_array,
    const ALocal_NodalBC * const &bc_part )
{
  const int nElem = alelem_ptr->get_nlocalele();
  const int loc_dof = dof * nLocBas;
  int loc_index, lrow_index;

  sol_a->GetLocalArray( array_a, node_ptr);
     
  for(int ee=0; ee<nElem; ++ee)
  {
    if( 1 ) //eptr_array[ee]->is_sizeNonzero()
    {

      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a);
  
      fnode_ptr->get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      lassem_ptr->Assem_Mass_Residual(local_a,// local_b, local_c,
  				      eptr_array[ee], ectrl_x, ectrl_y,
  				      ectrl_z, quad); 
      
      for(int i=0; i<nLocBas; ++i)
      {
        loc_index = IEN_e[i];
  
        for(int m=0; m<dof; ++m)
        {

          lrow_index = bc_part->get_LID(m, loc_index);
  
          row_index[dof * i + m] = dof * lrow_index + m;
          col_index[dof * i + m] = dof * lrow_index + m;
        }
      }
      MatSetValues(K, loc_dof, row_index, loc_dof, col_index,
          lassem_ptr->Tangent, ADD_VALUES);
  
      VecSetValues(G, loc_dof, row_index, lassem_ptr->Residual, ADD_VALUES);
    }
  }
  
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
  
  for(int fie = 0; fie<dof; ++fie)
    EssBC_KG( bc_part, fie);
  
  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);

}


void PGAssem_EP::Assem_mass_residual(const PDNSolution * const &sol_a,
				     const PDNTimeStep * const &time_info,
				     const ALocal_Elem * const &alelem_ptr,
				     std::vector<IPLocAssem *> &lassem_array, 
				     const ALocal_IEN_Mixed * const &lien_ptr,
				     const APart_Node * const &node_ptr,
				     const FEANode * const &fnode_ptr,
				     const std::vector< IQuadPts * > &quad_array,
				     std::vector<FEAElement*> &eptr_array,
				     const ALocal_NodalBC * const &bc_part )
{

  const int nElem = alelem_ptr->get_nlocalele();
  int loc_dof;
  int loc_index, lrow_index;
  int nlocbas_ee;

  sol_a->GetLocalArray( array_a, node_ptr);
     
  for(int ee=0; ee<nElem; ++ee)
  {
    nlocbas_ee = lien_ptr->get_nLocBas_loc(ee);
    loc_dof= dof * nlocbas_ee;
    IEN_e = new int [nlocbas_ee];
    local_a = new double [dof * nlocbas_ee];
    ectrl_x = new double [nlocbas_ee];
    ectrl_y = new double [nlocbas_ee];
    ectrl_z = new double [nlocbas_ee];
    row_index = new PetscInt [dof * nlocbas_ee];
    col_index = new PetscInt [dof * nlocbas_ee];
    
    if( 1 ) //eptr_array[ee]->is_sizeNonzero()
    {
      
      lien_ptr->get_LIEN_e(ee, IEN_e);
      GetLocal(array_a, IEN_e, local_a, nlocbas_ee);

      std::cout << "local_a : " <<"\n";
      for (int i = 0; i < nlocbas_ee; ++i)
	std::cout << local_a[i] << "\t";
      std::cout << std::endl;

      fnode_ptr->get_ctrlPts_xyz(nlocbas_ee, IEN_e, ectrl_x, ectrl_y, ectrl_z);

      std::cout << "element " << ee<< " coords : " <<"\n"
		<< "x \t y \t z \n";
      for (int i=0; i<nlocbas_ee; ++i){
	std::cout << ectrl_x[i] << "\t" << ectrl_y[i] <<"\t"<< ectrl_z[i] << "\n";
      }

      (lassem_array.at(ee))->Assem_Mass_Residual(local_a,eptr_array.at(ee), ectrl_x,
						 ectrl_y, ectrl_z, quad_array.at(ee)); 
      
      for(int i=0; i<nlocbas_ee; ++i)
      {
        loc_index = IEN_e[i];
  
        for(int m=0; m<dof; ++m)
        {

          lrow_index = bc_part->get_LID(m, loc_index);
  
          row_index[dof * i + m] = dof * lrow_index + m;
          col_index[dof * i + m] = dof * lrow_index + m;
        }
      }
      MatSetValues(K, loc_dof, row_index, loc_dof, col_index,
      		   (lassem_array.at(ee))->Tangent, ADD_VALUES);
      VecSetValues(G, loc_dof, row_index, (lassem_array.at(ee))->Residual, ADD_VALUES);

      std::cout <<"local tangent:" << std::endl;
      PetscScalarView(nlocbas_ee*nlocbas_ee, (lassem_array.at(ee))->Tangent,PETSC_VIEWER_STDOUT_WORLD);
      std::cout <<"local residual:" << std::endl;
      PetscScalarView(nlocbas_ee, (lassem_array.at(ee))->Residual,PETSC_VIEWER_STDOUT_WORLD);
    }

    delete IEN_e; IEN_e = nullptr;
    delete local_a; local_a = nullptr;
    delete ectrl_x; ectrl_x = nullptr;
    delete ectrl_y; ectrl_y = nullptr;
    delete ectrl_z; ectrl_z = nullptr;
    delete row_index; row_index = nullptr;
    delete col_index; col_index = nullptr;
  }
  
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
  
  for(int fie = 0; fie<dof; ++fie)
    EssBC_KG( bc_part, fie);
  
  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(G);
  VecAssemblyEnd(G);
}


void PGAssem_EP::Update_nodal_velo(const PDNSolution * const &sol_a, //disp
				   const PDNSolution * const &sol_b, //pre_hist
				   const double &t_n,
				   const double &dt,
				   const IonicModel * const &ionicmodel_ptr,
				   const ALocal_Elem * const &alelem_ptr,
				   const ALocal_IEN * const &lien_ptr,
				   const APart_Node * const &node_ptr,
				   const FEANode * const &fnode_ptr,
				   const ALocal_NodalBC * const &bc_part,
				   PDNSolution * const &sol_c, //new hist
				   PDNSolution * const &sol_d //new hist
				   )
{
  //  std::cout << "curr time: " << t_n << std::endl;
  int nElem = alelem_ptr->get_nlocalele();
  int loc_dof = dof * nLocBas;
  int loc_index, lrow_index; // lcol_index;

//  //int node_num = sol_a->get_nlocal();
  int node_num = node_ptr->get_nlocghonode();
//    
  sol_a->GetLocalArray( array_a, node_ptr );//V_in
  sol_b->GetLocalArray( array_b, node_ptr );//hist_old
  sol_c->GetLocalArray( array_c, node_ptr );//V_new
  sol_d->GetLocalArray( array_d, node_ptr );//hist_new
//
//  int node_locgho =node_ptr->get_nlocghonode();
//  int dof=  lassem_ptr->get_dof();
//  double *array_iion = new double [ dof * node_locgho ];
//  double *array_dphi = new double [ dof * node_locgho ];
//  
  double r_new, r_old, V_new, V_in;
  double ctrl_x, ctrl_y, ctrl_z;
  double t_n_half {t_n + dt/2.0};
  double t_n1 {t_n + dt};
  int num{1};
  std::vector<double> Istim;
  Istim.resize(3);
  
  for (int count{ 0 }; count < node_num; ++count)
    {
      V_in     = array_a[count];      
      r_old    = array_b[count];

      //global_idx = node_ptr -> get_local_to_global(count);
      fnode_ptr->get_ctrlPts_xyz(num, &count,
				 &ctrl_x, &ctrl_y, &ctrl_z);
      ionicmodel_ptr-> get_Istim (Istim.at(0), t_n,
				  ctrl_x, ctrl_y, ctrl_z);
      ionicmodel_ptr-> get_Istim (Istim.at(1), t_n_half,
				  ctrl_x, ctrl_y, ctrl_z);
      ionicmodel_ptr-> get_Istim (Istim.at(2), t_n1,
				  ctrl_x, ctrl_y, ctrl_z);

      //ionicmodel_ptr-> Forward_Euler(r_old, dt, V_in, Istim, r_new, V_new);
      ionicmodel_ptr-> Runge_Kutta_4(r_old, dt, V_in, Istim, r_new, V_new);

      array_c   [count] = V_new;
      array_d   [count] = r_new;
    }

  for( int ee=0; ee<nElem; ++ee )
  {
    if( 1 ) //eptr_array[ee]->is_sizeNonzero()
    {
      lien_ptr->get_LIEN_e(ee, IEN_e);
//      GetLocal(array_a, IEN_e, local_a); //V_in    
//      GetLocal(array_b, IEN_e, local_b); //hist_old
      GetLocal(array_c, IEN_e, local_c); //V_new   
      GetLocal(array_d, IEN_e, local_d); //hist_new

      for(int i=0; i<nLocBas; ++i)
      {
        loc_index = IEN_e[i];

        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);
  
          row_index[dof * i + m] = dof * lrow_index + m;
	}
      }
      
      VecSetValues(sol_c->solution, loc_dof, row_index, local_c, INSERT_VALUES);
      VecSetValues(sol_d->solution, loc_dof, row_index, local_d, INSERT_VALUES);
    }
  }

  VecAssemblyBegin(sol_c->solution);
  VecAssemblyEnd(sol_c->solution);
  VecAssemblyBegin(sol_d->solution);
  VecAssemblyEnd(sol_d->solution);
  
}


//this version is for the nonuniform element types
void PGAssem_EP::Update_nodal_velo(const PDNSolution * const &sol_a, //disp
				   const PDNSolution * const &sol_b, //pre_hist
				   const double &t_n,
				   const double &dt,
				   const std::vector< IonicModel * > &ionicmodel_array,
				   const ALocal_Elem * const &alelem_ptr,
				   const ALocal_IEN_Mixed * const &lien_ptr,
				   const APart_Node * const &node_ptr,
				   const FEANode * const &fnode_ptr,
				   const ALocal_NodalBC * const &bc_part,
				   PDNSolution * const &sol_c, //new hist
				   PDNSolution * const &sol_d //new hist
				   )
{
  //member variables  that are not assigned before and should be assigned in
  // element for loops :
  // nLocBas , row_index, col_index, local_a b c d, IEN_e, ectrl_xyz
  //std::cout << "new nodal velo update: " << std::endl;
  int nElem = alelem_ptr->get_nlocalele();
  int loc_dof, nlocbas_ee;
  int loc_index, lrow_index; // lcol_index;

  int node_num = node_ptr->get_nlocghonode();

  sol_a->GetLocalArray( array_a, node_ptr );//V_in
  sol_b->GetLocalArray( array_b, node_ptr );//hist_old
  sol_c->GetLocalArray( array_c, node_ptr );//V_new
  sol_d->GetLocalArray( array_d, node_ptr );//hist_new

  double r_new, r_old, V_new, V_in;
  double ctrl_x, ctrl_y, ctrl_z;
  double t_n_half {t_n + dt/2.0};
  double t_n1 {t_n + dt};
  int num{1};
  std::vector<int> ee;
  std::vector<double> Istim;
  Istim.resize(3);

  for (int count{ 0 }; count < node_num; ++count)
    {
      lien_ptr->get_node_to_elem(count, ee);
      //std::cout << "node : " << count << "is in elements \n";
      //VEC_T::print(ee);

      V_in     = array_a[count];      
      r_old    = array_b[count];

      fnode_ptr->get_ctrlPts_xyz(num, &count,
				 &ctrl_x, &ctrl_y, &ctrl_z);
      ionicmodel_array[ee[0]]-> get_Istim (Istim.at(0), t_n,
				  ctrl_x, ctrl_y, ctrl_z);
      ionicmodel_array[ee[0]]-> get_Istim (Istim.at(1), t_n_half,
				  ctrl_x, ctrl_y, ctrl_z);
      ionicmodel_array[ee[0]]-> get_Istim (Istim.at(2), t_n1,
				  ctrl_x, ctrl_y, ctrl_z);

      //ionicmodel_ptr-> Forward_Euler(r_old, dt, V_in, Istim, r_new, V_new);
      ionicmodel_array[ee[0]]-> Runge_Kutta_4(r_old, dt, V_in, Istim, r_new, V_new);

      array_c   [count] = V_new;
      array_d   [count] = r_new;
    }

  for( int ee=0; ee<nElem; ++ee )
  {
    if( 1 ) //eptr_array[ee]->is_sizeNonzero()
    {
      nlocbas_ee = lien_ptr->get_nLocBas_loc(ee);
      loc_dof= dof * nlocbas_ee;
      
      row_index = new PetscInt [dof * nlocbas_ee];
      local_c =     new double [dof * nlocbas_ee];  
      local_d =     new double [dof * nlocbas_ee];
      IEN_e =                new int [nlocbas_ee];

      //  col_index = new PetscInt [dof * nLocBas];
      //  local_a = new double [dof * nLocBas];
      //  local_b = new double [dof * nLocBas];

      lien_ptr->get_LIEN_e(ee, IEN_e);

      GetLocal(array_c, IEN_e, local_c, nlocbas_ee); //V_new   
      GetLocal(array_d, IEN_e, local_d, nlocbas_ee); //hist_new

      for(int i=0; i<nlocbas_ee; ++i)
      {
        loc_index = IEN_e[i];

        for(int m=0; m<dof; ++m)
        {
          lrow_index = bc_part->get_LID(m, loc_index);
  
          row_index[dof * i + m] = dof * lrow_index + m;
	}
      }
      
      VecSetValues(sol_c->solution, loc_dof, row_index, local_c, INSERT_VALUES);
      VecSetValues(sol_d->solution, loc_dof, row_index, local_d, INSERT_VALUES);

      delete row_index; row_index=nullptr;
      delete local_c  ; local_c  =nullptr;
      delete local_d  ; local_d  =nullptr;
      delete IEN_e    ; IEN_e    =nullptr;
    }
  }

  VecAssemblyBegin(sol_c->solution);
  VecAssemblyEnd(sol_c->solution);
  VecAssemblyBegin(sol_d->solution);
  VecAssemblyEnd(sol_d->solution);
  
}
				     


//EOF
