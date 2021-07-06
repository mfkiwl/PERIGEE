#include "VTK_Writer_EP_Mixed.hpp"

VTK_Writer_EP_Mixed::VTK_Writer_EP_Mixed(const IAGlobal_Mesh_Info * const &GMIptr,
					 const std::string &epart_file)
  :nElem(GMIptr->get_nElem())
{
  VIS_T::read_epart( epart_file, nElem, epart_map );
}


VTK_Writer_EP_Mixed::~VTK_Writer_EP_Mixed()
{
  VEC_T::clean( epart_map );
}



void VTK_Writer_EP_Mixed::writeOutput_compact(
    const IAGlobal_Mesh_Info * const &GMIptr,
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
    const bool &isXML )
{
  // Make sure to call element-type specific functions based on the elemType
  //int nqpts = quad->get_num_quadPts();
  int elemType;

  vtkUnstructuredGrid * gridData = vtkUnstructuredGrid::New();
  vtkPoints * points = vtkPoints::New();
  const int numDArrays = vdata_ptr->get_arrayCompSize();
  
  if(numDArrays != 1) SYS_T::print_fatal("Error: vdata size numDArrays != 1.\n");

  vtkDoubleArray ** dataVecs = new vtkDoubleArray * [numDArrays];
  for(int ii=0; ii<numDArrays; ++ii)
    {
      dataVecs[ii] = vtkDoubleArray::New();
      dataVecs[ii] -> SetNumberOfComponents( vdata_ptr->get_arraySizes(ii) );
      std::string temp_name = vdata_ptr->get_arrayNames(ii);
      dataVecs[ii] -> SetName( temp_name.c_str() );
      dataVecs[ii] -> SetNumberOfTuples( num_of_nodes );
    }

  vtkIntArray * anaprocId = vtkIntArray::New();
  anaprocId -> SetName("Analysis_Partition");
  anaprocId -> SetNumberOfComponents(1);

  vtkDoubleArray * vtk_fiber_ori = vtkDoubleArray::New();
  vtk_fiber_ori -> SetName("Fiber_Orientation");
  vtk_fiber_ori -> SetNumberOfComponents(3);
  vtk_fiber_ori -> SetNumberOfTuples(lelem_ptr->get_nlocalele());
  //

  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)  {
    
    nLocBas= GMIptr->get_nLocBas(lelem_ptr->get_elem_loc(ee));
    elemType= elemArray.at(ee)->get_Type();

    Interpolater intep {nLocBas, true};
			  
    lien_ptr -> get_LIEN_e(ee, IEN_e);
    fnode_ptr -> get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);

    //std::cout << "elemtype : " << elemType <<  "\n"
    //      	      << "nlocbas :  " << nLocBas  << std::endl;
    
    elemArray.at(ee)->buildBasis( quadArray.at(ee), ectrl_x, ectrl_y, ectrl_z );
    
    intep.interpolateVTKPts(IEN_e, ectrl_x, ectrl_y, ectrl_z, elemArray.at(ee), points );
    
    // Interpolate the transmembrane voltafe scalar
    std::vector<double> inputInfo; inputInfo.clear();
    int asize = vdata_ptr->get_arraySizes( 0 ); // 1
    for(int jj=0; jj<nLocBas; ++jj) {
      int pt_index = IEN_e[jj];
      for(int kk=0; kk<asize; ++kk) {
	inputInfo.push_back( pointArrays[0][pt_index * asize + kk ] );
      }
    }
    intep.interpolateVTKData( asize, IEN_e, inputInfo, elemArray.at(ee), dataVecs[0] );
    
    // Set mesh connectivity
    
    if (elemType == 512) { //line element
      VIS_T::setLineelem( IEN_e[0], IEN_e[1], gridData );
    }
    else if(elemType == 501) {//tet4 element
      VIS_T::setTetraelem( IEN_e[0], IEN_e[1], IEN_e[2], IEN_e[3], gridData );
    }
    else {
      SYS_T::print_fatal("Error VTK_Writer_EP_Mixed: this element type is not implemented.\n");
    }	
        
    // Analysis mesh partition 
    const int e_global = lelem_ptr->get_elem_loc(ee);
    std::vector<double> temp;
    lelem_ptr->get_fiber_ori_e(temp, ee);

    //vorticity -> SetNumberOfComponents( 9 );
    //vorticity -> SetName("Velocity Gradient");
    //vorticity -> SetNumberOfTuples( num_of_nodes );
    //vorticity -> InsertComponent(ii, 0, ux[ii] / num_adj_cell[ii] );
    
    anaprocId->InsertNextValue( epart_map[e_global] );
    vtk_fiber_ori->InsertComponent(ee, 0, temp.at(0) );
    vtk_fiber_ori->InsertComponent(ee, 1, temp.at(1) );
    vtk_fiber_ori->InsertComponent(ee, 2, temp.at(2) );
  }
  
    gridData -> SetPoints( points );
    points -> Delete();

  for(int ii=0; ii<numDArrays; ++ii)
    {
      gridData->GetPointData()->AddArray( dataVecs[ii] );
      dataVecs[ii]->Delete();
    }

  delete [] dataVecs;

  // Add cell data
  gridData->GetCellData()->AddArray(anaprocId);
  gridData->GetCellData()->AddArray(vtk_fiber_ori);
  anaprocId->Delete(); 
  vtk_fiber_ori->Delete();

  // If postprocess is parallel, record its partition
  if(size > 1)
    {
      int numCells = gridData->GetNumberOfCells();
      vtkIntArray * procId = vtkIntArray::New();
      procId->SetName("PostProcess_ID");
      procId->SetNumberOfComponents(1);
      for(int ii=0; ii<numCells; ++ii){
	procId->InsertComponent(ii, 0 , rank);
      }
      gridData->GetCellData()->AddArray(procId);
      procId->Delete();
    }

  // Write gridData
  VIS_T::writeVisFile( gridData, outputBName, outputName,
		       rank, size, sol_time, isXML );

  // Clean gridData
  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object ... \n");
  gridData->Delete();
  PetscPrintf(PETSC_COMM_WORLD, "-- End of vtk writer ... \n");
}




//void VTK_Writer_EP_Mixed::writeOutput(
//    const FEANode * const &fnode_ptr,
//    const ALocal_IEN * const &lien_ptr,
//    const ALocal_Elem * const &lelem_ptr,
//    const IVisDataPrep * const &vdata_ptr,
//    FEAElement * const &elemptr,
//    const IQuadPts * const &quad,
//    const double * const * const &pointArrays,
//    const int &rank, const int &size,
//    const double &sol_time,
//    const std::string &basename,
//    const std::string &outputBName,
//    const std::string &outputName,
//    const bool &isXML )
//{
//  vtkUnstructuredGrid * gridData = vtkUnstructuredGrid::New();
//  vtkPoints * points = vtkPoints::New();
//  const int numDArrays = vdata_ptr->get_arrayCompSize();
//  if(numDArrays != 2) SYS_T::print_fatal("Error: vdata size numDArrays != 2. \n");
//
//  vtkDoubleArray ** dataVecs = new vtkDoubleArray * [numDArrays];
//  for(int ii=0; ii<numDArrays; ++ii)
//  {
//    dataVecs[ii] = vtkDoubleArray::New();
//    dataVecs[ii] -> SetNumberOfComponents( vdata_ptr->get_arraySizes(ii) );
//    std::string temp_name = vdata_ptr->get_arrayNames(ii);
//    dataVecs[ii] -> SetName( temp_name.c_str() );
//  }
//
//  vtkIntArray * anaprocId = vtkIntArray::New();
//  anaprocId -> SetName("Analysis_Partition");
//  anaprocId -> SetNumberOfComponents(1);
//
//  int ptOffset = 0;
//  for(int ee=0; ee<lelem_ptr->get_nlocalele(); ++ee)
//  {
//    lien_ptr -> get_LIEN_e(ee, IEN_e);
//    fnode_ptr -> get_ctrlPts_xyz(nLocBas, IEN_e, ectrl_x, ectrl_y, ectrl_z);
//    elemptr->buildBasis( quad, ectrl_x, ectrl_y, ectrl_z );
//
//    intep.interpolateVTKPts(ptOffset, ectrl_x, ectrl_y, ectrl_z,
//        elemptr, points );
//
//    // Interpolate the pressure scalar
//    std::vector<double> inputInfo; inputInfo.clear();
//    int asize = vdata_ptr->get_arraySizes( 0 );
//    for(int jj=0; jj<nLocBas; ++jj)
//    {
//      int pt_index = IEN_e[jj];
//      for(int kk=0; kk<asize; ++kk)
//        inputInfo.push_back( pointArrays[0][pt_index * asize + kk ] );
//    }
//    intep.interpolateVTKData( asize, ptOffset, &inputInfo[0],
//        elemptr, dataVecs[0] );
//
//    // Interpolate velo 
//    asize = vdata_ptr -> get_arraySizes( 1 );
//    inputInfo.clear();
//    for(int jj=0; jj<nLocBas; ++jj)
//    {
//      int pt_index = IEN_e[jj];
//      for(int kk=0; kk<asize; ++kk)
//        inputInfo.push_back( pointArrays[1][pt_index * asize + kk] );
//    }
//    intep.interpolateVTKData( asize, ptOffset, &inputInfo[0],
//        elemptr, dataVecs[1] );
//
//    // Set mesh connectivity
//    VIS_T::setTetraelem( ptOffset, gridData );
//
//    // Analysis mesh partition 
//    int e_global = lelem_ptr->get_elem_loc(ee);
//    anaprocId->InsertNextValue( epart_map[e_global] );
//
//    // update offset
//    ptOffset += 4;
//  }
//
//  gridData -> SetPoints( points );
//  points -> Delete();
//
//  for(int ii=0; ii<numDArrays; ++ii)
//  {
//    gridData->GetPointData()->AddArray( dataVecs[ii] );
//    dataVecs[ii]->Delete();
//  }
//
//  gridData->GetCellData()->AddArray(anaprocId);
//
//  delete [] dataVecs;
//  anaprocId->Delete();
//
//  // If postprocess is parallel, record its partition
//  if(size > 1)
//  {
//    int numCells = gridData->GetNumberOfCells();
//    vtkIntArray * procId = vtkIntArray::New();
//    procId->SetName("PostProcess_ID");
//    procId->SetNumberOfComponents(1);
//    for(int ii=0; ii<numCells; ++ii)
//      procId->InsertComponent(ii, 0 , rank);
//    gridData->GetCellData()->AddArray(procId);
//    procId->Delete();
//  }
//
//  // Write gridData
//  VIS_T::writeVisFile( gridData, outputBName, outputName,
//      rank, size, sol_time, isXML );
//
//  // Clean gridData
//  PetscPrintf(PETSC_COMM_WORLD, "-- Clean gridData object ... \n");
//  gridData->Delete();
//}
//
//


// EOF
