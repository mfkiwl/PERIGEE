#include "ElemBC_3D_tet_wall.hpp"

ElemBC_3D_tet_wall::ElemBC_3D_tet_wall(
    const std::vector<std::string> &vtkfileList,
    const std::vector<double> &thickness_to_radius,
    const std::vector<double> &youngsmod_alpha,
    const std::vector<double> &youngsmod_beta,
    const std::string &centerlineFile,
    const int &fluid_density,
    const int &elemtype )
: ElemBC_3D_tet( vtkfileList, elemtype )
{
  // Check inputs
  SYS_T::print_fatal_if( thickness_to_radius.size() != vtkfileList.size(),
      "Error: thickness_to_radius length does not match that of the vtkfileList.\n");

  SYS_T::print_fatal_if( youngsmod_alpha.size() != vtkfileList.size(),
      "Error: youngsmod_alpha length does not match that of the vtkfileList.\n");

  SYS_T::print_fatal_if( youngsmod_beta.size() != vtkfileList.size(),
      "Error: youngsmod_beta length does not match that of the vtkfileList.\n");

  SYS_T::print_fatal_if( elemtype != 501,
      "Error: unsupported element type.\n");

  radius.resize(num_ebc);
  thickness.resize(num_ebc);
  youngsmod.resize(num_ebc);
  for(int ii=0; ii<num_ebc; ++ii) 
  {
    radius[ii].resize( num_node[ii] );
    thickness[ii].resize( num_node[ii] );
    youngsmod[ii].resize( num_node[ii] );
  }

  vtkXMLPolyDataReader * reader = vtkXMLPolyDataReader::New();
  reader -> SetFileName( centerlineFile.c_str() );
  reader -> Update();

  vtkPolyData * centerlineData = reader -> GetOutput();
  
  vtkPointLocator * locator = vtkPointLocator::New();
  locator -> Initialize();
  locator -> SetDataSet( centerlineData );
  locator -> BuildLocator();

  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
  {
    const double rho_alpha2 = fluid_density * youngsmod_alpha[ebc_id] * youngsmod_alpha[ebc_id];
    const double beta_exp   = 2.0 * youngsmod_beta[ebc_id] - 1.0; 

    for(int ii=0; ii<num_node[ebc_id]; ++ii)
    {
      const double pt[3] = {pt_xyz[ebc_id][3*ii], pt_xyz[ebc_id][3*ii+1], pt_xyz[ebc_id][3*ii+2]};

      const int closest_id = locator -> FindClosestPoint(&pt[0]);

      const double * cl_pt = centerlineData -> GetPoints() -> GetPoint(closest_id);

      radius[ebc_id][ii] = MATH_T::norm2(cl_pt[0] - pt[0], cl_pt[1] - pt[1], cl_pt[2] - pt[2]);
   
      thickness[ebc_id][ii] = radius[ebc_id][ii] * thickness_to_radius[ebc_id]; 

      youngsmod[ebc_id][ii] = rho_alpha2 / ( thickness[ebc_id][ii] * pow( 2.0*radius[ebc_id][ii], beta_exp ) );
    }
  }

  // Write out vtp's with wall properties
  for(int ebc_id = 0; ebc_id < num_ebc; ++ebc_id)
    write_wall_prop(ebc_id, "varwallprop_" + SYS_T::to_string(ebc_id));

  // clean memory
  locator -> Delete();
  reader -> Delete();
}


ElemBC_3D_tet_wall::~ElemBC_3D_tet_wall()
{
  for(int ii=0; ii<num_ebc; ++ii)
  {
    VEC_T::clean( radius[ii]    );
    VEC_T::clean( thickness[ii] );
    VEC_T::clean( youngsmod[ii] );
  }

  VEC_T::clean( radius    );
  VEC_T::clean( thickness );
  VEC_T::clean( youngsmod );
}


void ElemBC_3D_tet_wall::print_info() const
{
  ElemBC_3D_tet::print_info();

  for(int face=0; face<num_ebc; ++face)
  {
    VEC_T::print( radius[face],    "wall_id_" + SYS_T::to_string(face) + "_radius.txt",    '\n');
    VEC_T::print( thickness[face], "wall_id_" + SYS_T::to_string(face) + "_thickness.txt", '\n');
    VEC_T::print( youngsmod[face], "wall_id_" + SYS_T::to_string(face) + "_youngsmod.txt", '\n');
  }
}


void ElemBC_3D_tet_wall::write_wall_prop( const int &ebc_id, const std::string &filename ) const
{
    vtkPolyData * grid_w = vtkPolyData::New();
    TET_T::gen_triangle_grid( grid_w, 
        num_node[ebc_id], num_cell[ebc_id],
        pt_xyz[ebc_id], tri_ien[ebc_id], 
        global_node[ebc_id], global_cell[ebc_id] );    
    
    // Add thickness
    vtkDoubleArray * prop0 = vtkDoubleArray::New();
    prop0 -> SetName("Thickness");
    prop0 -> SetNumberOfComponents(1);
    for(int ii=0; ii<num_node[ebc_id]; ++ii)
    {
      prop0 -> InsertNextValue( thickness[ebc_id][ii] );
    }
    grid_w -> GetPointData() -> AddArray( prop0 );
    prop0->Delete();

    // Add Young's modulus
    vtkDoubleArray * prop1 = vtkDoubleArray::New();
    prop1 -> SetName("YoungsModulus");
    prop1 -> SetNumberOfComponents(1);
    for(int ii=0; ii<num_node[ebc_id]; ++ii)
      prop1 -> InsertNextValue( youngsmod[ebc_id][ii] );
    grid_w -> GetPointData() -> AddArray( prop1 );
    prop1->Delete();

    // write vtp
    TET_T::write_vtkXMLPolyData(filename, grid_w);

    // Clean memory
    grid_w->Delete();

}


// EOF
