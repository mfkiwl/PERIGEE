#include "ElemBC_3D_tet_wall.hpp"

ElemBC_3D_tet_wall::ElemBC_3D_tet_wall(
    const std::vector<std::string> &walls_combined,
    const std::string &centerlines_combined,
    const double &thickness2radius_combined,
    const double &fluid_density,
    const int &elemtype )
: ElemBC_3D_tet( walls_combined, elemtype ),
  fluid_density( fluid_density )
{
  // Check inputs
  SYS_T::print_fatal_if( walls_combined.size() != 1,
      "Error: walls_combined does not have a length of 1.\n");

  SYS_T::print_fatal_if( elemtype != 501,
      "Error: unsupported element type.\n");

  // num_ebc = 1 for wall elem bc
  const int ebc_id = 0;

  radius.resize(    num_node[ebc_id] );
  thickness.resize( num_node[ebc_id] );
  youngsmod.resize( num_node[ebc_id] );

  vtkXMLPolyDataReader * reader = vtkXMLPolyDataReader::New();
  reader -> SetFileName( centerlines_combined.c_str() );
  reader -> Update();
  vtkPolyData * centerlineData = reader -> GetOutput();

  vtkPointLocator * locator = vtkPointLocator::New();
  locator -> Initialize();
  locator -> SetDataSet( centerlineData );
  locator -> BuildLocator();

  // const double rho_alpha2 = fluid_density * youngsmod_alpha * youngsmod_alpha;
  // const double beta_exp   = 2.0 * youngsmod_beta - 1.0;

  for(int ii=0; ii<num_node[ebc_id]; ++ii)
  {
    const double pt[3] = {pt_xyz[ebc_id][3*ii], pt_xyz[ebc_id][3*ii+1], pt_xyz[ebc_id][3*ii+2]};

    const int closest_id = locator -> FindClosestPoint(&pt[0]);

    const double * cl_pt = centerlineData -> GetPoints() -> GetPoint(closest_id);

    radius[ii] = MATH_T::norm2(cl_pt[0] - pt[0], cl_pt[1] - pt[1], cl_pt[2] - pt[2]);
  
    thickness[ii] = radius[ii] * thickness2radius_combined; 

    youngsmod[ii] = 0.0;  // to be changed!

    // youngsmod[ii] = rho_alpha2 / ( thickness[ii] * pow( 2.0*radius[ii], beta_exp ) );
  }
 
  // clean memory
  locator -> Delete();
  reader -> Delete();

  // Write out vtp's with wall properties
  write_vtk(ebc_id, "varwallprop");
}


ElemBC_3D_tet_wall::ElemBC_3D_tet_wall(
    const std::vector<std::string> &walls_combined,
    const std::string &centerlines_combined,
    const double &thickness2radius_combined,
    const std::vector<std::string> &wallsList,
    const std::vector<std::string> &centerlinesList,
    const std::vector<double> &thickness2radiusList,
    const double &fluid_density,
    const int &elemtype )
: ElemBC_3D_tet_wall( walls_combined, centerlines_combined, thickness2radius_combined,
                      fluid_density, elemtype)
{
  // Check inputs
  SYS_T::print_fatal_if( centerlinesList.size() != wallsList.size(),
      "Error: centerlinesList length does not match that of wallsList.\n");

  SYS_T::print_fatal_if( thickness2radiusList.size() != wallsList.size(),
      "Error: thickness2radiusList length does not match that of wallsList.\n");
}


ElemBC_3D_tet_wall::~ElemBC_3D_tet_wall()
{
  VEC_T::clean( radius    );
  VEC_T::clean( thickness );
  VEC_T::clean( youngsmod );
}


void ElemBC_3D_tet_wall::print_info() const
{
  ElemBC_3D_tet::print_info();

  VEC_T::print( radius,    "wall_radius.txt",    '\n');
  VEC_T::print( thickness, "wall_thickness.txt", '\n');
  VEC_T::print( youngsmod, "wall_youngsmod.txt", '\n');
}


void ElemBC_3D_tet_wall::write_vtk( const int &ebc_id, const std::string &filename ) const
{
  vtkPolyData * grid_w = vtkPolyData::New();
  
  TET_T::gen_triangle_grid( grid_w, 
      num_node[ebc_id], num_cell[ebc_id],
      pt_xyz[ebc_id], tri_ien[ebc_id] );    

  // Add nodal indices
  TET_T::add_int_PointData( grid_w, global_node[ebc_id], "GlobalNodeID" );

  // Add cell indices
  TET_T::add_int_CellData( grid_w, global_cell[ebc_id], "GlobalElementID");

  // Add thickness
  TET_T::add_double_PointData( grid_w, thickness, "Thickness" );

  // Add Young's modulus
  TET_T::add_double_PointData( grid_w, youngsmod, "YoungsModulus" );

  // write vtp
  TET_T::write_vtkPointSet(filename, grid_w);

  // Clean memory
  grid_w->Delete();
}

// EOF
