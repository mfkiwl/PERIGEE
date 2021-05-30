#include "EBC_Partition_vtp_wall.hpp"

EBC_Partition_vtp_wall::EBC_Partition_vtp_wall( 
    const IPart * const &part,
    const Map_Node_Index * const &mnindex,
    const ElemBC * const &ebc )
: EBC_Partition_vtp(part, mnindex, ebc),
  fluid_density( ebc-> get_fluid_density() )
{
  const int ebc_id = 0;
  
  part_thickness.clear();
  part_youngsmod.clear();

  // wall has only one surface per the assumption in wall ebc  
  if(num_ebc == 0)
  {}
  else if(num_ebc == 1)
  {
    if( num_local_cell_node[ebc_id] > 0 )
    {
      // access wall properties of the whole surface
      std::vector<double> temp_th = ebc -> get_wall_thickness();
      std::vector<double> temp_E  = ebc -> get_wall_youngsmod();

      // save wall properties only belonging to nodes in the partition
      for( int ii=0; ii<num_local_cell_node[ebc_id]; ++ii )
      {
        part_thickness.push_back( temp_th[ local_cell_node[ebc_id][ii] ] );
        part_youngsmod.push_back( temp_E[  local_cell_node[ebc_id][ii] ] );
      }
    }
  }
  else
    SYS_T::print_fatal("Error: the num_ebc in EBC_Partition_vtp_wall should be 0 or 1. \n");

  // For wall surface, we keep a copy of the nodes belonging to this CPU's
  // subdomain
  local_node_pos.clear();
  for(int ii=0; ii<ebc->get_num_node(ebc_id); ++ii)
  {
    const int node_id = mnindex -> get_old2new( ebc->get_global_node(ebc_id, ii) );
    if( part -> isNodeInPart( node_id ) )
      local_node_pos.push_back( part -> get_nodeLocGhoIndex( node_id ) );
  }

  num_local_node = static_cast<int>( local_node_pos.size() );

}

EBC_Partition_vtp_wall::~EBC_Partition_vtp_wall()
{
  VEC_T::clean( part_thickness );
  VEC_T::clean( part_youngsmod );
  VEC_T::clean( local_node_pos );
}

void EBC_Partition_vtp_wall::write_hdf5( const char * FileName ) const
{
  // Call base class writer to write base class data at ebc_wall folder
  // certain wall info are in ebc_wall/ebc_0 sub-folder
  EBC_Partition_vtp::write_hdf5( FileName, "ebc_wall" );

  const std::string input_fName(FileName);
  const std::string fName = SYS_T::gen_partfile_name( input_fName, cpu_rank );

  // re-open the file
  hid_t file_id = H5Fopen(fName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  // open the folder at fName/ebc_wall again to append additional data 
  hid_t g_id = H5Gopen( file_id, "ebc_wall", H5P_DEFAULT );

  HDF5_Writer * h5w = new HDF5_Writer( file_id );

  h5w -> write_doubleScalar( g_id, "fluid_density", fluid_density );

  h5w -> write_intScalar( g_id, "num_local_node", num_local_node );

  // num_ebc = 1 for wall elem bc, and ebc_id is 0
  const int ebc_id = 0;
  if( num_local_cell[ebc_id] > 0 )
  {
    hid_t group_id = H5Gopen( g_id, "ebcid_0", H5P_DEFAULT );

    h5w->write_doubleVector( group_id, "thickness", part_thickness );
    h5w->write_doubleVector( group_id, "youngsmod", part_youngsmod );

    H5Gclose( group_id );
  }

  if( num_local_node > 0 )
  {
    hid_t group_id = H5Gopen( g_id, "ebcid_0", H5P_DEFAULT );

    h5w -> write_intVector( group_id, "local_node_pos", local_node_pos );

    H5Gclose( group_id );
  }

  delete h5w; H5Gclose( g_id ); H5Fclose( file_id );
}

void EBC_Partition_vtp_wall::write_hdf5( const char * FileName, const char * GroupName ) const
{
  // This function is NOT allowed. If the user uses the same groupname for the
  // wall and outlets, the wall data and the 0-th outlet data will be mixed.
  // We enforce the users to call the default write_hdf5, where the groupname
  // for the wall is fixed to be ebc_wall.
  SYS_T::print_fatal("Error: EBC_Partition_vtp_wall, write_hdf5 with groupname is not allowed.\n");
}

// EOF
