#include "ALocal_Elem_Fiber.hpp"

ALocal_Elem_Fiber::ALocal_Elem_Fiber(const std::string &fileBaseName,
				     const int &cpu_rank)
  : ALocal_Elem {fileBaseName, cpu_rank}
{
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  const std::string gname("/Local_Elem");
  
  std::vector<double> fiber_ori_vec;
  int row_num, col_num;
  h5r->read_doubleMatrix( gname.c_str(), "fiber_ori_loc_vec",
			  fiber_ori_vec, row_num, col_num );

  std::vector<double>::iterator it = fiber_ori_vec.begin();
  for (int ii=0; ii < row_num ; ++ii){
    fiber_ori_loc.push_back( std::vector<double> (it+ ii*col_num,
						  it+((ii+1)*col_num) ) );
  }

  //read in the physical tag of each local element.
  h5r->read_intVector( gname.c_str(), "phy_tag_loc", phy_tag_loc);

  delete h5r;
  H5Fclose( file_id );
}


ALocal_Elem_Fiber::~ALocal_Elem_Fiber()
{
  
}


//void ALocal_Elem_Fiber::print_info() const
//{
//  std::cout<<"nlocalelem: "<<nlocalele<<std::endl;
//  std::cout<<"elem_loc: \n";
//  VEC_T::print(elem_loc);
//}

// EOF
