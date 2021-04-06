#include "ALocal_IEN_Mixed.hpp"

ALocal_IEN_Mixed::ALocal_IEN_Mixed( const std::string &fileBaseName, const int &cpu_rank )
{
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );

  nlocalele = h5r -> read_intScalar("Local_Elem", "nlocalele");

  //nLocBas = h5r -> read_intScalar("Global_Mesh_Info", "nLocBas");
  h5r -> read_intVector("Local_Elem", "nLocBas_loc", nLocBas_loc);
  int nElemXnLocBas_loc;
  nElemXnLocBas_loc = std::accumulate(nLocBas_loc.begin(), nLocBas_loc.end(), 0);

  stride.resize(nlocalele);
  std::partial_sum(nLocBas_loc.begin(), nLocBas_loc.end(), stride.begin()+1);
  
  //std::vector<int> row_LIEN;
  h5r -> read_intVector("LIEN", "LIEN", LIEN);
  if(int(LIEN.size()) != nElemXnLocBas_loc )
    SYS_T::print_fatal("Error: ALocal_IEN_Mixed::LIEN is in wrong format.\n");

  delete h5r; H5Fclose( file_id );
}


ALocal_IEN_Mixed::~ALocal_IEN_Mixed()
{}


void ALocal_IEN_Mixed::print_info() const
{
  std::cout<<"ALocal_IEN_Mixed: \n";

  std::cout<<"nlocalele = "<<nlocalele<<'\n';
  std::cout<<"nLocBas_loc = "<<'\n';
  VEC_T::print(nLocBas_loc);
  std::cout<<"stride = "<<'\n';
  VEC_T::print(stride);

  std::cout<<"LIEN = "<<'\n';
  for(int ee = 0; ee<nlocalele; ++ee)
  {
    std::vector<int> temp;
    get_LIEN_e(ee, temp);
    std::cout<<ee<<" : \t";
    VEC_T::print(temp);
    std::cout<<'\n';
  }
}

// EOF
