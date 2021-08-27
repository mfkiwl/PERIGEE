#include "ALocal_IEN_Mixed.hpp"
#include <numeric>

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

  stride.resize(nlocalele+1);
  std::partial_sum(nLocBas_loc.begin(), nLocBas_loc.end(), stride.begin()+1);
  
  //std::vector<int> row_LIEN;
  h5r -> read_intVector("LIEN", "LIEN", LIEN);
  if(int(LIEN.size()) != nElemXnLocBas_loc )
    SYS_T::print_fatal("Error: ALocal_IEN_Mixed::LIEN is in wrong format.\n");

  //const int nlocalnode = h5r -> read_intScalar("Local_Node", "nlocalnode");
  const int nlocghonode = h5r -> read_intScalar("Local_Node", "nlocghonode");
  const int nlocalnode = h5r -> read_intScalar("Local_Node", "nlocalnode");
  //std::vector<int> node_loc;
  //h5r -> read_intVector("Local_Node", "node_loc", node_loc);
  
  ////now find for each local node, which local elements it belongs
  //// (only the first element, that is, myocardium preferred)
  // below with a loop on elements

  // modify this in future: below code piece check node neighborhood only with local
  // elements. for a more robust implementation it should check global elements.
  // for example when a purkinje-myocardium junction node is a ghost node,  then
  // this check will see only purkinje (or myocardium) element.
  // alternatively,: node_to_elemet should be in nlocalnode instead of nlocghonode.
  node_to_element.resize(nlocghonode);

  std::vector<int> nodes;
  for (int ee=0; ee<nlocalele; ee++) {
    nodes.clear();
    nodes.assign(LIEN.begin()+stride.at(ee), LIEN.begin()+stride.at(ee+1));
        
    for (auto it=nodes.begin(); it!=nodes.end(); it++){
      node_to_element.at(*it).push_back(ee);
    }
  }
  node_to_element.erase(node_to_element.begin()+nlocalnode,
			node_to_element.end() );
  //std::cout << "node_to_element size:" << node_to_element.size() << std::endl;

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
  std::cout<<"node_to_element = "<<'\n';
  for( auto it = node_to_element.begin(); it != node_to_element.end(); ++it )
    VEC_T::print(*it);
  std::cout<<'\n';

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
