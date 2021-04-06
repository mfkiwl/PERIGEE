#include "AGlobal_Mesh_Info_Mixed.hpp"

AGlobal_Mesh_Info_Mixed::AGlobal_Mesh_Info_Mixed( 
    const std::string &fileBaseName,
    const int &cpu_rank )
{
  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );

  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

  HDF5_Reader * h5r = new HDF5_Reader( file_id );
  
  std::vector<int> stu_degree;
  int num_row, num_col;
  h5r -> read_intMatrix("Global_Mesh_Info", "degree",
			stu_degree, num_row, num_col);

  xdegree.resize(num_row);
  ydegree.resize(num_row);
  zdegree.resize(num_row);
  
  for (int ii=0; ii < num_row ; ++ii){
	xdegree.at(ii)= stu_degree.at(ii*num_col + 0);
	ydegree.at(ii)= stu_degree.at(ii*num_col + 1);
	zdegree.at(ii)= stu_degree.at(ii*num_col + 2);
  }

  nElem    = h5r -> read_intScalar("Global_Mesh_Info", "nElem");
  nFunc    = h5r -> read_intScalar("Global_Mesh_Info", "nFunc");
  probDim  = h5r -> read_intScalar("Global_Mesh_Info", "probDim");
  h5r -> read_intVector("Global_Mesh_Info", "nLocBas", nLocBas);
  h5r -> read_intVector("Global_Mesh_Info", "elemType", elemType);

  if(num_col != 3)
    SYS_T::print_fatal("Error: s-t-u degree is not given for 3 dims. \n");
  if(num_row != nElem)
    SYS_T::print_fatal("Error: s-t-u degree is not given for all elements. \n");

  delete h5r;
  H5Fclose( file_id );
}


AGlobal_Mesh_Info_Mixed::AGlobal_Mesh_Info_Mixed( 
    const std::string &fileBaseName,
    const int &cpu_rank, const int &num_enrich_node )
{
//  std::string fName = SYS_T::gen_partfile_name( fileBaseName, cpu_rank );
//
//  hid_t file_id = H5Fopen( fName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
//
//  HDF5_Reader * h5r = new HDF5_Reader( file_id );
//  
//  std::vector<int> vdeg;
//
//  h5r -> read_intVector("Global_Mesh_Info", "degree", vdeg);
//  xdegree = vdeg[0];
//  ydegree = vdeg[1];
//  zdegree = vdeg[2];
//
//  nElem    = h5r -> read_intScalar("Global_Mesh_Info", "nElem");
//  nFunc    = h5r -> read_intScalar("Global_Mesh_Info", "nFunc");
//  nLocBas  = h5r -> read_intScalar("Global_Mesh_Info", "nLocBas");
//  probDim  = h5r -> read_intScalar("Global_Mesh_Info", "probDim");
//  elemType = h5r -> read_intScalar("Global_Mesh_Info", "elemType");
//
//  delete h5r; H5Fclose( file_id );
//  
//  if(num_enrich_node < 0)
  SYS_T::print_fatal("Error: AGlobal_Mesh_Info_Mixed::AGlobal_Mesh_Info_Mixed enrichment nodes is not implemented . \n");

  //nFunc = nFunc + nElem * num_enrich_node;
  //nLocBas = nLocBas + num_enrich_node;
  //elemType = elemType + 10;
}


AGlobal_Mesh_Info_Mixed::~AGlobal_Mesh_Info_Mixed()
{}


void AGlobal_Mesh_Info_Mixed::print_info() const
{
  std::cout<<"AGlobal_Mesh_Info_Mixed:"<<std::endl;
  std::cout<<"degree: \n" ;
  for (int i=0; i < nElem ; ++i){
    std::cout<< xdegree.at(i)<<'\t'<<ydegree.at(i)<<'\t'<<zdegree.at(i)<<'\n';
  }
  std::cout<<"nElem: "<<nElem<<'\n';
  std::cout<<"nFunc: "<<nFunc<<'\n';
  std::cout<<"nLocBas: "<<std::endl;
  VEC_T::print(nLocBas);
  std::cout<<"probDim: "<<probDim<<std::endl;
  std::cout<<"elemType: "<<std::endl;
  VEC_T::print(elemType);
  
}


// EOF
