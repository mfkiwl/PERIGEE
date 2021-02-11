// ==================================================================
// preprocess_tet_h5.cpp
// ------------------------------------------------------------------
// This preprocess code is used to handle 3D geometry discretized
// by arbitrary order tetrahedral elements, which is stored in a h5
// file.
//
// Date created: Dec. 6 2017
// Author: Ju Liu
// ==================================================================
#include "IEN_Gmsh.hpp"
#include "Mesh_FEM.hpp"
#include "Mesh_Tet4.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Map_Node_Index.hpp"
#include "NodalBC_3D.hpp"
#include "ElemBC_3D.hpp"
#include "Part_FEM.hpp"
#include "NBC_Partition_3D.hpp"
#include "EBC_Partition_FEM.hpp"

#include "NodalBC_3D_vtp.hpp"
#include "NodalBC_3D_vtu.hpp"
#include "IEN_Tetra_P1.hpp"
#include "Part_Tet.hpp"
#include "ElemBC_3D_tet4.hpp"
#include "EBC_Partition_vtp.hpp"

using std::cout;
using std::endl;

int main( int argc, char * argv[] )
{
  int sysret = system("rm -rf part_p*.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf preprocessor_cmd.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf NumLocalNode.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");

  const int probDim = 3;
  const int dofNum = 7;
  const int dofMat = 4;
  const int elemType = 531; // first order simplicial element

  const std::string part_file("part");

  std::string geo_file("/home/oguz/PERIGEE/examples/tet4_beam_bending/ventricle_inflation/mesh/coarse/mesh-complete/mesh-complete.mesh.vtu");
  std::string geo_base("/home/oguz/PERIGEE/examples/tet4_beam_bending/ventricle_inflation/mesh/coarse/mesh-complete/mesh-surfaces/base.vtp");
  std::string geo_endo("/home/oguz/PERIGEE/examples/tet4_beam_bending/ventricle_inflation/mesh/coarse/mesh-complete/mesh-surfaces/endo.vtp");
  std::string geo_epi ("/home/oguz/PERIGEE/examples/tet4_beam_bending/ventricle_inflation/mesh/coarse/mesh-complete/mesh-surfaces/epi.vtp");
  int cpu_size = 1;
  int in_ncommon = 2;
  bool isDualGraph = true;

  PetscMPIInt rank, size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  SYS_T::print_fatal_if(size!= 1,"ERROR: preprocessor is a serial program!\n");

  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionString("-geo_file", geo_file);
  SYS_T::GetOptionString("-geo_base", geo_base);
  SYS_T::GetOptionString("-geo_endo", geo_endo);
  SYS_T::GetOptionString("-geo_epi" , geo_epi );

  std::cout<<"==== Command Line Arguments ===="<<std::endl;
  std::cout<<" -geo_file: "<<geo_file<<std::endl;
  std::cout<<" -geo_base: "<<geo_base<<std::endl;
  std::cout<<" -geo_endo: "<<geo_endo<<std::endl;
  std::cout<<" -geo_epi: " <<geo_epi <<std::endl;
  std::cout<<" -cpu_size: "<<cpu_size<<std::endl;
  std::cout<<" -in_ncommon: "<<in_ncommon<<std::endl;
  std::cout<<" -isDualGraph: true \n";

  // ----------------------------------------------------------------
  // ------------- Write the input argument into a HDF5 file
  SYS_T::file_check(geo_file); std::cout<<geo_file<<" found. \n";
  SYS_T::file_check(geo_base); std::cout<<geo_base<<" found. \n";
  SYS_T::file_check(geo_endo); std::cout<<geo_endo<<" found. \n";
  SYS_T::file_check(geo_epi ); std::cout<<geo_epi <<" found. \n";


  hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5",
      H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

  cmdh5w->write_intScalar("cpu_size", cpu_size);
  cmdh5w->write_intScalar("in_ncommon", in_ncommon);
  cmdh5w->write_intScalar("dofNum", dofNum);
  cmdh5w->write_intScalar("dofMat", dofMat);
  cmdh5w->write_intScalar("elemType", elemType);
  cmdh5w->write_string("geo_file", geo_file);
  cmdh5w->write_string("geo_base", geo_base);
  cmdh5w->write_string("geo_endo", geo_endo);
  cmdh5w->write_string("geo_epi" , geo_epi );  
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; H5Fclose(cmd_file_id);
  // ----- Finish writing

  // Read the geometry file for the whole FSI domain
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<int> phy_tag;
  std::vector<double> ctrlPts;

  TET_T::read_vtu_grid(geo_file.c_str(), nFunc, nElem, ctrlPts, vecIEN);

  // Generate IEN
  IIEN * IEN = new IEN_Tetra_P1(nElem, vecIEN);
  VEC_T::clean( vecIEN );

  IMesh * mesh = new Mesh_Tet4(nFunc, nElem);
  mesh -> print_mesh_info();

  // Call METIS for partitioning
  IGlobal_Part * global_part;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
					 isDualGraph, mesh, IEN, "epart", "npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( mesh, "epart", "npart" );
  else
    {
      std::cerr<<"ERROR: wrong cpu_size: "<<cpu_size<<std::endl;
      exit(EXIT_FAILURE);
    }

  // Generate nodal index mapping
  Map_Node_Index * mnindex = new Map_Node_Index(global_part, cpu_size, mesh->get_nFunc());
  mnindex->write_hdf5("node_mapping");
 
  // Nodal BC 
  std::vector<INodalBC *> NBC_list; NBC_list.clear();
  NBC_list.resize( dofMat );
 
  std::vector<std::string> dir_list;
  dir_list.push_back( geo_base );
  // dir_list.push_back( geo_endo );
  // dir_list.push_back( geo_epi );
 
  NBC_list[0] = new NodalBC_3D_vtp( nFunc );
  NBC_list[1] = new NodalBC_3D_vtp(dir_list, nFunc );
  NBC_list[2] = new NodalBC_3D_vtp(dir_list, nFunc );
  NBC_list[3] = new NodalBC_3D_vtp(dir_list, nFunc );


  // Element BC
  std::vector<std::string> ebc_vtp_list; ebc_vtp_list.clear();
  ebc_vtp_list.push_back(geo_endo);
 
  ElemBC * ebc = new ElemBC_3D_tet4( ebc_vtp_list );

  // ----------------------------------------------------------------
  // Partition the mesh & BC
  const bool isPrintPartInfo = true;
  const int proc_size = cpu_size;
  std::vector<int> list_nlocalnode, list_nghostnode, list_ntotalnode, list_nbadnode;
  std::vector<double> list_ratio_g2l;
  int sum_nghostnode = 0;
  SYS_T::Timer * mytimer = new SYS_T::Timer();

  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
    {
      mytimer->Reset(); mytimer->Start();
      IPart * part = new Part_Tet( mesh, global_part, mnindex, IEN,
				   ctrlPts, proc_rank, proc_size, dofNum, dofMat,
				   elemType, isPrintPartInfo );
      mytimer->Stop();
      cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

      part -> write( part_file.c_str() );

      part -> print_part_loadbalance_edgecut();

      list_nlocalnode.push_back(part->get_nlocalnode());
      list_nghostnode.push_back(part->get_nghostnode());
      list_ntotalnode.push_back(part->get_ntotalnode());
      list_nbadnode.push_back(part->get_nbadnode());
      list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());
      sum_nghostnode += part->get_nghostnode();
  
      INBC_Partition * nbcpart = new NBC_Partition_3D(part, mnindex, NBC_list);
      nbcpart -> write_hdf5(part_file.c_str());
  
      //IEBC_Partition * ebcpart = new EBC_Partition_FEM(part, mnindex, ebc);
      IEBC_Partition * ebcpart = new EBC_Partition_vtp(part, mnindex, ebc);
      ebcpart -> write_hdf5(part_file.c_str());

      delete part;
      delete nbcpart;
      delete ebcpart;
    }
  // ----------------------------------------------------------------
  //VEC_T::write_int_h5("NumLocalNode","nln", list_nlocalnode);

  cout<<"\n===> Partition Quality: "<<endl;
  cout<<"The largest ghost / local node ratio is: ";
  cout<<*std::max_element(&list_ratio_g2l[0], &list_ratio_g2l[cpu_size-1])<<endl;
  cout<<"The smallest ghost / local node ratio is: ";
  cout<<*std::min_element(&list_ratio_g2l[0], &list_ratio_g2l[cpu_size-1])<<endl;
  cout<<"The summation of the number of ghost nodes is: "<<sum_nghostnode<<endl;
  cout<<"The maximum badnode number is: ";
  cout<<*std::max_element(&list_nbadnode[0], &list_nbadnode[cpu_size-1])<<endl;
  const int maxpart_nlocalnode = *std::max_element(&list_nlocalnode[0],
						   &list_nlocalnode[cpu_size-1]);
  const int minpart_nlocalnode = *std::min_element(&list_nlocalnode[0],
						   &list_nlocalnode[cpu_size-1]);
  cout<<"The maximum and minimum local node numbers are ";
  cout<<maxpart_nlocalnode<<"\t";
  cout<<minpart_nlocalnode<<endl;
  cout<<"The maximum / minimum of local node is: ";
  cout<<(double) maxpart_nlocalnode / (double) minpart_nlocalnode<<endl;

  // clean memory
  delete mytimer;
  delete ebc;
 
  std::vector<INodalBC *>::iterator it_nbc;
  for(it_nbc=NBC_list.begin(); it_nbc != NBC_list.end(); ++it_nbc) delete *it_nbc;
  delete IEN; delete mesh; delete global_part; delete mnindex;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
