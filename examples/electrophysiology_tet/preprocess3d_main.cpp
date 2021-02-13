// ==================================================================
// preprocess_3d_main.cpp for purkinje line model
// ------------------------------------------------------------------
// This preprocess code is used for handling the 3D geometry described
// by tetrahdral elements 
//
//?? Users should call the sv_fsi_converter to convert the node and
//?? element indices before calling this routine.
//
// The users are also responsible for providing the proper file names
// for this routine to handle.
//
// Author: Ju Liu
// ModifiedL: Oguz Ziya Tikenogullari
// Date: November 2019
// ==================================================================
#include "Math_Tools.hpp"
#include "Mesh_Line_3D.hpp"
#include "IEN_Line_P1.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_Line.hpp"
#include "NodalBC_Line_3D_vtp.hpp"
#include "NodalBC_3D_vtu.hpp"
#include "NodalBC_Line_3D_stimulus.hpp"
#include "ElemBC_3D_Line.hpp"
#include "NBC_Partition_3D.hpp"
#include "NBC_Partition_3D_inflow.hpp"
#include "EBC_Partition_vtp_outflow.hpp"

int main( int argc, char * argv[] )
{
  // Remove previously existing hdf5 files
  int sysret = system("rm -rf part_p*.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf preprocessor_cmd.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  sysret = system("rm -rf NumLocalNode.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");

  // Define basic settings
  const int dofNum = 1; // degree-of-freedom for the physical problem
  const int dofMat = 1; // degree-of-freedom in the matrix problem
  const int elemType = 512; //2-node line element in 3d. check if this
  //element type number coincides with another element type number.
  //because I gave this number to this element.


  // Input files
  // volume
  std::string geo_file("./purkinje.vtu");

//  // faces purkinje mesh  
  std::string sur_file_tip0("./tip0_curve.vtp");
  std::string sur_file_tip1("./tip1_curve.vtp");
  
//  // faces for HLHS mesh  
//  std::string sur_file_Base("./Base_Vol.vtp");
//  std::string sur_file_Epi ("./Epi_Vol.vtp");
//  std::string sur_file_RV  ("./RV_Vol.vtp");
//  std::string sur_file_LV  ("./LV_Vol.vtp");
//  
//  // faces for cube/beam mesh  
//  //std::string sur_file_top("./top_vol.vtp");
//  //std::string sur_file_bot("./bot_vol.vtp");
//  //std::string sur_file_lef("./lef_vol.vtp");
//  //std::string sur_file_rig("./rig_vol.vtp");
//  //std::string sur_file_fro("./fro_vol.vtp");
//  //std::string sur_file_bac("./bac_vol.vtp");
//
  const std::string part_file("part");

  int cpu_size = 1;
  int in_ncommon = 1;
  const bool isDualGraph = true;

  PetscMPIInt size, rank;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // Enforce serial run.
  if(size != 1)
  {
    cout<<"ERROR: Given processor number is greater than 1."; 
    cout<<"Preprocess code has to be serial!"<<endl;
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }

  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionString("-geo_file", geo_file);
  SYS_T::GetOptionString("-sur_file_tip0", sur_file_tip0);
  SYS_T::GetOptionString("-sur_file_tip1", sur_file_tip1 );


  std::cout<<"==== /Command Line Arguments ===="<<std::endl;
  std::cout<<" -geo_file: "<<geo_file<<std::endl;
  std::cout<<" -sur_file_tip0" <<sur_file_tip0<<std::endl;
  std::cout<<" -sur_file_tip1" <<sur_file_tip1<<std::endl;


  std::cout<<" -part_file: "<<part_file<<std::endl;
  std::cout<<" -cpu_size: "<<cpu_size<<std::endl;
  std::cout<<" -in_ncommon: "<<in_ncommon<<std::endl;
  std::cout<<" -isDualGraph: true \n";
  std::cout<<"----------------------------------\n";
  std::cout<<" dofNum: "<<dofNum<<std::endl;
  std::cout<<" dofMat: "<<dofMat<<std::endl;
  std::cout<<" elemType: "<<elemType<<std::endl;
  std::cout<<"====  Command Line Arguments/ ===="<<std::endl;

  // Check if the geometrical file exist on disk
  SYS_T::file_check(geo_file); std::cout<<geo_file<<" found. \n";
  SYS_T::file_check(sur_file_tip0); std::cout<<sur_file_tip0<<" found. \n";
  SYS_T::file_check(sur_file_tip1); std::cout<<sur_file_tip1<<" found. \n";


// ----- Write the input argument into a HDF5 file
  hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5",
				H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

  cmdh5w->write_intScalar("cpu_size", cpu_size);
  cmdh5w->write_intScalar("in_ncommon", in_ncommon);
  cmdh5w->write_intScalar("dofNum", dofNum);
  cmdh5w->write_intScalar("dofMat", dofMat);
  cmdh5w->write_intScalar("elemType", elemType);
  cmdh5w->write_string("geo_file", geo_file);
  cmdh5w->write_string("sur_file_tip0", sur_file_tip0);
  cmdh5w->write_string("sur_file_tip1", sur_file_tip1);

  
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; H5Fclose(cmd_file_id);
  // ----- Finish writing

  // Read the geometry file for the whole FSI domain
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<int> phy_tag;
  std::vector<double> ctrlPts;

  TET_T::read_purkinje_lines(geo_file.c_str(), nFunc, nElem, ctrlPts, vecIEN, phy_tag);

  //for(unsigned int ii=0; ii<phy_tag.size(); ++ii)
  //{
  //  if(phy_tag[ii] != 0 && phy_tag[ii] != 1) SYS_T::print_fatal("Error: FSI problem, the physical tag for element should be 0 (fluid domain) or 1 (solid domain).\n");
  //}

  // Generate IEN
  IIEN * IEN = new IEN_Line_P1(nElem, vecIEN);
  //IEN->print_IEN();

  if(elemType == 512)
  {
    SYS_T::print_fatal_if(vecIEN.size() / nElem != 2, "Error: the mesh connectivity array size does not match with the element type 512. \n");
  }
  else
  {
    SYS_T::print_fatal_if(1, "Error: this script doesn't support this element type. \n");
  }
  
  VEC_T::clean( vecIEN );

  //  // Generate the list of nodes for fluid and solid
  //  std::vector<int> node_f, node_s; node_f.clear(); node_s.clear();
  //  for(int ee=0; ee<nElem; ++ee)
  //  {
  //    if( phy_tag[ee] == 0 )
  //    {
  //      for(int ii=0; ii<4; ++ii) node_f.push_back( IEN->get_IEN(ee, ii) );
  //    }
  //    else
  //    {
  //      for(int ii=0; ii<4; ++ii) node_s.push_back( IEN->get_IEN(ee, ii) );
  //    }
  //  }
  //  
  //  VEC_T::sort_unique_resize( node_f );
  //  VEC_T::sort_unique_resize( node_s );
  //
  //  std::cout<<'\n'<<"Fluid domain number of nodes: "<<node_f.size()<<'\n';
  //  std::cout<<"Solid domain number of nodes: "<<node_s.size()<<'\n';
  //
  
  // Check the mesh: I don't check line elements like tet elements
  // because aspect ratios are not critical 
  //TET_T::tetmesh_check(ctrlPts, IEN, nElem, 3.5);


  // Generate the mesh
  IMesh * mesh = new Mesh_Line_3D(nFunc, nElem);
  mesh -> print_mesh_info();

  // Partition
  IGlobal_Part * global_part;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
					 isDualGraph, mesh, IEN,
					 "epart", "npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( mesh, "epart", "npart" );
  else
    {
      std::cerr<<"ERROR: wrong cpu_size: "<<cpu_size<<std::endl;
      exit(EXIT_FAILURE);
    }

  Map_Node_Index * mnindex =
    new Map_Node_Index(global_part, cpu_size, mesh->get_nFunc());
  mnindex->write_hdf5("node_mapping");

  // ----------------------------------------------------------------
  // Setup boundary condition
  std::vector<INodalBC *> NBC_list;
  NBC_list.clear();
  NBC_list.resize( dofMat );

  // Use dir list to set Dirichlet BCs on these vtp surfaces
  //std::vector<std::string> dir_list;
  //dir_list.push_back(sur_file_top);
  
  NBC_list[0] = new NodalBC_Line_3D_vtp( nFunc );
  //NBC_list[0]->print_info();
  //NBC_list[1] = new NodalBC_3D_vtp( dir_list, nFunc );
  //NBC_list[2] = new NodalBC_3D_vtp( dir_list, nFunc );
  //NBC_list[3] = new NodalBC_3D_vtp( dir_list, nFunc );

  //std::vector<double> inflow_outward_normal;
  //inflow_outward_normal.push_back(0.0);
  //inflow_outward_normal.push_back(0.0);
  //inflow_outward_normal.push_back(1.0);
  INodalBC * InFBC = new NodalBC_Line_3D_stimulus(  nFunc );

  std::vector<std::string> ebclist;
  ebclist.clear();
  ElemBC * ebc = new ElemBC_3D_Line( ebclist );
  //ebc->print_info();

  const bool isPrintPartInfo = true;
  const int proc_size = cpu_size;

  std::vector<int> list_nlocalnode, list_nghostnode, list_ntotalnode, list_nbadnode;
  std::vector<double> list_ratio_g2l;

  int sum_nghostnode = 0;

  SYS_T::Timer * mytimer = new SYS_T::Timer();
  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    mytimer->Reset();
    mytimer->Start();
    IPart * part = new Part_Line( mesh, global_part, mnindex, IEN,
        ctrlPts, proc_rank, proc_size, dofNum, dofMat, elemType,
        isPrintPartInfo );
    mytimer->Stop();
    //part->print_part_ele() ;
    //part->print_part_node();
    //part->print_part_ghost_node() ;
    //part->print_part_local_to_global() ;
    //part->print_part_LIEN() ;
    //part->print_part_loadbalance_edgecut() ;
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";

    part -> write( part_file.c_str() );

    part -> print_part_loadbalance_edgecut();

    // Partition Nodal BC
    INBC_Partition * nbcpart = new NBC_Partition_3D(part, mnindex, NBC_list);
    nbcpart -> write_hdf5(part_file.c_str());

    // Partition Inflow BC
    INBC_Partition * infpart
      = new NBC_Partition_3D_inflow(part, mnindex, InFBC);
    infpart->write_hdf5( part_file.c_str() );

    // Partition Elem BC
    IEBC_Partition * ebcpart = new EBC_Partition_vtp(part, mnindex, ebc);
    ebcpart -> write_hdf5(part_file.c_str());

    list_nlocalnode.push_back(part->get_nlocalnode());
    list_nghostnode.push_back(part->get_nghostnode());
    list_ntotalnode.push_back(part->get_ntotalnode());
    list_nbadnode.push_back(part->get_nbadnode());
    list_ratio_g2l.push_back((double)part->get_nghostnode()/(double) part->get_nlocalnode());

    sum_nghostnode += part->get_nghostnode();
    delete part; delete nbcpart; delete ebcpart; delete infpart;
  }

  VEC_T::write_int_h5("NumLocalNode","nln", list_nlocalnode);

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

  // Free memory
  delete ebc; delete InFBC;
  std::vector<INodalBC *>::iterator it_nbc;
  for(it_nbc=NBC_list.begin(); it_nbc != NBC_list.end(); ++it_nbc) delete *it_nbc;

  delete mnindex; delete global_part; delete mesh; delete IEN; 
  delete mytimer;
  PetscFinalize();
  return EXIT_SUCCESS;
}

// EOF
