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
#include "Mesh_Tet4.hpp"
#include "Mesh_Mixed.hpp"
#include "IEN_Line_P1.hpp"
#include "IEN_Tetra_P1.hpp"
#include "IEN_Mixed.hpp"
#include "Global_Part_METIS_Mixed.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_Line.hpp"
#include "Part_Mixed_Mesh.hpp"
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
  const int elemType_myo = 501; //first order (4 node) tet element.
  const int elemType_LVpur = 512; //2-node line element in 3d. check if this
  const int elemType_RVpur = 512; //2-node line element in 3d. check if this
  const int phy_tag_myo    = 1;   // tag number to be assigned to myocardium  
  const int phy_tag_LVpur  = 2;   // tag number to be assigned to the LV
  const int phy_tag_RVpur  = 3;   // and RV purkinje cells
  // check if element type number 512 coincides with another element type number.
  //because I gave this number to this element.
  //WARNING:append myo and pur with same order everytime
  std::vector< int > elemType_list
    {elemType_myo, elemType_LVpur, elemType_RVpur}; 
  std::vector< int > phy_tag_list {phy_tag_myo, phy_tag_LVpur, phy_tag_RVpur}; 

  // Input files
  char * char_home_dir = getenv("HOME");
  std::string home_dir (char_home_dir);
  
  //test mesh 
  //std::string geo_file_myo("./myo.vtu");
  std::string geo_file_myo
    (home_dir+"/PERIGEE/examples/EP-FEA/mesh/growth_test.vtu");  
  std::string geo_file_LVpur("./pur1.vtu");
  std::string geo_file_RVpur("./pur2.vtu");
  std::string LVendnodes_file
    (home_dir+"/PERIGEE/examples/EP-FEA/mesh/endnodes.txt");
  std::string RVendnodes_file
    (home_dir+"/PERIGEE/examples/EP-FEA/mesh/endnodes.txt");
  //criteria (distance) for matching purkinje junction nodes to myocardium 
  const double LV_tol= 0.1;
  const double RV_tol= 0.1;
  
  // //heart mesh 
  // std::string geo_file_myo
  //   (home_dir+"/PERIGEE/examples/EP-FEA/mesh/HLHS_fibers.vtu");
  // std::string geo_file_LVpur
  //   (home_dir+"/PERIGEE/examples/EP-FEA/mesh/LV-purkinje.vtu");
  // std::string geo_file_RVpur
  //   (home_dir+"/PERIGEE/examples/EP-FEA/mesh/RV-purkinje.vtu");
  // //warning: check that the first node of purkinje is not in the endnodes list.
  // std::string LVendnodes_file
  //   (home_dir+"/PERIGEE/examples/EP-FEA/mesh/LV_endnodes-picked.txt");
  // std::string RVendnodes_file
  //   (home_dir+"/PERIGEE/examples/EP-FEA/mesh/RV_endnodes-picked.txt");
  // //criteria (distance) for matching purkinje junction nodes to myocardium 
  // const double LV_tol= 1.1;
  // const double RV_tol= 1.0;

  const std::string part_file("part");

  int cpu_size = 1; 
  int in_ncommon = 1;
  const bool isDualGraph = true;

  PetscMPIInt size, rank;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionString("-geo_file_myo", geo_file_myo);  
  SYS_T::GetOptionString("-geo_file_LVpur", geo_file_LVpur);
  SYS_T::GetOptionString("-geo_file_RVpur", geo_file_RVpur);
  SYS_T::GetOptionString("-LVendnodes_file", LVendnodes_file);
  SYS_T::GetOptionString("-RVendnodes_file", RVendnodes_file);

  std::cout<<"==== /Command Line Arguments ===="<<std::endl;
  std::cout<<" -geo_file_myo: "<<geo_file_myo<<std::endl;
  std::cout<<" -geo_file_LVpur: "<<geo_file_LVpur<<std::endl;
  std::cout<<" -geo_file_RVpur: "<<geo_file_RVpur<<std::endl;
  std::cout<<" -LVendnodes_file: "<<LVendnodes_file<<std::endl;
  std::cout<<" -RVendnodes_file: "<<RVendnodes_file<<std::endl;

  std::cout<<" -part_file: "<<part_file<<std::endl;
  std::cout<<" -cpu_size: "<<cpu_size<<std::endl;
  std::cout<<" -in_ncommon: "<<in_ncommon<<std::endl;
  std::cout<<" -isDualGraph: true \n";
  std::cout<<"----------------------------------\n";
  std::cout<<" dofNum: "<<dofNum<<std::endl;
  std::cout<<" dofMat: "<<dofMat<<std::endl;
  std::cout<<" elemType_myo: "<<elemType_myo<<std::endl;
  std::cout<<" elemType_LVpur: "<<elemType_LVpur<<std::endl;
  std::cout<<" elemType_RVpur: "<<elemType_RVpur<<std::endl;
  std::cout<<"====  Command Line Arguments/ ===="<<std::endl;

  // Check if the geometrical file exist on disk
  SYS_T::file_check(geo_file_myo); std::cout<<geo_file_myo<<" found. \n";
  SYS_T::file_check(geo_file_LVpur); std::cout<<geo_file_LVpur<<" found. \n";
  SYS_T::file_check(geo_file_RVpur); std::cout<<geo_file_RVpur<<" found. \n";
  SYS_T::file_check(LVendnodes_file); std::cout<<LVendnodes_file<<" found. \n";
  SYS_T::file_check(RVendnodes_file); std::cout<<RVendnodes_file<<" found. \n";


  // ----- Write the input argument into a HDF5 file
  std::cout<< "// ----- Write the input argument into a HDF5 file" << std::endl;
  hid_t cmd_file_id = H5Fcreate("preprocessor_cmd.h5",
				H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  HDF5_Writer * cmdh5w = new HDF5_Writer(cmd_file_id);

  cmdh5w->write_intScalar("cpu_size", cpu_size);
  cmdh5w->write_intScalar("in_ncommon", in_ncommon);
  cmdh5w->write_intScalar("dofNum", dofNum);
  cmdh5w->write_intScalar("dofMat", dofMat);
  cmdh5w->write_intScalar("elemType_myo", elemType_myo);
  cmdh5w->write_intScalar("elemType_LVpur", elemType_LVpur);
  cmdh5w->write_intScalar("elemType_RVpur", elemType_RVpur);  
  cmdh5w->write_string("geo_file_myo", geo_file_myo);
  cmdh5w->write_string("geo_file_LVpur", geo_file_LVpur);
  cmdh5w->write_string("geo_file_RVpur", geo_file_RVpur);  
  cmdh5w->write_string("part_file", part_file);

  delete cmdh5w; H5Fclose(cmd_file_id);
  std::cout<< "// ----- Finish writing" << std::endl;

  // Read the geometry files
  int nFunc_LVpur, nElem_LVpur, nFunc_RVpur, nElem_RVpur, nFunc_myo, nElem_myo;
  std::vector<int> vecIEN_LVpur, vecIEN_RVpur, vecIEN_myo;
  std::vector<double> ctrlPts_LVpur, ctrlPts_RVpur, ctrlPts_myo,
    ctrlPts_combined, displacement_myo, displacement_LVpur, displacement_RVpur;
  std::vector< std::vector< double > > myo_fiber;
  
  TET_T::read_vtu_grid(geo_file_myo.c_str(), nFunc_myo, nElem_myo,
		       ctrlPts_myo, vecIEN_myo, myo_fiber, displacement_myo);
  TET_T::read_purkinje_lines(geo_file_LVpur.c_str(), 
			     nFunc_LVpur, nElem_LVpur, 
			     ctrlPts_LVpur, vecIEN_LVpur);
  TET_T::read_purkinje_lines(geo_file_RVpur.c_str(), 
			     nFunc_RVpur, nElem_RVpur, 
			     ctrlPts_RVpur, vecIEN_RVpur);

  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
  //std::cout << "nelem myo" << nElem_myo <<"\n"
  //	    << "nFunc myo" << nFunc_myo <<"\n"
  //	    << "ctrl pts myo" ;
  //for (auto i = ctrlPts_myo.begin(); i != ctrlPts_myo.end(); ++i)
  //  std::cout << *i << ' ';
  //std::cout << "\n" << "vec IEN myo" ;
  //for (auto i = vecIEN_myo.begin(); i != vecIEN_myo.end(); ++i)
  //  std::cout << *i << ' ';
  //std::cout << "\n" << std::endl;
  //
  //  std::cout << "nelem RVpur" << nElem_RVpur <<"\n"
  //  	    << "nFunc RVpur" << nFunc_RVpur <<"\n"
  //  	    << "ctrl pts RVpur" ;
  //  for (auto i = ctrlPts_RVpur.begin(); i != ctrlPts_RVpur.end(); ++i)
  //    std::cout << *i << ' ';
  //  std::cout << "\n" << "vec IEN RVpur" ;
  //  for (auto i = vecIEN_RVpur.begin(); i != vecIEN_RVpur.end(); ++i)
  //    std::cout << *i << ' ';
  //  std::cout << "\n" << std::endl;
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
  
  // Generate IEN
  IIEN * IEN_myo = new IEN_Tetra_P1(nElem_myo, vecIEN_myo);
  //std::cout << "IEN myo" << std::endl;
  //IEN_myo->print_IEN();

  IIEN * IEN_LVpur = new IEN_Line_P1(nElem_LVpur, vecIEN_LVpur);
  IIEN * IEN_RVpur = new IEN_Line_P1(nElem_RVpur, vecIEN_RVpur);  
  //std::cout << "IEN LVpur" << std::endl;
  //IEN_LVpur->print_IEN();
  //std::cout << "IEN RVpur" << std::endl;
  //IEN_RVpur->print_IEN();  

  if(elemType_myo == 501)  {
    SYS_T::print_fatal_if(vecIEN_myo.size() / nElem_myo != 4, "Error: the mesh connectivity array size does not match with the element type 501. \n");
  }
  else  {
    SYS_T::print_fatal_if(1, "Error: this script doesn't support this element type. \n");
  }

  if(elemType_LVpur == 512)  {
    SYS_T::print_fatal_if(vecIEN_LVpur.size() / nElem_LVpur != 2, "Error: the mesh connectivity array size does not match with the element type 512. \n");
  }
  else  {
    SYS_T::print_fatal_if(1, "Error: this script doesn't support this element type. \n");
  }
  
  if(elemType_RVpur == 512)  {
    SYS_T::print_fatal_if(vecIEN_RVpur.size() / nElem_RVpur != 2, "Error: the mesh connectivity array size does not match with the element type 512. \n");
  }
  else  {
    SYS_T::print_fatal_if(1, "Error: this script doesn't support this element type. \n");
  }
  
  // Check the tet mesh for aspect ratio
  TET_T::tetmesh_check(ctrlPts_myo, IEN_myo, nElem_myo, 3.5);

  // Generate the mesh
  IMesh * mesh_myo = new Mesh_Tet4(nFunc_myo, nElem_myo);
  std::cout << "Myocardium Mesh:" << std::endl;
  mesh_myo -> print_mesh_info();

  IMesh * mesh_LVpur = new Mesh_Line_3D(nFunc_LVpur, nElem_LVpur);
  std::cout << "LV-Purkinje Mesh:" << std::endl;
  mesh_LVpur -> print_mesh_info();

  IMesh * mesh_RVpur = new Mesh_Line_3D(nFunc_RVpur, nElem_RVpur);
  std::cout << "RV-Purkinje Mesh:" << std::endl;
  mesh_RVpur -> print_mesh_info();

  //WARNING:APPEND MYO AND PUR WITH SAME ORDER EVERYTIME
  std::vector<IMesh *> mesh_list;
  mesh_list.clear();
  mesh_list.push_back(mesh_myo);
  mesh_list.push_back(mesh_LVpur);
  mesh_list.push_back(mesh_RVpur);  

  //std::vector<IIEN *> IEN_list;
  std::vector< std::vector< int > > IEN_list;
  IEN_list.clear();
  IEN_list.push_back(vecIEN_myo);
  IEN_list.push_back(vecIEN_LVpur);
  IEN_list.push_back(vecIEN_RVpur);  

  std::vector< std::vector<double> > ctrlPts_list;
  ctrlPts_list.clear();
  ctrlPts_list.push_back(ctrlPts_myo);
  ctrlPts_list.push_back(ctrlPts_LVpur);
  ctrlPts_list.push_back(ctrlPts_RVpur);

  std::vector< std::vector<double> > displacement_list;
  displacement_LVpur.resize(nFunc_LVpur*3);
  displacement_RVpur.resize(nFunc_RVpur*3);
  displacement_list.clear();
  displacement_list.push_back(displacement_myo);
  displacement_list.push_back(displacement_LVpur);
  displacement_list.push_back(displacement_RVpur);  
  
  // std::cout << "disp LVpur:" << std::endl;
  // VEC_T::print( displacement_LVpur );
  // std::cout << "disp RVpur:" << std::endl;
  // VEC_T::print( displacement_RVpur );
  // std::cout << "disp myo:" << std::endl;
  // VEC_T::print( displacement_myo );
  
  IIEN * IEN_combined= new IEN_Mixed ( IEN_list, mesh_list, 
   				       ctrlPts_list,
   				       displacement_list,
   				       LVendnodes_file.c_str(),
   				       RVendnodes_file.c_str(),
   				       ctrlPts_combined, LV_tol, RV_tol);

  //  std::cout << "ctrlpts combined:" << std::endl;
  //  VEC_T::print( ctrlPts_combined );
  //  
  //  VEC_T::clean( vecIEN_myo );
  //  VEC_T::clean( vecIEN_pur );
  
  //IMesh * mesh_combined = new Mesh_Mixed(nFunc_tot, nElem_tot);
  IMesh * mesh_combined = new Mesh_Mixed(mesh_list, elemType_list,
					 IEN_combined, myo_fiber, phy_tag_list); 
  mesh_combined -> print_mesh_info();  
  
  // Partition
  IGlobal_Part * global_part;
  if(cpu_size > 1) {
    global_part = new Global_Part_METIS_Mixed( cpu_size, in_ncommon,
					       isDualGraph, mesh_combined,
					       IEN_combined,"epart", "npart");
  }
  
  else if(cpu_size == 1) {
    global_part = new Global_Part_Serial( mesh_combined, "epart", "npart" );
  }
  else  {
    std::cerr<<"ERROR: wrong cpu_size: "<<cpu_size<<std::endl;
    exit(EXIT_FAILURE);
  }

  Map_Node_Index * mnindex =
    new Map_Node_Index(global_part, cpu_size, mesh_combined->get_nFunc());
  //  mnindex->print_info();
  mnindex->write_hdf5("node_mapping");

  
  //WARNING: need re-implementation of boundary conditions 
  // if there will be any present
  // ----------------------------------------------------------------
  // Setup boundary condition
  std::vector<INodalBC *> NBC_list;
  NBC_list.clear();
  NBC_list.resize( dofMat );
  
  // Use dir list to set Dirichlet BCs on these vtp surfaces
  //std::vector<std::string> dir_list;
  //dir_list.push_back(sur_file_top);

  //WARNING: need re-implementation of nodal boundary conditions implementation
  // if there will be any present. (dir_list is populated)
  NBC_list[0] = new NodalBC_Line_3D_vtp( mesh_combined->get_nFunc() );
  //NBC_list[0]->print_info();
  /////NBC_list[1] = new NodalBC_3D_vtp( dir_list, nFunc );
  /////NBC_list[2] = new NodalBC_3D_vtp( dir_list, nFunc );
  /////NBC_list[3] = new NodalBC_3D_vtp( dir_list, nFunc );
  
  //WARNING: may need re-implementation of nodal boundary conditions implementation
  // if there will be any present. (di
  //std::vector<double> inflow_outward_normal;
  //inflow_outward_normal.push_back(0.0);
  //inflow_outward_normal.push_back(0.0);
  //inflow_outward_normal.push_back(1.0);
  INodalBC * InFBC = new NodalBC_Line_3D_stimulus( mesh_combined->get_nFunc());

  //WARNING: may need re-implementation of nodal boundary conditions implementation
  // if there will be any present. (di
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
    //IPart * part = new Part_Line( mesh, global_part, mnindex, IEN,
    //				  ctrlPts, proc_rank, proc_size, dofNum,
    //				  dofMat, elemType, isPrintPartInfo );
    IPart * part = new Part_Mixed_Mesh( mesh_combined, global_part, mnindex,
					IEN_combined, ctrlPts_combined,
					proc_rank, proc_size, dofNum, dofMat, 
					isPrintPartInfo );       
    mytimer->Stop();
    //    if (proc_rank ==0) {
    //cout<<"-- proc "<<proc_rank<<"  part_mixed_mesh info. \n";
    //part->print_part_ele() ;
    //part->print_part_node();
    //part->print_part_ghost_node() ;
    //part->print_part_local_to_global() ;
    //part->print_part_LIEN(mesh_combined) ;
    //part->print_part_loadbalance_edgecut() ;
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";
    //}

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
  delete ebc;
  delete InFBC;
  std::vector<INodalBC *>::iterator it_nbc;
  for(it_nbc=NBC_list.begin(); it_nbc != NBC_list.end(); ++it_nbc) delete *it_nbc;

  delete mnindex;
  delete global_part;
  delete mesh_combined;
  delete mesh_myo;
  delete mesh_LVpur;
  delete mesh_RVpur;
  delete IEN_combined;
  delete IEN_myo;
  delete IEN_LVpur;
  delete IEN_RVpur;   
  delete mytimer;

  PetscFinalize();

  return EXIT_SUCCESS;
}

// EOF
