// ==================================================================
// Pre_postprocess.cpp for purkinje model with line elements
//
// ------------------------------------------------------------------
// Objective:
// This routine provides a mesh partition for all postprocessers. This
// routine is needed because the postprocessors should be run in 
// parallel. But we do not require the postprocessor use the same mesh
// partition as the analysis part, since most of the time, postprocess
// requires less cpu.
//
// Date: Dec. 10th 2013
// Author: Ju Liu
// Modified: Oguz Ziya Tikenogullari
// ==================================================================
#include "HDF5_Reader.hpp"
#include "Tet_Tools.hpp"
#include "Mesh_Mixed.hpp"
#include "Mesh_Tet4.hpp"
#include "Mesh_Line_3D.hpp"
#include "IEN_Line_P1.hpp"
#include "IEN_Tetra_P1.hpp"
#include "IEN_Mixed.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_METIS_Mixed.hpp"
#include "Global_Part_Serial.hpp"
#include "Part_Line.hpp"
#include "Part_Mixed_Mesh.hpp"

//#include "NURBS_Tools.hpp"
//#include "kRefinement.hpp"
//#include "Mesh_NURBS_1Patch_3D.hpp"
//#include "IEN_NURBS_1Patch_3D.hpp"
//#include "Global_Part_METIS.hpp"
//#include "Global_Part_Serial.hpp"
//#include "Map_Node_Index.hpp"
//#include "Part_NURBS_1Patch_3D_METIS.hpp"
//#include "HDF5_PartReader.hpp"

using namespace std;

int main(int argc, char *argv[])
{

  //warning: check that the first node of purkinje is not in the endnodes list.
  char * char_home_dir = getenv("HOME");
  std::string home_dir (char_home_dir);

  ////test mesh endnodes
  //std::string LVendnodes_file
  //  (home_dir+"/PERIGEE/examples/EP-FEA/mesh/endnodes.txt");
  //std::string RVendnodes_file
  //  (home_dir+"/PERIGEE/examples/EP-FEA/mesh/endnodes.txt");
  ////criteria (distance) for matching purkinje junction nodes to myocardium 
  //const double LV_tol= 0.1;
  //const double RV_tol= 0.1;
  
  //heart mesh endnodes.
  std::string LVendnodes_file
    (home_dir+"/PERIGEE/examples/EP-FEA/mesh/LV_endnodes-picked.txt");
  std::string RVendnodes_file
    (home_dir+"/PERIGEE/examples/EP-FEA/mesh/RV_endnodes-picked.txt");
  //  criteria (distance) for matching purkinje junction nodes to myocardium 
  const double LV_tol= 1.0;
  const double RV_tol= 1.1;
  
  int sysret = system("rm -rf postpart_p*.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");

  std::string geo_file_myo, geo_file_LVpur, geo_file_RVpur;
  int dofNum, dofMat, elemType_myo, elemType_LVpur, elemType_RVpur, in_ncommon, probDim;

  std::string part_file("postpart");
  int cpu_size = 6;
  bool isDualGraph = true;
  bool isread_part = true;

  PetscMPIInt size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  SYS_T::print_fatal_if(size!=1, "ERROR: preprocessor is a serial program! \n");

  // Read preprocessor command-line arguements recorded in the .h5 file
  hid_t prepcmd_file
    = H5Fopen("preprocessor_cmd.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  HDF5_Reader * cmd_h5r = new HDF5_Reader( prepcmd_file );
  
  cmd_h5r -> read_string("/", "geo_file_myo", geo_file_myo);
  cmd_h5r -> read_string("/", "geo_file_LVpur", geo_file_LVpur);
  cmd_h5r -> read_string("/", "geo_file_RVpur", geo_file_RVpur);
  elemType_myo = cmd_h5r -> read_intScalar("/","elemType_myo");
  elemType_LVpur = cmd_h5r -> read_intScalar("/","elemType_LVpur");
  elemType_RVpur = cmd_h5r -> read_intScalar("/","elemType_RVpur");
  dofNum = cmd_h5r -> read_intScalar("/","dofNum");
  dofMat = cmd_h5r -> read_intScalar("/","dofMat");
  //probDim = cmd_h5r -> read_intScalar("/","probDim");
  in_ncommon = cmd_h5r -> read_intScalar("/","in_ncommon");
  
  delete cmd_h5r; H5Fclose(prepcmd_file);
  
  // The user can specify the new mesh partition options
  SYS_T::GetOptionInt("-cpu_size", cpu_size);
  SYS_T::GetOptionInt("-in_ncommon", in_ncommon);
  SYS_T::GetOptionBool("-METIS_isDualGraph", isDualGraph);
  SYS_T::GetOptionBool("-isread_part", isread_part); 
  
  cout<<"==== Command Line Arguments ===="<<endl;
  cout<<" -part_file: "<<part_file<<endl;
  cout<<" -cpu_size: "<<cpu_size<<endl;
  cout<<" -in_ncommon: "<<in_ncommon<<endl;
  if(isDualGraph) cout<<" -METIS_isDualGraph: true \n";
  else cout<<" -METIS_isDualGraph: false \n";
  if(isread_part) cout<<" -isread_part: true \n";
  else cout<<" -isread_part: false \n";
  cout<<"----------------------------------\n";
  cout<<"geo_file_myo: "<<geo_file_myo<<endl;
  cout<<"geo_file_LVpur: "<<geo_file_LVpur<<endl;
  cout<<"geo_file_RVpur: "<<geo_file_RVpur<<endl;  
  //cout<<"probDim: "<<probDim<<endl;
  cout<<"dofNum: "<<dofNum<<endl;
  cout<<"elemType_myo: "<<elemType_myo<<endl;
  cout<<"elemType_LVpur: "<<elemType_LVpur<<endl;
  cout<<"elemType_RVpur: "<<elemType_RVpur<<endl;  
  cout<<"==== Command Line Arguments ===="<<endl;
  
  // Read the geo_file
  int nFunc_LVpur,  nFunc_RVpur, nElem_LVpur, nElem_RVpur, nFunc_myo, nElem_myo;
  std::vector<int> vecIEN_LVpur, vecIEN_RVpur, vecIEN_myo;
  std::vector<int> phy_tag_LVpur, phy_tag_RVpur, phy_tag_myo;
  std::vector<double> ctrlPts_LVpur, ctrlPts_RVpur, ctrlPts_myo, ctrlPts_combined;
  std::vector< std::vector< double > > myo_fiber;
  std::vector< int > elemType_list {elemType_myo, elemType_LVpur, elemType_RVpur};
  
  //int nFunc, nElem;
  //std::vector<int> vecIEN;
  //std::vector<double> ctrlPts;
  //std::vector<int> phy_tag;
  
  // Check if the given geo file exist
  SYS_T::file_exist_check( geo_file_myo.c_str() );
  SYS_T::file_exist_check( geo_file_LVpur.c_str() );
  SYS_T::file_exist_check( geo_file_RVpur.c_str() );
  
  // Warning: this function returns phy_tag as 1 only, for now.
  //TET_T::read_purkinje_lines(geo_file.c_str(), nFunc, nElem, ctrlPts, vecIEN, phy_tag);
  TET_T::read_vtu_grid(geo_file_myo.c_str(), nFunc_myo, nElem_myo,
		       ctrlPts_myo, vecIEN_myo, myo_fiber);
  TET_T::read_purkinje_lines(geo_file_LVpur.c_str(),nFunc_LVpur, nElem_LVpur, 
			     ctrlPts_LVpur, vecIEN_LVpur, phy_tag_LVpur);
  TET_T::read_purkinje_lines(geo_file_RVpur.c_str(),nFunc_RVpur, nElem_RVpur, 
			     ctrlPts_RVpur, vecIEN_RVpur, phy_tag_RVpur);  

  //std::cout<<"elemType_myo: "<<elemType_myo<<std::endl;
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

  
  if(int(ctrlPts_LVpur.size()) != nFunc_LVpur * 3) SYS_T::print_fatal("Error: the ctrlPts from geo_file does not match the given number of nodes. \n");
  if(int(ctrlPts_RVpur.size()) != nFunc_RVpur * 3) SYS_T::print_fatal("Error: the ctrlPts from geo_file does not match the given number of nodes. \n");
  if(int(ctrlPts_myo.size()) != nFunc_myo * 3) SYS_T::print_fatal("Error: the ctrlPts from geo_file does not match the given number of nodes. \n");
  
  //std::cout<<"nElem_myo: "<<nElem_myo<<std::endl;
  //std::cout<<"nElem_LVpur: "<<nElem_LVpur<<std::endl;
  //std::cout<<"nElem_RVpur: "<<nElem_RVpur<<std::endl;  
  //std::cout<<"nFunc_myo: "<<nFunc_myo<<std::endl;
  //std::cout<<"nFunc_LVpur: "<<nFunc_LVpur<<std::endl;
  //std::cout<<"nFunc_RVpur: "<<nFunc_RVpur<<std::endl;  

  // Generate IEN
  IIEN * IEN_myo = new IEN_Tetra_P1(nElem_myo, vecIEN_myo);
  //std::cout << "IEN myo" << std::endl;
  //IEN_myo->print_IEN();

  IIEN * IEN_LVpur = new IEN_Line_P1(nElem_LVpur, vecIEN_LVpur);
  //std::cout << "IEN LVpur" << std::endl;
  //IEN_LVpur->print_IEN();

  IIEN * IEN_RVpur = new IEN_Line_P1(nElem_RVpur, vecIEN_RVpur);
  //std::cout << "IEN RVpur" << std::endl;
  //IEN_RVpur->print_IEN();

  // Generate the mesh
  IMesh * mesh_myo = new Mesh_Tet4(nFunc_myo, nElem_myo);
  std::cout << "mesh myo:" << std::endl;
  mesh_myo -> print_mesh_info();

  IMesh * mesh_LVpur = new Mesh_Line_3D(nFunc_LVpur, nElem_LVpur);
  std::cout << "mesh LVpur:" << std::endl;
  mesh_LVpur -> print_mesh_info();

  IMesh * mesh_RVpur = new Mesh_Line_3D(nFunc_RVpur, nElem_RVpur);
  std::cout << "mesh RVpur:" << std::endl;
  mesh_RVpur -> print_mesh_info();

  //WARNING:append myo and pur with same order everytime
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

  //
  
  IIEN * IEN_combined= new IEN_Mixed ( IEN_list, mesh_list, 
				       ctrlPts_list,
				       LVendnodes_file.c_str(),
				       RVendnodes_file.c_str(),
				       ctrlPts_combined, LV_tol, RV_tol);

  
  //std::cout << "ctrlpts combined:" << std::endl;
  //VEC_T::print( ctrlPts_combined );
  ////std::cout << "elemType combined:" << std::endl;
  ////VEC_T::print( elemType_combined );
  //VEC_T::clean( vecIEN_myo );
  //VEC_T::clean( vecIEN_pur );

  IMesh * mesh_combined = new Mesh_Mixed(mesh_list, elemType_list,
					 IEN_combined, myo_fiber);
  mesh_combined -> print_mesh_info();

  //int db_elem = 47395;
  //std::cout << "nodes of element :" << db_elem << "\n" 
  //	    << IEN_combined->get_IEN(db_elem, 0) << ", "
  //	    << IEN_combined->get_IEN(db_elem, 1) <<  std::endl;
  
  
  // Partition
  IGlobal_Part * global_part;
  if(cpu_size > 1) {
    //global_part = new Global_Part_METIS( cpu_size, in_ncommon,
    //					 isDualGraph, mesh, IEN,
    //					 "epart", "npart" );
    global_part = new Global_Part_METIS_Mixed( cpu_size, in_ncommon,
					       isDualGraph, mesh_combined,
					       IEN_combined,"epart", "npart");
  }
  
  else if(cpu_size == 1) {
    global_part = new Global_Part_Serial( mesh_combined, "epart", "npart" );
    //std::cerr<<"ERROR: complete serial parting: "<<cpu_size<<std::endl;
  }
  else
    {
      std::cerr<<"ERROR: wrong cpu_size: "<<cpu_size<<std::endl;
      exit(EXIT_FAILURE);
    }

  Map_Node_Index * mnindex =
    new Map_Node_Index(global_part, cpu_size, mesh_combined->get_nFunc());
  //mnindex->print_info();
  mnindex->write_hdf5("post_node_mapping");
  
  cout<<"\n=== Start Partition ... \n";
  int proc_size = cpu_size; bool isPrintPartInfo = true;
  SYS_T::Timer * mytimer = new SYS_T::Timer();
  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)    {
    mytimer->Reset(); mytimer->Start();
    
    //IPart * part = new Part_Line( mesh, global_part, mnindex, IEN,
    //				  ctrlPts, proc_rank, proc_size, dofNum, elemType,
    //				  isPrintPartInfo );
    
    IPart * part = new Part_Mixed_Mesh( mesh_combined, global_part, mnindex,
					IEN_combined, ctrlPts_combined, proc_rank,
					proc_size, dofNum, dofMat, 
					isPrintPartInfo );
    
    part->write(part_file.c_str());
    mytimer->Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";
    delete part;
  }
  
  // Clean memory
  cout<<"\n=== Clean memory. \n";
  delete mnindex; delete global_part; delete mytimer;
  delete mesh_myo;delete mesh_LVpur; delete mesh_RVpur; delete mesh_combined;
  delete IEN_myo;delete IEN_LVpur;delete IEN_RVpur;delete IEN_combined;
  PetscFinalize();
  return EXIT_SUCCESS;

}
