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
//#include "Mesh_Tet4.hpp"
#include "Mesh_Line_3D.hpp"
//#include "IEN_Tetra_P1.hpp"
#include "IEN_Line_P1.hpp"
#include "Global_Part_METIS.hpp"
#include "Global_Part_Serial.hpp"
//#include "Part_Tet.hpp"
#include "Part_Line.hpp"

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
  int sysret = system("rm -rf postpart_p*.h5");
  SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");

  std::string geo_file;
  int dofNum, elemType, in_ncommon, probDim;

  std::string part_file("postpart");
  int cpu_size = 1;
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

  cmd_h5r -> read_string("/", "geo_file", geo_file);
  elemType = cmd_h5r -> read_intScalar("/","elemType");
  dofNum = cmd_h5r -> read_intScalar("/","dofNum");
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
  cout<<"geo_file: "<<geo_file<<endl;
  //cout<<"probDim: "<<probDim<<endl;
  cout<<"dofNum: "<<dofNum<<endl;
  cout<<"elemType: "<<elemType<<endl;
  cout<<"==== Command Line Arguments ===="<<endl;

  // Read the geo_file
  int nFunc, nElem;
  std::vector<int> vecIEN;
  std::vector<double> ctrlPts;
  std::vector<int> phy_tag;

  // Check if the given geo file exist
  SYS_T::file_exist_check( geo_file.c_str() );

  TET_T::read_purkinje_lines(geo_file.c_str(), nFunc, nElem, ctrlPts, vecIEN, phy_tag);

    if(elemType == 512)
  {
    SYS_T::print_fatal_if(vecIEN.size() / nElem != 2, "Error: the IEN from geo_file does not match the given number of element. \n");
  }
  else
  {
    SYS_T::print_fatal_if(1, "Error: this script doesn't support this element type. \n");
  }

  if(int(ctrlPts.size()) != nFunc * 3) SYS_T::print_fatal("Error: the ctrlPts from geo_file does not match the given number of nodes. \n");

  std::cout<<"nElem: "<<nElem<<std::endl;
  std::cout<<"nFunc: "<<nFunc<<std::endl;

  IIEN * IEN = new IEN_Line_P1(nElem, vecIEN);
  //IEN->print_IEN();
  
  VEC_T::clean(vecIEN);

  IMesh * mesh = new Mesh_Line_3D(nFunc, nElem);
  mesh -> print_mesh_info();

  IGlobal_Part * global_part;
  if(cpu_size > 1)
    global_part = new Global_Part_METIS( cpu_size, in_ncommon,
        isDualGraph, mesh, IEN, "post_epart", "post_npart" );
  else if(cpu_size == 1)
    global_part = new Global_Part_Serial( mesh, "post_epart", "post_npart" );
  else
  {
    cerr<<"ERROR: wrong cpu_size: "<<cpu_size<<endl;
    exit(EXIT_FAILURE);
  }

  Map_Node_Index * mnindex
    = new Map_Node_Index(global_part, cpu_size, mesh->get_nFunc());
  mnindex->write_hdf5("post_node_mapping");

  
  cout<<"\n=== Start Partition ... \n";
  int proc_size = cpu_size; bool isPrintPartInfo = true;
  SYS_T::Timer * mytimer = new SYS_T::Timer();
  for(int proc_rank = 0; proc_rank < proc_size; ++proc_rank)
  {
    mytimer->Reset(); mytimer->Start();
    IPart * part = new Part_Line( mesh, global_part, mnindex, IEN,
        ctrlPts, proc_rank, proc_size, dofNum, elemType,
        isPrintPartInfo );
    part->write(part_file.c_str());
    mytimer->Stop();
    cout<<"-- proc "<<proc_rank<<" Time taken: "<<mytimer->get_sec()<<" sec. \n";
    delete part;
  }

  // Clean memory
  cout<<"\n=== Clean memory. \n";
  delete mnindex; delete global_part; delete mesh; delete IEN;
  PetscFinalize();
  return EXIT_SUCCESS;
}
