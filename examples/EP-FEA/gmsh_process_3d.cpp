// ==================================================================
// gmsh_process.cpp
//
// Code that handles the gmsh file and write data into vtk files.
//
// Author: Ju Liu
// Modified: Oguz Ziya Tikenogullari
// ==================================================================
#include "Gmsh_FileIO.hpp"

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  char * char_home_dir = getenv("HOME");

  std::string gmshVol (char_home_dir);
  gmshVol.append("/PERIGEE/examples/EP-FEA/mesh/tets_cube.msh");
  //gmshVol.append("/PERIGEE/examples/EP-FEA/mesh/single_tet.msh");
  //gmshVol.append("/niederer_meshes/niederer-01mm.msh");
  //gmshVol.append("/niederer_meshes/niederer-02mm.msh");
  //gmshVol.append("/niederer_meshes/niederer-05mm.msh");
  //gmshVol.append("/PERIGEE/examples/EP-FEA/mesh/beam.msh");
  SYS_T::GetOptionString("-gmsh_Vol", gmshVol);
  std::cout<<" -gmsh_Vol: "<<gmshVol<<std::endl;
  Gmsh_FileIO * GIOVol = new Gmsh_FileIO( gmshVol );
  GIOVol -> print_info();

  std::string gmshVol2 (char_home_dir);
  gmshVol2.append("/PERIGEE/examples/EP-FEA/mesh/tets_cube_grown.msh");
  //gmshVol2.append("/PERIGEE/examples/EP-FEA/mesh/single_tet_grown.msh");
  //gmshVol2.append("/niederer_meshes/niederer-01mm.msh");
  //gmshVol2.append("/niederer_meshes/niederer-02mm.msh");
  //gmshVol2.append("/niederer_meshes/niederer-05mm.msh");
  //gmshVol2.append("/PERIGEE/examples/EP-FEA/mesh/beam.msh");
  SYS_T::GetOptionString("-gmsh_Vol2", gmshVol2);
  std::cout<<" -gmsh_Vol2: "<<gmshVol2<<std::endl;
  Gmsh_FileIO * GIOVol2 = new Gmsh_FileIO( gmshVol2 );
  GIOVol2 -> print_info();
  
  std::string gmshPur1 (char_home_dir);
  //gmshPur1.append("/PERIGEE/examples/EP-FEA/mesh/longline1.msh");
  //gmshPur1.append("/PERIGEE/examples/EP-FEA/mesh/twolines1.msh");
  //gmshPur1.append("/PERIGEE/examples/EP-FEA/mesh/threelines.msh");
  gmshPur1.append("/PERIGEE/examples/EP-FEA/mesh/oneline1.msh");
  SYS_T::GetOptionString("-gmsh_Pur1", gmshPur1);
  std::cout<<" -gmsh_Pur1: "<<gmshPur1<<std::endl;
  Gmsh_FileIO * GIOPur1 = new Gmsh_FileIO( gmshPur1 );
  GIOPur1 -> print_info();

  std::string gmshPur2 (char_home_dir);
  //gmshPur2.append("/PERIGEE/examples/EP-FEA/mesh/longline2.msh");
  //gmshPur2.append("/PERIGEE/examples/EP-FEA/mesh/twolines2.msh");
  gmshPur2.append("/PERIGEE/examples/EP-FEA/mesh/oneline2.msh");
  SYS_T::GetOptionString("-gmsh_Pur2", gmshPur2);
  std::cout<<" -gmsh_Pur2: "<<gmshPur2<<std::endl;
  Gmsh_FileIO * GIOPur2 = new Gmsh_FileIO( gmshPur2 );
  GIOPur2 -> print_info();

  //std::string gmshPur3 (char_home_dir);
  //gmshPur3.append("/PERIGEE/examples/EP-FEA/mesh/oneline1_grown.msh");
  //SYS_T::GetOptionString("-gmsh_Pur3", gmshPur3);
  //std::cout<<" -gmsh_Pur3: "<<gmshPur3<<std::endl;
  //Gmsh_FileIO * GIOPur3 = new Gmsh_FileIO( gmshPur3 );
  //GIOPur3 -> print_info();
  //
  //std::string gmshPur4 (char_home_dir);
  //gmshPur4.append("/PERIGEE/examples/EP-FEA/mesh/twolines2_grown.msh");
  //SYS_T::GetOptionString("-gmsh_Pur4", gmshPur4);
  //std::cout<<" -gmsh_Pur4: "<<gmshPur4<<std::endl;
  //Gmsh_FileIO * GIOPur4 = new Gmsh_FileIO( gmshPur4 );
  //GIOPur4 -> print_info();

  ////three lines
  //GIO -> write_vtp_purkinje(0,0,true);//tip0
  //GIO -> write_vtp_purkinje(1,0,true);//tip1

  //  //tet and line mesh  
  //  GIO -> write_vtp(0,0);//tip
  //  GIO -> write_vtp(1,0);//pur
  //  GIO -> write_vtp(2,0);//diag
  //  GIO -> write_vtp(3,0);//vol
  //  
  //  ////faces for HLHS myocardium mesh
  //  //GIO -> write_vtp(0,0);//RV
  //  //GIO -> write_vtp(1,0);//Base
  //  //GIO -> write_vtp(2,0);//Epi
  //  //GIO -> write_vtp(3,0);//LV
  //
  //  //faces for Beam or cube mesh
  //  //GIO -> write_vtp(0,0);//top
  //  //GIO -> write_vtp(1,0);//bot
  //  //GIO -> write_vtp(2,0);//lef
  //  //GIO -> write_vtp(3,0);//rig
  //  //GIO -> write_vtp(4,0);//fro
  //  //GIO -> write_vtp(5,0);//bac
  //
  //  //faces for single tet mesh 
  //  //GIO -> write_vtp(0,0);//back
  //  //GIO -> write_vtp(1,0);//left
  //  //GIO -> write_vtp(2,0);//diagface
  //  //GIO -> write_vtp(3,0);//bottom

  const bool isXML = true;
  
  const std::string wmname_vol("myo");
  GIOVol -> write_vtu( wmname_vol, isXML );
  delete GIOVol;

  const std::string wmname_vol2("myo2");
  GIOVol2 -> write_vtu( wmname_vol2, isXML );
  delete GIOVol2;
  
  const std::string wmname_pur1("pur1");
  GIOPur1 -> write_vtu_purkinje( wmname_pur1, isXML );
  delete GIOPur1;
  const std::string wmname_pur2("pur2");
  GIOPur2 -> write_vtu_purkinje( wmname_pur2, isXML );
  delete GIOPur2;
  //const std::string wmname_pur3("pur3");
  //GIOPur3 -> write_vtu_purkinje( wmname_pur3, isXML );
  //delete GIOPur3;
  //const std::string wmname_pur4("pur4");
  //GIOPur4 -> write_vtu_purkinje( wmname_pur4, isXML );
  //delete GIOPur4;  

  PetscFinalize();
  return 0;
}

// EOF
