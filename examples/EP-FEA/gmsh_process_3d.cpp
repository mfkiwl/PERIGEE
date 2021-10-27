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
  //gmshVol.append("/PERIGEE/examples/EP-FEA/mesh/HLHS-coarse.msh");
  //gmshVol.append("/PERIGEE/examples/EP-FEA/mesh/tets_cube.msh");
  // gmshVol.append("/PERIGEE/examples/EP-FEA/mesh/niederer-1mm.msh");
  // gmshVol.append("/PERIGEE/examples/EP-FEA/mesh/niederer-05mm.msh");
  // gmshVol.append("/PERIGEE/examples/EP-FEA/mesh/niederer-025mm.msh");
  gmshVol.append("/PERIGEE/examples/EP-FEA/mesh/tets_cube.msh");
  SYS_T::GetOptionString("-gmsh_Vol", gmshVol);
  std::cout<<" -gmsh_Vol: "<<gmshVol<<std::endl;
  Gmsh_FileIO * GIOVol = new Gmsh_FileIO( gmshVol );
  GIOVol -> print_info();
  
  std::string gmshPur1 (char_home_dir);
  //gmshPur1.append("/PERIGEE/examples/EP-FEA/mesh/longline1.msh");
  //gmshPur1.append("/PERIGEE/examples/EP-FEA/mesh/twolines1.msh");
  gmshPur1.append("/PERIGEE/examples/EP-FEA/mesh/oneline1.msh");
  SYS_T::GetOptionString("-gmsh_Pur1", gmshPur1);
  std::cout<<" -gmsh_Pur1: "<<gmshPur1<<std::endl;
  Gmsh_FileIO * GIOPur1 = new Gmsh_FileIO( gmshPur1 );
  GIOPur1 -> print_info();

  std::string gmshPur2 (char_home_dir);
  //  gmshPur2.append("/PERIGEE/examples/EP-FEA/mesh/longline2.msh");
  //gmshPur2.append("/PERIGEE/examples/EP-FEA/mesh/twolines2.msh");
  gmshPur2.append("/PERIGEE/examples/EP-FEA/mesh/oneline2.msh");
  SYS_T::GetOptionString("-gmsh_Pur2", gmshPur2);
  std::cout<<" -gmsh_Pur2: "<<gmshPur2<<std::endl;
  Gmsh_FileIO * GIOPur2 = new Gmsh_FileIO( gmshPur2 );
  GIOPur2 -> print_info();

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
  
  const std::string wmname_pur1("pur1");
  GIOPur1 -> write_vtu_purkinje( wmname_pur1, isXML );
  delete GIOPur1;
  const std::string wmname_pur2("pur2");
  GIOPur2 -> write_vtu_purkinje( wmname_pur2, isXML );
  delete GIOPur2;

  PetscFinalize();
  return 0;
}

// EOF
