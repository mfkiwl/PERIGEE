// ==================================================================
// gmsh_process.cpp
//
// Code that handles the gmsh file and write data into vtk files.
//
// Date Created: Dec. 6 2017
// Author: Ju Liu
// Modified: Oguz Ziya Tikenogullari
// ==================================================================
#include "Gmsh_FileIO.hpp"

int main( int argc, char * argv[] )
{
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  char * char_home_dir = getenv("HOME");

  std::string gmshVol (char_home_dir);
  gmshVol.append("/PERIGEE/examples/electrophysiology_tet/mesh/HLHS-coarse.msh");
  //gmshVol.append("/PERIGEE/examples/electrophysiology_tet/mesh/tets_cube.msh");
  SYS_T::GetOptionString("-gmsh_Vol", gmshVol);
  std::cout<<" -gmsh_Vol: "<<gmshVol<<std::endl;
  Gmsh_FileIO * GIOVol = new Gmsh_FileIO( gmshVol );
  GIOVol -> print_info();
  
  //std::string gmshPur (char_home_dir);
  //gmshPur.append("/PERIGEE/examples/electrophysiology_tet/mesh/3lines.msh");
  //SYS_T::GetOptionString("-gmsh_Pur", gmshPur);
  //std::cout<<" -gmsh_Pur: "<<gmshPur<<std::endl;
  //Gmsh_FileIO * GIOPur = new Gmsh_FileIO( gmshPur );
  //GIOPur -> print_info();

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
  
  //const std::string wmname_pur("purkinje");
  //GIOPur -> write_vtu_purkinje( wmname_pur, isXML );
  //delete GIOPur;

  PetscFinalize();
  return 0;
}

// EOF
