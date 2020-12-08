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
  char * char_home_dir = getenv("HOME");
  std::string gmshFile (char_home_dir);
  gmshFile.append("/PERIGEE/examples/electrophysiology_tet/mesh/cube.msh");
  //gmshFile.append("/PERIGEE/examples/electrophysiology_tet/single_tet.msh");  

  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  SYS_T::GetOptionString("-gmsh_file", gmshFile);
  std::cout<<" -gmsh_file: "<<gmshFile<<std::endl;

  Gmsh_FileIO * GIO = new Gmsh_FileIO( gmshFile );

  GIO -> print_info();

  GIO -> write_vtp(0,0);//top
  GIO -> write_vtp(1,0);//bot
  GIO -> write_vtp(2,0);//lef
  GIO -> write_vtp(3,0);//rig
  GIO -> write_vtp(4,0);//fro
  GIO -> write_vtp(5,0);//bac

  //GIO -> write_vtp(0,0);//back
  //GIO -> write_vtp(1,0);//left
  //GIO -> write_vtp(2,0);//diagface
  //GIO -> write_vtp(3,0);//bottom

  const std::string wmname("whole_vol");
  const bool isXML = true;
  GIO -> write_vtu( wmname, isXML );

  delete GIO; 
  PetscFinalize();
  return 0;
}

// EOF
