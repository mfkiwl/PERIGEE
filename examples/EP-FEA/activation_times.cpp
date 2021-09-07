// ==================================================================
// activation_times.cpp 
// ------------------------------------------------------------------
// This calculates activation times and EP dissynchrony index
//
// Date: Aug 2021
// Author: Ju Liu
// Modified: Oguz Ziya Tikenogullari
// ==================================================================
//#include <cmath>
//#include <sstream>

#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "APart_Basic_Info.hpp"
#include "APart_Node.hpp"
#include "QuadPts_vis_tet4.hpp"
#include "QuadPts_vis.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Line2_3D_der1.hpp"
#include "VisDataPrep_3DEP.hpp"
#include "VTK_Writer_EP_Mixed.hpp"

#include "AGlobal_Mesh_Info_Mixed.hpp"
#include "ALocal_Elem_Fiber.hpp"
#include "ALocal_IEN_Mixed.hpp"


using namespace std;

int main( int argc, char * argv[] )
{
  // say that tissue is activated beyond this value.
  double act_threshold= -40;
  
  int nqp_line = 2;
  int nqp_tet = 4;
  
  string element_part_file = "epart.h5";
  string anode_mapping_file = "node_mapping.h5"; 
  string pnode_mapping_file = "post_node_mapping.h5";
  
  const std::string part_file("postpart");
  
  std::string sol_bname("SOL_");

  std::string out_bname("ACT_");

  const int dof = 1;

  int time_start = 0;
  int time_step = 50;
  int time_end = 10000;
  double dt = 0.1;

  bool isXML = true;
  bool is_write_fibers = false;
  bool isRestart = false;

  PetscMPIInt rank, size;
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionReal("-dt", dt);
  SYS_T::GetOptionString("-sol_bname", sol_bname);
  SYS_T::GetOptionString("-out_bname", out_bname);
  SYS_T::GetOptionBool("-xml", isXML);
  SYS_T::GetOptionBool("-write_fibers", is_write_fibers);
  SYS_T::GetOptionBool("-restart", isRestart);

  SYS_T::cmdPrint("-sol_bname:", sol_bname);
  SYS_T::cmdPrint("-out_bname:", out_bname);
  SYS_T::cmdPrint("-time_start:", time_start);
  SYS_T::cmdPrint("-time_step:", time_step);
  SYS_T::cmdPrint("-time_end:", time_end);
  SYS_T::cmdPrint("-dt:",dt);
//
  if(isXML) PetscPrintf(PETSC_COMM_WORLD, "-xml: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-xml: false \n");
  if(is_write_fibers) PetscPrintf(PETSC_COMM_WORLD, "-write_fibers: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-write_fibers: false \n");

  if(isRestart) PetscPrintf(PETSC_COMM_WORLD, "-restart: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-restart: false \n");

  // If this is not a restart run, clean all previous visualization files
  if( !isRestart )
  {
    int sysret = system("rm -rf *ACT_*.vtu");
    SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
    sysret = system("rm -rf *ACT_*.pvtu");
    SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
    sysret = system("rm -rf *ACT_*.pvd");
    SYS_T::print_fatal_if(sysret != 0, "Error: system call failed. \n");
  }

  SYS_T::commPrint("===> Reading mesh files ... ");
  FEANode * fNode = new FEANode(part_file, rank);
  ALocal_IEN_Mixed * locIEN = new ALocal_IEN_Mixed(part_file, rank);
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_Mixed(part_file,rank);
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
  ALocal_Elem * locElem = new ALocal_Elem_Fiber(part_file, rank);
  APart_Node * pNode = new APart_Node(part_file, rank);
  SYS_T::commPrint("Done! \n");

  if(size != PartBasic->get_cpu_size()) SYS_T::print_fatal(
      "Error: number of processors does not match with prepost! \n");

  PetscPrintf(PETSC_COMM_WORLD,
      "\n===> %d processor(s) are assigned for:", size);
  PetscPrintf(PETSC_COMM_WORLD, "Postprocessing - activation times.\n");

  SYS_T::commPrint("===> Build sampling points.");

  //int nLocBas = GMIptr->get_nLocBas();

  // ===== Quadrature rules =====
  SYS_T::commPrint("===> Build quadrature rules. \n");
  SYS_T::commPrint("===> and Setup element container. \n");

  std::vector<IQuadPts *> quadArray; 
  quadArray.resize(locElem->get_nlocalele());
  std::vector<FEAElement *> elemArray; 
  elemArray.resize(locElem->get_nlocalele());
  int local_elemType;

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee) {

    local_elemType = GMIptr->get_elemType(locElem->get_elem_loc(ee));

    if ( local_elemType == 512 ){
      quadArray[ee] = new QuadPts_vis( nqp_line );
      elemArray[ee] = new FEAElement_Line2_3D_der1( nqp_line ); 
    }
    else if ( local_elemType == 501 ){
      quadArray[ee] = new QuadPts_vis_tet4( );
      elemArray[ee] = new FEAElement_Tet4( nqp_tet ); 
    } else {
      SYS_T::print_fatal("Error: Element type not supported.\n");      
    }
  }
  
  IVisDataPrep * visprep = new VisDataPrep_3DEP("activation time");
  //visprep->print_info();

  // Allocate memory for pointArrays and activation times
  double ** pointArrays = new double * [visprep->get_arrayCompSize()];
  double ** activation_times = new double * [visprep->get_arrayCompSize()];

  for(int ii=0; ii<visprep->get_arrayCompSize(); ++ii){
    pointArrays[ii] = new double [pNode->get_nlocghonode() * visprep->get_arraySizes(ii)];
    activation_times[ii] = new double [pNode->get_nlocghonode() * visprep->get_arraySizes(ii)];
    //initialize activation times array to -100.
    for(int jj=0; jj< pNode->get_nlocghonode() * visprep->get_arraySizes(ii) ; ++jj) {
      activation_times[ii][jj] = -100.0;
    }
  }

  VTK_Writer_EP_Mixed * vtk_w = new VTK_Writer_EP_Mixed( GMIptr,  element_part_file);
  
  std::ostringstream time_index;

  for(int time = time_start; time<=time_end; time+= time_step)  {
    std::string name_to_read(sol_bname);
    time_index.str("");
    time_index<< 900000000 + time;
    name_to_read.append(time_index.str());

    PetscPrintf(PETSC_COMM_WORLD, "Time %d: Read %s \n",
		time, name_to_read.c_str() );
    
    visprep->get_pointArray(name_to_read, anode_mapping_file,
			    pnode_mapping_file, pNode, GMIptr,
			    dof, pointArrays);
    
    for (int ii=0 ; ii < pNode->get_nlocghonode() ; ++ii){
      //std::cout << "pointArrays " << pointArrays[0][ii] <<std::endl;
      
      if ((pointArrays[0][ii] > act_threshold) && (activation_times[0][ii] <0.0 )) {
	activation_times[0][ii] = time*dt;
      }
    }
    		
  }

  std::string name_to_write(out_bname);

  vtk_w->writeOutput_compact( GMIptr, fNode, locIEN, locElem,
			      visprep, elemArray, quadArray, activation_times, rank, size, 
			      pNode -> get_ntotalnode(), 0.0, sol_bname, 
			      out_bname, name_to_write, isXML, is_write_fibers);

  MPI_Barrier(PETSC_COMM_WORLD);

  // Finalize
  for(int ii=0; ii<visprep->get_ptarray_size(); ++ii)
    delete [] pointArrays[ii];

  std::vector<FEAElement *>::iterator it_elema;
  for(it_elema = elemArray.begin(); it_elema != elemArray.end(); ++it_elema) {
    delete *it_elema;
  }

  std::vector<IQuadPts *>::iterator it_quada;
  for(it_quada = quadArray.begin(); it_quada != quadArray.end(); ++it_quada) {
    delete *it_quada;
  }
  
  delete [] pointArrays;
  delete [] activation_times;
  delete vtk_w;
  delete visprep;
  delete pNode; delete locElem; delete PartBasic; delete GMIptr;
  delete locIEN; delete fNode;

  SYS_T::commPrint("===> FIN. \n");
  PetscFinalize();
  return 0;
}

// EOF
