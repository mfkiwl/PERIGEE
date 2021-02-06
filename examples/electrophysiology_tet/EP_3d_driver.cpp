// ==================================================================
// Input:
// partition file in hdf5 format.
//
// Output:
// finite element solution vector.
//
// Note:
// This code relies on the PETSc library.
//
// Date: November 2020
// Author: Ju Liu
// Modified: Oguz Ziya Tikenogullari
// ==================================================================
#include <cmath>
#include <iomanip>
#include "Vec_Tools.hpp"
//#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss.hpp"
//#include "QuadPts_Gauss_Triangle.hpp"
#include "HDF5_PartReader.hpp"
#include "FEANode.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
#include "FEAElement_Triangle6_3D_der0.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "ALocal_NodalBC.hpp"
#include "ALocal_Inflow_NodalBC.hpp"
#include "ALocal_BC_3D.hpp"
#include "AInt_Weight.hpp"
#include "APart_Node.hpp"
#include "APart_Basic_Info.hpp"
#include "PDNTimeStep.hpp"
#include "PDNSolution_EP.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "PLocAssem_EP_3D.hpp"
#include "PGAssem_EP.hpp"
#include "PTime_Solver_EP_OperatorSplit.hpp"
#include "IonicModel_AP.hpp"

int main(int argc, char *argv[])
{
  // Number of quadrature points for tets and triangles
  // Use: 1 / 1 for linear, 2 / 1 for quadratic
  int nqp_line = 1, nqp_vertex = 0;
  //Note: nqp_vertex=1 is redundant fix this

  // Estimate of the nonzero per row for the sparse matrix 
  int nz_estimate = 300;

  // partition file base name
  std::string part_file("part");

  // Nonlinear solver parameters
  double nl_rtol = 1.0e-10;
  double nl_atol = 1.0e-10;
  double nl_dtol = 0.9;
  int nl_maxits = 20;
  int nl_refreq = 4;
  //  int nl_threshold = 4;
  
  // Time step initailization
  double initial_time = 0.0;
  double initial_step = 1.0;
  int initial_index = 0;
  double final_time = 1000.0;

  // Time solver parameters
  std::string sol_bName("SOL_");
  int ttan_renew_freq = 1;
  int sol_record_freq = 5;

  //// Restart options
  //bool is_restart = false;
  //int restart_index = 0;
  //double restart_time = 0.0;
  //double restart_step = 1.0e-3;
  //std::string restart_name = "SOL_";

  PetscMPIInt rank, size;
  // ======= PETSc Initialize =======
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // ======= Read Command Line Arguments =======
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");

  SYS_T::GetOptionInt("-nqp_line", nqp_line);
  SYS_T::GetOptionInt("-nqp_vertex", nqp_vertex);
  SYS_T::GetOptionString("-part_file", part_file);
  SYS_T::GetOptionReal("-nl_rtol", nl_rtol);
  SYS_T::GetOptionReal("-nl_atol", nl_atol);
  SYS_T::GetOptionReal("-nl_dtol", nl_dtol);
  SYS_T::GetOptionInt("-nl_maxits", nl_maxits);
  SYS_T::GetOptionInt("-nl_refreq", nl_refreq);
  SYS_T::GetOptionReal("-init_time", initial_time);
  SYS_T::GetOptionReal("-fina_time", final_time);
  SYS_T::GetOptionReal("-init_step", initial_step);
  SYS_T::GetOptionInt("-init_index", initial_index);
  SYS_T::GetOptionInt("-ttan_freq", ttan_renew_freq);
  SYS_T::GetOptionInt("-sol_rec_freq", sol_record_freq);
  SYS_T::GetOptionString("-sol_name", sol_bName);

  SYS_T::cmdPrint("-part_file:", part_file);
  SYS_T::cmdPrint("-nqp_line:", nqp_line);
  SYS_T::cmdPrint("-nqp_vertex:", nqp_vertex);
  SYS_T::cmdPrint("-nl_atol:", nl_atol); 
  SYS_T::cmdPrint("-nl_dtol:", nl_dtol); 
  SYS_T::cmdPrint("-nl_maxits:", nl_maxits);
  SYS_T::cmdPrint("-nl_refreq:", nl_refreq); 
  SYS_T::cmdPrint("-init_time:", initial_time); 
  SYS_T::cmdPrint("-init_step:", initial_step); 
  SYS_T::cmdPrint("-init_index:", initial_index); 
  SYS_T::cmdPrint("-fina_time:", final_time); 
  SYS_T::cmdPrint("-ttan_freq:", ttan_renew_freq); 
  SYS_T::cmdPrint("-sol_rec_freq:", sol_record_freq); 
  SYS_T::cmdPrint("-sol_name:", sol_bName); 

  // ======= Generate Main Data Structure =======
  SYS_T::commPrint("===> Reading mesh files ... \n");
  //HDF5_PartReader * h5reader = new HDF5_PartReader(part_file, rank);
  
  //// 1.1 Get points' coordinates 
  FEANode * fNode = new FEANode(part_file, rank);
  fNode->print_info();
  
//
//  // 1.2 Get mesh size for each local element
//  IALocal_meshSize * locmSize = new ALocal_meshSize_3D_NURBS(h5reader);
//
//  // 1.3 Get extraction operator for each local elements
//  IAExtractor * fExt = new AExtractor_3D_NURBS_xyz(h5reader);
  
  
  // 1.4 Get LIEN for each local elements
  ALocal_IEN * locIEN = new ALocal_IEN(part_file, rank);
  locIEN->print_info();  

  // 1.5 Get Global Mesh Info
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_FEM_3D(part_file,rank);
  GMIptr->print_info();
  
  // 1.6 Get partition info
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
  PartBasic->print_info();
  
  // 1.7 Get local element info
  ALocal_Elem * locElem = new ALocal_Elem(part_file, rank);
  locElem->print_info();

  // 1.8 Get local BC info
  ALocal_NodalBC * locbc = new ALocal_NodalBC(part_file, rank);
  locbc->print_info();

  ALocal_EBC * locebc = new ALocal_EBC(part_file, rank);
  locebc->print_info();

  APart_Node * pNode = new APart_Node(part_file, rank);
  pNode->print_info();

  if(size != PartBasic->get_cpu_size())
  {
    PetscPrintf(PETSC_COMM_WORLD,
        "Error: Assigned CPU number does not match the partition number. \n");
    MPI_Abort(PETSC_COMM_WORLD,1);
  }

  PetscPrintf(PETSC_COMM_WORLD, "\n===> %d processor(s) are assigned for FEM analysis. ", size);
  
  // ======= Generate Finite Element =======
  // ===== Quadrature rules =====
  SYS_T::commPrint("===> Build quadrature rules. \n");
  IQuadPts * quad_line   = new QuadPts_Gauss( nqp_line );
  //IQuadPts * quad_vertex = new QuadPts_Gauss( nqp_vertex );
  quad_line->print_info();
  //quad_vertex->print_info();

  SYS_T::commPrint("===> Build quadrature weight ... \n");
  AInt_Weight * Int_w_vol = new AInt_Weight(quad_line);
  Int_w_vol->print_info();
//  
//  // ===== Finite Element Container =====
//  SYS_T::commPrint("===> Setup element container. \n");
//
//  std::vector<FEAElement *> elemArray; 
//  elemArray.resize(locElem->get_nlocalele());
//  std::vector<FEAElement *> elemsArray; 
//  elemsArray.resize(locElem->get_nlocalele());
//  double feaelement_memsize = 0.0; clock_t elem_timer = clock();
//
//  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
//    {
//      //FEAElement * elementv	= nullptr; 
//      //FEAElement * elements	= nullptr; 
//      
//      if( GMIptr->get_elemType() == 501 ){
//	elemArray[ee] = new FEAElement_Tet4( nqp_tet ); // elem type 501
//	feaelement_memsize += elemArray[ee]->get_memory_usage(); 
//	elemsArray[ee] = new FEAElement_Triangle3_3D_der0( nqp_tri ); 
//      }
//      else if( GMIptr->get_elemType() == 502 ){
//	SYS_T::print_fatal_if( nqp_tet < 29,
//			       "Error: not enough quadrature points for tets.\n" );
//	SYS_T::print_fatal_if( nqp_tri < 13,
//			       "Error: not enough quadrature points for triangles.\n" );
//	elemArray[ee] = new FEAElement_Tet10_v2( nqp_tet ); // elem type 502
//	feaelement_memsize += elemArray[ee]->get_memory_usage(); 
//	elemsArray[ee] = new FEAElement_Triangle6_3D_der0( nqp_tri ); 
//      }
//      else SYS_T::print_fatal("Error: Element type not supported.\n");
//    }
//  
//  elem_timer = clock() - elem_timer;
//  MPI_Barrier(PETSC_COMM_WORLD);
//  SYS_T::synPrintElementInfo(locElem->get_nlocalele(), feaelement_memsize,
//      (double)elem_timer/(double)CLOCKS_PER_SEC, rank);
//  // ---------------------------------------
//
//  // ===== Generate Generalized-alpha method
//  SYS_T::commPrint("===> Genereate the Generalized-alpha time scheme ... \n");
//  //this was 0.5
//  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha(1.0,1.0,1.0);
//  tm_galpha_ptr->print_info();
//
//  //====== Ionic model setup
//  SYS_T::commPrint("===> Generate Ionic Model ... \n");
//  IonicModel * ionicmodel_ptr = new IonicModel_AP () ;
//  ionicmodel_ptr -> print_info();
//
//  //====== Local assembly pointer
//  SYS_T::commPrint("===> Initialize local assembly routine ... \n");
//  IPLocAssem * locAssem_ptr =
//    new PLocAssem_EP_3D(tm_galpha_ptr, ionicmodel_ptr,
//			GMIptr->get_nLocBas(), quadv->get_num_quadPts() );
//  
//  //std::cout << "gmi locbas" << GMIptr->get_nLocBas() << std::endl;
//  //std::cout << "elements locbas " << elements->get_nLocBas() << std::endl;
//  //std::cout << "quadv num " << quadv->get_num_quadPts() << std::endl;  
//
//  // ---------------------------------------
//
//  // ======= Solution Initialization =======
//  PDNSolution * disp = new PDNSolution_EP(pNode, fNode, locbc, 2); // make it 2
//  PDNSolution * velo = new PDNSolution_EP(pNode, fNode, locbc, 0);
//  PDNSolution * hist = new PDNSolution_EP(pNode, fNode, locbc, 0);
//
//  //if( is_restart )
//  //{
//  //  initial_index = restart_index;
//  //  initial_time  = restart_time;
//  //  initial_step  = restart_step;
//  //
//  //  SYS_T::file_exist_check(restart_name.c_str());
//  //
//  //  sol->ReadBinary(restart_name.c_str());
//  //  PetscPrintf(PETSC_COMM_WORLD, "===> Read sol from disk as a restart run... \n");
//  //  PetscPrintf(PETSC_COMM_WORLD, "     restart_name: %s \n", restart_name.c_str());
//  //  PetscPrintf(PETSC_COMM_WORLD, "     restart_time: %e \n", restart_time);
//  //  PetscPrintf(PETSC_COMM_WORLD, "     restart_index: %d \n", restart_index);
//  //  PetscPrintf(PETSC_COMM_WORLD, "     restart_step: %e \n", restart_step);
//  //}
//  
//  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);
//
//  // ============= Global Assembly pointer
//  int vpetsc_type = 0; // petsc version controller
//  PGAssem_EP * gloAssem_ptr
//    = new PGAssem_EP(locAssem_ptr, GMIptr, pNode, vpetsc_type);
//  // ============= Estimate the matrix structure
//  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_ptr, locIEN, pNode, locbc );
//  
//  gloAssem_ptr->Fix_nonzero_err_str();
//  SYS_T::commPrint("===> Matrix nonzero structure fixed ... \n");
//
//  // 2.5 Setup linear solver context
//  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc();
//  SYS_T::commPrint("===> PETSc linear solver setted up:\n");
//  SYS_T::commPrint("----------------------------------------------------------- \n");
//  lsolver->Info();
//  SYS_T::commPrint("----------------------------------------------------------- \n");
//  //PC preproc; lsolver->GetPC(&preproc);
//  //PCSetType( preproc, PCASM );
//
//  // 2.6 Assembly mass matrix and solve for consistent initial solution
//  gloAssem_ptr->Clear_KG();
//  gloAssem_ptr->Assem_mass_residual( disp,// hist,
//				     timeinfo, //ionicmodel_ptr,
//				     locElem, locAssem_ptr, locIEN, pNode,
//				     fNode, quadv, elemArray, locbc );
//
//    
//  lsolver->Solve( gloAssem_ptr->K, gloAssem_ptr->G, velo); 
//  SYS_T::commPrint("initial solution's time derivative obtained. \n"); 
//
//  // 2.7 Setup nonlinear solver context
//  PNonlinear_Solver_EP * nsolver
//    = new PNonlinear_Solver_EP(nl_rtol, nl_atol,
//					    nl_dtol, nl_maxits, nl_refreq);
//  SYS_T::commPrint("===> Nonlinear solver setted up:\n");
//  nsolver->Info();
//
//  // 2.8 Setup time marching context
//  PTime_Solver_EP_OperatorSplit * tsolver =
//    new PTime_Solver_EP_OperatorSplit( sol_bName, sol_record_freq,
//				      ttan_renew_freq, final_time );
//
//  SYS_T::commPrint("===> Time marching solver setted up:\n");
//  tsolver->Info();
//
//  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
//    tsolver->TM_generalized_alpha(
//      velo, disp, hist, timeinfo, tm_galpha_ptr, locElem, locIEN, pNode,
//      fNode, locbc, quadv, elemArray, ionicmodel_ptr, locAssem_ptr,
//      gloAssem_ptr, lsolver, nsolver );
//
//  // ======= PETSc Finalize =======
//  SYS_T::commPrint("\n===> Clean memory ... \n");
//  delete fNode;
//  delete GMIptr;
//  delete PartBasic;
//  delete locIEN;
//  delete locElem;
//  delete locbc;
//  delete quadv; delete quads;
//  delete pNode;
//  delete Int_w_vol;
//  delete disp;
//  delete velo;
//  delete hist;  
//  delete timeinfo;
//  delete tm_galpha_ptr;
//  delete ionicmodel_ptr;
//  delete locAssem_ptr;
//  delete gloAssem_ptr;
//  delete lsolver;
//  delete nsolver;
//  delete tsolver;
//  std::vector<FEAElement *>::iterator it_elema;
//  for(it_elema = elemArray.begin(); it_elema != elemArray.end(); ++it_elema) {
//    delete *it_elema;
//  }
//  ////---------------------------------------
  SYS_T::commPrint("===> Exit program. \n");
  PetscFinalize();
  return 0;
}
//EOF
