// ==================================================================
// !!!WARNING: set up the linear solver in a better way, initial velo
// calculation isn't converging due to mass assembly, and solver
// convergence is taking a lot of iterations through time stepping
//
// Input:
// partition file in hdf5 format.
//
// Output:
// finite element solution vector.
//
// Note:
// This code relies on the PETSc library.
//
// IMPROVEMENT:
// individual local assembly instantiation is not necessary for each element.
// create one for each mesh type: line and tet elements.
// later assign these to the corresponding entries in the locassem_array .
// take the quadarray as an example. Elemarray needs to contain  individial
// instatntiations for each element.
//
// Date: November 2020
// Author: Ju Liu
// Modified: Oguz Ziya Tikenogullari
// ==================================================================
#include <cmath>
#include <iomanip>
#include <numeric>

#include "Vec_Tools.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss.hpp"
#include "QuadPts_Gauss_Triangle.hpp"
#include "HDF5_PartReader.hpp"
#include "FEANode.hpp"
#include "FEAElement_Tet4.hpp"
//#include "FEAElement_Tet10_v2.hpp"
#include "FEAElement_Triangle3_3D_der0.hpp"
//#include "FEAElement_Triangle6_3D_der0.hpp"
#include "FEAElement_Line2_3D_der1.hpp"
#include "AGlobal_Mesh_Info_FEM_3D.hpp"
#include "AGlobal_Mesh_Info_Mixed.hpp"
#include "ALocal_Elem_Fiber.hpp"
#include "ALocal_IEN.hpp"
#include "ALocal_IEN_Mixed.hpp"
#include "ALocal_NodalBC.hpp"
#include "ALocal_Inflow_NodalBC.hpp"
#include "ALocal_BC_3D.hpp"
//#include "AInt_Weight.hpp"
#include "APart_Node.hpp"
#include "APart_Basic_Info.hpp"
#include "PDNTimeStep.hpp"
#include "PDNSolution_EP.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "PLocAssem_EP_3D.hpp"
#include "PGAssem_EP.hpp"
#include "PTime_Solver_EP_OperatorSplit.hpp"
#include "IonicModel_AP.hpp"
#include "IonicModel_TTP.hpp"
#include "IonicModel_Purkinje.hpp"
#include "IonicModel_Passive.hpp"

int main(int argc, char *argv[])
{
  // Number of quadrature points for tets and triangles
  // Use: 1 / 1 for linear, 2 / 1 for quadratic
  int nqp_line = 1, nqp_vertex = 0;
  int nqp_tet = 4, nqp_tri = 3;
  //Note: nqp_vertex=1 is redundant fix this

  //// Estimate of the nonzero per row for the sparse matrix 
  //int nz_estimate = 300;

  // partition file base name
  std::string part_file("part");

  // Nonlinear solver parameters
  double nl_rtol = 1.0e-10;
  double nl_atol = 1.0e-10;
  double nl_dtol = 0.9;
  int nl_maxits = 20;
  int nl_refreq = 5;
  //  int nl_threshold = 4;
  
  // Time step initailization
  double initial_time = 0.0;
  double initial_step = 0.1;
  int initial_index = 0;
  double final_time = 100.0;

  // Time solver parameters
  std::string sol_bName("SOL_");
  int ttan_renew_freq = 1;
  int sol_record_freq = 1;

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
  SYS_T::GetOptionInt("-nqp_tet", nqp_tet);
  SYS_T::GetOptionInt("-nqp_tri", nqp_tri);
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
  SYS_T::cmdPrint("-nqp_tet:", nqp_tet);
  SYS_T::cmdPrint("-nqp_tri:", nqp_tri);
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

  //// 1.1 Get points' coordinates
  SYS_T::commPrint("===> FEANode ... \n");
  FEANode * fNode = new FEANode(part_file, rank);
  //if (rank==1)
  //  fNode->print_info();
    
  // 1.4 Get LIEN for each local elements
  SYS_T::commPrint("===> ALocal_IEN_Mixed ... \n");
  ALocal_IEN_Mixed * locIEN = new ALocal_IEN_Mixed(part_file, rank);
  //if (rank==1)
  //  locIEN->print_info();
  
  // 1.5 Get Global Mesh Info
  SYS_T::commPrint("===> AGlobal_Mesh_Info_Mixed ... \n");
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_Mixed(part_file,rank);
  //if (rank==1)
  //  GMIptr->print_info();

  // 1.6 Get partition info
  SYS_T::commPrint("===> APart_Basic_Info ... \n");
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
  //if (rank==1)
  //  PartBasic->print_info();

  
  // 1.7 Get local element info
  SYS_T::commPrint("===> ALocal_Elem ... \n");
  ALocal_Elem * locElem = new ALocal_Elem_Fiber(part_file, rank);
  //  if (rank==1)
  //locElem->print_info();

  // 1.8 Get local BC info
  SYS_T::commPrint("===> ALocal_NodalBC ... \n");
  ALocal_NodalBC * locbc = new ALocal_NodalBC(part_file, rank);
  //locbc->print_info();
  
  //ALocal_EBC * locebc = new ALocal_EBC(part_file, rank);
  ////locebc->print_info();

  SYS_T::commPrint("===> APart_Node ... \n");
  APart_Node * pNode = new APart_Node(part_file, rank);
  //if (rank==1)
  //  pNode->print_info();
    
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
  //quad_line->print_info();
  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );
  IQuadPts * quads = new QuadPts_Gauss_Triangle( nqp_tri );
  //quadv->print_info();
  //quads->print_info();

  //SYS_T::commPrint("===> Build quadrature weight ... \n");
  //AInt_Weight * Int_w_line = new AInt_Weight(quad_line);
  //Int_w_line->print_info();
  //AInt_Weight * Int_w_v = new AInt_Weight(quadv);
  //Int_w_v->print_info();

  std::vector<IQuadPts *> quadArray; 
  quadArray.resize(locElem->get_nlocalele());
  int local_elemType;

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee) {

    local_elemType = GMIptr->get_elemType(locElem->get_elem_loc(ee));

    if ( local_elemType == 512 ){
      quadArray[ee] =quad_line;
    }
    else if ( local_elemType == 501 ){
      quadArray[ee] =quadv;
    }
  }
    
  // ===== Finite Element Container =====
  SYS_T::commPrint("===> Setup element container. \n");

  std::vector<FEAElement *> elemArray; 
  elemArray.resize(locElem->get_nlocalele());
  double feaelement_memsize = 0.0; clock_t elem_timer = clock();

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee)
    {
      //FEAElement * elementv	= nullptr; 
      //FEAElement * elements	= nullptr;
     
      local_elemType = GMIptr->get_elemType(locElem->get_elem_loc(ee));
     
      if ( local_elemType == 512 ){
	elemArray[ee] = new FEAElement_Line2_3D_der1( nqp_line ); // elem type 512
	feaelement_memsize += elemArray[ee]->get_memory_usage(); 

      } else if( local_elemType == 501 ){
	elemArray[ee] = new FEAElement_Tet4( nqp_tet ); // elem type 501
	feaelement_memsize += elemArray[ee]->get_memory_usage(); 

      }
      else std::cout<<"Error: Element type not supported.\n";
    }


  //check if fiber orientations are nonzero for tet, zero for purkinje mesh
  std::vector<double> fiber_ori_e;
  int nElem = locElem->get_nlocalele();
  for(int ee=0; ee<nElem; ++ee) {
    
    locElem->get_fiber_ori_e(fiber_ori_e, ee);
	  
    if (((elemArray.at(ee))->get_elemDim() ==3) &&
	(std::accumulate(fiber_ori_e.begin(), fiber_ori_e.end(), 0.0)==0.0)){
      std::cout << "element no: " << ee << std::endl;
      std::cout << "element fiber ori: " <<  std::endl;
      VEC_T::print(fiber_ori_e);
      SYS_T::print_exit("Error: myocardium element received a zero fiber orientation \n");
      
    } else if (((elemArray.at(ee))->get_elemDim() ==1) &&
	       (std::accumulate(fiber_ori_e.begin(), fiber_ori_e.end(), 0.0)!=0.0)){
      SYS_T::print_exit("Error: purkinje element received a non-zero fiber orientation \n");
    }
  }
  //end of check 


  elem_timer = clock() - elem_timer;
  MPI_Barrier(PETSC_COMM_WORLD);
  SYS_T::synPrintElementInfo(locElem->get_nlocalele(), feaelement_memsize,
			     (double)elem_timer/(double)CLOCKS_PER_SEC, rank);
  // ---------------------------------------

  // ===== Generate Generalized-alpha method
  SYS_T::commPrint("===> Genereate the Generalized-alpha time scheme ... \n");
  //this was 0.5
  TimeMethod_GenAlpha * tm_galpha_ptr = new TimeMethod_GenAlpha(1.0,1.0,1.0);
  tm_galpha_ptr->print_info();

  //====== Ionic model setup
  //warning: I allocate space in pdnsolution for the internal variables,
  // using the largest number of internal variables required by the assigned
  // ionic models. for example if you have ttp and ap assigned as ionic models,
  // pdnsolution is allocated for
  //       number_of_int_vars_for_ttp x number_of_total_nodes.
  // for a more efficient memory use, take the nonuniformity in number of
  // internal varibles into accout. but this will complicate retrieveing the
  // values of local and element internal variables.
  // (because a differing sized vectro will need to be retrieved each time and 
  //  that size will need to be tracked.)
  // Idea: implement stride vector in PDNsolution class. 
  SYS_T::commPrint("===> Generate Ionic Models of myocardium and purkinje ... \n");
  IonicModel * ionicmodel_myo = new IonicModel_TTP () ;
  //IonicModel * ionicmodel_myo = new IonicModel_AP () ;
  //IonicModel * ionicmodel_pur = new IonicModel_TTP () ;
  //IonicModel * ionicmodel_pur = new IonicModel_Purkinje () ;
  //IonicModel * ionicmodel_pur = new IonicModel_AP () ;
  IonicModel * ionicmodel_pur = new IonicModel_Passive () ;
  
  int ionicmodel_dof ; 
  ionicmodel_dof =   std::max( ionicmodel_myo->get_n_int_vars(),
			       ionicmodel_pur->get_n_int_vars() );

  ionicmodel_myo -> print_info();
  ionicmodel_pur -> print_info();
  
  //====== Local assembly pointer
  SYS_T::commPrint("===> Initialize local assembly routine ... \n");

  std::vector<IPLocAssem *> locAssem_array; 
  locAssem_array.resize(locElem->get_nlocalele());
  std::vector<IonicModel *> ionicmodel_array; 
  ionicmodel_array.resize(locElem->get_nlocalele());

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee) {

    local_elemType = GMIptr->get_elemType(locElem->get_elem_loc(ee));

    if ( local_elemType == 512 ){
      locAssem_array[ee] =
    	new PLocAssem_EP_3D(tm_galpha_ptr, ionicmodel_pur,
    			    GMIptr->get_nLocBas(locElem->get_elem_loc(ee)),
			    quad_line->get_num_quadPts() );
      ionicmodel_array[ee]= ionicmodel_pur;
    }
    else if ( local_elemType == 501 ){
      locAssem_array[ee] =
	new PLocAssem_EP_3D(tm_galpha_ptr, ionicmodel_myo,
			    GMIptr->get_nLocBas(locElem->get_elem_loc(ee)),
			    quadv->get_num_quadPts() );
      ionicmodel_array[ee]= ionicmodel_myo;
    }
  }

  // ---------------------------------------

  // ======= Solution Initialization =======
  PDNSolution * disp = new PDNSolution_EP(pNode, fNode, locbc,  3); 
  PDNSolution * velo = new PDNSolution_EP(pNode, fNode, locbc,  0);
  PDNSolution * hist = new PDNSolution_EP(pNode, ionicmodel_dof,
					  fNode,  locbc, 0, locIEN, ionicmodel_array);

  //std::cout << "initial solution: " << std::endl;
  //hist->PrintNoGhost();

  //if( is_restart )
  //{
  //  initial_index = restart_index;
  //  initial_time  = restart_time;
  //  initial_step  = restart_step;
  //
  //  SYS_T::file_exist_check(restart_name.c_str());
  //
  //  sol->ReadBinary(restart_name.c_str());
  //  PetscPrintf(PETSC_COMM_WORLD, "===> Read sol from disk as a restart run... \n");
  //  PetscPrintf(PETSC_COMM_WORLD, "     restart_name: %s \n", restart_name.c_str());
  //  PetscPrintf(PETSC_COMM_WORLD, "     restart_time: %e \n", restart_time);
  //  PetscPrintf(PETSC_COMM_WORLD, "     restart_index: %d \n", restart_index);
  //  PetscPrintf(PETSC_COMM_WORLD, "     restart_step: %e \n", restart_step);
  //}
  
  PDNTimeStep * timeinfo = new PDNTimeStep(initial_index, initial_time, initial_step);

  // ============= Global Assembly pointer
  int vpetsc_type = 0; // petsc version controller
  PGAssem_EP * gloAssem_ptr
    = new PGAssem_EP(locAssem_array, GMIptr, locElem, locIEN, pNode, locbc, rank);
  // ============= Estimate the matrix structure
  gloAssem_ptr->Assem_nonzero_estimate( locElem, locAssem_array,
					locIEN, pNode, locbc );
  
  gloAssem_ptr->Fix_nonzero_err_str();
  SYS_T::commPrint("===> Matrix nonzero structure fixed ... \n");

  // 2.5 Setup linear solver context
  PLinear_Solver_PETSc * lsolver = new PLinear_Solver_PETSc(PCJACOBI);
  SYS_T::commPrint("===> PETSc linear solver setted up:\n");
  SYS_T::commPrint("----------------------------------------------------------- \n");
  lsolver->Info();
  SYS_T::commPrint("----------------------------------------------------------- \n");

  //PC preproc; lsolver->GetPC(&preproc);
  //PCSetType( preproc, PCASM );

  // 2.6 Assembly mass matrix and solve for consistent initial solution
  SYS_T::commPrint("Mass residual: \n"); 
  gloAssem_ptr->Clear_KG();
  gloAssem_ptr->Assem_mass_residual( disp, timeinfo,
				     locElem, locAssem_array, locIEN, pNode,
				     fNode, quadArray, elemArray, locbc );
  //gloAssem_ptr->Print_G();
  //SYS_T::commPrint("Mass tangent: \n"); 
  //MatView(gloAssem_ptr->K, PETSC_VIEWER_STDOUT_WORLD);

  lsolver->Solve( gloAssem_ptr->K, gloAssem_ptr->G, velo); 

  SYS_T::commPrint("initial solution's time derivative obtained. \n");
  //SYS_T::commPrint("Velo: \n"); 
  //velo->PrintNoGhost();
  //SYS_T::commPrint("Disp. \n"); 
  //disp->PrintNoGhost();  

  // 2.7 Setup nonlinear solver context
  PNonlinear_Solver_EP * nsolver
    = new PNonlinear_Solver_EP(nl_rtol, nl_atol,
			       nl_dtol, nl_maxits, nl_refreq);
  SYS_T::commPrint("===> Nonlinear solver setted up:\n");
  nsolver->Info();

  // 2.8 Setup time marching context
  PTime_Solver_EP_OperatorSplit * tsolver =
    new PTime_Solver_EP_OperatorSplit( sol_bName, sol_record_freq,
				       ttan_renew_freq, final_time );

  SYS_T::commPrint("===> Time marching solver setted up:\n");
  tsolver->Info();

  SYS_T::commPrint("===> Start Finite Element Analysis:\n");
  tsolver->TM_generalized_alpha(velo, disp, hist, timeinfo, tm_galpha_ptr,
				locElem, locIEN, pNode,
				fNode, locbc, quadArray, elemArray,
				ionicmodel_array, locAssem_array,
				gloAssem_ptr, lsolver, nsolver );

  // ======= PETSc Finalize =======
  SYS_T::commPrint("\n===> Clean memory ... \n");
  delete fNode;
  delete GMIptr;
  delete PartBasic;
  delete locIEN;
  delete locElem;
  delete locbc;
  delete quad_line; delete quadv; delete quads;
  delete pNode;
  delete disp;
  delete velo;
  delete hist;  
  delete timeinfo;
  delete gloAssem_ptr;
  delete tm_galpha_ptr;
  delete ionicmodel_pur;
  delete ionicmodel_myo;
  delete lsolver;
  delete nsolver;
  delete tsolver;
  std::vector<FEAElement *>::iterator it_elema;
  for(it_elema = elemArray.begin(); it_elema != elemArray.end(); ++it_elema) {
    delete *it_elema;
  }
  std::vector<IPLocAssem *>::iterator it_loca;
  for(it_loca = locAssem_array.begin(); it_loca != locAssem_array.end(); ++it_loca) {
    delete *it_loca;
  }
  
  //---------------------------------------
  SYS_T::commPrint("===> Exit program. \n");
  PetscFinalize();
  return 0;
}
//EOF
