// ==================================================================
// compute_ECG.cpp
// ------------------------------------------------------------------
// This is a script to calculate ECG signals from simulation results.
// 
// We follow the method given in the paper below. it doesn't account for
// organs etc. and it's a simplistic way of calculating ECG.
//
// Sahli Costabal, Francisco, Jiang Yao, and Ellen Kuhl.
// "Predicting drug‐induced arrhythmias by multiscale modeling."
// International journal for numerical methods in biomedical engineering
// 34.5 (2018): e2964.
//
// Date: Jul 2021
// Author: Oguz Ziya Tikenogullari, Ju Liu
// ==================================================================
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include "Sys_Tools.hpp"
#include "QuadPts_Gauss_Tet.hpp"
#include "QuadPts_Gauss.hpp"
#include "FEANode.hpp"
#include "FEAElement_Tet4.hpp"
#include "FEAElement_Line2_3D_der1.hpp"
#include "APart_Node.hpp"
#include "APart_Basic_Info.hpp"
#include "PostVectSolution.hpp"
#include "Post_tools.hpp"
#include "ALocal_Elem_Fiber.hpp"
#include "ALocal_IEN_Mixed.hpp"
#include "AGlobal_Mesh_Info_Mixed.hpp"

int main( int argc, char * argv[] )
{
  //remove previously existing ECG recording
  char ecg_out_file[]="ecg_recording.csv";

  std::string ecg_file_to_remove {ecg_out_file};
  std::string rm_command ("rm -rf " + ecg_file_to_remove);
  int sysret = system(rm_command.c_str());
  SYS_T::print_fatal_if(sysret != 0, "Error: removing the ecg recording failed. \n");
  
  // Number of quadrature points
  int nqp_line = 1;
  int nqp_tet = 1;

  std::string element_part_file = "epart.h5";
  std::string anode_mapping_file = "node_mapping.h5"; 
  std::string pnode_mapping_file = "post_node_mapping.h5";
  
  const std::string part_file("postpart");
  
  std::string sol_bname("SOL_");
  std::string out_bname = sol_bname;

  const int dof = 1;

  int time_start = 0;
  int time_step = 2;
  int time_end = 200;
  double dt = 0.5;

  bool isXML = true;

  // // electrode location for ECG computation.
  // //3cm away from the left ventricle. 
  //double electrode_x = -152.2731;
  //double electrode_y = -327.1169;  
  //double electrode_z =  248.9918;
  //  //opposite direction
  //double electrode_x = -96.99;
  //double electrode_y = -303.80;  
  //double electrode_z =  248.30;
  ////v2 lead
  //double electrode_x =-81.6568;
  //double electrode_y =-261.7974;  
  //double electrode_z =226.3804;
  //for slab example:
  double v2coor_x = 0.00;
  double v2coor_y = 10.00;
  double v2coor_z = 10.00;
  double v6coor_x = 20.00;
  double v6coor_y = 10.00;
  double v6coor_z = 10.00;

  PetscMPIInt rank, size;
  // ====== PETSc Initialize =====
  PetscInitialize(&argc, &argv, (char *)0, PETSC_NULL);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // ====== Read command line =====
  SYS_T::commPrint("===> Reading arguments from Command line ... \n");
  
  SYS_T::GetOptionInt("-time_start", time_start);
  SYS_T::GetOptionInt("-time_step", time_step);
  SYS_T::GetOptionInt("-time_end", time_end);
  SYS_T::GetOptionReal("-dt", dt);
  SYS_T::GetOptionString("-sol_bname", sol_bname);
  //SYS_T::GetOptionString("-out_bname", out_bname);
  SYS_T::GetOptionBool("-xml", isXML);
  SYS_T::GetOptionReal("-v2coor_x", v2coor_x);
  SYS_T::GetOptionReal("-v2coor_y", v2coor_y);
  SYS_T::GetOptionReal("-v2coor_z", v2coor_z);
  SYS_T::GetOptionReal("-v6coor_x", v6coor_x);
  SYS_T::GetOptionReal("-v6coor_y", v6coor_y);
  SYS_T::GetOptionReal("-v6coor_z", v6coor_z);

  SYS_T::cmdPrint("-sol_bname:", sol_bname);
  //  SYS_T::cmdPrint("-out_bname:", out_bname);
  SYS_T::cmdPrint("-time_start:", time_start);
  SYS_T::cmdPrint("-time_step:", time_step);
  SYS_T::cmdPrint("-time_end:", time_end);
  SYS_T::cmdPrint("-dt:",dt);
  SYS_T::cmdPrint("-v2coor_x:",v2coor_x);
  SYS_T::cmdPrint("-v2coor_y:",v2coor_y);
  SYS_T::cmdPrint("-v2coor_z:",v2coor_z);
  SYS_T::cmdPrint("-v6coor_x:",v6coor_x);
  SYS_T::cmdPrint("-v6coor_y:",v6coor_y);
  SYS_T::cmdPrint("-v6coor_z:",v6coor_z);

  double n_electrodes = 2;
  std::vector< std::vector<double> > electrode_coors(n_electrodes);
  electrode_coors.at(0)={v2coor_x, v2coor_y, v2coor_z};
  electrode_coors.at(1)={v6coor_x, v6coor_y, v6coor_z};  

  if(isXML) PetscPrintf(PETSC_COMM_WORLD, "-xml: true \n");
  else PetscPrintf(PETSC_COMM_WORLD, "-xml: false \n");

  // ===== Read Partition file =====
  SYS_T::commPrint("===> Reading mesh files ... \n");
  FEANode * fNode = new FEANode(part_file, rank);
  ALocal_IEN_Mixed * locIEN = new ALocal_IEN_Mixed(part_file, rank);
  IAGlobal_Mesh_Info * GMIptr = new AGlobal_Mesh_Info_Mixed(part_file,rank);
  APart_Basic_Info * PartBasic = new APart_Basic_Info(part_file, rank);
  ALocal_Elem * locElem = new ALocal_Elem_Fiber(part_file, rank);
  APart_Node * pNode = new APart_Node(part_file, rank);
  SYS_T::commPrint("Done! \n");

  if(size != PartBasic->get_cpu_size())
    SYS_T::print_fatal("Error: number of processors does not match with prepost! \n");

  PetscPrintf(PETSC_COMM_WORLD,
	      "\n===> %d processor(s) are assigned for:", size);
  PetscPrintf(PETSC_COMM_WORLD, "Postprocessing - visualization.\n");

  // ==========build quadrature rules and element container  ===========
  SYS_T::commPrint("\n===> Build quadrature rules and element container ... \n");

  IQuadPts * quad_line   = new QuadPts_Gauss( nqp_line );
  IQuadPts * quadv = new QuadPts_Gauss_Tet( nqp_tet );
  
  std::vector<IQuadPts *> quadArray; 
  quadArray.resize(locElem->get_nlocalele());
  std::vector<FEAElement *> elemArray; 
  elemArray.resize(locElem->get_nlocalele());
  int local_elemType;

  for(int ee=0; ee<locElem->get_nlocalele(); ++ee) {

    local_elemType = GMIptr->get_elemType(locElem->get_elem_loc(ee));

    if ( local_elemType == 512 ){
      quadArray[ee] = quad_line;
      elemArray[ee] = new FEAElement_Line2_3D_der1( nqp_line ); 
    }
    else if ( local_elemType == 501 ){
      quadArray[ee] = quadv;
      elemArray[ee] = new FEAElement_Tet4( nqp_tet ); 
    } else {
      SYS_T::print_fatal("Error: Element type not supported.\n");      
    }
  }

  // ============ start time stepping ============
  int * IEN_e;
  double * ectrl_x; double * ectrl_y; double * ectrl_z;
  double * esol;
  std::vector<double> ecg_el(n_electrodes), ecg_tn_proc(n_electrodes),
    ecg_tn_tot(n_electrodes) ;
  PostVectSolution * pSolu;
  FILE *fp;
  
  //std::ofstream out_str{ ecg_out_file };
  std::ostringstream time_index;
  PetscFOpen(PETSC_COMM_WORLD, ecg_out_file, "w", &fp);
  PetscFPrintf(PETSC_COMM_WORLD, fp, "Time");
  for (int ii=0; ii<n_electrodes; ++ii){
    PetscFPrintf(PETSC_COMM_WORLD, fp, ",\t Lead %d",ii);
  }
  PetscFPrintf(PETSC_COMM_WORLD, fp, " \n");
  
  for(int time = time_start; time<=time_end; time+= time_step)  {

    std::fill(ecg_tn_proc.begin(), ecg_tn_proc.end(), 0.0);
    
    std::string name_to_read(sol_bname);
    std::string name_to_write(out_bname);
    time_index.str("");
    time_index<< 900000000 + time;
    name_to_read.append(time_index.str());
    name_to_write.append(time_index.str());

    //PetscPrintf(PETSC_COMM_WORLD, "Time %d: Read %s  \n",
    //		time, name_to_read.c_str() );
    
    // ===== Manage Solution for postprocessing
    pSolu =
      new PostVectSolution( name_to_read, anode_mapping_file,
			    pnode_mapping_file , pNode, GMIptr, dof );    
    
    // ========= start element loop ==================

    for(int ee=0; ee<locElem->get_nlocalele(); ++ee)  {

      //PetscPrintf(PETSC_COMM_WORLD, " local Element : %d \n", ee);

      std::fill(ecg_el.begin(), ecg_el.end(), 0.0);
      
      int nlocbas_ee= locIEN->get_nLocBas_loc(ee);
      IEN_e   = new int [nlocbas_ee];
      ectrl_x = new double [nlocbas_ee];
      ectrl_y = new double [nlocbas_ee];
      ectrl_z = new double [nlocbas_ee];
      esol    = new double [nlocbas_ee];
      
      locIEN -> get_LIEN_e(ee, IEN_e);
      fNode -> get_ctrlPts_xyz( nlocbas_ee, IEN_e, ectrl_x, ectrl_y, ectrl_z);
      elemArray[ee] -> buildBasis( quadArray[ee], ectrl_x, ectrl_y, ectrl_z );
      
      pSolu -> get_esol(0, nlocbas_ee, IEN_e, esol);

      if ((elemArray[ee]->get_elemDim())  ==3) {
	//only 3 dimensional elements contribute to ecg
	for (int ii=0; ii<n_electrodes; ++ii){
	  ecg_el.at(ii) = POST_T::calculate_ecg(esol, elemArray[ee],
						quadArray[ee], ectrl_x, ectrl_y,
						ectrl_z, electrode_coors.at(ii),
						nlocbas_ee, time); 
	}
      }

      //if(v2_el > 1e-10) {
      //std::cout << "v2_el= " << v2_el <<std::endl;
      //}
      for (int ii=0; ii<n_electrodes; ++ii){
	ecg_tn_proc.at(ii) = ecg_tn_proc.at(ii) + ecg_el.at(ii);
      }
      
      delete [] IEN_e    ; IEN_e    =nullptr;
      delete [] ectrl_x  ; ectrl_x  =nullptr;
      delete [] ectrl_y  ; ectrl_y  =nullptr;
      delete [] ectrl_z  ; ectrl_z  =nullptr;
      delete [] esol     ; esol     =nullptr;
    }

    std::fill(ecg_tn_tot.begin(), ecg_tn_tot.end(), 0.0);

    PetscBarrier(NULL);

    for (int ii=0; ii<n_electrodes; ++ii){
      MPI_Reduce(&(ecg_tn_proc.at(ii)), &(ecg_tn_tot.at(ii)), 1, MPI_DOUBLE,
		 MPI_SUM, 0,  PETSC_COMM_WORLD);
    }
    
    ////write the ecg signal at the current time point.
    //std::cout << "time: " << time*dt << "\t"
    //	      << "v2_tn_proc= " << v2_tn_proc <<std::endl;
    ////out_str << time*dt << "\t" << ecg_tn  << '\n';

    PetscFPrintf(PETSC_COMM_WORLD, fp, "%e", time*dt);
    PetscPrintf(PETSC_COMM_WORLD, "Time: %e", time*dt);
    for (int ii=0; ii<n_electrodes; ++ii){
      PetscFPrintf(PETSC_COMM_WORLD, fp, ", \t %e", ecg_tn_tot.at(ii));
      PetscPrintf(PETSC_COMM_WORLD, " \t Lead%d: %e", ii, ecg_tn_tot.at(ii));
    }
    PetscFPrintf(PETSC_COMM_WORLD, fp, " \n");
    PetscPrintf(PETSC_COMM_WORLD, " \n");	
    delete pSolu;
  }
  
  PetscFClose(PETSC_COMM_WORLD, fp);

  // ===== PETSc Finalize =====
  delete pNode;
  delete GMIptr;
  delete fNode;
  delete locIEN;
  delete locElem;
  delete PartBasic;
  delete quad_line; delete quadv; 



  std::vector<FEAElement *>::iterator it_elema;
  for(it_elema = elemArray.begin(); it_elema != elemArray.end(); ++it_elema) {
    delete *it_elema;
  }

  PetscFinalize();
  return 0;
}

// EOF
