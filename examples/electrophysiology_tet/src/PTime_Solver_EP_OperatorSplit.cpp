#include "PTime_Solver_EP_OperatorSplit.hpp"

PTime_Solver_EP_OperatorSplit::PTime_Solver_EP_OperatorSplit(
    const std::string &input_name, const int &input_record_freq,
    const int &input_renew_tang_freq, const double &input_final_time )
: final_time(input_final_time), sol_record_freq(input_record_freq),
  renew_tang_freq(input_renew_tang_freq), pb_name(input_name)
{
}


PTime_Solver_EP_OperatorSplit::~PTime_Solver_EP_OperatorSplit()
{}


std::string PTime_Solver_EP_OperatorSplit::Name_Generator(const int &counter) const
{
  int aux = 900000000 + counter;
  std::ostringstream temp;
  temp<<aux;

  std::string out_name(pb_name);
  out_name.append(temp.str());
  return out_name;
}

void PTime_Solver_EP_OperatorSplit::Info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "final time: %e \n", final_time);
  PetscPrintf(PETSC_COMM_WORLD, "solution record frequency : %d \n", sol_record_freq);
  PetscPrintf(PETSC_COMM_WORLD, "tangent update frequency over time steps: %d \n", renew_tang_freq);
  PetscPrintf(PETSC_COMM_WORLD, "solution base name: %s \n", pb_name.c_str());
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


//generalized alpha time marching with history variables
void PTime_Solver_EP_OperatorSplit::TM_generalized_alpha(
    const PDNSolution * const &init_velo,
    const PDNSolution * const &init_disp,
    const PDNSolution * const &init_hist,    
    PDNTimeStep * const &time_info,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &bc_part,
    //const AInt_Weight * const &wei_ptr,
    const IQuadPts * const &quad,
    std::vector<FEAElement *> &ele_ptr,
    const IonicModel * const &ionicmodel_ptr,
    IPLocAssem * const &lassem_ptr,
    PGAssem_EP * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_Solver_EP * const &nsolver_ptr
    ) const
{
  PDNSolution * pre_disp = new PDNSolution(*init_disp);
  PDNSolution * cur_disp = new PDNSolution(*init_disp);
  PDNSolution * tmp_disp = new PDNSolution(*init_disp);
  PDNSolution * tmp2_disp = new PDNSolution(*init_disp);
  PDNSolution * pre_velo = new PDNSolution(*init_velo);
  PDNSolution * cur_velo = new PDNSolution(*init_velo);
  PDNSolution * pre_hist = new PDNSolution(*init_hist);//history var.
  PDNSolution * cur_hist = new PDNSolution(*init_hist);//history var.
  PDNSolution * tmp_hist = new PDNSolution(*init_hist);//history var.

  double time, dt, dt_half, time_half;

  // save the initial solution
  std::string sol_name = Name_Generator(time_info->get_index());
  std::string hist_sol_name = "hist_" + sol_name;
  cur_disp->WriteBinary(sol_name.c_str());
  cur_hist->WriteBinary(hist_sol_name.c_str()); // does it overwrite?

  bool conv_flag, renew_flag;
  int nl_counter;
  while( time_info->get_time() < final_time )
    {
      if(time_info->get_index() % renew_tang_freq == 0)
	renew_flag = true;
      else
	renew_flag = false;

      dt        = time_info->get_step();
      dt_half   = dt/2.0;
      time      = time_info->get_time();
      time_half =time+dt_half;

     // //Manuf soln no diffusion
     // gassem_ptr->Update_nodal_velo(pre_disp, pre_hist,
     // 			    time, dt,
     // 			    ionicmodel_ptr, alelem_ptr, lien_ptr,
     // 			    anode_ptr, feanode_ptr, ele_ptr, bc_part,
     //             	    cur_disp, cur_hist);

      //Step 1 - Update of ionic variables
      gassem_ptr->Update_nodal_velo(pre_disp, pre_hist,
      			    time, dt_half,
      			    ionicmodel_ptr, alelem_ptr, lien_ptr,
			    anode_ptr, feanode_ptr, bc_part,
      			    tmp_disp, tmp_hist   );
      
      //Step 2 - gen_alpha_solve  for the diffusion problem 
      nsolver_ptr
      	->Gen_alpha_solve(renew_flag, time, dt, 
      			  pre_velo, tmp_disp,
      			  tmga_ptr, alelem_ptr, lien_ptr,
      			  anode_ptr, feanode_ptr, bc_part,
      			  quad, ele_ptr, 
      			  lassem_ptr, gassem_ptr, 
      			  lsolver_ptr, cur_velo, tmp2_disp,
      			  conv_flag, nl_counter);
      
      //Step 3 - Update of ionic variables
      gassem_ptr->Update_nodal_velo(tmp2_disp, tmp_hist,
      	         	    time_half, dt_half,
      			    ionicmodel_ptr, alelem_ptr, lien_ptr,
      			    anode_ptr, feanode_ptr, bc_part,
      			    cur_disp, cur_hist   ); 

      time_info->TimeIncrement();

      PetscPrintf(PETSC_COMM_WORLD, "Time = %e, dt = %e, index = %d \n",
		  time_info->get_time(), time_info->get_step(),
		  time_info->get_index());

      if(time_info->get_index()%sol_record_freq == 0)
	{
	  sol_name = Name_Generator( time_info->get_index() );
	  hist_sol_name = "hist_" + sol_name;
	  cur_disp->WriteBinary(sol_name.c_str());
	  cur_hist->WriteBinary(hist_sol_name.c_str()); // does it overwrite?
	}

      pre_disp->Copy(*cur_disp);
      pre_velo->Copy(*cur_velo);
      pre_hist->Copy(*cur_hist);//update history
    }

  delete pre_disp; delete cur_disp; delete tmp_disp; delete tmp2_disp;
  delete pre_velo; delete cur_velo; 
  delete pre_hist; delete cur_hist; delete tmp_hist;//clear memory alloc
}



// ! This version is for mesh with nonuniform element type
void PTime_Solver_EP_OperatorSplit::TM_generalized_alpha(
    const PDNSolution * const &init_velo,
    const PDNSolution * const &init_disp,
    const PDNSolution * const &init_hist,    
    PDNTimeStep * const &time_info,
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN_Mixed * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &bc_part,
    //const AInt_Weight * const &wei_ptr,
    const std::vector< IQuadPts * > &quad_array,
    std::vector<FEAElement *> &ele_array,
    const std::vector< IonicModel * > &ionicmodel_array,
    std::vector< IPLocAssem * > &lassem_array,
    PGAssem_EP * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PNonlinear_Solver_EP * const &nsolver_ptr
    ) const
{
  PDNSolution * pre_disp = new PDNSolution(*init_disp);
  PDNSolution * cur_disp = new PDNSolution(*init_disp);
  PDNSolution * tmp_disp = new PDNSolution(*init_disp);
  PDNSolution * tmp2_disp= new PDNSolution(*init_disp);
  PDNSolution * pre_velo = new PDNSolution(*init_velo);
  PDNSolution * cur_velo = new PDNSolution(*init_velo);
  PDNSolution * pre_hist = new PDNSolution(*init_hist);//history var.
  PDNSolution * cur_hist = new PDNSolution(*init_hist);//history var.
  PDNSolution * tmp_hist = new PDNSolution(*init_hist);//history var.

  double time, dt, dt_half, time_half;

  // save the initial solution
  std::string sol_name = Name_Generator(time_info->get_index());
  std::string hist_sol_name = "hist_" + sol_name;
  cur_disp->WriteBinary(sol_name.c_str());
  cur_hist->WriteBinary(hist_sol_name.c_str()); // does it overwrite?

  bool conv_flag, renew_flag;
  int nl_counter;
  while( time_info->get_time() < final_time )
    {
      if(time_info->get_index() % renew_tang_freq == 0)
	renew_flag = true;
      else
	renew_flag = false;

      dt        = time_info->get_step();
      dt_half   = dt/2.0;
      time      = time_info->get_time();
      time_half =time+dt_half;

      ////Manuf soln no diffusion
      //gassem_ptr->Update_nodal_velo(pre_disp, pre_hist,
      //				    time, dt,
      //				    ionicmodel_array, alelem_ptr, lien_ptr,
      //				    anode_ptr, feanode_ptr, bc_part,
      //				    cur_disp, cur_hist);
      //
      //Step 1 - Update of ionic variables
      gassem_ptr->Update_nodal_velo(pre_disp, pre_hist,
      			    time, dt_half,
      			    ionicmodel_array, alelem_ptr, lien_ptr,
      			    anode_ptr, feanode_ptr, bc_part,
      			    tmp_disp, tmp_hist );
      
      //Step 2 - gen_alpha_solve  for the diffusion problem 
      nsolver_ptr
      	->Gen_alpha_solve(renew_flag, time, dt, 
      			  pre_velo, tmp_disp,
      			  tmga_ptr, alelem_ptr, lien_ptr,
      			  anode_ptr, feanode_ptr, bc_part,
      			  quad_array, ele_array, 
      			  lassem_array, gassem_ptr, 
      			  lsolver_ptr, cur_velo, tmp2_disp,
      			  conv_flag, nl_counter);
      
      //Step 3 - Update of ionic variables
      gassem_ptr->Update_nodal_velo(tmp2_disp, tmp_hist,
      	         	    time_half, dt_half,
      			    ionicmodel_array, alelem_ptr, lien_ptr,
      			    anode_ptr, feanode_ptr, bc_part,
      			    cur_disp, cur_hist   ); 

      time_info->TimeIncrement();

      PetscPrintf(PETSC_COMM_WORLD, "Time = %e, dt = %e, index = %d \n",
		  time_info->get_time(), time_info->get_step(),
		  time_info->get_index());

      SYS_T::commPrint("Disp. \n"); 
      cur_disp->PrintNoGhost();  

      if(time_info->get_index()%sol_record_freq == 0)
	{
	  sol_name = Name_Generator( time_info->get_index() );
	  hist_sol_name = "hist_" + sol_name;
	  cur_disp->WriteBinary(sol_name.c_str());
	  cur_hist->WriteBinary(hist_sol_name.c_str()); // does it overwrite?
	}

      pre_disp->Copy(*cur_disp);
      pre_velo->Copy(*cur_velo);
      pre_hist->Copy(*cur_hist);//update history
    }

  delete pre_disp; delete cur_disp; delete tmp_disp; delete tmp2_disp;
  delete pre_velo; delete cur_velo; 
  delete pre_hist; delete cur_hist; delete tmp_hist;//clear memory alloc
}
