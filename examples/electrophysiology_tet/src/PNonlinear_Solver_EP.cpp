#include "PNonlinear_Solver_EP.hpp"

PNonlinear_Solver_EP::PNonlinear_Solver_EP(
    const double &input_nrtol, const double &input_natol,
    const double &input_ndtol,
    const int &input_max_iteration, const int &input_renew_freq )
: nr_tol(input_nrtol), na_tol(input_natol), nd_tol(input_ndtol),
  nmaxits(input_max_iteration), nrenew_freq(input_renew_freq)
{
}


PNonlinear_Solver_EP::~PNonlinear_Solver_EP()
{}


void PNonlinear_Solver_EP::Info() const
{
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
  PetscPrintf(PETSC_COMM_WORLD, "relative tolerance: %e \n", nr_tol);
  PetscPrintf(PETSC_COMM_WORLD, "absolute tolerance: %e \n", na_tol);
  PetscPrintf(PETSC_COMM_WORLD, "divergence tolerance: %e \n", nd_tol);
  PetscPrintf(PETSC_COMM_WORLD, "maximum iteration: %d \n", nmaxits);
  PetscPrintf(PETSC_COMM_WORLD, "tangent matrix renew frequency: %d \n", nrenew_freq);
  PetscPrintf(PETSC_COMM_WORLD, "----------------------------------------------------------- \n");
}


void PNonlinear_Solver_EP::Gen_alpha_solve(
    const bool &new_tangent_flag,
    const double &curr_time,
    const double &dt,
    const PDNSolution * const &pre_velo,
    const PDNSolution * const &pre_disp,
    //const PDNSolution * const &pre_hist,    
    const TimeMethod_GenAlpha * const &tmga_ptr,
    const ALocal_Elem * const &alelem_ptr,
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &anode_ptr,
    const FEANode * const &feanode_ptr,
    const ALocal_NodalBC * const &bc_part,
    const IQuadPts * const &quad,
    //const AInt_Weight * const &wei_ptr,
    std::vector<FEAElement *> &ele_ptr,
    //const IonicModel * const &ionicmodel_ptr,
    IPLocAssem * const &lassem_ptr,
    PGAssem_EP * const &gassem_ptr,
    PLinear_Solver_PETSc * const &lsolver_ptr,
    PDNSolution * const &velo,
    PDNSolution * const &disp,
    //PDNSolution * const &hist,
    bool &conv_flag,
    int &nl_counter ) const
{
  nl_counter = 0;
  double ksp_its_num;
  double ksp_max_its_num = (double) lsolver_ptr->get_ksp_maxits();
  double ksp_its_check = 0.3 * ksp_max_its_num;
  double residual_norm = 0.0;
  double initial_norm = 0.0;
  double relative_error = 0.0;

  const double gamma   = tmga_ptr->get_gamma();
  const double alpha_m = tmga_ptr->get_alpha_m();
  const double alpha_f = tmga_ptr->get_alpha_f();

  // predictor
  disp->ScaleValue(0.0);
  velo->ScaleValue(0.0);
  // do not need to predict hist_new,
  // that's going to be predicted and calculated in ionic model
  // no need: hist->ScaleValue(0.0);
  
  PDNSolution velo_step(*velo);
  
  disp->PlusAX(*pre_disp, 1.0); //iga book p/205 eq7.40 
  //hist->PlusAX(*pre_hist, 1.0); 
  velo->PlusAX(*pre_velo, (gamma-1.0)/gamma);
  //hist_dot->PlusAX(*pre_hist_dot, (gamma-1.0)/gamma);  

  // define displacement and history at alpha_f and velocity at alpha_m
  PDNSolution disp_alpha(*pre_disp);
  PDNSolution velo_alpha(*pre_velo);
  
  disp_alpha.ScaleValue((1.0 - alpha_f));
  velo_alpha.ScaleValue((1.0 - alpha_m));
  
  disp_alpha.PlusAX(*disp, alpha_f);
  velo_alpha.PlusAX(*velo, alpha_m);

  //  double dt_alp = alpha_f*dt ;

//========================================================
  
// if new tanget flag is true, update the tangent matrix,
// otherwise, keep using the tangent matrix from the previous
// time step
  if( new_tangent_flag )
    {
      gassem_ptr->Clear_KG();

      gassem_ptr->Assem_tangent_residual
	( &velo_alpha, &disp_alpha, //pre_hist, hist, //here hist is hist_alpha
	  curr_time, dt, /// dt_alp, ionicmodel_ptr,
	  alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
	  feanode_ptr, quad, ele_ptr, bc_part );
      
      lsolver_ptr->SetOperator(gassem_ptr->K);
    }
  else
    {
      gassem_ptr->Clear_G();

      gassem_ptr->Assem_residual
	( &velo_alpha, &disp_alpha, //  pre_hist, hist, //here hist is hist_alpha
	  curr_time, dt, //dt_alp, ionicmodel_ptr,
	  alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
	  feanode_ptr, quad, ele_ptr, bc_part );
    }
  
  VecNorm(gassem_ptr->G, NORM_2, &initial_norm);
  PetscPrintf(PETSC_COMM_WORLD, "  Initial residual 2-norm: %e \n", initial_norm);
  
  do
  {
    lsolver_ptr->Solve( gassem_ptr->G, &velo_step );
      
    nl_counter += 1;

    ksp_its_num = (double) lsolver_ptr->get_ksp_it_num();
   
    if(ksp_its_num > ksp_its_check)
    {
      PetscPrintf(PETSC_COMM_WORLD, "Warning: linear solver converges slowly! \n");
      break;
    }

    // corrector
    disp->PlusAX( velo_step, -1.0 * gamma * dt );
    velo->PlusAX( velo_step, -1.0 );

    disp_alpha.PlusAX(velo_step, -1.0 * alpha_f * gamma * dt);
    velo_alpha.PlusAX(velo_step, -1.0 * alpha_m);

//
    if(nl_counter % nrenew_freq == 0)
    {
      gassem_ptr->Clear_KG();
      gassem_ptr->Assem_tangent_residual
	( &velo_alpha, &disp_alpha,// pre_hist, hist, //here hist is hist_alpha
	  curr_time, dt, // dt_alp, ionicmodel_ptr,
	  alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
	  feanode_ptr, quad, ele_ptr, bc_part );
    
      lsolver_ptr->SetOperator(gassem_ptr->K);
    }
    else
    {
      gassem_ptr->Clear_G();
      gassem_ptr->Assem_residual
	( &velo_alpha, &disp_alpha, //  pre_hist, hist, //here hist is hist_alpha
	  curr_time, dt, //dt_alp, ionicmodel_ptr,
	  alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
	  feanode_ptr, quad, ele_ptr, bc_part );

    }

    VecNorm(gassem_ptr->G, NORM_2, &residual_norm);
    PetscPrintf(PETSC_COMM_WORLD, "  --- residual norm: %e \n",
                residual_norm);

    relative_error = residual_norm / initial_norm;

    if( relative_error >= nd_tol )
    {
      PetscPrintf(PETSC_COMM_WORLD, "Warning: nonlinear solver is diverging with error %e \n",
          relative_error);
      break;
    }

  }while(nl_counter < nmaxits && relative_error > nr_tol && 
	 residual_norm > na_tol );

  ////========================================================
  //// calculate hist_new again at t_n+1, at new time step
  ////find new hist variables and ionic currents, and derivatives
  //
  //gassem_ptr->Assem_residual
  //  ( &velo_alpha, &disp_alpha,  pre_hist, hist, // now hist is hist new
  //    curr_time, dt, dt, ionicmodel_ptr, alelem_ptr, lassem_ptr, lien_ptr, anode_ptr,
  //    feanode_ptr, wei_ptr, ele_ptr, bc_part );

  //std::cout<< "hist_new at the end of NLin solver" << std::endl;
  //hist->PrintNoGhost();
  //========================================================
  
  Print_convergence_info(nl_counter, relative_error, residual_norm);

  if(relative_error <= nr_tol || residual_norm <= na_tol)
    conv_flag = true;
  else
    conv_flag = false;

}


// EOF
