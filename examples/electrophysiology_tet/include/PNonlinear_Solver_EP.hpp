#ifndef PNONLINEAR_SOLVER_EP_HPP
#define PNONLINEAR_SOLVER_EP_HPP
// ==================================================================
// PNonlinear_Solver_EP.hpp
//
// This class implements the nonlinear solver procedure.
//
// Date: Oct 4 2020
// Author: Ju Liu, liujuy@gmail.com
// Modified: Oguz Ziya Tikenogullari, o.z.tikenogullari@gmail.com
// ==================================================================
#include "TimeMethod_GenAlpha.hpp"
#include "ALocal_NodalBC.hpp"
#include "IPLocAssem.hpp"
#include "PGAssem_EP.hpp"
#include "IPGAssem.hpp"
#include "PLinear_Solver_PETSc.hpp"
#include "Matrix_PETSc.hpp"
//#include "IonicModel.hpp"

class PNonlinear_Solver_EP
{
  public:
    PNonlinear_Solver_EP(
        const double &input_nrtol, const double &input_natol,
        const double &input_ndtol,
        const int &input_max_iteration, const int &input_renew_freq );

    ~PNonlinear_Solver_EP();

    int get_non_max_its() const {return nmaxits;}

    void Info() const;

    //! Generalized-alpha nonlinear solver
    //  which is used in each time step of generalized-alpha time 
    //  scheme
    //  -------------------------------------------------------------
    //  \para bool new_tangent_flag: true->assembly tangent matrix
    //  \para double curr_time: t
    //  \para double dt: dt
    //  \para PDNSolution pre_disp: d_n
    //  \para PDNSolution pre_velo: v_n
    //  \para TimeMethod_GenAlpha tmga_ptr: pointer to time method
    //  \para ALocal_Elem alelem_ptr: local element info
    //  \para ALocal_IEN lien_ptr: IEN arrays
    //  \para APart_Node anode_ptr: local_to_global info
    //  \para FEANode feanode_ptr: control points
    //  \para ALocal_NodalBC bc_part: boundary conditions
    //  \para AInt_Weight wei_ptr: quadrature weights
    //  \para vector<FEAElement*>: basis function at quadrature pts
    //  \para IPLocAssem lassem_ptr: interface for local assembly 
    //  \para PGAssem_NLHeat_GenAlpha gassem_ptr: global assembly method
    //  \para PLinear_Solver_PETSc: linear solver method
    //  \para output PDNSolution disp: d_n+1
    //  \para output PDNSolution velo: v_n+1
    //  \para output bool conv_flag: true if nl iteration converged
    //  \para output int non_ite_counter: number of nl iterations 
    //  -------------------------------------------------------------

    //generalized alpha solver when history variables exist.
    void Gen_alpha_solve(
        const bool &new_tangent_flag,
        const double &curr_time,
        const double &dt,
        const PDNSolution * const &pre_velo,
        const PDNSolution * const &pre_disp,
	//const PDNSolution * const &pre_hist_dot,
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
	//PDNSolution * const &hist_dot,	
	//PDNSolution * const &hist,
        bool &conv_flag,
        int &nl_counter ) const;

  //this version for mixed mesh 
  void Gen_alpha_solve(const bool &new_tangent_flag,
		       const double &curr_time,
		       const double &dt,
		       const PDNSolution * const &pre_velo,
		       const PDNSolution * const &pre_disp,
		       //const PDNSolution * const &pre_hist,    
		       const TimeMethod_GenAlpha * const &tmga_ptr,
		       const ALocal_Elem * const &alelem_ptr,
		       const ALocal_IEN_Mixed * const &lien_ptr,
		       const APart_Node * const &anode_ptr,
		       const FEANode * const &feanode_ptr,
		       const ALocal_NodalBC * const &bc_part,
		       const std::vector< IQuadPts * > &quad_array,
		       ///const IQuadPts * const &quad,
		       //const AInt_Weight * const &wei_ptr,
		       std::vector<FEAElement *> &ele_array,
		       //const IonicModel * const &ionicmodel_ptr,
		       std::vector< IPLocAssem * > &lassem_array,
		       //IPLocAssem * const &lassem_ptr,
		       PGAssem_EP * const &gassem_ptr,
		       PLinear_Solver_PETSc * const &lsolver_ptr,
		       PDNSolution * const &velo,
		       PDNSolution * const &disp,
		       //PDNSolution * const &hist,
		       bool &conv_flag,
		       int &nl_counter ) const;

  private:
    const double nr_tol;
    const double na_tol;
    const double nd_tol;
    const int nmaxits;
    const int nrenew_freq;

    void Print_convergence_info( const int &count, const double rel_err,
        const double abs_err ) const
    {PetscPrintf(PETSC_COMM_WORLD,
        "  === NR ite: %d, r_error: %e, a_error: %e \n", 
        count, rel_err, abs_err);}

};

#endif
