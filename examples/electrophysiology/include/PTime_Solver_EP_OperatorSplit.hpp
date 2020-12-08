#ifndef PTIME_SOLVER_EP_OPERATORSPLIT_HPP
#define PTIME_SOLVER_EP_OPERATORSPLIT_HPP
// ==================================================================
// PTime_Solver_EP_OperatorSplit.hpp
// 
// Parallel Time marching solver. We implement
// EP solver with operator splitting in this class.
//
// refer to:
// Mathematical Cardiac Electrophysiology book
// Qu, Garfinkel 1999 and Sundnes et al 2005  papers 
//
// Date: Oct 4 2020
// Author: Ju Liu, liujuy@gmail.com
// Modified: Oguz Ziya Tikenogullari, o.z.tikenogullari@gmail.com
// ==================================================================
#include "PDNTimeStep.hpp"
#include "PNonlinear_Solver_EP.hpp"
#include "IonicModel.hpp"

class PTime_Solver_EP_OperatorSplit
{
  public:
    PTime_Solver_EP_OperatorSplit(const std::string &input_name, const int &input_record_freq,
        const int &input_renew_tang_freq, const double &input_final_time );

    ~PTime_Solver_EP_OperatorSplit();

    // ! print the key parameters of the time solver on screen
    void Info() const;

   // ! Perform time marching with history variables
    void TM_generalized_alpha(
        const PDNSolution * const &init_velo,
        const PDNSolution * const &init_disp,
	const PDNSolution * const &init_hist,
        PDNTimeStep * const &time_info,
        const TimeMethod_GenAlpha * const &tmga_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &anode_ptr,
        const FEANode * const &feanode_ptr,
        const IALocal_BC * const &bc_part,
        const AInt_Weight * const &wei_ptr,
        const std::vector<FEAElement *> &ele_ptr,
	const IonicModel * const &ionicmodel_ptr,	
        IPLocAssem * const &lassem_ptr,
	PGAssem_EP * const &gassem_ptr,
        PLinear_Solver_PETSc * const &lsolver_ptr,
        PNonlinear_Solver_EP * const &nsolver_ptr
        ) const;
    
    
  private:
    const double final_time;
    const int sol_record_freq; // the frequency for writing solutions
    const int renew_tang_freq; // the frequency for renewing tangents
    const std::string pb_name; // the problem base name for the solution

    std::string Name_Generator( const int &counter ) const;
};

#endif
