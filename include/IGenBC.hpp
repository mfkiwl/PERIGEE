#ifndef IGENBC_HPP
#define IGENBC_HPP
// ==================================================================
// IGenBC.hpp
//
// Interface for GenBC classes
//
// This is a pure virtual class for General BC for outflow boundary
// conditions.
//
// Author: Ju Liu
// Date: Aug 21 2019
// Reference: J. Liu, et al. The nested block preconditioning technique
// for the incompressible Navier-Stokes equations with emphasis on
// hemodynamic simulations.
// ==================================================================
#include "Sys_Tools.hpp"

class IGenBC
{
  public:
    bool UserLPM_Dirichlet_flag=false;
    IGenBC(){};

    virtual ~IGenBC(){};

    virtual void print_info() const = 0;

    // --------------------------------------------------------------
    // return the number of faces with elemental boundary condtiions.
    // This value is assumed to match the num_ebc in ALocal_EBC.
    // --------------------------------------------------------------
    virtual int get_num_ebc() const = 0;

    // --------------------------------------------------------------
    // Get the dP/dQ for surface ii
    // for implicit BC's, this value is defined via a difference quotient
    // that is, (get_P(Q+epsilon) - get_P(Q)) / epsilon
    // for simple models, e.g. resistance type bc, this value is just
    // the resistance value on this bc.
    // --------------------------------------------------------------
    virtual double get_m( const int &ii, const double &dot_Q,
       const double &Q )

    {
      SYS_T::print_fatal("Error: IGenBC::get_m is not implemented.\n");
      return 0;
    }

    // --------------------------------------------------------------
    // Get the dP/d(dot_Q) for surface ii
    // for implicit BC's, this value is defined via a difference quotient
    // that is, (get_P(dot_Q+epsilon) - get_P(dot_Q)) / epsilon
    // for simple models, e.g. inductance type bc, this value is just
    // the inductance value on this bc.
    // --------------------------------------------------------------
    virtual double get_n( const int &ii, const double &dot_Q,
       const double &Q )
    {
      SYS_T::print_fatal("Error: IGenBC::get_n is not implemented.\n");
      return 0;
    }

    // --------------------------------------------------------------
    // Get the P value for surface ii, the traction on the surface is
    // modeled as h = P I
    // for resistance bc, for example, this value is
    // Resistance x Q + P_offset
    // --------------------------------------------------------------
    virtual double get_P( const int &ii, const double &dot_Q,
       const double &Q )
    {
      SYS_T::print_fatal("Error: IGenBC::get_P is not implemented.\n");
      return 0;
    }

    // --------------------------------------------------------------
    // Return the pressure at the time step n, which is used as the
    // initial value for ODE integration.
    // For Resistance bc, it is Resistance x Q_previous + P_offset
    // For RCR, it is get_P(ii, Q_previous), which is also stored
    // as Pi_0 + Q_0 x Rp
    // --------------------------------------------------------------
    virtual double get_P0( const int &ii ) const
    {
      SYS_T::print_fatal("Error: IGenBC::get_P0 is not implemented.\n");
      return 0;
    }

    // --------------------------------------------------------------
    // Record solution values as initial conditions for the next time step
    // para ii : the outlet face id, ranging from 0 to num_ebc - 1
    // in_Q_0 : the initial value for the ODE integration of flow rate
    // in_P_0 : the initial value for the ODE integration of averaged Pressure
    // For problems like RCR, it is often convenient to integrate the
    // ODE for Pi, the pressure over the capacitor. So, there can be
    // more data to be initialized. Check the details of each class
    // implementation.
    // --------------------------------------------------------------
    virtual void reset_initial_sol( const int &ii, const double &in_Q_0,
        const double &in_P_0, const double &curr_time )
    {
      SYS_T::print_fatal("Error: IGenBC::reset_initial_sol is not implemented.\n");
    }
    virtual void get_m( double * const &dot_Q, double * const &Q, double * const &m) const
    {
    // --------------------------------------------------------------
    // Get the dP/d(dot_Q) for all outlet surfaces
    // for implicit BC's, this value is defined via a difference quotient
    // that is, (get_P(dot_Q+epsilon) - get_P(dot_Q)) / epsilon
    // for simple models, e.g. inductance type bc, this value is just
    // the inductance value on this bc.
    // --------------------------------------------------------------
      SYS_T::print_fatal("Error: IGenBC::get_m is not implemented.\n");
    }

    virtual void get_n( double * const &dot_Q, double * const &Q, double * const &m ) const
    {
    // --------------------------------------------------------------
    // Get the P value for all outlet surfaces, the traction on the surface is
    // modeled as h = P I
    // for resistance bc, for example, this value is
    // Resistance x Q + P_offset
    // --------------------------------------------------------------

    SYS_T::print_fatal("Error: IGenBC::get_n is not implemented.\n");
    }
    virtual void get_P(double * const &dot_Q, double * const &Q, double * const &P)const
    {
    // --------------------------------------------------------------
    // Return the outlet pressures at the time step n, which is used as the
    // initial value for ODE integration.
    // For Resistance bc, it is Resistance x Q_previous + P_offset
    SYS_T::print_fatal("Error: IGenBC::get_P is not implemented.\n");
    }
    virtual void get_P0( double * const &P0 )const
    {

    // --------------------------------------------------------------
    // Record solution values as initial conditions for the next time step
    // para ii : the outlet face id, ranging from 0 to num_ebc - 1
    // in_Q_0 : the initial value for the ODE integration of flow rate
    // in_P_0 : the initial value for the ODE integration of averaged Pressure
    // For problems like RCR, it is often convenient to integrate the
    // ODE for Pi, the pressure over the capacitor. So, there can be
    // more data to be initialized. Check the details of each class
    // implementation.
    // --------------------------------------------------------------
    SYS_T::print_fatal("Error: IGenBC::get_P0 is not implemented.\n");
    }
    virtual void reset_initial_sol( double * const &in_Q_0,
       double * const &in_P_0, const double &curr_time )
    {
      SYS_T::print_fatal("Error: IGenBC::reset_initial_sol is not implemented.\n");
    }

    virtual void reset_initial_sol( double * const &in_Q_0_Neumann,
       double * const &in_P_0_Neumann, double * const &in_Q_0_Dirichlet,
       double * const &in_P_0_Dirichlet, const double &curr_time )
    {
      SYS_T::print_fatal("Error: IGenBC::reset_initial_sol is not implemented.\n");
    }


    virtual void get_Q0(double * const &Qn) const
    {
      SYS_T::print_fatal("Error: IGenBC::get_Q0 is not implemented.\n");
    }
    virtual double get_Q0( const int &ii ) const
    {
      SYS_T::print_fatal("Error: IGenBC::get_Q0 is not implemented.\n");
      return 0;
    }

    virtual void get_m( double * const &in_dot_Q, double * const &in_Q, double * const &in_P,
      double * const &m ) const
    {
      SYS_T::print_fatal("Error: IGenBC::get_m is not implemented.\n");
    }

    virtual void get_P_Q( double * const &in_dot_Q,
       double * const &in_Q, double * const &in_P, double * const &P_Neumann, double * const &Q_Dirichlet)const
    {
      SYS_T::print_fatal("Error: IGenBC::get_GenBC_P_Q is not implemented.\n");
    }

    virtual double get_curr_P(const int &ii )const
    {
     SYS_T::print_fatal("Error: IGenBC::get_curr_P is not implemented.\n");
     return 0;
    }


    virtual double get_curr_Q(const int &ii)const
    {
     SYS_T::print_fatal("Error: IGenBC::get_curr_Q is not implemented.\n");
     return 0;
    }

    virtual double get_curr_m(const int &ii)const
    {
     SYS_T::print_fatal("Error: IGenBC::get_curr_m is not implemented.\n");
     return 0;
    }

    virtual double get_curr_n(const int &ii)const
    {
     SYS_T::print_fatal("Error: IGenBC::get_curr_n is not implemented.\n");
     return 0;
    }

    virtual void set_curr_P(const int & ii, const double & Pi){
     SYS_T::print_fatal("Error: IGenBC::set_curr_P is not implemented.\n");
    }
    virtual void set_curr_Q(const int & ii, const double & Qi){
      SYS_T::print_fatal("Error: IGenBC::set_curr_Q is not implemented.\n");
    }
    virtual void set_curr_m(const int & ii,const double & mi){
      SYS_T::print_fatal("Error: IGenBC::set_curr_m is not implemented.\n");
    }
    virtual void set_curr_n(const int &ii, const double &ni){
      SYS_T::print_fatal("Error: IGenBC::set_curr_n is not implemented.\n");
    }
    virtual int get_num_Dirichlet_faces()const
    {
      return 0;
    }

};

#endif
