#ifndef IMATERIALMODEL_HPP
#define IMATERIALMODEL_HPP
// ==================================================================
// IMaterialModel.hpp
// 
// Base class for electroactive material models .
// differential eqn : dV/dt = -1/C_m * (I_ion +1/chi*Istim)
// implementation   : V_new = V_old - dt/C_m*Iion
// Iion includes the stimulus current Istim.

// Date: May 25 2020
// Author: Oguz Ziya Tikenogullari, o.z.tikenogullari@gmail.com
// ==================================================================
#include "Sys_Tools.hpp"
#include "Tensor4_3D.hpp"
#include "Math_Tools.hpp"
#include "Vec_Tools.hpp"

class IonicModel
{
public:
  IonicModel(const double d_iso_in,
	     const double d_ani_in,
	     const double chi_in,
	     const double C_m_in,
	     const double n_int_vars_in);
  
  virtual ~IonicModel();
  
  virtual void print_info() const; //const =0 ;

  double get_diso() const;

  double get_dani() const;

  double get_chi() const;

  double get_C_m() const;

  double get_n_int_vars() const;

  virtual void get_int_vars( double* val) const;

  //only one I_stim at t_n
  void Forward_Euler(const std::vector<double> &r_old_in,
		     const double &dt_in,
		     const double &V_in,
		     const std::vector<double> &I_stim,
		     std::vector<double> &r_new,
		     double &V_new) const ;

  //for this function we need 3 values of Istim,
  //at t_n , at  t_{n+1/2} and t_{n+1}
  void Runge_Kutta_4(const std::vector<double> &r_old_in,
		     const double &dt_in,
		     const double &V_in,
		     const std::vector<double> &I_stim,
		     std::vector<double> &r_new,
		     double &V_new) const ;

  virtual void get_Istim(double &Istim,
			 const double &t,
			 const double &x,
			 const double &y,
			 const double &z ) const;

protected:
//  virtual void get_Iion(const std::vector<double> &r_old_in,
//			const double &V_in,
//			const double &I_stim,
//			std::vector<double> &f_r,
//			double &Iion ) const;


  virtual void run_ionic(const std::vector<double> &r_old_in,
			 const double &V_old_in,
			 const double &I_stim,
			 const double &dt_in,
			 std::vector<double> &r_new,
			 double &V_new) const;

  const double d_iso, d_ani, chi, C_m, n_int_vars;

};

#endif
