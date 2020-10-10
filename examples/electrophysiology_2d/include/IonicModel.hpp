#ifndef IMATERIALMODEL_HPP
#define IMATERIALMODEL_HPP
// ==================================================================
// IMaterialModel.hpp
// 
// Base class for electroactive material models .
// differential eqn : dV/dt = -1/C_m * I_ion
// implementation   : V_new = V_old - dt/C_m

// Date: May 25 2020
// Author: Oguz Ziya Tikenogullari, Ju Liu
// Contact: o.z.tikenogullari@gmail.com, liujuy@gmail.com
// ==================================================================
#include "Sys_Tools.hpp"
#include "Tensor4_3D.hpp"
#include "Math_Tools.hpp"

class IonicModel
{
public:
  IonicModel(const double d_iso_in,
	     const double d_ani_in,
	     const double chi_in,
	     const double C_m_in );
  
  virtual ~IonicModel();
  
  virtual void print_info() const; //const =0 ;

  double get_diso() const;

  double get_dani() const;

  double get_chi() const;

  double get_C_m() const;
  
  void Forward_Euler(const double &r_old_in,
			      const double &dt_in,
			      const double &V_in,
			      double &r_new,
			      double &V_new) const ;

  void Runge_Kutta_4(const double &r_old_in,
		     const double &dt_in,
		     const double &V_in,
		     double &r_new,
		     double &V_new) const ;

protected:
  virtual void get_Iion(const double &r_old_in,
			const double &V_in,
			double &f_r,
			double &Iion ) const;

private:
  const double d_iso, d_ani, chi, C_m;

};

#endif
