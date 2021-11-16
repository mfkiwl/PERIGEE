#ifndef IONICMODEL_PURKINJE_HPP
#define IONICMODEL_PURKINJE_HPP
// ==================================================================
// IonicModel_Purkinje.hpp
// 
// Aliev-Panfilov model from Goktepe&Kuhl,2009
//
// Date: May 25 2020
// Author: Oguz Ziya Tikenogullari
// Contact: o.z.tikenogullari@gmail.com
// ==================================================================
#include "Sys_Tools.hpp"
#include "Tensor4_3D.hpp"
#include "Math_Tools.hpp"
#include "IonicModel.hpp"

class IonicModel_Purkinje : public IonicModel
{
public:
  IonicModel_Purkinje(const double &cond_scale,
		      const double &LV_delay_in,
		      const double &RV_delay_in);

  IonicModel_Purkinje(const double &cond_scale);
  
  virtual ~IonicModel_Purkinje();
  
  virtual void print_info() const; //const =0 ;

  virtual void get_Istim(double &Istim,
			 const double &t,
			 const double &x,
			 const double &y,
			 const double &z ) const;

  virtual void get_int_vars( double* val) const;

protected:

  virtual void run_ionic(const std::vector<double> &r_old_in,
			 const double &V_old_in,
			 const double &I_stim,
			 const double &dt_in,
			 std::vector<double> &r_new,
			 double &V_new) const;

private:
  const double ap_1, ap_2, ap_3, m1, m2, alpha, gamma, b, c,
    LV_pur_delay, RV_pur_delay; 
};

#endif
