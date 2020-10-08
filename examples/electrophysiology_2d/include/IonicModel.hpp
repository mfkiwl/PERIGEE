#ifndef IMATERIALMODEL_HPP
#define IMATERIALMODEL_HPP
// ==================================================================
// IMaterialModel.hpp
// 
// Interface for electroactive material models .
// Aliev-Panfilov model from Goktepe&Kuhl,2009
//
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
  IonicModel();
  
  virtual ~IonicModel();
  
  virtual void print_info() const; //const =0 ;

  virtual double get_diso() const;

  virtual double get_dani() const;

  virtual double get_chi() const;

  virtual double get_C_m() const;
  
  virtual void run_model(const double &r_old_in,
			 const double &dt_in,
			 const double &Phi_in,
			 double &r_new,
			 double &Phi_new      ) const;
private:
  //  const double ap_1, ap_2, ap_3, m1, m2, alpha, gamma,
  //b, c, d_iso, d_ani, tol, chi, C_m;

};

#endif
