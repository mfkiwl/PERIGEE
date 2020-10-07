#ifndef IONICMODEL_HPP
#define IONICMODEL_HPP
// ==================================================================
// IonicModel.hpp
// 
// Interface for ionic models of electroactive materials.
//
// Date: May 25 2020
// Author: Oguz Ziya Tikenogullari
// Contact: o.z.tikenogullari@gmail.com
// ==================================================================
#include "Sys_Tools.hpp"
#include "Tensor4_3D.hpp"
#include "Math_Tools.hpp"

class IonicModel
{
public:
  IonicModel();
  
  virtual ~IonicModel();
  
  virtual void print_info() const; // =0

  virtual void run_model(const double &r_old,
			 const double &dt,
			 const double &Phi,
			 double &f_Phi,
			 double &dP_fP,
			 double &r_new) const; // =0;

  virtual double get_diso() const;

  virtual double get_dani() const;
  
  virtual double get_chi() const;

  virtual double get_C_m() const;

private:
  //conduction parameters
  const double d_iso, d_ani, tol, chi, C_m;
};

#endif
