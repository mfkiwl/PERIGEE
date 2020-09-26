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
  
  virtual void material_routine(const double &r_old,
				const double &dt,
				const double &Phi,
				double &f_Phi,
				double &dP_fP,
				double &r_new) const;
private:
  const double fh_1, fh_2, fh_3, alpha,
    a, b, c, d_iso, d_ani, tol, chi, C_m;

};

#endif
