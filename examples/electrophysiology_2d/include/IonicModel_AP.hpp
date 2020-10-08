#ifndef IONICMODEL_ALIEVPANFILOV_HPP
#define IONICMODEL_ALIEVPANFILOV_HPP
// ==================================================================
// IonicModel_AP.hpp
// 
// Interface for electroactive material models .
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

class IonicModel_AP : public IonicModel
{
public:
  IonicModel_AP();
  
  virtual ~IonicModel_AP();
  
  virtual void print_info() const; //const =0 ;

  virtual double get_diso() const;

  virtual double get_dani() const;

  virtual double get_chi() const;

  virtual double get_C_m() const;
  
  virtual void run_model(const double &r_old_in,
			 const double &dt_in,
			 const double &Phi_in,
			 double &r_new,
			 double &Phi_new       ) const;
private:
  const double ap_1, ap_2, ap_3, m1, m2, alpha, gamma,
    b, c, d_iso, d_ani, tol, chi, C_m;

};

#endif
