#ifndef IONICMODEL_AP_HPP
#define IONICMODEL_AP_HPP
// ==================================================================
// IonicModel_AP.hpp
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

class IonicModel_AP : public IonicModel
{
public:
  IonicModel_AP();
  
  virtual ~IonicModel_AP();
  
  virtual void print_info() const; //const =0 ;

  //virtual double get_diso() const;
  //
  //virtual double get_dani() const;
  //
  //virtual double get_chi() const;
  //
  //virtual double get_C_m() const;

protected:
  virtual void get_Iion(const double &r_old_in,
			const double &Phi_in,
			const double &I_stim,
			double &f_r,
			double &Iion) const;
private:
  const double ap_1, ap_2, ap_3, m1, m2, alpha, gamma, b, c;

};

#endif
