#ifndef IONICMODEL_FHN_HPP
#define IONICMODEL_FHN_HPP
// ==================================================================
// IonicModel_FHN.hpp
// 
// FitzHugh-Nagumo model from Goktepe&Kuhl,2009
//
// Date: May 25 2020
// Author: Oguz Ziya Tikenogullari
// Contact: o.z.tikenogullari@gmail.com
// ==================================================================
#include "Sys_Tools.hpp"
#include "Tensor4_3D.hpp"
#include "Math_Tools.hpp"
#include "IonicModel.hpp"

class IonicModel_FHN : public IonicModel
{
public:
  IonicModel_FHN();
  
  virtual ~IonicModel_FHN();
  
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
  const double fh_1, fh_2, fh_3, alpha, a, b, c;

};

#endif
