#ifndef IONICMODEL_TTP_HPP
#define IONICMODEL_TTP_HPP
// ==================================================================
// IonicModel_TTP.hpp
// 
// TenTusscher-Panfilov model
//
// this model has 7 ionic currents and 12 gating variables. one of the
// ionic currents is V and it's a field varible in FE analysis.
// Therefore we have 18 internal variables in total (6+12)
//
// Date: Aug 2021
// Author: Oguz Ziya Tikenogullari
// Contact: o.z.tikenogullari@gmail.com
// ==================================================================
#include "Sys_Tools.hpp"
#include "Tensor4_3D.hpp"
#include "Math_Tools.hpp"
#include "IonicModel.hpp"

class IonicModel_TTP : public IonicModel
{
public:
  IonicModel_TTP();
  
  virtual ~IonicModel_TTP();
  
  virtual void print_info() const; //const =0 ;

  //virtual double get_diso() const;
  //
  //virtual double get_dani() const;
  //
  //virtual double get_chi() const;
  //
  //virtual double get_C_m() const;

  virtual void get_int_vars( double* val) const;

protected:
  virtual void get_Iion(const double &r_old_in,
			const double &Phi_in,
			const double &I_stim,
			double &f_r,
			double &Iion) const;
private:
  //const double a;

};

#endif
