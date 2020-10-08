#include "IonicModel.hpp"

IonicModel::IonicModel()
{
  SYS_T::commPrint("IonicModel constructor. \n");
}

IonicModel::~IonicModel()
{
  SYS_T::commPrint("IonicModel destructor. \n");
}

void IonicModel::print_info () const
{
  SYS_T::commPrint("IonicModel print info. \n");
}

void IonicModel::run_model(const double &r_old_in,
			   const double &dt_in,
			   const double &Phi_in,
			   double &r_new,
			   double &Phi_new) const
{
  SYS_T::commPrint("IonicModel run_model. \n");
  r_new   = 1.0+r_old_in;
  Phi_new = 1.0+Phi_in;
}

double IonicModel::get_diso() const
{
  SYS_T::commPrint("get_diso not implemented . \n");
  return 0.0;
}

double IonicModel::get_dani() const
{
  SYS_T::commPrint("get_dani not implemented . \n");
  return 0.0;
}

double IonicModel::get_chi() const
{
  SYS_T::commPrint("get_chi not implemented . \n");
  return 0.0;
}

double IonicModel::get_C_m() const
{
  SYS_T::commPrint("get_cm not implemented . \n");
  return 0.0;
}




// EOF
