#include "IonicModel.hpp"

IonicModel::IonicModel()
  :fh_1{65.0}, fh_2{35.0}, fh_3{200.0}, alpha{-0.5},
   a{0.0}, b{-0.6}, c{10.0}, d_iso{0.1}, d_ani{0.0}, tol{1e-8},
   chi{140}, C_m{1}
{
  SYS_T::commPrint("IonicModel constructor. \n");
};

IonicModel::~IonicModel()
{
  SYS_T::commPrint("IonicModel destructor. \n");
};

void IonicModel::print_info () const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  Test ionic model: \n");
};

double IonicModel::get_diso() const
{
  return d_iso;
};

double IonicModel::get_dani() const
{
  return d_ani;
};

double IonicModel::get_chi() const
{
  return chi;
};

double IonicModel::get_C_m() const
{
  return C_m;
};

void IonicModel::material_routine(const double &r_old_in,
				  const double &dt_in,
				  const double &Phi_in,
				  double &f_Phi,
				  double &dP_fP,
				  double &r_new) const
{
  f_Phi=1.0;
  dP_fP=0.0;
  r_new=0.0;

}

// EOF
