#include "IonicModel.hpp"

IonicModel::IonicModel()
  :d_iso{0.176}, d_ani{1.158}, tol{1e-8}, chi{1400}, C_m{1.0}
{
  SYS_T::commPrint("IonicModel constructor. \n");
}

IonicModel::~IonicModel()
{
  SYS_T::commPrint("IonicModel destructor. \n");
}

void IonicModel::print_info () const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  Conduction parameters: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  d_iso = %e \n", get_diso());
}

void IonicModel::run_model(const double &r_old_in,
			   const double &dt_in,
			   const double &Phi_in,
			   double &f_Phi,
			   double &dP_fP,
			   double &r_new) const
{
  r_new= 1.0;
  f_Phi= 0.0;
  dP_fP= 0.0;
}

double IonicModel::get_diso() const
{
  return d_iso;
}

double IonicModel::get_dani() const
{
  return d_ani;
}

double IonicModel::get_chi() const
{
  return chi;
}

double IonicModel::get_C_m() const
{
  return C_m;
}




// EOF
