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
  PetscPrintf(PETSC_COMM_WORLD, "\t  FitzHugh Nagumo EP: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  fh_1 = %e \n", fh_1);
  PetscPrintf(PETSC_COMM_WORLD, "\t  fh_2 = %e \n", fh_2);
  PetscPrintf(PETSC_COMM_WORLD, "\t  fh_3 = %e \n", fh_3);
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
//  double dP_dr=0;
//  //non dimensionalize 
  const double Phi_nd { (Phi_in+fh_2)/fh_1};
  const double dt_nd { dt_in/fh_3};
//
//
  r_new = (dt_nd*(Phi_nd + a)+r_old_in) / (b*dt_nd+1.0) ;
  f_Phi = c*(-std::pow(Phi_nd,3.0) + std::pow(Phi_nd,2.0)*(alpha+1.0)
	     -Phi_nd*alpha - r_new);
  dP_fP = c*(-3.0*std::pow(Phi_nd,2.0) + 2.0*(alpha+1.0)*Phi_nd-alpha)
    -(c*dt_nd)/(b*dt_nd+1);
  //redimensionalize
  f_Phi=fh_1*f_Phi/fh_3;
  dP_fP=dP_fP/fh_3;

//  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  //for manufactured solution:
//  r_new=dt_in;
//  f_Phi=10;
//  dP_fP=0;
  
}


// EOF
