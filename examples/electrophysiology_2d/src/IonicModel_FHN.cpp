#include "IonicModel_FHN.hpp"

IonicModel_FHN::IonicModel_FHN()
  //           d_iso, d_ani, chi,  C_m
  : IonicModel(0.1, 0.0, 140.0, 1.0),
    fh_1{65.0}, fh_2{35.0}, fh_3{200.0}, 
    alpha{-0.5}, a{0.0}, b{-0.6}, c{50.0}
{
  //SYS_T::commPrint("FHN constructor. \n");
};

IonicModel_FHN::~IonicModel_FHN()
{
  //SYS_T::commPrint("FHN destructor. \n");
};

void IonicModel_FHN::print_info () const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  FitzHugh-Nagumo EP: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  fhn_1 = %e \n", fh_1);
  PetscPrintf(PETSC_COMM_WORLD, "\t  fhn_2 = %e \n", fh_2);
  PetscPrintf(PETSC_COMM_WORLD, "\t  fhn_3 = %e \n", fh_3);
  PetscPrintf(PETSC_COMM_WORLD, "\t  alpha = %e \n", alpha);
  PetscPrintf(PETSC_COMM_WORLD, "\t  a     = %e \n", a);
  PetscPrintf(PETSC_COMM_WORLD, "\t  b     = %e \n", b);
  PetscPrintf(PETSC_COMM_WORLD, "\t  c     = %e \n", c);
};

//implicit solving for FHN model is very cheap 
//I recommend using that instead. 

void IonicModel_FHN::get_Iion(const double &r_old_in,
			      const double &V_in,
			      const double &I_stim,
			      double &f_r,
			      double &Iion) const
{

  //non dimensionalize 
  const double V_old { (V_in+fh_2)/fh_1};
  const double r_old {r_old_in};

  f_r  = V_old + a - b*r_old; 
  
  Iion = -c*(-std::pow(V_old,3.0) + std::pow(V_old,2.0)*(alpha+1.0)
	     -V_old*alpha - r_old);

  //redimensionalize to be multiplied with dt
  // and add Istim/chi to the ionic current
  Iion = Iion*fh_1/fh_3 + I_stim/chi;
  f_r  = f_r/fh_3;
}


// EOF
