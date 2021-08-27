#include "IonicModel_TTP.hpp"

IonicModel_TTP::IonicModel_TTP()
  //           d_iso,     d_ani,     chi,   C_m, n_int_vars
  : IonicModel(0.012*4.0, 0.078*4.0, 140.0, 0.1, 18)
{
  //SYS_T::commPrint("AP constructor. \n");
};

IonicModel_TTP::~IonicModel_TTP()
{
  //SYS_T::commPrint("AP destructor. \n");
};

void IonicModel_TTP::print_info () const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  Ten-Tusscher Panfilov EP: \n");
//  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_1 = %e \n", ap_1);
//  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_2 = %e \n", ap_2);
//  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_3 = %e \n", ap_3);
};

void IonicModel_TTP::get_Iion(const double &r_old_in,
			     const double &V_in,
			     const double &I_stim,
			     double &f_r,
			     double &Iion) const
{

//  //non dimensionalize 
//  const double V_old { (V_in+ap_2)/ap_1};
//  const double r_old {r_old_in};
//  
//  Iion = -(c*V_old*(V_old-alpha)*(1-V_old)- r_old*V_old );
//
//  f_r  = (gamma+(m1*r_old)/(m2+V_old))
//                          * (-r_old - c*V_old*(V_old-b-1.0));
//
//  //redimensionalize to be multiplied with dt
  Iion = 0.0;
  f_r  = 0.0;
}

void IonicModel_TTP::get_int_vars(double* val) const
{
  SYS_T::commPrint("TTP get int vars. \n"); 
  val[ 0 ]=  0 ;
  val[ 1 ]=  1 ;
  val[ 2 ]=  2 ;
  val[ 3 ]=  3 ;
  val[ 4 ]=  4 ;
  val[ 5 ]=  5 ;
  val[ 6 ]=  6 ;
  val[ 7 ]=  7 ;
  val[ 8 ]=  8 ;
  val[ 9 ]=  9 ;
  val[10 ]= 10 ;
  val[11 ]= 11 ;
  val[12 ]= 12 ;
  val[13 ]= 13 ;
  val[14 ]= 14 ;
  val[15 ]= 15 ;
  val[16 ]= 16 ;
  val[17 ]= 17 ;
}

// EOF
