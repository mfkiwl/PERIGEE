#include "IonicModel_AP.hpp"

IonicModel_AP::IonicModel_AP()
  //           d_iso, d_ani, chi,  C_m
  : IonicModel(0.0176, 0.0, 140.0, 0.1),
    ap_1{100}, ap_2{80}, ap_3{12.9}, m1{0.2},
    m2{0.3}, alpha{0.01}, gamma{0.002}, b{0.15}, c{8.0}
{
  //SYS_T::commPrint("AP constructor. \n");
};

IonicModel_AP::~IonicModel_AP()
{
  //SYS_T::commPrint("AP destructor. \n");
};

void IonicModel_AP::print_info () const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  Aliev-Panfilov EP: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_1 = %e \n", ap_1);
  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_2 = %e \n", ap_2);
  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_3 = %e \n", ap_3);
};

//double IonicModel_AP::get_diso() const
//{
//  return d_iso;
//};
//
//double IonicModel_AP::get_dani() const
//{
//  return d_ani;
//};
//
//double IonicModel_AP::get_chi() const
//{
//  return chi;
//};
//
//double IonicModel_AP::get_C_m() const
//{
//  return C_m;
//};

void IonicModel_AP::get_Iion(const double &r_old_in,
			     const double &V_in,
			     const double &I_stim,
			     double &f_r,
			     double &Iion) const
{

  //non dimensionalize 
  const double V_old { (V_in+ap_2)/ap_1};
  const double r_old {r_old_in};
  
  Iion = -(c*V_old*(V_old-alpha)*(1-V_old)- r_old*V_old );

  f_r  = (gamma+(m1*r_old)/(m2+V_old))
                          * (-r_old - c*V_old*(V_old-b-1.0));

  //redimensionalize to be multiplied with dt
  Iion = Iion*ap_1/ap_3 + I_stim/chi;;
  f_r  = f_r/ap_3;
}


// EOF
