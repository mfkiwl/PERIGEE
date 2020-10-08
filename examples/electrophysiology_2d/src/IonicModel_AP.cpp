#include "IonicModel_AP.hpp"

IonicModel_AP::IonicModel_AP()
  :ap_1{100}, ap_2{80}, ap_3{12.9}, m1{0.2}, m2{0.3},alpha{0.01},
   gamma{0.002}, b{0.15}, c{8.0}, d_iso{0.176}, d_ani{1.158}, tol{1e-8},
   chi{1400}, C_m{1.0}
{
  SYS_T::commPrint("AP constructor. \n");
};

IonicModel_AP::~IonicModel_AP()
{
  SYS_T::commPrint("AP destructor. \n");
};

void IonicModel_AP::print_info () const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  Aliev Panfilov EP: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_1 = %e \n", ap_1);
  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_2 = %e \n", ap_2);
  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_3 = %e \n", ap_3);
};

double IonicModel_AP::get_diso() const
{
  return d_iso;
};

double IonicModel_AP::get_dani() const
{
  return d_ani;
};

double IonicModel_AP::get_chi() const
{
  return chi;
};

double IonicModel_AP::get_C_m() const
{
  return C_m;
};

void IonicModel_AP::run_model(const double &r_old_in,
			   const double &dt_in,
			   const double &Phi_in,
			   double &r_new,
			   double &Phi_new) const
{

  //non dimensionalize 
  double Phi_old { (Phi_in+ap_2)/ap_1};
  double r_old {r_old_in};
  const double dt { dt_in/ap_3};
  int   time_steps= 10;
  double  dt_i  = dt/time_steps;

  for(int i=0; i<time_steps; ++i)
    {
      Phi_new = Phi_old
	         + dt_i * (c*Phi_old*(Phi_old-alpha)*(1-Phi_old)
			   - r_old*Phi_old );

      r_new   = r_old + dt_i * (gamma+(m1*r_old)/(m2+Phi_old))
                          	* (-r_old - c*Phi_old*(Phi_old-b-1.0));

      Phi_old=Phi_new; 
      r_old=r_new;
      
    }

  //redimensionalize
  Phi_new = Phi_new*ap_1-ap_2;
}


// EOF
