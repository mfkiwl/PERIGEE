#include "IonicModel_AP.hpp"

IonicModel_AP::IonicModel_AP()
  //           d_iso,     d_ani,     chi,   C_m, n_int_vars
  : IonicModel(0.012*1.0, 0.078*1.0, 140.0, 0.1, 1),
    ap_1{100.0}, ap_2{80.0}, ap_3{12.9}, m1{0.2},
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


void IonicModel_AP::run_ionic(const std::vector<double> &r_old_in,
			      const double &V_old_in,
			      const double &I_stim,
			      const double &dt_in,
			      std::vector<double> &r_new,
			      double &V_new) const
{
  //non dimensionalize 
  const double V_old { (V_old_in+ap_2)/ap_1};
  const double r_old= r_old_in.at(0);
  //const double dt= dt_in/ap_3;
  double Iion, fr; 

  Iion = -(c*V_old*(V_old-alpha)*(1-V_old)- r_old*V_old );
  fr   = (gamma+(m1*r_old)/(m2+V_old))
	  * (-r_old - c*V_old*(V_old-b-1.0));

  Iion=Iion*ap_1/ap_3;
  fr  = fr/ap_3;

  V_new       = V_old_in - Iion * dt_in/C_m - I_stim/chi * dt_in/C_m;
  r_new.at(0) = r_old_in.at(0)+ dt_in * fr;
  
  // V_new = V_old
  //   + (dt/C_m)*(c*V_old*(V_old-alpha)*(1-V_old)- r_old*V_old );

  // r_new.at(0) = r_old
  //   + dt*((gamma+(m1*r_old)/(m2+V_old))
  // 	  * (-r_old - c*V_old*(V_old-b-1.0)));

  // V_new = V_new * ap_1 - ap_2; 
  // V_new = V_new -  I_stim/chi * dt_in/C_m;
}


void IonicModel_AP::get_Istim(double &Istim,
			       const double &t,
			       const double &x,
			       const double &y,
			       const double &z ) const
{
  if (( x <= 0.1 ) && ( y <= 0.1 ) && ( z <= 0.1 )) {
    if(t <= 2.0){
      Istim = -10.0;
    }
  }
  else {
    Istim = 0.0 ;
  }
};

void IonicModel_AP::get_int_vars(double* val) const
{
  // initiate the internal variable (recovery variable r) 
  val[ 0 ]=  0 ;
}

// EOF
