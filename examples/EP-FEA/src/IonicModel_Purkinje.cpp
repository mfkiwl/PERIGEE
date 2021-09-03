#include "IonicModel_Purkinje.hpp"

IonicModel_Purkinje::IonicModel_Purkinje()
  //           d_iso,   d_ani, chi,   C_m  n_int_var
  : IonicModel(5.0*4.0, 0.0,   140.0, 0.1, 1),
    ap_1{80.0}, ap_2{80.0}, ap_3{3.0}, m1{0.2},
    m2{0.3}, alpha{0.01}, gamma{0.002}, b{0.15}, c{8.0}
{
  //SYS_T::commPrint("IonicModel_Purkinje constructor. \n");
};

IonicModel_Purkinje::~IonicModel_Purkinje()
{
  //SYS_T::commPrint("IonicModel_Purkinje destructor. \n");
};

void IonicModel_Purkinje::print_info () const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  Purkinje ionic model: \n");
  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_1 = %e \n", ap_1);
  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_2 = %e \n", ap_2);
  PetscPrintf(PETSC_COMM_WORLD, "\t  ap_3 = %e \n", ap_3);
};

void IonicModel_Purkinje::get_Istim(double &Istim,
				    const double &t,
				    const double &x,
				    const double &y,
				    const double &z ) const
{
  //if (( x >= 1.0 ) || ( y >= 1.0 )) {
  //excite 1st node of purkinje network:
  if (( std::sqrt(  std::pow(x-(-106.7), 2.0)
   		    + std::pow(y-(-301.9), 2.0)
   		    + std::pow(z-( 248.2), 2.0)  ) <= 1.0 )
      || ( std::sqrt( std::pow(x-(-104.9), 2.0)
   		      + std::pow(y-(-304.0), 2.0)
   		      + std::pow(z-( 234.0), 2.0)  ) <= 1.0 ) ){
    if(t <= 2.0){
      Istim = -300.0;
    }
  }
  else {
    Istim = 0.0 ;
  }
};

void IonicModel_Purkinje::run_ionic(const std::vector<double> &r_old_in,
				    const double &V_old_in,
				    const double &I_stim,
				    const double &dt_in,
				    std::vector<double> &r_new,
				    double &V_new) const
{
  //non dimensionalize 
  const double V_old { (V_old_in+ap_2)/ap_1};
  const double r_old= r_old_in.at(0);
  const double dt= dt_in/ap_3;

  V_new = V_old
    + (dt/C_m)*(c*V_old*(V_old-alpha)*(1-V_old)- r_old*V_old );
  
  r_new.at(0) = r_old
    + dt*((gamma+(m1*r_old)/(m2+V_old))
	  * (-r_old - c*V_old*(V_old-b-1.0)));

  V_new = V_new * ap_1 - ap_2; 
  V_new = V_new -  I_stim/chi * dt_in/C_m;
}


void IonicModel_Purkinje::get_int_vars(double* val) const
{
  //SYS_T::commPrint("Purkinje ionic,  get int vars. \n"); 
  val[ 0 ]=  0 ;
};

// EOF
