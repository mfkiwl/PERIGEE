#include "IonicModel_Purkinje.hpp"

IonicModel_Purkinje::IonicModel_Purkinje()
  //           d_iso,   d_ani, chi,   C_m  n_int_var
  : IonicModel(5.0*4.0, 0.0,   140.0, 0.1, 1),
    ap_1{100}, ap_2{80}, ap_3{12.9}, m1{0.2},
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
  //excite 1st node of purkinje network:
  if (( std::sqrt(  std::pow(x-(-106.7), 2.0)
		    + std::pow(y-(-301.9), 2.0)
		    + std::pow(z-( 248.2), 2.0)  ) <= 1.0 )
      || ( std::sqrt( std::pow(x-(-104.9), 2.0)
		      + std::pow(y-(-304.0), 2.0)
		      + std::pow(z-( 234.0), 2.0)  ) <= 1.0 ) ){
    //  if ( x<=-1.0 ) {
    if(t <= 1.0){
      Istim = -300.0;
    }
  }
  else {
    Istim = 0.0 ;
  }
};

//double IonicModel_Purkinje::get_diso() const
//{
//  return d_iso;
//};
//
//double IonicModel_Purkinje::get_dani() const
//{
//  return d_ani;
//};
//
//double IonicModel_Purkinje::get_chi() const
//{
//  return chi;
//};
//
//double IonicModel_Purkinje::get_C_m() const
//{
//  return C_m;
//};

void IonicModel_Purkinje::get_Iion(const double &r_old_in,
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

void IonicModel_Purkinje::get_int_vars(double* val) const
{
  SYS_T::commPrint("Purkinje ionic,  get int vars. \n"); 
  val[ 0 ]=  -1 ;
}


// EOF
