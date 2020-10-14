#include "IonicModel.hpp"

IonicModel::IonicModel(const double d_iso_in,
		       const double d_ani_in,
		       const double chi_in,
		       const double C_m_in )
  : d_iso(d_iso_in), d_ani(d_ani_in), chi(chi_in), C_m(C_m_in)
{
  //SYS_T::commPrint("IonicModel constructor. \n");
}


IonicModel::~IonicModel()
{
  //SYS_T::commPrint("IonicModel destructor. \n");
}


void IonicModel::print_info () const
{
  SYS_T::commPrint("IonicModel print info. \n");
}


void IonicModel::get_Iion(const double &r_old_in,
			  const double &V_in,
			  const double &I_stim,
			  double &f_r,
			  double &Iion) const
{
  SYS_T::commPrint("IonicModel get_Iion. \n");
}

void IonicModel::get_Istim(double &Istim,
			   const double &time,
			   const double &ctrl_x,
			   const double &ctrl_y,
			   const double &ctrl_z ) const
{
  Istim = 1.0;
}


void IonicModel::Forward_Euler(const double &r_old_in,
			       const double &dt_in,
			       const double &V_in,
			       const std::vector<double> &I_stim,
			       double &r_new,
			       double &V_new) const
{
  double V_old = V_in;
  double r_old   = r_old_in;
  int time_steps = 10;
  double dt_i    = dt_in / time_steps;
  double Iion, f_r;

  //std::cout << "begin forw euelr" <<std::endl; 
  for(int i=0; i<time_steps; ++i)
    {
      //I stim is sent to ionic model and included
      // in Iion and state variable evolution.
      get_Iion( r_old, V_old, I_stim.at(0), f_r, Iion);
      
      std::cout << "Istim[0]: " << I_stim.at(0)  <<std::endl;
      std::cout << "Iion    : " << Iion        <<std::endl;

      V_new   = V_old - dt_i/C_m * Iion;
      r_new   = r_old   + dt_i * f_r;

      //std::cout << "V_new : " << V_new <<std::endl;
      //std::cout << "r_new : " << r_new <<std::endl; 

      V_old=V_new; 
      r_old=r_new;
    }
}


void IonicModel::Runge_Kutta_4(const double &r_old_in,
			       const double &dt_in,
			       const double &V_in,
			       const std::vector<double> &I_stim,
			       double &r_new,
			       double &V_new) const
{
  double V_old = V_in;
  double r_old   = r_old_in;
  int time_steps = 10;
  double dt_i    = dt_in / time_steps;
  double Iion, f_r, K1r, K2r, K3r, K4r, K1V, K2V, K3V, K4V;
  
  for(int i=0; i<time_steps; ++i)
    {
      //1st pass
      get_Iion( r_old, V_old, I_stim.at(0), f_r, Iion);
      K1V   = - dt_i/C_m * Iion;
      K1r   =   dt_i * f_r;

      //2nd pass
      get_Iion( r_old+0.5*K1r, V_old+0.5*K1V,
		I_stim.at(1), f_r, Iion);
      K2V   = - dt_i/C_m * Iion;
      K2r   =   dt_i * f_r;

      //3rd pass
      get_Iion( r_old+0.5*K2r, V_old+0.5*K2V,
		I_stim.at(1), f_r, Iion);
      K3V   = - dt_i/C_m * Iion;
      K3r   =   dt_i * f_r;
	

      //4th pass
      get_Iion( r_old+K3r, V_old+K3V,
		I_stim.at(2),  f_r, Iion);
      K4V   = - dt_i/C_m * Iion;
      K4r   =   dt_i * f_r;

      V_new = V_old + (1.0/6.0)*(K1V + 2.0*K2V + 2.0*K3V + K4V);
      r_new = r_old + (1.0/6.0)*(K1r + 2.0*K2r + 2.0*K3r + K4r);

      //update before the next time step 
      V_old=V_new; 
      r_old=r_new;
    }
}

double IonicModel::get_diso() const
{
  //SYS_T::commPrint("IonicModel get_diso %e \n", d_iso);
  return d_iso;
}

double IonicModel::get_dani() const
{
  //  SYS_T::commPrint("IonicModel get_ani  %e \n", d_ani);
  return d_ani;
}

double IonicModel::get_chi() const
{
  //  SYS_T::commPrint("IonicModel get_chi %e \n", chi);
  return chi;
}

double IonicModel::get_C_m() const
{
  //  SYS_T::commPrint("IonicModel C_m %e \n", C_m);
  return C_m;
}




// EOF
