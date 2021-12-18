#ifndef IONICMODEL_TTPPURKINJE_HPP
#define IONICMODEL_TTPPURKINJE_HPP
// ==================================================================
// IonicModel_TTPPurkinje.hpp
// 
// TenTusscher-Panfilov model for purkinje fibers
// From the paper : Ten Tusscher and Panfilov 2007, Modelling of the ventricular
// conduction system
//
// Differences compared to ventricle model are:
/// 0.35*G_Ks , 2.94*G_Na, 8.0*conductivity
//
// Date: 2021
// Author: Oguz Ziya Tikenogullari
// Contact: o.z.tikenogullari@gmail.com
// ==================================================================
#include "Sys_Tools.hpp"
#include "Tensor4_3D.hpp"
#include "Math_Tools.hpp"
#include "IonicModel.hpp"

class IonicModel_TTPPurkinje : public IonicModel
{
public:
  IonicModel_TTPPurkinje(const double &cond_scale,
		 const double &delay_in);

  IonicModel_TTPPurkinje(const double &cond_scale);
  
  virtual ~IonicModel_TTPPurkinje();
  
  virtual void print_info() const; //const =0 ;

  virtual void get_int_vars( double* val) const;

protected:

  virtual void run_ionic(const std::vector<double> &r_old_in,
			 const double &V_old_in,
			 const double &I_stim,
			 const double &dt_in,
			 std::vector<double> &r_new,
			 double &V_new) const;  

  virtual void get_Istim(double &Istim,
			 const double &t,
			 const double &x,
			 const double &y,
			 const double &z ) const;
private:

  //Cm and sV are the member variables of the parent class,
  //IonicModel.cpp(.hpp)  called C_m and chi. so we assign them in the
  //constructor.
  
  //Cm: Cell capacitance per unit surface area
  //const double Cm =0.185e0;;         // units: uF/cm^{2}  
  ////sV: Surface to volume ratio
  //const double sV =0.2e0;            // units: um^{-1}

  //R: Gas constant
  const double Rc ;       // units: J/mol/K
  //T: Temperature
  const double Tc ;          // units: K
  //F: Faraday constant
  const double Fc ;          // units: C/mmol
    //rho: Cellular resistivity
  const double rho ;        // units: \Omega-cm
  //V_c: Cytoplasmic volume
  const double V_c ;      // units: um^{3}
  //V_sr: Sacroplasmic reticulum volume
  const double V_sr ;      // units: um^{3}
  //V_ss: Subspace volume
  const double V_ss ;      // units: um^{3}
  //K_o: Extracellular K concentration
  const double K_o ;          // units: mM
  //Na_o: Extracellular Na concentration
  const double Na_o ;       // units: mM
  //Ca_o: Extracellular Ca concentration
  const double Ca_o ;         // units: mM
  //G_Na: Maximal I_Na conductance
  const double G_Na ;      // units: nS/pF
  //G_K1: Maximal I_K1 conductance
  const double G_K1 ;       // units: nS/pF
  //G_to: Maximal epicardial I_to conductance, units: nS/pF
  const std::vector<double> G_to ;
  //G_Kr: Maximal I_Kr conductance
  const double G_Kr ;       // units: nS/pF
  ////////G_Kr for spiral wave breakup
  ////// G_Kr = 0.172e0       // units: nS/pF
  //G_Ks: Maximal epicardial I_Ks conductance, units: nS/pF
  const std::vector<double> G_Ks ;
  ////////G_Ks for spiral wave breakup (epi)
  ////// G_Ks(3) = (/0.441e0, 0.392e0, 0.098e0/)
  //p_KNa: Relative I_Ks permeability to Na
  const double p_KNa ;       // dimensionless
  //G_CaL: Maximal I_CaL conductance
  const double G_CaL ;      // units: cm^{3}/uF/ms
  //K_NaCa: Maximal I_NaCa
  const double K_NaCa ;    // units: pA/pF
  //gamma: Voltage dependent parameter of I_NaCa
  const double gamma ;       // dimensionless
  //K_mCa: Ca_i half-saturation constant for I_NaCa
  const double K_mCa ;       // units: mM
  //K_mNai: Na_i half-saturation constant for I_NaCa
  const double K_mNai ;      // units: mM
  //K_sat: Saturation factor for I_NaCa
  const double K_sat ;        // dimensionless
  //alpha: Factor enhancing outward nature of I_NaCa
  const double alpha ;        // dimensionless
  //p_NaK: Maximal I_NaK
  const double P_NaK ;      // units: pA/pF
  //K_mK: K_o half-saturation constant of I_NaK
  const double K_mK ;         // units: mM
  //K_mNa: Na_i half-saturation constant of I_NaK
  const double K_mNa ;       // units: mM
  //G_pK: Maximal I_pK conductance
  const double G_pK ;       // units: nS/pF
  ////////G_pK for spiral wave breakup
  ////// G_pK = 2.19e-3       // units: nS/pF
  //G_pCa: Maximal I_pCa conductance
  const double G_pCa ;     // units: pA/pF
  ////////G_pCa for spiral wave breakup
  ////// G_pCa = 0.8666e0     // units: pA/pF
  //K_pCa: Half-saturation constant of I_pCa
  const double K_pCa ;       // units: mM
  //G_bNa: Maximal I_bNa conductance
  const double G_bNa ;       // units: nS/pF
  //G_bCa: Maximal I_bCa conductance
  const double G_bCa ;      // units: nS/pF
  //Vmax_up: Maximal I_up conductance
  const double Vmax_up ;   // units: mM/ms
  //K_up: Half-saturation constant of I_up
  const double K_up ;        // units: mM
  //V_rel: Maximal I_rel conductance
  const double V_rel ;      // units: mM/ms
  //k1p: R to O and RI to I, I_rel transition rate
  const double k1p ;         // units: mM^{-2}/ms
  //k2p: O to I and R to RI, I_rel transition rate
  const double k2p ;        // units: mM^{-1}/ms
  //k3: O to R and I to RI, I_rel transition rate
  const double k3 ;         // units: ms^{-1}
  //k4: I to O and Ri to I, I_rel transition rate
  const double k4 ;         // units: ms^{-1}
  //EC: Ca_sr half-saturation constant of k_casr
  const double EC ;          // units: mM
  //max_sr: Maximum value of k_casr
  const double max_sr ;       // dimensionless
  //min_sr: Minimum value of k_casr
  const double min_sr ;       // dimensionless
  //V_leak: Maximal I_leak conductance
  const double V_leak ;      // units: mM/ms
  //V_xfer: Maximal I_xfer conductance
  const double V_xfer ;      // units: mM/ms
  //Buf_c: Total cytoplasmic buffer concentration
  const double Buf_c ;        // units: mM
  //K_bufc: Ca_i half-saturation constant for cytplasmic buffer
  const double K_bufc ;      // units: mM
  //Buf_sr: Total sacroplasmic buffer concentration
  const double Buf_sr ;      // units: mM
  //K_bufsr: Ca_sr half-saturation constant for subspace buffer
  const double K_bufsr ;      // units: mM
  //Buf_ss: Total subspace buffer concentration
  const double Buf_ss ;       // units: mM
  //K_bufss: Ca_ss half-saturation constant for subspace buffer
  const double K_bufss ;     // units: mM
  //Resting potential
  //   Vrest = don't set it here 

  const double pur_delay; 

};

#endif


