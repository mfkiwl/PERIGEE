#include "IonicModel_TTP.hpp"

IonicModel_TTP::IonicModel_TTP()
  //           d_iso,     d_ani,     chi,   C_m, n_int_vars
  : IonicModel(0.012*1.0, 0.078*1.0, 140.0, 0.01, 18),
    Rc{8.314}, Tc{310.0e0}, Fc{96.485},
    rho{162.0e0}, V_c{16.404e-3}, V_sr{1.094e-3}, V_ss{5.468e-5},
    K_o{5.4e0}, Na_o{140.0e0}, Ca_o{2.0e0}, G_Na{14.838e0},
    G_K1{5.405e0}, G_to{0.073e0, 0.294e0, 0.294e0}, G_Kr{0.153e0},
    G_Ks{0.392e0, 0.392e0, 0.098e0}, p_KNa{0.03e0}, G_CaL{3.98e-5},
    K_NaCa{1000.0e0}, gamma{0.35e0}, K_mCa{1.38e0}, K_mNai{87.5e0},
    K_sat{0.1e0}, alpha{2.5e0}, P_NaK{2.724e0}, K_mK{1.0e0},
    K_mNa{40.0e0}, G_pK{1.46e-2}, G_pCa{0.1238e0}, K_pCa{5.0e-4},
    G_bNa{2.9e-4}, G_bCa{5.92e-4}, Vmax_up{6.375e-3}, K_up{2.5e-4},
    V_rel{0.102e0}, k1p{0.15e0}, k2p{0.045e0}, k3{0.06e0}, k4{5.0e-3},
    EC{1.5e0}, max_sr{2.5e0}, min_sr{1.0e0}, V_leak{3.6e-4},
    V_xfer{3.8e-3}, Buf_c{0.2e0}, K_bufc{1.0e-3}, Buf_sr{10.0e0},
    K_bufsr{0.3e0}, Buf_ss{0.4e0}, K_bufss{2.5e-4}
{
  //SYS_T::commPrint("AP constructor. \n");
};

IonicModel_TTP::~IonicModel_TTP()
{
  //SYS_T::commPrint("AP destructor. \n");
};

void IonicModel_TTP::print_info () const
{
  PetscPrintf(PETSC_COMM_WORLD, "\t  Ten-Tusscher Panfilov EP \n");
  //PetscPrintf(PETSC_COMM_WORLD, "\t  Rc = %.2e \n", Rc);
  //PetscPrintf(PETSC_COMM_WORLD, "\t  G_to = [%.2e, %.2e, %.2e] \n",
  //                  G_to[0], G_to[1],G_to[2] ); 
};

void IonicModel_TTP::run_ionic(const std::vector<double> &r_old_in,
			       const double &V_old_in,
			       const double &I_stim,
			       const double &dt_in,
			       std::vector<double> &r_new,
			       double &V_new) const
{

  int i=2; // i=0 for endo, 1 for myo, 2 for epi
  
  double V, K_i, Na_i, Ca_i, Ca_ss, Ca_sr, R_bar, xr1, xr2, xs, m, h, j, d, f,
    f2, fcass, s, r, xr1i, a, b, tau, xr2i, xsi, mi, hi, ji, di, c, fi, f2i,
    fcassi, si, ri, RT, E_K, E_Na, E_Ca, E_Ks, I_Na, I_to, e1, e2, e3, e4,
    alpK1, betK1, xK1i, sq5, I_K1, I_Kr, I_Ks, I_CaL, n1, d1, d2, d3, I_NaCa,
    I_NaK, I_pCa, I_pK, I_bCa, I_bNa, I_leak, I_up, k_casr, k1, O,
    I_rel, I_xfer, I_ion, dV, dK_i, dNa_i, n2, dCa_i, dCa_ss, dCa_sr,
    k2, dR_bar;
  
  // Local copies of state variables
  V     = V_old_in        ;
  K_i   = r_old_in.at( 0 );
  Na_i  = r_old_in.at( 1 );
  Ca_i  = r_old_in.at( 2 );
  Ca_ss = r_old_in.at( 3 );
  Ca_sr = r_old_in.at( 4 );
  R_bar = r_old_in.at( 5 );
  // Local copies of gating variables
  xr1   = r_old_in.at(  6 );
  xr2   = r_old_in.at(  7 );
  xs    = r_old_in.at(  8 );
  m     = r_old_in.at(  9 );
  h     = r_old_in.at( 10 );
  j     = r_old_in.at( 11 );
  d     = r_old_in.at( 12 );
  f     = r_old_in.at( 13 );
  f2    = r_old_in.at( 14 );
  fcass = r_old_in.at( 15 );
  s     = r_old_in.at( 16 );
  r     = r_old_in.at( 17 );

  
  //Calculate and update the gating variables.
  //xr1: activation gate for I_Kr
  xr1i   = 1.0/(1.0 + exp(-(26.0+V)/7.0));
  a      = 450.0/(1.0 + exp(-(45.0+V)/10.0));
  b      = 6.0/(1.0 + exp((30.0+V)/11.5));
  tau    = a*b;
  xr1    = xr1i - (xr1i - xr1)*exp(-dt_in/tau);

  //xr2: inactivation gate for I_Kr
  xr2i   = 1.0 /(1.0 + exp((88.0+V)/24.0));
  a      = 3.0 /(1.0 + exp(-(60.0+V)/20.0));
  b      = 1.12/(1.0 + exp(-(60.0-V)/20.0));
  tau    = a*b;
  xr2    = xr2i - (xr2i - xr2)*exp(-dt_in/tau);

  //xs: activation gate for I_Ks
  xsi    = 1.0/(1.0 + exp(-(5.0+V)/14.0));
  a      = 1400.0/sqrt(1.0 + exp((5.0-V)/6.0));
  b      = 1.0/(1.0 + exp((V-35.0)/15.0));
  tau    = a*b + 80.0;
  xs     = xsi - (xsi - xs)*exp(-dt_in/tau);

  //m: activation gate for I_Na
  mi     = 1.0/( pow(1.0 + exp(-(56.86+V)/9.03),2.0) );
  a      = 1.0/(1.0 + exp(-(60.0+V)/5.0));
  b      = 0.1/(1.0 + exp((35.0+V)/5.0))
    + 0.1/(1.0 + exp((V-50.0)/200.0));
  tau    = a*b;
  m      = mi - (mi - m)*exp(-dt_in/tau);
  // m= m+ dt_in*(a*(1-m) - b*m);

  //h: fast inactivation gate for I_Na
  hi  = 1.0/( pow(1.0 + exp((71.55+V)/7.43),2.0) );
  if (V >= -40.0) {
    a   = 0.0;
    b   = 0.77/(0.13*(1.0 + exp(-(10.66+V)/11.1)));
  }
  else{
    a   = 5.7e-2*exp(-(80.0+V)/6.8);
    b   = 2.7*exp(0.079*V) + 310000.0*exp(0.3485*V);
  }
 
  tau = 1.0 / (a + b);
  h   = hi - (hi - h)*exp(-dt_in/tau);

  //j: slow inactivation gate for I_Na
  ji  = 1.0/( pow(1.0 + exp((71.55+V)/7.43),2.0) );
  if (V >= -40.0) {
    a   = 0.0;
    b   = 6e-1*exp(5.7e-2*V)/(1.0 + exp(-0.1*(V+32.0)));
  }
  else{
    a   = -(25428.0*exp(0.2444*V) + 6.948e-6*exp(-0.04391*V))
      * (V+37.78) / (1.0 + exp(0.311*(79.23+V)));
    b   = 0.02424*exp(-0.01052*V) / (1.0 + exp(-0.1378*(40.14+V))); 
  }
  tau = 1.0 / (a + b);
  j   = ji - (ji - j)*exp(-dt_in/tau);

  //d: activation gate for I_CaL
  di     = 1.0/(1.0 + exp(-(8.0+V)/7.5));
  a      = 1.4/(1.0 + exp(-(35.0+V)/13.0)) + 0.25;
  b      = 1.4/(1.0 + exp((5.0+V)/5.0));
  c      = 1.0/(1.0 + exp((50.0-V)/20.0));
  tau    = a*b + c;
  d      = di - (di - d)*exp(-dt_in/tau);

  //f: slow inactivation gate for I_CaL
  fi     = 1.0/(1.0 + exp((20.0+V)/7.0));
  a      = 1102.5*exp(-(pow(V+27.0,2.0))/225.0);
  b      = 200.0/(1.0 + exp((13.0-V)/10.0));
  c      = 180.0/(1.0 + exp((30.0+V)/10.0)) + 20.0;
  tau    = a + b + c;
  // % c!     for spiral wave breakup
  // % c      IF (V .GT. 0) tau = tau*2.0
  f      = fi - (fi - f)*exp(-dt_in/tau);

  //f2: fast inactivation gate for I_CaL
  f2i    = 0.67/(1.0 + exp((35.0+V)/7.0)) + 0.33;
  a      = 562.0*exp(-(pow(27.0+V,2.0)) /240.0);
  b      = 31.0/(1.0 + exp((25.0-V)/10.0));
  c      = 80.0/(1.0 + exp((30.0+V)/10.0));
  tau    = a + b + c;
  f2     = f2i - (f2i - f2)*exp(-dt_in/tau);

  //fCass: inactivation gate for I_CaL into subspace
  c      = 1.0/(1.0 + pow(Ca_ss/0.05,2.0));
  fcassi = 0.6*c  + 0.4;
  tau    = 80.0*c + 2.0;
  fcass  = fcassi - (fcassi - fcass)*exp(-dt_in/tau);
	      
  // s: inactivation gate for I_to
  if ((i==0)||(i==2)) {
    si  = 1.0/(1.0 + exp((20.0+V)/5.0));
    tau = 85.0*exp(-(pow(V+45.0,2.0)) /320.0) 
      + 5.0/(1.0+exp((V-20.0)/5.0)) + 3.0;
  }
  else if (i == 1){
    si  = 1.0/(1.0 + exp((28.0+V)/5.0));
    tau = 1000.0*exp(-(pow(V+67.0,2.0)) /1000.0) + 8.0;
  } else {
    SYS_T::print_exit("IonicModel_TTP: i should be  1, 2 or 3");
  }
  
  s   = si - (si - s)*exp(-dt_in/tau);

  // r: activation gate for I_to
  ri     = 1.0/(1.0 + exp((20.0-V)/6.0));
  tau    = 9.5*exp(-(pow(V+40.0,2.0)) /1800.0) + 0.8;
  r      = ri - (ri - r)*exp(-dt_in/tau);

  //UPDATE GATING VARIABLES
  r_new.at(  6 ) = xr1  ;
  r_new.at(  7 ) = xr2  ;
  r_new.at(  8 ) = xs   ;
  r_new.at(  9 ) = m    ;
  r_new.at( 10 ) = h    ;
  r_new.at( 11 ) = j    ;
  r_new.at( 12 ) = d    ;
  r_new.at( 13 ) = f    ;
  r_new.at( 14 ) = f2   ;
  r_new.at( 15 ) = fcass;
  r_new.at( 16 ) = s    ;
  r_new.at( 17 ) = r    ;


  // START CALCULATING CURRENTS AND TIME RATES
  RT   = Rc * Tc / Fc;
  E_K  = RT * log(K_o/K_i);
  E_Na = RT * log(Na_o/Na_i);
  E_Ca = 0.5 * RT * log(Ca_o/Ca_i);
  E_Ks = RT * log( (K_o + p_KNa*Na_o)/(K_i + p_KNa*Na_i) );

  //I_Na: Fast sodium current
  I_Na = G_Na * pow(m,3.0) * h * j * (V - E_Na);

  //I_to: transient outward current
  I_to = G_to.at(i) * r * s * (V - E_K);

  //I_K1: inward rectifier outward current
  e1   = exp(0.06*(V - E_K - 200.0));
  e2   = exp(2.0e-4*(V - E_K + 100.0));
  e3   = exp(0.1*(V - E_K - 10.0));
  e4   = exp(-0.5*(V - E_K));
  alpK1 = 0.1/(1.0 + e1);
  betK1 = (3.0*e2 + e3) / (1.0 + e4);
  xK1i  = alpK1 / (alpK1 + betK1);
  sq5  = sqrt(K_o/5.4);
  I_K1 = G_K1 * sq5 * xK1i * (V - E_K);

  //I_Kr: rapid delayed rectifier current
  I_Kr = G_Kr * sq5 * xr1 * xr2 * (V - E_K);

  //I_Ks: slow delayed rectifier current 
  // this is different between '04 and '06 models
  I_Ks = G_Ks.at(i) * pow(xs,2.0) * (V - E_Ks);

  //I_CaL: L-type Ca current
  a     = 2.0*(V-15.0)/RT;
  b     = 2.0*a*Fc * (0.25*Ca_ss*exp(a) - Ca_o) / (exp(a)-1.0);
  I_CaL = G_CaL * d * f * f2 * fcass * b;

  //I_NaCa: Na-Ca exchanger current
  e1     = exp(gamma*V/RT);
  e2     = exp((gamma-1.0)*V/RT);
  n1     = e1*pow(Na_i,3.0)*Ca_o - e2*pow(Na_o,3.0)*Ca_i*alpha;
  d1     = pow(K_mNai,3.0) + pow(Na_o,3.0);
  d2     = K_mCa + Ca_o;
  d3     = 1.0 + K_sat*e2;
  I_NaCa = K_NaCa * n1 / (d1*d2*d3);

  //I_NaK: Na-K pump current
  e1    = exp(-0.1*V/RT);
  e2    = exp(-V/RT);
  n1    = P_NaK * K_o * Na_i;
  d1    = K_o + K_mK;
  d2    = Na_i + K_mNa;
  d3    = 1.0 + 0.1245*e1 + 0.0353*e2;
  I_NaK = n1 / (d1*d2*d3);

  //I_pCa: plateau Ca current
  I_pCa = G_pCa * Ca_i / (K_pCa + Ca_i);

  //I_pK: plateau K current
  I_pK  = G_pK * (V-E_K) / (1.0 + exp((25.0-V)/5.98));

  //I_bCa: background Ca current
  I_bCa = G_bCa * (V - E_Ca);

  //I_bNa: background Na current
  I_bNa = G_bNa * (V - E_Na);

  //I_leak: Sacroplasmic Reticulum Ca leak current
  I_leak = V_leak * (Ca_sr - Ca_i);

  //I_up: Sacroplasmic Reticulum Ca pump current
  I_up  = Vmax_up / (1.0 + pow((K_up/Ca_i),2.0));

  //I_rel: Ca induced Ca current (CICR)
  k_casr = max_sr - ((max_sr-min_sr)/(1.0 + pow((EC/Ca_sr),2.0)) );
  k1     = k1p / k_casr;
  O      = k1 * R_bar * pow(Ca_ss,2.0) / (k3 + k1*pow(Ca_ss,2.0));
  I_rel  = V_rel * O * (Ca_sr - Ca_ss);

  //I_xfer: diffusive Ca current between Ca subspae and cytoplasm
  I_xfer = V_xfer * (Ca_ss - Ca_i);

  I_ion = -( I_stim  +  //no I_sac
	     I_K1   + 
	     I_Na   + 
	     I_to   + 
	     I_Kr   + 
	     I_Ks   + 
	     I_CaL  + 
	     I_NaCa + 
	     I_NaK  + 
	     I_pCa  + 
	     I_pK   + 
	     I_bCa  + 
	     I_bNa         );
  
  //-----------------------------------------------------------------------
  //Now compute time derivatives

  //dV/dt
  dV    = I_ion;

  //dK_i/dt
  dK_i  = -(C_m/(V_c*Fc))
    * (I_K1+ I_to + I_Kr + I_Ks + I_pK -2.0*I_NaK + I_stim);

  //dNa_i/dt
  dNa_i  = -(C_m/(V_c*Fc)) * (I_Na + I_bNa + 3.0*(I_NaK + I_NaCa));

  //dCa_i/dt
  n1     = (I_leak - I_up)*V_sr/V_c + I_xfer;
  n2     = -(C_m/(V_c*Fc)) * (I_bCa + I_pCa - 2.0*I_NaCa) /2.0;
  d1     = 1.0 + K_bufc*Buf_c/pow((Ca_i + K_bufc),2.0);
  dCa_i  = (n1 + n2)/d1;

  //dCa_ss: rate of change of Ca_ss
  n1     = (-I_CaL*C_m/(2.0*Fc) + I_rel*V_sr - V_c*I_xfer) / V_ss;
  d1     = 1.0 + K_bufss*Buf_ss/(Ca_ss + K_bufss)*2.0;
  dCa_ss = n1 / d1;

  //dCa_sr: rate of change of Ca_sr
  n1     = I_up - I_leak - I_rel;
  d1     = 1.0 + K_bufsr*Buf_sr/pow((Ca_sr + K_bufsr),2.0);
  dCa_sr = n1 / d1;

  //Rbar: ryanodine receptor
  k2     = k2p * k_casr;
  dR_bar = -k2*Ca_ss*R_bar + k4*(1.0 - R_bar);

  //Update the state variables from the derivatives.
  V_new         = V     + dt_in*dV     ;
  r_new.at( 0 ) = K_i   + dt_in*dK_i   ;
  r_new.at( 1 ) = Na_i  + dt_in*dNa_i  ;
  r_new.at( 2 ) = Ca_i  + dt_in*dCa_i  ;
  r_new.at( 3 ) = Ca_ss + dt_in*dCa_ss ;
  r_new.at( 4 ) = Ca_sr + dt_in*dCa_sr ;
  r_new.at( 5 ) = R_bar + dt_in*dR_bar ;
 
};


void IonicModel_TTP::get_Istim(double &Istim,
			       const double &t,
			       const double &x,
			       const double &y,
			       const double &z ) const
{
  //if ( ( x>1.5 ) || ( y>1.5 ) ) {
  //  if ((t >= 0.0) && (t <= 1.0)) {
  //    Istim = -52.0;
  //  }
  //}
  //else {
  //  Istim = 0.0 ;
  //}
  Istim = 0.0 ;
};

void IonicModel_TTP::get_int_vars(double* val) const
{
  //SYS_T::commPrint("TTP ionic, get int vars. \n"); 
  val[ 0 ]= 138.4    ;//K_i    
  val[ 1 ]= 10.355   ;//Na_i   
  val[ 2 ]= 1.3e-4   ;//Ca_i   
  val[ 3 ]= 3.6e-4   ;//Ca_ss  
  val[ 4 ]= 3.715    ;//Ca_sr  
  val[ 5 ]= 9.068e-1 ;//R_bar  
  val[ 6 ]= 4.48e-3  ;//x_r1   
  val[ 7 ]= 4.76e-1  ;//x_r2   
  val[ 8 ]= 8.7e-3   ;//x_s    
  val[ 9 ]= 1.55e-3  ;//m      
  val[10 ]= 7.573e-1 ;//h      
  val[11 ]= 7.225e-1 ;//j      
  val[12 ]= 3.164e-5 ;//d      
  val[13 ]= 8.009e-1 ;//f      
  val[14 ]= 9.778e-1 ;//f_2    
  val[15 ]= 9.953e-1 ;//f_cass 
  val[16 ]= 3.212e-1 ;//s      
  val[17 ]= 2.235e-8 ;//r      
};

// EOF
