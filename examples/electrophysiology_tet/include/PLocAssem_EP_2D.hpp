#ifndef PLOCASSEM_EP_2D_HPP
#define PLOCASSEM_EP_2D_HPP
// ==================================================================
// PLocAssem_EP_2D.hpp
// This is the local assembly routine for EP equation in
// two dimension, using generalized alpha method as time marching
// scheme. See J.A. Cottrell, et. al Isogeometric Analysis, pp 198
// for details.
// Author: Ju Liu, liujuy@gmail.com
// Modified: Oguz Ziya Tikenogullari, o.z.tikenogullari@gmail.com
//
// ==================================================================
#include "Sys_Tools.hpp"
#include "Math_Tools.hpp"
#include "TimeMethod_GenAlpha.hpp"
#include "AInt_Weight.hpp"
#include "IPLocAssem.hpp"
#include "FEAElement.hpp"
#include "IonicModel.hpp"

class PLocAssem_EP_2D : public IPLocAssem
{
  public:
    PLocAssem_EP_2D(
        const class TimeMethod_GenAlpha * const &tm_gAlpha,
	const class IonicModel * const &ionicmodel,
        const int &in_nlocbas, const int &in_nqp
        );
    virtual ~PLocAssem_EP_2D();

    virtual int get_dof() const {return dof_per_node;}

    virtual void Zero_Tangent_Residual();
    virtual void Zero_Residual();

    virtual void Assem_Estimate();
    
    virtual void Assem_Residual(
        double t_n, double dt,
        const double * const &vec_a,
        const double * const &vec_b,
	//const double * const &vec_c,
        //const double * const &vec_d,
        const class FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const class AInt_Weight * const &weight );


    virtual void Assem_Tangent_Residual(
        double t_n, double dt,
        const double * const &vec_a,
        const double * const &vec_b,
	//const double * const &vec_c,
        //const double * const &vec_d,
        const class FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const class AInt_Weight * const &weight );

    
    virtual void Assem_Mass(
        const class FEAElement * const &element,
        const class AInt_Weight * const &weight );


    virtual void Assem_Mass_Residual(
        const double * const &vec_a,
	//const double * const &vec_b,
	//const double * const &vec_c,    
        const class FEAElement * const &element,
        const double * const &eleCtrlPts_x,
        const double * const &eleCtrlPts_y,
        const double * const &eleCtrlPts_z,
        const class AInt_Weight * const &weight );

  private:
    // ! Necessary data for assembly:
    // generalized-alpha method
    double alpha_f, alpha_m, gamma;
  
    //  conduction coefficients
    double d_iso, d_ani, chi, C_m;
  
    // vec_size = nLocBas * dof_per_node
    int vec_size, nLocBas, dof_per_node;

    // number of quadrature points
    int nqp;

    // ! Dynamic array allocations
    double * R;
    double * dR_dx;
    double * dR_dy;

    // ! Material properties and external functions:
    // ! define the (nonlinear) conductivity tensor
    void get_k( const double &u, const double &x, const double &y,
        double &k11, double &k12, double &k21, double &k22 ) const
    {
      //assume 1 fiber direction is x
      k11 = d_ani+d_iso;    k12 = 0.0;
      k21 = 0.0;            k22 = d_iso;
    }

    // ! define the derivative the conductivity tensor w.r.t. u
    void get_dk_du( const double &u, const double &x, const double &y, 
        double &dk11, double &dk12, double &dk21, double &dk22 ) const
    {
      dk11 = 0.0;  dk12 = 0.0;
      dk21 = 0.0;  dk22 = 0.0;
    }


    // ! define the external  source
    double get_f( const double &x, const double &y, const double &t ) const
    { 
      const double pi = MATH_T::PI;
      //const double a = sin(pi*x);
      //const double b = sin(pi*y);
      //const double pi2 = pi * pi;
      //double val;
      //
      //if (t<=0.5) {
      //	val = ( pi * cos(pi*t) + 2.0*pi2*sin(pi*t) )*a*b ;
      //}
      //else {
      //	val = 0;
      //}
      //	  
      //return val;
      double val=0.0;

      ////MANUFACTURED SOLUTION, constant ionic output:
      //val =140.0*(cos(2.0*pi*x)/4.0 - 1.0/4.0)*(cos(2.0*pi*y) - 1.0)
      //+ (2.0*t*std::pow(pi,2.0)*cos(2.0*pi*y)*(cos(2.0*pi*x)/4.0 - 1.0/4.0))/5.0
      //+ (t*std::pow(pi,2.0)*cos(2.0*pi*x)*(cos(2.0*pi*y) - 1.0))/10.0 - 1400.0;

      ////MANUFACTURED SOLUTION, FitzHugh Nagumo Model:
      //val = (43750.0*cos(2.0*pi*x))/9.0 - (175.0*t)/12.0 + (43750.0*cos(2.0*pi*y))/9.0
      //	- (175.0*exp((3.0*t)/1000.0)
      //	   *(250.0*cos(2.0*pi*x) + 250.0*cos(2.0*pi*y) - 250.0*cos(2.0*pi*x)*cos(2.0*pi*y) - 115.0))/9.0
      //	+ (175.0*t*cos(2.0*pi*x))/12.0 + (175.0*t*cos(2.0*pi*y))/12.0 - (43750.0*cos(2.0*pi*x)*cos(2.0*pi*y))/9.0
      //	+ 140.0*(cos(2.0*pi*x)/4.0 - 1.0/4.0)*(cos(2.0*pi*y) - 1.0)
      //	+ 2275.0*((t*(cos(2.0*pi*x)/4.0 - 1.0/4.0)*(cos(2.0*pi*y) - 1.0))/65.0 - 9.0/13.0)
      //	*((t*(cos(2.0*pi*x)/4.0 - 1.0/4.0)*(cos(2.0*pi*y) - 1.0))/65.0 - 5.0/26.0)
      //	*((t*(cos(2.0*pi*x)/4.0 - 1.0/4.0)*(cos(2.0*pi*y) - 1.0))/65.0 - 22.0/13.0)
      //	- (175.0*t*cos(2.0*pi*x)*cos(2.0*pi*y))/12.0 + (2.0*t*pi*pi*cos(2.0*pi*y)*(cos(2.0*pi*x)/4.0 - 1.0/4.0))/5.0
      //	+ (t*pi*pi*cos(2.0*pi*x)*(cos(2.0*pi*y) - 1.0))/10.0 - 20125.0/9.0;

      ////MANUFACTURED SOLUTION, FitzHugh Nagumo Model wiht c=10 instead of 50:      
      //      val=  (8750.0*cos(2.0*pi*x))/9.0 - (35.0*t)/12.0 + (8750.0*cos(2.0*pi*y))/9.0
      //	- (35.0*exp((3.0*t)/1000.0)*(250.0*cos(2.0*pi*x) + 250.0*cos(2.0*pi*y)
      //				     - 250.0*cos(2.0*pi*x)*cos(2.0*pi*y) - 115.0))/9.0
      //	+ (35.0*t*cos(2.0*pi*x))/12.0 + (35.0*t*cos(2.0*pi*y))/12.0
      //	- (8750.0*cos(2.0*pi*x)*cos(2.0*pi*y))/9.0 + 140.0*(cos(2.0*pi*x)/4.0 - 1.0/4.0)
      //	*(cos(2.0*pi*y) - 1.0) + 455.0*((t*(cos(2.0*pi*x)/4.0 - 1.0/4.0)
      //					 *(cos(2.0*pi*y) - 1.0))/65.0 - 9.0/13.0)
      //	*((t*(cos(2.0*pi*x)/4.0 - 1.0/4.0)*(cos(2.0*pi*y) - 1.0))/65.0 - 5.0/26.0)
      //	*((t*(cos(2.0*pi*x)/4.0 - 1.0/4.0)*(cos(2.0*pi*y) - 1.0))/65.0 - 22.0/13.0)
      //	- (35.0*t*cos(2.0*pi*x)*cos(2.0*pi*y))/12.0
      //	+ (2.0*t*pi*pi*cos(2.0*pi*y)*(cos(2.0*pi*x)/4.0 - 1.0/4.0))/5.0
      //	+ (t*pi*pi*cos(2.0*pi*x)*(cos(2.0*pi*y) - 1.0))/10.0 - 4025.0/9.0;

      return val; 
    } 
};

#endif
