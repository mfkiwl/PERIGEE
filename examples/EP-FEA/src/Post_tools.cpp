#include "Post_tools.hpp"

double POST_T::exact_scalar(const double &x, const double &y, const double &z,
		const double &t)
{
	const double pi = MATH_T::PI;
  
  //const double beta = 0.001 * pi;

  //const double u = t*t*( x * cos(beta*z) - y * sin(beta*z) - x );
  //const double u = t*t*( x * sin(beta*z) + y * cos(beta*z) - y );
  //const double u = 0.0;

  //const double J = 2*t*t*cos(beta*z) - 2*t*t*t*t*cos(beta*z) - 2*t*t + 2*t*t*t*t + 1.0;

  //const double p = 0.5 * 4.009433333333775e+06 * (J - 1.0 / J);

  return t*t*sin(pi*x)*sin(pi*y)*sin(pi*z);
}


double POST_T::exact_scalar(const double &x, const double &y, const double &t)
{
  const double pi = MATH_T::PI;
  
  const double V = 1.0/4.0*(1-cos(2.0*pi*x))*(1-cos(2.0*pi*y))*t-80;

  return V;
}


void POST_T::exact_3dvector(const double &x, const double &y, const double &z,
		const double &t, double &val_x, double &val_y, double &val_z )
{
  const double pi = MATH_T::PI;

  val_x = pi*t*t*cos(pi*x)*sin(pi*y)*sin(pi*z);
  val_y = pi*t*t*cos(pi*y)*sin(pi*x)*sin(pi*z);
  val_z = pi*t*t*cos(pi*z)*sin(pi*x)*sin(pi*y);
}


void POST_T::exact_2dvector(const double &x, const double &y,
    const double &t, double &val_x, double &val_y )
{
  const double pi = MATH_T::PI;

  const double V_x = 1.0/4.0*(1-cos(2.0*pi*y))*t*sin(2*pi*x)*2*pi;
  const double V_y = 1.0/4.0*(1-cos(2.0*pi*x))*t*sin(2*pi*y)*2*pi;

  val_x = V_x;
  val_y = V_y;
}


double POST_T::get_manu_scalar_l2_error(
    const double * const &sol,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const AInt_Weight * const &weight,
    double * const &R,
    const double &curr )
{
  double error = 0.0;

  double exa_qua, sol_qua;
  double coor_x, coor_y, coor_z;

  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  for(int qua=0; qua<nqp; ++qua)
  {
    double gwts = element->get_detJac(qua) * weight->get_weight(qua);

    sol_qua = 0.0;
    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    element->get_R(qua, R);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];
      sol_qua += sol[ii] * R[ii];
    } 

    exa_qua = exact_scalar(coor_x, coor_y, coor_z, curr);

    error += (sol_qua - exa_qua)*(sol_qua - exa_qua) * gwts;
  }

  return error;
}


double POST_T::get_manu_scalar_l2_error(
    const double * const &sol,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const AInt_Weight * const &weight,
    double * const &R,
    const double &curr )
{
  double error = 0.0;

  double exa_qua, sol_qua;
  double coor_x, coor_y;

  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  for(int qua=0; qua<nqp; ++qua)
  {
    double gwts = element->get_detJac(qua) * weight->get_weight(qua);

    sol_qua = 0.0;
    coor_x = 0.0; coor_y = 0.0;

    element->get_R(qua, R);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      sol_qua += sol[ii] * R[ii];
    } 

    exa_qua = exact_scalar(coor_x, coor_y, curr);

    error += (sol_qua - exa_qua)*(sol_qua - exa_qua) * gwts;
  }

  return error;
}


double POST_T:: get_manu_scalar_h1_error(
    const double * const &solu,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const double * const &ectrlPts_z,
    const AInt_Weight * const &weight, 
    double * const &R,
    double * const &dR_dx,
    double * const &dR_dy,
    double * const &dR_dz,
    const double &curr )
{
  double error = 0.0;
  double exa_x, exa_y, exa_z, sol_x, sol_y, sol_z, exa, sol;
  double coor_x, coor_y, coor_z;
  double gwts;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  for(int qua=0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    sol = 0.0; sol_x = 0.0; sol_y = 0.0; sol_z = 0.0;
    exa = 0.0; exa_x = 0.0; exa_y = 0.0; exa_z = 0.0;
    coor_x = 0.0; coor_y = 0.0; coor_z = 0.0;

    element->get_R(qua, R);
    element->get_gradR(qua, dR_dx, dR_dy, dR_dz);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];
      sol    += solu[ii] * R[ii];
      sol_x  += solu[ii] * dR_dx[ii];
      sol_y  += solu[ii] * dR_dy[ii];
      sol_z  += solu[ii] * dR_dz[ii];
    }

    exa = exact_scalar(coor_x, coor_y, coor_z, curr);
    exact_3dvector(coor_x, coor_y, coor_z, curr, exa_x, exa_y, exa_z);

    error += (sol - exa) * (sol - exa) * gwts;
    error += (sol_x - exa_x) * (sol_x - exa_x) * gwts;
    error += (sol_y - exa_y) * (sol_y - exa_y) * gwts;
    error += (sol_z - exa_z) * (sol_z - exa_z) * gwts;
  }

  return error;
}


double POST_T:: get_manu_scalar_h1_error(
    const double * const &solu,
    const FEAElement * const &element,
    const double * const &ectrlPts_x,
    const double * const &ectrlPts_y,
    const AInt_Weight * const &weight, 
    double * const &R,
    double * const &dR_dx,
    double * const &dR_dy,
    const double &curr )
{
  double error = 0.0;
  double exa_x, exa_y, sol_x, sol_y, exa, sol;
  double coor_x, coor_y;
  double gwts;
  int nqp = weight->get_num();
  int nLocBas = element->get_nLocBas();

  for(int qua=0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    sol = 0.0; sol_x = 0.0; sol_y = 0.0;
    exa = 0.0; exa_x = 0.0; exa_y = 0.0;
    coor_x = 0.0; coor_y = 0.0;

    element->get_R_gradR(qua, R, dR_dx, dR_dy);

    for(int ii=0; ii<nLocBas; ++ii)
    {
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      sol    += solu[ii] * R[ii];
      sol_x  += solu[ii] * dR_dx[ii];
      sol_y  += solu[ii] * dR_dy[ii];
    }

    exa = exact_scalar(coor_x, coor_y, curr);
    exact_2dvector(coor_x, coor_y, curr, exa_x, exa_y);

    error += (sol - exa) * (sol - exa) * gwts;
    error += (sol_x - exa_x) * (sol_x - exa_x) * gwts;
    error += (sol_y - exa_y) * (sol_y - exa_y) * gwts;
  }

  return error;
}

double POST_T::calculate_ecg(const double * const &solu,
			     const FEAElement * const &element,
			     const IQuadPts * const &quadPtr,
			     const double * const &ectrlPts_x,
			     const double * const &ectrlPts_y,
			     const double * const &ectrlPts_z,
			     const std::vector<double> &xe,
			     //double * const &R,
			     //double * const &dR_dx,
			     //double * const &dR_dy,
			     //double * const &dR_dz,
			     const int &nlocbas_ee,
			     const double &curr ) {

  int nqp = quadPtr->get_num_quadPts();
  //std::cout << "nqp=" << nqp << std::endl;

  double JxW;
  double solQ;
  double Integral;
  double dsolQ_x=0;  double  dsolQ_y=0;   double dsolQ_z=0;
  double coor_x =0;  double  coor_y =0;   double coor_z =0;
  double R[nlocbas_ee];
  double Rx[nlocbas_ee];
  double Ry[nlocbas_ee];
  double Rz[nlocbas_ee];
  
  for(int qua=0; qua<nqp; ++qua)    {

    solQ=0;
    dsolQ_x=0;    dsolQ_y=0;    dsolQ_z=0;
    coor_x =0;    coor_y =0;    coor_z =0;

    element->get_R_gradR(qua, R, Rx, Ry, Rz);

    // Inner product with basis functions
    // ! This loop may be speed up by calling the sdot function in cblas
//    if (qua==0) {
//      std::cout << "ectrlpts_x: \n" ;
//      for (int ii=0; ii<nlocbas_ee; ++ii)	{
//	std::cout << ectrlPts_x[ii] << "\t" ;
//      }
//      std::cout << std::endl;
//
//      std::cout << "ectrlpts_y: \n" ;
//      for (int ii=0; ii<nlocbas_ee; ++ii)	{
//	std::cout << ectrlPts_y[ii] << "\t" ;
//      }
//      std::cout << std::endl;
//
//      std::cout << "ectrlpts_z: \n" ;
//      for (int ii=0; ii<nlocbas_ee; ++ii)	{
//	std::cout << ectrlPts_z[ii] << "\t" ;
//      }
//      std::cout << std::endl;
//
//      std::cout << "R: \n" ;
//      for (int ii=0; ii<nlocbas_ee; ++ii)	{
//	std::cout << R[ii] << "\t" ;
//      }
//      std::cout << std::endl;
//
//      std::cout << "Rx: \n" ;
//      for (int ii=0; ii<nlocbas_ee; ++ii)	{
//	std::cout << Rx[ii] << "\t" ;
//      }
//      std::cout << std::endl;
//
//      std::cout << "Ry: \n" ;
//      for (int ii=0; ii<nlocbas_ee; ++ii)	{
//	std::cout << Ry[ii] << "\t" ;
//      }
//      std::cout << std::endl;
//
//      std::cout << "Rz: \n" ;
//      for (int ii=0; ii<nlocbas_ee; ++ii)	{
//	std::cout << Rz[ii] << "\t" ;
//      }
//      std::cout << std::endl;
//    }
    
    for(int ii=0; ii<nlocbas_ee; ++ii)	{

      solQ   += solu[ii] * R[ii];
      
      dsolQ_x += solu[ii] * Rx[ii];
      dsolQ_y += solu[ii] * Ry[ii];
      dsolQ_z += solu[ii] * Rz[ii];
      coor_x += ectrlPts_x[ii] * R[ii];
      coor_y += ectrlPts_y[ii] * R[ii];
      coor_z += ectrlPts_z[ii] * R[ii];
    }

    double r_sq = (xe[0]*xe[0] + xe[1]*xe[1] + xe[2]*xe[2]
		   - 2*(xe[0]*coor_x + xe[1]*coor_y + xe[2]*coor_z)
		   + coor_x*coor_x + coor_y*coor_y + coor_z*coor_z ) ;

    double drinv_x = -1./2.* std::pow(r_sq, -3./2.)
      * (-2.0 * xe[0] + 2.0 * coor_x);
    double drinv_y = -1./2.* std::pow(r_sq, -3./2.)
      * (-2.0 * xe[1] + 2.0 * coor_y);
    double drinv_z = -1./2.* std::pow(r_sq, -3./2.)
      * (-2.0 * xe[2] + 2.0 * coor_z);
	    
    
    JxW = element->get_detJac(qua) * quadPtr->get_qw(qua);

    //    get_k(d, coor_x, coor_y, coor_z, fiber_ori_e,
    //	  k11, k12, k13, k21, k22, k23, k31, k32, k33);


    Integral = 0;
    Integral += JxW * ( -dsolQ_x * drinv_x
			-dsolQ_y * drinv_y
			-dsolQ_z * drinv_z );
      
    //  //	*(chi * C_m * R[ii] * v + k11 * d_x * dR_dx[ii]
    //  //	 + k12 * d_y * dR_dx[ii] + k13 * d_z * dR_dx[ii]
    //  //	 + k21 * d_x * dR_dy[ii] + k22 * d_y * dR_dy[ii]
    //  //	 + k23 * d_z * dR_dy[ii]
    //  //	 + k31 * d_x * dR_dz[ii]
    //  //	 + k32 * d_y * dR_dz[ii] + k33 * d_z * dR_dz[ii] );
    //}

    
  }

  return Integral;
}

// EOF
