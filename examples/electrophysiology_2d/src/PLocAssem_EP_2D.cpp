#include "PLocAssem_NLHeat_2D_GenAlpha.hpp"

PLocAssem_NLHeat_2D_GenAlpha::PLocAssem_NLHeat_2D_GenAlpha(
    const class TimeMethod_GenAlpha * const &tm_gAlpha,
    const class IonicModel * const &ionicmodel,
    const int &in_locbas, const int &in_nqp )
{
  alpha_m = tm_gAlpha->get_alpha_m();
  alpha_f = tm_gAlpha->get_alpha_f();
  gamma   = tm_gAlpha->get_gamma();
  d_iso   =ionicmodel->get_diso();
  d_ani   =ionicmodel->get_dani();
  chi     =ionicmodel->get_chi();
  C_m     =ionicmodel->get_C_m();

  nLocBas = in_locbas;
  
  dof_per_node = 1;
  
  vec_size = nLocBas * dof_per_node;
  
  nqp = in_nqp;

  Tangent  = new PetscScalar[vec_size * vec_size];
  Residual = new PetscScalar[vec_size];

  for(int ii=0; ii<vec_size; ++ii)
    Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size * vec_size; ++ii)
    Tangent[ii] = 0.0;

  R     = new double [nLocBas];
  dR_dx = new double [nLocBas];
  dR_dy = new double [nLocBas];
}

PLocAssem_NLHeat_2D_GenAlpha::~PLocAssem_NLHeat_2D_GenAlpha()
{
  delete [] Tangent;
  delete [] Residual;
  delete [] R;
  delete [] dR_dx;
  delete [] dR_dy;
}

void PLocAssem_NLHeat_2D_GenAlpha::Zero_Tangent_Residual()
{
  for(int ii=0; ii<vec_size; ++ii)
    Residual[ii] = 0.0;
  for(int ii=0; ii<vec_size * vec_size; ++ii)
    Tangent[ii] = 0.0;
}


void PLocAssem_NLHeat_2D_GenAlpha::Zero_Residual()
{
  for(int ii=0; ii<vec_size; ++ii)
    Residual[ii] = 0.0;
}


void PLocAssem_NLHeat_2D_GenAlpha::Assem_Estimate()
{
  for(int ii=0; ii<vec_size * vec_size; ++ii)
    Tangent[ii] = 1.0;
}

void PLocAssem_NLHeat_2D_GenAlpha::Assem_Residual(
    double time, double dt,
    const double * const &velo,
    const double * const &disp,
    const double * const &Iion,
    const double * const &dPhi_Iion,    
    const class FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const class AInt_Weight * const &weight )
{
  double v, d, d_x, d_y, Im, Istim;
  double gwts; // quadrature weights
  int qua; // quadrature index
  double coor_x, coor_y; // xy coor of the quad pts
  double  k11, k12, k21, k22; // material
  std::vector<double> E_Im(nLocBas);// E_ is for element
  std::vector<double> dIm(nLocBas);
  std::vector<double> dIstim(nLocBas);
  std::vector<double> E_dIm(nLocBas);  
  double curr = time + alpha_f * dt;

  int ii; // iterator

  for (ii=0; ii<nLocBas; ++ii)
    {
      E_Im[ii] = 0 - chi * Iion[ii]; //E_Istim=0
      E_dIm[ii] = 0 - chi * dPhi_Iion[ii];//E_dIstim=0
    } //is this ok for parallel proc

  Zero_Residual(); // zero all values for assembly

  for(qua=0; qua<nqp; ++qua)
  {
    v = 0.0; d = 0.0; Im= 0.0; d_x = 0.0; d_y = 0.0;
    coor_x = 0.0; coor_y = 0.0;
    
    element->get_R_gradR(qua, R, dR_dx, dR_dy);

    // Inner product with basis functions
    // ! This loop may be speed up by calling the sdot function in cblas 
    for(ii=0; ii<nLocBas; ++ii)
    {
      v   += velo[ii] * R[ii];
      d   += disp[ii] * R[ii];
      Im  += E_Im[ii] * R[ii];     //NOTE!: update I_m with Istim below
      dIm[ii] = E_dIm[ii] * R[ii]; //d_Im is a vector
      
      d_x += disp[ii] * dR_dx[ii];
      d_y += disp[ii] * dR_dy[ii];
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
    }
    
    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    //get stimulus for the volume (not the boundary)
    Istim = get_f(coor_x, coor_y, curr);
    //dIstim= 0; already
    
    Im = Im + Istim;
    for(ii=0; ii<nLocBas; ++ii)
      {      
	dIm[ii] = dIm[ii] +  dIstim[ii] ;
      }

    get_k(d, coor_x, coor_y, k11, k12, k21, k22);
    
    for(ii=0; ii<nLocBas; ++ii)
    {
      Residual[ii] += gwts * (chi * C_m * R[ii] * v
          + k11 * d_x * dR_dx[ii] + k12 * d_y * dR_dx[ii]
          + k21 * d_x * dR_dy[ii] + k22 * d_y * dR_dy[ii]
          - Im * R[ii] );
    }
  }
}

void PLocAssem_NLHeat_2D_GenAlpha::Assem_Tangent_Residual(
    double time, double dt,
    const double * const &velo,
    const double * const &disp,
    const double * const &Iion,
    const double * const &dPhi_Iion,    
    const class FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const class AInt_Weight * const &weight )
{
  int ii, qua, A, B, index; // iterator
  double v, d, d_x, d_y, Im, Istim;
  double gwts;
  double coor_x, coor_y;
  double k11, k12, k21, k22;
  double dk11, dk12, dk21, dk22;
  std::vector<double> E_Im(nLocBas);// E_ is for element
  std::vector<double> dIm(nLocBas);
  std::vector<double> dIstim(nLocBas);  
  std::vector<double> E_dIm(nLocBas);  
  double curr = time + alpha_f * dt;

  //calculate Im and dIm from Iion&Istim and dIion&dIstim
  // at nodes! for the local element
  for (ii=0; ii<nLocBas; ++ii)
    {
      E_Im[ii] = 0 - chi * Iion[ii]; //E_Istim=0
      E_dIm[ii] = 0 - chi * dPhi_Iion[ii];//E_dIstim=0
    }


    
  Zero_Tangent_Residual(); // zero all values for assembly
  
  for(qua =0; qua<nqp; ++qua)
  {
    v = 0.0; d = 0.0;  d_x = 0.0; d_y = 0.0;
    coor_x = 0.0; coor_y = 0.0;
    Im= 0.0;

    element->get_R_gradR(qua, R, dR_dx, dR_dy);
   
    // ! scalar product
    for(ii=0; ii<nLocBas; ++ii)
    {
      v   += velo[ii] * R[ii];
      d   += disp[ii] * R[ii];
      Im  += E_Im[ii] * R[ii];     //NOTE!: update I_m with Istim below
      dIm[ii] = E_dIm[ii] * R[ii]; //d_Im is a vector
      
      d_x += disp[ii] * dR_dx[ii];
      d_y += disp[ii] * dR_dy[ii];
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
    }
  
    gwts = element->get_detJac(qua) * weight->get_weight(qua);

    //get stimulus for the volume (not the boundary)
    Istim = get_f(coor_x, coor_y, curr);
    //dIstim= 0; already
    
    Im = Im + Istim;
    for(ii=0; ii<nLocBas; ++ii)
      {      
	dIm[ii] = dIm[ii] +  dIstim[ii] ;
      }
    
    get_k(d, coor_x, coor_y, k11, k12, k21, k22);
    get_dk_du(d, coor_x, coor_y, dk11, dk12, dk21, dk22);

    for(A=0; A<nLocBas; ++A)
    {
      //Residual[A] += gwts * (R[A] * v
      //    + k11 * d_x * dR_dx[A] + k12 * d_y * dR_dx[A]
      //    + k21 * d_x * dR_dy[A] + k22 * d_y * dR_dy[A]
      //    - f * R[A] );
      Residual[A] += gwts * (chi * C_m * R[A] * v
      	     + k11 * d_x * dR_dx[A] + k12 * d_y * dR_dx[A]
      	     + k21 * d_x * dR_dy[A] + k22 * d_y * dR_dy[A]
      	     - Im * R[A] );

      for(B=0; B<nLocBas; ++B)
      {
        index = A * nLocBas + B;

        //Tangent[index] += gwts * (R[A] * R[B] * alpha_m);
	Tangent[index] += gwts * (chi * C_m * R[A] * R[B] * alpha_m);

        Tangent[index] += gwts * alpha_f * gamma * dt * 
          ( k11 * dR_dx[A] * dR_dx[B] + k12 * dR_dx[A] * dR_dy[B]
          + k21 * dR_dy[A] * dR_dx[B] + k22 * dR_dy[A] * dR_dy[B]
          );

        Tangent[index] += gwts * alpha_f * gamma * dt *
          ( dk11 * d_x * dR_dx[A] + dk12 * d_y * dR_dx[A]
          + dk21 * d_x * dR_dy[A] + dk22 * d_y * dR_dy[A]
           ) * R[B];

	Tangent[index] += gwts * alpha_f * gamma *dt *
	  ( - R[A] * dIm[B]);
      }
    }
  }
}

void PLocAssem_NLHeat_2D_GenAlpha::Assem_Mass(
    const class FEAElement * const &element,
    const class AInt_Weight * const &weight )
{
  int qua, A, B, index; // iterator
  double gwts;

  for(A=0; A<vec_size*vec_size; ++A)
    Tangent[A] = 0.0;

  for(qua = 0; qua<nqp; ++qua)
  {
    gwts = element->get_detJac(qua) * weight->get_weight(qua);
    element->get_R(qua, R);
    
    for(A =0; A<nLocBas; ++A)
    {
      for(B =0; B<nLocBas; ++B)
      {
        index = A * nLocBas + B;
        Tangent[index] += gwts * chi * C_m * R[A] * R[B];
      }
    }
  }
}

void PLocAssem_NLHeat_2D_GenAlpha::Assem_Mass_Residual(
    const double * const &disp,
    const double * const &Iion,
    const double * const &dPhi_Iion,    
    const class FEAElement * const &element,
    const double * const &eleCtrlPts_x,
    const double * const &eleCtrlPts_y,
    const double * const &eleCtrlPts_z,
    const class AInt_Weight * const &weight )
{
  int ii, qua, A, B, index; // iterator
  double d, d_x, d_y, Im, Istim;
  double gwts;
  double coor_x, coor_y;
  double k11, k12, k21, k22;
  double dk11, dk12, dk21, dk22;
  std::vector<double> E_Im(nLocBas);// E_ is for element
  std::vector<double> dIm(nLocBas);
  std::vector<double> dIstim(nLocBas);  
  std::vector<double> E_dIm(nLocBas);  

  double curr = 0.0;

  //calculate Im and dIm from Iion&Istim and dIion&dIstim
  // at nodes! for the local element
  for (ii=0; ii<nLocBas; ++ii)
    {
      E_Im[ii] = 0 - chi * Iion[ii]; //E_Istim=0
      E_dIm[ii] = 0 - chi * dPhi_Iion[ii];//E_dIstim=0
    }

  Zero_Tangent_Residual(); // zero all values for assembly
  
  for(qua =0; qua<nqp; ++qua)
  {
    d = 0.0; d_x = 0.0; d_y = 0.0;
    coor_x = 0.0; coor_y = 0.0;

    element->get_R_gradR(qua, R, dR_dx, dR_dy);
   
    // ! scalar product
    for(ii=0; ii<nLocBas; ++ii)
    {
      d   += disp[ii] * R[ii];
      Im  += E_Im[ii] * R[ii];//NOTE!: update I_m with Istim below
      dIm[ii] = E_dIm[ii] * R[ii]; //d_Im is a vector
     
      d_x += disp[ii] * dR_dx[ii];
      d_y += disp[ii] * dR_dy[ii];
      coor_x += eleCtrlPts_x[ii] * R[ii];
      coor_y += eleCtrlPts_y[ii] * R[ii];
    }

    gwts = element->get_detJac(qua) * weight->get_weight(qua);


    //get stimulus for the volume (not the boundary)
    Istim = get_f(coor_x, coor_y, curr);
    //dIstim= 0; already
    
    Im = Im + Istim;
    for(ii=0; ii<nLocBas; ++ii)
      {      
	dIm[ii] = dIm[ii] +  dIstim[ii] ;
      }
    
    get_k(d, coor_x, coor_y, k11, k12, k21, k22);
    get_dk_du(d, coor_x, coor_y, dk11, dk12, dk21, dk22);

    for(A=0; A<nLocBas; ++A)
    {
      Residual[A] -= gwts * (
            k11 * d_x * dR_dx[A] + k12 * d_y * dR_dx[A]
          + k21 * d_x * dR_dy[A] + k22 * d_y * dR_dy[A]
          - Im * R[A] );

      for(B=0; B<nLocBas; ++B)
      {
        index = A * nLocBas + B;
        Tangent[index] += gwts * chi * C_m * R[A] * R[B];
      }
    }
  }
}

//EOF
