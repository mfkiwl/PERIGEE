#ifndef FEAELEMENT_LINE2_3D_DER1_HPP
#define FEAELEMENT_LINE2_3D_DER1_HPP
// ==================================================================
// FEAElement_Line2_3D_der1.hpp
//
// This is an implementation of the element routine for 1D P1 
// element with basis function and its gradients  evaluated.
// This routine is modified from FEAElement_Line2_3D_der0 which 
// is only for boundary conditions
// Note: this element returns only the x-component of gradient
//       because this is a 1-d elmeent 
//
// Node index :  0 --------------- 1
//
// Author: Ju Liu
// modified: Oguz Ziya Tikenogullari
// Date created: Feb 2021
// ==================================================================
#include "FEAElement.hpp"

class FEAElement_Line2_3D_der1 : public FEAElement
{
  public:
    FEAElement_Line2_3D_der1( const int &in_nqua );

    virtual ~FEAElement_Line2_3D_der1();

    virtual int get_elemDim() const {return 1;}

    virtual int get_Type() const {return 512;}

    virtual int get_numType() const {return 1;}

    virtual int get_numQuapts() const {return numQuapts;}

    virtual int get_nLocBas() const {return nLocBas;}

    virtual void print() const;

    virtual double get_memory_usage() const;

    virtual void buildBasis( const IQuadPts * const &quad_rule,
        const double * const &ctrl_x, const double * const &ctrl_y,
        const double * const &ctrl_z );

    virtual void get_R( const int &quaindex, double * const &basis ) const;

    virtual void get_R_gradR( const int &quaindex, 
			      double * const &basis,
			      double * const &basis_x) const;
  
    virtual void get_normal_out( const int &quaindex,
        const double * const &ctrl_x, const double * const &ctrl_y,
        const double * const &ctrl_z,
        const double &intpt_x, const double &intpt_y, const double &intpt_z,
        double &nx, double &ny, double &nz, double &area ) const;

    virtual double get_detJac(const int &quaindex) const 
    {return detJac;};

  private:
    const int nLocBas, numQuapts;

    // length nLocBas x numQuapts = 2 x numQuapts
    double * R;

    double dR_dx[2];
  
    // length should be numQuapts, 
    // here for linear element, they are all constant
    double dx_dr, dy_dr, dz_dr, detJac;
};

#endif
