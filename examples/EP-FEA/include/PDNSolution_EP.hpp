#ifndef PDNSOLUTION_EP_HPP
#define PDNSOLUTION_EP_HPP
// ==================================================================
// PDNSolution_EP.hpp
// This is the class that we define initial solution for EP
// Author: Ju Liu, liujuy@gmail.com
// Modified: Oguz Ziya Tikenogullari, o.z.tikenogullari@gmail.com
// ==================================================================
#include "Sys_Tools.hpp"
#include "PDNSolution.hpp"
#include "APart_Node.hpp"
#include "ALocal_NodalBC.hpp"
#include "FEANode.hpp"
#include "ALocal_IEN_Mixed.hpp"
#include "IonicModel.hpp"
#include <petscvec.h>
#include <petscsys.h>
#include <petscis.h>
#include <petscviewer.h>

class PDNSolution_EP : public PDNSolution
{
public:
  PDNSolution_EP( const class APart_Node * const &pNode,
		  const FEANode * const &fNode,
		  const class ALocal_NodalBC * const &lbc,
		  int type );

  //this constructor is for initializing EP internal variables.
  // it asks the ionicmodel class to set the internal variable
  // at each node. so this is suitable only with mass lumping
  // used with with reaction/ionic part of the EP problem.
  PDNSolution_EP(const class APart_Node * const &pNode,
		 const int &input_dof_num,
		 const FEANode * const &fNode,
		 const class ALocal_NodalBC * const &locbc ,
		 const int type,
		 const ALocal_IEN_Mixed * const &lien_ptr,
		 const std::vector< IonicModel * > &ionicmodel_array);
  
  virtual ~PDNSolution_EP();

  // Initial solution setting
  // case 0: all zero for initial temperature
  void Init_ZeroTemp( const class ALocal_NodalBC * const &LBC );

  // case 1: interior values 1.0, boundary values 0.0
  void Init_OneTemp( const class ALocal_NodalBC * const &LBC );
  
  // case 2: values -80.0 and 0, partially
  void Init_Partial( const class APart_Node * const &pNode,
		     const FEANode * const &fNode,
		     const class ALocal_NodalBC * const &LBC);
  // case 3: Whole domain rest at -80.0mV
  void Init_Rest( const class ALocal_NodalBC * const &LBC );
  // get the global size of solution vector
  virtual int GetSize() const;
};
#endif
