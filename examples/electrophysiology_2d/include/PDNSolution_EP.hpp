#ifndef PDNSOLUTION_EP_HPP
#define PDNSOLUTION_EP_HPP
// ==================================================================
// PDNSolution_EP.hpp
// This is the class that we define initial solution for EP
// ==================================================================
#include "Sys_Tools.hpp"
#include "PDNSolution.hpp"
#include "APart_Node.hpp"
#include "IALocal_BC.hpp"
#include "FEANode.hpp"

class PDNSolution_EP : public PDNSolution
{
  public:
    PDNSolution_EP( const class APart_Node * const &pNode,
			 const FEANode * const &fNode,
			 const class IALocal_BC * const &lbc,
			 int type );
  
    virtual ~PDNSolution_EP();

    // Initial solution setting
    // case 0: all zero for initial temperature
    void Init_ZeroTemp( const class IALocal_BC * const &LBC );

    // case 1: interior values 1.0, boundary values 0.0
    void Init_OneTemp( const class IALocal_BC * const &LBC );
  
    // case 2: values -80.0 and 0, partially
    void Init_Partial( const class APart_Node * const &pNode,
		       const FEANode * const &fNode,
		       const class IALocal_BC * const &LBC);
    // case 3: Whole domain rest at -80.0mV
    void Init_Rest( const class IALocal_BC * const &LBC );
    // get the global size of solution vector
    virtual int GetSize() const;
};
#endif
