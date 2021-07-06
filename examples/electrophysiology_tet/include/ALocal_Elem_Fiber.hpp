#ifndef ALOCAL_ELEM_FIBER_HPP
#define ALOCAL_ELEM_FIBER_HPP
// ==================================================================
// ALocal_Elem_Fiber.hpp
// Analysis use Local Element class. This class records the local 
// partition domain's element index and total number of local elements.
//
// Author: Ju Liu
// Date: Nov. 10 2013
// ==================================================================
#include "HDF5_Reader.hpp"
#include "ALocal_Elem.hpp"

class ALocal_Elem_Fiber : public ALocal_Elem
{
public:
  // Constructor : read h5 file by giving the part file base name and rank
  ALocal_Elem_Fiber(const std::string &fileBaseName, const int &cpu_rank);

  virtual ~ALocal_Elem_Fiber();

  virtual void get_fiber_ori_e( std::vector<double> &fiber_ori_e, const int &ee)
    const {fiber_ori_e = fiber_ori_loc.at(ee);}
  //    // Assess the data
  //    virtual int get_elem_loc(const int &index) const {return elem_loc[index];}
  //    
  //    virtual int get_nlocalele() const {return nlocalele;}
  //
  //    virtual void print_info() const;
  //
  //    // This is a virtual function for multiphysics simulations. A tag
  //    // is attached to an element to identify different physical domain,
  //    // for example, fluid subdomain and solid subdomain.
  //    // For single domain problem, this function is NOT needed.
  //    virtual int get_elem_tag(const int &index) const
  //    {
  //      SYS_T::print_fatal("Error: ALocal_Elem_Fiber::get_elem_tag is not implemented.\n");
  //      return -1;
  //    }

private:
  // A vector of vector storing the local partition's fiber orientations
  //for each element.
  std::vector<std::vector<double>> fiber_ori_loc;

};

#endif
