#ifndef MESH_MIXED_HPP
#define MESH_MIXED_HPP
// ==================================================================
// Mesh_Mixed.hpp
//
// This is the instantiation of the IMesh class for mixed element types
// we assume element types of linear tet(501) and linear line (512) :
//               search for 501 and 512 in cpp file.
//
// Date:   2021
// ==================================================================
#include "IMesh.hpp"
#include "IIEN.hpp"
#include "Vec_Tools.hpp"

class Mesh_Mixed : public IMesh
{
public:
  Mesh_Mixed(const std::vector< IMesh * > &mesh_list,
	     const std::vector< int > &elemType_list,
	     const IIEN * const &ien_ptr,
	     const std::vector< std::vector<double> > &myo_fiber,
	     const std::vector< int > &phy_tag_list);

  virtual ~Mesh_Mixed();

  virtual void print_mesh_info() const;
  
  virtual void get_fiber_ori_loc(std::vector< std::vector< double > > &fiber_ori_loc,
				 const std::vector<int> &elem_loc) const ;
  
  virtual void get_phy_tag_loc( std::vector< int > &phy_tag_loc,
				const std::vector<int> &elem_loc) const ;

  virtual int get_s_degree()
    const {SYS_T::print_exit("Error: Mesh.get_s_degree is not implemented. \n");return 0;}
  virtual int get_t_degree()							
    const {SYS_T::print_exit("Error: Mesh.get_s_degree is not implemented. \n");return 0;}
  virtual int get_u_degree()							
    const {SYS_T::print_exit("Error: Mesh.get_s_degree is not implemented. \n");return 0;}
  
  virtual void get_stu_degree( std::vector <int >  &stu_of_elem, const int &ee) const { stu_of_elem = stu_degrees.at(ee);}

  virtual void get_stu_deg_vec(std::vector< std::vector< int > > &stu_deg_vec_in) const {stu_deg_vec_in = stu_degrees;}

  virtual void get_nLocBas_vec(std::vector< int >  &nLocBas_vec_in) const {nLocBas_vec_in = nLocBas;}
  
  virtual void get_elemType_vec(std::vector< int >  &elemType_vec_in) const {elemType_vec_in = elemType;}
  
  virtual int get_nFunc() const {return nFunc;}

  virtual int get_nElem() const {return nElem;}

  virtual int get_nElemXnLocBas() const {return nElemXnLocBas;}

  virtual int get_nLocBas() 
    const {SYS_T::print_exit("Error: you need to provide element number because mesh is mixed with element types. \n");return 0;}

  virtual int get_nLocBas( const int &ee) const {return nLocBas.at(ee);}

  virtual int get_elemType( const int &ee) const {return elemType.at(ee);}
  
private:
  
  int nFunc, nElem, nElemXnLocBas;
  std::vector < int > nLocBas; // nlocbas per element
  std::vector < int > elemType; //elemTypes per element
  std::vector < int > phy_tag; //phy_tag per element
  std::vector < std::vector < int > > stu_degrees;
  std::vector < std::vector < double > > fiber_ori;
};

#endif
