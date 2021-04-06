#ifndef MESH_MIXED_HPP
#define MESH_MIXED_HPP
// ==================================================================
// Mesh_Mixed.hpp
//
// This is the instantiation of the IMesh class mixed element types
//
// Date:   2021
// ==================================================================
#include "IMesh.hpp"
#include "IIEN.hpp"
#include "Vec_Tools.hpp"

class Mesh_Mixed : public IMesh
{
public:
  Mesh_Mixed(const std::vector< IMesh * > mesh_list,
	     const std::vector< int > &elemType_list,
	     const IIEN * const &ien_ptr);
  //Mesh_Mixed(const int &in_nfunc, const int &in_nelem);

  virtual ~Mesh_Mixed();

  virtual void print_mesh_info() const;


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
  std::vector < std::vector < int > > stu_degrees;
};

#endif
