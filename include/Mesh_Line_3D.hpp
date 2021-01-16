#ifndef MESH_LINE_3D_HPP
#define MESH_LINE_3D_HPP
// ==================================================================
// Mesh_Line_3D.hpp
//
// This is the instantiation of the IMesh class for 2-node line
// element,  unstructured mesh.
//
// Date: Jan 12  2021
// ==================================================================
#include "IMesh.hpp"

class Mesh_Line_3D : public IMesh
{
  public:
    Mesh_Line_3D(const int &in_nFunc, const int &in_nElem);

    virtual ~Mesh_Line_3D();

    virtual void print_mesh_info() const;

    virtual int get_s_degree() const {return 1;}
    virtual int get_t_degree() const {SYS_T::print_exit("Error: line element has no shape functions in t-direction. \n"); return 0;}
    virtual int get_u_degree() const {SYS_T::print_exit("Error: line element has no shape functions in u-direction. \n"); return 0;}

    virtual int get_nFunc() const {return nFunc;}

    virtual int get_nElem() const {return nElem;}

    virtual int get_nLocBas() const {return 2;}

  private:
    const int nFunc, nElem;
};

#endif
