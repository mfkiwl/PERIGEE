#include "Mesh_Line_3D.hpp"

Mesh_Line_3D::Mesh_Line_3D(const int &in_nfunc, const int &in_nelem)
  : nFunc(in_nfunc), nElem(in_nelem)
{}


Mesh_Line_3D::~Mesh_Line_3D()
{}


void Mesh_Line_3D::print_mesh_info() const
{
  std::cout<<'\n';
  std::cout<<"======= Mesh_Line_3D ======="<<std::endl;
  std::cout<<"Degree: "<<get_s_degree()<<'\t'
	   <<"and no t and u degrees"<<std::endl;
  std::cout<<"Total Elem: "<<get_nElem()<<std::endl;
  std::cout<<"Total Func: "<<get_nFunc()<<std::endl;
  std::cout<<"Local Basis #: "<<get_nLocBas()<<std::endl;
  std::cout<<"========================="<<std::endl;
}

// EOF
