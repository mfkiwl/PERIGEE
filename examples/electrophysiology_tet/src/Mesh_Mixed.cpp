#include "Mesh_Mixed.hpp"

Mesh_Mixed::Mesh_Mixed(const std::vector< IMesh * > mesh_list,
		       const std::vector< int > &elemType_list,
		       const IIEN * const &ien_ptr)
{
  nFunc = ien_ptr->get_nFunc_tot();
  nElem = ien_ptr->get_nElem_tot();
  nElemXnLocBas = 0;
  nLocBas.clear();
  elemType.clear();
  stu_degrees.clear();
  int distance; 
  
  for (auto it = mesh_list.begin(); it != mesh_list.end(); ++it){
    distance= std::distance(mesh_list.begin(), it);
    
    nElemXnLocBas= nElemXnLocBas
      + ((*it)->get_nElem())*((*it)->get_nLocBas());
    nLocBas.insert(nLocBas.end(),
		   (*it)->get_nElem(), (*it)->get_nLocBas());
    elemType.insert(elemType.end(),
		    (*it)->get_nElem(),
		    elemType_list.at(distance));
    stu_degrees.insert( stu_degrees.end(), (*it)->get_nElem(),
			std::vector<int> {(*it)->get_s_degree(),
			    		  (*it)->get_t_degree(),
			                  (*it)->get_u_degree() });
  }
}

//Mesh_Mixed::Mesh_Mixed(const int &in_nfunc, const int &in_nelem)
//  : nFunc(in_nfunc), nElem(in_nelem)
//{
//}


Mesh_Mixed::~Mesh_Mixed()
{}

  
void Mesh_Mixed::print_mesh_info() const
{
  std::cout<<'\n';
  std::cout<<"======= Mesh_Mixed ======="<<std::endl;
  std::cout<<"Degree S-T-U: "<<std::endl;
  for( auto it=stu_degrees.begin(); it!=stu_degrees.end(); ++it ){
    VEC_T::print(*it);
  }
  std::cout<<"Total Elem: "<<get_nElem()<<std::endl;
  std::cout<<"Total Func: "<<get_nFunc()<<std::endl;
  std::cout<<"Local Basis times element: "<<get_nElemXnLocBas()<<std::endl;
  std::cout<<"Local Basis per element: "<<std::endl;
  VEC_T::print(nLocBas);
  std::cout<<"Element type per element: "<<std::endl;
  VEC_T::print(elemType);
  std::cout<<"========================="<<std::endl;
}

// EOF
