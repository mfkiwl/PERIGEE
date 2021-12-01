#include "Mesh_Mixed.hpp"

Mesh_Mixed::Mesh_Mixed(const std::vector< IMesh * > &mesh_list,
		       const std::vector< int > &elemType_list,
		       const IIEN * const &ien_ptr,
		       const std::vector< std::vector<double> > &myo_fiber,
		       const std::vector< int > &phy_tag_list)
{
  nFunc = ien_ptr->get_nFunc_tot();
  nElem = ien_ptr->get_nElem_tot();
  nElemXnLocBas = 0;
  nLocBas.clear();
  elemType.clear();
  phy_tag.clear();
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
    phy_tag.insert(phy_tag.end(),
		   (*it)->get_nElem(),
		   phy_tag_list.at(distance));
    stu_degrees.insert( stu_degrees.end(), (*it)->get_nElem(),
			std::vector<int> {(*it)->get_s_degree(),
			    		  (*it)->get_t_degree(),
			                  (*it)->get_u_degree() });

    if(elemType_list.at(distance) == 501){
      fiber_ori.insert( fiber_ori.end(), myo_fiber.begin(), myo_fiber.end());
    }else if (elemType_list.at(distance) == 512){
      fiber_ori.insert( fiber_ori.end(), (*it)->get_nElem(),
			std::vector<double> {0.0, 0.0, 0.0} );
    }else{
      SYS_T::print_exit("Mesh Mixer constructor: element type is not implemented. \n");
    }
  }
}

//Mesh_Mixed::Mesh_Mixed(const int &in_nfunc, const int &in_nelem)
//  : nFunc(in_nfunc), nElem(in_nelem)
//{
//}


Mesh_Mixed::~Mesh_Mixed()
{}

void Mesh_Mixed::get_fiber_ori_loc(std::vector<std::vector<double>> &fiber_ori_loc,
				   const std::vector<int> &elem_loc) const
{
  fiber_ori_loc.clear();
  for (size_t i = 0; i < elem_loc.size(); ++i ) {
    fiber_ori_loc.push_back(fiber_ori.at(elem_loc[i]));
  }
}

void Mesh_Mixed::get_phy_tag_loc(std::vector<int> &phy_tag_loc,
				 const std::vector<int> &elem_loc) const
{
  phy_tag_loc.clear();
  for (size_t i = 0; i < elem_loc.size(); ++i ) {
    phy_tag_loc.push_back(phy_tag.at(elem_loc[i]));
  }
}

  
void Mesh_Mixed::print_mesh_info() const
{
  std::cout<<'\n';
  std::cout<<"======= Mesh_Mixed ======="<<std::endl;
  //std::cout<<"Degree S-T-U: "<<std::endl;
  //for( auto it=stu_degrees.begin(); it!=stu_degrees.end(); ++it ){
  //  VEC_T::print(*it);
  //}
  std::cout<<"Total Elem: "<<get_nElem()<<std::endl;
  std::cout<<"Total Func: "<<get_nFunc()<<std::endl;
  std::cout<<"Local Basis times element: "<<get_nElemXnLocBas()<<std::endl;
  //std::cout<<"Local Basis per element: "<<std::endl;
  //VEC_T::print(nLocBas);
  //std::cout<<"Element type per element: "<<std::endl;
  //VEC_T::print(elemType);
  //std::cout<<"Phy tag per element: "<<std::endl;
  //VEC_T::print(phy_tag);
  std::cout<<"========================="<<std::endl;
}

// EOF
