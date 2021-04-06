#include "IEN_Mixed.hpp"


IEN_Mixed::IEN_Mixed(const std::vector< std::vector<int> > &IEN_list,
		     const std::vector< IMesh * >  &mesh_list,
		     //const std::vector< int > &elemType_list,
		     const std::vector< std::vector<double> > &ctrlPts_list,
		     const std::string &endnodes_filename,
		     //std::vector< int > &elemType_combined,
		     std::vector<double> &ctrlPts_combined)
		     //int &nFunc_tot,
		     //int &nElem_tot)
{
  if(mesh_list.size() != IEN_list.size())
    {
      std::cerr<<"ERROR: IEN_list and mesh_list sizes don't match. \n";
      exit(1);
    }

  auto vecIEN_1= IEN_list.at(0);
  auto vecIEN_2= IEN_list.at(1);
  IMesh *Mesh_1  = mesh_list.at(0);
  IMesh *Mesh_2  = mesh_list.at(1);
  int nnode1= Mesh_1->get_nFunc();
  int nnode2= Mesh_2->get_nFunc();
  int nelem1= Mesh_1->get_nElem();
  int nelem2= Mesh_2->get_nElem();
  int nLocBas1= Mesh_1->get_nLocBas();
  int nLocBas2= Mesh_2->get_nLocBas();
  
  //1- read endnodes file contents
  std::vector<int> endnodes; 
  std::ifstream endnode_file(endnodes_filename);
  if (!endnode_file)
    {
      // Print an error and exit
      std::cerr << endnodes_filename
		<< "could not be opened for reading!" << std::endl;
      exit(1);
    }
  int number ;
  while (endnode_file>>number) {
    endnodes.push_back(number);
  }
  VEC_T::sort_unique_resize(endnodes);
  //  std::cout<< "endnodes:" <<std::endl;
  //VEC_T::print(endnodes);

  //2-Now find the closest nodes on Mesh1 to the nodes on Mesh2
  std::vector<double> ctrlPts_1= ctrlPts_list.at(0);
  std::vector<double> ctrlPts_2= ctrlPts_list.at(1);
  std::vector<int> node1_list;
  //node1_list.resize(endnodes.size());
  //std::cout<< "node1s list size:" << node1_list.size() <<std::endl;

  double distance, x1, x2, y1, y2, z1, z2;
  

  //std::cout<< "nnode1:" << nnode1 <<std::endl;
  //std::cout<< "nnode2:" << nnode2 <<std::endl;
  int node1, node2;
  
  for (auto it2=endnodes.begin(); it2!=endnodes.end(); ++it2) {

    node2=*it2;
    x2=ctrlPts_2.at(3*node2);
    y2=ctrlPts_2.at(3*node2+1);
    z2=ctrlPts_2.at(3*node2+2);
    std::cout<< "node 2 coords" <<std::endl;
    std::cout<< x2 <<","
    	     << y2 <<","
    	     << z2 <<std::endl;

    for (int ii=0; ii<nnode1; ii++) {
      x1=ctrlPts_1.at(3*ii+0);
      y1=ctrlPts_1.at(3*ii+1);
      z1=ctrlPts_1.at(3*ii+2);

      std::cout<< "node 1 coords" <<std::endl;
      std::cout<< x1 <<","
	       << y1 <<","
	       << z1 <<std::endl;
      distance = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
      std::cout<< "distance:" << distance <<std::endl;
      if (distance < 1e-2) {
      	node1= ii;
      	break; 
      } else if((ii+1)==nnode1) { //is this else robust? 
	std::cerr<<"ERROR: couldn't find a matching node for purkinje node: "
		 << node2 << " \n";
	exit(1);
      }
    }
    node1_list.push_back(node1);
    std::cout<< "node1-node2 couple is:" << node1<< "," << node2 <<std::endl;
  }
  std::cout<< "node1-list:" <<std::endl;
  VEC_T::print(node1_list);
  
  //3- find locations of nodes in the IEN2 
  std::vector<int>::iterator location;
  std::vector< std::vector<int> > node2_locations;
  node2_locations.resize(nnode2);//resize the 1st dimension only

  for (int ii=0; ii<nnode2 ; ii++) {
    location= std::find(vecIEN_2.begin(), vecIEN_2.end(), ii);
    
    while (location!=vecIEN_2.end()){
      (node2_locations.at(ii)).push_back(*location);
      location= std::find(location+1, vecIEN_2.end(), ii);
    }
  }

  std::cout<< "node2-locations" << std::endl;
  for( auto it = node2_locations.begin(); it != node2_locations.end(); it++ ){
    for( auto it2 = (*it).begin(); it2 != (*it).end(); it2++ ){
      std::cout<< *it2 <<std::endl;
    }
  }

  //4- subtract nodes of endnodes from IEN2
  //   and keep the rest of the nodes separate
  std::vector<int> endnodes_reverse;
  endnodes_reverse.resize(endnodes.size());
  std::reverse_copy(endnodes.begin(), endnodes.end(), endnodes_reverse.begin());
  
  std::vector< std::vector<int> > locs2_replace, locs2_reorder;
  locs2_replace.resize(endnodes.size());
  //locs2_reorder.resize(nnode2-endnodes.size());
  locs2_reorder=node2_locations;
  
  int ii=0;
  for( auto it = endnodes.begin(); it != endnodes.end(); it++ ){
    //(locs2_replace.at(ii)).resize((node2_locations.at(*it)).size());
    locs2_replace.at(ii)= node2_locations.at(*it);
    ii++;
  }
  for( auto it = endnodes_reverse.begin(); it != endnodes_reverse.end(); it++ ){
    locs2_reorder.erase(locs2_reorder.begin() + (*it));
    ctrlPts_2.erase(ctrlPts_2.begin()+(*it), ctrlPts_2.begin()+(*it)+3);
  }
  
  std::cout<< "locs 2 replace" << std::endl;
  for( auto it = locs2_replace.begin(); it != locs2_replace.end(); it++ ){
    for( auto it2 = (*it).begin(); it2 != (*it).end(); it2++ ){
      std::cout<< *it2 <<std::endl;
    }
  }
  std::cout<< "locs 2 reorder" << std::endl;
  for( auto it = locs2_reorder.begin(); it != locs2_reorder.end(); it++ ){
    for( auto it2 = (*it).begin(); it2 != (*it).end(); it2++ ){
      std::cout<< *it2 <<std::endl;
    }
  }

  //5- merge the numbering of endnodes of IEN2 to IEN1
  //   and renumber rest of the nodes
  ii=0;
  for( auto it = locs2_replace.begin(); it != locs2_replace.end(); it++ ){
    for( auto it2 = (*it).begin(); it2 != (*it).end(); it2++ ){
      vecIEN_2.at(*it2) = node1_list.at(ii);
    }
    ii++;
  }

  ii=nnode1;
  for( auto it = locs2_reorder.begin(); it != locs2_reorder.end(); it++ ){
    for( auto it2 = (*it).begin(); it2 != (*it).end(); it2++ ){
      vecIEN_2.at(*it2) = ii;
    }
    ii++;
  }
  nFunc_tot= ii;

  //6- store IENs as 2d vectors
  std::vector< std::vector <int> > IEN1, IEN2;
  IEN1.resize(nelem1);  IEN2.resize(nelem2);
  for( auto it = IEN1.begin(); it != IEN1.end(); it++ ){
    (*it).resize(nLocBas1);
  }
  for (int jj=0; jj<vecIEN_1.size() ; jj++){
    int row= jj/nLocBas1;
    int col= jj%nLocBas1;
    (IEN1.at(row)).at(col)=vecIEN_1.at(jj);
  }
     
  for( auto it = IEN2.begin(); it != IEN2.end(); it++ ){
    (*it).resize(nLocBas2);
  }
  for (int jj=0; jj<vecIEN_2.size() ; jj++){
    int row= jj/nLocBas2;
    int col= jj%nLocBas2;
    (IEN2.at(row)).at(col)=vecIEN_2.at(jj);
  }
  
  IEN1.insert(IEN1.end(), IEN2.begin(), IEN2.end());


  //copy IEN1 to IEN
  IEN=IEN1;
  std::cout<< "ien  " << std::endl;
  for( auto it = IEN.begin(); it != IEN.end(); it++ ){
    VEC_T::print(*it);
  }

  //assign output parameters
  nElem_tot = IEN.size();
  //nFunc = ii  assigned above

  ctrlPts_combined=ctrlPts_1;
  ctrlPts_combined.insert(ctrlPts_combined.end(),
			  ctrlPts_2.begin(),ctrlPts_2.end());

  //elemType_combined=std::vector<int> (nelem1, elemType_list.at(0));
  //elemType_combined.insert(elemType_combined.end(),
  //                         nelem2, elemType_list.at(1));
  

  std::cout << "nElem_tot= " << nElem_tot << std::endl;
  std::cout << "nFunc_tot= " << nFunc_tot << std::endl;

}


IEN_Mixed::~IEN_Mixed()
{

}


int IEN_Mixed::get_IEN( const int &ee, const int &l_node ) const
{
  return (IEN.at(ee)).at(l_node);
}


void IEN_Mixed::print_IEN() const
{
  std::cout<<std::endl;
  std::cout<<"====== IEN ====== \n";
  for( auto it = IEN.begin(); it != IEN.end(); it++ ){
    VEC_T::print(*it);
  }
  std::cout<<"================= \n";
}

// EOF
