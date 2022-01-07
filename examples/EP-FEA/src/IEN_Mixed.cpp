#include "IEN_Mixed.hpp"


IEN_Mixed::IEN_Mixed(const std::vector< std::vector<int> > &IEN_list,
		     const std::vector< IMesh * >  &mesh_list,
		     const std::vector< std::vector<double> > &ctrlPts_list,
		     const std::vector< std::vector< double > > &displacement_list,
		     const std::string &LVendnodes_filename,
		     const std::string &RVendnodes_filename,
		     std::vector<double> &ctrlPts_combined,
		     const double &LV_tol,
		     const double &RV_tol)
{
  if(mesh_list.size() != IEN_list.size())    {
    std::cerr<<"ERROR: IEN_list and mesh_list sizes don't match. \n";
    exit(1);
  }

  auto vecIEN_1= IEN_list.at(0);
  auto vecIEN_2= IEN_list.at(1);
  auto vecIEN_3= IEN_list.at(2);
  IMesh *Mesh_1  = mesh_list.at(0);
  IMesh *Mesh_2  = mesh_list.at(1);
  IMesh *Mesh_3  = mesh_list.at(2);
  int nnode1= Mesh_1->get_nFunc();
  int nnode2= Mesh_2->get_nFunc();
  int nnode3= Mesh_3->get_nFunc();
  int nelem1= Mesh_1->get_nElem();
  int nelem2= Mesh_2->get_nElem();
  int nelem3= Mesh_3->get_nElem();
  int nLocBas1= Mesh_1->get_nLocBas();
  int nLocBas2= Mesh_2->get_nLocBas();
  int nLocBas3= Mesh_3->get_nLocBas();
  std::vector<double> displacement_1= displacement_list.at(0);
  std::vector<double> displacement_2= displacement_list.at(1);
  std::vector<double> displacement_3= displacement_list.at(2);

  
  //1- read endnodes file contents
  std::vector<int> LVendnodes; 
  std::ifstream LVendnode_file(LVendnodes_filename);
  if (!LVendnode_file)
    {
      // Print an error and exit
      std::cerr << LVendnodes_filename
		<< "LV endnodes could not be opened for reading!" << std::endl;
      exit(1);
    }
  int number2 ;
  while (LVendnode_file>>number2) {
    LVendnodes.push_back(number2);
  }
  VEC_T::sort_unique_resize(LVendnodes);
  
  std::vector<int> RVendnodes; 
  std::ifstream RVendnode_file(RVendnodes_filename);
  if (!RVendnode_file)
    {
      // Print an error and exit
      std::cerr << RVendnodes_filename
		<< "RV endnodes could not be opened for reading!" << std::endl;
      exit(1);
    }
  int number3 ;
  while (RVendnode_file>>number3) {
    RVendnodes.push_back(number3);
  }
  VEC_T::sort_unique_resize(RVendnodes);

  //2-Now find the closest nodes on Mesh1 to the nodes on Mesh2
  std::vector<double> ctrlPts_1= ctrlPts_list.at(0);
  std::vector<double> ctrlPts_2= ctrlPts_list.at(1);
  std::vector<double> ctrlPts_3= ctrlPts_list.at(2);
  std::vector<int> node12_list, node13_list;
  //node1_list.resize(endnodes.size());
  //std::cout<< "node1s list size:" << node1_list.size() <<std::endl;

  double distance, x1, x2, x3, y1, y2, y3, z1, z2, z3;
  
  //std::cout<< "nnode1:" << nnode1 <<std::endl;
  //std::cout<< "nnode2:" << nnode2 <<std::endl;
  //std::cout<< "nnode3:" << nnode3 <<std::endl;
  int node12, node2, node13, node3;

  //Match myocardium with LV purkinje 
  for (auto it2=LVendnodes.begin(); it2!=LVendnodes.end(); ++it2) {

    node2=*it2;
    x2=ctrlPts_2.at(3*node2);
    y2=ctrlPts_2.at(3*node2+1);
    z2=ctrlPts_2.at(3*node2+2);
    //std::cout<< "node 2 coords" <<std::endl;
    //std::cout<< x2 <<","
    //	     << y2 <<","
    //	     << z2 <<std::endl;

    for (int ii=0; ii<nnode1; ii++) {
      x1=ctrlPts_1.at(3*ii+0);
      y1=ctrlPts_1.at(3*ii+1);
      z1=ctrlPts_1.at(3*ii+2);

      //std::cout<< "node 1 coords" <<std::endl;
      //std::cout<< x1 <<","
      //	       << y1 <<","
      //	       << z1 <<std::endl;
      distance = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
      //std::cout<< "distance:" << distance <<std::endl;
      if (distance < LV_tol) { // turn this tolerance value into a user param
      	node12= ii;
	//std::cout<< "node 1 coords" <<std::endl;
	//std::cout<< x1 <<","
	//	 << y1 <<","
	//	 << z1 <<std::endl;
      	break; 
      } else if((ii+1)==nnode1) { //is this else robust? 
	std::cerr<<"ERROR: couldn't find an LV matching node for purkinje node: "
		 << node2 << " \n";
	exit(1);
      }
    }
    node12_list.push_back(node12);
    //std::cout<< "node12-node2 couple is:" << node12<< "," << node2 
    //     << "(dist :" << distance << ")" <<std::endl;
  }

  //Match myocardium with RV purkinje 
  for (auto it2=RVendnodes.begin(); it2!=RVendnodes.end(); ++it2) {
    node3=*it2;
    x3=ctrlPts_3.at(3*node3);
    y3=ctrlPts_3.at(3*node3+1);
    z3=ctrlPts_3.at(3*node3+2);
    //std::cout<< "node 3 coords" <<std::endl;
    //std::cout<< x3 <<","
    //	     << y3 <<","
    //	     << z3 <<std::endl;

    for (int ii=0; ii<nnode1; ii++) {
      x1=ctrlPts_1.at(3*ii+0);
      y1=ctrlPts_1.at(3*ii+1);
      z1=ctrlPts_1.at(3*ii+2);

      //std::cout<< "node 1 coords" <<std::endl;
      //std::cout<< x1 <<","
      //	       << y1 <<","
      //	       << z1 <<std::endl;
      distance = sqrt((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3)+(z1-z3)*(z1-z3));
      //std::cout<< "distance:" << distance <<std::endl;
      if (distance < RV_tol) {
      	node13= ii;
	//std::cout<< "node 1 coords" <<std::endl;
	//std::cout<< x1 <<","
	//	 << y1 <<","
	//	 << z1 <<std::endl;
      	break; 
      } else if((ii+1)==nnode1) { //is this else robust? 
	std::cerr<<"ERROR: couldn't find a matching RV node for purkinje node: "
		 << node3 << " \n";
	exit(1);
      }
    }
    node13_list.push_back(node13);
    //std::cout<< "node13-node3 couple is:" << node13<< "," << node3 
    //     << "(dist :" << distance << ")" <<std::endl;
  }

  std::cout<< "IEN_Mixed: Successfully matched every purkinje endnode." <<std::endl;
  //std::cout<< "node12-list:" <<std::endl;
  //VEC_T::print(node12_list);
  //std::cout<< "node13-list:" <<std::endl;
  //VEC_T::print(node13_list);

  //3- find locations of nodes in the IEN2  and IEN3
  std::vector<int>::iterator location;
  std::vector< std::vector<int> > node2_locations;
  node2_locations.resize(nnode2);//resize the 1st dimension only

  for (int ii=0; ii<nnode2 ; ii++) {
    location= std::find(vecIEN_2.begin(), vecIEN_2.end(), ii);
    
    while (location!=vecIEN_2.end()){
      //(node2_locations.at(ii)).push_back(*location);
      (node2_locations.at(ii)).push_back(std::distance(vecIEN_2.begin(),location));
      location= std::find(location+1, vecIEN_2.end(), ii);
    }
  }
   
  std::vector< std::vector<int> > node3_locations;
  node3_locations.resize(nnode3);//resize the 1st dimension only

  for (int ii=0; ii<nnode3 ; ii++) {
    location= std::find(vecIEN_3.begin(), vecIEN_3.end(), ii);
    
    while (location!=vecIEN_3.end()){
      //(node3_locations.at(ii)).push_back(*location);
      (node3_locations.at(ii)).push_back(std::distance(vecIEN_3.begin(),location));
      location= std::find(location+1, vecIEN_3.end(), ii);
    }
  }

  //std::cout<< "node2-locations -first 10" << std::endl;
  //for( auto it = node2_locations.begin(); it != node2_locations.begin()+10; it++ ){
  //  VEC_T::print(*it);
  //}
  
  
  //4- subtract nodes of endnodes from IEN2 and IEN3
  //   and keep the rest of the nodes separate
  std::vector<int> LVendnodes_reverse;
  LVendnodes_reverse.resize(LVendnodes.size());
  std::reverse_copy(LVendnodes.begin(), LVendnodes.end(), LVendnodes_reverse.begin());
  
  std::vector< std::vector<int> > locs2_replace, locs2_reorder;
  locs2_replace.resize(LVendnodes.size());
  //locs2_reorder.resize(nnode2-LVendnodes.size());
  locs2_reorder=node2_locations;
  
  int ii=0;
  for( auto it = LVendnodes.begin(); it != LVendnodes.end(); it++ ){
    //(locs2_replace.at(ii)).resize((node2_locations.at(*it)).size());
    locs2_replace.at(ii)= node2_locations.at(*it);
    ii++;
  }

  //  std::cout<< "LVendnodes reverse erases ctrlpts and displacements" << std::endl;
  for( auto it = LVendnodes_reverse.begin(); it != LVendnodes_reverse.end(); it++ ){
    //std::cout << "LVendnode: " << *it << "\n";
    locs2_reorder.erase(locs2_reorder.begin() + (*it));
    ctrlPts_2.erase(ctrlPts_2.begin()+(*it)*3, ctrlPts_2.begin()+(*it)*3+3);
    displacement_2.erase(displacement_2.begin()+(*it)*3,
			 displacement_2.begin()+(*it)*3+3); 
  }

  std::vector<int> RVendnodes_reverse;
  RVendnodes_reverse.resize(RVendnodes.size());
  std::reverse_copy(RVendnodes.begin(), RVendnodes.end(), RVendnodes_reverse.begin());
  
  std::vector< std::vector<int> > locs3_replace, locs3_reorder;
  locs3_replace.resize(RVendnodes.size());
  //locs3_reorder.resize(nnode3-RVendnodes.size());
  locs3_reorder=node3_locations;
  
  ii=0;
  for( auto it = RVendnodes.begin(); it != RVendnodes.end(); it++ ){
    //(locs3_replace.at(ii)).resize((node3_locations.at(*it)).size());
    locs3_replace.at(ii)= node3_locations.at(*it);
    ii++;
  }

  //  std::cout<< "RVendnodes reverse erases ctrlpts" << std::endl;
  for( auto it = RVendnodes_reverse.begin(); it != RVendnodes_reverse.end(); it++ ){
    //std::cout << "RVendnode: " << *it << "\n";
    locs3_reorder.erase(locs3_reorder.begin() + (*it));
    ctrlPts_3.erase(ctrlPts_3.begin()+(*it)*3, ctrlPts_3.begin()+(*it)*3+3);
    displacement_3.erase(displacement_3.begin()+(*it)*3,
			 displacement_3.begin()+(*it)*3+3); 
  }

  //4.5- update the control point coordinates according to the displacement
  //values
  if ( (ctrlPts_1.size() != displacement_1.size())
       || (ctrlPts_2.size() != displacement_2.size())
       || (ctrlPts_3.size() != displacement_3.size()) )  {
    std::cerr<<"ERROR: control points and displacement vector sizes don't match. \n";
    exit(1);
  }
  std::transform(ctrlPts_1.begin(), ctrlPts_1.end(), displacement_1.begin(),
		 ctrlPts_1.begin(),std::plus<double>()); 
  std::transform(ctrlPts_2.begin(), ctrlPts_2.end(),displacement_2.begin(),
		 ctrlPts_2.begin(),std::plus<double>());
  std::transform(ctrlPts_3.begin(), ctrlPts_3.end(),displacement_3.begin(),
		 ctrlPts_3.begin(), std::plus<double>()); 
		 
  
  //std::cout<< "ctrlpts2 before replace- first 10" << std::endl;
  //for( auto it = ctrlPts_2.begin(); it != ctrlPts_2.begin()+30; it=it+3 ){
  //  std::cout<< *(it) <<","
  //	     << *(it+1) <<","
  //	     << *(it+2) <<std::endl;
  //}
  //
  //std::cout<< "locs 2 replace" << std::endl;
  //for( auto it = locs2_replace.begin(); it != locs2_replace.end(); it++ ){
  //  VEC_T::print(*it);
  //}
  //
  //std::cout<< "locs 2 reorder - first 10" << std::endl;
  //for( auto it = locs2_reorder.begin(); it != locs2_reorder.begin()+10; it++ ){
  //  VEC_T::print(*it);
  //}
  //
  //std::cout<< "ctrlpts2 after replace- first 10" << std::endl;
  //for( auto it = ctrlPts_2.begin(); it != ctrlPts_2.begin()+30; it=it+3 ){
  //  std::cout<< *(it) <<","
  //	     << *(it+1) <<","
  //	     << *(it+2) <<std::endl;
  //}

  //5- merge the numbering of endnodes of IEN2 and IEN3 to IEN1
  //   and renumber rest of the nodes
  ii=0;
  for( auto it = locs2_replace.begin(); it != locs2_replace.end(); it++ ){
    for( auto it2 = (*it).begin(); it2 != (*it).end(); it2++ ){
      vecIEN_2.at(*it2) = node12_list.at(ii);
    }
    ii++;
  }

  ii=0;
  for( auto it = locs3_replace.begin(); it != locs3_replace.end(); it++ ){
    for( auto it2 = (*it).begin(); it2 != (*it).end(); it2++ ){
      vecIEN_3.at(*it2) = node13_list.at(ii);
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
  for( auto it = locs3_reorder.begin(); it != locs3_reorder.end(); it++ ){
    for( auto it2 = (*it).begin(); it2 != (*it).end(); it2++ ){
      vecIEN_3.at(*it2) = ii;
    }
    ii++;
  }
  nFunc_tot= ii;

  //6- store IENs as 2d vectors
  std::vector< std::vector <int> > IEN1, IEN2, IEN3;
  IEN1.resize(nelem1);  IEN2.resize(nelem2); IEN3.resize(nelem3); 
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
  for( auto it = IEN3.begin(); it != IEN3.end(); it++ ){
    (*it).resize(nLocBas3);
  }
  for (int jj=0; jj<vecIEN_3.size() ; jj++){
    int row= jj/nLocBas3;
    int col= jj%nLocBas3;
    (IEN3.at(row)).at(col)=vecIEN_3.at(jj);
  }

  IEN1.insert(IEN1.end(), IEN2.begin(), IEN2.end());
  IEN1.insert(IEN1.end(), IEN3.begin(), IEN3.end());

  //copy IEN1 to IEN
  IEN=IEN1;
  //std::cout<< "ien  " << std::endl;
  //for( auto it = IEN.begin(); it != IEN.end(); it++ ){
  //  VEC_T::print(*it);
  //}

  //assign output parameters
  nElem_tot = IEN.size();
  //nFunc = ii  assigned above

  ctrlPts_combined=ctrlPts_1;
  ctrlPts_combined.insert(ctrlPts_combined.end(),
			  ctrlPts_2.begin(),ctrlPts_2.end());
  ctrlPts_combined.insert(ctrlPts_combined.end(),
			  ctrlPts_3.begin(),ctrlPts_3.end());

  //elemType_combined=std::vector<int> (nelem1, elemType_list.at(0));
  //elemType_combined.insert(elemType_combined.end(),
  //                         nelem2, elemType_list.at(1));
  
  std::cout << "nElem_tot= " << nElem_tot << std::endl;
  std::cout << "nFunc_tot= " << nFunc_tot << std::endl;
  
}

IEN_Mixed::IEN_Mixed(const std::vector< std::vector<int> > &IEN_list,
		     const std::vector< std::vector<int> > &IEN_list_G,
		     const std::vector<IMesh *>  &mesh_list,
		     const std::vector<IMesh *>  &mesh_list_G,
		     const std::vector< std::vector<double> > &ctrlPts_list,
		     const std::vector< std::vector<double> > &ctrlPts_list_G,  
		     const std::string &LVendnodes_filename,
		     const std::string &RVendnodes_filename,
		     std::vector<double> &ctrlPts_combined,
		     const double &LV_tol,
		     const double &RV_tol)
{
  if(mesh_list.size() != IEN_list.size())
    {
      std::cerr<<"ERROR: IEN_list and mesh_list sizes don't match. \n";
      exit(1);
    }

  auto vecIEN_1= IEN_list.at(0);
  auto vecIEN_2= IEN_list.at(1);
  auto vecIEN_3= IEN_list.at(2);
  IMesh *Mesh_1  = mesh_list.at(0);
  IMesh *Mesh_2  = mesh_list.at(1);
  IMesh *Mesh_3  = mesh_list.at(2);
  int nnode1= Mesh_1->get_nFunc();
  int nnode2= Mesh_2->get_nFunc();
  int nnode3= Mesh_3->get_nFunc();
  int nelem1= Mesh_1->get_nElem();
  int nelem2= Mesh_2->get_nElem();
  int nelem3= Mesh_3->get_nElem();
  int nLocBas1= Mesh_1->get_nLocBas();
  int nLocBas2= Mesh_2->get_nLocBas();
  int nLocBas3= Mesh_3->get_nLocBas();  
  
  //1- read endnodes file contents
  std::vector<int> LVendnodes; 
  std::ifstream LVendnode_file(LVendnodes_filename);
  if (!LVendnode_file)
    {
      // Print an error and exit
      std::cerr << LVendnodes_filename
		<< "LV endnodes could not be opened for reading!" << std::endl;
      exit(1);
    }
  int number2 ;
  while (LVendnode_file>>number2) {
    LVendnodes.push_back(number2);
  }
  VEC_T::sort_unique_resize(LVendnodes);
  
  std::vector<int> RVendnodes; 
  std::ifstream RVendnode_file(RVendnodes_filename);
  if (!RVendnode_file)
    {
      // Print an error and exit
      std::cerr << RVendnodes_filename
		<< "RV endnodes could not be opened for reading!" << std::endl;
      exit(1);
    }
  int number3 ;
  while (RVendnode_file>>number3) {
    RVendnodes.push_back(number3);
  }
  VEC_T::sort_unique_resize(RVendnodes);

  //2-Now find the closest nodes on Mesh1 to the nodes on Mesh2
  std::vector<double> ctrlPts_1= ctrlPts_list.at(0);
  std::vector<double> ctrlPts_2= ctrlPts_list.at(1);
  std::vector<double> ctrlPts_3= ctrlPts_list.at(2);
  std::vector<double> ctrlPts_1_G= ctrlPts_list_G.at(0);
  std::vector<double> ctrlPts_2_G= ctrlPts_list_G.at(1);
  std::vector<double> ctrlPts_3_G= ctrlPts_list_G.at(2);
  std::vector<int> node12_list, node13_list;
  //node1_list.resize(endnodes.size());
  //std::cout<< "node1s list size:" << node1_list.size() <<std::endl;

  double distance, x1, x2, x3, y1, y2, y3, z1, z2, z3;
  
  //std::cout<< "nnode1:" << nnode1 <<std::endl;
  //std::cout<< "nnode2:" << nnode2 <<std::endl;
  //std::cout<< "nnode3:" << nnode3 <<std::endl;
  int node12, node2, node13, node3;

  //Match myocardium with LV purkinje 
  for (auto it2=LVendnodes.begin(); it2!=LVendnodes.end(); ++it2) {

    node2=*it2;
    x2=ctrlPts_2.at(3*node2);
    y2=ctrlPts_2.at(3*node2+1);
    z2=ctrlPts_2.at(3*node2+2);
    //std::cout<< "node 2 coords" <<std::endl;
    //std::cout<< x2 <<","
    //	     << y2 <<","
    //	     << z2 <<std::endl;

    for (int ii=0; ii<nnode1; ii++) {
      x1=ctrlPts_1.at(3*ii+0);
      y1=ctrlPts_1.at(3*ii+1);
      z1=ctrlPts_1.at(3*ii+2);

      //std::cout<< "node 1 coords" <<std::endl;
      //std::cout<< x1 <<","
      //	       << y1 <<","
      //	       << z1 <<std::endl;
      distance = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
      //std::cout<< "distance:" << distance <<std::endl;
      if (distance < LV_tol) { // turn this tolerance value into a user param
      	node12= ii;
	//std::cout<< "node 1 coords" <<std::endl;
	//std::cout<< x1 <<","
	//	 << y1 <<","
	//	 << z1 <<std::endl;
      	break; 
      } else if((ii+1)==nnode1) { //is this else robust? 
	std::cerr<<"ERROR: couldn't find an LV matching node for purkinje node: "
		 << node2 << " \n";
	exit(1);
      }
    }
    node12_list.push_back(node12);
    //std::cout<< "node12-node2 couple is:" << node12<< "," << node2 
    //     << "(dist :" << distance << ")" <<std::endl;
  }

  //Match myocardium with RV purkinje 
  for (auto it2=RVendnodes.begin(); it2!=RVendnodes.end(); ++it2) {
    node3=*it2;
    x3=ctrlPts_3.at(3*node3);
    y3=ctrlPts_3.at(3*node3+1);
    z3=ctrlPts_3.at(3*node3+2);
    //std::cout<< "node 3 coords" <<std::endl;
    //std::cout<< x3 <<","
    //	     << y3 <<","
    //	     << z3 <<std::endl;

    for (int ii=0; ii<nnode1; ii++) {
      x1=ctrlPts_1.at(3*ii+0);
      y1=ctrlPts_1.at(3*ii+1);
      z1=ctrlPts_1.at(3*ii+2);

      //std::cout<< "node 1 coords" <<std::endl;
      //std::cout<< x1 <<","
      //	       << y1 <<","
      //	       << z1 <<std::endl;
      distance = sqrt((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3)+(z1-z3)*(z1-z3));
      //std::cout<< "distance:" << distance <<std::endl;
      if (distance < RV_tol) {
      	node13= ii;
	//std::cout<< "node 1 coords" <<std::endl;
	//std::cout<< x1 <<","
	//	 << y1 <<","
	//	 << z1 <<std::endl;
      	break; 
      } else if((ii+1)==nnode1) { //is this else robust? 
	std::cerr<<"ERROR: couldn't find a matching RV node for purkinje node: "
		 << node3 << " \n";
	exit(1);
      }
    }
    node13_list.push_back(node13);
    //std::cout<< "node13-node3 couple is:" << node13<< "," << node3 
    //     << "(dist :" << distance << ")" <<std::endl;
  }

  std::cout<< "IEN_Mixed: Successfully matched every purkinje endnode." <<std::endl;
  //std::cout<< "node12-list:" <<std::endl;
  //VEC_T::print(node12_list);
  //std::cout<< "node13-list:" <<std::endl;
  //VEC_T::print(node13_list);

  //3- find locations of nodes in the IEN2  and IEN3
  std::vector<int>::iterator location;
  std::vector< std::vector<int> > node2_locations;
  node2_locations.resize(nnode2);//resize the 1st dimension only

  for (int ii=0; ii<nnode2 ; ii++) {
    location= std::find(vecIEN_2.begin(), vecIEN_2.end(), ii);
    
    while (location!=vecIEN_2.end()){
      //(node2_locations.at(ii)).push_back(*location);
      (node2_locations.at(ii)).push_back(std::distance(vecIEN_2.begin(),location));
      location= std::find(location+1, vecIEN_2.end(), ii);
    }
  }
   
  std::vector< std::vector<int> > node3_locations;
  node3_locations.resize(nnode3);//resize the 1st dimension only

  for (int ii=0; ii<nnode3 ; ii++) {
    location= std::find(vecIEN_3.begin(), vecIEN_3.end(), ii);
    
    while (location!=vecIEN_3.end()){
      //(node3_locations.at(ii)).push_back(*location);
      (node3_locations.at(ii)).push_back(std::distance(vecIEN_3.begin(),location));
      location= std::find(location+1, vecIEN_3.end(), ii);
    }
  }

  //std::cout<< "node2-locations -first 10" << std::endl;
  //for( auto it = node2_locations.begin(); it != node2_locations.begin()+10; it++ ){
  //  VEC_T::print(*it);
  //}
  
  
  //4- subtract nodes of endnodes from IEN2 and IEN3
  //   and keep the rest of the nodes separate
  std::vector<int> LVendnodes_reverse;
  LVendnodes_reverse.resize(LVendnodes.size());
  std::reverse_copy(LVendnodes.begin(), LVendnodes.end(), LVendnodes_reverse.begin());
  
  std::vector< std::vector<int> > locs2_replace, locs2_reorder;
  locs2_replace.resize(LVendnodes.size());
  //locs2_reorder.resize(nnode2-LVendnodes.size());
  locs2_reorder=node2_locations;
  
  int ii=0;
  for( auto it = LVendnodes.begin(); it != LVendnodes.end(); it++ ){
    //(locs2_replace.at(ii)).resize((node2_locations.at(*it)).size());
    locs2_replace.at(ii)= node2_locations.at(*it);
    ii++;
  }

  //  std::cout<< "LVendnodes reverse erases ctrlpts" << std::endl;
  for( auto it = LVendnodes_reverse.begin(); it != LVendnodes_reverse.end(); it++ ){
    //std::cout << "LVendnode: " << *it << "\n";
    locs2_reorder.erase(locs2_reorder.begin() + (*it));
    ctrlPts_2_G.erase(ctrlPts_2_G.begin()+(*it)*3, ctrlPts_2_G.begin()+(*it)*3+3);
  }

  std::vector<int> RVendnodes_reverse;
  RVendnodes_reverse.resize(RVendnodes.size());
  std::reverse_copy(RVendnodes.begin(), RVendnodes.end(), RVendnodes_reverse.begin());
  
  std::vector< std::vector<int> > locs3_replace, locs3_reorder;
  locs3_replace.resize(RVendnodes.size());
  //locs3_reorder.resize(nnode3-RVendnodes.size());
  locs3_reorder=node3_locations;
  
  ii=0;
  for( auto it = RVendnodes.begin(); it != RVendnodes.end(); it++ ){
    //(locs3_replace.at(ii)).resize((node3_locations.at(*it)).size());
    locs3_replace.at(ii)= node3_locations.at(*it);
    ii++;
  }

  //  std::cout<< "RVendnodes reverse erases ctrlpts" << std::endl;
  for( auto it = RVendnodes_reverse.begin(); it != RVendnodes_reverse.end(); it++ ){
    //std::cout << "RVendnode: " << *it << "\n";
    locs3_reorder.erase(locs3_reorder.begin() + (*it));
    ctrlPts_3_G.erase(ctrlPts_3_G.begin()+(*it)*3, ctrlPts_3_G.begin()+(*it)*3+3);
  }
  
  //std::cout<< "ctrlpts2 before replace- first 10" << std::endl;
  //for( auto it = ctrlPts_2.begin(); it != ctrlPts_2.begin()+30; it=it+3 ){
  //  std::cout<< *(it) <<","
  //	     << *(it+1) <<","
  //	     << *(it+2) <<std::endl;
  //}
  //
  //std::cout<< "locs 2 replace" << std::endl;
  //for( auto it = locs2_replace.begin(); it != locs2_replace.end(); it++ ){
  //  VEC_T::print(*it);
  //}
  //
  //std::cout<< "locs 2 reorder - first 10" << std::endl;
  //for( auto it = locs2_reorder.begin(); it != locs2_reorder.begin()+10; it++ ){
  //  VEC_T::print(*it);
  //}
  //
  //std::cout<< "ctrlpts2 after replace- first 10" << std::endl;
  //for( auto it = ctrlPts_2.begin(); it != ctrlPts_2.begin()+30; it=it+3 ){
  //  std::cout<< *(it) <<","
  //	     << *(it+1) <<","
  //	     << *(it+2) <<std::endl;
  //}

  //5- merge the numbering of endnodes of IEN2 and IEN3 to IEN1
  //   and renumber rest of the nodes
  ii=0;
  for( auto it = locs2_replace.begin(); it != locs2_replace.end(); it++ ){
    for( auto it2 = (*it).begin(); it2 != (*it).end(); it2++ ){
      vecIEN_2.at(*it2) = node12_list.at(ii);
    }
    ii++;
  }

  ii=0;
  for( auto it = locs3_replace.begin(); it != locs3_replace.end(); it++ ){
    for( auto it2 = (*it).begin(); it2 != (*it).end(); it2++ ){
      vecIEN_3.at(*it2) = node13_list.at(ii);
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
  for( auto it = locs3_reorder.begin(); it != locs3_reorder.end(); it++ ){
    for( auto it2 = (*it).begin(); it2 != (*it).end(); it2++ ){
      vecIEN_3.at(*it2) = ii;
    }
    ii++;
  }
  nFunc_tot= ii;

  //6- store IENs as 2d vectors
  std::vector< std::vector <int> > IEN1, IEN2, IEN3;
  IEN1.resize(nelem1);  IEN2.resize(nelem2); IEN3.resize(nelem3); 
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
  for( auto it = IEN3.begin(); it != IEN3.end(); it++ ){
    (*it).resize(nLocBas3);
  }
  for (int jj=0; jj<vecIEN_3.size() ; jj++){
    int row= jj/nLocBas3;
    int col= jj%nLocBas3;
    (IEN3.at(row)).at(col)=vecIEN_3.at(jj);
  }

  IEN1.insert(IEN1.end(), IEN2.begin(), IEN2.end());
  IEN1.insert(IEN1.end(), IEN3.begin(), IEN3.end());


  //copy IEN1 to IEN
  IEN=IEN1;
  //std::cout<< "ien  " << std::endl;
  //for( auto it = IEN.begin(); it != IEN.end(); it++ ){
  //  VEC_T::print(*it);
  //}

  //assign output parameters
  nElem_tot = IEN.size();
  //nFunc = ii  assigned above

  ctrlPts_combined=ctrlPts_1_G;
  ctrlPts_combined.insert(ctrlPts_combined.end(),
			  ctrlPts_2_G.begin(),ctrlPts_2_G.end());
  ctrlPts_combined.insert(ctrlPts_combined.end(),
			  ctrlPts_3_G.begin(),ctrlPts_3_G.end());

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
