#ifndef IEN_MIXED_HPP
#define IEN_MIXED_HPP
// ==================================================================
// IEN_Mixed.hpp
//
// This objects gives the IEN array for a mesh that consists of
// elements of non-uniform type
//
// initial intention is to create a single IEN structure for
// electrophysiology, which involves 1d (line) and 3d(tet) elements
//
// Author: Oguz Ziya Tikenogullari 
// from the original of: IEN_Tetra_P1 by Ju Liu
// Year: 2021
// ==================================================================
#include <vector>
#include <fstream>
#include "IIEN.hpp"
#include "IMesh.hpp"
#include "Vec_Tools.hpp"


class IEN_Mixed : public IIEN
{
  public:
    // Constructor: take IENs created by other readers 
//  IEN_Mixed(const IIEN * const &IEN1,
//	    const IIEN * const &IEN2 );
  IEN_Mixed(const std::vector< std::vector<int> > &IEN_list,
	    const std::vector<IMesh *>  &mesh_list,
	    const std::vector< std::vector<double> > &ctrlPts_list,	    
	    const std::string &LVendnodes_filename,
	    const std::string &RVendnodes_filename,
	    std::vector<double> &ctrlPts_combined,
	    const double &LV_tol,
	    const double &RV_tol);

  IEN_Mixed(const std::vector< std::vector<int> > &IEN_list,
	    const std::vector< std::vector<int> > &IEN_list_G,
	    const std::vector<IMesh *>  &mesh_list,
	    const std::vector<IMesh *>  &mesh_list_G,
	    const std::vector< std::vector<double> > &ctrlPts_list,
	    const std::vector< std::vector<double> > &ctrlPts_list_G,  
	    const std::string &LVendnodes_filename,
	    const std::string &RVendnodes_filename,
	    std::vector<double> &ctrlPts_combined,
	    const double &LV_tol,
	    const double &RV_tol);
    
  ~IEN_Mixed();

  virtual int get_IEN( const int &ee, const int &l_node ) const;

  virtual void print_IEN() const;

  //int get_size_IEN () const {return size_IEN;}
  virtual int get_nFunc_tot () const {return nFunc_tot;}
  virtual int get_nElem_tot () const {return nElem_tot;}

private:
  std::vector<std::vector<int>> IEN;

  int nFunc_tot, nElem_tot;
  //int size_IEN;
};

#endif
