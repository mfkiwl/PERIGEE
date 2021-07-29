#ifndef AGLOBAL_MESH_INFO_MIXED_HPP
#define AGLOBAL_MESH_INFO_MIXED_HPP
// ==================================================================
// AGlobal_Mesh_Info_Mixed.hpp
//
// This is the instantiation of global mesh info class which is used
// for mesh with element type non uniform.
//
// Date Created: Jan 21 2017
// Modified: 2021 Oguz Ziya Tikenogullari
// ==================================================================
#include "IAGlobal_Mesh_Info.hpp"
#include "HDF5_Reader.hpp"
#include "Vec_Tools.hpp"

class AGlobal_Mesh_Info_Mixed : public IAGlobal_Mesh_Info
{
  public:
    AGlobal_Mesh_Info_Mixed( const std::string &fileBaseName,
        const int &cpu_rank );

    // Construct a global mesh info based on a mesh by enriching it
    // in each cell with additional bubble nodes
    // Input:  num_enrich_node : the number of enriched nodes
    // Output: nFunc = original nFunc + nElem * num_enrich_node
    //         nLocBas = nLocBas + num_enrich_node
    //         elemType = elemType + 10
    AGlobal_Mesh_Info_Mixed( const std::string &fileBaseName,
        const int &cpu_rank, const int &num_enrich_node );
    
    virtual ~AGlobal_Mesh_Info_Mixed();

    virtual int get_xdegree() const {SYS_T::print_fatal("Error: AGLobal_Mesh_Info_Mixed::get_xdegree is not implemented. \n"); return -1;}
    virtual int get_ydegree() const {SYS_T::print_fatal("Error: AGLobal_Mesh_Info_Mixed::get_ydegree is not implemented. \n"); return -1;}
    virtual int get_zdegree() const {SYS_T::print_fatal("Error: AGLobal_Mesh_Info_Mixed::get_zdegree is not implemented. \n"); return -1;}
    virtual int get_xdegree(const int &ee) const {return xdegree.at(ee);}
    virtual int get_ydegree(const int &ee) const {return xdegree.at(ee);}
    virtual int get_zdegree(const int &ee) const {return xdegree.at(ee);}

    virtual double get_max_hx() const
    {SYS_T::print_fatal("Error: AGlobal_Mesh_Info_Mixed::get_max_hx is not implemented. \n"); return 0.0;}

    virtual double get_max_hy() const
    {SYS_T::print_fatal("Error: AGlobal_Mesh_Info_Mixed::get_max_hy is not implemented. \n"); return 0.0;}

    virtual double get_max_hz() const
    {SYS_T::print_fatal("Error: AGlobal_Mesh_Info_Mixed::get_max_hz is not implemented. \n"); return 0.0;}
    
    virtual int get_nElem() const {return nElem;}

    virtual int get_nFunc() const {return nFunc;}

    virtual int get_nLocBas() const {SYS_T::print_fatal("Error: AGLobal_Mesh_Info_Mixed::get_nLocBas is not implemnted. \n"); return -1;}
  
    virtual int get_nLocBas(const int &ee) const {return nLocBas.at(ee);}

    virtual int get_probDim() const {return probDim;}

    virtual int get_elemType() const {SYS_T::print_fatal("Error: AGLobal_Mesh_Info_Mixed::get_elemType is not implemented. \n"); return -1;}

    virtual int get_elemType(const int &ee) const {return elemType.at(ee);}

    virtual void print_info() const;

  private:
  std::vector< int >  xdegree, ydegree, zdegree;
  std::vector< int >  elemType, nLocBas;
  int nElem, nFunc, probDim ;
};

#endif
