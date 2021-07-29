#ifndef PART_MIXED_MESH_HPP
#define PART_MIXED_MESH_HPP
// ==================================================================
// Part_Mixed_Mesh.hpp
//
// Object: Partition mixed mesh into subdomains,
//         and record each subdomain in hdf5 file format.
//
// Date: Jan. 12 2017
// Modified on May 22 2017: I added a new constructor to handle the 
//         cases for dofNum is not the dof for the matrix problem in
//         the implicit solver.
// Modified 2021: derived partition function for nonuniform elements
// ==================================================================
#include "IMesh.hpp"
#include "IPart.hpp"

class Part_Mixed_Mesh : public IPart
{
  public:

  Part_Mixed_Mesh( const IMesh * const &mesh,
		   const IGlobal_Part * const &gpart,
		   const Map_Node_Index * const &mnindex,
		   const IIEN * const &IEN,
		   const std::vector<double> &ctrlPts,	
		   //const std::vector<std::vector<double>> &myo_fiber,	     
		   const int &in_cpu_rank, const int &in_cpu_size,
		   const int &in_dofNum, const int &in_dofMat,
		   const bool isPrintInfo );

    virtual ~Part_Mixed_Mesh();

    virtual void write( const char * inputFileName ) const;
    
    virtual bool isElemInPart(int gloindex) const
    {return VEC_T::is_invec(elem_loc, gloindex);}
    
    virtual bool isNodeInPart(int gloindex) const
    {return VEC_T::is_invec(node_loc, gloindex);}
   
    // Determine the position of a given index in the elem_loc array 
    virtual int get_elemLocIndex(const int &gloindex) const;

    // Determine the position of a given index in the local_to_global
    // array
    virtual int get_nodeLocGhoIndex(const int &gloindex) const
    {return VEC_T::get_pos(local_to_global, gloindex);}

    virtual void print_part_ele() const;

    virtual void print_part_node() const;

    virtual void print_part_ghost_node() const;

    virtual void print_part_local_to_global() const;

    virtual void print_part_LIEN(const IMesh * const &mesh) const;

    virtual void print_part_loadbalance_edgecut() const;

    virtual int get_elem_loc(int pos) const {return elem_loc[pos];}
    virtual int get_nlocalele() const {return nlocalele;}
    virtual int get_node_loc(int pos) const {return node_loc[pos];}
    virtual int get_node_loc_original(int pos) const {return node_loc_original[pos];}
    virtual int get_node_ghost(int pos) const {return node_ghost[pos];}
    virtual int get_local_to_global(int pos) const {return local_to_global[pos];}
    virtual int get_nlocalnode() const {return nlocalnode;}
    virtual int get_nghostnode() const {return nghostnode;}
    virtual int get_ntotalnode() const {return ntotalnode;}
    virtual int get_nbadnode() const {return nbadnode;}
    virtual int get_nlocghonode() const {return nlocghonode;}
    virtual int get_cpu_rank() const {return cpu_rank;}
    virtual int get_cpu_size() const {return cpu_size;}
    virtual bool get_part_isDual() const {return part_isdual;}
    virtual bool get_isMETIS() const {return isMETIS;}
    virtual int get_dual_edge_ncommon() const {return dual_edge_ncommon;}

    virtual int get_nElem() const {return nElem;}
    virtual int get_nFunc() const {return nFunc;}
    virtual int get_sDegree() const
      {std::cerr<<"ERROR: part.get_degree "<<std::endl;return 0;}
    virtual int get_tDegree() const
      {std::cerr<<"ERROR: part.get_degree "<<std::endl;return 0;}
    virtual int get_uDegree() const
      {std::cerr<<"ERROR: part.get_degree "<<std::endl;return 0;}
    virtual int get_nLocBas() const
      {std::cerr<<"ERROR: part.get_nlocbas "<<std::endl;return 0;}
    virtual int get_LIEN(int ee, int ii) const {return LIEN[ee][ii];}

    virtual double get_ctrlPts_x_loc(int pos) const {return ctrlPts_x_loc[pos];}
    virtual double get_ctrlPts_y_loc(int pos) const {return ctrlPts_y_loc[pos];}
    virtual double get_ctrlPts_z_loc(int pos) const {return ctrlPts_z_loc[pos];}

  protected:
  // 1. local element
  std::vector<int> elem_loc;
  int nlocalele;
  std::vector< std::vector< double > > fiber_ori_loc; 

  // 2. local node
  std::vector<int> node_loc;
  std::vector<int> node_loc_original;
  std::vector<int> node_ghost;
  std::vector<int> local_to_global;

  int nlocalnode;
  int nghostnode;
  int ntotalnode;
  int nbadnode;
  int nlocghonode;

  // 3. CPU info and partition parameters
  int cpu_rank, cpu_size;
  bool isMETIS, part_isdual;
  int dual_edge_ncommon;

  // 4. global mesh info
  const int nElem, nFunc;
  std::vector< std::vector< int > >  stu_Degree; //stu_Degree_loc;
  std::vector<int>  nLocBas, elemType, nLocBas_loc; // elemType_loc;
  const int probDim, dofNum, dofMat;

  // 5. LIEN
  int ** LIEN;
  //    std::vector< std::vector< int > > LIEN;

  // 6. local point coordinates (i.e. control point geometry)
  std::vector<double> ctrlPts_x_loc, ctrlPts_y_loc, ctrlPts_z_loc;

  void Generate_Partition( const IMesh * const &mesh,
			   const IGlobal_Part * const &gpart,
			   const Map_Node_Index * const &mnindex,
			   const IIEN * const &IEN,
			   const std::vector<double> &ctrlPts,
			   //const std::vector<std::vector<double>> &myo_fiber,
			   const bool &isPrintinfo );
};

#endif
