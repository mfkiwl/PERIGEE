#ifndef NODALBC_3D_VTP_HPP
#define NODALBC_3D_VTP_HPP
// ==================================================================
// NodalBC_3D_vtp.hpp
//
// This is an instantiation of INodalbc for 3D problems by reading the
// dirichlet nodes from the vtp file.
//
// This class is designed to handle the mesh for unstructural tetrahedral
// mesh generated by automatic mesher like tetgen.
// 
// The data contained in this class include:
// dir_nodes : the nodal indices for the Dirichlet nodes;
// num_dir_nodes : the number of the Dirichlet nodes, i.e., the length
//                 of the dir_nodes array;
// per_slave_nodes : the nodal indices for the slave nodes;
// per_master_nodes : the nodal indices for the master nodes;
// num_per_nodes : the number of periodic-type nodes, i.e., the length 
//                 of the per_slave_nodes / per_master_nodes.
//
// ID : the vector for the ID array, which is generated based on the 
//      nFunc, the number of total basis functions, and the dir_nodes
//      and per_slave_nodes/per_master_nodes. Once the dir_nodes, per_
//      xxx_nodes, and nFunc are given, the ID array will be generated
//      by the function create_ID.
//
// Date: Jan. 6 2017
// Author: Ju Liu
// ==================================================================
#include "INodalBC.hpp"
#include "Tet_Tools.hpp"
#include "Vector_3.hpp"

class NodalBC_3D_vtp : public INodalBC
{
  public:
    // --------------------------------------------------------------
    // Default constructor: clear the dir_nodes, per_slave_nodes,
    // per_master_nodes; set num_dir_nodes, num_per_nodes to be zero;
    // set ID based on the above "no-nodal bc" setting.
    // --------------------------------------------------------------
    NodalBC_3D_vtp( const int &nFunc );
   
    // --------------------------------------------------------------
    // Specify the Dirichlet nodes for CMM. This includes all inlet
    // nodes and the outline (`ring') nodes for each outlet surface.
    // Used for inlet & outlet clamping.
    // --------------------------------------------------------------
    NodalBC_3D_vtp( const std::string &inflow_vtp_file,
        const std::string &wall_vtp_file,
        const std::vector<std::string> &outflow_vtp_files,
        const int &nFunc );
     
    // --------------------------------------------------------------
    // Specify the Dirichlet nodes for CMM. This includes all interior
    // inlet nodes. For each inlet/outlet, ring nodes are also included
    // for the velocity dof corresponding to the unit normal's dominant
    // component. Used for inlet & outlet in-plane or purely radial motion.
    //     \para type: 0 for in-plane motion. 1 for purely radial motion.
    //     \para comp: velocity component. 0, 1, or 2.
    // --------------------------------------------------------------
    NodalBC_3D_vtp( const std::string &inflow_vtp_file,
        const std::vector<double> &inflow_outward_vec,
        const std::string &wall_vtp_file,
        const std::vector<std::string> &outflow_vtp_files,
        const std::vector< std::vector<double> > &outflow_outward_vec,
        const int &type, const int &comp, const int &nFunc );
     
    // --------------------------------------------------------------
    // The vtp file specifies the Dirichlet nodes. No periodical BC.
    // --------------------------------------------------------------
    NodalBC_3D_vtp( const std::string &vtpfileName, const int &nFunc );

    // --------------------------------------------------------------
    // The list of vtp files specifies the Dirichlet nodes. 
    // No periodical type BC nodes.
    // --------------------------------------------------------------
    NodalBC_3D_vtp( const std::vector<std::string> &vtpfileList, 
        const int &nFunc );

    virtual ~NodalBC_3D_vtp();

  private:
    NodalBC_3D_vtp() {};

    // Compute centroid coordinates given a cap's nodal coordinates
    void compute_cap_centroid( const std::vector<double> &pts, Vector_3 &centroid ) const;

    // Return the dominant component index of a ring node's unit tangential vector
    // \para outvec  : corresponding cap's unit outward normal
    // \para centroid: corresponding cap's centroidal coordinates
    // \para pt_x, pt_y, pt_z: nodal coordinates
    int compute_tangential( const Vector_3 &outvec, const Vector_3 &centroid,
        const int &pt_x, const int &pt_y, const int &pt_z );
};

#endif
