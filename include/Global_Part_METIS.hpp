#ifndef GLOBAL_PART_METIS_HPP
#define GLOBAL_PART_METIS_HPP
// ==================================================================
// Global_Part_METIS.hpp
// Object:
// Use METIS routien to get the global mesh partition stored as:
// element partition information: epart
// node partition information: npart.
//
// Date Created: Oct 2nd 2013
// ==================================================================
#include "IMesh.hpp"
#include "IIEN.hpp"
#include "Vec_Tools.hpp"
#include "IGlobal_Part.hpp"

class Global_Part_METIS : public IGlobal_Part
{
  public:
    // Constructor:
    // It will create eptr and eind arrays and call METIS_PartMeshDual or
    // METIS_PartMeshNodal for mesh partition
    Global_Part_METIS( const int &cpu_size,
        const int &in_ncommon, const bool &isDualGraph,
        const IMesh * const &mesh,
        const IIEN * const &IEN,
        const char * const &element_part_name,
        const char * const &node_part_name );

  Global_Part_METIS( const int &cpu_size,
		     const int &in_ncommon, const bool &isDualGraph,
		     const std::vector<IMesh *> &mesh_list,
		     const std::vector<IIEN *> &IEN_list,
		     const char * const &element_part_name,
		     const char * const &node_part_name );

    virtual ~Global_Part_METIS();

    virtual idx_t get_epart( const int &ee ) const {return epart[ee];}
    
    virtual idx_t get_npart( const int &nn ) const {return npart[nn];}

    virtual bool get_isMETIS() const {return isMETIS;};
    
    virtual bool get_isDual() const {return isDual;};
    
    virtual int get_dual_edge_ncommon() const {return dual_edge_ncommon;}
  
  private:
    const bool isMETIS, isDual;
    const int dual_edge_ncommon;

    idx_t * epart;
    idx_t * npart;

    virtual void write_part_hdf5( const char * const &fileName, 
        const idx_t * const &part_in,
        const int &part_size, const int &cpu_size,
        const bool &part_isdual, const int &in_ncommon,
        const bool &isMETIS ) const;

    // --------------------------------------------------------------
    // This function will write the data of part_in in 64bit HDF5 format. 
    // This function should be called if idx_t is the int64_t.
    // --------------------------------------------------------------
    virtual void write_part_hdf5_64bit( const char * const &fileName, 
        const int64_t * const &part_in,
        const int64_t &part_size, const int &cpu_size,
        const bool &part_isdual, const int &in_ncommon,
        const bool &isMETIS ) const;
};

#endif
