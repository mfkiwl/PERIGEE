#ifndef PGASSEM_EP_HPP
#define PGASSEM_EP_HPP
// ==================================================================
// PGAssem_EP.hpp
// Parallel Global Assembly routine.
//
// This class assembles the global tangent matrix / mass matrix, and
// residual vector, with essential boundary conditions (including 
// dirichlet bc, strongly enforced periodic bc).
// 
// Date:
// Author: Ju Liu, liujuy@gmail.com
// Modified: Oguz Ziya Tikenogullari, o.z.tikenogullari@gmail.com
// ==================================================================
#include "Sys_Tools.hpp"
#include "Vec_Tools.hpp"
#include "APart_Node.hpp"
#include "ALocal_Elem.hpp"
#include "ALocal_IEN.hpp"
#include "IAGlobal_Mesh_Info.hpp"
#include "IALocal_BC.hpp"
#include "IPLocAssem.hpp"
#include "FEANode.hpp"
#include "PDNSolution.hpp"
#include "IALocal_meshSize.hpp"
#include "ALocal_meshSize_3D_NURBS.hpp"
#include "BernsteinBasis_Array.hpp"
#include "IAExtractor.hpp"
#include "AExtractor_3D_NURBS_xyz.hpp"
#include "FEAElement.hpp"
#include "AInt_Weight.hpp"

#include "IonicModel.hpp"
#include "PDNTimeStep.hpp"

class PGAssem_EP
{
  public:
    // K is the tangent matrix
    Mat K;
    
    // G is the residual vector
    Vec G;

    PGAssem_EP( const IPLocAssem * const &locassem_ptr,
        const IAGlobal_Mesh_Info * const &agmi_ptr,
        const APart_Node * const &pnode_ptr,
        const int &petsc_version_type );
      
    // ------------------------------------------------------------------------
    // ! Constructor that does not need to specify the petsc version
    // ! for matrix creation. 
    // ------------------------------------------------------------------------
    PGAssem_EP( const IPLocAssem * const &locassem_ptr,
        const IAGlobal_Mesh_Info * const &agmi_ptr,
        const APart_Node * const &pnode_ptr );


    // ------------------------------------------------------------------------
    // ! Constructor that calls Get_dnz_onz to get a more accurate estimate
    //   nonzero structure of the sparse matrix.
    // ------------------------------------------------------------------------
    PGAssem_EP( const IPLocAssem * const &locassem_ptr,
        const IAGlobal_Mesh_Info * const &agmi_ptr,
        const ALocal_Elem * const &alelem_ptr,
        const ALocal_IEN * const &aien_ptr,
        const APart_Node * const &pnode_ptr,
        const IALocal_BC * const &part_bc );


    // ------------------------------------------------------------------------
    // ! Destructor
    // ------------------------------------------------------------------------
    virtual ~PGAssem_EP();

    // ------------------------------------------------------------------------
    // ! Flag : Fix nonzero structure
    // !        Add or insert in new location is ignored
    // ------------------------------------------------------------------------
    void Fix_nonzero_str() 
    {MatSetOption(K, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);}
    
   
    // ------------------------------------------------------------------------
    // ! Flag : New allocation error 
    //          Add or insert in new location will generate error
    // ------------------------------------------------------------------------
    void Fix_nonzero_err_str()
    {MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);}

    
    // ------------------------------------------------------------------------
    // ! Flag : Ignore new allocation 
    //          Add or insert in new location will NOT generate error
    // ------------------------------------------------------------------------
    void Release_nonzero_err_str()
    {MatSetOption(K, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);}
   
   
    // ------------------------------------------------------------------------
    // ! Flag : Keep nonzero pattern for the matrix K 
    // ------------------------------------------------------------------------
    void Keep_nonzero_pattern()
    {MatSetOption(K, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);} 
   
    // ------------------------------------------------------------------------
    // ! clear K G  to be zero
    // ------------------------------------------------------------------------
    void Clear_KG()
    {
      MatZeroEntries(K);
      VecSet(G, 0.0);
    }
   
    // ------------------------------------------------------------------------
    // ! clear G to be zero
    // ------------------------------------------------------------------------
    void Clear_G()
    {VecSet(G, 0.0);}

    // ------------------------------------------------------------------------
    // ! Assembly estimated nonzero structure
    // ------------------------------------------------------------------------
    void Assem_nonzero_estimate(
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr, 
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const IALocal_BC * const &bc_part );



    // ------------------------------------------------------------------------
    // ! Assembly tangent matrix and residual vector,
    //   if history variables exist.
    //       Needs the basis function quadrature info
    //       sol_a : disp /or sol_n+1
    //       sol_b : velo /or sol_n
    //       sol_c : hist
    //       sol_d : hist_dot
    // ------------------------------------------------------------------------
    void Assem_tangent_residual(
        const PDNSolution * const &sol_a,
        const PDNSolution * const &sol_b,
	//const PDNSolution * const &sol_c,//pre_hist
	//PDNSolution * const &sol_d,//new hist
        const double &t_n,
        const double &dt,
        //const double &dt_ion,
	//const IonicModel * const &ionicmodel_ptr,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr, 
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const AInt_Weight * const &wei_ptr,
        const std::vector<FEAElement*> &eptr_array,
        const IALocal_BC * const &bc_part );

    
    // ------------------------------------------------------------------------
    // ! Assembly tangent matrix and residual vector together.
    //     The quadrature info is computed within the assembly routine.
    //     This function is for 3D problems
    //     sol_a : disp or sol_{n+1}
    //     sol_b : velo or sol_{n}
    //   -----------------
    //   The following input are needed for element quadrature evaluation
    //     IALocal_meshSize : object saving the mesh size
    //     BernsteinBasis : pre-evaluated Bernstein polynomial at Gauss 
    //                      quadrature points for volumetric integration
    //     IAExtractor : Bezier extraction operator
    // ------------------------------------------------------------------------

    // ------------------------------------------------------------------------
    // ! Assembly residual vector when history variables are present
    //  sol_c is hist_dot, sol_d is hist
    // ------------------------------------------------------------------------
    void Assem_residual(
    const PDNSolution * const &sol_a, //velo
    const PDNSolution * const &sol_b, //disp
    //const PDNSolution * const &sol_c, //pre_hist
    //PDNSolution * const &sol_d, //new hist
    const double &t_n,
    const double &dt,
    //const double &dt_ion,
    //const IonicModel * const &ionicmodel_ptr,
    const ALocal_Elem * const &alelem_ptr,
    IPLocAssem * const &lassem_ptr, 
    const ALocal_IEN * const &lien_ptr,
    const APart_Node * const &node_ptr,
    const FEANode * const &fnode_ptr,
    const AInt_Weight * const &wei_ptr,
    const std::vector<FEAElement*> &eptr_array,
    const IALocal_BC * const &bc_part );

    

    // ------------------------------------------------------------------------
    // ! Assembly mass matrix and its residual vector
    // ------------------------------------------------------------------------
    void Assem_mass_residual(
        const PDNSolution * const &sol_a,
	//const PDNSolution * const &sol_b, //hist
	const PDNTimeStep * const &time_info,
	//const IonicModel * const &ionicmodel_ptr,
        const ALocal_Elem * const &alelem_ptr,
        IPLocAssem * const &lassem_ptr, 
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const FEANode * const &fnode_ptr,
        const AInt_Weight * const &wei_ptr,
        const std::vector<FEAElement*> &eptr_array,
        const IALocal_BC * const &bc_part );


    // ------------------------------------------------------------------------
    // ! Assembly mass matrix and its residual vector for 3D problems
    //   without cached quadrature info.
    // ! -----------------
    // ! The following input are needed for element quadrature evaluation
    // ! IALocal_meshSize : object saving the mesh size
    // ! BernsteinBasis : pre-evaluated Bernstein polynomial at Gauss 
    // !                  quadrature points for volumetric integration
    // ! IAExtractor : Bezier extraction operator
    // ------------------------------------------------------------------------
  void Update_nodal_velo(const PDNSolution * const &sol_a, //disp
			 //const PDNSolution * const &sol_d, //velo
			 const PDNSolution * const &sol_b, //pre_hist
			 const double &t_n,
			 const double &dt,
			 const IonicModel * const &ionicmodel_ptr,
			 const ALocal_Elem * const &alelem_ptr,
			 //IPLocAssem * const &lassem_ptr, 
			 const ALocal_IEN * const &lien_ptr,
			 const APart_Node * const &node_ptr,
			 const FEANode * const &fnode_ptr,
			 //const AInt_Weight * const &wei_ptr,
			 const std::vector<FEAElement*> &eptr_array,
			 const IALocal_BC * const &bc_part,
			 PDNSolution * const &sol_c, //new hist
			 PDNSolution * const &sol_d //new hist
			 );  

    // ------------------------------------------------------------------------
    // ! Print_G : print the residual vector G on screen
    // ------------------------------------------------------------------------
    void Print_G() const
    {VecView(G, PETSC_VIEWER_STDOUT_WORLD);}


  private:
    // ------------------------------------------------------------------------
    // ! Initialization for MPIAIJ matrix using PETSc-3.2-px
    // ------------------------------------------------------------------------
    void Init_petsc_32(const int &nonzero_per_row, const int &num_loc_row);


    // ------------------------------------------------------------------------
    // ! Initialization for MPIAIJ matrix using PETSc-3.5
    // ------------------------------------------------------------------------
    void Init_petsc_35(const int &nonzero_per_row, const int &num_loc_row);


    // ------------------------------------------------------------------------
    // ! Initialization for MPIAIJ matrix given the estimate of nonzero
    //   per row in diagonal and off-diagonal portion.
    // ------------------------------------------------------------------------
    void Init_petsc_35(const int &dnz, const int &onz, const int &num_loc_row);


    // ------------------------------------------------------------------------
    // ! Initialization for MPIAIJ matrix using PETSc-3.5
    //   The nonzero of diagonal portion and off-diagonal portion are given
    //   by the d_nnz and o_nnz PetscInt arrays with length nlocalnode*dof
    //   = nlocrow. These arrays can be obtained by calling Get_dnz_onz function
    // ------------------------------------------------------------------------
    void Init_petsc_35(const PetscInt * const &dnz, 
        const PetscInt * const &onz, const int &num_loc_row );


    // ------------------------------------------------------------------------
    // ! Enforce essential bc for tangent matrix and residual vector
    // Note: This function use ADD_VALUES for matrix K. Hence, this function
    // should be used with LID in assemly, which avoids putting entries in 
    // essential bc nodes. The advantage is that only one MatAssemblyBegin/End
    // pair needs to be called.
    // ------------------------------------------------------------------------
    void EssBC_KG( const IALocal_BC * const &bc_part, const int &field );

    // ------------------------------------------------------------------------
    // ! Enforce essential bc for the residual vector
    // ------------------------------------------------------------------------
    void EssBC_G( const IALocal_BC * const &bc_part, const int &field );


    // ------------------------------------------------------------------------
    // ! Get local element's array.
    //   the long array should be listed compatible with the LIEN
    //   structure, for example, local_to_global() ~ LIEN[][].
    //   user is responsible for creating and deleting local_array pointer
    // ------------------------------------------------------------------------
    void GetLocal(const double * const &array, const int * const &IEN,
        double * const &local_array) const
    {
      int offset1, offset2;
      for(int ii=0; ii<nLocBas; ++ii)
      {
        offset1 = ii * dof;
        offset2 = IEN[ii] * dof;
        for(int jj=0; jj<dof; ++jj)
          local_array[offset1 + jj] = array[offset2 + jj];
      }
    }


    // ------------------------------------------------------------------------
    // ! Get_dnz_onz
    //   This function gives an estimate for the number of the diagonal
    //   nonzeros and the off-diagonal nonzeros. The results are returned
    //   as two arrays with length dof * nlocalnode.
    //   This function is used to give an accurate matrix preallocation
    //   for the PETSc matrix Mat.
    //   Output : 
    //   \para dnz : the number of nonzeros in diagonal portion
    //   \para onz : the number of nonzeros in off-diagonal portion
    // ------------------------------------------------------------------------
    void Get_dnz_onz( const int &nElem,
        const ALocal_IEN * const &lien_ptr,
        const APart_Node * const &node_ptr,
        const IALocal_BC * const &bc_part,
        PetscInt * const &dnz, PetscInt * const &onz ) const;


    // ------------------------------------------------------------------------
    // ! Get_part_id
    //   Get the partition id for given node index from the given list of nodal
    //   partition list. nlist has to start as nlist[0] = 0;
    //   nlist = {0, 3, 5, 8} means that partition assigns 3 cpus.
    //   cpu 0 owns node 0 to node 2.
    //   cpu 1 owns node 3 to node 4.
    //   cpu 2 owns node 5 to node 8.
    //   if val = 1, returns 0.
    //   if val > *(nlist.end()), return -1.
    // ------------------------------------------------------------------------
    int Get_part_id( const std::vector<unsigned int> &nlist, const unsigned int &val ) const
    {
      for(unsigned int ii=1; ii<nlist.size(); ++ii)
      {
        if(val < nlist[ii])
          return int(ii-1);
      }
      return -1;
    }

    // ------------------------------------------------------------------------
    // ! Private data structures
    // ------------------------------------------------------------------------
    //   nLocBas: number of basis per element
    int nLocBas;

    //   dof: degree of freedom of this problem
    int dof;

    //   row_index and col_index 
    //     They are used in the  assembly routines to store element's LID 
    //     and local_to_global node indeices.
    PetscInt * row_index;
    PetscInt * col_index;

    //   array_a & array_b are used to store this partition's nodes'
    //   solution vector
    double * array_a;
    double * array_b;
    double * array_c;
    double * array_d;

    //   local_a & local_b are element solution vector, i.e. with size
    //   nLocBas * dof
    double * local_a;
    double * local_b;
    double * local_c;
    double * local_d;  

    //   IEN_e: LIEN indices in one particular element
    int * IEN_e;

    //   ectrl_x/y/z: element's control points' xyz coordinates
    double * ectrl_x;
    double * ectrl_y;
    double * ectrl_z;
};

#endif
