#include "ElemBC_3D_turbulence_wall_model.hpp"

ElemBC_3D_turbulence_wall_model::ElemBC_3D_turbulence_wall_model( 
    const std::vector<std::string> &vtkfileList,
    const int &in_wall_model_type, const IIEN * const &VIEN, 
    const int &elemtype )
: ElemBC_3D ( vtkfileList, elemtype ), 
  wall_model_type {in_wall_model_type}
{
  SYS_T::print_fatal_if(VEC_T::get_size(vtkfileList) > 1, "Error, ElemBC_3D_turbulence_wall_model: The number of wall file should not be more than 1.\n");

  if(VEC_T::get_size(vtkfileList) == 1)
  {
    face_id.resize(num_cell[0]);

    if(elem_type == 501 || elem_type == 502)
    {
      TET_T::Tet4 * tetcell = new TET_T::Tet4();

      for(int ee{0}; ee < num_cell[0]; ++ee)
      {
        const int node_t[3] { get_ien(0, ee, 0), get_ien(0, ee, 1), get_ien(0, ee, 2) };

        const std::array<int,3> node_t_gi {{ get_global_node(0, node_t[0]),
                                             get_global_node(0, node_t[1]),
                                             get_global_node(0, node_t[2]) }};

        const int cell_gi = get_global_cell(0, ee);

        const std::array<int,4> tet_n {{ VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                             VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3) }};
        
        tetcell->reset( tet_n );

        face_id[ee] = tetcell->get_face_id( node_t_gi );
      }

      delete tetcell;
    }
    else if(elem_type == 601 || elem_type == 602)
    {
      HEX_T::Hex8 * hexcell = new HEX_T::Hex8();

      for(int ee{0}; ee < num_cell[0]; ++ee)
      {
        const int node_q[4] { get_ien(0, ee, 0), get_ien(0, ee, 1),
                              get_ien(0, ee, 2), get_ien(0, ee, 3) };
        
        const std::array<int,4> node_q_gi {{ get_global_node(0, node_q[0]),
                                             get_global_node(0, node_q[1]),
                                             get_global_node(0, node_q[2]),
                                             get_global_node(0, node_q[3]) }};
        
        const int cell_gi = get_global_cell(0, ee);

        const std::array<int,8> hex_n {{ VIEN->get_IEN(cell_gi, 0), VIEN->get_IEN(cell_gi, 1),
                                         VIEN->get_IEN(cell_gi, 2), VIEN->get_IEN(cell_gi, 3),
                                         VIEN->get_IEN(cell_gi, 4), VIEN->get_IEN(cell_gi, 5),
                                         VIEN->get_IEN(cell_gi, 6), VIEN->get_IEN(cell_gi, 7) }};
        
        hexcell->reset( hex_n );

        face_id[ee] = hexcell->get_face_id( node_q_gi );
      }

      delete hexcell;
    }
    else
      SYS_T::print_fatal("Error: ElemBC_3D_turbulence_wall_model, unknown element type.\n");
  }
  else
    ; // The weak_list is empty and the wall file was put in the dir_list. Strongly enforced Dirichlet BC will be appiled.
}

// EOF
