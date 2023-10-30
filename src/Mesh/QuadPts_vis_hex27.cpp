#include "QuadPts_vis_hex27.hpp"

QuadPts_vis_hex27::QuadPts_vis_hex27()
{
  qp[0 ] = 0.0; qp[1 ] = 0.0; qp[2 ] = 0.0; 
  qp[3 ] = 1.0; qp[4 ] = 0.0; qp[5 ] = 0.0; 
  qp[6 ] = 1.0; qp[7 ] = 1.0; qp[8 ] = 0.0;
  qp[9 ] = 0.0; qp[10] = 1.0; qp[11] = 0.0; 
  qp[12] = 0.0; qp[13] = 0.0; qp[14] = 1.0;
  qp[15] = 1.0; qp[16] = 0.0; qp[17] = 1.0;
  qp[18] = 1.0; qp[19] = 1.0; qp[20] = 1.0;
  qp[21] = 0.0; qp[22] = 1.0; qp[23] = 1.0;
  qp[24] = 0.5; qp[25] = 0.0; qp[26] = 0.0; 
  qp[27] = 1.0; qp[28] = 0.5; qp[29] = 0.0; 
  qp[30] = 0.5; qp[31] = 1.0; qp[32] = 0.0;
  qp[33] = 0.0; qp[34] = 0.5; qp[35] = 0.0; 
  qp[36] = 0.5; qp[37] = 0.0; qp[38] = 1.0;
  qp[39] = 1.0; qp[40] = 0.5; qp[41] = 1.0;
  qp[42] = 0.5; qp[43] = 1.0; qp[44] = 1.0;
  qp[45] = 0.0; qp[46] = 0.5; qp[47] = 1.0;
  qp[48] = 0.0; qp[49] = 0.0; qp[50] = 0.5; 
  qp[51] = 1.0; qp[52] = 0.0; qp[53] = 0.5; 
  qp[54] = 1.0; qp[55] = 1.0; qp[56] = 0.5;
  qp[57] = 0.0; qp[58] = 1.0; qp[59] = 0.5; 
  qp[60] = 0.0; qp[61] = 0.5; qp[62] = 0.5;
  qp[63] = 1.0; qp[64] = 0.5; qp[65] = 0.5;
  qp[66] = 0.5; qp[67] = 0.0; qp[68] = 0.5;
  qp[69] = 0.5; qp[70] = 1.0; qp[71] = 0.5; 
  qp[72] = 0.5; qp[73] = 0.5; qp[74] = 0.0; 
  qp[75] = 0.5; qp[76] = 0.5; qp[77] = 1.0; 
  qp[78] = 0.5; qp[79] = 0.5; qp[80] = 0.5;

  qw[0 ] = 0.5; qw[1 ] = 0.5;
  qw[2 ] = 0.5; qw[3 ] = 0.5;
  qw[4 ] = 0.5; qw[5 ] = 0.5; 
  qw[6 ] = 0.5; qw[7 ] = 0.5;

  qw[8 ] = 0.5; qw[9 ] = 0.5;
  qw[10] = 0.5; qw[11] = 0.5;
  qw[12] = 0.5; qw[13] = 0.5; 
  qw[14] = 0.5; qw[15] = 0.5;

  qw[16] = 0.5; qw[17] = 0.5;
  qw[18] = 0.5; qw[19] = 0.5;
  qw[20] = 0.5; qw[21] = 0.5; 
  qw[22] = 0.5; qw[23] = 0.5;

  qw[24] = 0.5; qw[25] = 0.5;
  qw[26] = 0.5;
}

QuadPts_vis_hex27::~QuadPts_vis_hex27()
{}

void QuadPts_vis_hex27::print_info() const
{
  SYS_T::commPrint("\n===== Visualization Points for Hex27 ===== \n");
  for(int ii=0; ii<27; ++ii)
    SYS_T::commPrint("%e, %e, %e, %e \n", 
        qw[ii], qp[3*ii], qp[3*ii+1], qp[3*ii+2]);
  SYS_T::commPrint("========================================= \n");
}

// EOF
