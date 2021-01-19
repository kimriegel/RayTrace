/*******************************************************************************
 *
 *  VisualizeMatrix.cpp
 *
 *  written by:     Sang Ik Cho
 *                  stc142@psu.edu 
 *  last modified:  7 / 22 / 2011
 *     
 *    This module saves the pressure data to output files of 
 *  LoBoomFDTD3D_MultiBldg program.  The snapshots of pressure field are saved 
 *  in a format ready to be post-processed on Tecplot and time history of 
 *  pressure is also recorded in columns of numeric values to be processed on 
 *  Matlab or other software.
 *    
 *******************************************************************************/     

#include "Headers.h"

void VisualizeMatrix(Grid3D &U, int n)
{
  int i, j, k, l;
  
  // Highlight the building interior points for visualization
  for (i=1; i<=I_MAX+1; i++) {
    for (j=1; j<=J_MAX+1; j++) {
      for (k=1; k<=K_MAX+1; k++) {
        if (U[i][j][k].InOrOut==1) {
          U[i][j][k].px = 100;
          U[i][j][k].py = 100;
          U[i][j][k].pz = 100;
        }
      }
    }
  }

  // Snapshots of pressure field recorded on slices in three directions
  if (n%n_skip_yz==0)
    for (k=k_range_yz1; k<k_range_yz2; k++)
      for (j=j_range_yz1; j<j_range_yz2; j++)
        p_2DYZ<<U[i_slice_yz][j][k].px+U[i_slice_yz][j][k].py+U[i_slice_yz][j][k].pz<<endl;
        
  if (n%n_skip_xz==0)
    for (k=k_range_xz1; k<k_range_xz2; k++)
      for (i=i_range_xz1; i<i_range_xz2; i++)
        p_2DXZ<<U[i][j_slice_xz][k].px+U[i][j_slice_xz][k].py+U[i][j_slice_xz][k].pz<<endl;
        
  if (n%n_skip_xy==0)
    for (j=j_range_xy1; j<j_range_xy2; j++)
      for (i=i_range_xy1; i<i_range_xy2; i++)
        p_2DXY<<U[i][j][k_slice_xy].px+U[i][j][k_slice_xy].py+U[i][j][k_slice_xy].pz<<endl;

  // Time series of pressure recorded at Microphone locations
  p_mic << n*dt << " ";
  for (l=0; l<num_mics; l++)
    p_mic << U[mics[l].i][mics[l].j][mics[l].k].px+U[mics[l].i][mics[l].j][mics[l].k].py+U[mics[l].i][mics[l].j][mics[l].k].pz << " ";
  p_mic << endl;
  
  return;
}
