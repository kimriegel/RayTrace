/*******************************************************************************
 *
 *  VisualizeMatrix.cpp
 *
 *  written by:     Sang Ik Cho
 *                  stc142@psu.edu 
 *  last modified:  7 / 22 / 2011
 *     
 *    This module saves the pressure data to output files of 
 *  LoBoomFDTD3D_SingleBldg program.  The snapshots of pressure field are saved 
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
  for (i=wall_loc_i1; i<=wall_loc_i2; i++) {
    for (j=wall_loc_j1; j<=wall_loc_j2; j++) {
      for (k=wall_loc_k1; k<=wall_loc_k2; k++) {
        U[i][j][k].px = 100;
        U[i][j][k].py = 100;
        U[i][j][k].pz = 100;
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
  
  // Time history of pressure recorded at microphone locations specified by user
  p_mic << n*dt << " ";
  for (l=0; l<num_mics; l++)
    p_mic << U[mics[l].i][mics[l].j][mics[l].k].px+U[mics[l].i][mics[l].j][mics[l].k].py+U[mics[l].i][mics[l].j][mics[l].k].pz << " ";
  p_mic << endl;
  
  // Pressure distribution on exterior surfaces
  if (n%n_skip_surface==0) {
    // Roof
    p_surface<<"ZONE T=\"Roof\"\nSTRANDID=0, SOLUTIONTIME="<<n*dt<<endl;
    p_surface<<"I="<<wall_loc_i2-wall_loc_i1+3<<", J="<<wall_loc_j2-wall_loc_j1+3<<", K=1, ZONETYPE=Ordered\nDATAPACKING=POINT"<<endl;
    p_surface<<"DT=(SINGLE SINGLE SINGLE SINGLE )"<<endl;
    for (int j=wall_loc_j1-1;j<=wall_loc_j2+1;j++)
      for (int i=wall_loc_i1-1;i<=wall_loc_i2+1;i++)
        p_surface<<(i-(int)(wall_loc_i1+wall_loc_i2)/2)*dx<<"\t"<<(j-(int)(wall_loc_j1+wall_loc_j2)/2)*dx<<"\t"<<(wall_loc_k2+1)*dx<<"\t"<<U[i][j][wall_loc_k2+1].px+U[i][j][wall_loc_k2+1].py+U[i][j][wall_loc_k2+1].pz<<endl;
    // Front Wall
    p_surface<<"ZONE T=\"Front Wall\"\nSTRANDID=0, SOLUTIONTIME="<<n*dt<<endl;
    p_surface<<"I=1, J="<<wall_loc_j2-wall_loc_j1+3<<", K="<<wall_loc_k2-wall_loc_k1+3<<", ZONETYPE=Ordered\nDATAPACKING=POINT"<<endl;
    p_surface<<"DT=(SINGLE SINGLE SINGLE SINGLE )"<<endl;
    for (int k=wall_loc_k1-1;k<=wall_loc_k2+1;k++)
      for (int j=wall_loc_j1-1;j<=wall_loc_j2+1;j++)
        p_surface<<(wall_loc_i1-(int)(wall_loc_i1+wall_loc_i2)/2-1)*dx<<"\t"<<(j-(int)(wall_loc_j1+wall_loc_j2)/2)*dx<<"\t"<<k*dx<<"\t"<<U[wall_loc_i1-1][j][k].px+U[wall_loc_i1-1][j][k].py+U[wall_loc_i1-1][j][k].pz<<endl;
    // Back Wall
    p_surface<<"ZONE T=\"Back Wall\"\nSTRANDID=0, SOLUTIONTIME="<<n*dt<<endl;
    p_surface<<"I=1, J="<<wall_loc_j2-wall_loc_j1+3<<", K="<<wall_loc_k2-wall_loc_k1+3<<", ZONETYPE=Ordered\nDATAPACKING=POINT"<<endl;
    p_surface<<"DT=(SINGLE SINGLE SINGLE SINGLE )"<<endl;
    for (int k=wall_loc_k1-1;k<=wall_loc_k2+1;k++)
      for (int j=wall_loc_j1-1;j<=wall_loc_j2+1;j++)
        p_surface<<(wall_loc_i2-(int)(wall_loc_i1+wall_loc_i2)/2+1)*dx<<"\t"<<(j-(int)(wall_loc_j1+wall_loc_j2)/2)*dx<<"\t"<<k*dx<<"\t"<<U[wall_loc_i2+1][j][k].px+U[wall_loc_i2+1][j][k].py+U[wall_loc_i2+1][j][k].pz<<endl;
    // Left Wall
    p_surface<<"ZONE T=\"Left Wall\"\nSTRANDID=0, SOLUTIONTIME="<<n*dt<<endl;
    p_surface<<"I="<<wall_loc_i2-wall_loc_i1+3<<", J=1, K="<<wall_loc_k2-wall_loc_k1+3<<", ZONETYPE=Ordered\nDATAPACKING=POINT"<<endl;
    p_surface<<"DT=(SINGLE SINGLE SINGLE SINGLE )"<<endl;
    for (int k=wall_loc_k1-1;k<=wall_loc_k2+1;k++)
      for (int i=wall_loc_i1-1;i<=wall_loc_i2+1;i++)
        p_surface<<(i-(int)(wall_loc_i1+wall_loc_i2)/2)*dx<<"\t"<<(wall_loc_j1-(int)(wall_loc_j1+wall_loc_j2)/2-1)*dx<<"\t"<<k*dx<<"\t"<<U[i][wall_loc_j1-1][k].px+U[i][wall_loc_j1-1][k].py+U[i][wall_loc_j1-1][k].pz<<endl;
    // Right Wall
    p_surface<<"ZONE T=\"Right Wall\"\nSTRANDID=0, SOLUTIONTIME="<<n*dt<<endl;
    p_surface<<"I="<<wall_loc_i2-wall_loc_i1+3<<", J=1, K="<<wall_loc_k2-wall_loc_k1+3<<", ZONETYPE=Ordered\nDATAPACKING=POINT"<<endl;
    p_surface<<"DT=(SINGLE SINGLE SINGLE SINGLE )"<<endl;
    for (int k=wall_loc_k1-1;k<=wall_loc_k2+1;k++)
      for (int i=wall_loc_i1-1;i<=wall_loc_i2+1;i++)
        p_surface<<(i-(int)(wall_loc_i1+wall_loc_i2)/2)*dx<<"\t"<<(wall_loc_j2-(int)(wall_loc_j1+wall_loc_j2)/2+1)*dx<<"\t"<<k*dx<<"\t"<<U[i][wall_loc_j2+1][k].px+U[i][wall_loc_j2+1][k].py+U[i][wall_loc_j2+1][k].pz<<endl;
  }

  return;
}
