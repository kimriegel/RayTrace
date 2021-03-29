/*******************************************************************************
 *
 *  InitializeDataFiles.cpp
 *
 *  written by:     Sang Ik Cho
 *                  stc142@psu.edu 
 *  last modified:  7 / 22 / 2011
 *     
 *    This module initializes the output files of LoBoomFDTD3D_MultiBldg 
 *  program with their proper headings to be read into Tecplot.  The output 
 *  files include the pressure field over the entire computational domain at the
 *  initial stage of the simulation, the snapshots of pressure in three 
 *  orthogonal slices over the duration of the simulation, and the pressure time
 *  history at the microphone locations.  The numerical precision for 
 *  pressure field data should be set low in order to keep the data file as 
 *  small as possible, but not so low as to lose the detail in the contour.  The 
 *  precision for the pressure time history data is set to be relatively high to
 *  maintain the fine detail of the waveforms.
 *    
 *******************************************************************************/     

#include "Headers.h"

void InitializeDataFiles()
{
  p_3D_init.open("p_3D_initial.dat");
  p_2DYZ.open("p_2D_slices_YZ.dat");
  p_2DXZ.open("p_2D_slices_XZ.dat");
  p_2DXY.open("p_2D_slices_XY.dat");
  p_mic.open("p_mics.dat");
  
  // Initial pressure field over the entire computational domain
  p_3D_init<<"TITLE = \"Initial Pressure in 3D\""<<endl;
  p_3D_init<<"VARIABLES = \"P\""<<endl;
  p_3D_init<<"ZONE DATAPACKING=BLOCK I="<<I_MAX<<", J="<<J_MAX<<", K="<<K_MAX<<endl;
  p_3D_init.setf(ios::fixed,ios::floatfield);
  p_3D_init.precision(2);

  // Pressure snapshots in YZ plane over the duration of the simulation
  p_2DYZ<<"TITLE = \"Propagating Pressure Pulse, in YZ plane\""<<endl;
  p_2DYZ<<"VARIABLES = \"P\""<<endl;
  p_2DYZ<<"ZONE DATAPACKING=BLOCK I="<<(j_range_yz2-j_range_yz1)<<", J="<<(k_range_yz2-k_range_yz1)<<", K="<<num_frame_yz<<endl;
  p_2DYZ.setf(ios::fixed,ios::floatfield);
  p_2DYZ.precision(2);

  // Pressure snapshots in XZ plane over the duration of the simulation
  p_2DXZ<<"TITLE = \"Propagating Pressure Pulse, in XZ plane\""<<endl;
  p_2DXZ<<"VARIABLES = \"P\""<<endl;
  p_2DXZ<<"ZONE DATAPACKING=BLOCK I="<<(i_range_xz2-i_range_xz1)<<", J="<<(k_range_xz2-k_range_xz1)<<", K="<<num_frame_xz<<endl;
  p_2DXZ.setf(ios::fixed,ios::floatfield);
  p_2DXZ.precision(2);

  // Pressure snapshots in XY plane over the duration of the simulation
  p_2DXY<<"TITLE = \"Propagating Pressure Pulse, in XY plane\""<<endl;
  p_2DXY<<"VARIABLES = \"P\""<<endl;
  p_2DXY<<"ZONE DATAPACKING=BLOCK I="<<(i_range_xy2-i_range_xy1)<<", J="<<(j_range_xy2-j_range_xy1)<<", K="<<num_frame_xy<<endl;
  p_2DXY.setf(ios::fixed,ios::floatfield);
  p_2DXY.precision(2);
  
  // Pressure time history at microphone locations
  p_mic.setf(ios::fixed,ios::floatfield);
  p_mic.precision(5);

  return;
}
