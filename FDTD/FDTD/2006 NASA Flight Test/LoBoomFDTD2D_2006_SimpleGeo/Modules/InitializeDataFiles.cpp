/*******************************************************************************
 *
 *  InitializeDataFiles.cpp
 *
 *  written by:     Sang Ik Cho
 *                  stc142@psu.edu 
 *  last modified:  7 / 15 / 2011
 *     
 *    This module initializes the output files of LoBoomFDTD2D_2006_simple
 *  program with their proper headings to be read into Tecplot.  The output 
 *  files include the pressure field over the entire computational domain at the 
 *  initial and final stages of the simulation, the snapshots of pressure in a 
 *  selected region of domain over the duration of the simulation, and the 
 *  pressure time history at the microphone locations.  The numerical precision 
 *  for pressure field data should be set low in order to keep the data file as 
 *  small as possible, but not so low as to lose the detail in the contour.  
 *  The precision for the pressure time history data is set to be relatively 
 *  high to maintain the fine detail of the waveforms.     
 *    
 *******************************************************************************/     

#include "Headers.h"

void InitializeDataFiles()
{
  int x_indrange = x_range_ind_hi-x_range_ind_lo;
  
  p_2D_initial.open("p_2D_initial.dat");
  p_2D_final.open("p_2D_final.dat");
  p_domain.open("p_domain.dat");
  p_mic.open("p_mics.dat");
  
  // Initial pressure field over the entire computational domain
  p_2D_initial << "TITLE = \"Initial Pressure\"" << endl;
  p_2D_initial << "VARIABLES = \"P\"" << endl ;
  p_2D_initial << "ZONE DATAPACKING=BLOCK I=" << x_indmax << ", J=" << z_indmax << endl;
  
  // Final pressure field over the entire computational domain
  p_2D_final << "TITLE = \"Final Pressure\"" << endl;
  p_2D_final << "VARIABLES = \"P\"" << endl;
  p_2D_final << "ZONE DATAPACKING=BLOCK I=" << x_indmax << ", J=" << z_indmax << endl;
  
  // Pressure snapshots over the duration of the simulation
  p_domain << "TITLE = \"Propagating Pressure Pulse\"" << endl;
  p_domain << "VARIABLES = \"P\"" << endl;
  p_domain << "ZONE DATAPACKING=BLOCK I=" << x_indrange << ", J=" << z_range_ind_hi << ", K=" << num_frame << endl;
  
  // Pressure time history at microphone locations
  p_mic.setf(ios::fixed,ios::floatfield);
  p_mic.precision(5);

  return;
}
