/*******************************************************************************
 *
 *  InitializeMatrix.cpp
 *
 *  written by:     Sang Ik Cho
 *                  stc142@psu.edu 
 *  last modified:  7 / 22 / 2011
 *     
 *    This module initializes the 3-D computational domain with an incoming
 *  sonic boom for LoBoomFDTD3D_MultiBldg program.  Both the incident wave and 
 *  its reflection form the ground are included.  The boundary constants for PML 
 *  boundaries are also calculated.
 *               
 *******************************************************************************/

#include "Headers.h"

void InitializeMatrix(Grid3D &U, Array &sigma, Array &boom_shape)
{
  int i,j,k,prop_time_inc,prop_time_ref;
  float a_inc,b_inc,c_inc,d_inc,a_ref,b_ref,c_ref,d_ref;
  int input_count = 0;
  
  // Calculate PML boundary constants
  for (int i=0; i<D; i++)
    sigma[i] = sigma_m*pow(((D-i)/(float)D),B);

  // Read the input boom signiture
  ifstream f_input;
  f_input.open(str_input_name.c_str());
  if(!f_input) {
    cout << "ERROR: The program could not locate an input file named: " << str_input_name << endl;
    exit(1);
  } 
  else {
    while(!f_input.eof()) {
      f_input >> scientific >> boom_shape[input_count++];
      if (input_count>input_length_i+1) {
        cout << "ERROR: The input file is longer than specified by input_length_i." << endl;
        exit(1);
      }
    }
    f_input.close();
  }

  // Fill in the domain with incident & reflected waves
  a_inc = -cos(theta);
  b_inc = 0;
  c_inc = -sin(theta);
  d_inc = -a_inc*i_intersect;
  a_ref = -cos(theta);
  b_ref = 0;
  c_ref = sin(theta);
  d_ref = -a_ref*i_intersect;
  
  for (k=1; k<K_MAX+1; k++) {
    for (j=1; j<J_MAX+1; j++) {
      for (i=1; i<I_MAX+1; i++) {
        if (U[i][j][k].InOrOut==1) {
          U[i][j][k].px = 100;
          U[i][j][k].py = 100;
          U[i][j][k].pz = 100;
        }
        else {
          // Time for the planewave wavefront to get to a point, measured in time index
          prop_time_inc = -(int)floor(fs*dx*(a_inc*i+b_inc*j+c_inc*k+d_inc)/c0);
          prop_time_ref = -(int)floor(fs*dx*(a_ref*i+b_ref*j+c_ref*k+d_ref)/c0);
  
          // Orient the incident boom in the 3-D domain
          if (prop_time_inc>=0 && prop_time_inc<input_count) {
            U[i][j][k].px = boom_shape[prop_time_inc]/3;
            U[i][j][k].py = boom_shape[prop_time_inc]/3;
            U[i][j][k].pz = boom_shape[prop_time_inc]/3;
            U[i][j][k].u = boom_shape[prop_time_inc]/z0*a_inc;
            U[i][j][k].v = boom_shape[prop_time_inc]/z0*b_inc;
            U[i][j][k].w = boom_shape[prop_time_inc]/z0*c_inc;
          }
          else {
            U[i][j][k].px = 0;
            U[i][j][k].py = 0;
            U[i][j][k].pz = 0;
            U[i][j][k].u = 0;
            U[i][j][k].v = 0;
            U[i][j][k].w = 0;
          }
          // Orient the reflected boom in the 3-D domain
          if (prop_time_ref>=0 && prop_time_ref<input_count) {
            U[i][j][k].px = U[i][j][k].px+boom_shape[prop_time_ref]/3;
            U[i][j][k].py = U[i][j][k].py+boom_shape[prop_time_ref]/3;
            U[i][j][k].pz = U[i][j][k].pz+boom_shape[prop_time_ref]/3;
            U[i][j][k].u = U[i][j][k].u+boom_shape[prop_time_ref]/z0*a_ref;
            U[i][j][k].v = U[i][j][k].v+boom_shape[prop_time_ref]/z0*b_ref;
            U[i][j][k].w = U[i][j][k].w+boom_shape[prop_time_ref]/z0*c_ref;
          }
        }
        // Store the initial pressure field to an output file
        p_3D_init << U[i][j][k].px+U[i][j][k].py+U[i][j][k].pz << endl;
      }
    }
  }

  p_3D_init.close();
  
  return;
}
