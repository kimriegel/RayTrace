/*******************************************************************************
 *
 *  ReadParameter_mod.cpp
 *
 *  written by:     Sang Ik Cho
 *                  stc142@psu.edu 
 *  last modified:  7 / 15 / 2011
 *     
 *    This module reads in the input files "In_parameters_sim.txt", 
 *  "In_parameters_ani.txt" and "In_parameters_mic.txt" and assigns the values
 *  to the appropriate parameters for LoBoomFDTD2D_2006_simple program.  
 *  
 *******************************************************************************/

#include "Headers.h"
  
void ReadParameters()
{
  /*****************************************************************************
   *  Simulation Parameters
   *****************************************************************************/
  // Read parameters from input file
  ifstream f_par_sim("In_parameters_sim.txt");
  if (!f_par_sim) {
    cout << "ERROR: Could not open In_parameters_sim.txt" <<endl;
    exit(1);
  }
  f_par_sim.ignore(256,'=');
  f_par_sim >> input_name;      // Input boom signature file
  f_par_sim.ignore(256,'=');
  f_par_sim >> input_length_i;  // Length of the input boom file [samples]
  f_par_sim.ignore(256,'=');
  f_par_sim >> input_fs;        // Sampling frequency of the input wave time series [Hz]
  f_par_sim.ignore(256,'=');
  f_par_sim >> c0;              // Speed of sound in the medium [m/s]
  f_par_sim.ignore(256,'=');
  f_par_sim >> theta_inc;       // Elevation angle of the incident boom direction [deg]
  f_par_sim.ignore(256,'=');
  f_par_sim >> num_floors;      // Number of floors of the house
  f_par_sim.ignore(256,'=');
  f_par_sim >> x_max;           // Domain size in horizontal direction[m]
  f_par_sim.ignore(256,'=');
  f_par_sim >> z_max;           // Domain size in vertical direction[m]
  f_par_sim.ignore(256,'=');
  f_par_sim >> x_intersect;     // x-coordinate of the initial wavefront intersection with the ground [m]
  f_par_sim.ignore(256,'=');
  f_par_sim >> x0;              // x-coordinate of the building wall on the incident side [m]
  f_par_sim.ignore(256,'=');
  f_par_sim >> dx;              // Grid spacing of a uniform grid [m]
  f_par_sim.ignore(256,'=');
  f_par_sim >> tmax;            // Time elapsed from the beginning to the end of the simulation [s]
  f_par_sim.close();

  // Simulation variables specification & calculation
  x_indmax = (int)(x_max/dx);
  z_indmax = (int)(z_max/dx);  
  dt = dx/(1.415*c0);
  t_indmax = (int)(tmax/dt);
  cXcXdtXdt = c0*c0*dt*dt;
  dxXdx = dx*dx;
  theta_ref = -theta_inc;
  tan14 = tan(14./180*PI);
  z1 = -(x0+4.4229)*tan14+3.2865+(num_floors-1)*2.5;
  z2 = -(x0+4.4229)*tan14+3.4389+(num_floors-1)*2.5;
  z3 = (x0+4.4229)*tan14+3.2865+(num_floors-1)*2.5;
  z4 = (x0+4.4229)*tan14+3.4389+(num_floors-1)*2.5;


  /*****************************************************************************
   *  Animation Visualization Parameters
   *****************************************************************************/
  // Read parameters from input file
  ifstream f_par_ani("In_parameters_ani.txt");
  if (!f_par_ani) {
    cout << "ERROR: Could not open In_parameters_ani.txt" <<endl;
    exit(1);
  }
  f_par_ani.ignore(256,'=');
  f_par_ani >> num_frame;     // Total number of frames for animation
  f_par_ani.ignore(256,'=');
  f_par_ani >> x_range_lo;    // Starting x-coordinate of animation snapshot [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> x_range_hi;    // Ending x-coordinate of animation snapshot [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> z_range_hi;    // Starting z-coordinate of animation snapshot [m]
  f_par_ani.close();
  
  
  // Animation & Visualization Parameter Specification & Variable Calculation   
  x_range_ind_lo = (int)(x_range_lo/dx);       
  x_range_ind_hi = (int)(x_range_hi/dx);
  z_range_ind_hi = (int)(z_range_hi/dx);

  mov_skip = t_indmax/num_frame; // Number of time steps to skip between each animation frame
  
  // Preventing "floating point exception" and wrong number of frames for short simulations
  if (mov_skip == 0) {
    mov_skip = 1;
    num_frame = t_indmax;
  }


  /*****************************************************************************
   *  Microphone selection
   *****************************************************************************/
  ifstream f_par_mic("In_parameters_mic.txt");
  if (!f_par_mic) {
    cout << "ERROR: Could not open In_parameters_mic.txt" <<endl;
    exit(1);
  }
  Mic mic_temp;
  string line;
  while (f_par_mic.good()) {
    getline (f_par_mic,line);
    if (line=="196") {
      mic_temp.i = (int)ceil(x0/dx)-1;
      mic_temp.j = 1;
      mics.push_back(mic_temp);
    }
    else if (line=="197") {
      mic_temp.i = (int)ceil(x0/dx)-1;
      mic_temp.j = (int)floor(1.2192/dx);
      mics.push_back(mic_temp);
    }
    else if (line=="179") {
      mic_temp.i = (int)floor((x0-0.7732)/dx);
      mic_temp.j = (int)ceil((2.1796+(num_floors-1)*2.5)/dx);
      mics.push_back(mic_temp);
    }
    else if (line=="180") {
      mic_temp.i = (int)floor((x0+1.6662)/dx);
      mic_temp.j = (int)ceil(2.7878/dx);
      mics.push_back(mic_temp);
    }
    else if (line=="181") {
      mic_temp.i = (int)floor((x0+4.2782)/dx);
      mic_temp.j = (int)ceil((3.4389+(num_floors-1)*2.5)/dx);
      mics.push_back(mic_temp);
    }
    else if (line=="182") {
      mic_temp.i = (int)ceil((x0+7.1054)/dx);
      mic_temp.j = (int)ceil(2.8062/dx);
      mics.push_back(mic_temp);
    }
    else if (line=="183") {
      mic_temp.i = (int)ceil((x0+9.6558)/dx);
      mic_temp.j = (int)ceil((2.1703+(num_floors-1)*2.5)/dx);
      mics.push_back(mic_temp);
    }
    else if (line=="184") {
      mic_temp.i = (int)ceil((x0+8.8458)/dx);
      mic_temp.j = (int)floor(1.2192/dx);
      mics.push_back(mic_temp);
    }
    else if (line=="198") {
      mic_temp.i = (int)ceil((x0+8.8458)/dx);
      mic_temp.j = 1;
      mics.push_back(mic_temp);
    }
  }
  f_par_mic.close();
  num_mics = mics.size();   // Total number of microphones recording
  

  return;
}
