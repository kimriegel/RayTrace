/*******************************************************************************
 *
 *  ReadParameter_mod.cpp
 *
 *  written by:     Sang Ik Cho
 *                  stc142@psu.edu 
 *  last modified:  7 / 22 / 2011
 *     
 *    This module reads in the input files "In_parameters_sim.txt", 
 *  "In_parameters_ani.txt" and "In_parameters_dom.txt" and assigns the values
 *  to the appropriate parameters for LoBoomFDTD2D_2006 program.  
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
  f_par_sim >> str_input_name;  // Input boom signature file
  f_par_sim.ignore(256,'=');
  f_par_sim >> fs;              // Sampling frequency of the input wave time series [Hz]
  f_par_sim.ignore(256,'=');
  f_par_sim >> input_length_i;  // Length of the input boom file [samples]
  f_par_sim.ignore(256,'=');
  f_par_sim >> boom_length_i;   // Length of the boom portion of the input boom file [samples]
  f_par_sim.ignore(256,'=');
  f_par_sim >> rho;             // Volumetric density of the medium [kg/m^3]
  f_par_sim.ignore(256,'=');
  f_par_sim >> c0;              // Speed of sound in the medium [m/s]
  f_par_sim.ignore(256,'=');
  f_par_sim >> dx;              // Spatial spacing of a uniform grid [m]
  f_par_sim.ignore(256,'=');
  f_par_sim >> CFL;             // Courant-Fredrich-Lewy number
  f_par_sim.ignore(256,'=');
  f_par_sim >> t_max;           // Time elapsed from the beginning to the end of the simulation [s]
  f_par_sim.ignore(256,'=');
  f_par_sim >> sigma_m;         // Damping coefficient for PML boundaries
  f_par_sim.ignore(256,'=');
  f_par_sim >> B;               // Exponential factor for setting damping constants for PML boundaries
  f_par_sim.ignore(256,'=');
  f_par_sim >> D;               // Thickness of PML boundaries
  f_par_sim.ignore(256,'=');
  f_par_sim >> angle_elv;       // Elevation angle of the incident boom direction [deg]
  f_par_sim.ignore(256,'=');
  f_par_sim >> angle_azm;       // Azimuthal angle of the incident boom direction, measured from True North [deg]
  f_par_sim.close();

  // Simulation variables specification & calculation
  z0 = rho*c0;                      // Characteristic impedance of the medium [Pa-s/m]
  dt = CFL*dx/c0;                   // Timestep size for the simulation [s]
  n_max = (int)round(t_max/dt);     // Length of the simulation in timesteps
  theta = angle_elv/180*PI;         // Elevation angle of the incident boom direction [rad]
  phi = (270-angle_azm)/180*PI;     // Azimuthal angle of the incident boom direction, measured from True North [rad]

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
  f_par_ani >> num_frame_yz;   // Total number of frames for YZ slice visualization 
  f_par_ani.ignore(256,'=');
  f_par_ani >> x_slice_yz;     // x-coordinate of the slice [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> y_range_yz1;    // Starting y-coordinate of YZ slice [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> y_range_yz2;    // Ending y-coordinate of YZ slice [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> z_range_yz1;    // Starting z-coordinate of YZ slice [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> z_range_yz2;    // Ending z-coordinate of YZ slice [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> num_frame_xz;   // Total number of frames for XZ slice visualization
  f_par_ani.ignore(256,'=');
  f_par_ani >> y_slice_xz;     // y-coordinate of the slice [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> x_range_xz1;    // Starting x-coordinate of XZ slice [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> x_range_xz2;    // Ending x-coordinate of XZ slice [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> z_range_xz1;    // Starting z-coordinate of XZ slice [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> z_range_xz2;    // Ending z-coordinate of XZ slice [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> num_frame_xy;   // Total number of frames for XY slice visualization
  f_par_ani.ignore(256,'=');
  f_par_ani >> z_slice_xy;     // Distance from the ground in z-direction to be cut [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> x_range_xy1;    // Starting x-coordinate of XY slice [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> x_range_xy2;    // Ending x-coordinate of XY slice [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> y_range_xy1;    // Starting y-coordinate of XY slice [m]
  f_par_ani.ignore(256,'=');
  f_par_ani >> y_range_xy2;    // Ending y-coordinate of XY slice [m]
  f_par_ani.close();
  
  // Animation & Visualization Parameter Specification & Variable Calculation   
  i_slice_yz = int(round(x_slice_yz/dx))+D;
  j_range_yz1 = int(round(y_range_yz1/dx))+D;
  j_range_yz2 = int(round(y_range_yz2/dx))+D;
  k_range_yz1 = int(round(z_range_yz1/dx));
  k_range_yz2 = int(round(z_range_yz2/dx));
  j_slice_xz = int(round(y_slice_xz/dx))+D;
  i_range_xz1 = int(round(x_range_xz1/dx))+D;
  i_range_xz2 = int(round(x_range_xz2/dx))+D;
  k_range_xz1 = int(round(z_range_xz1/dx));
  k_range_xz2 = int(round(z_range_xz2/dx));
  k_slice_xy = int(round(z_slice_xy/dx))+1;
  i_range_xy1 = int(round(x_range_xy1/dx))+D;
  i_range_xy2 = int(round(x_range_xy2/dx))+D;
  j_range_xy1 = int(round(y_range_xy1/dx))+D;
  j_range_xy2 = int(round(y_range_xy2/dx))+D;
  n_skip_yz = n_max/num_frame_yz; // Number of time steps to skip between each movie frame 
  n_skip_xz = n_max/num_frame_xz; // Number of time steps to skip between each movie frame
  n_skip_xy = n_max/num_frame_xy; // Number of time steps to skip between each movie frame
  
  // Preventing "floating point exception" and wrong number of frames for short simulations
  if (n_skip_yz == 0) {
    n_skip_yz = 1;
    num_frame_yz = n_max;
  }
  if (n_skip_xz == 0) {
    n_skip_xz = 1;
    num_frame_xz = n_max;
  }
  if (n_skip_xy == 0) {
    n_skip_xy = 1;
    num_frame_xy = n_max;
  }


  /*****************************************************************************
   *  Computational Domain Parameters
   *****************************************************************************/
  // Read parameters from input file
  ifstream f_par_dom("In_parameters_dom.txt");
  if (!f_par_dom) {
    cout << "ERROR: Could not open In_parameters_dom.txt" <<endl;
    exit(1);
  }
  f_par_dom.ignore(256,'=');
  f_par_dom >> extra_dim_x;
  f_par_dom.ignore(256,'=');
  f_par_dom >> extra_dim_y;
  f_par_dom.ignore(256,'=');
  f_par_dom >> buffer_dim_x;
  f_par_dom.ignore(256,'=');
  f_par_dom >> buffer_dim_z;
  f_par_dom.ignore(256,'=');
  f_par_dom >> wall_height;
  f_par_dom.ignore(256,'=');
  f_par_dom >> wall_thickness;  
  f_par_dom.close();
  

  return;
}
