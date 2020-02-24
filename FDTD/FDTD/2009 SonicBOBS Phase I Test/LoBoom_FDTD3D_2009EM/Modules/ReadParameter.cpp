/*******************************************************************************
 *
 *  ReadParameter.cpp
 *
 *  written by:     Sang Ik Cho
 *                  stc142@psu.edu 
 *  last modified:  9 / 28 / 2011
 *     
 *    This module reads in the input files "In_parameters_sim.txt" and 
 *  "In_parameters_vis.txt" and assigns the values to the appropriate parameters
 *  for LoBoomFDTD2D_2006 program.  
 *  
 *******************************************************************************/
  
#include "Headers.h"

#define PI 3.14159265

void ReadParameters()
{
  /*****************************************************************************
   *  Simulation Parameters
   *****************************************************************************/
  // Read parameters from input files
  ifstream f_par_sim("In_parameters_sim.txt");
  if (!f_par_sim) {
    cout << "ERROR: Could not open In_parameters_sim.txt" <<endl;
    exit(1);
  }
  f_par_sim.ignore(256,'=');
  f_par_sim >> str_input_name; // Input boom signature file
  f_par_sim.ignore(256,'=');
  f_par_sim >> fs;             // Sampling frequency of the input wave time series [Hz]
  f_par_sim.ignore(256,'=');
  f_par_sim >> input_length_i; // Length of the input boom file [samples]
  f_par_sim.ignore(256,'=');
  f_par_sim >> boom_length_i;  // Length of the boom portion of the input boom file [samples]
  f_par_sim.ignore(256,'=');
  f_par_sim >> rho;            // Volumetric density of the medium [kg/m^3]
  f_par_sim.ignore(256,'=');
  f_par_sim >> c0;             // Speed of sound in the medium [m/s]
  f_par_sim.ignore(256,'=');
  f_par_sim >> dx;             // Spatial spacing of a uniform grid [m]
  f_par_sim.ignore(256,'=');
  f_par_sim >> CFL;            // Courant-Fredrich-Lewy number
  f_par_sim.ignore(256,'=');
  f_par_sim >> t_max;          // Time elapsed from the beginning to the end of the simulation [s]
  f_par_sim.ignore(256,'=');
  f_par_sim >> sigma_m;        // Damping coefficient for PML boundaries
  f_par_sim.ignore(256,'=');
  f_par_sim >> B;              // Exponential factor for setting damping constants for PML boundaries
  f_par_sim.ignore(256,'=');
  f_par_sim >> D;              // Thickness of PML boundaries
  f_par_sim.ignore(256,'=');
  f_par_sim >> angle_elv;      // Elevation angle of the incident boom direction [deg]
  f_par_sim.ignore(256,'=');
  f_par_sim >> angle_azm;      // Azimuthal angle of the incident boom direction, measured from negative x-direction, aligned with the left building wall [deg]
  f_par_sim.ignore(256,'=');
  f_par_sim >> bldg_dim_x;     // Building dimension in x [m]
  f_par_sim.ignore(256,'=');
  f_par_sim >> bldg_dim_y;     // Building dimension in y [m]
  f_par_sim.ignore(256,'=');
  f_par_sim >> bldg_dim_z;     // Building dimension in z [m]
  f_par_sim.ignore(256,'=');
  f_par_sim >> bldg_offset_x;  // Building offset from the PML boundary in x [m]
  f_par_sim.ignore(256,'=');
  f_par_sim >> bldg_offset_y;  // Building offset from the PML boundary in y [m]
  f_par_sim.ignore(256,'=');
  f_par_sim >> buffer_dim_x;   // Extra space to avoid input wave getting distorted towards the back in x [m]
  f_par_sim.ignore(256,'=');
  f_par_sim >> buffer_dim_y;   // Extra space to avoid input wave getting distorted towards the back in y [m]
  f_par_sim.ignore(256,'=');
  f_par_sim >> buffer_dim_z;   // Extra space to avoid input wave getting distorted towards the back in z [m]
  f_par_sim.close();
  
  // Simulation variables specification & calculation
  z0 = rho*c0;                      // Characteristic impedance of the medium [Pa-s/m]
  dt = CFL*dx/c0;                   // Timestep size for the simulation [s]
  n_max = (int)round(t_max/dt);     // Length of the simulation in timesteps
  theta = angle_elv/180*PI;         // Elevation angle of the incoming wavefront [rad]
  phi = angle_azm/180*PI;           // Azimuthal angle of the incoming wavefront [rad]
  boom_duration = boom_length_i/fs;   // Time duration of the useful boom signature [s]
  x_intersect = bldg_dim_z*sin(theta)+bldg_dim_y*sin(phi)+bldg_dim_x; // x-coordinate of the initial wavefront intersection with the ground [m]
  i_intersect = (int)round((bldg_offset_x+x_intersect)/dx)+D;         // Index at x_intersect
  
  // Domain size calculation
  dom_dim_x = x_intersect+boom_duration*c0/(cos(phi)*cos(theta))+buffer_dim_x;  // Domain dimension in x [m]
  dom_dim_y = bldg_offset_y+bldg_dim_y+(boom_duration*c0/(cos(phi)*cos(theta))+bldg_dim_x)*cos(phi)*sin(phi)+buffer_dim_y;  // Domain dimension in y [m]
  dom_dim_z = bldg_dim_z+(boom_duration*c0/(cos(phi)*cos(theta))+bldg_dim_x)*cos(theta)*sin(theta)+buffer_dim_z; // Domain dimension in z [m]
  I_MAX = (int)round(dom_dim_x/dx)+2*D;           // Index maximum in x-direction
  J_MAX = (int)round(dom_dim_y/dx)+2*D;           // Index maximum in y-direction
  K_MAX = (int)round(dom_dim_z/dx)+2*D;           // Index maximum in z-direction

  // Building size and location
  bldg_dim_i = (int)round(bldg_dim_x/dx);         // Building dimension in index i
  bldg_dim_j = (int)round(bldg_dim_y/dx);         // Building dimension in index j
  bldg_dim_k = (int)round(bldg_dim_z/dx);         // Building dimension in index k
  bldg_offset_i = (int)round(bldg_offset_x/dx);   // Building offset from boundary in index i
  bldg_offset_j = (int)round(bldg_offset_y/dx);   // Building offset from boundary in index j
  wall_loc_i1 = bldg_offset_i+D;                  // West wall location index i
  wall_loc_i2 = wall_loc_i1+bldg_dim_i-1;         // East wall location index i
  wall_loc_j1 = bldg_offset_j+D;                  // South wall location index j
  wall_loc_j2 = wall_loc_j1+bldg_dim_j-1;         // North wall location index j
  wall_loc_k1 = 1;                                // "Floor" location index k (ground)
  wall_loc_k2 = wall_loc_k1+bldg_dim_k-1;         // Roof location index k
  

  /*****************************************************************************
   *  Visalization Parameters
   *****************************************************************************/
  // Read parameters from input files
  ifstream f_par_vis("In_parameters_vis.txt");
  if (!f_par_vis) {
    cout << "ERROR: Could not open In_parameters_vis.txt" <<endl;
    exit(1);
  }
  f_par_vis.ignore(256,'=');
  f_par_vis >> num_frame_yz;   // Total number of frames for YZ slice visualization 
  f_par_vis.ignore(256,'=');
  f_par_vis >> x_slice_yz;     // x-coordinate of the slice [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> y_range_yz1;    // Starting y-coordinate of YZ slice [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> y_range_yz2;    // Ending y-coordinate of YZ slice [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> z_range_yz1;    // Starting z-coordinate of YZ slice [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> z_range_yz2;    // Ending z-coordinate of YZ slice [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> num_frame_xz;   // Total number of frames for XZ slice visualization
  f_par_vis.ignore(256,'=');
  f_par_vis >> y_slice_xz;     // y-coordinate of the slice [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> x_range_xz1;    // Starting x-coordinate of XZ slice [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> x_range_xz2;    // Ending x-coordinate of XZ slice [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> z_range_xz1;    // Starting z-coordinate of XZ slice [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> z_range_xz2;    // Ending z-coordinate of XZ slice [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> num_frame_xy;   // Total number of frames for XY slice visualization
  f_par_vis.ignore(256,'=');
  f_par_vis >> z_slice_xy;     // Distance from the ground in z-direction to be cut [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> x_range_xy1;    // Starting x-coordinate of XY slice [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> x_range_xy2;    // Ending x-coordinate of XY slice [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> y_range_xy1;    // Starting y-coordinate of XY slice [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> y_range_xy2;    // Ending y-coordinate of XY slice [m]
  f_par_vis.ignore(256,'=');
  f_par_vis >> n_frame_surface;  // Number of frames for external pressure loading visualization
  f_par_vis.ignore(256,'=');
  f_par_vis >> num_mics;       // Number of microphones to record the pressure time history
  mics = new Mic[num_mics];
  float mic_x,mic_y,mic_z;
  for (int l=0; l<num_mics; l++) {
    f_par_vis.ignore(256,'=');
    f_par_vis >> mic_x;        // x-coordinate of a mic [m]
    f_par_vis.ignore(256,'=');
    f_par_vis >> mic_y;        // y-coordinate of a mic [m]
    f_par_vis.ignore(256,'=');
    f_par_vis >> mic_z;        // z-coordinate of a mic [m]
    mics[l].i = (int)round((mic_x+bldg_offset_x)/dx)+D;
    mics[l].j = (int)round((mic_y+bldg_offset_y)/dx)+D;
    mics[l].k = (int)floor(mic_z/dx)+1;
  }
  f_par_vis.close();

  // Adjust microphone location if it coincides with a building interior point
  for (int l=0; l<num_mics; l++) {
    if (mics[l].i==wall_loc_i1 && mics[l].j>=wall_loc_j1 && mics[l].j<=wall_loc_j2 && mics[l].k<=wall_loc_k2)
      mics[l].i--;
    if (mics[l].i==wall_loc_i2 && mics[l].j>=wall_loc_j1 && mics[l].j<=wall_loc_j2 && mics[l].k<=wall_loc_k2)
      mics[l].i++;
    if (mics[l].j==wall_loc_j1 && mics[l].i>=wall_loc_i1 && mics[l].i<=wall_loc_i2 && mics[l].k<=wall_loc_k2)
      mics[l].j--;
    if (mics[l].j==wall_loc_j2 && mics[l].i>=wall_loc_i1 && mics[l].i<=wall_loc_i2 && mics[l].k<=wall_loc_k2)
      mics[l].j++;
    if (mics[l].k==wall_loc_k2 && mics[l].i>=wall_loc_i1 && mics[l].i<=wall_loc_i2 && mics[l].j>=wall_loc_j1 && mics[l].j<=wall_loc_j2)
      mics[l].k++;
  }

  // Simulation & Visualization Parameter Specification & Variable Calculation   
  i_slice_yz = (int)round((x_slice_yz+bldg_offset_x)/dx)+D;
  if (i_slice_yz<D || i_slice_yz>I_MAX-D+1) {
    cout << "ERROR: Invalid x_slice_yz value" <<endl;
    exit(1);
  }
  j_range_yz1 = (int)round((y_range_yz1+bldg_offset_y)/dx)+D;
  if (j_range_yz1<0) {
    cout << "ERROR: Invalid y_range_yz1 value" <<endl;
    exit(1);
  }
  j_range_yz2 = (int)round((y_range_yz2+bldg_offset_y)/dx)+D;
  if (j_range_yz2>J_MAX+2 || j_range_yz2<j_range_yz1) {
    cout << "ERROR: Invalid y_range_yz2 value" <<endl;
    exit(1);
  }
  k_range_yz1 = (int)round(z_range_yz1/dx);
  if (k_range_yz1<0) {
    cout << "ERROR: Invalid z_range_yz1 value" <<endl;
    exit(1);
  }
  k_range_yz2 = (int)round(z_range_yz2/dx);
  if (k_range_yz2>K_MAX+2 || k_range_yz2<k_range_yz1) {
    cout << "ERROR: Invalid z_range_yz2 value" <<endl;
    exit(1);
  }
  j_slice_xz = (int)round((y_slice_xz+bldg_offset_y)/dx)+D;
  if (j_slice_xz<D || j_slice_xz>J_MAX-D+1) {
    cout << "ERROR: Invalid j_slice_xz value" <<endl;
    exit(1);
  }
  i_range_xz1 = (int)round((x_range_xz1+bldg_offset_x)/dx)+D;
  if (i_range_xz1<0) {
    cout << "ERROR: Invalid x_range_xz1 value" <<endl;
    exit(1);
  }
  i_range_xz2 = (int)round((x_range_xz2+bldg_offset_x)/dx)+D;
  if (i_range_xz2>I_MAX+2 || i_range_xz2<i_range_xz1) {
    cout << "ERROR: Invalid x_range_xz2 value" <<endl;
    exit(1);
  }
  k_range_xz1 = (int)round(z_range_xz1/dx);
  if (k_range_xz1<0) {
    cout << "ERROR: Invalid z_range_xz1 value" <<endl;
    exit(1);
  }
  k_range_xz2 = (int)round(z_range_xz2/dx);
  if (k_range_xz2>K_MAX+2 || k_range_xz2<k_range_xz1) {
    cout << "ERROR: Invalid z_range_xz2 value" <<endl;
    exit(1);
  }
  k_slice_xy = (int)round(z_slice_xy/dx);
  if (k_slice_xy<0 || k_slice_xy>K_MAX-D+1) {
    cout << "ERROR: Invalid z_slice_xy value" <<endl;
    exit(1);
  }
  else if (k_slice_xy==0) k_slice_xy++;
  i_range_xy1 = (int)round((x_range_xy1+bldg_offset_x)/dx)+D;
  if (i_range_xy1<0) {
    cout << "ERROR: Invalid x_range_xy1 value" <<endl;
    exit(1);
  }
  i_range_xy2 = (int)round((x_range_xy2+bldg_offset_x)/dx)+D;
  if (i_range_xy2>I_MAX+2 || i_range_xy2<i_range_xy1) {
    cout << "ERROR: Invalid x_range_xy2 value" <<endl;
    exit(1);
  }
  j_range_xy1 = (int)round((y_range_xy1+bldg_offset_y)/dx)+D;
  if (j_range_xy1<0) {
    cout << "ERROR: Invalid y_range_xz1 value" <<endl;
    exit(1);
  }
  j_range_xy2 = (int)round((y_range_xy2+bldg_offset_y)/dx)+D;
  if (j_range_xy2>J_MAX+2 || j_range_xy2<j_range_xy1) {
    cout << "ERROR: Invalid y_range_xy2 value" <<endl;
    exit(1);
  }

  n_skip_yz = n_max/num_frame_yz; // Number of time steps to skip between each 2D YZ movie frame 
  n_skip_xz = n_max/num_frame_xz; // Number of time steps to skip between each 2D XZ movie frame
  n_skip_xy = n_max/num_frame_xy; // Number of time steps to skip between each 2D XY movie frame
  n_skip_surface = n_max/n_frame_surface; // Number of time steps to skip between each surface movie frame
  
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
  if (n_skip_surface == 0) {
    n_skip_surface = 1;
    n_frame_surface = n_max;
  }

  // Display the domain size 
  long long TotalNumPt = (long long)I_MAX*(long long)J_MAX*(long long)K_MAX; 
  cout << endl;
  if (boom_length_i > input_length_i)
    cout << "WARNING: boom_length_i > input_length_i" << endl; 
  cout << "Using the input file\t\t: " << str_input_name << endl;
  cout << "Domain size in x-direction\t= " << dom_dim_x << " m" << endl;
  cout << "Domain size in y-direction\t= " << dom_dim_y << " m" << endl;
  cout << "Domain size in z-direction\t= " << dom_dim_z << " m" << endl;
  cout << "Grid points in x-direction\t= " << I_MAX << endl;
  cout << "Grid points in y-direction\t= " << J_MAX << endl;
  cout << "Grid points in z-direction\t= " << K_MAX << endl;
  cout << "Total number of grid points\t= " << TotalNumPt << endl;
  
  return;
}
