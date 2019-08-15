/*******************************************************************************
 *
 *  LoBoomFDTD2D_2006.cpp
 *   
 *  written by:     Sang Ik Cho
 *                  stc142@psu.edu 
 *  last modified:  7 / 5 / 11
 *   
 *    This program uses a 2-dimensional FDTD method to simulate low-boom 
 *  propagation and diffraction around the house structure used in 2006 NASA 
 *  flight test.  A second-order leapfrog scheme is applied to a 2-D linear
 *  wave equation.
 *    The house geometry is modeled as a rigid structure with the shape of a 
 *  cross-section of the house.  An experimentally recorded low-boom pressure 
 *  time series from the flight test is taken as an input and is used to 
 *  initialize the computational domain with an incident boom and its reflection
 *  off the rigid ground.  The external boundaries of the computational domain 
 *  is rigid, so the size of the domain must be set with care so as not to avoid
 *  the result being contaminated from unphysical reflections from these 
 *  boundaries.  The ground is modeled as a flat, perfectly rigid surface.
 *    The user specifies a number of simulation parameters including the domain 
 *  size, input boom, grid size, and total simulation time.  "num_floors" is a 
 *  parameter added to specify how tall the building is by specifying how many 
 *  floors the building has.  The height of the house geometry is extended by 
 *  2.5 m per each additional floor.  This enables testing of multi-story 
 *  building cases with the same simulation configurations.  The size of the 
 *  domain and the initial location of the boom wavefront must be carefully set
 *  for different building heights in order to avoid the wavefront touching the
 *  building geometry and to ensure the building is exposed to the full extent 
 *  of the boom energy.  The user can also change the visualization parameters 
 *  for the output file intended for generating animations on Tecplot.   
 *    The program can be run in parallel using OpenMP using a user-specified 
 *  number of threads. The OpenMP #pragma statements are ignored if OpenMP 
 *  option is not enabled at the time of compilation.
 *
 *******************************************************************************/

#include "Headers.h"

string input_name;
int input_length_i, num_floors, x_indmax, z_indmax, t_indmax;
double input_fs, x_max, z_max, x_intersect, x0, dx, tmax, c0, theta_inc, dt, cXcXdtXdt, dxXdx, theta_ref;
int num_frame, mov_skip, x_range_ind_lo, x_range_ind_hi, z_range_ind_hi;
double x_range_lo, x_range_hi, z_range_hi;
double sin_inc, sin_ref, cos_inc, cos_ref, tan14, z1, z2, z3, z4;
int prop_time_inc, prop_time_ref, input_count;
int num_mics;
vector<Mic> mics;

ofstream p_2D_initial, p_domain, p_mic, p_2D_final;

int main(int args, char *argv[])
{
  int i, j, n, l;
  
  cout << endl << "Program running... " << endl;
  /*****************************************************************************
   *  Read Parameters from Input Files
   *****************************************************************************/
  cout << "ReadParameter()... ";
  ReadParameters();
  cout << "DONE" << endl;


  /*****************************************************************************
   *  Initialization of Computational Domain & Output Files    
   *****************************************************************************/
  // Initialize output data files
  cout << "InitializeDataFiles()... ";
  InitializeDataFiles();
  cout << "DONE" << endl;

  
  // Read in the boom signature from the input boom file
  double boom_shape[input_length_i];
  input_count = 0;
  ifstream f_input;
  f_input.open(input_name.c_str());
  if(!f_input) {
    cout << "ERROR: The program could not locate an input file named: " << input_name << endl;
    exit(1);
  } 
  else {
    while(!f_input.eof()) {
      f_input >> scientific >> boom_shape[input_count++];
      if (input_count>input_length_i+1) {
        cout << "ERROR: The input file is longer than specified by input_length." << endl;
        exit(1);
      }
    }
    f_input.close();
  }
  
  // Setup the computational grid
  double x_pos[x_indmax], z_pos[z_indmax];
  double p[x_indmax][z_indmax], p_prev[x_indmax][z_indmax], p_next[x_indmax][z_indmax];
  for ( i = 0; i < x_indmax; i++ )
    x_pos[i] = dx*i;
  for ( j = 0; j < z_indmax; j++ )
    z_pos[j] = dx*j;
    
  // The computational domain should be initialized with the incoming boom for 
  // two time steps, because the second-order leapfrog scheme uses two previous 
  // time steps to calculate the next.
  sin_inc = sin(theta_inc/180*PI);
  cos_inc = cos(theta_inc/180*PI);
  sin_ref = sin(theta_ref/180*PI);
  cos_ref = cos(theta_ref/180*PI);

  #pragma omp parallel for default(shared) private(i,j)
  for ( j = 0; j < z_indmax; j++ ) {
    for ( i = 0; i < x_indmax; i++ ) {
      // Initializing grid at current time step
      prop_time_inc = input_fs*(-cos_inc*x_pos[i]+sin_inc*z_pos[j]+cos_inc*x_intersect)/c0;  // Time for the incident wavefront to get to a grid point, measured in time index
      prop_time_ref = input_fs*(-cos_ref*x_pos[i]+sin_ref*z_pos[j]+cos_ref*x_intersect)/c0;  // Time for the reflected wavefront to get to a grid point, measured in time index
      if ( prop_time_inc >= 0 && prop_time_inc < input_count-1 )
        p[i][j] = boom_shape[prop_time_inc];
      else
        p[i][j] = 0;
      if ( prop_time_ref > 0 && prop_time_ref < input_count-1 )
        p[i][j] = p[i][j]+boom_shape[prop_time_ref];

      // Initializing grid at  previous time step
      prop_time_inc = input_fs*((-cos_inc*x_pos[i]+sin_inc*z_pos[j]+cos_inc*x_intersect)/c0-dt);
      prop_time_ref = input_fs*((-cos_ref*x_pos[i]+sin_ref*z_pos[j]+cos_ref*x_intersect)/c0-dt);
      if ( prop_time_inc >= 0 && prop_time_inc < input_count-1 )
        p_prev[i][j] = boom_shape[prop_time_inc];
      else
        p_prev[i][j] = 0;
      if ( prop_time_ref >= 0 && prop_time_ref < input_count-1 )
        p_prev[i][j] = p_prev[i][j]+boom_shape[prop_time_ref];
      
      // Grid points inside the building are given higher values, only for visualization purpose.
      if ( (x_pos[i] >= x0 && x_pos[i] <= x0+8.8458 && z_pos[j] <= tan14*x_pos[i]+z1 && z_pos[j] <= -tan14*x_pos[i]+z3) || (x_pos[i] >= x0-0.972 && x_pos[i] <= x0+9.818 && z_pos[j] <= tan14*x_pos[i]+z2 && z_pos[j] <= -tan14*x_pos[i]+z4 && (z_pos[j] >= tan14*x_pos[i]+z1 || z_pos[j] >= -tan14*x_pos[i]+z3)) )
        p[i][j] = 100;
    }
  }
  for ( j = 0; j < z_indmax; j++ )
    for ( i = 0; i < x_indmax; i++ )
      p_2D_initial << p[i][j] << endl;


  /*****************************************************************************
   *  FDTD Calculation
   *****************************************************************************/
  // Advancing solution through time steps
  cout << "FDTD calculation running" << endl;
  for ( n = 0; n < t_indmax; n++ ) {
    cout << "Time propagated = " << n*dt << endl; // To observe the simulation progress  
    
    // Rigid boundary condition is applied to the building structure by setting the acoustic 
    // pressure at the adjacent grid points equal. 
    #pragma omp parallel for default(shared) private(i,j)
    for ( j = 0; j < z_indmax; j++ ) {
      for ( i = 0; i < x_indmax; i++ ) {
        if ( i == (int)ceil(x0/dx) && z_pos[j] < tan14*x_pos[i]+z1 )
          p[i][j] = p[i-1][j];
        if ( i == (int)floor((x0+8.8458)/dx) && z_pos[j] < -tan14*x_pos[i]+z3 )
          p[i][j] = p[i+1][j];
        if ( i == (int)ceil((x0-0.9718)/dx) && z_pos[j] < tan14*x_pos[i]+z2 && z_pos[j] > tan14*x_pos[i]+z1 )
          p[i][j] = p[i-1][j];
        if ( i == (int)floor((x0+9.818)/dx) && z_pos[j] < -tan14*x_pos[i]+z4 && z_pos[j] > -tan14*x_pos[i]+z3 )
          p[i][j] = p[i+1][j];
        if ( j == (int)floor((tan14*x_pos[i]+z2)/dx) && x_pos[i] > x0-0.9718 && x_pos[i] < x0+4.4229 )
          p[i][j] = p[i][j+1];
        if ( j == (int)floor((-tan14*x_pos[i]+z4)/dx) && x_pos[i] > x0+4.4229 && x_pos[i] < x0+9.818 )
          p[i][j] = p[i][j+1];          
        if ( j == (int)ceil((tan14*x_pos[i]+z1)/dx) && x_pos[i] > x0-0.9718 && x_pos[i] < x0 )
          p[i][j] = p[i][j-1];
        if ( j == (int)ceil((-tan14*x_pos[i]+z3)/dx) && x_pos[i] > x0+8.846 && x_pos[i] < x0+9.818 )
          p[i][j] = p[i][j-1];
      }
    }

    // Rigid boundary condition for the domain boundary
    for ( i = 1; i < x_indmax-1; i++ ) {
      p[i][z_indmax-1] = p[i][z_indmax-2];
      p[i][0] = p[i][1];
    }
    for ( j = 0; j < z_indmax; j++ ) {
      p[x_indmax-1][j] = p[x_indmax-2][j];
      p[0][j] = p[1][j];
    }
    
    // Finite difference calculation
    #pragma omp parallel for default(shared) private(i,j)
    for ( j = 1; j < z_indmax-1; j++ ) {
      for ( i = 1; i < x_indmax-1; i++ ) {
        // Interior points are skipped
        if ( x_pos[i] >= x0 && x_pos[i] <= x0+8.8458 && z_pos[j] <= tan14*x_pos[i]+z1 && z_pos[j] <= -tan14*x_pos[i]+z3 )
          p_next[i][j] = 0;
        else if ( x_pos[i] >= x0-0.972 && x_pos[i] <= x0+9.818 && z_pos[j] <= tan14*x_pos[i]+z2 && z_pos[j] <= -tan14*x_pos[i]+z4 && (z_pos[j] >= tan14*x_pos[i]+z1 || z_pos[j] >= -tan14*x_pos[i]+z3) )
          p_next[i][j] = 0;          
        // FDTD calculation is performed at exterior points
        else
  	      p_next[i][j] = 2*p[i][j]-p_prev[i][j]+cXcXdtXdt*((p[i+1][j]-2*p[i][j]+p[i-1][j])/dxXdx+(p[i][j+1]-2*p[i][j]+p[i][j-1])/dxXdx);
  	  }
    }

    // Advance in time and store pressure data
    for ( j = 0; j < z_indmax; j++ ) {
      for ( i = 0; i < x_indmax; i++ ) {
        p_prev[i][j] = p[i][j];
        p[i][j] = p_next[i][j];
        if ( (x_pos[i] >= x0 && x_pos[i] <= x0+8.8458 && z_pos[j] <= tan14*x_pos[i]+z1 && z_pos[j] <= -tan14*x_pos[i]+z3) || (x_pos[i] >= x0-0.972 && x_pos[i] <= x0+9.818 && z_pos[j] <= tan14*x_pos[i]+z2 && z_pos[j] <= -tan14*x_pos[i]+z4 && (z_pos[j] >= tan14*x_pos[i]+z1 || z_pos[j] >= -tan14*x_pos[i]+z3)) )   // To display the geometry of the structure
          p[i][j] = 100;
        if ( (n%mov_skip) == 0 && i >= x_range_ind_lo && i < x_range_ind_hi && j < z_range_ind_hi )
          p_domain << p[i][j] << endl;
      }
    }
    // Time history of pressure recorded at microphone locations specified by user
    p_mic << n*dt << " ";
    for (l=0; l<num_mics; l++)
      p_mic << p[mics[l].i][mics[l].j]<< " ";
    p_mic << endl;
  }

  // Plot final state of the domain to check for contamination of the pressure field from unphysical reflections
  for ( j = 0; j < z_indmax; j++ )
    for ( i = 0; i < x_indmax; i++ )
      p_2D_final << p[i][j] << endl;

  p_2D_initial.close();
  p_domain.close();
  p_mic.close();
  p_2D_final.close();

  cout << "\n***********************************************\n\n";
  cout << "\tThis program terminated normally.\n";
  cout << "\n***********************************************\n";

  return 0;
}