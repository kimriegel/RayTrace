/*******************************************************************************
 *
 *  LoBoomFDTD3D_SingleBldg.cpp
 *  
 *  written by:     Sang Ik Cho
 *                  stc142@psu.edu 
 *  last modified:  4 / 30 / 2012
 *  
 *    This program simulates low-boom propagation and diffraction around a
 *  single rectangular building structure using a 3-dimensional FDTD method.  
 *  A simple second-order centered-difference scheme in spatial domain and a 
 *  4-stage low-dispersion low-dissipation Runge-Kutta time integration method 
 *  were applied to the linearized Euler equations coupling acoustic pressure 
 *  and particle velocity.
 *    The orientation of the computational domain is aligned with the building 
 *  such that the two side walls are parallel with the x-axis and the incident &
 *  back walls are parallel with the y-axis.  An experimentally recorded 
 *  low-boom pressure time series is taken as an input and is used to initialize
 *  the computational domain with an incident boom and its reflection according 
 *  to its incoming azimuth and elevation angles.  The size of the computational
 *  domain is automatically calculated based on the location of the building and 
 *  the incoming angle of the boom.  A simple rigid boundary condition is used 
 *  for the ground and building and a perfectly matched layer (PML) absorbing 
 *  boundary is used at the other boundaries of the domain to let the out-going 
 *  wave energy radiate without unphysical reflections.
 *    The user specifies a number of simulation parameters including the input 
 *  boom, grid size, total simulation time, PML parameters, building dimensions,
 *  and parameters for determining the domain size.  Although the calculation of
 *  the domain size is automated, the "buffer" in each coordinate direction 
 *  should be carefully set to ensure the building is exposed to the
 *  undistorted, full extent of the incoming boom.  The user can also change the
 *  visualization parameters for the output files intended for generating 
 *  animations on Tecplot.  The number of microphones to record the pressure
 *  time history and their locations are also set by the user. 
 *    The program can be run in parallel using OpenMP using a user-specified 
 *  number of threads. The OpenMP #pragma statements are ignored if OpenMP 
 *  option is not enabled at the time of compilation.
 *   
 *******************************************************************************/

#include "Headers.h"

string str_input_name;
float fs,rho,c0,z0,dx,CFL,dt,t_max,sigma_m,B,angle_elv,theta,angle_azm,phi;
int input_length_i,boom_length_i,D,n_max;
float bldg_dim_x,bldg_dim_y,bldg_dim_z,bldg_offset_x,bldg_offset_y;
float buffer_dim_x,buffer_dim_y,buffer_dim_z;
float boom_duration,x_intersect;
float dom_dim_x,dom_dim_y,dom_dim_z;
int i_intersect,I_MAX,J_MAX,K_MAX;
int bldg_dim_i,bldg_dim_j,bldg_dim_k,bldg_offset_i,bldg_offset_j;
int wall_loc_i1,wall_loc_i2,wall_loc_j1,wall_loc_j2,wall_loc_k1,wall_loc_k2;
int num_frame_yz,num_frame_xz,num_frame_xy,n_skip_yz,n_skip_xz,n_skip_xy;
float x_slice_yz,y_range_yz1,y_range_yz2,z_range_yz1,z_range_yz2;
float y_slice_xz,x_range_xz1,x_range_xz2,z_range_xz1,z_range_xz2;
float z_slice_xy,x_range_xy1,x_range_xy2,y_range_xy1,y_range_xy2;
int i_slice_yz,j_range_yz1,j_range_yz2,k_range_yz1,k_range_yz2;
int j_slice_xz,i_range_xz1,i_range_xz2,k_range_xz1,k_range_xz2;
int k_slice_xy,i_range_xy1,i_range_xy2,j_range_xy1,j_range_xy2;
int n_frame_surface,n_skip_surface,num_mics;
Mic* mics;

ofstream p_3D_init, p_2DYZ, p_2DXZ, p_2DXY, p_mic, p_surface;

int main(int args, char *argv[])
{  
  int i,j,k;
  
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
  
  // Declare & initialize arrays & matrices
  cout << "Declaring Arrays & Matrices... ";
  Grid3D U;
  U.resize(I_MAX+2);
  for (i=0; i<I_MAX+2; i++) {
    U[i].resize(J_MAX+2);
    for (j=0; j<J_MAX+2; j++) {
      U[i][j].resize(K_MAX+2);
    }
  }
  Grid3D Um(U),Km(U);
  Array sigma(D), boom_shape(input_length_i);
  cout << "DONE" << endl;
  
  cout << "InitializeMatrix()... ";
  InitializeMatrix(U, sigma, boom_shape);
  cout << "DONE" << endl << endl;


  /*****************************************************************************
   *  FDTD Calculation
   *****************************************************************************/
    // Coefficients for LDDRK4 scheme
    double c_LDDRK4[4] = {1, 0.5, 0.162997, 0.0407574};
    double beta_LDDRK4[4] = {0, c_LDDRK4[3]/c_LDDRK4[2], c_LDDRK4[2]/c_LDDRK4[1], c_LDDRK4[1]};

  cout << "FDTD calculation running" << endl;
  // The centered-difference spatial finite differencing module that includes
  // the rigid boundary condition for the ground and building is called for each
  // level of Runge Kutta calculation.
  for (int n=0; n<n_max; n++) {
    cout << "Time propagated = " << n*dt << endl;
    #pragma omp parallel default(shared) private(i,j,k)
    {
      for (int m=0; m<4; m++) {
        #pragma omp for
        for (i=1; i<I_MAX+1; i++) {
          for (j=1; j<J_MAX+1; j++) {
            for (k=1; k<K_MAX+1; k++) {
              Um[i][j][k].px = U[i][j][k].px+beta_LDDRK4[m]*Km[i][j][k].px;
              Um[i][j][k].py = U[i][j][k].py+beta_LDDRK4[m]*Km[i][j][k].py;
              Um[i][j][k].pz = U[i][j][k].pz+beta_LDDRK4[m]*Km[i][j][k].pz;
              Um[i][j][k].u = U[i][j][k].u+beta_LDDRK4[m]*Km[i][j][k].u;
              Um[i][j][k].v = U[i][j][k].v+beta_LDDRK4[m]*Km[i][j][k].v;
              Um[i][j][k].w = U[i][j][k].w+beta_LDDRK4[m]*Km[i][j][k].w;
            }
          }
        }
        EulerPML3D(Um,Km,sigma);
      }
      #pragma omp for
      for (i=1; i<I_MAX+1; i++) {
        for (j=1; j<J_MAX+1; j++) {
          for (k=1; k<K_MAX+1; k++) {
            U[i][j][k].px = U[i][j][k].px+Km[i][j][k].px;
            U[i][j][k].py = U[i][j][k].py+Km[i][j][k].py;
            U[i][j][k].pz = U[i][j][k].pz+Km[i][j][k].pz;
            U[i][j][k].u = U[i][j][k].u+Km[i][j][k].u;
            U[i][j][k].v = U[i][j][k].v+Km[i][j][k].v;
            U[i][j][k].w = U[i][j][k].w+Km[i][j][k].w;
          }
        }
      }        
    }
    
    // Pressure data is stored at the end of every time step
    VisualizeMatrix(U,n);
  }

  p_2DYZ.close();
  p_2DXZ.close();
  p_2DXY.close();
  p_mic.close();
  p_surface.close();
  
  cout << endl << "***********************************************" << endl;
  cout << endl << "\tThis program terminated normally." << endl;
  cout << endl << "***********************************************" << endl;

  return 0;
}
