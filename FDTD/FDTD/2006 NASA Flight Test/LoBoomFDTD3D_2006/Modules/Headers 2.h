/******************************************************************************
 *
 *  Headers.h
 *
 *  written by:     Sang Ik Cho
 *                  stc142@psu.edu 
 *  last modified:  5 / 1 / 2012
 *     
 *    This headear module contains all of the library headers, variable type 
 *  definition, and declaration of global variables, file input/output streams, 
 *  and function prototypes for LoBoomFDTD3D_MultiBldg program.  This header 
 *  must be included in the main function and every member modules for correct 
 *  compilation of the program.
 *  
 *  *** This program is written for the 2006 test house geometry ***
 *    
 ******************************************************************************/  

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <omp.h>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

#define PI 3.14159265

// Definition of Custom Variable Types
struct GridPoint {
  bool InOrOut;
  float px;
  float py;
  float pz;
  float u;
  float v;
  float w;
};

struct Mic {
  int i,j,k;
};

typedef vector<vector<vector<GridPoint> > > Grid3D;
typedef vector<float> Array;

// Variables
extern string str_input_name;
extern float fs,rho,c0,z0,dx,CFL,dt,t_max,sigma_m,B,angle_elv,theta;
extern int input_length_i,boom_length_i,D,n_max;
extern int i_intersect,I_MAX,J_MAX,K_MAX;
extern int num_frame_yz,num_frame_xz,num_frame_xy,n_skip_yz,n_skip_xz,n_skip_xy;
extern float x_slice_yz,y_range_yz1,y_range_yz2,z_range_yz1,z_range_yz2;
extern float y_slice_xz,x_range_xz1,x_range_xz2,z_range_xz1,z_range_xz2;
extern float z_slice_xy,x_range_xy1,x_range_xy2,y_range_xy1,y_range_xy2;
extern int i_slice_yz,j_range_yz1,j_range_yz2,k_range_yz1,k_range_yz2;
extern int j_slice_xz,i_range_xz1,i_range_xz2,k_range_xz1,k_range_xz2;
extern int k_slice_xy,i_range_xy1,i_range_xy2,j_range_xy1,j_range_xy2;
extern float extra_dim_x,extra_dim_y,buffer_dim_x,buffer_dim_z;
extern int num_mics;
extern vector<Mic> mics;

// File I/O stream
extern ofstream p_3D_init, p_2DYZ, p_2DXZ, p_2DXY, p_mic;

// Function Prototypes
extern void ReadParameters();
extern void PrepareDomain(Grid3D &);
extern void InitializeDataFiles();
extern void InitializeMatrix(Grid3D &, Array &, Array &);
extern void EulerPML3D(Grid3D &, Grid3D &, Array &);
extern void VisualizeMatrix(Grid3D &, int);